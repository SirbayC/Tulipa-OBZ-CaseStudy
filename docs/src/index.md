# OBZ case study

Tutorial for the OBZ case study as an example of the full workflow of Tulipa.

We are basing ourselves on the Tulipa [data pipeline/workflow](#TODO).
To help us navigate this workflow, we'll reproduce the diagram from the link above here.
For more details on the steps of the workflow, check the original link, or follow the tutorial.

![Tulipa Workflow. Textual explanation below.](./assets/tulipa-workflow.jpg)

## External source

We begin by inspecting our external data source.
This data is specific to OBZ, but it might be useful.

TODO: Expand this exploration.

```@example obz
user_input_dir = joinpath(@__DIR__, "..", "..", "user-input-files")
readdir(user_input_dir)
```

For the Tulipa workflow, we will need to transform some of this data into a specific format.
This can be done externally in whatever tools you are already comfortable with,
or through Julia via DuckDB and Tulipa's convenience functions.

## Create connection

Once we are done manipulating the data externally, it is time to create a DuckDB connection.

You can create a connection storing the DB locally, or keep everything in-memory only.
Let's assume you want to store the DB, otherwise you can just remove this argument.

```@example obz
using DuckDB: DBInterface, DuckDB

connection = DBInterface.connect(DuckDB.DB, "obz.db")
# Manually cleaning everything
for row in DuckDB.query(connection, "SHOW TABLES")
    try
        DuckDB.query(connection, "DROP TABLE $(row.name)")
    catch
        DuckDB.query(connection, "DROP VIEW $(row.name)")
    end
end
```

We will be performing various queries with DuckDB. To format them nicely, we can wrap the results in a `DataFrame`:

```@example obz
using DataFrames: DataFrame

nice_query(str) = DataFrame(DuckDB.query(connection, str))
```

## Load data

Once we are done manipulating the data externally, it is time to load it into the DuckDB connection.
We can load them manually with `DuckDB`, but we also have a convenience function:

```@example obz
using TulipaIO: TulipaIO

TulipaIO.read_csv_folder(
    connection, 
    user_input_dir,
    table_name_prefix = "raw_", 
    replace_if_exists = true,
    skip = 1, # Data in this folder has an extra header line
)
TulipaIO.read_csv_folder(
    connection, 
    joinpath(user_input_dir, "profiles"), 
    table_name_prefix = "raw_", 
    replace_if_exists = true,
)

nice_query("SHOW TABLES")
```

## Data processing for scenarios with DuckDB/TulipaIO

As we mentioned before, you can process your data externally and then load it.
But you can also use Julia and DuckDB to process the data.

This step is required to prepare the data for TulipaClustering for the
clustering of the profile data.

We need a single profiles table with 4 columns:

- `profile_name`
- `year`
- `timestep`
- `value`

And we have 2 profiles tables:

```@example obz
nice_query("FROM raw_profiles LIMIT 10")
```

```@example obz
nice_query("FROM raw_min_max_reservoir_levels LIMIT 10")
```

Notice that these are all hourly profiles for the whole year:

```@example obz
nice_query("SELECT year, MAX(timestep) FROM raw_profiles GROUP BY year")
```

```@example obz
nice_query("SELECT year, MAX(timestep) FROM raw_min_max_reservoir_levels GROUP BY year")
```

So we will transform both these tables to long format and stack them:

```@example obz
using TulipaClustering: TulipaClustering 

TulipaClustering.transform_wide_to_long!(connection, "raw_profiles", "pivot_profiles")
TulipaClustering.transform_wide_to_long!(connection, "raw_min_max_reservoir_levels", "pivot_min_max_reservoir_levels")

DuckDB.query(
    connection,
    "CREATE OR REPLACE TABLE input_profiles AS
    FROM pivot_profiles
    UNION
    FROM pivot_min_max_reservoir_levels
    ORDER BY profile_name, year, timestep
    "
)

nice_query("SELECT COUNT(*) FROM input_profiles")
```

The expected number of rows is

```@example obz
8760 * (165 - 2 + 14 - 2)
```

We can get some specific profile the following way:

```@example obz
using Plots

profile_name = "NL_Solar"
year = 2050

# Notice that we don't use DataFrame, because we don't need the conversion, and
# it's cheaper to not do it.
subtable = DuckDB.query(
    connection, 
    "SELECT timestep, value FROM input_profiles 
    WHERE profile_name='$profile_name' AND year=$year
    ORDER BY timestep
    ",
)
df = DataFrame(subtable)
plot(df.timestep, df.value, lab=profile_name)
```

## Cluster into representative periods using TulipaClustering

```@example obz
using Distances: SqEuclidean

clustering_params = (
    ## Data for clustering
    n_rp = 3,                         # number of representative periods
    period_duration = div(8760, 365), # hours of the representative period
    method = :k_means,
    distance = SqEuclidean(),
    ## Data for weight fitting
    weight_type = :convex,
    tol = 1e-2,
    ## Data for projected subgradient
    niters = 100,
    learning_rate = 0.001,
    adaptive_grad = false,
)

TulipaClustering.cluster!(
    connection,
    clustering_params.period_duration,
    clustering_params.n_rp;
    clustering_params.method,
    clustering_params.distance,
    clustering_params.weight_type,
    clustering_params.tol,
    clustering_params.niters,
    clustering_params.learning_rate,
    clustering_params.adaptive_grad,
)
```

## Prepare data for TulipaEnergyModel's format

Now the hard part starts. We need to create several files for Tulipa using our files.
Again, we remind you that you can create most of these files externally, i.e.,
you don't have to use DuckDB to join them here.

### Assets

First, let's join all basic assets data. I will use code to generate the list of relevant tables, but you can do it manually too.

```@example obz
table_names = [
    row.table_name
    for row in DuckDB.query(connection, "FROM duckdb_tables WHERE table_name LIKE 'raw_assets_%_basic_data'")
]
```

Now, we stack all these tables.

```@example obz
query_string = join(
    ["FROM $t" for t in table_names],
    " UNION BY NAME ",
)
DuckDB.query(
    connection, 
    "CREATE TABLE input_asset AS 
    SELECT
        name AS asset,
        type,
        COALESCE(capacity, 0) AS capacity,
        COALESCE(capacity_storage_energy, 0) AS capacity_storage_energy,
        COALESCE(is_seasonal, false) AS is_seasonal,
    FROM ($query_string)
    ORDER BY asset
    ",
)

nice_query("FROM input_asset ORDER BY random() LIMIT 10")
```

Similarly, we join the `raw_assets_%_yearly_data` tables:

```@example obz
table_names = [
    row.table_name
    for row in DuckDB.query(connection, "FROM duckdb_tables WHERE table_name LIKE 'raw_assets_%_yearly_data'")
]

query_string = join(
    ["FROM $t" for t in table_names],
    " UNION BY NAME ",
)

DuckDB.query(
    connection, 
    "CREATE TABLE t_asset_yearly AS 
    FROM ($query_string)
    ",
)

nice_query("FROM t_asset_yearly ORDER BY random() LIMIT 10")
```

This `t_asset_yearly` tables is used to create three other asset tables.

```@example obz
DuckDB.query(
    connection,
    "CREATE TABLE input_asset_commission AS 
    SELECT
        name AS asset,
        year AS commission_year,
    FROM t_asset_yearly
    ORDER by asset
    "
)

DuckDB.query(
    connection,
    "CREATE TABLE input_asset_milestone AS 
    SELECT
        name AS asset,
        year AS milestone_year,
        coalesce(0, peak_demand) AS peak_demand,
        coalesce(0, initial_storage_level) AS initial_storage_level,
        coalesce(0, storage_inflows) AS storage_inflows,
    FROM t_asset_yearly
    ORDER by asset
    "
)

DuckDB.query(
    connection,
    "CREATE TABLE input_asset_both AS 
    SELECT
        name AS asset,
        year AS milestone_year,
        year AS commission_year,
        coalesce(0, initial_units) AS initial_units,
        coalesce(0, initial_storage_units) AS initial_storage_units,
    FROM t_asset_yearly
    ORDER by asset
    "
)
```

### Flows

We repeat the steps above for flows:

```@example obz
table_names = [
    row.table_name
    for row in DuckDB.query(connection, "FROM duckdb_tables WHERE table_name LIKE 'raw_flows_%_basic_data'")
]
query_string = join(
    ["FROM $t" for t in table_names],
    " UNION BY NAME ",
)
DuckDB.query(
    connection, 
    "CREATE TABLE input_flow AS 
    SELECT
        from_asset,
        to_asset,
        carrier,
        COALESCE(capacity, 0) AS capacity,
        COALESCE(is_transport, false) AS is_transport,
    FROM ($query_string)
    ORDER BY from_asset, to_asset
    ",
)

table_names = [
    row.table_name
    for row in DuckDB.query(connection, "FROM duckdb_tables WHERE table_name LIKE 'raw_flows_%_yearly_data'")
]

query_string = join(
    ["FROM $t" for t in table_names],
    " UNION BY NAME ",
)

DuckDB.query(
    connection, 
    "CREATE TABLE t_flow_yearly AS 
    FROM ($query_string)
    ",
)

DuckDB.query(
    connection,
    "CREATE TABLE input_flow_commission AS 
    SELECT
        from_asset,
        to_asset,
        year AS commission_year,
        coalesce(0, efficiency) AS efficiency,
    FROM t_flow_yearly
    ORDER by from_asset, to_asset
    "
)

DuckDB.query(
    connection,
    "CREATE TABLE input_flow_milestone AS 
    SELECT
        from_asset, 
        to_asset,
        year AS milestone_year,
        coalesce(0, variable_cost) AS variable_cost,
    FROM t_flow_yearly
    ORDER by from_asset, to_asset
    "
)

DuckDB.query(
    connection,
    "CREATE TABLE input_flow_both AS 
    SELECT
        from_asset, 
        to_asset,
        year AS milestone_year,
        year AS commission_year,
        coalesce(0, initial_export_units) AS initial_export_units,
        coalesce(0, initial_import_units) AS initial_import_units,
    FROM t_flow_yearly
    ORDER by from_asset, to_asset
    "
)
```

```@example obz
1
```

```@example obz
# DEBUG
using TulipaEnergyModel: TulipaEnergyModel as TEM

expected_folder = joinpath(@__DIR__, "..", "..", "tulipa-energy-model-files") 
TulipaIO.read_csv_folder(
    connection, 
    expected_folder, 
    schemas = TEM.schema_per_table_name, 
    table_name_prefix = "expected_",
)
nice_query("FROM expected_assets_profiles ORDER BY asset LIMIT 10")

# Copy remaining tables
copied_over = String[]
for row in DuckDB.query(connection, "FROM duckdb_tables WHERE table_name LIKE 'expected_%'")
    expected_table_name = row.table_name
    table_name = replace(expected_table_name, "expected_" => "")
    if table_name in ("rep_periods_data", "rep_periods_mapping", "profiles_rep_periods")
        continue
    end
    input_table_name = "input_" * table_name
    exitflag = DuckDB.query(
        connection, 
        "CREATE TABLE IF NOT EXISTS $input_table_name AS FROM $expected_table_name",
    )
    if length(collect(exitflag)) > 0
        push!(copied_over, table_name)
    end
end
copied_over
```

```@example obz
TEM.populate_with_defaults!(connection)
```

```@example obz
try
    energy_problem = TEM.EnergyProblem(connection)
catch ex
    @show ex
end
# nice_query("SHOW TABLES")
```

```@example obz
# TEM.create_model!(energy_problem)
nice_query("from input_assets_rep_periods_partitions LIMIT 10")
```

## Create internal tables for the model indices

## Create model

## Solve model

## Store primal and dual solution

## Data processing for plots and dashboard

## Create plots

## Export solution

## [TODO](@id TODO)

- [ ] Link to Tulipa data pipeline

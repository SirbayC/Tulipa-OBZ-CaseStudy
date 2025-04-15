# TODO/Abel: Transform into documentation
# Step 0: Preparation
@info "Activating the environment"
cd(@__DIR__)
using Pkg: Pkg
Pkg.activate(".")
#Pkg.update()
Pkg.instantiate()

@info "Loading the packages"
# Data
using CSV
using DataFrames
using DuckDB

# Plots
using Plots
using Plots.PlotMeasures
using StatsPlots

# Misc
using Glob
using Statistics

# Tulipa
using Distances
using HiGHS
using JuMP
using TulipaClustering
using TulipaEnergyModel
using TulipaIO

# TODO/Abel: review these functions
@info "Including the helper functions"
include("utils/functions.jl")

@info "Defining the directories"
user_input_dir = joinpath(@__DIR__, "user-input-files")
# TODO/Abel: Shouldn't need tulipa_files_dir
tulipa_files_dir = joinpath(@__DIR__, "tulipa-energy-model-files")
output_dir = joinpath(@__DIR__, "outputs")

## Create or reset the tulipa files and outputs
if isdir(tulipa_files_dir)
end
rm(tulipa_files_dir; force = true, recursive = true)
mkdir(tulipa_files_dir)

if !isdir(output_dir)
    mkdir(output_dir)
end

@info "Defining default values"
# TODO/Abel: delete this
default_values = get_default_values(; default_year = 2050)

@info "Defining TulipaClustering data"
clustering_params = (
    ## Data for clustering
    n_rp = 3,               # number of representative periods
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

@info "Pre-processing the profiles"
include("utils/preprocess-profiles.jl")

@info "Pre-processing the user inputs"
include("utils/preprocess-user-inputs.jl")

@info "Setting up the solver"
optimizer = HiGHS.Optimizer
optimizer_parameters =
    Dict("output_flag" => true, "mip_rel_gap" => 0.0, "mip_feasibility_tolerance" => 1e-5)

# using Gurobi
# optimizer = Gurobi.Optimizer
# optimizer_parameters = Dict("OutputFlag" => 1, "MIPGap" => 0.0, "FeasibilityTol" => 1e-5)

@info "Reading Tulipa input files"
connection = DBInterface.connect(DuckDB.DB, "obz-script.db")
# Manually cleaning everything
for row in DuckDB.query(connection, "show tables")
    try
        DuckDB.query(connection, "DROP TABLE $(row.name)")
    catch
        DuckDB.query(connection, "DROP VIEW $(row.name)")
    end
end
read_csv_folder(
    connection,
    tulipa_files_dir;
    schemas = TulipaEnergyModel.schema_per_table_name,
    replace_if_exists = true,
)

@info "Solving the optimization problem"
error("NO NEED")
energy_problem = run_scenario(
    connection;
    output_folder = output_dir,
    optimizer = optimizer,
    optimizer_parameters = optimizer_parameters,
    #model_file_name = "model.lp",
    show_log = false,
    log_file = "log_file.log",
    #enable_names = false,
)

if energy_problem.termination_status == INFEASIBLE
    @info "Getting the IIS model to find infeasibilities"
    compute_conflict!(energy_problem.model)
    iis_model, reference_map = copy_conflict(energy_problem.model)
    print(iis_model)
end

@info "Creating auxiliary file to post-process the results"
assets_bidding_zone_tecnology_file = "assets-bidding-zone-tecnology-data.csv"
df_assets_basic_data = create_one_file_for_assets_basic_info(
    assets_bidding_zone_tecnology_file,
    user_input_dir,
    output_dir,
    default_values,
)

@info "Post-processig of prices"
prices = get_prices_dataframe(connection, energy_problem)

@info "Post-processig of storage levels"
intra_storage_levels = get_intra_storage_levels_dataframe(connection)

@info "Post-processig of balance per balancing asset"
balances = get_balance_per_asset(connection, energy_problem, df_assets_basic_data)

@info "Saving the results to CSV files"
prices_file_name = joinpath(output_dir, "eu-case-prices.csv")
CSV.write(prices_file_name, unstack(prices, :asset, :price))

intra_storage_levels_file_name = joinpath(output_dir, "eu-case-intra-storage-levels.csv")
CSV.write(intra_storage_levels_file_name, unstack(intra_storage_levels, :asset, :SoC))

balance_file_name = joinpath(output_dir, "eu-case-balance-per-bidding-zone.csv")
CSV.write(balance_file_name, unstack(balances, :technology, :solution; fill = 0))

@info "Plotting the results"
prices_plot = plot_prices(
    prices;
    assets = ["NL_E_Balance", "UK_E_Balance", "OBZLL_E_Balance"],
    #rep_periods = [1, 2],
    #plots_args = (xlims = (8760 / 2, 8760 / 2 + 168), ylims = (0, 100)),
    plots_args = (xticks = 0:730:8760, ylim = (0, 100)),
    duration_curve = true,
)
prices_plot_name = joinpath(output_dir, "eu-case-price-duration-curve.png")
savefig(prices_plot, prices_plot_name)

batteries_storage_levels_plot = plot_intra_storage_levels(
    intra_storage_levels;
    assets = ["NL_Battery", "UK_Battery"],
    #rep_periods = [1, 2],
    plots_args = (xlims = (8760 / 2, 8760 / 2 + 168), xticks = 0:12:8760, ylims = (0, 1)),
)
batteries_storage_levels_plot_name = joinpath(output_dir, "eu-case-batteries-storage-levels.png")
savefig(batteries_storage_levels_plot, batteries_storage_levels_plot_name)

if n_rp > 1
    hydro_storage_levels_plot = plot_inter_storage_levels(
        connection;
        assets = ["ES_Hydro_Reservoir", "NO_Hydro_Reservoir", "FR_Hydro_Reservoir"],
        #plots_args = (xticks = 0:730:8760, ylims = (0, 1)),
    )
else
    hydro_storage_levels_plot = plot_intra_storage_levels(
        intra_storage_levels;
        assets = ["ES_Hydro_Reservoir", "NO_Hydro_Reservoir", "FR_Hydro_Reservoir"],
        #rep_periods = [1, 2],
        plots_args = (xticks = 0:730:8760, ylims = (0, 1)),
    )
end
hydro_storage_levels_plot_name = joinpath(output_dir, "eu-case-hydro-storage-levels.png")
savefig(hydro_storage_levels_plot, hydro_storage_levels_plot_name)

asset = "NL_E_Balance" # Any hub or consumer is a valid assets
balance_plot = plot_asset_balance(
    balances;
    asset = asset,
    year = 2050,
    rep_period = 1,
    plots_args = (xlims = (8760 / 2, 8760 / 2 + 168), xticks = 0:6:8760),
)
balance_plot_name = joinpath(output_dir, "eu-case-balance-$asset.png")
savefig(balance_plot, balance_plot_name)

from_asset = "OBZLL_E_Balance"
to_asset = "NL_E_Balance"
year = 2050
rep_period = 1
flow_plot = plot_flow(
    connection,
    from_asset,
    to_asset,
    year,
    rep_period;
    plots_args = (xlims = (8760 / 2, 8760 / 2 + 168), xticks = 0:12:8760),
)
flow_plot_name = joinpath(output_dir, "flows-$from_asset-$to_asset.png")
savefig(flow_plot, flow_plot_name)

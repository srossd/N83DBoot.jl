"""
CSV caching system and Mathematica interface for functional computations.

Provides:
- Filename generation based on functional parameters
- Interface to call external Mathematica scripts for computing functionals
- CSV file reading/writing with BigFloat support
- Utilities for managing cached functionals

## Cache directory layout

```
cache_dir/
├── Crossing/
│   └── cross_1_m_4_n_0_precision_1024_rOrder_100/
│       ├── .cache_metadata          # creation timestamp
│       ├── identity.csv             # single-row: (value)
│       ├── Bp0020.csv               # single-row: (value)
│       ├── Bp0040.csv               # single-row: (value)
│       ├── Btwo0200.csv             # single-row: (value)
│       ├── semishort_Ap0020_J_0.csv # single-row: (value)
│       ├── semishort_Ap0020_J_2.csv
│       ├── semishort_Atwo0100_J_1.csv
│       ├── spin_0.csv               # (tau, value) rows for LongOperator J=0
│       ├── spin_2.csv
│       └── ...
└── Integral/
    └── b_0.0_precision_256_rOrder_20/
        ├── .cache_metadata
        ├── identity.csv
        ├── ...
        └── spin_0.csv
```
"""

using CSV
using DataFrames
using Dates
using Downloads
using JLD2
using SHA

# Path to the Mathematica scripts bundled with this package
const _SCRIPTS_DIR = joinpath(pkgdir(@__MODULE__), "scripts")

# ============================================================================
# Functional type name dispatch
# ============================================================================

"""
    functional_type_name(functional::Functional) -> String

Return the canonical type-name string used for cache subdirectories.
"""
functional_type_name(::CrossingFunctional) = "Crossing"
functional_type_name(::IntegralFunctional) = "Integral"

# ============================================================================
# Cache filename templates
# ============================================================================

"""
    get_cache_filename(functional::Functional) -> String

Return the cache subdirectory name for a given functional.
"""
function get_cache_filename(f::CrossingFunctional)
    return "cross_$(f.cross)_m_$(f.m)_n_$(f.n)_precision_$(f.precision)_rOrder_$(f.r_order)"
end

function get_cache_filename(f::IntegralFunctional)
    return "b_$(f.b)_precision_$(f.precision)_rOrder_$(f.r_order)"
end

# ============================================================================
# Operator filename dispatch
# ============================================================================

"""
    operator_cache_filename(operator::Operator) -> String

Return the CSV filename (without path) for a given operator's cache entry.
"""
operator_cache_filename(::IdentityOperator)    = "identity.csv"
operator_cache_filename(::Bp0020Operator)      = "Bp0020.csv"
operator_cache_filename(::Bp0040Operator)      = "Bp0040.csv"
operator_cache_filename(::Btwo0200Operator)    = "Btwo0200.csv"
operator_cache_filename(op::SemishortAp0020)   = "semishort_Ap0020_J_$(op.J).csv"
operator_cache_filename(op::SemishortAtwo0100) = "semishort_Atwo0100_J_$(op.J).csv"
operator_cache_filename(op::LongOperator)      = "spin_$(op.J).csv"

# ============================================================================
# Cache path helpers
# ============================================================================

"""
    get_cache_path(functional::Functional) -> String

Return the full path to the cache directory for `functional`.
"""
function get_cache_path(functional::Functional)
    subdir   = functional_type_name(functional)
    filename = get_cache_filename(functional)
    return joinpath(get_cache_dir(), subdir, filename)
end

# ============================================================================
# Cache directory management
# ============================================================================

"""
    ensure_cache_dir(functional::Functional)

Create the cache subdirectory for `functional` if it does not exist,
and write a metadata timestamp file.
"""
function ensure_cache_dir(functional::Functional)
    full_dir = joinpath(get_cache_dir(), functional_type_name(functional))
    mkpath(full_dir)

    func_dir = get_cache_path(functional)
    if !isdir(func_dir)
        mkpath(func_dir)
        write_cache_metadata(func_dir)
    end
end

# ============================================================================
# Cache metadata helpers
# ============================================================================

const CACHE_METADATA_FILENAME = ".cache_metadata"

"""
    write_cache_metadata(dir::String)

Write a `.cache_metadata` file in `dir` containing the current UTC timestamp.
"""
function write_cache_metadata(dir::String)
    open(joinpath(dir, CACHE_METADATA_FILENAME), "w") do f
        write(f, string(now(UTC)))
    end
end

"""
    read_cache_metadata(dir::String) -> Union{DateTime,Nothing}

Read the UTC timestamp from a `.cache_metadata` file, or `nothing` on failure.
"""
function read_cache_metadata(dir::String)
    path = joinpath(dir, CACHE_METADATA_FILENAME)
    isfile(path) || return nothing
    try
        return DateTime(strip(read(path, String)))
    catch
        return nothing
    end
end

# ============================================================================
# Mathematica parameter builders
# ============================================================================

"""
    build_mathematica_params(functional::Functional) -> String

Return a Mathematica-style parameter string for the functional's fields,
including the `type` key consumed by `compute_functional.m`.
"""
function build_mathematica_params(functional::CrossingFunctional)
    return "type -> \"Crossing\", cross -> $(functional.cross), m -> $(functional.m), n -> $(functional.n), precision -> $(functional.precision), rOrder -> $(functional.r_order)"
end

function build_mathematica_params(functional::IntegralFunctional)
    return "type -> \"Integral\", b -> $(functional.b), precision -> $(functional.precision), rOrder -> $(functional.r_order)"
end

# ============================================================================
# Operator argument builders
# ============================================================================

"""Build the Mathematica `operators -> {...}` argument for fixed operators."""
function _build_operators_arg(fixed_ops::Vector)
    names = [_fixed_op_mathematica_name(op) for op in fixed_ops]
    return "{$(join(names, ", "))}"
end

"""Build the Mathematica `semishorts -> {{family, {J...}}, ...}` argument."""
function _build_semishorts_arg(ap0020_ops::Vector, atwo0100_ops::Vector)
    entries = String[]
    if !isempty(ap0020_ops)
        spins = join([string(op.J) for op in ap0020_ops], ", ")
        push!(entries, "{\"Ap0020\", {$spins}}")
    end
    if !isempty(atwo0100_ops)
        spins = join([string(op.J) for op in atwo0100_ops], ", ")
        push!(entries, "{\"Atwo0100\", {$spins}}")
    end
    return "{$(join(entries, ", "))}"
end

"""Build the Mathematica `long -> {{J, {tau...}}, ...}` argument."""
function _build_long_arg(long_ops::Vector)
    spins_dict = Dict{Int,Vector{Float64}}()
    for op in long_ops
        push!(get!(spins_dict, op.J, Float64[]), op.tau)
    end
    entries = ["{$spin, {$(join(string.(twists), ", "))}}"
               for (spin, twists) in sort(collect(spins_dict))]
    return "{$(join(entries, ", "))}"
end

# ============================================================================
# Slurm
# ============================================================================

"""
    SlurmConfig

Configuration for Slurm job submission.

# Fields
- `job_name::String`: Job name (default: `"n83dboot_job"`)
- `time::String`: Time limit (default: `"6:00:00"`)
- `nodes::Int`: Number of nodes (default: 1)
- `ntasks::Int`: Number of tasks (default: 1)
- `cpus_per_task::Int`: CPUs per task (default: 1)
- `mem::String`: Memory request (default: `"4G"`)
- `partition::String`: Partition name (default: `"cpu"`)
- `account::String`: Account name (default: `""`)
- `output_dir`: Log output directory; `nothing` uses the global config default
- `array`: Slurm array spec, e.g. `"0-150:2"` (default: `nothing`)
- `extra_options::Dict{String,String}`: Additional `--key=value` SBATCH options
"""
Base.@kwdef mutable struct SlurmConfig
    job_name::String = "n83dboot_job"
    time::String = "6:00:00"
    nodes::Int = 1
    ntasks::Int = 1
    cpus_per_task::Int = 1
    mem::String = "4G"
    partition::String = "cpu"
    account::String = ""
    output_dir::Union{String,Nothing} = nothing
    array::Union{String,Nothing} = nothing
    extra_options::Dict{String,String} = Dict{String,String}()
end

"""
    generate_slurm_script(config::SlurmConfig, command::String;
                          post_command=nothing) -> String

Generate a Slurm batch script string.
`post_command`, if given, is appended after the main command.
"""
function generate_slurm_script(config::SlurmConfig, command::String;
                                post_command::Union{String,Nothing}=nothing)
    output_dir = isnothing(config.output_dir) ? get_config().slurm_output_dir : config.output_dir
    isdir(output_dir) || mkpath(output_dir)

    output_path = isnothing(config.array) ?
        joinpath(output_dir, "$(config.job_name).%j.out") :
        joinpath(output_dir, "$(config.job_name).%A.%a.out")

    script  = "#!/bin/bash\n\n"
    script *= "#SBATCH --job-name=$(config.job_name)\n"
    isnothing(config.array) || (script *= "#SBATCH --array=$(config.array)\n")
    script *= "#SBATCH --time=$(config.time)\n"
    script *= "#SBATCH --nodes=$(config.nodes)\n"
    script *= "#SBATCH --ntasks=$(config.ntasks)\n"
    script *= "#SBATCH --cpus-per-task=$(config.cpus_per_task)\n"
    script *= "#SBATCH --mem=$(config.mem)\n"
    script *= "#SBATCH --partition=$(config.partition)\n"
    isempty(config.account) || (script *= "#SBATCH --account=$(config.account)\n")
    script *= "#SBATCH --output=$(output_path)\n"
    for (key, value) in config.extra_options
        script *= "#SBATCH --$(key)=$(value)\n"
    end
    script *= "\n" * command * "\n"
    if !isnothing(post_command)
        script *= "\n# Post-processing\n" * post_command * "\n"
    end
    return script
end

"""
    submit_slurm_job(config::SlurmConfig, command::String;
                     script_dir=nothing, wait=false, post_command=nothing) -> String

Write a Slurm batch script to `script_dir`, submit it, and return the job ID.
"""
function submit_slurm_job(config::SlurmConfig, command::String;
                          script_dir::Union{String,Nothing}=nothing,
                          wait::Bool=false,
                          post_command::Union{String,Nothing}=nothing)
    output_dir = isnothing(config.output_dir) ? get_config().slurm_output_dir : config.output_dir
    script_dir = isnothing(script_dir) ? joinpath(get_cache_dir(), "slurm_scripts") : script_dir
    isdir(script_dir) || mkpath(script_dir)

    script_content = generate_slurm_script(config, command; post_command=post_command)
    script_path    = joinpath(script_dir, "$(config.job_name)_$(time()).sh")

    open(script_path, "w") do f; write(f, script_content); end
    chmod(script_path, 0o755)

    submit_cmd = wait ? `sbatch --wait $(script_path)` : `sbatch $(script_path)`
    output     = read(submit_cmd, String)

    m = match(r"Submitted batch job (\d+)", output)
    isnothing(m) && error("Failed to extract job ID from sbatch output: $output")
    job_id = m.captures[1]

    @info "Submitted Slurm job $(config.job_name) with ID: $job_id"
    @info "Script saved to: $script_path"
    @info "Output will be written to: $(joinpath(output_dir, config.job_name)).$job_id.out"

    return job_id
end

# ============================================================================
# compute_functional
# ============================================================================

"""
    compute_functional(functional, operators; kwargs...) -> (cache_path, job_id)

Compute functional values at `operators` by invoking `compute_functional.m`.
All operator families (fixed, semishort, long) are bundled into a single
Mathematica call. Submits a Slurm job unless `slurm_config=nothing`.

# Arguments
- `functional::Functional`: The functional to compute
- `operators::Vector{<:Operator}`: Operators to evaluate at
- `math_command::String`: Mathematica executable (default: `"math"`)
- `slurm_config`: `SlurmConfig` to submit via Slurm, or `nothing` for direct run
- `wait::Bool`: Wait for Slurm job completion (default: `false`)

# Returns
`(cache_path::String, job_id::Union{String,Nothing})`
"""
function compute_functional(functional::Functional, operators::Vector{<:Operator};
                            math_command::String="math",
                            slurm_config::Union{SlurmConfig,Nothing}=SlurmConfig(),
                            wait::Bool=false)
    ensure_cache_dir(functional)
    cache_path   = get_cache_path(functional)
    fixed_ops    = filter(op -> op isa FixedOperator, operators)
    ap0020_ops   = filter(op -> op isa SemishortAp0020, operators)
    atwo0100_ops = filter(op -> op isa SemishortAtwo0100, operators)
    long_ops     = filter(op -> op isa LongOperator, operators)

    params      = build_mathematica_params(functional)
    ops_arg     = _build_operators_arg(fixed_ops)
    ss_arg      = _build_semishorts_arg(ap0020_ops, atwo0100_ops)
    long_arg    = _build_long_arg(long_ops)
    script_path = joinpath(_SCRIPTS_DIR, "compute_functional.m")

    cmd = """$(math_command) -script $(script_path) '{$(params), """ *
          """operators -> $(ops_arg), semishorts -> $(ss_arg), long -> $(long_arg), """ *
          """filename -> "$(cache_path)", dataDir -> "$(get_config().data_dir)"}'"""

    return _run_or_submit(cmd, slurm_config, functional, "all"; wait=wait)
end

"""
    _fixed_op_mathematica_name(op::FixedOperator) -> String

Return the Mathematica symbol name for a fixed-OPE operator.
"""
_fixed_op_mathematica_name(::IdentityOperator)  = "\"Id\""
_fixed_op_mathematica_name(::Bp0020Operator)    = "\"Bp0020\""
_fixed_op_mathematica_name(::Bp0040Operator)    = "\"Bp0040\""
_fixed_op_mathematica_name(::Btwo0200Operator)  = "\"Btwo0200\""

"""
    _run_or_submit(cmd, slurm_config, functional, label; wait)

Run a Mathematica command directly or submit it as a Slurm job.
Returns (cache_path, job_id).
"""
function _run_or_submit(cmd::String,
                        slurm_config::Union{SlurmConfig,Nothing},
                        functional::Functional,
                        label::String;
                        wait::Bool=false)
    cache_path = get_cache_path(functional)
    if !isnothing(slurm_config)
        cfg = deepcopy(slurm_config)
        cfg.job_name = "$(cfg.job_name)_$(label)"
        job_id = submit_slurm_job(cfg, cmd; wait=wait)
        return (cache_path, job_id)
    else
        run(`sh -c $cmd`)
        return (cache_path, nothing)
    end
end

# ============================================================================
# Cache existence
# ============================================================================

"""
    cache_exists(functional::Functional) -> Bool

Return `true` if the cache directory for `functional` exists.
"""
function cache_exists(functional::Functional)
    return isdir(get_cache_path(functional))
end

"""
    get_max_spin(operators::Vector{<:Operator}) -> Union{Int,Nothing}

Return the maximum spin among all `LongOperator`s in `operators`,
or `nothing` if there are none.
"""
function get_max_spin(operators::Vector{<:Operator})
    long_ops = filter(op -> op isa LongOperator, operators)
    isempty(long_ops) && return nothing
    return maximum(op.J for op in long_ops)
end

# ============================================================================
# RAM cache
# ============================================================================

"""Global RAM cache: cache_directory_path → Dict(filename → DataFrame)."""
const RAM_CACHE = Dict{String, Dict{String, DataFrame}}()

"""Global operator values cache: (functional, operator) → parsed BigFloat value."""
const OPERATOR_VALUES_CACHE = Dict{Tuple{Functional,Operator}, BigFloat}()

"""
    clear_operator_values_cache!()

Clear the in-memory cache of (functional, operator) → BigFloat values.
Call this if functional CSV data on disk has changed.
"""
function clear_operator_values_cache!()
    empty!(OPERATOR_VALUES_CACHE)
end

"""
    load_from_cache!(functional; override_ram_cache=false)

Load the functional's cached CSV files into `RAM_CACHE` and return the
resulting `Dict{String,DataFrame}`.

Pass `override_ram_cache=true` to force a re-read from disk.
"""
function load_from_cache!(functional::Functional;
                          override_ram_cache::Bool=false,
                          verbose::Bool=false)
    cache_dir_path = get_cache_path(functional)

    if haskey(RAM_CACHE, cache_dir_path) && !override_ram_cache
        return RAM_CACHE[cache_dir_path]
    end

    if !isdir(cache_dir_path)
        @warn "Cache directory does not exist at $cache_dir_path"
        return nothing
    end

    cache_data = _load_cache_from_directory(cache_dir_path; verbose=verbose)
    isempty(cache_data) && return nothing

    RAM_CACHE[cache_dir_path] = cache_data
    return cache_data
end

"""
    _load_cache_from_directory(cache_dir) -> Dict{String,DataFrame}

Read all `.csv` files in `cache_dir`, assigning column names based on file type:
- `spin_J.csv` → columns [:tau, :value]
- `identity.csv`, `Bp0020.csv`, `Bp0040.csv`, `Btwo0200.csv` → column [:value]
- `semishort_*.csv` → column [:value]
"""
function _load_cache_from_directory(cache_dir::String; verbose::Bool=false)
    cache_data = Dict{String,DataFrame}()

    for file in readdir(cache_dir)
        endswith(file, ".csv") || continue
        filepath = joinpath(cache_dir, file)
        verbose && @info "Loading cache file: $filepath"
        try
            df = if startswith(file, "spin_")
                # LongOperator: (tau, value) pairs
                d = CSV.read(filepath, DataFrame, header=false, types=[Float64, String])
                rename!(d, [:tau, :value])
            elseif file in ("identity.csv", "Bp0020.csv", "Bp0040.csv", "Btwo0200.csv") ||
                   startswith(file, "semishort_")
                # Single-value operators: one scalar
                d = CSV.read(filepath, DataFrame, header=false, types=[String])
                rename!(d, [:value])
            else
                CSV.read(filepath, DataFrame, header=false, types=String)
            end
            cache_data[file] = df
        catch e
            if e isa DimensionMismatch
                verbose && @warn "Empty file: $file at $filepath, skipping."
            else
                verbose && @warn "Failed to load $file at $filepath: $e"
            end
        end
    end

    return cache_data
end

# ============================================================================
# get_functional_value
# ============================================================================

"""
    get_functional_value(functional, operator) -> Union{BigFloat,Nothing}

Return the cached value of `functional` at `operator`, or `nothing` if not found.
Loads from CSV into RAM cache on first call.
"""
function get_functional_value(functional::Functional, operator::LongOperator; kwargs...)
    cache_data = load_from_cache!(functional)
    isnothing(cache_data) && return nothing

    spin_file = "spin_$(operator.J).csv"
    haskey(cache_data, spin_file) || return nothing

    df = cache_data[spin_file]
    ("tau" in names(df)) || error("Cache file $spin_file must contain 'tau' column")

    rows = findall(df.tau .≈ operator.tau)
    isempty(rows) && return nothing

    "value" in names(df) || error("Cache file $spin_file must contain 'value' column")
    return parse(BigFloat, string(df[first(rows), :value]))
end

function get_functional_value(functional::Functional, operator::FixedOperator; kwargs...)
    cache_data = load_from_cache!(functional)
    isnothing(cache_data) && return nothing

    fname = operator_cache_filename(operator)
    haskey(cache_data, fname) || return nothing

    df = cache_data[fname]
    nrow(df) == 0 && return nothing
    "value" in names(df) || error("Cache file $fname must contain 'value' column")
    return parse(BigFloat, string(df[1, :value]))
end

function get_functional_value(functional::Functional, operator::SemishortAp0020; kwargs...)
    _get_semishort_value(functional, operator)
end

function get_functional_value(functional::Functional, operator::SemishortAtwo0100; kwargs...)
    _get_semishort_value(functional, operator)
end

function _get_semishort_value(functional::Functional, operator::FreeOperator)
    cache_data = load_from_cache!(functional)
    isnothing(cache_data) && return nothing

    fname = operator_cache_filename(operator)
    haskey(cache_data, fname) || return nothing

    df = cache_data[fname]
    nrow(df) == 0 && return nothing
    "value" in names(df) || error("Cache file $fname must contain 'value' column")
    return parse(BigFloat, string(df[1, :value]))
end

# ============================================================================
# LinearCombinationObjective — no cache, values come directly from the dict
# ============================================================================

# load_from_cache! is a no-op: there is nothing to load from disk.
load_from_cache!(::LinearCombinationObjective; kwargs...) = nothing

# get_functional_value dispatches are defined for each concrete operator type to
# avoid method ambiguity with the existing (Functional, ConcreteOp) methods.
get_functional_value(f::LinearCombinationObjective, ::IdentityOperator;  kwargs...) = BigFloat(0)
get_functional_value(f::LinearCombinationObjective, ::Bp0020Operator;    kwargs...) = BigFloat(0)
get_functional_value(f::LinearCombinationObjective, ::Bp0040Operator;    kwargs...) = BigFloat(0)
get_functional_value(f::LinearCombinationObjective, ::Btwo0200Operator;  kwargs...) = BigFloat(0)

function get_functional_value(f::LinearCombinationObjective, op::SemishortAp0020; kwargs...)
    return BigFloat(get(f.operator_coefficients, op, 0.0))
end
function get_functional_value(f::LinearCombinationObjective, op::SemishortAtwo0100; kwargs...)
    return BigFloat(get(f.operator_coefficients, op, 0.0))
end
function get_functional_value(f::LinearCombinationObjective, op::LongOperator; kwargs...)
    return BigFloat(get(f.operator_coefficients, op, 0.0))
end

# ============================================================================
# RAM cache serialization
# ============================================================================

"""
    dump_ram_cache(path::String)

Save the in-memory RAM cache to a JLD2 file for fast reloading.
"""
function dump_ram_cache(path::String)
    jldsave(path; ram_cache=RAM_CACHE)
    @info "Saved $(length(RAM_CACHE)) cache entries to $path"
end

"""
    load_ram_cache!(path::String; replace::Bool=false)

Populate the in-memory RAM cache from a JLD2 snapshot created by `dump_ram_cache`.
Skips keys already present unless `replace=true`.
"""
function load_ram_cache!(path::String; replace::Bool=false)
    loaded = jldopen(path, "r") do f; f["ram_cache"]; end
    n_new  = 0
    for (k, v) in loaded
        if replace || !haskey(RAM_CACHE, k)
            RAM_CACHE[k] = v
            n_new += 1
        end
    end
    @info "Loaded $n_new new cache entries from $path ($(length(loaded)) total in file)"
    return RAM_CACHE
end

# ============================================================================
# Cache completeness checks
# ============================================================================

"""
    check_cache_complete(functional, operators) -> Bool

Return `true` if every operator in `operators` has a cached value.
"""
function check_cache_complete(functional::Functional, operators::Vector{<:Operator})
    cache_exists(functional) || return false
    return all(!isnothing(get_functional_value(functional, op)) for op in operators)
end

"""
    find_missing_functionals(functionals, operators; verbose=true) -> Vector{Functional}

Return the subset of `functionals` that have incomplete or missing caches for `operators`.
"""
function find_missing_functionals(functionals::Vector{<:Functional},
                                  operators::Vector{<:Operator};
                                  verbose::Bool=true)
    missing_functionals = Functional[]

    if verbose
        println("\n" * "="^70)
        println("Checking cache completeness for $(length(functionals)) functionals")
        println("="^70)
    end

    for (i, functional) in enumerate(functionals)
        verbose && print("[$i/$(length(functionals))] Checking $(functional)... ")
        if check_cache_complete(functional, operators)
            verbose && println("complete")
        else
            verbose && println("missing/incomplete")
            push!(missing_functionals, functional)
        end
    end

    if verbose
        println("="^70)
        isempty(missing_functionals) ?
            println("All functionals cached!") :
            println("$(length(missing_functionals))/$(length(functionals)) functionals need computation")
        println("="^70)
        println()
    end

    return missing_functionals
end

"""
    find_missing_operators(functional, operators; verbose=true) -> Vector{Operator}

Return the operators in `operators` that are missing from the cache for `functional`.
"""
function find_missing_operators(functional::Functional,
                                operators::Vector{<:Operator};
                                verbose::Bool=true)
    missing_operators = Operator[]
    verbose && println("\nChecking operators for $(functional):")

    for op in operators
        if isnothing(get_functional_value(functional, op))
            push!(missing_operators, op)
            verbose && println("  Missing: $op")
        else
            verbose && println("  Cached: $op")
        end
    end

    if verbose
        isempty(missing_operators) ?
            println("  All operators cached!") :
            println("  $(length(missing_operators))/$(length(operators)) operators missing")
    end

    return missing_operators
end

# ============================================================================
# Cache bundle packaging and distribution
# ============================================================================

"""
    package_cache(functionals; output_path="cache_bundle.tar.gz") -> String

Create a `.tar.gz` archive containing the cache directories for all given
functionals and save it to `output_path`.

Returns the absolute path to the created archive.
"""
function package_cache(functionals::Vector{<:Functional};
                       output_path::String="cache_bundle.tar.gz")
    cache_dir  = get_cache_dir()
    output_abs = isabspath(output_path) ? output_path : joinpath(pwd(), output_path)
    tar_path   = endswith(output_abs, ".gz") ? output_abs[1:end-3] : output_abs * ".tar"

    isfile(tar_path) && rm(tar_path)

    n_included = 0
    cd(cache_dir) do
        for functional in functionals
            func_dir = get_cache_path(functional)
            if !isdir(func_dir)
                @warn "Cache directory missing, skipping: $func_dir"
                continue
            end
            rel_path = relpath(func_dir, cache_dir)
            flag = n_included == 0 ? "cf" : "rf"
            run(pipeline(`tar $flag $tar_path $rel_path`, stdout=devnull, stderr=devnull))
            n_included += 1
        end
    end

    n_included == 0 && error("No valid cache directories found to package.")

    run(`gzip -f $tar_path`)
    checksum = bytes2hex(open(sha256, output_abs))
    @info "Cache bundle written" output=output_abs sha256=checksum n_functionals=n_included

    return output_abs
end

"""
    download_cache(url_or_filename; cache_dir=get_cache_dir(), expected_sha256=nothing)

Download a `.tar.gz` cache bundle and extract it into `cache_dir`.
The RAM cache is cleared after extraction.
"""
function download_cache(url_or_filename::String;
                        cache_dir::String=get_cache_dir(),
                        expected_sha256::Union{String,Nothing}=nothing)
    _GITHUB_RELEASES_BASE = "https://github.com/srossd/N83DBoot.jl/releases/latest/download/"
    url = startswith(url_or_filename, r"https?://") ? url_or_filename :
          _GITHUB_RELEASES_BASE * url_or_filename

    isdir(cache_dir) || mkpath(cache_dir)

    tmp = tempname() * ".tar.gz"
    try
        @info "Downloading cache bundle..." url
        Downloads.download(url, tmp)

        if !isnothing(expected_sha256)
            actual = bytes2hex(open(sha256, tmp))
            if actual != lowercase(expected_sha256)
                error("SHA256 mismatch: expected $(expected_sha256), got $actual")
            end
            @info "Checksum verified" sha256=actual
        end

        run(`tar xzf $tmp -C $cache_dir`)

        # Invalidate overlapping RAM cache entries
        cache_dir_norm = rstrip(cache_dir, '/')
        for key in collect(keys(RAM_CACHE))
            if startswith(rstrip(key, '/'), cache_dir_norm)
                delete!(RAM_CACHE, key)
            end
        end
        empty!(OPERATOR_VALUES_CACHE)

        @info "Cache bundle extracted" destination=cache_dir
    finally
        isfile(tmp) && rm(tmp)
    end
end

"""
    list_cached_functionals(subdir::String) -> Vector{Functional}

Return all functionals whose cache directories exist under `cache_dir/subdir`.
"""
function list_cached_functionals(subdir::String)
    full_dir = joinpath(get_cache_dir(), subdir)
    isdir(full_dir) || return Functional[]

    functionals = Functional[]
    for entry in readdir(full_dir)
        path = joinpath(full_dir, entry)
        isdir(path) || continue
        try
            push!(functionals, _functional_from_cache_path(subdir, entry))
        catch e
            @warn "Failed to parse functional from path $path: $e"
        end
    end
    return functionals
end

"""
    _functional_from_cache_path(type_name, dirname) -> Functional

Reconstruct a `Functional` from its cache directory name.
"""
function _functional_from_cache_path(type_name::String, dirname::String)
    parts = split(dirname, "_")
    if type_name == "Crossing"
        # cross_1_m_4_n_0_precision_1024_rOrder_100
        cross     = parse(Int, parts[2])
        m         = parse(Int, parts[4])
        n         = parse(Int, parts[6])
        precision = parse(Int, parts[8])
        r_order   = parse(Int, parts[10])
        return CrossingFunctional(m, n, cross; precision=precision, r_order=r_order)
    elseif type_name == "Integral"
        # b_0.0_precision_256_rOrder_20
        b         = parse(Float64, parts[2])
        precision = parse(Int, parts[4])
        r_order   = parse(Int, parts[6])
        return IntegralFunctional(b; precision=precision, r_order=r_order)
    else
        error("Unknown functional type: $type_name")
    end
end

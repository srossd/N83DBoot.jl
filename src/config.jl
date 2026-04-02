"""
Configuration system for N83DBoot (3D N=8 SCFT bootstrap).

Configuration priority (highest to lowest):
1. Environment variables
2. Runtime values set via `set_config!()`
3. Saved preferences in `LocalPreferences.toml` (written by `save_config!()`)
4. Hardcoded fallback defaults

Use `save_config!()` to persist your settings across Julia sessions.

# Environment Variables
- `N83DBOOT_CACHE_DIR`: Cache directory for functional values
- `N83DBOOT_SLURM_OUT`: Slurm job output directory
- `N83DBOOT_SDPB_BUILD`: Path to SDPB build directory, or `"docker"` to run via Docker
- `N83DBOOT_SDPB_DOCKER_IMAGE`: Docker image to use when `sdpb_build_dir = "docker"`
- `N83DBOOT_SDPB_WORK`: SDPB working directory
- `N83DBOOT_SDPB_RESULTS`: SDPB results directory
- `N83DBOOT_LOCALIZATION_DIR`: Directory containing localization CSV files
- `N83DBOOT_DATA_DIR`: Directory containing ABJM input data files for Mathematica scripts
"""

using Preferences

"""
    N83DBootConfig

Configuration for N83DBoot.

# Fields
- `cache_dir::String`: Directory for cached functional values
- `slurm_output_dir::String`: Directory for Slurm job output logs
- `sdpb_build_dir::String`: Path to SDPB build directory containing `sdpb` and `sdp2input`
  executables, **or** the special value `"docker"` to run via Docker (see `sdpb_docker_image`).
- `sdpb_docker_image::String`: Docker image used when `sdpb_build_dir = "docker"`.
  Defaults to `"bootstrapcollaboration/sdpb:master"`. The working directory is
  automatically mounted into the container; all paths are remapped accordingly.
- `sdpb_work_dir::String`: Working directory for SDPB computations
- `sdpb_results_dir::String`: Directory for storing SDPB results
- `localization_data_dir::String`: Directory containing localization CSV files
  (short_coeffs_k1.csv, short_coeffs_k2.csv, short_coeffs_k2_M1.csv, rhs_k*.csv)
- `data_dir::String`: Directory containing ABJM input data files used by Mathematica scripts
  (superblocks.m, scalar_blocks.m, block-expansion-Lmax70-20.m, blocks/, gint_cache/).
  Defaults to the package scripts directory if empty.
"""
Base.@kwdef mutable struct N83DBootConfig
    cache_dir::String             = ""
    slurm_output_dir::String      = ""
    sdpb_build_dir::String        = ""
    sdpb_docker_image::String     = "bootstrapcollaboration/sdpb:master"
    sdpb_work_dir::String         = ""
    sdpb_results_dir::String      = ""
    localization_data_dir::String = ""
    data_dir::String              = ""
end

"""
    _load_preferences() -> N83DBootConfig

Build a `N83DBootConfig` from values saved in `LocalPreferences.toml`,
falling back to hardcoded defaults for any key that has not been saved yet.
"""
function _load_preferences()
    return N83DBootConfig(
        cache_dir             = @load_preference("cache_dir",             ""),
        slurm_output_dir      = @load_preference("slurm_output_dir",      ""),
        sdpb_build_dir        = @load_preference("sdpb_build_dir",        ""),
        sdpb_docker_image     = @load_preference("sdpb_docker_image",     "bootstrapcollaboration/sdpb:master"),
        sdpb_work_dir         = @load_preference("sdpb_work_dir",         ""),
        sdpb_results_dir      = @load_preference("sdpb_results_dir",      ""),
        localization_data_dir = @load_preference("localization_data_dir", ""),
        data_dir              = @load_preference("data_dir",              ""),
    )
end

"""Global configuration instance, initialized from saved preferences."""
const CONFIG = Ref{N83DBootConfig}(_load_preferences())

"""
    get_config() -> N83DBootConfig

Return the current configuration with environment variable overrides applied.
"""
function get_config()
    config = CONFIG[]

    return N83DBootConfig(
        cache_dir             = get(ENV, "N83DBOOT_CACHE_DIR",          config.cache_dir),
        slurm_output_dir      = get(ENV, "N83DBOOT_SLURM_OUT",          config.slurm_output_dir),
        sdpb_build_dir        = get(ENV, "N83DBOOT_SDPB_BUILD",          config.sdpb_build_dir),
        sdpb_docker_image     = get(ENV, "N83DBOOT_SDPB_DOCKER_IMAGE",   config.sdpb_docker_image),
        sdpb_work_dir         = get(ENV, "N83DBOOT_SDPB_WORK",           config.sdpb_work_dir),
        sdpb_results_dir      = get(ENV, "N83DBOOT_SDPB_RESULTS",        config.sdpb_results_dir),
        localization_data_dir = get(ENV, "N83DBOOT_LOCALIZATION_DIR",    config.localization_data_dir),
        data_dir              = get(ENV, "N83DBOOT_DATA_DIR",            config.data_dir),
    )
end

"""
    set_config!(; kwargs...)

Set configuration values for the current Julia session. Keyword arguments
match `N83DBootConfig` field names. Changes are lost when Julia exits —
use `save_config!()` to persist them.

# Examples
```julia
set_config!(cache_dir="/my/cache", sdpb_build_dir="/my/sdpb/build")
```
"""
function set_config!(; kwargs...)
    config = CONFIG[]

    env_map = [
        (:cache_dir,             "N83DBOOT_CACHE_DIR"),
        (:slurm_output_dir,      "N83DBOOT_SLURM_OUT"),
        (:sdpb_build_dir,        "N83DBOOT_SDPB_BUILD"),
        (:sdpb_docker_image,     "N83DBOOT_SDPB_DOCKER_IMAGE"),
        (:sdpb_work_dir,         "N83DBOOT_SDPB_WORK"),
        (:sdpb_results_dir,      "N83DBOOT_SDPB_RESULTS"),
        (:localization_data_dir, "N83DBOOT_LOCALIZATION_DIR"),
        (:data_dir,              "N83DBOOT_DATA_DIR"),
    ]

    env_overrides = [
        "$key (overridden by $env_var=$(ENV[env_var]))"
        for (key, env_var) in env_map
        if haskey(kwargs, key) && haskey(ENV, env_var)
    ]

    if !isempty(env_overrides)
        @warn "The following settings will be overridden by environment variables:\n  " *
              join(env_overrides, "\n  ")
    end

    for (key, value) in kwargs
        if hasfield(N83DBootConfig, key)
            setfield!(config, key, value)
        else
            @warn "Unknown configuration key: $key"
        end
    end

    CONFIG[] = config
end

"""
    save_config!(; kwargs...)

Persist the current configuration to `LocalPreferences.toml` in the active
Julia project, making the settings the default for future Julia sessions.
Accepts the same keyword arguments as `set_config!()` — if any are provided,
they are applied before saving.

Call `reset_config!()` (or restart Julia) to reload the saved values as the
session defaults.

# Examples
```julia
save_config!(
    cache_dir              = "/scratch/myuser/cache",
    localization_data_dir  = "/scratch/myuser/localization",
    sdpb_build_dir         = "/home/myuser/sdpb/build",
)
```
"""
function save_config!(; kwargs...)
    isempty(kwargs) || set_config!(; kwargs...)
    config = get_config()
    @set_preferences!(
        "cache_dir"             => config.cache_dir,
        "slurm_output_dir"      => config.slurm_output_dir,
        "sdpb_build_dir"        => config.sdpb_build_dir,
        "sdpb_docker_image"     => config.sdpb_docker_image,
        "sdpb_work_dir"         => config.sdpb_work_dir,
        "sdpb_results_dir"      => config.sdpb_results_dir,
        "localization_data_dir" => config.localization_data_dir,
        "data_dir"              => config.data_dir,
    )
    @info "Configuration saved to LocalPreferences.toml. Call reset_config!() or restart Julia to apply as defaults."
end

"""
    reset_config!()

Reset the session configuration to the saved preferences (or hardcoded
defaults if no preferences have been saved yet).
"""
function reset_config!()
    CONFIG[] = _load_preferences()
end

"""
    set_cache_dir!(path::String)

Set the cache directory path (convenience wrapper for `set_config!(cache_dir=path)`).
"""
function set_cache_dir!(path::String)
    set_config!(cache_dir=path)
end

"""
    get_cache_dir() -> String

Return the active cache directory path.
"""
function get_cache_dir()
    return get_config().cache_dir
end

const _FIELD_LABELS = Dict{Symbol,String}(
    :cache_dir             => "Cache directory (for computed functional values)",
    :sdpb_build_dir        => "SDPB build directory (path to dir containing sdpb/sdp2input, or \"docker\" to run via Docker)",
    :sdpb_docker_image     => "Docker image (used when sdpb_build_dir=\"docker\"; default: bootstrapcollaboration/sdpb:master)",
    :sdpb_work_dir         => "SDPB working directory (for intermediate files)",
    :sdpb_results_dir      => "SDPB results directory",
    :localization_data_dir => "Localization data directory (short_coeffs_*.csv, rhs_*.csv files)",
    :data_dir              => "ABJM data directory (block expansion files, blocks/, etc.)",
    :slurm_output_dir      => "Slurm output directory",
)

"""
    configure_interactively(fields=nothing)

Prompt the user to set any unset configuration values, then save them to
`LocalPreferences.toml`. If `fields` is provided (a vector of `Symbol`s
matching `N83DBootConfig` field names), only those fields are considered;
otherwise all fields are checked.

Only fields whose current value is empty are prompted — already-configured
fields are left untouched.

# Examples
```julia
# Prompt only for the fields needed to run an SDP locally
configure_interactively([:cache_dir, :sdpb_build_dir, :sdpb_work_dir, :sdpb_results_dir])
```
"""
function configure_interactively(fields=nothing)
    all_fields = [
        :cache_dir, :slurm_output_dir, :sdpb_build_dir,
        :sdpb_work_dir, :sdpb_results_dir, :localization_data_dir, :data_dir, :script_dir,
    ]
    fields = isnothing(fields) ? all_fields : fields

    config = get_config()
    unconfigured = filter(fields) do f
        val = getfield(config, f)
        val isa String ? isempty(val) : false
    end

    isempty(unconfigured) && return

    println("The following configuration values are not set.")
    println("Enter a path for each, then press Enter. Values will be saved to LocalPreferences.toml.\n")

    kwargs = Dict{Symbol,Any}()
    for field in unconfigured
        label = get(_FIELD_LABELS, field, string(field))
        print("  $label: ")
        input = strip(readline())
        if !isempty(input)
            kwargs[field] = String(input)
        end
    end

    if !isempty(kwargs)
        save_config!(; kwargs...)
        println()
    end
end

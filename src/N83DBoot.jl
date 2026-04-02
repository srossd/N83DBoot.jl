module N83DBoot

"""
    N83DBoot

Julia package for conformal bootstrap calculations in 3D N=8 SCFTs (ABJM theory).

Provides type hierarchies for operators and functionals, a CSV caching system
for computed functional values, an interface to external Mathematica computation
scripts, and SDP construction / SDPB solver integration.

The operator spectrum contains:
- Fixed-OPE operators: Identity, Bp0020, Bp0040, Btwo0200
  (OPE coefficients determined by localization via nosaka.m)
- Semishort families: SemishortAp0020(J) [even J], SemishortAtwo0100(J) [odd J]
  (free OPE coefficients, dimension pinned to J+1)
- Long multiplets: LongOperator(tau, J) (free OPE coefficients, τ = Δ−J)

# Examples

```julia
using N83DBoot

# Configure directories
set_cache_dir!("/path/to/cache")
set_config!(localization_data_dir="/path/to/localization_data")

# Load localization data for k=1 ABJM at rank N=4
loc = load_localization_data(4, 1)

# Build operators and constraints
operators = operator_grid(long_spins=collect(0:20), max_twist=40.0)
constraints = vcat(
    crossing_constraints(Lambda=9),
    integral_constraints(loc)
)

# Solve
sdp = SDPProblem(CrossingFunctional(4, 0, 1), constraints, operators, loc)
result = solve_sdp(sdp)
```
"""

# Re-export LaTeXStrings so users can write L"..." labels without a separate import
using LaTeXStrings: LaTeXString, @L_str
export LaTeXString, @L_str

# Include source files
include("config.jl")
include("operators.jl")
include("functionals.jl")
include("localization.jl")
include("cache.jl")
include("plot.jl")
include("sdp.jl")
include("util.jl")

# ── Operator types ──────────────────────────────────────────────────────────
export Operator, FixedOperator, FreeOperator
export IdentityOperator, Bp0020Operator, Bp0040Operator, Btwo0200Operator
export Identity, Bp0020, Bp0040, Btwo0200
export SemishortAp0020, SemishortAtwo0100
export LongOperator

# ── Functional types ─────────────────────────────────────────────────────────
export Functional, CrossingFunctional, IntegralFunctional, LinearCombinationObjective

# ── Cache functions ───────────────────────────────────────────────────────────
export set_cache_dir!, get_cache_dir, compute_functional
export cache_exists, load_from_cache!, get_functional_value
export list_cached_functionals, get_max_spin
export check_cache_complete, find_missing_functionals, find_missing_operators
export dump_ram_cache, load_ram_cache!
export package_cache, download_cache
export clear_operator_values_cache!

# ── Slurm ────────────────────────────────────────────────────────────────────
export SlurmConfig, generate_slurm_script, submit_slurm_job

# ── SDP types ─────────────────────────────────────────────────────────────────
export ComparisonOperator, LessThanOrEqual, Equal, GreaterThanOrEqual
export SDPConstraint, SDPProblem, DualizedSDP, SDPBResult

# ── SDP functions ─────────────────────────────────────────────────────────────
export dualize_sdp, write_sdp_json, write_sdp_mathematica
export parse_sdpb_output
export build_sdp2input_command, build_sdpb_command
export solve_sdp
export extract_NkM, extract_Lambda

# ── Plot ──────────────────────────────────────────────────────────────────────
export Metadata, Objective, Plot
export JointBoundsPlot
export get_data, generate_plot

# ── Util ──────────────────────────────────────────────────────────────────────
export operator_grid, twist_spacing, crossing_constraints, integral_constraints
export default_twist_spacing

# ── Configuration ─────────────────────────────────────────────────────────────
export DN8BootConfig, get_config, set_config!, save_config!, reset_config!, configure_interactively

# ── Localization ──────────────────────────────────────────────────────────────
export LocalizationData, load_localization_data, get_localization_dir
export package_localization_data, download_localization_data

end

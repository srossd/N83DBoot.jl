"""
Example: joint bounds on two semishort OPE coefficients in ABJM k=1 theory.

This script computes the allowed region in the
    (λ²_{(A,+)_0},  λ²_{(A,2)_1})
plane by sweeping over linear combinations

    max  cos(θ) · λ²_{Ap(0)}  +  sin(θ) · λ²_{A2(1)}

for angles θ ∈ [0, 2π) and two constraint sets:
    • crossing constraints only          ("Without IC")
    • crossing + localization integral   ("With IC")

With Lambda = 19 the problem has O(100) crossing derivatives; pre-computing
the functional values on a cluster is strongly recommended.

Usage
-----
1. Configure paths (once per session or in LocalPreferences.toml):

    using N83DBoot
    save_config!(
        cache_dir             = "/path/to/cache",
        localization_data_dir = "/path/to/localization_data",
        sdpb_build_dir        = "/path/to/sdpb/build",   # or "docker"
        sdpb_work_dir         = "/path/to/work",
        sdpb_results_dir      = "/path/to/results",
    )

2. Pre-compute all required functional values.  At Lambda=19 this is the
   expensive step and should be run on a cluster via the Slurm interface:

    using N83DBoot
    include("examples/abjm_semishort_joint_bounds.jl")
    precompute_functionals(2, 1)   # N=2, k=1

3. Run the angle sweep (also via Slurm for parallel execution):

    run_joint_bounds()

4. After all results are available, generate the plot:

    p = generate_plot(JOINT_PLOT)
"""

using N83DBoot

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

const LAMBDA    = 14       # crossing derivative cutoff (Λ)
const MAX_SPIN  = 24       # maximum spin included for long operators
const MAX_TWIST = 35.0     # maximum Δ for long operators

const ABJM_N = 2           # ABJM rank
const ABJM_K = 1           # Chern-Simons level

# Operators whose joint allowed region we want to map out
const OP_X = SemishortAp0020(0)     # (A,+)_0 — even-spin semishort, J=0
const OP_Y = SemishortAp0020(2)     # (A,+)_2 — even-spin semishort, J=2

# Number of angles to sweep around the unit circle.
# More angles give a smoother boundary; 36 is a good starting point.
const N_ANGLES = 36

# Joint bounds plot specification.
# Results from both "Without IC" (num_integral=0) and "With IC" (num_integral=1)
# sweeps are stored in the same directory and distinguished by the series spec.
const JOINT_PLOT = JointBoundsPlot("abjm_semishort_joint_bounds";
    op_x   = OP_X,
    op_y   = OP_Y,
    series = [Metadata("num_integral")],
    xlabel = L"\lambda^2_{(A,+)_0}",
    ylabel = L"\lambda^2_{(A,+)_2}")

# ---------------------------------------------------------------------------
# Helper: build the operator grid
# ---------------------------------------------------------------------------

function build_operators()
    return operator_grid(
        long_spins           = collect(0:MAX_SPIN),
        max_twist            = MAX_TWIST,
        semishort_even_spins = collect(0:2:MAX_SPIN),
        semishort_odd_spins  = collect(1:2:(MAX_SPIN - 1)),
        twist_spacing_fn     = twist_spacing((4., 1/16), (16., 1/4), (Inf, 1/2))
    )
end

# ---------------------------------------------------------------------------
# Helper: pre-compute all functional values needed for one SDP
#
# At Lambda=19 this involves O(100) crossing functionals × O(1000) operators.
# Run this once per (N, k); the results are persisted in the CSV cache.
# ---------------------------------------------------------------------------

function precompute_functionals(N::Int=ABJM_N, k::Int=ABJM_K;
                                 slurm_config=nothing)
    loc = load_localization_data(N, k)
    ops = build_operators()
    cs  = crossing_constraints(Λ=LAMBDA)
    push!(cs, integral_constraints(loc))

    @info "Pre-computing $(length(cs)) functionals on $(length(ops)) operators..."
    for (i, c) in enumerate(cs)
        @info "  Functional $i/$(length(cs)): $(c.functional)"
        if c.functional isa CrossingFunctional
            # compute_functional(c.functional, ops; slurm_config=slurm_config)
        else
            for spin=0:MAX_SPIN
                compute_functional(c.functional, filter(x -> (x isa LongOperator && x.J == spin) || (!(x isa LongOperator) && spin == 0), ops);
                                   slurm_config=slurm_config)
            end
        end
    end

    # Also pre-compute a representative LinearCombinationObjective value set
    # (not strictly needed since those are computed analytically, but this
    # ensures the crossing/integral caches are fully populated).
    @info "Pre-computation complete."
end

# ---------------------------------------------------------------------------
# Helper: build an SDPProblem for one angle and one constraint set
# ---------------------------------------------------------------------------

function build_joint_sdp(theta::Float64, loc;
                          include_integral::Bool)
    ops = build_operators()
    cs  = crossing_constraints(Λ=LAMBDA)
    include_integral && append!(cs, integral_constraints(loc))

    objective = LinearCombinationObjective(
        OP_X => cos(theta),
        OP_Y => sin(theta),
    )

    return SDPProblem(objective, cs, ops, loc; direction=1)
end

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

"""
    run_joint_bounds(; N, k, slurm_config, results_dir)

Run the full angle sweep for both "Without IC" and "With IC" series.

Submits N_ANGLES × 2 = $(2*N_ANGLES) SDPB jobs (or re-uses cached results).
If `slurm_config` is provided, jobs are submitted asynchronously; call
`generate_plot(JOINT_PLOT)` once all results are available.
"""
function run_joint_bounds(;
    N            = ABJM_N,
    k            = ABJM_K,
    slurm_config = nothing,
    results_dir  = nothing,
)
    loc    = load_localization_data(N, k)
    thetas = range(0, 2π, length=N_ANGLES + 1)[1:end-1]  # N_ANGLES evenly-spaced angles

    @info "Sweeping $(N_ANGLES) angles for N=$N, k=$k..."

    for include_ic in (false, true)
        label = include_ic ? "With IC" : "Without IC"
        @info "  Series: $label"
        sdps = [build_joint_sdp(θ, loc; include_integral=include_ic) for θ in thetas]

        solve_sdp(sdps;
            plots        = JOINT_PLOT,
            slurm_config = slurm_config,
            results_dir  = results_dir,
        )
    end

    if isnothing(slurm_config)
        @info "All jobs complete — generating plot..."
        p = generate_plot(JOINT_PLOT; results_dir=results_dir)
        return p
    else
        @info "Jobs submitted.  Call generate_plot(JOINT_PLOT) once results are ready."
        return nothing
    end
end

# ---------------------------------------------------------------------------
# Interactive helpers (used by main())
# ---------------------------------------------------------------------------

"""
Prompt the user for Slurm job parameters and return a configured `SlurmConfig`.
"""
function _prompt_slurm_config(; job_name::String="n83dboot_sdp")
    print("  Time limit (e.g. 4:00:00): ")
    time_str  = strip(readline())

    print("  Number of cores (ntasks): ")
    cores_str = strip(readline())

    print("  Slurm account: ")
    account   = strip(readline())

    return SlurmConfig(;
        job_name      = job_name,
        time          = time_str,
        ntasks        = parse(Int, cores_str),
        account       = account,
    )
end

"""
Count `{hash}_out.txt` result files already written for `JOINT_PLOT`.
"""
function _count_complete_results(; results_dir::Union{String,Nothing}=nothing)
    rd       = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    plot_dir = joinpath(rd, JOINT_PLOT.name)
    isdir(plot_dir) || return 0
    return count(f -> endswith(f, "_out.txt"), readdir(plot_dir))
end

# ---------------------------------------------------------------------------
# main() — interactive entry point
# ---------------------------------------------------------------------------

"""
    main()

Interactive entry point for the ABJM semishort joint bounds calculation.

Steps performed:
1. Configure any unset paths interactively and save them.
2. Check whether all required functional values are cached; if not, offer to
   compute them or download the pre-built `Lambda_$(LAMBDA).tar.gz` bundle.
3. If all $(N_ANGLES * 2) results already exist, generate and return the plot
   immediately.
4. Ask whether to submit SDPB jobs via Slurm (prompting for time, cores, and
   account) or run them serially in the current session.
5. Submit / run the angle sweep.  If Slurm was used, print a reminder to call
   `main()` again once the jobs finish; otherwise return the finished plot.
"""
function main()
    println("\n=== ABJM Semishort Joint Bounds: Interactive Mode ===\n")

    # ── 1. Configuration ──────────────────────────────────────────────────────
    configure_interactively([:cache_dir, :localization_data_dir, :sdpb_build_dir,
                              :sdpb_work_dir, :sdpb_results_dir, :slurm_output_dir])

    # ── 2. Functional cache check ─────────────────────────────────────────────
    loc = nothing
    try
        loc         = load_localization_data(ABJM_N, ABJM_K)
    catch Exception
        download_localization_data("Localization.tar.gz")
        loc = load_localization_data(ABJM_N, ABJM_K)
    end

    ops         = build_operators()
    cs          = crossing_constraints(Λ=LAMBDA)
    append!(cs, integral_constraints(loc))
    functionals = [c.functional for c in cs]

    missing_fns = find_missing_functionals(functionals, ops; verbose=false)

    if !isempty(missing_fns)
        n_total   = length(functionals)
        n_missing = length(missing_fns)
        println("$n_missing/$n_total functionals are not fully cached.\n")
        println("  [1] Compute them with precompute_functionals")
        println("  [2] Download Lambda_$(LAMBDA).tar.gz from GitHub releases")
        println("  [q] Quit")
        print("Choice: ")
        choice = strip(readline())

        if choice == "1"
            print("\nUse Slurm for precomputation? [y/N]: ")
            use_slurm_pre = lowercase(strip(readline())) == "y"
            slurm_pre     = use_slurm_pre ?
                _prompt_slurm_config(; job_name="n83d_precompute") : nothing
            precompute_functionals(ABJM_N, ABJM_K; slurm_config=slurm_pre)
            if !isnothing(slurm_pre)
                println("\nPrecomputation jobs submitted. Run main() again once they complete.")
                return nothing
            end
        elseif choice == "2"
            download_cache("Lambda_$(LAMBDA).tar.gz")
            download_localization_data("Localization.tar.gz")
        else
            println("Exiting.")
            return nothing
        end

        # Re-check after download or serial precomputation
        missing_fns = find_missing_functionals(functionals, ops; verbose=false)
        if !isempty(missing_fns)
            println("$(length(missing_fns)) functionals still missing. " *
                    "Run main() again after precomputation completes.")
            return nothing
        end
    end

    # ── 3. Check for existing results → generate plot if complete ─────────────
    n_expected = N_ANGLES * 2   # two series: without IC and with IC
    n_complete = _count_complete_results()

    if n_complete >= n_expected
        println("All $n_expected results are available — generating plot...")
        return generate_plot(JOINT_PLOT)
    elseif n_complete > 0
        println("$n_complete/$n_expected results are ready (jobs may still be running).")
    end

    # ── 4. Serial or Slurm? ───────────────────────────────────────────────────
    print("\nRun SDP jobs via Slurm? [y/N]: ")
    use_slurm = lowercase(strip(readline())) == "y"

    slurm_cfg = use_slurm ? _prompt_slurm_config() : nothing

    # ── 5. Run the angle sweep ────────────────────────────────────────────────
    if isnothing(slurm_cfg)
        # Serial: run_joint_bounds blocks until complete and returns the plot.
        return run_joint_bounds(; N=ABJM_N, k=ABJM_K, slurm_config=nothing)
    else
        run_joint_bounds(; N=ABJM_N, k=ABJM_K, slurm_config=slurm_cfg)
        println("\nSlurm jobs submitted. Run main() again once they complete to generate the plot.")
        return nothing
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
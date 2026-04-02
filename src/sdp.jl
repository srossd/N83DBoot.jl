# ============================================================================
# SDP (Semidefinite Program) Construction and SDPB Solver Integration
# ============================================================================

using JSON
using Logging
using SHA
using Dates
using Printf

# ============================================================================
# Type Definitions
# ============================================================================

"""
    ComparisonOperator

Enumeration for constraint comparison operators.

# Values
- `LessThanOrEqual`: ≤ constraint
- `Equal`: = constraint
- `GreaterThanOrEqual`: ≥ constraint
"""
@enum ComparisonOperator begin
    LessThanOrEqual
    Equal
    GreaterThanOrEqual
end

"""
    SDPConstraint

A single constraint in an SDP problem.

Represents a constraint of the form: functional(operators) ⊙ rhs
where ⊙ is one of ≤, =, or ≥.

# Fields
- `functional::Functional`: The functional that defines this constraint
- `rhs::BigFloat`: Right-hand side value
- `comparison::ComparisonOperator`: Comparison operator

# Examples
```julia
f = CrossingFunctional(4, 0, 1)
c1 = SDPConstraint(f)                              # f = 0 (default)
c2 = SDPConstraint(f, rhs=1.5, comparison=GreaterThanOrEqual)
```
"""
struct SDPConstraint
    functional::Functional
    rhs::BigFloat
    comparison::ComparisonOperator

    function SDPConstraint(functional::Functional;
                          rhs::Real=0,
                          comparison::ComparisonOperator=Equal)
        @assert functional isa Functional "functional must be a Functional type"
        new(functional, BigFloat(rhs), comparison)
    end
end

function Base.show(io::IO, c::SDPConstraint)
    comp_symbol = c.comparison == Equal ? "=" :
                  c.comparison == LessThanOrEqual ? "≤" : "≥"
    print(io, "SDPConstraint($(c.functional) $(comp_symbol) $(c.rhs))")
end

"""
    SDPProblem

A semidefinite program in user-friendly form (before dualization).

# Fields
- `objective::Functional`: The objective functional to maximize
- `constraints::Vector{SDPConstraint}`: List of constraint functionals
- `operators::Vector{Operator}`: Operators at which to evaluate functionals
- `localization::LocalizationData`: Localization data (OPE² for fixed operators, integral RHS)
- `direction::Int`: Optimization direction (+1 = maximize, -1 = minimize)
- `precision::Int`: Precision in bits for numerical computations

The `operators` vector must contain exactly one `IdentityOperator`, which anchors
the normalization.  The other `FixedOperator` instances (`Bp0020`, `Bp0040`,
`Btwo0200`) are looked up from `localization`; they do **not** need to be
listed separately in `operators` — they are always included internally.

# Examples
```julia
loc = load_localization_data(4, 1)
objective = CrossingFunctional(4, 0, 1)
constraints = [SDPConstraint(IntegralFunctional(0.0), comparison=GreaterThanOrEqual)]
operators = [Identity, Bp0020, Bp0040, Btwo0200,
             SemishortAp0020(0), LongOperator(2.0, 0)]
sdp = SDPProblem(objective, constraints, operators, loc)
```
"""
struct SDPProblem
    objective::Functional
    constraints::Vector{SDPConstraint}
    operators::Vector{Operator}
    localization::LocalizationData
    direction::Int
    precision::Int

    function SDPProblem(objective::Functional,
                       constraints::Vector{SDPConstraint},
                       operators::Vector{Operator},
                       localization::LocalizationData;
                       precision::Int=1024, direction::Int=1)
        @assert !isempty(constraints) "Must have at least one constraint"
        @assert !isempty(operators) "Must have at least one operator"
        @assert precision > 0 "Precision must be positive"
        @assert direction ∈ (1, -1) "direction must be +1 or -1"

        # Verify there is exactly one IdentityOperator (required for normalization)
        id_ops = filter(op -> op isa IdentityOperator, operators)
        @assert length(id_ops) == 1 "Must have exactly one IdentityOperator (found $(length(id_ops)))"

        new(objective, constraints, operators, localization, direction, precision)
    end
end

function Base.show(io::IO, sdp::SDPProblem)
    loc = sdp.localization
    print(io, "SDPProblem(objective=$(sdp.objective), " *
              "$(length(sdp.constraints)) constraints, " *
              "$(length(sdp.operators)) operators, " *
              "N=$(loc.N), k=$(loc.k), M=$(loc.M), " *
              "precision=$(sdp.precision), direction=$(sdp.direction))")
end

function Base.getproperty(sdp::SDPProblem, name::Symbol)
    if name == :functionals
        return [c.functional for c in sdp.constraints]
    elseif name == :normalization
        # normalization[j] = F_j(Id)
        #                   + λ²_Bp0020   · F_j(Bp0020)
        #                   + λ²_Bp0040   · F_j(Bp0040)
        #                   + λ²_Btwo0200 · F_j(Btwo0200)
        #                   - rhs_j
        # normalization[end] = 1  (α₀ coefficient)
        loc = getfield(sdp, :localization)
        normalization = zeros(BigFloat, length(sdp.constraints) + 1)
        for (idx, constraint) in enumerate(sdp.constraints)
            functional = constraint.functional
            rhs = constraint.rhs

            id_val = get_functional_value(functional, IdentityOperator())
            bp0020_val  = get_functional_value(functional, Bp0020Operator())
            bp0040_val  = get_functional_value(functional, Bp0040Operator())
            btwo0200_val = get_functional_value(functional, Btwo0200Operator())

            fixed_sum = (something(id_val,  BigFloat(0))
                       + loc.lambda2_Bp0020   * something(bp0020_val,  BigFloat(0))
                       + loc.lambda2_Bp0040   * something(bp0040_val,  BigFloat(0))
                       + loc.lambda2_Btwo0200 * something(btwo0200_val, BigFloat(0)))
            normalization[idx] = fixed_sum - rhs
        end
        normalization[end] = BigFloat(1)
        return normalization
    else
        return getfield(sdp, name)
    end
end

"""
    DualizedSDP

An SDP problem converted to SDPB's required form:
maximize alpha_0 subject to (alpha_0 - dual_obj) = 1 and PSD constraints ≥ 0

# Fields
- `alpha_index::Int`: Position of alpha_0 in the decision variable vector
- `normalization::Vector{BigFloat}`: Coefficients for normalization constraint
- `objective::Vector{BigFloat}`: Coefficients for objective function
- `constraint_matrices::Vector`: PSD constraint structures for SDPB
- `operator_values::Dict{Tuple{Functional,Operator},BigFloat}`: Cached functional values
- `precision::Int`: Precision in bits
"""
struct DualizedSDP
    alpha_index::Int
    normalization::Vector{BigFloat}
    objective::Vector{BigFloat}
    constraint_matrices::Vector
    operator_values::Dict{Tuple{Functional,Operator},BigFloat}
    precision::Int

    function DualizedSDP(alpha_index::Int,
                        normalization::Vector{BigFloat},
                        objective::Vector{BigFloat},
                        constraint_matrices::Vector,
                        operator_values::Dict{Tuple{Functional,Operator},BigFloat};
                        precision::Int=1024)
        @assert alpha_index >= 1 "alpha_index must be positive"
        @assert length(normalization) == length(objective) "normalization and objective must have same length"
        @assert !isempty(constraint_matrices) "Must have at least one constraint matrix"
        new(alpha_index, normalization, objective, constraint_matrices, operator_values, precision)
    end
end

function Base.show(io::IO, d::DualizedSDP)
    print(io, "DualizedSDP($(length(d.objective)) variables, $(length(d.constraint_matrices)) constraint matrices, precision=$(d.precision))")
end

"""
    SDPBResult

Result from solving an SDP with SDPB.

# Fields
- `status::String`: Termination status
- `objective_value::Union{BigFloat,Nothing}`: Raw alpha_0 value from SDPB
- `true_objective_value::Union{BigFloat,Nothing}`: Corrected dual objective
- `primal_solution::Union{Vector{BigFloat},Nothing}`: Primal solution y vector
- `dual_solution::Union{Vector{BigFloat},Nothing}`: Dual solution (if available)
- `output_dir::String`: Directory containing full SDPB output
- `solve_time::Float64`: Wall-clock solve time in seconds
- `z_solution::Union{Vector{BigFloat},Nothing}`: Full z solution vector
"""
struct SDPBResult
    status::String
    objective_value::Union{BigFloat,Nothing}
    true_objective_value::Union{BigFloat,Nothing}
    primal_solution::Union{Vector{BigFloat},Nothing}
    dual_solution::Union{Vector{BigFloat},Nothing}
    output_dir::String
    solve_time::Float64
    z_solution::Union{Vector{BigFloat},Nothing}
end

function Base.show(io::IO, r::SDPBResult)
    z_info = isnothing(r.z_solution) ? "no z" : "$(length(r.z_solution))-element z"
    print(io, "SDPBResult(status=\"$(r.status)\", true_objective=$(r.true_objective_value), time=$(r.solve_time)s, $z_info)")
end

# ============================================================================
# Dualization
# ============================================================================

"""
    dualize_sdp(sdp::SDPProblem) -> DualizedSDP

Convert an SDP problem to SDPB's required form.

All `FixedOperator` instances (Identity, Bp0020, Bp0040, Btwo0200) enter the
normalization vector and are skipped when building PSD constraints.
All `FreeOperator` instances (semishort and long) get their own 1×1 PSD constraint.
"""
function dualize_sdp(sdp::SDPProblem)
    @info "Loading functional values from cache..."
    operator_values = Dict{Tuple{Functional,Operator},BigFloat}()

    all_functionals = [sdp.objective; [c.functional for c in sdp.constraints]]

    for functional in all_functionals
        load_from_cache!(functional)
    end

    for functional in all_functionals
        for operator in sdp.operators
            key = (functional, operator)
            value = if haskey(OPERATOR_VALUES_CACHE, key)
                OPERATOR_VALUES_CACHE[key]
            else
                v = get_functional_value(functional, operator)
                if isnothing(v)
                    error("Functional value not cached for ($functional, $operator). " *
                          "Run compute_functional first.")
                end
                OPERATOR_VALUES_CACHE[key] = v
                v
            end
            operator_values[key] = value
        end
    end

    @info "Loaded $(length(operator_values)) functional values"

    n_constraints = length(sdp.constraints)
    n_vars = n_constraints + 1
    alpha_index = n_vars

    @info "Problem has $(n_constraints) constraints, $(n_vars) total decision variables (dual)"

    normalization = sdp.normalization

    objective = zeros(BigFloat, n_vars)
    objective[alpha_index] = BigFloat(-sdp.direction)

    constraint_matrices = []

    for operator in sdp.operators
        if operator isa FixedOperator
            # All fixed operators (Identity, Bp0020, Bp0040, Btwo0200) go into
            # the normalization vector; they do not generate PSD constraints.
            continue
        end

        push!(constraint_matrices, build_dual_operator_constraint(
            operator, sdp.constraints, sdp.objective, operator_values, normalization, sdp.direction))
    end

    @info "Built $(length(constraint_matrices)) constraint matrices"

    return DualizedSDP(alpha_index, normalization, objective, constraint_matrices,
                      operator_values, precision=sdp.precision)
end

"""
    build_dual_operator_constraint(operator, constraints, objective, operator_values,
                                   normalization, direction) -> Dict

Build a 1×1 dual feasibility PSD constraint for a single `FreeOperator`.
"""
function build_dual_operator_constraint(operator::Operator,
                                       constraints::Vector{SDPConstraint},
                                       objective::Functional,
                                       operator_values::Dict,
                                       normalization::Vector{BigFloat},
                                       direction::Int)
    coefficients = []

    for constraint in constraints
        functional = constraint.functional
        value = operator_values[(functional, operator)]
        push!(coefficients, value)
    end
    push!(coefficients, 0)

    obj_value = operator_values[(objective, operator)]
    coefficients = direction .* (coefficients .- (obj_value .* normalization))

    polynomials = [[[string.(coefficients)]]]

    return Dict(
        "DampedRational" => Dict(
            "base" => "0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527684399520824697579279012900862665358949409878309219436737733811504863899112514561634498772",
            "constant" => "1",
            "poles" => []
        ),
        "polynomials" => polynomials
    )
end

# ============================================================================
# File Writing
# ============================================================================

"""
    write_sdp_json(dualized::DualizedSDP, output_path::String)

Write a dualized SDP to a JSON file in SDPB's PMP format.
"""
function write_sdp_json(dualized::DualizedSDP, output_path::String)
    @assert endswith(output_path, ".json") "Output path must end with .json"

    json_dict = Dict(
        "objective" => [string(x) for x in dualized.objective],
        "normalization" => [string(x) for x in dualized.normalization],
        "PositiveMatrixWithPrefactorArray" => dualized.constraint_matrices
    )

    open(output_path, "w") do f
        JSON.print(f, json_dict, 4)
    end

    @info "Wrote SDP to $output_path"
end

format_fixed(x::Real, prec::Int = ceil(precision(BigFloat(x)))) = begin
    bx = BigFloat(x)
    if bx == 0
        return "0"
    end
    digits = max(1, Int(ceil(prec * log10(2))))
    return @sprintf("%.*f", digits, bx)
end

"""
    write_sdp_mathematica(dualized::DualizedSDP, output_path::String)

Write a dualized SDP to a Mathematica (.m) file in SDPB's format.
"""
function write_sdp_mathematica(dualized::DualizedSDP, output_path::String)
    @assert endswith(output_path, ".m") "Output path must end with .m"

    open(output_path, "w") do f
        print(f, "SDP[")

        print(f, "{")
        for (i, val) in enumerate(dualized.objective)
            print(f, format_fixed(val, dualized.precision))
            if i < length(dualized.objective)
                print(f, ", ")
            end
        end
        print(f, "}, ")

        print(f, "{")
        for (i, val) in enumerate(dualized.normalization)
            print(f, format_fixed(val, dualized.precision))
            if i < length(dualized.normalization)
                print(f, ", ")
            end
        end
        print(f, "}, ")
        println(f)

        println(f, "  {")
        for (idx, matrix) in enumerate(dualized.constraint_matrices)
            write_mathematica_matrix(f, matrix, indent="    ", prec=dualized.precision)
            if idx < length(dualized.constraint_matrices)
                println(f, ",")
            else
                println(f)
            end
        end
        println(f, "  }")
        println(f, "]")
    end

    @info "Wrote SDP to $output_path"
end

"""
    write_mathematica_matrix(io::IO, matrix::Dict; indent::String="", prec::Int=256)

Helper function to write a constraint matrix in Mathematica format.
"""
function write_mathematica_matrix(io::IO, matrix::Dict; indent::String="", prec::Int=256)
    damped_rational = matrix["DampedRational"]
    polynomials = matrix["polynomials"]

    println(io, indent, "PositiveMatrixWithPrefactor[")

    print(io, indent, "  DampedRational[")
    print(io, damped_rational["constant"], ", ")
    print(io, "{")
    poles = damped_rational["poles"]
    for (i, pole) in enumerate(poles)
        print(io, pole)
        if i < length(poles)
            print(io, ", ")
        end
    end
    print(io, "}, ")
    print(io, damped_rational["base"], ", x]")
    println(io, ",")

    for (row_idx, row) in enumerate(polynomials)
        print(io, indent, "    {")
        for (col_idx, col) in enumerate(row)
            print(io, "{")
            for (sample_idx, sample) in enumerate(col)
                print(io, "{")
                for (coeff_idx, coeff) in enumerate(sample)
                    print(io, format_fixed(parse(BigFloat, coeff), prec))
                    if coeff_idx < length(sample)
                        print(io, ", ")
                    end
                end
                print(io, "}")
                if sample_idx < length(col)
                    print(io, ", ")
                end
            end
            print(io, "}")
            if col_idx < length(row)
                print(io, ", ")
            end
        end
        print(io, "}")
        if row_idx < length(polynomials)
            println(io, ",")
        else
            println(io)
        end
    end
    print(io, "]")
end

# ============================================================================
# SDPB Integration — Command Builders
# ============================================================================

"""
Build a cleanup command that removes `full_work_dir`.

When using Docker, files in the work directory are owned by root, so
we run `rm` inside the container using the same volume mount.
"""
function _build_clean_cmd(full_work_dir::String,
                          sdpb_build_dir::String,
                          sdpb_docker_image::String)
    if sdpb_build_dir == "docker"
        work_dir  = dirname(full_work_dir)
        hash_name = basename(full_work_dir)
        return "docker run -v $work_dir:/usr/local/share/sdpb $sdpb_docker_image rm -rf /usr/local/share/sdpb/$hash_name"
    else
        return "rm -rf $full_work_dir"
    end
end

"""
    build_sdp2input_command(pmp_path, sdp_dir, precision; ...) -> String

Build the command string for running the sdp2input preprocessor.
Supports both local build and Docker execution.
"""
function build_sdp2input_command(pmp_path::String, sdp_dir::String, precision::Int;
                                 sdpb_build_dir::Union{String,Nothing}=nothing,
                                 sdpb_docker_image::Union{String,Nothing}=nothing,
                                 use_slurm::Bool=false)
    config = get_config()
    sdpb_build_dir    = isnothing(sdpb_build_dir)    ? config.sdpb_build_dir    : sdpb_build_dir
    sdpb_docker_image = isnothing(sdpb_docker_image) ? config.sdpb_docker_image : sdpb_docker_image

    if sdpb_build_dir == "docker"
        local_dir = dirname(pmp_path)
        remap(p)  = joinpath("/usr/local/share/sdpb", relpath(p, local_dir))
        return join([
            "docker", "run",
            "-v", "$local_dir:/usr/local/share/sdpb",
            sdpb_docker_image,
            "sdp2input",
            "--precision=$precision",
            "--input=$(remap(pmp_path))",
            "--output=$(remap(sdp_dir))",
        ], " ")
    else
        exe = joinpath(sdpb_build_dir, "sdp2input")
        mpi_parts = use_slurm ? ["srun"] : ["mpirun", "-n", 1]
        return join([mpi_parts..., exe, "--precision=$precision", "--input=$pmp_path", "--output=$sdp_dir"], " ")
    end
end

"""
    build_sdpb_command(sdp_dir, output_dir, precision; ...) -> String

Build the command string for running the SDPB solver.
Supports both local build and Docker execution.
"""
function build_sdpb_command(sdp_dir::String, output_dir::String, precision::Int;
                            checkpoint_dir::Union{String,Nothing}=nothing,
                            sdpb_build_dir::Union{String,Nothing}=nothing,
                            sdpb_docker_image::Union{String,Nothing}=nothing,
                            use_slurm::Bool=false,
                            n_cores::Union{Int,Nothing}=nothing,
                            extra_args::Vector{String}=String[])
    config = get_config()
    sdpb_build_dir    = isnothing(sdpb_build_dir)    ? config.sdpb_build_dir    : sdpb_build_dir
    sdpb_docker_image = isnothing(sdpb_docker_image) ? config.sdpb_docker_image : sdpb_docker_image

    if sdpb_build_dir == "docker"
        local_dir = dirname(sdp_dir)
        remap(p)  = joinpath("/usr/local/share/sdpb", relpath(p, local_dir))

        mpi_parts = ["mpirun", "--allow-run-as-root", "-n", isnothing(n_cores) ? "1" : string(n_cores)]

        cmd_parts = [
            "docker", "run",
            "-v", "$local_dir:/usr/local/share/sdpb",
            sdpb_docker_image,
            mpi_parts...,
            "sdpb",
            "--precision=$precision",
            "-s $(remap(sdp_dir))",
            "-o $(remap(output_dir))",
            "--writeSolution=z",
        ]

        if !isnothing(checkpoint_dir)
            push!(cmd_parts, "-c $(remap(checkpoint_dir))")
        end
        append!(cmd_parts, extra_args)

        if use_slurm
            cmd_parts = ["srun", "--export=ALL", cmd_parts...]
        end
    else
        mpi_parts = use_slurm ? [] : ["mpirun", "-n", isnothing(n_cores) ? "1" : string(n_cores)]
        cmd_parts = [
            mpi_parts...,
            joinpath(sdpb_build_dir, "sdpb"),
            "--precision=$precision",
            "-s $sdp_dir",
            "-o $output_dir",
            "--writeSolution=z",
        ]

        if use_slurm
            cmd_parts = ["srun", "--export=ALL", cmd_parts...]
        end
        if !isnothing(checkpoint_dir)
            push!(cmd_parts, "-c $checkpoint_dir")
        end
        append!(cmd_parts, extra_args)
    end

    return join(filter(!isempty, cmd_parts), " ")
end

# ============================================================================
# Result Parsing
# ============================================================================

"""
    parse_z_txt(z_file::String) -> Vector{BigFloat}

Parse SDPB's `z.txt` solution file into a BigFloat vector.
"""
function parse_z_txt(z_file::String)
    isfile(z_file) || error("z.txt not found: $z_file")
    lines = filter(!isempty, strip.(split(read(z_file, String), '\n')))
    isempty(lines) && error("z.txt is empty: $z_file")
    header = split(lines[1])
    height = parse(Int, header[1])
    values = BigFloat[]
    for line in @view lines[2:end]
        for tok in split(line)
            push!(values, parse(BigFloat, tok))
        end
    end
    length(values) == height || error("z.txt: expected $height values, got $(length(values))")
    return values
end

"""
    parse_out_txt(output_file, direction, fixed_contribution)

Low-level parser for a single SDPB `out.txt` file.

`fixed_contribution` is the value of the objective functional summed over all
fixed operators (weighted by their localization OPE² coefficients), which is
excluded from the SDP and must be added by hand to recover the true dual bound.
"""
function parse_out_txt(output_file::String, direction::Int, fixed_contribution::Real)
    isfile(output_file) || error("Output file not found: $output_file")
    content = read(output_file, String)

    status_match = match(r"terminateReason = \"([^\"]+)\"", content)
    status = isnothing(status_match) ? "unknown" : status_match.captures[1]

    primal_obj_match = match(r"primalObjective = ([^;]+);", content)
    objective_value = isnothing(primal_obj_match) ? nothing :
        parse(BigFloat, strip(primal_obj_match.captures[1]))

    true_objective_value = isnothing(objective_value) ? nothing :
        objective_value * (-direction) - BigFloat(1) + BigFloat(fixed_contribution)

    y_match = match(r"y = \{([^}]+)\}", content)
    primal_solution = isnothing(y_match) ? nothing :
        [parse(BigFloat, strip(s)) for s in split(y_match.captures[1], ",")]

    time_match = match(r"Solver runtime  = ([0-9.]+)", content)
    solve_time = isnothing(time_match) ? 0.0 : parse(Float64, time_match.captures[1])

    return (status=status, objective_value=objective_value,
            true_objective_value=true_objective_value,
            primal_solution=primal_solution, solve_time=solve_time)
end

"""
    parse_sdpb_output(output_dir::String, sdp::SDPProblem; read_z::Bool=false) -> SDPBResult

Parse SDPB output files to extract the solution.

The true dual objective is:
    alpha_0 · (-direction) - 1
    + F_obj(Id)
    + λ²_Bp0020   · F_obj(Bp0020)
    + λ²_Bp0040   · F_obj(Bp0040)
    + λ²_Btwo0200 · F_obj(Btwo0200)
"""
function parse_sdpb_output(output_dir::String, sdp::SDPProblem; read_z::Bool=false)
    loc = sdp.localization

    id_val       = get_functional_value(sdp.objective, IdentityOperator())
    bp0020_val   = get_functional_value(sdp.objective, Bp0020Operator())
    bp0040_val   = get_functional_value(sdp.objective, Bp0040Operator())
    btwo0200_val = get_functional_value(sdp.objective, Btwo0200Operator())

    fixed_contribution = (something(id_val, BigFloat(0))
                        + loc.lambda2_Bp0020   * something(bp0020_val,   BigFloat(0))
                        + loc.lambda2_Bp0040   * something(bp0040_val,   BigFloat(0))
                        + loc.lambda2_Btwo0200 * something(btwo0200_val, BigFloat(0)))

    parsed = parse_out_txt(joinpath(output_dir, "out.txt"), sdp.direction, fixed_contribution)
    @info "Parsed SDPB output: status=$(parsed.status), true_objective=$(parsed.true_objective_value)"

    z_solution = if read_z
        z_file = joinpath(output_dir, "z.txt")
        if isfile(z_file)
            try
                parse_z_txt(z_file)
            catch e
                @warn "Failed to parse z.txt: $e"
                nothing
            end
        else
            @warn "read_z=true but z.txt not found in $output_dir"
            nothing
        end
    else
        nothing
    end

    return SDPBResult(
        parsed.status,
        parsed.objective_value,
        parsed.true_objective_value,
        parsed.primal_solution,
        nothing,
        output_dir,
        parsed.solve_time,
        z_solution
    )
end

# ============================================================================
# Problem Metadata Helpers
# ============================================================================

"""
    extract_NkM(sdp::SDPProblem) -> Tuple{Int,Int,Int}

Extract the (N, k, M) ABJM theory point from the problem's localization data.
"""
function extract_NkM(sdp::SDPProblem)
    loc = sdp.localization
    return (loc.N, loc.k, loc.M)
end

"""
    extract_Lambda(sdp::SDPProblem) -> Union{Int,Nothing}

Infer the crossing cutoff Λ from the crossing constraints in the problem.

For ABJM, Λ is defined so that the constraint set equals exactly the union of:
- Channel 1: {(m, 0) : m even, 0 ≤ m ≤ Λ}
- Channel 2: {(m, n) : m+n even, 0 ≤ n ≤ m, m+n ≤ Λ}

Returns `nothing` if no crossing constraints are present or if the set does
not form a complete block up to any Λ.
"""
function extract_Lambda(sdp::SDPProblem)
    crossing_pairs_ch1 = Set{Tuple{Int,Int}}()
    crossing_pairs_ch2 = Set{Tuple{Int,Int}}()

    for c in sdp.constraints
        if c.functional isa CrossingFunctional
            f = c.functional
            if f.cross == 1
                push!(crossing_pairs_ch1, (f.m, f.n))
            else
                push!(crossing_pairs_ch2, (f.m, f.n))
            end
        end
    end

    isempty(crossing_pairs_ch1) && isempty(crossing_pairs_ch2) && return nothing

    all_pairs = union(crossing_pairs_ch1, crossing_pairs_ch2)
    Λ_candidate = maximum(m + n for (m, n) in all_pairs)

    # Channel 1: (m, 0) with m even, m ∈ [0, Λ]
    expected_ch1 = Set{Tuple{Int,Int}}()
    for m in 0:2:Λ_candidate
        push!(expected_ch1, (m, 0))
    end

    # Channel 2: (m, n) with m+n even, 0 ≤ n ≤ m, m+n ≤ Λ
    expected_ch2 = Set{Tuple{Int,Int}}()
    for m in 0:Λ_candidate
        for n in 0:min(m, Λ_candidate - m)
            iseven(m + n) && push!(expected_ch2, (m, n))
        end
    end

    return (crossing_pairs_ch1 == expected_ch1 && crossing_pairs_ch2 == expected_ch2) ?
           Λ_candidate : nothing
end

"""
    functional_type_key(f::Functional) -> String

Short snake_case label for a functional type.
"""
function functional_type_key(f::Functional)
    if f isa CrossingFunctional;              return "crossing"
    elseif f isa IntegralFunctional;          return "integral"
    elseif f isa LinearCombinationObjective;  return "linear_combination"
    else;                                     return "other"
    end
end

"""
    operator_key_string(op::FreeOperator) -> String

Short identifier for a free operator, used as a metadata key suffix.
"""
operator_key_string(op::SemishortAp0020)   = "SemishortAp0020_J$(op.J)"
operator_key_string(op::SemishortAtwo0100) = "SemishortAtwo0100_J$(op.J)"
operator_key_string(op::LongOperator)      = "Long_tau$(op.tau)_J$(op.J)"

"""
    functional_properties(f::Functional) -> Dict{String,Any}

Return a dict of the defining (non-precision) parameters of a functional.
"""
function functional_properties(f::Functional)
    props = Dict{String, Any}()
    if f isa CrossingFunctional
        props["m"]     = f.m
        props["n"]     = f.n
        props["cross"] = f.cross
    elseif f isa IntegralFunctional
        props["b"] = f.b
    elseif f isa LinearCombinationObjective
        for (op, coeff) in f.operator_coefficients
            props["coeff_" * operator_key_string(op)] = coeff
        end
    end
    return props
end

"""
    generate_problem_hash(sdp::SDPProblem) -> String

Generate a unique 16-character hex hash for an SDP problem.
"""
function generate_problem_hash(sdp::SDPProblem)
    parts = String[]

    push!(parts, "objective:$(sdp.objective)")

    constraint_strs = [string(c.functional, "_", c.comparison, "_", c.rhs) for c in sdp.constraints]
    sort!(constraint_strs)
    push!(parts, "constraints:" * join(constraint_strs, ","))

    operator_strs = String[]
    for op in sdp.operators
        if op isa IdentityOperator
            push!(operator_strs, "Identity")
        elseif op isa Bp0020Operator
            push!(operator_strs, "Bp0020")
        elseif op isa Bp0040Operator
            push!(operator_strs, "Bp0040")
        elseif op isa Btwo0200Operator
            push!(operator_strs, "Btwo0200")
        elseif op isa SemishortAp0020
            push!(operator_strs, "SemishortAp0020_J$(op.J)")
        elseif op isa SemishortAtwo0100
            push!(operator_strs, "SemishortAtwo0100_J$(op.J)")
        elseif op isa LongOperator
            push!(operator_strs, "Long_tau$(op.tau)_J$(op.J)")
        end
    end
    sort!(operator_strs)
    push!(parts, "operators:" * join(operator_strs, ","))

    loc = sdp.localization
    push!(parts, "localization:N$(loc.N)_k$(loc.k)_M$(loc.M)")

    push!(parts, "precision:$(sdp.precision)")
    push!(parts, "direction:$(sdp.direction)")

    full_string = join(parts, "|")
    hash_bytes = sha256(full_string)
    hash_hex = bytes2hex(hash_bytes)
    return hash_hex[1:16]
end

"""
    write_problem_metadata(sdp::SDPProblem, filepath::String)

Write a JSON metadata file describing the SDP problem.
"""
function write_problem_metadata(sdp::SDPProblem, filepath::String)
    N, k, M = extract_NkM(sdp)
    Λ_val   = extract_Lambda(sdp)
    loc     = sdp.localization

    # Compute the fixed-operator contribution to the objective (excluded from
    # the SDP; re-added when computing the true dual bound).
    id_val       = get_functional_value(sdp.objective, IdentityOperator())
    bp0020_val   = get_functional_value(sdp.objective, Bp0020Operator())
    bp0040_val   = get_functional_value(sdp.objective, Bp0040Operator())
    btwo0200_val = get_functional_value(sdp.objective, Btwo0200Operator())

    fixed_contribution = Float64(
        something(id_val, BigFloat(0))
        + loc.lambda2_Bp0020   * something(bp0020_val,   BigFloat(0))
        + loc.lambda2_Bp0040   * something(bp0040_val,   BigFloat(0))
        + loc.lambda2_Btwo0200 * something(btwo0200_val, BigFloat(0))
    )

    type_counts = Dict{String, Int}()
    for c in sdp.constraints
        key = "num_" * functional_type_key(c.functional)
        type_counts[key] = get(type_counts, key, 0) + 1
    end

    metadata = Dict{String, Any}(
        "timestamp"             => string(now()),
        "precision"             => sdp.precision,
        "direction"             => sdp.direction,
        "fixed_contribution"    => fixed_contribution,
        "N"                     => N,
        "k"                     => k,
        "M"                     => M,
        "cT"                    => Float64(loc.cT),
        "Lambda"                => Λ_val,
        "objective"             => string(sdp.objective),
        "objective_properties"  => functional_properties(sdp.objective),
        "num_constraints"       => length(sdp.constraints),
        "num_operators"         => length(sdp.operators),
        type_counts...,
        "constraints" => [
            Dict(
                "functional" => string(c.functional),
                "comparison" => string(c.comparison),
                "rhs"        => string(c.rhs)
            ) for c in sdp.constraints
        ],
        "operators" => [
            if op isa IdentityOperator
                Dict("type" => "Identity")
            elseif op isa Bp0020Operator
                Dict("type" => "Bp0020")
            elseif op isa Bp0040Operator
                Dict("type" => "Bp0040")
            elseif op isa Btwo0200Operator
                Dict("type" => "Btwo0200")
            elseif op isa SemishortAp0020
                Dict("type" => "SemishortAp0020", "J" => op.J)
            elseif op isa SemishortAtwo0100
                Dict("type" => "SemishortAtwo0100", "J" => op.J)
            elseif op isa LongOperator
                Dict("type" => "Long", "tau" => op.tau, "J" => op.J)
            else
                Dict("type" => "Unknown")
            end
            for op in sort(sdp.operators)
        ]
    )

    open(filepath, "w") do f
        JSON.print(f, metadata, 4)
    end

    @info "Wrote problem metadata to $filepath"
end

# ============================================================================
# Result Storage Helpers
# ============================================================================

"""
    get_result_paths(problem_hash, results_dir=nothing; plot_name=nothing)

Return `(result_path, metadata_path)` for a given problem hash.
"""
function get_result_paths(problem_hash::String,
                          results_dir::Union{String,Nothing}=nothing;
                          plot_name::Union{String,Nothing}=nothing)
    results_dir = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    dir = isnothing(plot_name) ? results_dir : joinpath(results_dir, plot_name)
    result_path   = joinpath(dir, "$(problem_hash)_out.txt")
    metadata_path = joinpath(dir, "$(problem_hash)_metadata.json")
    return (result_path, metadata_path)
end

"""
    get_z_result_path(problem_hash, results_dir=nothing; plot_name=nothing)

Return the path for the saved `z.txt` solution file.
"""
function get_z_result_path(problem_hash::String,
                            results_dir::Union{String,Nothing}=nothing;
                            plot_name::Union{String,Nothing}=nothing)
    results_dir = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    dir = isnothing(plot_name) ? results_dir : joinpath(results_dir, plot_name)
    return joinpath(dir, "$(problem_hash)_z.txt")
end

"""
    check_result_exists(problem_hash, results_dir=nothing; plot_name=nothing) -> Bool

Check whether a result file already exists for this problem hash.
"""
function check_result_exists(problem_hash::String,
                              results_dir::Union{String,Nothing}=nothing;
                              plot_name::Union{String,Nothing}=nothing)
    result_path, _ = get_result_paths(problem_hash, results_dir; plot_name=plot_name)
    return isfile(result_path)
end

# ============================================================================
# High-Level Workflow
# ============================================================================

"""
    solve_sdp(sdp::SDPProblem; ...) -> Union{SDPBResult, Nothing}

Complete workflow: check cache → dualize → write PMP → sdp2input → sdpb → parse.

If `slurm_config` is provided, the job is submitted asynchronously and
`nothing` is returned. If results are already cached, they are loaded and
returned immediately.

# Keyword Arguments
- `plots`: `Plot`, vector of `Plot`s, or `nothing` for flat storage
- `work_dir`: Working directory for intermediate files (default: from config)
- `sdpb_build_dir`: SDPB build directory or `"docker"` (default: from config)
- `sdpb_docker_image`: Docker image when `sdpb_build_dir="docker"` (default: from config)
- `slurm_config`: `SlurmConfig` or `nothing` for direct execution (default: `SlurmConfig()`)
- `keep_intermediate`: Retain intermediate files after solving (default: false)
- `save_result`: Copy result to `results_dir` (default: true)
- `extra_sdpb_args`: Additional arguments passed to sdpb (default: `[]`)
- `results_dir`: Permanent results directory (default: from config)
- `n_cores`: Number of MPI ranks (default: nothing → SDPB chooses)
"""
function solve_sdp(sdp::SDPProblem;
                   plots::Union{Plot, JointBoundsPlot, Vector, Nothing}=nothing,
                   work_dir::Union{String,Nothing}=nothing,
                   sdpb_build_dir::Union{String,Nothing}=nothing,
                   sdpb_docker_image::Union{String,Nothing}=nothing,
                   slurm_config::Union{SlurmConfig,Nothing}=SlurmConfig(),
                   keep_intermediate::Bool=false,
                   save_result::Bool=true,
                   extra_sdpb_args::Vector{String}=String[],
                   results_dir::Union{String,Nothing}=nothing,
                   n_cores::Union{Int,Nothing}=nothing)
    config = get_config()
    work_dir          = isnothing(work_dir)          ? config.sdpb_work_dir     : work_dir
    sdpb_build_dir    = isnothing(sdpb_build_dir)    ? config.sdpb_build_dir    : sdpb_build_dir
    sdpb_docker_image = isnothing(sdpb_docker_image) ? config.sdpb_docker_image : sdpb_docker_image
    results_dir       = isnothing(results_dir)       ? config.sdpb_results_dir  : results_dir

    plot_names = if isnothing(plots)
        [nothing]
    elseif plots isa Union{Plot, JointBoundsPlot}
        [plots.name]
    else
        [p.name for p in plots]
    end

    @info "Starting SDP solve workflow..."
    problem_hash  = generate_problem_hash(sdp)
    full_work_dir = joinpath(work_dir, problem_hash)
    @info "Problem hash: $problem_hash"
    @info "Working directory: $full_work_dir"

    # ------------------------------------------------------------------ #
    # Step 0: Check for cached results                                     #
    # ------------------------------------------------------------------ #
    @info "Step 0/5: Checking for cached results..."

    cached_result_path = nothing
    for pn in plot_names
        rp, _ = get_result_paths(problem_hash, results_dir; plot_name=pn)
        if isfile(rp) && !isempty(read(rp, String))
            cached_result_path = rp
            break
        end
    end

    if !isnothing(cached_result_path)
        temp_result_dir = joinpath(full_work_dir, "cached_result")
        mkpath(temp_result_dir)
        cp(cached_result_path, joinpath(temp_result_dir, "out.txt"), force=true)
        result = parse_sdpb_output(temp_result_dir, sdp, read_z=false)

        if result.status != "found primal-dual optimal solution"
            @warn "Cached result at $cached_result_path has status $(result.status). Re-submitting."
            cached_result_path = nothing
            rm(temp_result_dir, force=true, recursive=true)
        else
            @info "✓ Found cached result at $cached_result_path"

            for pn in plot_names
                zp = get_z_result_path(problem_hash, results_dir; plot_name=pn)
                if isfile(zp)
                    cp(zp, joinpath(temp_result_dir, "z.txt"), force=true)
                    break
                end
            end
            result = parse_sdpb_output(temp_result_dir, sdp, read_z=true)
            rm(temp_result_dir, force=true, recursive=true)

            if save_result
                for pn in plot_names
                    rp, mp = get_result_paths(problem_hash, results_dir; plot_name=pn)
                    d = dirname(rp)
                    isdir(d) || mkpath(d)
                    isfile(rp) || cp(cached_result_path, rp, force=true)
                    isfile(mp) || write_problem_metadata(sdp, mp)
                end
            end

            @info "Loaded cached result: status=$(result.status), objective=$(result.true_objective_value)"
            return result
        end
    else
        @info "No cached result found, will submit new job"
    end

    # ------------------------------------------------------------------ #
    # Create work directory; write metadata                                #
    # ------------------------------------------------------------------ #
    mkpath(full_work_dir)

    if save_result
        for pn in plot_names
            rp, mp = get_result_paths(problem_hash, results_dir; plot_name=pn)
            d = dirname(rp)
            isdir(d) || mkpath(d)
            @info "Writing metadata to $mp..."
            write_problem_metadata(sdp, mp)
        end
    end

    # ------------------------------------------------------------------ #
    # Step 1: Dualize                                                       #
    # ------------------------------------------------------------------ #
    @info "Step 1/5: Dualizing SDP problem..."
    dualized = dualize_sdp(sdp)

    # ------------------------------------------------------------------ #
    # Step 2: Write PMP                                                    #
    # ------------------------------------------------------------------ #
    pmp_path = joinpath(full_work_dir, "sdp.m")
    @info "Step 2/5: Writing PMP to $pmp_path..."
    write_sdp_mathematica(dualized, pmp_path)

    # ------------------------------------------------------------------ #
    # Step 3: Build commands                                               #
    # ------------------------------------------------------------------ #
    sdp_dir        = joinpath(full_work_dir, "sdp_preprocessed")
    output_dir     = joinpath(full_work_dir, "output")
    checkpoint_dir = joinpath(full_work_dir, "checkpoint")

    @info "Step 3/5: Preparing combined sdp2input + sdpb job..."

    sdp2input_cmd = build_sdp2input_command(pmp_path, sdp_dir, sdp.precision,
                                            sdpb_build_dir=sdpb_build_dir,
                                            sdpb_docker_image=sdpb_docker_image)
    mkpath(checkpoint_dir)

    sdpb_cmd = build_sdpb_command(sdp_dir, output_dir, sdp.precision,
                                  checkpoint_dir=checkpoint_dir,
                                  sdpb_build_dir=sdpb_build_dir,
                                  sdpb_docker_image=sdpb_docker_image,
                                  use_slurm=(!isnothing(slurm_config)),
                                  extra_args=extra_sdpb_args,
                                  n_cores=n_cores)

    clean_cmd = _build_clean_cmd(full_work_dir, sdpb_build_dir, sdpb_docker_image)

    out_txt_path = joinpath(output_dir, "out.txt")
    z_txt_path   = joinpath(output_dir, "z.txt")
    post_cmds = String[]
    if save_result
        for pn in plot_names
            rp, _ = get_result_paths(problem_hash, results_dir; plot_name=pn)
            zp    = get_z_result_path(problem_hash, results_dir; plot_name=pn)
            push!(post_cmds, "cp $out_txt_path $rp")
            push!(post_cmds, "cp $z_txt_path $zp 2>/dev/null || true")
        end
    end

    combined_cmd = "$sdp2input_cmd && $sdpb_cmd"
    if !isempty(post_cmds)
        combined_cmd *= " && " * join(post_cmds, " && ")
    end

    # ------------------------------------------------------------------ #
    # Step 4: Submit / execute                                             #
    # ------------------------------------------------------------------ #
    @info "Step 4/5: Submitting combined sdp2input + sdpb job..."

    if !isnothing(slurm_config)
        post_cmd = join(post_cmds, " && ") * (!keep_intermediate ? "\n$clean_cmd" : "")
        job_id = submit_slurm_job(slurm_config, combined_cmd, wait=false, post_command=post_cmd)
        @info "Job submitted with ID: $job_id"
        for pn in plot_names
            rp, mp = get_result_paths(problem_hash, results_dir; plot_name=pn)
            zp     = get_z_result_path(problem_hash, results_dir; plot_name=pn)
            @info "Results will be copied to: $rp"
            @info "z solution will be copied to: $zp"
            @info "Metadata available at: $mp"
        end
        return nothing
    else
        @info "Running sdp2input + sdpb directly..."
        run(`sh -c $combined_cmd`)
        result = parse_sdpb_output(output_dir, sdp, read_z=true)
        if !keep_intermediate
            run(`sh -c $clean_cmd`)
        end
        return result
    end
end

"""
    solve_sdp(sdps::Vector{<:SDPProblem}; ...) -> Vector{Union{SDPBResult,Nothing}}

Batch solve multiple SDP problems, reusing cached results where available.
"""
function solve_sdp(sdps::Vector{<:SDPProblem};
                   plots::Union{Plot, JointBoundsPlot, Vector, Nothing}=nothing,
                   work_dir::Union{String,Nothing}=nothing,
                   sdpb_build_dir::Union{String,Nothing}=nothing,
                   sdpb_docker_image::Union{String,Nothing}=nothing,
                   slurm_config::Union{SlurmConfig,Nothing}=SlurmConfig(),
                   keep_intermediate::Bool=false,
                   save_result::Bool=true,
                   extra_sdpb_args::Vector{String}=String[],
                   results_dir::Union{String,Nothing}=nothing,
                   n_cores::Union{Int,Nothing}=nothing,
                   batch_size::Int=1,
                   verbosity::Int=2)
    config            = get_config()
    work_dir          = isnothing(work_dir)          ? config.sdpb_work_dir     : work_dir
    sdpb_build_dir    = isnothing(sdpb_build_dir)    ? config.sdpb_build_dir    : sdpb_build_dir
    sdpb_docker_image = isnothing(sdpb_docker_image) ? config.sdpb_docker_image : sdpb_docker_image
    results_dir       = isnothing(results_dir)       ? config.sdpb_results_dir  : results_dir

    plot_names = if isnothing(plots)
        [nothing]
    elseif plots isa Union{Plot, JointBoundsPlot}
        [plots.name]
    else
        [p.name for p in plots]
    end

    results       = Vector{Union{Nothing, SDPBResult}}(nothing, length(sdps))
    to_solve_indices = Int[]
    output_dirs   = String[]
    combined_cmds = String[]
    clean_cmds    = String[]

    for (i, sdp) in enumerate(sdps)
        problem_hash  = generate_problem_hash(sdp)
        full_work_dir = joinpath(work_dir, problem_hash)

        cached_result_path = nothing
        for pn in plot_names
            rp, _ = get_result_paths(problem_hash, results_dir; plot_name=pn)
            if isfile(rp) && !isempty(read(rp, String))
                cached_result_path = rp
                break
            end
        end

        if !isnothing(cached_result_path)
            verbosity >= 2 && @info "Problem $i: ✓ Found cached result at $cached_result_path"
            temp_result_dir = joinpath(full_work_dir, "cached_result")
            mkpath(temp_result_dir)
            cp(cached_result_path, joinpath(temp_result_dir, "out.txt"), force=true)
            for pn in plot_names
                zp = get_z_result_path(problem_hash, results_dir; plot_name=pn)
                if isfile(zp)
                    cp(zp, joinpath(temp_result_dir, "z.txt"), force=true)
                    break
                end
            end
            results[i] = parse_sdpb_output(temp_result_dir, sdp, read_z=true)
            rm(temp_result_dir, force=true, recursive=true)
            if save_result
                for pn in plot_names
                    rp, mp = get_result_paths(problem_hash, results_dir; plot_name=pn)
                    d = dirname(rp)
                    isdir(d) || mkpath(d)
                    isfile(rp) || cp(cached_result_path, rp, force=true)
                    isfile(mp) || write_problem_metadata(sdp, mp)
                end
            end
            continue
        end

        verbosity >= 2 && @info "Problem $i: No cached result found, will solve"
        mkpath(full_work_dir)
        if save_result
            for pn in plot_names
                rp, mp = get_result_paths(problem_hash, results_dir; plot_name=pn)
                d = dirname(rp)
                isdir(d) || mkpath(d)
                write_problem_metadata(sdp, mp)
            end
        end

        dualized = dualize_sdp(sdp)
        pmp_path  = joinpath(full_work_dir, "sdp.m")
        write_sdp_mathematica(dualized, pmp_path)

        sdp_dir        = joinpath(full_work_dir, "sdp_preprocessed")
        output_dir_i   = joinpath(full_work_dir, "output")
        checkpoint_dir = joinpath(full_work_dir, "checkpoint")

        sdp2input_cmd = build_sdp2input_command(pmp_path, sdp_dir, sdp.precision,
                                                sdpb_build_dir=sdpb_build_dir,
                                                sdpb_docker_image=sdpb_docker_image)
        mkpath(checkpoint_dir)
        sdpb_cmd = build_sdpb_command(sdp_dir, output_dir_i, sdp.precision,
                                      checkpoint_dir=checkpoint_dir,
                                      sdpb_build_dir=sdpb_build_dir,
                                      sdpb_docker_image=sdpb_docker_image,
                                      use_slurm=(!isnothing(slurm_config)),
                                      extra_args=extra_sdpb_args,
                                      n_cores=n_cores)

        out_txt_path_i = joinpath(output_dir_i, "out.txt")
        z_txt_path_i   = joinpath(output_dir_i, "z.txt")
        post_cmds_i    = String[]
        if save_result
            for pn in plot_names
                rp, _ = get_result_paths(problem_hash, results_dir; plot_name=pn)
                zp    = get_z_result_path(problem_hash, results_dir; plot_name=pn)
                push!(post_cmds_i, "cp $out_txt_path_i $rp")
                push!(post_cmds_i, "cp $z_txt_path_i $zp 2>/dev/null || true")
            end
        end

        combined_cmd_i = "$sdp2input_cmd && $sdpb_cmd"
        if !isempty(post_cmds_i)
            combined_cmd_i *= " && " * join(post_cmds_i, " && ")
        end

        push!(to_solve_indices, i)
        push!(output_dirs, output_dir_i)
        push!(combined_cmds, combined_cmd_i)
        push!(clean_cmds, _build_clean_cmd(full_work_dir, sdpb_build_dir, sdpb_docker_image))
    end

    isempty(to_solve_indices) && return results

    n = length(to_solve_indices)

    if !isnothing(slurm_config)
        chunks = Iterators.partition(1:n, batch_size)
        n_jobs = length(chunks)
        for (chunk_idx, chunk) in enumerate(chunks)
            chunk_cmds  = combined_cmds[chunk]
            chunk_clean = clean_cmds[chunk]
            cmds_str    = join(chunk_cmds, "\n")
            post_cmd    = !keep_intermediate ? join(chunk_clean, "\n") : nothing
            job_id      = submit_slurm_job(slurm_config, cmds_str, wait=false,
                                           post_command=post_cmd)
            verbosity >= 2 && @info "Slurm job $chunk_idx/$n_jobs submitted with ID: $job_id ($(length(chunk)) problem(s))"
        end
    else
        verbosity >= 2 && @info "Running $n SDPs directly..."

        for (j, (i, cmd, out_dir)) in enumerate(zip(to_solve_indices, combined_cmds, output_dirs))
            if verbosity == 1
                print("\r  [$(lpad(j - 1, ndigits(n)))/$n] Solving...")
                flush(stdout)
            elseif verbosity >= 2
                @info "Running problem $j/$n..."
            end

            if verbosity <= 1
                with_logger(NullLogger()) do
                    run(pipeline(`sh -c $cmd`, stdout=devnull, stderr=devnull))
                    results[i] = parse_sdpb_output(out_dir, sdps[i], read_z=true)
                end
            else
                run(`sh -c $cmd`)
                results[i] = parse_sdpb_output(out_dir, sdps[i], read_z=true)
            end

            if !keep_intermediate
                run(`sh -c $(clean_cmds[j])`)
            end
        end

        if verbosity == 1
            println("\r  [$n/$n] Done.$(repeat(' ', 30))")
        end
    end

    return results
end

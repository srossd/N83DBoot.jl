"""
Utility helpers for building ABJM bootstrap problems.

- `operator_grid`: build the full operator list (fixed + semishort + long)
- `crossing_constraints`: generate all ABJM crossing constraints up to cutoff Λ
- `integral_constraints`: generate the localization integral constraint
- `twist_spacing`: factory for piecewise-constant twist spacing functions
"""

# ============================================================================
# Twist-spacing helpers
# ============================================================================

"""
    twist_spacing(pairs::Vararg{Tuple{Real,Real}}) -> Function

Factory for a piecewise-constant twist-spacing function.

Each pair `(τ_max, spacing)` says: use `spacing` while `τ < τ_max`.
The last pair's spacing applies for all larger values.

# Example
```julia
fn = twist_spacing((4.0, 0.5), (10.0, 1.0), (20.0, 2.0), (Inf, 4.0))
fn(1.5, 0)  # → 0.5
fn(6.0, 2)  # → 1.0
```
"""
function twist_spacing(pairs::Vararg{Tuple{<:Real,<:Real}})
    return function (tau::Float64, _::Int)
        for i in eachindex(pairs)
            tau < pairs[i][1] && return pairs[i][2]
        end
        return pairs[end][2]
    end
end

"""
    default_twist_spacing

Default twist-spacing function for `operator_grid`. Dense near the unitarity bound,
coarser at large twist.
"""
const default_twist_spacing = twist_spacing(
    (4.0, 0.5),
    (10.0, 1.0),
    (20.0, 2.0),
    (Inf, 4.0),
)

# ============================================================================
# operator_grid
# ============================================================================

"""
    operator_grid(; long_spins, max_twist, twist_spacing_fn,
                    semishort_even_spins, semishort_odd_spins,
                    include_fixed) -> Vector{Operator}

Build the full operator list for an ABJM bootstrap calculation.

The returned vector contains (in canonical sort order):
1. Fixed-OPE operators: `Identity`, `Bp0020`, `Bp0040`, `Btwo0200`
   (only if `include_fixed=true`)
2. Semishort multiplets: `SemishortAp0020(J)` for each even J in
   `semishort_even_spins` and `SemishortAtwo0100(J)` for each odd J in
   `semishort_odd_spins`
3. Long multiplets: `LongOperator(τ, J)` for each J in `long_spins`,
   with τ running from 1 (unitarity bound) up to `max_twist` in
   steps given by `twist_spacing_fn(τ, J)`

# Arguments
- `long_spins`: Spins for the long-operator grid (default: 0:30)
- `max_twist`: Upper cutoff on twist for long operators (default: 60.0)
- `twist_spacing_fn`: `(tau, J) → spacing` function (default: `default_twist_spacing`)
- `semishort_even_spins`: Even spins for `SemishortAp0020` (default: 0:2:30)
- `semishort_odd_spins`: Odd spins for `SemishortAtwo0100` (default: 1:2:29)
- `include_fixed`: Whether to include the four fixed operators (default: true)

# Examples
```julia
ops = operator_grid()
ops = operator_grid(long_spins=collect(0:20), max_twist=40.0)
ops = operator_grid(semishort_even_spins=collect(0:2:20),
                    semishort_odd_spins=collect(1:2:19))
```
"""
function operator_grid(;
    long_spins::Vector{Int}           = collect(0:30),
    max_twist::Float64                = 60.0,
    twist_spacing_fn::Function        = default_twist_spacing,
    semishort_even_spins::Vector{Int} = collect(0:2:30),
    semishort_odd_spins::Vector{Int}  = collect(1:2:29),
    include_fixed::Bool               = true,
)::Vector{Operator}

    operators = Operator[]

    # 1. Fixed operators (exactly one of each; Identity is required for normalization)
    if include_fixed
        push!(operators, Identity)
        push!(operators, Bp0020)
        push!(operators, Bp0040)
        push!(operators, Btwo0200)
    end

    # 2. Semishort multiplets
    for J in semishort_even_spins
        push!(operators, SemishortAp0020(J))
    end
    for J in semishort_odd_spins
        push!(operators, SemishortAtwo0100(J))
    end

    # 3. Long multiplets: τ from unitarity bound (τ=1) up to max_twist
    for J in long_spins
        tau = 1.0
        while tau <= max_twist
            push!(operators, LongOperator(tau, J))
            spacing = twist_spacing_fn(tau, J)
            tau += spacing
        end
    end

    return operators
end

# ============================================================================
# crossing_constraints
# ============================================================================

"""
    crossing_constraints(; Λ::Int, precision=1024, r_order=100) -> Vector{SDPConstraint}

Generate all ABJM crossing constraints up to derivative cutoff Λ.

ABJM has two crossing channels with parity condition m+n **even**:

- **Channel 1**: `CrossingFunctional(m, 0, 1)` for m = 0, 2, 4, …, Λ
  (even m, n=0)
- **Channel 2**: `CrossingFunctional(m, n, 2)` for m ∈ [0,Λ], n ∈ [0,m],
  m+n even (so both m,n even or both odd), m+n ≤ Λ

Each constraint has `rhs=0` and `comparison=Equal`.

# Arguments
- `Λ::Int`: Maximum total derivative order m+n (the bootstrap cutoff)
- `precision::Int`: BigFloat precision in bits (default: 1024)
- `r_order::Int`: Conformal-block expansion order in r (default: 100)

# Examples
```julia
cs = crossing_constraints(Λ=19)
cs = crossing_constraints(Λ=9, precision=512, r_order=50)
```
"""
function crossing_constraints(; Λ::Int, precision::Int=1024, r_order::Int=100)
    constraints = SDPConstraint[]

    # Channel 1: (m, 0), m even, 0 < m ≤ Λ
    for m in 2:2:Λ
        f = CrossingFunctional(m, 0, 1; precision=precision, r_order=r_order)
        push!(constraints, SDPConstraint(f))
    end

    # Channel 2: (m, n), m+n even, 0 ≤ n ≤ m, m+n ≤ Λ
    for m in 0:Λ
        for n in 0:min(m, Λ - m)
            iseven(m + n) || continue
            f = CrossingFunctional(m, n, 2; precision=precision, r_order=r_order)
            push!(constraints, SDPConstraint(f))
        end
    end

    return constraints
end

# ============================================================================
# integral_constraints
# ============================================================================

"""
    integral_constraints(loc::LocalizationData;
                         precision=64, r_order=20) -> SDPConstraint

Generate the localization integral constraint for a given ABJM theory point.

Wraps an `IntegralFunctional` (b=0) in an `SDPConstraint` whose right-hand side
is taken from `loc.rhs` (the nosaka.m matrix-model result for this (N, k, M) point).

# Arguments
- `loc::LocalizationData`: Localization data for the theory
- `precision::Int`: BigFloat precision in bits passed to `IntegralFunctional` (default: 64)
- `r_order::Int`: Block-expansion order passed to `IntegralFunctional` (default: 20)

# Examples
```julia
loc = load_localization_data(4, 1)
c   = integral_constraints(loc)
```
"""
function integral_constraints(loc::LocalizationData;
                              precision::Int=64,
                              r_order::Int=20)
    f = IntegralFunctional(0.0; precision=precision, r_order=r_order)
    return [SDPConstraint(f; rhs=loc.rhs, comparison=Equal)]
end

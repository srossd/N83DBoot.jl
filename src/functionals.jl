"""
Functional types for 3D N=8 SCFT (ABJM theory) conformal bootstrap calculations.

Two functional types cover the ABJM bootstrap:

- `CrossingFunctional(m, n, cross)`: derivative of the crossing equation at the
  crossing-symmetric point; two channels (cross ∈ {1,2}).
- `IntegralFunctional(b)`: integral of the crossing equation weighted by the
  ABJM localization measure; parameter `b` controls the integration boundary.
"""

"""
    Functional

Abstract base type for all functionals in the 3D N=8 SCFT bootstrap.
"""
abstract type Functional end

# ============================================================================
# CrossingFunctional
# ============================================================================

"""
    CrossingFunctional(m::Int, n::Int, cross::Int; precision::Int=1024, r_order::Int=100)

Crossing equation functional specified by derivative orders (m, n) and channel index.

ABJM has two crossing channels (from `abjm_defs.m`):
- Channel 1: (m, 0) with m even, m = 0, 2, …, Λ
- Channel 2: (m, n) with m+n even, m ∈ [0,Λ], n ∈ [0, min(m, Λ-m)]

Note: the parity condition is **m+n even** (opposite of 4D N=4 where m+n is odd).

# Arguments
- `m::Int`: z-derivative order, m ≥ 0
- `n::Int`: z̄-derivative order, n ≥ 0
- `cross::Int`: Crossing channel index, must be 1 or 2
- `precision::Int`: Precision in bits for extended arithmetic (default: 1024)
- `r_order::Int`: Order of the block expansion in r (default: 100)

# Constraints
- `cross ∈ {1, 2}`
- `m + n` must be even
- For channel 1: `n == 0` and `m` must be even

# Examples
```julia
f = CrossingFunctional(4, 0, 1)   # Channel 1, (4,0) derivative
f = CrossingFunctional(2, 2, 2)   # Channel 2, (2,2) derivative
f = CrossingFunctional(2, 2, 2; precision=512, r_order=50)
```
"""
struct CrossingFunctional <: Functional
    m::Int
    n::Int
    cross::Int
    precision::Int
    r_order::Int

    function CrossingFunctional(m::Int, n::Int, cross::Int;
                                precision::Int=1024, r_order::Int=100)
        @assert cross ∈ (1, 2) "cross must be 1 or 2, got cross=$cross"
        @assert m >= 0 && n >= 0 "m, n must be non-negative, got m=$m, n=$n"
        @assert iseven(m + n) "m + n must be even for ABJM crossing, got m+n=$(m+n)"
        if cross == 1
            @assert n == 0 "Channel 1 only uses n=0 derivatives, got n=$n"
            @assert iseven(m) "Channel 1 requires even m, got m=$m"
        end
        @assert precision > 0 "Precision must be positive, got precision=$precision"
        @assert r_order > 0 "r_order must be positive, got r_order=$r_order"
        new(m, n, cross, precision, r_order)
    end
end

Base.show(io::IO, f::CrossingFunctional) =
    print(io, "CrossingFunctional(m=$(f.m), n=$(f.n), cross=$(f.cross), prec=$(f.precision), r_order=$(f.r_order))")

# ============================================================================
# IntegralFunctional
# ============================================================================

"""
    IntegralFunctional(b::Float64=0.0; precision::Int=256, r_order::Int=20)

Integral functional that integrates the ABJM crossing equation against the
localization measure.

The Mathematica function `integral[superblock, b, prec]` in `abjm_defs.m`
combines the three prefactors (1/u, 1/u², v/u²) into a single numerical
value; this functional wraps that combined result. Different values of `b`
parametrize different integration domains.

The right-hand side of the integral constraint for a given theory is stored
in `LocalizationData.rhs` and comes from `nosaka.m`.

# Arguments
- `b::Float64`: Integration boundary parameter, b ≥ 0 (default: 0.0)
- `precision::Int`: Precision in bits (default: 256)
- `r_order::Int`: Order of the block expansion in r (default: 20)

# Examples
```julia
f = IntegralFunctional()              # b=0 (default)
f = IntegralFunctional(1.5)           # b=1.5
f = IntegralFunctional(0.0; precision=512, r_order=40)
```
"""
struct IntegralFunctional <: Functional
    b::Float64
    precision::Int
    r_order::Int

    function IntegralFunctional(b::Float64=0.0;
                                precision::Int=64, r_order::Int=20)
        @assert b >= 0.0 "b must be >= 0, got b=$b"
        @assert precision > 0 "Precision must be positive, got precision=$precision"
        @assert r_order > 0 "r_order must be positive, got r_order=$r_order"
        new(b, precision, r_order)
    end
end

Base.show(io::IO, f::IntegralFunctional) =
    print(io, "IntegralFunctional(b=$(f.b), prec=$(f.precision), r_order=$(f.r_order))")

# ============================================================================
# LinearCombinationObjective
# ============================================================================

"""
    LinearCombinationObjective(op1 => c1, op2 => c2, ...)

Objective that bounds a linear combination of free-operator OPE² coefficients.

Use this as the `objective` argument to `SDPProblem` when you want SDPB to
bound `c₁ λ²_{A} + c₂ λ²_{B}` rather than a single functional value.
By sweeping the coefficient vector `(c₁, c₂)` over many directions, you can
trace out the full boundary of the joint allowed region in `(λ²_A, λ²_B)` space.

Fixed operators (`IdentityOperator`, `Bp0020Operator`, etc.) are not included in
the combination — their contribution to the dual bound is zero.  No Mathematica
computation or cache lookup is required; values come directly from the coefficient
dictionary.

# Examples
```julia
# Bound cos(θ)·λ²_{Ap(0)} + sin(θ)·λ²_{A2(1)}
obj = LinearCombinationObjective(
    SemishortAp0020(0)   => cos(θ),
    SemishortAtwo0100(1) => sin(θ),
)
sdp = SDPProblem(obj, constraints, operators, loc; direction=1)
```
"""
struct LinearCombinationObjective <: Functional
    operator_coefficients::Dict{FreeOperator, Float64}

    function LinearCombinationObjective(pairs::Pair{<:FreeOperator, <:Real}...)
        new(Dict{FreeOperator, Float64}(op => Float64(c) for (op, c) in pairs))
    end
end

function Base.show(io::IO, f::LinearCombinationObjective)
    parts = sort(["$(round(c; digits=8)) * $(op)" for (op, c) in f.operator_coefficients])
    print(io, "LinearCombinationObjective(", join(parts, " + "), ")")
end

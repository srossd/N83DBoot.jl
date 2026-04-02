"""
Operator types for 3D N=8 SCFT (ABJM theory) conformal bootstrap calculations.

This module defines the type hierarchy for operators used in the bootstrap.

## Abstract types
- `Operator`: base type for all operators
- `FixedOperator <: Operator`: operators whose OPE coefficients are fixed
  (either = 1 by normalization, or by localization data)
- `FreeOperator <: Operator`: operators with free (≥ 0) OPE coefficients;
  these generate PSD constraints in the dual SDP

## Concrete types
### Fixed-OPE operators (singletons)
- `IdentityOperator`: Identity operator; OPE coeff = 1 by normalization convention
- `Bp0020Operator`: [0,2,0] short multiplet; OPE coeff² from LocalizationData
- `Bp0040Operator`: [0,4,0] short multiplet; OPE coeff² from LocalizationData
- `Btwo0200Operator`: twist-2 short multiplet; OPE coeff² from LocalizationData

Convenience singleton constants `Identity`, `Bp0020`, `Bp0040`, `Btwo0200` are exported.

### Semishort multiplets (FreeOperator)
- `SemishortAp0020(J)`: even-spin semishort; dimension pinned to J+1
- `SemishortAtwo0100(J)`: odd-spin semishort; dimension pinned to J+1

### Long multiplets (FreeOperator)
- `LongOperator(tau, J)`: generic operators, twist τ = Δ−J ≥ 1
"""

# ============================================================================
# Abstract type tree
# ============================================================================

"""
    Operator

Abstract base type for all operators in the 3D N=8 SCFT bootstrap.
"""
abstract type Operator end

"""
    FixedOperator <: Operator

Abstract type for operators whose OPE coefficients are fixed (not free SDP variables).
These operators contribute to the normalization vector in the dual problem:
- `IdentityOperator`: normalization anchor (OPE coeff = 1)
- `Bp0020Operator`, `Bp0040Operator`, `Btwo0200Operator`: localization-fixed OPE coefficients
"""
abstract type FixedOperator <: Operator end

"""
    FreeOperator <: Operator

Abstract type for operators with free (≥ 0) OPE coefficients.
These operators generate positive semidefinite constraints in the dual SDP.
"""
abstract type FreeOperator <: Operator end

# ============================================================================
# Fixed-OPE operators (singletons)
# ============================================================================

"""
    IdentityOperator()

The identity operator. Its OPE coefficient is fixed to 1 by the normalization
convention. Plays the role that `ShortOperator` plays in N4SYMBoot — exactly
one must appear in the `operators` vector of an `SDPProblem`.

# Examples
```julia
op = IdentityOperator()
op = Identity   # convenience singleton
```
"""
struct IdentityOperator <: FixedOperator end

"""
    Bp0020Operator()

The [0,2,0] short (BPS) multiplet. Its OPE coefficient squared is fixed by
the localization data in `LocalizationData.lambda2_Bp0020`.

# Examples
```julia
op = Bp0020Operator()
op = Bp0020   # convenience singleton
```
"""
struct Bp0020Operator <: FixedOperator end

"""
    Bp0040Operator()

The [0,4,0] short (BPS) multiplet. Its OPE coefficient squared is fixed by
the localization data in `LocalizationData.lambda2_Bp0040`.

# Examples
```julia
op = Bp0040Operator()
op = Bp0040   # convenience singleton
```
"""
struct Bp0040Operator <: FixedOperator end

"""
    Btwo0200Operator()

The twist-2 short multiplet. Its OPE coefficient squared is fixed by
the localization data in `LocalizationData.lambda2_Btwo0200`.

# Examples
```julia
op = Btwo0200Operator()
op = Btwo0200   # convenience singleton
```
"""
struct Btwo0200Operator <: FixedOperator end

"""Pre-constructed singleton instances for the four fixed operators."""
const Identity = IdentityOperator()
const Bp0020   = Bp0020Operator()
const Bp0040   = Bp0040Operator()
const Btwo0200 = Btwo0200Operator()

# ============================================================================
# Semishort operators
# ============================================================================

"""
    SemishortAp0020(J::Int)

Even-spin semishort multiplet from the Ap0020 family. The conformal dimension
is pinned to Δ = J + 1 (unitarity-saturating shortening condition), so only
the spin J is a free parameter.

# Arguments
- `J::Int`: Even spin, J ≥ 0

# Examples
```julia
op = SemishortAp0020(0)   # J=0 semishort
op = SemishortAp0020(4)   # J=4 semishort
```
"""
struct SemishortAp0020 <: FreeOperator
    J::Int

    function SemishortAp0020(J::Int)
        @assert iseven(J) && J >= 0 "Ap0020 requires even non-negative spin, got J=$J"
        new(J)
    end
end

"""
    SemishortAtwo0100(J::Int)

Odd-spin semishort multiplet from the Atwo0100 family. The conformal dimension
is pinned to Δ = J + 1 (unitarity-saturating shortening condition), so only
the spin J is a free parameter.

# Arguments
- `J::Int`: Odd spin, J ≥ 1

# Examples
```julia
op = SemishortAtwo0100(1)   # J=1 semishort
op = SemishortAtwo0100(3)   # J=3 semishort
```
"""
struct SemishortAtwo0100 <: FreeOperator
    J::Int

    function SemishortAtwo0100(J::Int)
        @assert isodd(J) && J >= 1 "Atwo0100 requires odd spin >= 1, got J=$J"
        new(J)
    end
end

# ============================================================================
# Long operator
# ============================================================================

"""
    LongOperator(tau::Float64, J::Int)

Generic (long) multiplet parametrized by twist τ = Δ − J and spin J.
The OPE coefficient is a free (≥ 0) variable in the bootstrap SDP.

# Arguments
- `tau::Float64`: Twist τ = Δ − J ≥ 1 (unitarity bound)
- `J::Int`: Spin J ≥ 0

# Examples
```julia
op = LongOperator(1.5, 0)   # τ=1.5, spin-0 long operator
op = LongOperator(2.0, 2)   # τ=2.0, spin-2 long operator (Δ=4.0)
```
"""
struct LongOperator <: FreeOperator
    tau::Float64
    J::Int

    function LongOperator(tau::Float64, J::Int)
        @assert J >= 0 "Spin must be non-negative, got J=$J"
        @assert tau >= 1.0 - 1e-10 "Must satisfy unitarity bound τ >= 1, got τ=$tau"
        new(tau, J)
    end
end

# ============================================================================
# Display
# ============================================================================

Base.show(io::IO, ::IdentityOperator)   = print(io, "IdentityOperator()")
Base.show(io::IO, ::Bp0020Operator)     = print(io, "Bp0020Operator()")
Base.show(io::IO, ::Bp0040Operator)     = print(io, "Bp0040Operator()")
Base.show(io::IO, ::Btwo0200Operator)   = print(io, "Btwo0200Operator()")
Base.show(io::IO, op::SemishortAp0020)    = print(io, "SemishortAp0020(J=$(op.J))")
Base.show(io::IO, op::SemishortAtwo0100)  = print(io, "SemishortAtwo0100(J=$(op.J))")
Base.show(io::IO, op::LongOperator)     = print(io, "LongOperator(τ=$(op.tau), J=$(op.J))")

# ============================================================================
# Sorting
# ============================================================================

# Fixed operators sort before free operators; within fixed they sort by "weight"
_fixed_order(::IdentityOperator)  = 0
_fixed_order(::Bp0020Operator)    = 1
_fixed_order(::Bp0040Operator)    = 2
_fixed_order(::Btwo0200Operator)  = 3

Base.isless(a::FixedOperator, b::FixedOperator)   = _fixed_order(a) < _fixed_order(b)
Base.isless(::FixedOperator, ::FreeOperator)       = true
Base.isless(::FreeOperator, ::FixedOperator)       = false

# Semishort operators sort before long operators, ordered by J
Base.isless(a::SemishortAp0020, b::SemishortAp0020)       = a.J < b.J
Base.isless(a::SemishortAtwo0100, b::SemishortAtwo0100)   = a.J < b.J
Base.isless(::SemishortAp0020, ::SemishortAtwo0100)       = true
Base.isless(::SemishortAtwo0100, ::SemishortAp0020)       = false
Base.isless(::SemishortAp0020,   ::LongOperator)          = true
Base.isless(::LongOperator,      ::SemishortAp0020)       = false
Base.isless(::SemishortAtwo0100, ::LongOperator)          = true
Base.isless(::LongOperator,      ::SemishortAtwo0100)     = false

# Long operators sort by J then tau
Base.isless(a::LongOperator, b::LongOperator) =
    a.J < b.J || (a.J == b.J && a.tau < b.tau)

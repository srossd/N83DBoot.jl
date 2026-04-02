using Test
using N83DBoot

# ============================================================================
# operators.jl
# ============================================================================

@testset "Operators" begin
    @testset "Fixed operators" begin
        @test Identity isa IdentityOperator
        @test Bp0020   isa Bp0020Operator
        @test Bp0040   isa Bp0040Operator
        @test Btwo0200 isa Btwo0200Operator

        @test Identity   isa FixedOperator
        @test Bp0020     isa FixedOperator
        @test Bp0040     isa FixedOperator
        @test Btwo0200   isa FixedOperator
    end

    @testset "SemishortAp0020 validation" begin
        @test SemishortAp0020(0) isa FreeOperator
        @test SemishortAp0020(4) isa FreeOperator
        @test_throws AssertionError SemishortAp0020(-1)   # negative spin
        @test_throws AssertionError SemishortAp0020(1)    # odd spin
        @test_throws AssertionError SemishortAp0020(3)    # odd spin
    end

    @testset "SemishortAtwo0100 validation" begin
        @test SemishortAtwo0100(1) isa FreeOperator
        @test SemishortAtwo0100(3) isa FreeOperator
        @test_throws AssertionError SemishortAtwo0100(0)   # even spin
        @test_throws AssertionError SemishortAtwo0100(2)   # even spin
        @test_throws AssertionError SemishortAtwo0100(-1)  # negative
    end

    @testset "LongOperator validation" begin
        @test LongOperator(1.0, 0) isa FreeOperator
        @test LongOperator(3.0, 2) isa FreeOperator
        @test_throws AssertionError LongOperator(0.5, 0)   # below unitarity bound
        @test_throws AssertionError LongOperator(2.0, 2)   # Δ < J+1
        @test_throws AssertionError LongOperator(1.0, -1)  # negative spin

        op = LongOperator(5.0, 2)
        @test op.delta == 5.0
        @test op.J     == 2
        @test op.spin  == 2
        @test op.twist ≈ 3.0
    end

    @testset "Operator sorting" begin
        ops = [LongOperator(3.0, 0), SemishortAtwo0100(1), Identity,
               SemishortAp0020(0), Bp0020, LongOperator(2.0, 1)]
        sorted = sort(ops)
        @test sorted[1] isa IdentityOperator
        @test sorted[2] isa Bp0020Operator
        # Fixed operators before free
        n_fixed = count(op -> op isa FixedOperator, sorted)
        @test all(sorted[1:n_fixed] isa FixedOperator for i in 1:n_fixed)
    end
end

# ============================================================================
# functionals.jl
# ============================================================================

@testset "Functionals" begin
    @testset "CrossingFunctional validation" begin
        @test CrossingFunctional(0, 0, 1) isa CrossingFunctional
        @test CrossingFunctional(4, 0, 1) isa CrossingFunctional
        @test CrossingFunctional(2, 2, 2) isa CrossingFunctional
        @test CrossingFunctional(0, 0, 2) isa CrossingFunctional

        # cross must be 1 or 2
        @test_throws AssertionError CrossingFunctional(2, 0, 3)
        # m+n must be even
        @test_throws AssertionError CrossingFunctional(1, 0, 2)  # 1+0 = 1 (odd)
        @test_throws AssertionError CrossingFunctional(3, 0, 1)  # 3+0 = 3 (odd)
        # Channel 1: n must be 0
        @test_throws AssertionError CrossingFunctional(2, 2, 1)
        # Channel 1: m must be even
        @test_throws AssertionError CrossingFunctional(3, 0, 1)  # odd m (also fails m+n even)
    end

    @testset "IntegralFunctional validation" begin
        @test IntegralFunctional(0.0) isa IntegralFunctional
        @test IntegralFunctional(1.5) isa IntegralFunctional
        @test_throws AssertionError IntegralFunctional(-0.1)   # b < 0
    end
end

# ============================================================================
# config.jl
# ============================================================================

@testset "Config" begin
    original_cache = get_config().cache_dir

    set_config!(cache_dir="/tmp/test_cache")
    @test get_config().cache_dir == "/tmp/test_cache"

    set_config!(cache_dir=original_cache)
    @test get_config().cache_dir == original_cache
end

# ============================================================================
# util.jl
# ============================================================================

@testset "Utilities" begin
    @testset "operator_grid" begin
        ops = operator_grid(
            long_spins           = [0, 2],
            max_delta            = 5.0,
            semishort_even_spins = [0, 2],
            semishort_odd_spins  = [1],
        )
        @test any(op -> op isa IdentityOperator,     ops)
        @test any(op -> op isa Bp0020Operator,        ops)
        @test any(op -> op isa Bp0040Operator,        ops)
        @test any(op -> op isa Btwo0200Operator,      ops)
        @test any(op -> op isa SemishortAp0020 && op.J == 0, ops)
        @test any(op -> op isa SemishortAp0020 && op.J == 2, ops)
        @test any(op -> op isa SemishortAtwo0100 && op.J == 1, ops)
        @test any(op -> op isa LongOperator, ops)

        # All long operators within bounds
        long_ops = filter(op -> op isa LongOperator, ops)
        @test all(op -> op.delta <= 5.0, long_ops)
        @test all(op -> op.delta >= op.J + 1.0 - 1e-10, long_ops)
    end

    @testset "operator_grid no fixed" begin
        ops = operator_grid(
            long_spins           = [0],
            max_delta            = 3.0,
            semishort_even_spins = Int[],
            semishort_odd_spins  = Int[],
            include_fixed        = false,
        )
        @test !any(op -> op isa FixedOperator, ops)
    end

    @testset "crossing_constraints" begin
        cs = crossing_constraints(Λ=4)
        @test !isempty(cs)
        @test all(c -> c.functional isa CrossingFunctional, cs)

        # All constraints should have rhs=0 and comparison=Equal
        @test all(c -> c.rhs == 0, cs)
        @test all(c -> c.comparison == Equal, cs)

        # ABJM parity: m+n must be even for all generated constraints
        for c in cs
            f = c.functional
            @test iseven(f.m + f.n)
        end

        # Channel 1: n=0 and m even
        ch1 = filter(c -> c.functional.cross == 1, cs)
        @test all(c -> c.functional.n == 0, ch1)
        @test all(c -> iseven(c.functional.m), ch1)

        # Channel 2: all valid (m,n) with m+n even and m+n ≤ Λ
        ch2 = filter(c -> c.functional.cross == 2, cs)
        @test !isempty(ch2)
    end

    @testset "delta_spacing" begin
        fn = delta_spacing((4.0, 0.5), (10.0, 1.0), (Inf, 2.0))
        @test fn(2.0, 0) ≈ 0.5   # twist = 2.0 < 4.0
        @test fn(6.0, 0) ≈ 1.0   # twist = 6.0 < 10.0
        @test fn(15.0, 0) ≈ 2.0  # twist = 15.0 ≥ 10.0
    end
end

# ============================================================================
# Cache path generation
# ============================================================================

@testset "Cache paths" begin
    f_cross    = CrossingFunctional(4, 0, 1; precision=1024, r_order=100)
    f_integral = IntegralFunctional(0.0; precision=256, r_order=20)

    @test N83DBoot.functional_type_name(f_cross)    == "Crossing"
    @test N83DBoot.functional_type_name(f_integral) == "Integral"

    name_cross = N83DBoot.get_cache_filename(f_cross)
    @test occursin("cross_1", name_cross)
    @test occursin("m_4",     name_cross)
    @test occursin("n_0",     name_cross)

    name_int = N83DBoot.get_cache_filename(f_integral)
    @test occursin("b_0.0",   name_int)

    # Operator-specific cache filenames
    @test N83DBoot.operator_cache_filename(Identity)            == "identity.csv"
    @test N83DBoot.operator_cache_filename(Bp0020)              == "Bp0020.csv"
    @test N83DBoot.operator_cache_filename(Bp0040)              == "Bp0040.csv"
    @test N83DBoot.operator_cache_filename(Btwo0200)            == "Btwo0200.csv"
    @test N83DBoot.operator_cache_filename(SemishortAp0020(4))   == "semishort_Ap0020_J_4.csv"
    @test N83DBoot.operator_cache_filename(SemishortAtwo0100(3)) == "semishort_Atwo0100_J_3.csv"
    @test N83DBoot.operator_cache_filename(LongOperator(3.0, 0)) == "spin_0.csv"
    @test N83DBoot.operator_cache_filename(LongOperator(4.5, 2)) == "spin_2.csv"
end

# ============================================================================
# SDPProblem construction (without cache — just structure)
# ============================================================================

@testset "SDPProblem structure" begin
    # Construct a minimal LocalizationData for testing.
    # Conventions: lambda2_Bp0020 = 256/cT(Mathematica)
    #              cT = 256 / lambda2_Bp0020
    lambda2_Bp0020   = BigFloat("3.35099165327933251")  # k=1, N=4
    lambda2_Btwo0200 = BigFloat("7.45176137921519754")
    lambda2_Bp0040   = (4 * lambda2_Bp0020 + lambda2_Btwo0200 + 16) / 5
    cT_val           = BigFloat(256) / lambda2_Bp0020
    loc = LocalizationData(4, 1, 0,
        cT_val,
        lambda2_Bp0020,
        lambda2_Bp0040,
        lambda2_Btwo0200,
        BigFloat("0.00649138483"),   # rhs (approximate)
    )

    objective   = CrossingFunctional(0, 0, 1)
    constraints = [SDPConstraint(CrossingFunctional(2, 0, 1))]
    operators   = [Identity, Bp0020, Bp0040, Btwo0200, SemishortAp0020(0), LongOperator(2.0, 0)]

    sdp = SDPProblem(objective, constraints, operators, loc)
    @test sdp.localization.N == 4
    @test sdp.localization.k == 1
    @test sdp.precision      == 1024
    @test sdp.direction      == 1

    N, k, M = extract_NkM(sdp)
    @test N == 4 && k == 1 && M == 0

    # Must have exactly one IdentityOperator
    @test_throws AssertionError SDPProblem(objective, constraints,
        [Bp0020, SemishortAp0020(0), LongOperator(2.0, 0)], loc)

    @test_throws AssertionError SDPProblem(objective, constraints,
        [Identity, Identity, LongOperator(2.0, 0)], loc)
end

# ============================================================================
# extract_Lambda
# ============================================================================

@testset "extract_Lambda" begin
    lam2 = BigFloat("3.35099165327933251")
    loc  = LocalizationData(4, 1, 0,
        BigFloat(256) / lam2, lam2,
        (4*lam2 + BigFloat("7.45176137921519754") + 16) / 5,
        BigFloat("7.45176137921519754"),
        BigFloat("0.00649138483"))
    operators   = [Identity, LongOperator(2.0, 0)]

    # Λ=4: channel 1 = {(0,0),(2,0),(4,0)}, channel 2 = {(0,0),(2,0),(4,0),(2,2),(0,0),...}
    cs = crossing_constraints(Λ=4)
    sdp = SDPProblem(cs[1].functional, cs[2:end], operators, loc)
    Λ = extract_Lambda(sdp)
    @test Λ == 4
end

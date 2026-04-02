"""
Localization-based spectral data for ABJM theory (3D N=8 SCFTs).

Loads nosaka.m localization outputs — OPE coefficients of protected operators
and the integral constraint right-hand side — for a given ABJM theory point
(N, k, M).

## File layout under `localization_data_dir`

Pre-computed CSV tables produced by `short_coeffs.m` and `rhs.m`:
```
short_coeffs_k1.csv    # k=1, M=0: columns [lambda2_Bp0020, lambda2_Bp0040, lambda2_Btwo0200], row i → N=i+1
short_coeffs_k2.csv    # k=2, M=0: same columns
short_coeffs_k2_M1.csv # k=2, M=1: same columns
rhs_k1.csv             # k=1, M=0: columns [rhs], row i → N=i+1
rhs_k2.csv             # k=2, M=0: columns [rhs], row i → N=i+1
rhs_k2_M1.csv          # k=2, M=1: columns [rhs], row i → N=i+1
```
Rows are 1-indexed with the first row corresponding to N=2.
OPE coefficient conventions (matching short_ope.m / InputLoader.cpp):
- `lambda2_Bp0020` = 256/cT(Mathematica) — stress-tensor OPE coefficient
- `lambda2_Bp0040` = λBp = (4*lambda2_Bp0020 + lambda2_Btwo0200 + 16)/5 — SUSY Ward identity
- `lambda2_Btwo0200` = λB2 — twist-2 OPE coefficient from localization

For N beyond the tabulated range, `load_localization_data` can call `nosaka.m`
via a Mathematica subprocess (the same pattern as `compute_functional` in `cache.jl`).
"""

using CSV
using DataFrames
using Downloads
using SHA

# ---------------------------------------------------------------------------
# Configuration helpers
# ---------------------------------------------------------------------------

"""
    get_localization_dir() -> String

Return the active localization data directory (`config.localization_data_dir`).
"""
function get_localization_dir()
    return get_config().localization_data_dir
end

# ---------------------------------------------------------------------------
# RAM cache:  (data_dir, k, M) → Dict(filename → DataFrame)
# ---------------------------------------------------------------------------

const _LOCALIZATION_RAM_CACHE = Dict{Tuple{String,Int,Int}, Dict{String,DataFrame}}()

# ---------------------------------------------------------------------------
# File naming helpers
# ---------------------------------------------------------------------------

function _localization_filename(k::Int, M::Int)
    if k == 1 && M == 0
        return ("short_coeffs_k1.csv", "rhs_k1.csv")
    elseif k == 2 && M == 0
        return ("short_coeffs_k2.csv", "rhs_k2.csv")
    elseif k == 2 && M == 1
        return ("short_coeffs_k2_M1.csv", "rhs_k2_M1.csv")
    else
        error("No pre-computed localization data for k=$k, M=$M. " *
              "Supported: (k=1,M=0), (k=2,M=0), (k=2,M=1).")
    end
end

# ---------------------------------------------------------------------------
# LocalizationData type
# ---------------------------------------------------------------------------

"""
    LocalizationData

Localization (matrix model) inputs for a specific ABJM theory point (N, k, M).

These values are computed by `nosaka.m` and describe the exact OPE coefficients
of the protected (short) multiplets appearing in the OPE of the stress-tensor
multiplet with itself.

# Fields
- `N::Int`: Rank (number of colors)
- `k::Int`: Chern-Simons level (1 or 2 for tabulated data)
- `M::Int`: Monopole sector (0 for vanilla ABJM, 1 for the M=1 sector at k=2)
- `cT::BigFloat`: Weyl anomaly coefficient (derived as 256/λ²_{Bp0020})
- `lambda2_Bp0020::BigFloat`: λ²_{Bp0020} — 256/cT(Mathematica), = `ope[0]` from short_ope.m
- `lambda2_Bp0040::BigFloat`: λ²_{Bp0040} — λBp = (4*λ²_Bp0020 + λ²_Btwo0200 + 16)/5, = `ope[1]`
- `lambda2_Btwo0200::BigFloat`: λ²_{Btwo0200} — λB2 from nosaka.m, = `ope[2]`
- `rhs::BigFloat`: Right-hand side of the integral constraint, from nosaka.m

Note on conventions: `lambda2_Bp0020` = 256/cT (in the Mathematica normalization), which is also
the stress-tensor multiplet OPE coefficient. `lambda2_Bp0040` is determined by the SUSY Ward
identity. Only `lambda2_Btwo0200` and `lambda2_Bp0020` are independently computed from localization.
"""
struct LocalizationData
    N::Int
    k::Int
    M::Int
    cT::BigFloat
    lambda2_Bp0020::BigFloat
    lambda2_Bp0040::BigFloat
    lambda2_Btwo0200::BigFloat
    rhs::BigFloat
end

function Base.show(io::IO, loc::LocalizationData)
    print(io, "LocalizationData(N=$(loc.N), k=$(loc.k), M=$(loc.M), " *
              "cT=$(Float64(loc.cT)), λ²_Bp0020=$(Float64(loc.lambda2_Bp0020)), " *
              "λ²_Bp0040=$(Float64(loc.lambda2_Bp0040)), λ²_Btwo0200=$(Float64(loc.lambda2_Btwo0200)), " *
              "rhs=$(Float64(loc.rhs)))")
end

# ---------------------------------------------------------------------------
# Internal CSV loading
# ---------------------------------------------------------------------------

"""
    _load_localization_csv!(data_dir, k, M) -> Dict{String, DataFrame}

Read the localization CSV files for (k, M) into a RAM-cached Dict.
"""
function _load_localization_csv!(data_dir::String, k::Int, M::Int)::Dict{String,DataFrame}
    cache_key = (data_dir, k, M)
    if haskey(_LOCALIZATION_RAM_CACHE, cache_key)
        return _LOCALIZATION_RAM_CACHE[cache_key]
    end

    coeffs_file, rhs_file = _localization_filename(k, M)
    result = Dict{String,DataFrame}()

    for fname in (coeffs_file, rhs_file)
        fpath = joinpath(data_dir, fname)
        if isfile(fpath)
            try
                df = CSV.read(fpath, DataFrame, header=false)
                result[fname] = df
            catch e
                @warn "Failed to load $fpath: $e"
            end
        else
            @warn "Localization data file not found: $fpath"
        end
    end

    if !isempty(result)
        _LOCALIZATION_RAM_CACHE[cache_key] = result
    end
    return result
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    load_localization_data(N::Int, k::Int, M::Int=0;
                           localization_dir=get_localization_dir(),
                           math_command="math",
                           compute_if_missing=true) -> LocalizationData

Load localization data for ABJM theory with rank N, Chern-Simons level k,
and monopole sector M.

First attempts to read from pre-computed CSV tables in `localization_dir`.
If N is beyond the tabulated range and `compute_if_missing=true`, calls
`nosaka.m` via a Mathematica subprocess to compute the exact values.

# Arguments
- `N::Int`: Rank (number of colors)
- `k::Int`: Chern-Simons level
- `M::Int`: Monopole sector (default: 0)
- `localization_dir::String`: Directory with pre-computed CSV files
- `math_command::String`: Mathematica executable name (default: `"math"`)
- `compute_if_missing::Bool`: Call nosaka.m if N is not in the CSV (default: true)

# Returns
`LocalizationData` with cT, OPE coefficients, and integral RHS.

# Examples
```julia
loc = load_localization_data(4, 1)      # k=1 ABJM at rank N=4
loc = load_localization_data(6, 2, 0)   # k=2 ABJM at rank N=6
```
"""
function load_localization_data(N::Int, k::Int, M::Int=0;
                                localization_dir::String=get_localization_dir(),
                                math_command::String="math",
                                compute_if_missing::Bool=true)::LocalizationData

    if !isdir(localization_dir) || isempty(localization_dir)
        if compute_if_missing
            return _compute_localization_data(N, k, M; math_command=math_command)
        end
        error("Localization data directory not configured. Call set_config!(localization_data_dir=...).")
    end

    data = _load_localization_csv!(localization_dir, k, M)

    coeffs_file, rhs_file = _localization_filename(k, M)

    if haskey(data, coeffs_file) && haskey(data, rhs_file)
        coeffs_df = data[coeffs_file]
        rhs_df    = data[rhs_file]

        # Rows are 1-indexed; the CSV starts at N=2, so row_index = N - 1.
        row_index = N - 1
        if row_index >= 1 && row_index <= nrow(coeffs_df) && row_index <= nrow(rhs_df)
            row_c = coeffs_df[row_index, :]
            row_r = rhs_df[row_index, :]

            # CSV columns: [lambda2_Bp0020, lambda2_Bp0040, lambda2_Btwo0200]
            # (matches short_ope.m output: {256/cT, λBp, λB2} = {ope[0], ope[1], ope[2]})
            # cT (Weyl anomaly) is derived as 256 / lambda2_Bp0020.
            lambda2_Bp0020   = BigFloat(row_c[1])
            lambda2_Bp0040   = BigFloat(row_c[2])
            lambda2_Btwo0200 = BigFloat(row_c[3])
            cT_val = BigFloat(256) / lambda2_Bp0020

            return LocalizationData(
                N,
                k,
                M,
                cT_val,
                lambda2_Bp0020,
                lambda2_Bp0040,
                lambda2_Btwo0200,
                BigFloat(row_r[1]),   # rhs
            )
        else
            if compute_if_missing
                @info "N=$N is beyond the tabulated range (N=2..$(nrow(coeffs_df)+1)). Computing via nosaka.m..."
                return _compute_localization_data(N, k, M; math_command=math_command)
            end
            error("N=$N is beyond the tabulated range (N=2..$(nrow(coeffs_df)+1)) and compute_if_missing=false.")
        end
    else
        if compute_if_missing
            @info "Localization CSV not found; computing via nosaka.m..."
            return _compute_localization_data(N, k, M; math_command=math_command)
        end
        error("Localization CSV files not found in $localization_dir for k=$k, M=$M.")
    end
end

"""
    _compute_localization_data(N, k, M; math_command="math") -> LocalizationData

Compute localization data by calling `nosaka.m` via a Mathematica subprocess.
Returns a `LocalizationData` with the resulting OPE coefficients and RHS.
"""
function _compute_localization_data(N::Int, k::Int, M::Int;
                                    math_command::String="math")::LocalizationData
    script_path = joinpath(pkgdir(@__MODULE__), "scripts", "nosaka_query.m")
    isfile(script_path) || error("nosaka_query.m script not found at $script_path")

    cmd = "$math_command -script $script_path $N $k $M"
    @info "Computing localization data via Mathematica..." cmd

    output = read(`sh -c $cmd`, String)
    lines  = filter(!isempty, strip.(split(output, '\n')))

    # Expected output: 5 whitespace-separated numbers on one line:
    #   cT  lambda2_Bp0020  lambda2_Bp0040  lambda2_Btwo0200  rhs
    data_line = last(filter(l -> occursin(r"^[0-9e+\-. ]+$", l), lines))
    vals = map(x -> parse(BigFloat, x), split(strip(data_line)))
    length(vals) == 5 || error("Unexpected output from nosaka_query.m: $output")

    return LocalizationData(N, k, M, vals[1], vals[2], vals[3], vals[4], vals[5])
end

# ---------------------------------------------------------------------------
# Bundle packaging and distribution
# ---------------------------------------------------------------------------

"""
    package_localization_data(; output_path="localization.tar.gz") -> String

Create a `.tar.gz` archive of all CSV files in `get_localization_dir()` and
save it to `output_path`. The SHA256 checksum is printed after writing.

Returns the absolute path to the created archive.
"""
function package_localization_data(; output_path::String="localization.tar.gz")
    data_dir   = get_localization_dir()
    output_abs = isabspath(output_path) ? output_path : joinpath(pwd(), output_path)

    isdir(data_dir) || error("Localization data directory not found: $data_dir")

    files = sort(filter(f -> endswith(f, ".csv"), readdir(data_dir)))
    isempty(files) && error("No CSV files found in $data_dir")

    n       = length(files)
    ndigits = Base.ndigits(n)
    tar_path = endswith(output_abs, ".gz") ? output_abs[1:end-3] : output_abs * ".tar"

    isfile(tar_path) && rm(tar_path)
    cd(data_dir) do
        for (j, filename) in enumerate(files)
            print("\r  [$(lpad(j, ndigits))/$n] $filename$(repeat(' ', 20))")
            flush(stdout)
            flag = j == 1 ? "cf" : "rf"
            run(pipeline(`tar $flag $tar_path $filename`, stdout=devnull, stderr=devnull))
        end
    end
    println("\r  [$n/$n] Compressing...$(repeat(' ', 40))")
    flush(stdout)

    run(`gzip -f $tar_path`)

    checksum = bytes2hex(open(sha256, output_abs))
    println("\r  Done.$(repeat(' ', 50))")
    @info "Localization data bundle written" output=output_abs sha256=checksum

    return output_abs
end

"""
    download_localization_data(url_or_filename;
                               data_dir=get_localization_dir(),
                               expected_sha256=nothing)

Download a `.tar.gz` localization data bundle and extract it into `data_dir`.
The in-memory cache is cleared after extraction.

If `url_or_filename` does not start with `http://` or `https://`, it is
treated as a bare filename and the GitHub releases URL is prepended automatically.

# Examples
```julia
download_localization_data("Localization.tar.gz")
download_localization_data("https://example.com/Localization.tar.gz")
```
"""
function download_localization_data(url_or_filename::String;
                                    data_dir::String=get_localization_dir(),
                                    expected_sha256::Union{String,Nothing}=nothing)
    _GITHUB_RELEASES_BASE = "https://github.com/srossd/N83DBoot.jl/releases/latest/download/"
    url = startswith(url_or_filename, r"https?://") ? url_or_filename :
          _GITHUB_RELEASES_BASE * url_or_filename

    isdir(data_dir) || mkpath(data_dir)

    tmp = tempname() * ".tar.gz"
    try
        @info "Downloading localization data bundle..." url
        Downloads.download(url, tmp)

        if !isnothing(expected_sha256)
            actual = bytes2hex(open(sha256, tmp))
            if actual != lowercase(expected_sha256)
                error("SHA256 mismatch: expected $(expected_sha256), got $actual")
            end
            @info "Checksum verified" sha256=actual
        end

        run(`tar xzf $tmp -C $data_dir`)

        # Invalidate overlapping RAM cache entries
        data_dir_norm = rstrip(data_dir, '/')
        for key in collect(keys(_LOCALIZATION_RAM_CACHE))
            if rstrip(key[1], '/') == data_dir_norm
                delete!(_LOCALIZATION_RAM_CACHE, key)
            end
        end

        @info "Localization data bundle extracted" destination=data_dir
    finally
        isfile(tmp) && rm(tmp)
    end
end

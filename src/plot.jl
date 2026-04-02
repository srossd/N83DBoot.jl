# ============================================================================
# Plot types and generation for bootstrap results
# ============================================================================

using Plots
using LaTeXStrings
using DataFrames

# Okabe-Ito colorblind-safe palette
const OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#000000"]

const GOLDEN_RATIO = (1 + sqrt(5)) / 2
const PANEL_WIDTH  = 600  # pixels per panel column

_numeric_key(v::Number) = (false, Float64(v), "")
_numeric_key(v)         = (true, 0.0, string(v))
_tuple_sort_key(t::Tuple) = map(_numeric_key, t)

# ============================================================================
# Axis/Series Specifier Types
# ============================================================================

"""
    Metadata(key::String)

Axis or series specifier that reads from the problem metadata JSON.

# Examples
```julia
Metadata("N")    # ABJM rank
Metadata("k")    # Chern-Simons level
Metadata("cT")   # stress-tensor two-point coefficient
```
"""
struct Metadata
    key::String
end

"""
    Objective(key::String)

Axis or series specifier that reads from the objective functional or the solve result.

Use `"OBJ"` for the optimized objective value. For other keys, reads from the
`objective_properties` dict saved in the metadata (i.e. fields of the objective
functional).

# Examples
```julia
Objective("OBJ")   # the optimal bound value
Objective("b")     # the b parameter of an IntegralFunctional objective
```
"""
struct Objective
    key::String
end

const AxisSpec = Union{Metadata, Objective}

Base.show(io::IO, s::Metadata)  = print(io, "Metadata(\"$(s.key)\")")
Base.show(io::IO, s::Objective) = print(io, "Objective(\"$(s.key)\")")

# ============================================================================
# Plot Struct
# ============================================================================

"""
    Plot(name::String; axes, series=[], panels=[], panel_layout=:row, xlabel=nothing, ylabel=nothing)

A plot specification for bootstrap bounds.

When `solve_sdp` is called with a `Plot`, results are saved into a subdirectory
named `plot.name` inside the results directory.  `generate_plot` then scans that
subdirectory and assembles a figure.

# Fields
- `name::String`: Unique name for the plot (used as the subdirectory name)
- `axes::Tuple{AxisSpec,AxisSpec}`: `(x_spec, y_spec)` — what goes on each axis
- `series::Vector{AxisSpec}`: specifiers whose distinct values define separate colored traces
- `panels::Vector{AxisSpec}`: specifiers whose distinct values define separate subplots
- `panel_layout::Union{Symbol,Tuple{Int,Int}}`: arrangement of panels — `:row` (default,
  all in one row), `:col` (all in one column), or `(nrows, ncols)` for an explicit grid
- `xlabel`: custom x-axis label (default: auto-generated from axis spec); supports `L"..."` LaTeX strings
- `ylabel`: custom y-axis label (default: auto-generated from axis spec); supports `L"..."` LaTeX strings

# Axis/series/panel specifiers
- `Metadata("key")` — reads `key` from the top-level problem metadata (e.g. `"N"`, `"k"`, `"Lambda"`)
- `Objective("OBJ")` — the optimized objective value
- `Objective("field")` — a field of the objective functional (e.g. `"b"` for `IntegralFunctional`)

# Examples
```julia
# Bound on the stress-tensor OPE coefficient vs 1/N at fixed k
plt = Plot("lambda_vs_1overN";
    axes   = (Metadata("N"), Objective("OBJ")),
    series = [Metadata("k")],
    xlabel = L"N",
    ylabel = L"\\lambda^2")

# Multiple panels by Chern-Simons level
plt2 = Plot("lambda_panels";
    axes         = (Metadata("N"), Objective("OBJ")),
    series       = [Metadata("Lambda")],
    panels       = [Metadata("k")],
    panel_layout = :row)
```
"""
struct Plot
    name::String
    axes::Tuple{AxisSpec, AxisSpec}
    series::Vector{AxisSpec}
    panels::Vector{AxisSpec}
    panel_layout::Union{Symbol, Tuple{Int,Int}}
    xlabel::Union{AbstractString, Nothing}
    ylabel::Union{AbstractString, Nothing}

    function Plot(name::String;
                  axes::Tuple{AxisSpec, AxisSpec},
                  series::AbstractVector = AxisSpec[],
                  panels::AbstractVector = AxisSpec[],
                  panel_layout::Union{Symbol, Tuple{Int,Int}} = :row,
                  xlabel::Union{AbstractString, Nothing} = nothing,
                  ylabel::Union{AbstractString, Nothing} = nothing)
        @assert !isempty(name) "Plot name must not be empty"
        @assert panel_layout isa Tuple || panel_layout in (:row, :col) "panel_layout must be :row, :col, or a (nrows, ncols) tuple"
        new(name, axes, collect(AxisSpec, series), collect(AxisSpec, panels), panel_layout, xlabel, ylabel)
    end
end

Base.show(io::IO, p::Plot) =
    print(io, "Plot(\"$(p.name)\"; axes=$(p.axes), series=$(p.series), panels=$(p.panels))")

# ============================================================================
# Value Extraction from Metadata
# ============================================================================

"""
    extract_axis_value(spec, metadata) -> Any

Extract the numeric value described by `spec` from a loaded metadata `Dict`.
Returns `nothing` if the key is absent.
"""
function extract_axis_value(spec::Metadata, metadata::Dict)
    haskey(metadata, spec.key) && return metadata[spec.key]
    return 0
end

function extract_axis_value(spec::Objective, metadata::Dict)
    if spec.key == "OBJ"
        result = get(metadata, "result", nothing)
        isnothing(result) && return nothing
        return get(result, "OBJ", nothing)
    else
        props = get(metadata, "objective_properties", Dict())
        return get(props, spec.key, nothing)
    end
end

function extract_axis_value(spec::Metadata, row::DataFrames.DataFrameRow)
    col = spec.key
    hasproperty(row, Symbol(col)) || return 0
    v = row[col]
    (ismissing(v) || isnothing(v)) && return 0
    return v
end

function extract_axis_value(spec::Objective, row::DataFrames.DataFrameRow)
    col = spec.key == "OBJ" ? "OBJ" : spec.key
    hasproperty(row, Symbol(col)) || return nothing
    v = row[col]
    (ismissing(v) || isnothing(v)) && return nothing
    return v
end

function axis_label(spec::Metadata)
    return spec.key
end

function axis_label(spec::Objective)
    return spec.key == "OBJ" ? "Objective" : spec.key
end

function series_label(spec::Metadata, value)
    return "$(spec.key)=$(value)"
end

function series_label(spec::Objective, value)
    return "$(spec.key)=$(value)"
end

# ============================================================================
# Data Aggregation
# ============================================================================

"""
    get_data(plt::Plot; results_dir=nothing, filters=[]) -> DataFrame

Scan `results_dir/plt.name/` for all SDP results and return the data as a
`DataFrame`.  Each row is one solved SDP result.  Columns come from the
flattened metadata (top-level scalar fields and `objective_properties`
sub-fields) plus `"OBJ"` (the computed objective value) and `"direction"`.

# Arguments
- `plt::Plot`: The plot specification (used for the subdirectory name)
- `results_dir`: Results directory (default: from `get_config()`)
- `filters`: Vector of `AxisSpec => value` pairs to keep only matching rows.

# Examples
```julia
get_data(plt; filters = [Metadata("k") => 1, Metadata("Lambda") => 19])
```
"""
function get_data(plt::Plot;
                  results_dir::Union{String,Nothing}=nothing,
                  filters::AbstractVector=Pair{AxisSpec,Any}[])
    results_dir = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    plot_dir    = joinpath(results_dir, plt.name)

    isdir(plot_dir) || error("Plot directory does not exist: $plot_dir")

    all_files      = readdir(plot_dir, join=true)
    metadata_files = filter(f -> endswith(f, "metadata.json"), all_files)
    isempty(metadata_files) && @warn "No metadata files found in $plot_dir"

    rows = Dict{String, Any}[]
    for metadata_file in metadata_files
        hash        = replace(basename(metadata_file), "_metadata.json" => "")
        result_file = joinpath(plot_dir, "$(hash)_out.txt")
        isfile(result_file) || continue

        metadata           = JSON.parsefile(metadata_file)
        direction          = get(metadata, "direction", 1)
        # 3DN8Boot uses "fixed_contribution" (sum over all 4 fixed operators)
        fixed_contribution = get(metadata, "fixed_contribution", 0.0)

        parsed  = parse_out_txt(result_file, direction, fixed_contribution)
        optimal = parsed.status == "found primal-dual optimal solution"
        true_obj = (optimal && !isnothing(parsed.true_objective_value)) ?
            Float64(parsed.true_objective_value) : missing

        row = Dict{String, Any}()
        for (k, v) in metadata
            if k == "objective_properties" && isa(v, Dict)
                for (pk, pv) in v
                    row[String(pk)] = pv
                end
            elseif k == "operators" && isa(v, AbstractVector)
                long_ops = filter(op -> isa(op, Dict) && get(op, "type", "") == "Long", v)
                if !isempty(long_ops)
                    row["max_tau"]  = maximum(Float64(op["tau"]) for op in long_ops)
                    row["max_spin"] = maximum(Int(op["J"]) for op in long_ops)
                end
            elseif !isa(v, Dict) && !isa(v, AbstractVector)
                row[String(k)] = v
            end
        end
        row["OBJ"]       = true_obj
        row["direction"] = direction
        row["hash"]      = hash
        push!(rows, row)
    end

    isempty(rows) && return DataFrame()

    all_keys = sort(collect(union(keys.(rows)...)))
    df = DataFrame(Dict(k => [let v = get(r, k, missing); isnothing(v) ? missing : v; end for r in rows] for k in all_keys))

    for (spec, value) in filters
        col = spec.key
        if ismissing(value)
            df = df[ismissing.(df[!, col]), :]
        else
            col in names(df) || continue
            df = df[coalesce.(df[!, col] .== value, false), :]
        end
    end

    return df
end

# ============================================================================
# Plot Generation
# ============================================================================

"""
    _resolve_panel_layout(layout, n) -> Tuple{Int,Int}

Convert a `panel_layout` specifier to a concrete `(nrows, ncols)` grid size for
`n` panels.
"""
_resolve_panel_layout(layout::Symbol, n::Int) =
    layout == :row ? (1, n) :
    layout == :col ? (n, 1) :
    error("panel_layout must be :row, :col, or a (nrows, ncols) tuple")
_resolve_panel_layout(layout::Tuple{Int,Int}, ::Int) = layout

"""
    _build_single_panel(df, plt; kwargs...) -> Plots.Plot

Internal helper. Build one panel using the series grouping logic.
"""
function _build_single_panel(df::DataFrame, plt::Plot; kwargs...)
    x_spec, y_spec = plt.axes

    dir_groups = Dict{Tuple{Any,Int}, Vector{Tuple{Any,Any}}}()

    for row in eachrow(df)
        x = extract_axis_value(x_spec, row)
        y = extract_axis_value(y_spec, row)
        (isnothing(x) || isnothing(y)) && continue

        series_key = if isempty(plt.series)
            ()
        else
            vals = [extract_axis_value(s, row) for s in plt.series]
            any(isnothing, vals) && continue
            tuple(vals...)
        end

        direction = Int(row["direction"])
        push!(get!(dir_groups, (series_key, direction), Tuple{Any,Any}[]), (x, y))
    end

    sorted_series_keys = sort(unique(k[1] for k in keys(dir_groups)), by=_tuple_sort_key)
    @info "Series groups found: $(sorted_series_keys)"
    paired_keys = Set(sk for sk in sorted_series_keys
                      if haskey(dir_groups, (sk, -1)) && haskey(dir_groups, (sk, 1)))

    p = plot(; kwargs...)

    for (i, series_key) in enumerate(sorted_series_keys)
        base_name = if isempty(plt.series)
            plt.name
        else
            join([series_label(s, k) for (s, k) in zip(plt.series, series_key)], ", ")
        end

        if series_key in paired_keys
            lower_pts = sort(dir_groups[(series_key, -1)], by=pt -> pt[1])
            upper_pts = sort(dir_groups[(series_key,  1)], by=pt -> pt[1])

            xs_l = [pt[1] for pt in lower_pts]
            ys_l = [Float64(pt[2]) for pt in lower_pts]
            xs_u = [pt[1] for pt in upper_pts]
            ys_u = [Float64(pt[2]) for pt in upper_pts]

            common_xs = sort(collect(intersect(Set(xs_l), Set(xs_u))))
            if !isempty(common_xs)
                l_map = Dict(zip(xs_l, ys_l))
                u_map = Dict(zip(xs_u, ys_u))
                plot!(p, common_xs, [l_map[x] for x in common_xs];
                      fillrange=[u_map[x] for x in common_xs], fillalpha=0.2,
                      linealpha=0, label=false, color=i)
            end

            plot!(p, xs_u, ys_u; label=base_name, color=i, marker=false)
            plot!(p, xs_l, ys_l; label=false,     color=i, marker=false)
        else
            for dir in (1, -1)
                haskey(dir_groups, (series_key, dir)) || continue
                pts = sort(dir_groups[(series_key, dir)], by=pt -> pt[1])
                xs = [pt[1] for pt in pts]
                ys = [Float64(pt[2]) for pt in pts]
                dir_label = isempty(paired_keys) ? base_name :
                            "$(base_name) ($(dir == 1 ? "upper" : "lower"))"
                plot!(p, xs, ys; label=dir_label, color=i, marker=:circle)
            end
        end
    end

    return p
end

"""
    generate_plot(plt::Plot; results_dir=nothing, output_path=nothing, kwargs...) -> Plots.Plot

Scan `results_dir/plt.name/` for all SDP results and generate a figure using
Plots.jl.

# Arguments
- `plt::Plot`: The plot specification
- `results_dir`: Results directory (default: from `get_config()`)
- `output_path`: Path for the output file (default: `results_dir/plt.name/plt.name.png`)
- `filters`: Vector of `AxisSpec => value` pairs passed through to `get_data`.
- `kwargs...`: Additional keyword arguments forwarded to `Plots.plot` (e.g.
  `title`, `xlabel`, `ylabel`, `xlims`, `ylims`, `legend`, `size`)

# Examples
```julia
plt = Plot("lambda_vs_N";
    axes   = (Metadata("N"), Objective("OBJ")),
    series = [Metadata("k")])
p = generate_plot(plt; filters=[Metadata("k") => 1], ylims=(0, 20))
```
"""
function generate_plot(plt::Plot;
                       results_dir::Union{String,Nothing}=nothing,
                       output_path::Union{String,Nothing}=nothing,
                       filters::AbstractVector=Pair{AxisSpec,Any}[],
                       plot_kwargs...)
    results_dir = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    plot_dir    = joinpath(results_dir, plt.name)

    x_spec, y_spec = plt.axes

    df  = get_data(plt; results_dir=results_dir, filters=filters)
    out = isnothing(output_path) ? joinpath(plot_dir, "$(plt.name).png") : output_path

    kw_nt      = NamedTuple(plot_kwargs)
    user_title = get(kw_nt, :title, nothing)
    user_size  = get(kw_nt, :size,  nothing)
    xl = something(get(kw_nt, :xlabel, nothing), plt.xlabel, axis_label(x_spec))
    yl = something(get(kw_nt, :ylabel, nothing), plt.ylabel, axis_label(y_spec))
    base_kwargs = merge(
        (xlabel=xl, ylabel=yl, palette=OKABE_ITO),
        NamedTuple(k => v for (k, v) in pairs(kw_nt) if k ∉ (:title, :size, :xlabel, :ylabel))
    )

    if isempty(df)
        @warn "No complete results found in $plot_dir (need both _metadata.json and _out.txt)"
        default_size = (PANEL_WIDTH, round(Int, PANEL_WIDTH / GOLDEN_RATIO))
        p = plot(; title=plt.name, size=isnothing(user_size) ? default_size : user_size, base_kwargs...)
        savefig(p, out)
        return p
    end

    @info "Found $(nrow(df)) data points for plot '$(plt.name)'"

    local p
    if isempty(plt.panels)
        panel_title  = isnothing(user_title) ? plt.name : user_title
        default_size = (PANEL_WIDTH, round(Int, PANEL_WIDTH / GOLDEN_RATIO))
        p = _build_single_panel(df, plt; title=panel_title,
                                size=isnothing(user_size) ? default_size : user_size,
                                base_kwargs...)
    else
        panel_row_indices = Dict{Any, Vector{Int}}()
        for (i, row) in enumerate(eachrow(df))
            vals = [extract_axis_value(s, row) for s in plt.panels]
            any(isnothing, vals) && continue
            pk = tuple(vals...)
            push!(get!(panel_row_indices, pk, Int[]), i)
        end

        sorted_panel_keys = sort(collect(keys(panel_row_indices)), by=_tuple_sort_key)
        @info "Panel groups found: $(sorted_panel_keys)"

        panel_plots = map(sorted_panel_keys) do pk
            panel_label = join([series_label(s, k) for (s, k) in zip(plt.panels, pk)], ", ")
            panel_df    = df[panel_row_indices[pk], :]
            _build_single_panel(panel_df, plt; title=panel_label, base_kwargs...)
        end

        nrows, ncols = _resolve_panel_layout(plt.panel_layout, length(panel_plots))
        default_size = (ncols * PANEL_WIDTH, round(Int, nrows * PANEL_WIDTH / GOLDEN_RATIO))
        p = plot(panel_plots...; layout=(nrows, ncols),
                 size=isnothing(user_size) ? default_size : user_size,
                 left_margin=10Plots.mm)
    end

    savefig(p, out)
    @info "Plot written to $out"

    return p
end

# ============================================================================
# JointBoundsPlot — boundary of joint OPE coefficient allowed region
# ============================================================================

"""
    JointBoundsPlot(name; op_x, op_y, series=[], panels=[], ...)

Plot specification for a **joint** bootstrap bound on two free-operator OPE²
coefficients `λ²_A` (x-axis) and `λ²_B` (y-axis).

Collect results from many `SDPProblem`s whose `objective` is a
`LinearCombinationObjective(op_x => cos θ, op_y => sin θ)` for a sweep
of angles `θ ∈ [0, 2π)`.  Each solve returns a supporting half-plane
`cos θ · λ²_A + sin θ · λ²_B ≤ bound`, and `generate_plot` reconstructs
the joint allowed region as the intersection of all such half-planes.

# Fields
- `name::String`: Subdirectory under `results_dir` where data are stored
- `op_x::FreeOperator`: Operator whose λ² is plotted on the x-axis
- `op_y::FreeOperator`: Operator whose λ² is plotted on the y-axis
- `series::Vector{AxisSpec}`: Specifiers that group results into separate
  boundary curves (e.g. `[Metadata("num_integral")]` for with / without IC)
- `panels`, `panel_layout`: Same as `Plot`
- `xlabel`, `ylabel`: Axis labels (default: string representation of the operators)

# Examples
```julia
plt = JointBoundsPlot("semishort_joint";
    op_x   = SemishortAp0020(0),
    op_y   = SemishortAtwo0100(1),
    series = [Metadata("num_integral")],
    xlabel = L"\\lambda^2_{(A,+)_0}",
    ylabel = L"\\lambda^2_{(A,2)_1}")
generate_plot(plt)
```
"""
struct JointBoundsPlot
    name::String
    op_x::FreeOperator
    op_y::FreeOperator
    series::Vector{AxisSpec}
    panels::Vector{AxisSpec}
    panel_layout::Union{Symbol, Tuple{Int,Int}}
    xlabel::Union{AbstractString, Nothing}
    ylabel::Union{AbstractString, Nothing}

    function JointBoundsPlot(name::String;
                             op_x::FreeOperator,
                             op_y::FreeOperator,
                             series::AbstractVector = AxisSpec[],
                             panels::AbstractVector = AxisSpec[],
                             panel_layout::Union{Symbol, Tuple{Int,Int}} = :row,
                             xlabel::Union{AbstractString, Nothing} = nothing,
                             ylabel::Union{AbstractString, Nothing} = nothing)
        @assert !isempty(name) "JointBoundsPlot name must not be empty"
        @assert panel_layout isa Tuple || panel_layout in (:row, :col) "panel_layout must be :row, :col, or (nrows,ncols)"
        new(name, op_x, op_y,
            collect(AxisSpec, series), collect(AxisSpec, panels),
            panel_layout, xlabel, ylabel)
    end
end

Base.show(io::IO, p::JointBoundsPlot) =
    print(io, "JointBoundsPlot(\"$(p.name)\"; op_x=$(p.op_x), op_y=$(p.op_y))")

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _operator_coeff_key(op) -> String

Column name in the DataFrame for the coefficient of `op` in a
`LinearCombinationObjective` (must match `operator_key_string` in sdp.jl).
"""
_operator_coeff_key(op::SemishortAp0020)   = "coeff_SemishortAp0020_J$(op.J)"
_operator_coeff_key(op::SemishortAtwo0100) = "coeff_SemishortAtwo0100_J$(op.J)"
_operator_coeff_key(op::LongOperator)      = "coeff_Long_tau$(op.tau)_J$(op.J)"

"""
    _intersect_supporting_lines(nx1, ny1, f1, nx2, ny2, f2)

Intersect two half-plane boundaries `nx1·x + ny1·y = f1` and
`nx2·x + ny2·y = f2`.  Returns `nothing` for nearly-parallel normals.
"""
function _intersect_supporting_lines(nx1::Float64, ny1::Float64, f1::Float64,
                                     nx2::Float64, ny2::Float64, f2::Float64)
    det = nx1 * ny2 - nx2 * ny1
    abs(det) < 1e-12 && return nothing
    x = (f1 * ny2 - f2 * ny1) / det
    y = (f2 * nx1 - f1 * nx2) / det
    return (x, y)
end

"""
    _reconstruct_boundary(normals, bounds) -> (xs, ys)

Reconstruct the boundary polygon of the joint allowed region from a collection
of supporting half-planes.  Each half-plane is `n_i · p ≤ bound_i` where
`n_i = (nx_i, ny_i)` is the outward normal direction.

Returns `(xs, ys)` — the x and y coordinates of the polygon vertices in
angle-sorted order, suitable for a `Plots.plot!(..., seriestype=:shape)` call.
"""
function _reconstruct_boundary(normals::Vector{Tuple{Float64,Float64}},
                                bounds::Vector{Float64})
    n = length(normals)
    n < 3 && return Float64[], Float64[]

    # Sort half-planes by angle of outward normal
    order   = sortperm(normals, by = t -> atan(t[2], t[1]))
    normals = normals[order]
    bounds  = bounds[order]

    tol = 1e-10
    is_feasible(pt) = all(normals[k][1] * pt[1] + normals[k][2] * pt[2] ≤ bounds[k] + tol
                          for k in 1:n)

    # Find the first non-redundant constraint to use as the starting point.
    # A redundant constraint has no feasible intersection with any other line.
    start = 0
    for i in 1:n
        for d in 1:n-1
            pt = _intersect_supporting_lines(normals[i]..., bounds[i],
                                             normals[mod1(i+d,n)]..., bounds[mod1(i+d,n)])
            if !isnothing(pt) && is_feasible(pt)
                start = i
                break
            end
        end
        start != 0 && break
    end
    start == 0 && return Float64[], Float64[]

    # Sweep around in angle order: from constraint i, find the nearest j (by
    # offset d = 1, 2, …) whose intersection with i is feasible.  That
    # intersection is the next polygon vertex; advance i = j and repeat until
    # we return to `start`.
    xs = Float64[]
    ys = Float64[]
    i  = start
    while true
        for d in 1:n-1
            j  = mod1(i + d, n)
            pt = _intersect_supporting_lines(normals[i]..., bounds[i],
                                             normals[j]..., bounds[j])
            if !isnothing(pt) && is_feasible(pt)
                push!(xs, pt[1])
                push!(ys, pt[2])
                i = j
                break
            end
        end
        i == start && break
    end

    return xs, ys
end

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

"""
    get_data(plt::JointBoundsPlot; results_dir=nothing, filters=[]) -> DataFrame

Load all joint-bounds results for `plt` from `results_dir/plt.name/`.
"""
function get_data(plt::JointBoundsPlot;
                  results_dir::Union{String,Nothing}=nothing,
                  filters::AbstractVector=Pair{AxisSpec,Any}[])
    # Reuse standard data loading by wrapping in a pass-through Plot
    proxy = Plot(plt.name;
                 axes   = (Metadata("_x"), Metadata("_y")),
                 series = plt.series,
                 panels = plt.panels)
    return get_data(proxy; results_dir=results_dir, filters=filters)
end

# ---------------------------------------------------------------------------
# Panel builder
# ---------------------------------------------------------------------------

"""
    _build_joint_panel(df, plt; kwargs...) -> Plots.Plot

Build one joint-bounds panel with (possibly multiple) series boundary curves.
"""
function _build_joint_panel(df::DataFrame, plt::JointBoundsPlot; kwargs...)
    key_x = _operator_coeff_key(plt.op_x)
    key_y = _operator_coeff_key(plt.op_y)

    # Group rows by series key
    series_groups = Dict{Any, Vector{DataFrames.DataFrameRow}}()
    for row in eachrow(df)
        sk = if isempty(plt.series)
            ()
        else
            vals = [extract_axis_value(s, row) for s in plt.series]
            any(isnothing, vals) && continue
            tuple(vals...)
        end
        push!(get!(series_groups, sk, DataFrames.DataFrameRow[]), row)
    end

    sorted_keys = sort(collect(keys(series_groups)), by = _tuple_sort_key)
    @info "JointBoundsPlot series found: $(sorted_keys)"

    p = plot(; kwargs...)

    for (i, sk) in enumerate(sorted_keys)
        rows    = series_groups[sk]
        normals = Tuple{Float64,Float64}[]
        bvals   = Float64[]

        for row in rows
            hasproperty(row, Symbol(key_x)) || continue
            hasproperty(row, Symbol(key_y)) || continue
            hasproperty(row, :OBJ)          || continue
            cx  = row[Symbol(key_x)]
            cy  = row[Symbol(key_y)]
            obj = row[:OBJ]
            (ismissing(cx) || ismissing(cy) || ismissing(obj)) && continue
            push!(normals, (Float64(cx), Float64(cy)))
            push!(bvals,   Float64(obj))
        end

        if isempty(normals)
            @warn "No valid data for series key $sk"
            continue
        end

        xs, ys = _reconstruct_boundary(normals, bvals)
        isempty(xs) && continue

        lbl = if isempty(plt.series)
            plt.name
        else
            join([series_label(s, k) for (s, k) in zip(plt.series, sk)], ", ")
        end

        color = OKABE_ITO[mod1(i, length(OKABE_ITO))]
        plot!(p, xs, ys;
              seriestype = :shape,
              fillalpha  = 0.35,
              linewidth  = 1.5,
              label      = lbl,
              color      = color)
    end

    return p
end

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

"""
    generate_plot(plt::JointBoundsPlot; results_dir=nothing, output_path=nothing, ...) -> Plots.Plot

Reconstruct and plot the joint allowed region in `(λ²_{op_x}, λ²_{op_y})` space
from a sweep of `LinearCombinationObjective` SDPB results.
"""
function generate_plot(plt::JointBoundsPlot;
                       results_dir::Union{String,Nothing}=nothing,
                       output_path::Union{String,Nothing}=nothing,
                       filters::AbstractVector=Pair{AxisSpec,Any}[],
                       plot_kwargs...)
    results_dir = isnothing(results_dir) ? get_config().sdpb_results_dir : results_dir
    plot_dir    = joinpath(results_dir, plt.name)

    df  = get_data(plt; results_dir=results_dir, filters=filters)
    out = isnothing(output_path) ? joinpath(plot_dir, "$(plt.name).png") : output_path

    kw_nt      = NamedTuple(plot_kwargs)
    user_title = get(kw_nt, :title, nothing)
    user_size  = get(kw_nt, :size,  nothing)
    xl = something(get(kw_nt, :xlabel, nothing), plt.xlabel, string(plt.op_x))
    yl = something(get(kw_nt, :ylabel, nothing), plt.ylabel, string(plt.op_y))
    base_kwargs = merge(
        (xlabel=xl, ylabel=yl, palette=OKABE_ITO, legend=:bottomright),
        NamedTuple(k => v for (k, v) in pairs(kw_nt) if k ∉ (:title, :size, :xlabel, :ylabel))
    )
    default_size = (PANEL_WIDTH, round(Int, PANEL_WIDTH / GOLDEN_RATIO))

    if isempty(df)
        @warn "No complete results found in $plot_dir"
        p = plot(; title=plt.name,
                 size=isnothing(user_size) ? default_size : user_size,
                 base_kwargs...)
        savefig(p, out)
        return p
    end

    @info "Found $(nrow(df)) data points for joint bounds plot '$(plt.name)'"

    local p
    if isempty(plt.panels)
        panel_title = isnothing(user_title) ? plt.name : user_title
        p = _build_joint_panel(df, plt;
                               title = panel_title,
                               size  = isnothing(user_size) ? default_size : user_size,
                               base_kwargs...)
    else
        panel_row_indices = Dict{Any, Vector{Int}}()
        for (i, row) in enumerate(eachrow(df))
            vals = [extract_axis_value(s, row) for s in plt.panels]
            any(isnothing, vals) && continue
            pk = tuple(vals...)
            push!(get!(panel_row_indices, pk, Int[]), i)
        end

        sorted_panel_keys = sort(collect(keys(panel_row_indices)), by = _tuple_sort_key)
        @info "JointBoundsPlot panels found: $(sorted_panel_keys)"

        panel_plots = map(sorted_panel_keys) do pk
            panel_label = join([series_label(s, k) for (s, k) in zip(plt.panels, pk)], ", ")
            panel_df    = df[panel_row_indices[pk], :]
            _build_joint_panel(panel_df, plt; title=panel_label, base_kwargs...)
        end

        nrows, ncols = _resolve_panel_layout(plt.panel_layout, length(panel_plots))
        p = plot(panel_plots...;
                 layout     = (nrows, ncols),
                 size       = isnothing(user_size) ?
                              (ncols * PANEL_WIDTH, round(Int, nrows * PANEL_WIDTH / GOLDEN_RATIO)) :
                              user_size,
                 left_margin = 10Plots.mm)
    end

    savefig(p, out)
    @info "Joint bounds plot written to $out"
    return p
end

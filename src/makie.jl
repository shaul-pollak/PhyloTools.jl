function make_coords(tr::Node{I,T}; u::Bool=false) where {I,T}
    d = treepositions(tr, upwards=u)
    toloop = prewalk(tr)
    v = Vector{Tuple}()
    for n in toloop
        if isroot(n)
            if isfinite(distance(n))
                (x, y) = d[id(n)]
                if !u
                    push!(v, (0.0, y))
                    push!(v, (x, y))
                elseif u
                    push!(v, (x, 0.0))
                    push!(v, (x, y))
                end
            end
        else
            (x1, y1) = d[id(n)]
            (x2, y2) = d[id(parent(n))]
            if !u
                push!(v, (x2, y2))
                push!(v, (x2, y1))
                push!(v, (x2, y1))
                push!(v, (x1, y1))
            elseif u
                push!(v, (x2, y2))
                push!(v, (x1, y2))
                push!(v, (x1, y2))
                push!(v, (x1, y1))
            end
        end
    end
    return v
end

function MakieCore.plot!(p::MakieCore.Plot(Phy.Node))
    tr = p[1][]
    up = get(p, :up, false)[]
    x = make_coords(tr, u=up)
    xs = [x[1] for x in x]
    ys = [x[2] for x in x]
    clr = get(p, :color, :black)
    lw = get(p, :linewidth, 1)
    MakieCore.linesegments!(p, xs, ys, color=clr, linewidth=lw)
    return p
end

function filter_treepos(tps::Dict{I,J}, d::Dict{<:AbstractString,K}, tr::Node{I,T}) where {I,J,K,T}
    lvs = filter(x -> name(x) ∈ keys(d), getleaves(tr))
    x, y, c = Vector{Float64}(), Vector{Float64}(), Vector{K}()
    for l in lvs
        p = tps[id(l)]
        push!(x, p[1])
        push!(y, p[2])
        push!(c, d[name(l)])
    end
    return x, y, c
end

MakieCore.@recipe(TipPoints, tr, d) do scene
    MakieCore.Theme(
        Axis=(
            # remove spines
            leftspinevisible=false,
            rightspinevisible=false,
            bottomspinevisible=false,
            topspinevisible=false,
            # remove y axis decorations
            ylabelvisible=false,
            yticklabelsvisible=false,
            yticksvisible=false,
            ygridvisible=false,
            yminorgridvisible=false,
            yminorticksvisible=false,
            # remove x axis decorations
            xlabelvisible=false,
            xticklabelsvisible=false,
            xticksvisible=false,
            xgridvisible=false,
            xminorgridvisible=false,
            xminorticksvisible=false
        )
    )
end

function MakieCore.plot!(p::TipPoints)
    tr = p[1][]
    d = p[2][]
    up = get(p, :up, false)[]
    cm = get(p, :colormap, :viridis)
    tps = treepositions(tr, upwards=up)
    x, y, c = filter_treepos(tps, d, tr)
    MakieCore.scatter!(p, x, y, color=c, colormap=cm)
    return p
end



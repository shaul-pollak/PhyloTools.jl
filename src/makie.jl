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

function MakieCore.plot!(p::MakieCore.Plot(Node))
    tr = p[1][]
    up = get(p, :up, false)
    x = make_coords(tr, u=up)
    xs = [x[1] for x in x]
    ys = [x[2] for x in x]
    clr = get(p, :color, :black)
    lw = get(p, :linewidth, 1)
    MakieCore.linesegments!(p, xs, ys, color=clr, linewidth=lw)
    return p
end

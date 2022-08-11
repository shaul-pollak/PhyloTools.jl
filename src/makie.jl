import MakieCore.convert_arguments
import MakieCore.plottype
import MakieCore.linesegments
import MakieCore.AbstractPlot

function make_coords(tr::Node{I,T}, u::Bool=false) where {I,T}
    d = treepositions(tr, upwards=u)
    toloop = prewalk(tr)
    v = Vector{Tuple}()
    for n in toloop
        if isroot(n)
            if isfinite(distance(n))
                (x, y) = d[id(n)]
                push!(v, (x, y - distance(n)))
                push!(v, (x, y))
            end
        else
            (x1, y1) = d[id(n)]
            (x2, y2) = d[id(parent(n))]
            push!(v, (x2, y2))
            push!(v, (x2, y1))
            push!(v, (x2, y1))
            push!(v, (x1, y1))
        end
    end
    return v
end

function convert_arguments(P::Type{<:AbstractPlot}, tr::Node)
    x = make_coords(tr)
    xs = [x[1] for x in x]
    ys = [x[2] for x in x]
    return convert_arguments(P, xs, ys)
end

plottype(::Node) = linesegments

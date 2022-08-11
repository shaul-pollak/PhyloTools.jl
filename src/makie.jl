function make_coords(tr::Node{I,T}, u::Bool=false) where {I,T}
    d = treepositions(tr, upwards=u)
    toloop = prewalk(tr)
    v = Vector{Point2}()
    for n in toloop
        if isroot(n)
            if isfinite(distance(n))
                (x, y) = d[id(n)]
                push!(v, Point(x, y - distance(n)))
                push!(v, Point(x, y))
            end
        else
            (x1, y1) = d[id(n)]
            (x2, y2) = d[id(parent(n))]
            push!(v, Point2(x2, y2))
            push!(v, Point2(x2, y1))
            push!(v, Point2(x2, y1))
            push!(v, Point2(x1, y1))
        end
    end
    return v
end

function Makie.convert_arguments(P::Type{<:AbstractPlot}, tr::Node)
    x = make_coords(tr)
    return Makie.convert_arguments(P, x)
end

Makie.plottype(::Node) = linesegments


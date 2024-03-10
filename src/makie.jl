using GeometryBasics: Point2f, Polygon


function treepositions(tr::Node{I,T}, transform=false; upwards=false) where {I,T}
  o = postwalk(tr);
  n = length(o);
  h = getheights(tr);
  maxh = maximum(values(h));
  nodepos = Dict{I,Tuple{Float64,Float64}}();
  i = 0.;
  for n in o
    if isleaf(n) 
      i += 1.;
      hn = transform ? maxh : h[id(n)];
    else
      hn = h[id(n)];
    end
    nodepos[id(n)] = upwards ? (i,hn) : (hn, i);
  end
  function walk(n)
    isleaf(n) && return nodepos[id(n)]
    cs = map(walk, children(n));
    xn = upwards ? nodepos[id(n)][2] : nodepos[id(n)][1];
    yn = (upwards ? sum(first.(cs)) : sum(last.(cs))) / length(cs);
    nodepos[id(n)] = upwards ? (yn,xn) : (xn, yn);
  end
  walk(tr)
  # if the root branch has non-zero length we add this:
  return nodepos
end


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

function MakieCore.plot!(p::MakieCore.Plot(PhyloTools.Node))
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


"""
    gene(x::Matrix{Any}, y::Number)

draws gene arrows in a given orientation at a specified y value
x is a matrix with a row for each gene.
the columns are: start, end, orientation (a string, "+" or "-"), and target gene designation (a Bool)

"""
function gene(x, y)
    yt = y + 0.35
    yb = y - 0.35
    if x[2] - x[1] > 100
        if x[3] == "+"
            p1 = Point2f(x[1], yt)
            p2 = Point2f(x[2] - 100, yt)
            p3 = Point2f(x[2], y)
            p4 = Point2f(x[2] - 100, yb)
            p5 = Point2f(x[1], yb)
            return Polygon([p1, p2, p3, p4, p5])
        else
            p1 = Point2f(x[1], y)
            p2 = Point2f(x[1] + 100, yt)
            p3 = Point2f(x[2], yt)
            p4 = Point2f(x[2], yb)
            p5 = Point2f(x[1] + 100, yb)
            return Polygon([p1, p2, p3, p4, p5])
        end
    else
        if x[3] == "+"
            p1 = Point2f(x[1], yt)
            p2 = Point2f(x[2], y)
            p3 = Point2f(x[1], yb)
            return Polygon([p1, p2, p3])
        else
            p1 = Point2f(x[1], y)
            p2 = Point2f(x[2], yt)
            p3 = Point2f(x[2], yb)
            return Polygon([p1, p2, p3])
        end
    end
end

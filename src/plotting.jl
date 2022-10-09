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

@recipe function f(tree::Node; showlabels=false, internal=true, fs=3, upwards=false, nodes="all")
  d = treepositions(tree,upwards=upwards)
  framestyle --> :none
  grid --> false
  legend --> false
  widen --> false
  if nodes == "leaves"
    toloop = getleaves(tree)
  elseif nodes == "nodes"
    toloop = filter(x -> !isleaf(x),prewalk(tree))
  else
    toloop = prewalk(tree)
  end
  for n in toloop
    if isroot(n) 
      if isfinite(distance(n))
        (x, y) = d[id(n)];
        @series begin
          seriestype --> :path
          seriescolor --> :black
          [(x, y-distance(n)), (x,y)]
        end
      end
    else
      @series begin
        seriestype --> :path
        seriescolor --> :black
        (x1, y1) = d[id(n)]
        if nodes == "all"
          (x2, y2) = d[id(parent(n))]
          upwards ? [(x2, y2), (x1, y2), (x1, y1)] : [(x2, y2), (x2, y1), (x1, y1)]
        else
          [(x1,y1)]
        end
      end
    end
  end
  if showlabels==true
    ns = internal ? prewalk(tree) : getleaves(tree)
    anns = [(d[id(n)]..., (" " * name(n), fs, :left)) for n in ns]
    @series begin
      annotations := anns
      [], []
    end
  end
end

@userplot TipPoint
@userplot NodePoint
@recipe function f(tp::TipPoint)
  nodes := "leaves";
  upwards --> false;
  seriestype := :scatter
  tree = tp.args[1];
  l = length(tp.args)
  if l==1
    tree
  elseif l==2
    colDict = tp.args[2]
    tree, colDict
  end
end
@recipe function f(np::NodePoint)
  tree = np.args[1];
  nodes := "nodes";
  upwards --> true;
  seriestype := :scatter
  tree
end

@recipe function f(tree::Node,colDict::Dict; upwards=false)
  d = treepositions(tree,upwards=upwards)
  framestyle --> :none
  grid --> false
  legend --> false
  if eltype(keys(colDict))==String
    toloop = filter(x -> name(x) ∈ keys(colDict),prewalk(tree))
    for n in toloop
      @series begin
        markercolor := colDict[name(n)]
        markerstrokecolor --> colDict[name(n)]
        (x1, y1) = d[id(n)]
        [(x1,y1)]
      end
    end
  elseif eltype(keys(colDict))==UInt32
    toloop = filter(x -> id(x) ∈ keys(colDict),prewalk(tree))
    for n in toloop
      @series begin
        marker_z := colDict[id(n)]
        markerstrokecolor := nothing
        (x1, y1) = d[id(n)]
        [(x1,y1)]
      end
    end
  end
end

# @userplot Hist
# @recipe function f(h::Hist)
#   x = h.args[1];
#   legend_position --> :outertopright
#   @series begin
#     seriestype := :stephist
#     seriescolor --> "black"
#     fillcolor --> "darkred"
#     fillrange := 0
#     ylims := (0,Inf)
#     x
#   end
# end

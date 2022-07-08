function treepositions(tr::N, transform=false; upwards=false) where N
    o = postwalk(tr);
    n = length(o);
    h = getheights(tr);
    maxh = maximum(values(h));
    nodepos = Dict{N,Tuple{Float64,Float64}}();
    i = 0.;
    for n in o
        if isleaf(n) 
            i += 1.;
            hn = transform ? maxh : h[id(n)];
        else
            hn = h[id(n)];
        end
        nodepos[n] = upwards ? (i,hn) : (hn, i);
    end
    function walk(n)
        isleaf(n) && return nodepos[n]
        cs = map(walk, children(n));
        xn = upwards ? nodepos[n][2] : nodepos[n][1];
        yn = (upwards ? sum(first.(cs)) : sum(last.(cs))) / length(cs);
        nodepos[n] = upwards ? (yn,xn) : (xn, yn);
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
        (x, y) = d[n];
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
        (x1, y1) = d[n]
        if nodes == "all"
          (x2, y2) = d[parent(n)]
          upwards ? [(x2, y2), (x1, y2), (x1, y1)] : [(x2, y2), (x2, y1), (x1, y1)]
        else
          [(x1,y1)]
        end
      end
    end
  end
  if showlabels==true
    ns = internal ? prewalk(tree) : getleaves(tree)
    anns = [(d[n]..., (" " * name(n), fs, :left)) for n in ns]
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
  markerstrokewidth := 0;
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
  markerstrokewidth := 0;
  seriestype := :scatter
  tree
end

@recipe function f(tree::Node,colDict::Dict; 
    showlabels=false, fs=3, upwards=false)
  d = treepositions(tree,upwards=upwards)
  framestyle --> :none
  grid --> false
  legend --> false
  toloop = filter(x -> name(x) ∈ keys(colDict),prewalk(tree))
  for n in toloop
    @series begin
      markercolor := colDict[name(n)]
      markerstrokecolor := colDict[name(n)]
      (x1, y1) = d[n]
      [(x1,y1)]
    end
  end
end



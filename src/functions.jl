leafnames(x::Node) = [name(a) for a in getleaves(x)];

function Ntip(x::Node)
	length(getleaves(x))
end

function coph!(out,orig,target,nd,ld)
  @inbounds begin
    pth = getpath(orig,target)[1]; 
    focal = [name(orig)];
    for a in 2:length(pth)
      k = pth[a];
      nfc = [x for x in children(k) if id(x)!=id(orig)];
      nonfocal = [isleaf(x) ? Set([name(x)]) : keys(nd[id(x)]) for x in nfc];
      nonfocal = reduce(union,nonfocal) |> collect;
      ii = [ld[x] for x in focal];
      jj = [ld[x] for x in nonfocal];
      for i in 1:length(ii)
        @simd for j in 1:length(jj)
          out[ii[i],jj[j]] = out[jj[j],ii[i]] = nd[id(k)][focal[i]] + nd[id(k)][nonfocal[j]];
        end
      end
      # go over subtrees
      if !isroot(k)
        for n in nfc
          if !isleaf(n) 
            coph!(out,getleaves(n)[1],n,nd,ld);
          end
        end
      end
      focal = union(focal,nonfocal);
      orig = k;
    end
  end
  return out
end

function cophenetic(x::Node{T,I}) where {T,I}
  out = zeros(Float64,Ntip(x),Ntip(x));
  leaves = leafnames(x);
  # leaf name to row/col index
  ld = Dict(leaves[i] => i for i in eachindex(leaves));
  # node dict
  nd = Dict{T,Dict{String,Float64}}();
  # populate node dict using postwalk traversal
  o = postwalk(x);
  @inbounds begin
    for k in 1:length(o)
      n = o[k];
      if !isleaf(n)
        lv = Vector{String}();
        vv = Vector{Float64}();
        for c in children(n)
          if isleaf(c)
            push!(lv,name(c));
            push!(vv,distance(c));
          else
            append!(lv,keys(nd[id(c)]));
            append!(vv,values(nd[id(c)]) .+ distance(c));
          end
        end
        # create entry in in
        nd[id(n)] = Dict(zip(lv,vv));
      end
    end
    ## populate output
    lvs = [getleaves(y) for y in children(x)];
    coph!(out,lvs[1][1],getroot(lvs[1][1]),nd,ld);
    coph!(out,lvs[2][1],getroot(lvs[2][1]),nd,ld);
  end
  return out
end

cart2lin(x::CartesianIndex,dims::Tuple) = dims[1]*(x[2]-1) + x[1];

make_independent_tree(x::Node) = readnw(nwstr(x));

function reroot!(node::Node{T,I}) where {T,I}
  curr_root = getroot(node);
  p = parent(node);
  # return original tree if node is already outgroup
  p == curr_root && length(children(curr_root)) == 2 && return p
  i = maximum(id, prewalk(curr_root)) + 1; # assign new root maximum(id) + 1
  newroot = Node(T(i), NewickData()); # create new root node
  alongtheway = [node];
  while !isnothing(p)
    push!(alongtheway,p);
    p = parent(p);
  end
  while length(alongtheway) > 1
    p = pop!(alongtheway);
    # this will become the parent of p
    a = filter(x->id(x) == id(alongtheway[end]), children(p))[1];
    d = distance(a);
    delete!(p, id(a));
    if node != a
      push!(a, p)
      setdistance!(p, d)
    else
      push!(newroot, p);
      push!(newroot, a);
      setdistance!(p, d/2)
      setdistance!(a, d/2)
      break
    end
  end 
  if length(children(curr_root)) == 1  
    # when the previous root was not multifurcating, we end up with a
    # superfluous node, so we delete it
    p = parent(curr_root)
    delete!(p, id(curr_root))
    push!(p, curr_root[1])
  end
  return newroot
end

function getdivs(x::Node{I,T}) where {I,T}
  o = postwalk(x);
  out = Dict(id(n) => 0.0 for n in o);
  for n in o
    if !isleaf(n)
      out[id(n)] = maximum([out[id(c)]+distance(c) for c in children(n)]);
    end
  end
  return out
end

function cluster_tree(tree::Node{I,T},f::Number,tp="color") where {I,T}
  lns = leafnames(tree);
  out = Dict(lns[i] => i for i in eachindex(lns));
  ds = getdivs(tree);
  tot = ds[id(tree)];
  th = tot/f;
  po = prewalk(tree);
  i = 1;
  while length(po)>0
    n = po[1];
    if ds[id(n)] <= th
      lvs = leafnames(n);
      for k in lvs
        out[k] = i;
      end
      setdiff!(po,postwalk(n));
      i += 1;
    else
      popfirst!(po);
    end
  end
  cv = Colors.distinguishable_colors(min(50,maximum(values(out))));
  out2 = Dict(l => cv[1] for l in keys(out));
  for k in keys(out)
    i = mod(out[k],length(cv))+1
    out2[k] = cv[i];
  end
  if tp=="color"
    return out2
  elseif tp=="int"
    return out
  else
    return out,out2
  end
end

function extract_with_tips(N::Node{T,I},x::Vector{F},nparents::Int=0) where {T,I,F<:AbstractString}
  checkall(p,x) = all([ln ∈ x for ln in leafnames(p)]);
  ls = getleaves(N);
  c = ls[name.(ls) .== x[1]][1];
  p = parent(c);
  flg=true;
  while flg
    if checkall(p,x)
      c = p;
      p = parent(c);
    else
      flg = false;
    end
  end
  if nparents>0
    for i in 1:nparents
      c = parent(c);
    end
  end
  return c
end

function extract_with_tips(N::Node{T,I},x::AbstractDict{K,V},nparents::Int=0) where {T,I,K<:AbstractString,V}
  x = collect(keys(x));
  extract_with_tips(N,x,nparents)
end

# create dictionary where keys are internal nodes and values are dictionaries
# where keys are tips and values are distance to that tip
function max_dist_to_descending_tip(tr::Node{I,T}) where {I,T}
  pw = postwalk(tr);
  d = Dict(id(x) => Dict{I,Float64}() for x in pw if !isleaf(x));
  for n in pw # leaves to root (find farthest descending node)
    if !isleaf(n) 
      for c in children(n)
        if isleaf(c)
          d[id(n)][id(c)] = distance(c);
        else
          for k in keys(d[id(c)])
            d[id(n)][k] = d[id(c)][k] + distance(c);
          end
        end
      end
    end
  end
  return d
end

function add_incoming_tips!(tr::Node{I,T}, d::Dict{I,Dict{I, Float64}}) where {I,T}
  pw = [x for x in prewalk(tr) if !isleaf(x) && !isroot(x)];
  for n in pw
    if isroot(n)
      cn = [x for x in children(n) if !isleaf(x)];
      for i in eachindex(cn)
        c = cn[i];
        for k in keys(d[id(c)])
          d[id(n)][k] = d[id(c)][k] + distance(c);
        end
      end
    else
      np = parent(n);
      for k in setdiff(keys(d[id(np)]),keys(d[id(n)]))
        d[id(n)][k] = d[id(np)][k] + distance(n);
      end
    end
  end
  return nothing
end

function find_midpoint_node(tr::Node{I,T}, d::Dict{I,Dict{I,Float64}}) where {I,T}
  m = Inf;
  i = I(0);
  for (k,vd) in d
    mv = maximum(values(vd));
    if mv < m
      m = mv;
      i = k;
    end
  end
  return filter(x -> id(x)==i,prewalk(tr))[1]
end

function midpoint_root(tr::Node{I,T}) where {I,T}
  all([isleaf(x) for x in children(tr)]) && return tr
  d = max_dist_to_descending_tip(tr)
  add_incoming_tips!(tr,d)
  n = find_midpoint_node(tr,d)
  tr2 = id(n)==id(tr) ? tr : reroot!(n)
  d = max_dist_to_descending_tip(tr2);
  cs = children(tr2);
  md = [isleaf(x) ? 0. : maximum(values(d[id(x)])) for x in cs];
  mdiff = diff(sort(md))[1];
  mdiff==0 && return tr2
  i = md[1]>md[2] ? 1 : 2;
  j = md[1]>md[2] ? 2 : 1;
  th = distance(cs[1]);
  newdist = mdiff/2 <= th ? mdiff/2 : th/2;
  setdistance!(cs[i],max(eps(),newdist))
  setdistance!(cs[j],distance(cs[j]) + newdist/2)
  return tr2
end

function whichmax(x::Matrix{T}) where T<:Number
  o = zeros(Int8,size(x,1))
  m = maximum(x,dims=2);
  for i in 1:size(x,1)
    o[i] = findfirst(x[i,:] .== m[i]);
  end
  return o
end

function whichmin(x::Matrix{T}) where T<:Number
  o = zeros(Int8,size(x,1))
  m = minimum(x,dims=2);
  for i in 1:size(x,1)
    o[i] = findfirst(x[i,:] .== m[i]);
  end
  return o
end

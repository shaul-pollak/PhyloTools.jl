function order_dict(x,y)
  # x is the gene tree
  # y is the species tree
  o1 = collect(values(x));
  kx = collect(keys(x));
  o2 = [y[split(k,'_')[1]] for k in kx];
  return o1,o2
end

function count(x)
  mx = maximum(x)
  out = zeros(Int64,mx);
  for a in x
    out[a] += 1
  end
  return out
end

function count(x,y)
  if length(x)!=length(y)
    error("both vectors must be the same length")
  end
  mx = maximum(x)
  my = maximum(y)
  out = zeros(Int64,mx,my)
  for i in 1:length(x) 
    out[x[i],y[i]] += 1
  end
  return out
end

function entropy(x)
  p = count(values(x));
  p = p./sum(p); # turn in probability
  -sum([x*log2(x) for x in p if x>0])
end

function entropy(vx,vy)
  pxy = count(vx,vy);
  pxy = pxy ./ sum(pxy); # turn into probability
  -sum([x*log2(x) for x in pxy if x>0])
end

function treeVI!(querytree::Node,stc::AbstractDict,out::Vector{Any};
    fv::UnitRange{Int64}=2:1.5:30)
  Xs = [cluster_tree(querytree,f,"int") for f in fv];
  Y = order_dict(Xs[1],stc)[2];
  nvi = zeros(Float64,length(Xs));
  Hx = [entropy(x) for x in Xs];
  Hy = entropy(Y);
  @inbounds begin
    for i in eachindex(Xs)
      Hxy = entropy(collect(values(Xs[i])),Y);
      vi = 2*Hxy - Hx[i] - Hy
      nvi[i] = 1-vi/Hxy;
    end
  end
  m = maximum(nvi);
  out[1] = m
  out[2] = findfirst(nvi .== m)
  return nothing
end

function treeVI!(querytree::Node,stc::Vector{Dict{String,Int64}},
    out::Vector{Vector{Any}}; fv=2:1.5:30)
  Xs = [cluster_tree(querytree,f,"int") for f in fv];
  Ys = [order_dict(Xs[1],stc[i])[2] for i in eachindex(stc)];
  nvi = zeros(Float64,length(Xs),length(stc));
  Hx = [entropy(x) for x in Xs];
  Hy = [entropy(y) for y in Ys];
  @inbounds begin
    for i in eachindex(Xs)
      for j in eachindex(Ys)
        Hxy = entropy(collect(values(Xs[i])),Ys[j])
        vi = 2*Hxy - Hx[i] - Hy[j]
        nvi[i,j] = 1-vi/Hxy;
      end
    end
  end
  m = maximum(nvi,dims=1);
  out[1] = m[:]
  out[2] = zeros(Float64,length(m));
  out[2] = [isfinite(m[i]) ? fv[findfirst(nvi[:,i] .== m[i])] : out[2][i] for i in 1:length(m)]
  return nothing
end

function duplication_score(x::Node{I,T},th::Float64) where {I,T}
  ng(x::Node) = [split(l,'_')[1] for l in leafnames(x)] |> unique |> length;
  tot = ng(x);
  any([ng(y)/tot for y in children(x)] .> th)
end



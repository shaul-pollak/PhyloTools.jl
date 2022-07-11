function order_dict(x,y)
  # x is the gene tree
  # y is the species tree
  o1 = collect(values(x));
  kx = collect(keys(x));
  o2 = [haskey(y,k[1]) ? y[k[1]] : NaN for k in split.(kx,'_')];
  return o1,o2
end

function count(x)
  mx = maximum((x[isfinite.(x)]));
  out = zeros(Int64,mx);
  for a in x[isfinite.(x)]
    out[a] += 1
  end
  return out
end

function count(x,y)
  if length(x)!=length(y)
    error("both vectors must be the same length")
  end
  I = @. isfinite(x) & isfinite(y);
  mx = maximum(x[I])
  my = maximum(y[I])
  out = zeros(Int64,mx,my)
  for i in findall(I)
    out[x[i],y[i]] += 1
  end
  return out
end

function entropy(x::Dict{X,Y}) where {X,Y<:Number}
  p = count(collect(values(x)));
  p = p./sum(p); # turn in probability
  max(0.,-sum([x*log2(x) for x in p if x>0]))
end

function entropy(x::Vector{T}) where T<:Number
  p = count(x);
  p = p./sum(p); # turn in probability
  max(0.,-sum([x*log2(x) for x in p if x>0]))
end

function entropy(vx,vy)
  pxy = count(vx,vy);
  pxy = pxy ./ sum(pxy); # turn into probability
  max(0.,-sum([x*log2(x) for x in pxy if x>0]))
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
        if Hxy != 0
          vi = 2*Hxy - Hx[i] - Hy[j]
          nvi[i,j] = 1-vi/Hxy;
        end
      end
    end
  end
  m = maximum(nvi,dims=1);
  out[1] = m[:]
  out[2] = zeros(Float64,length(m));
  out[2] = [isfinite(m[i]) ? fv[findfirst(nvi[:,i] .== m[i])] : out[2][i] for i in 1:length(m)]
  return nothing
end

function duplication_score(x::Node{I,T}) where {I,T}
  ng(x::Node) = [split(l,'_')[1] for l in leafnames(x)] |> unique |> length;
  tot = ng(x);
  sum([ng(y)/tot for y in children(x)])/length(children(x)) # mean
end

function treeVI(trees::Vector{Node{T,I}}, stc::Vector{Dict{String,K}}; 
    fv::Vector{F}=[2,10,30,50],dupth::Float64=.5) where {T,I,K<:Integer,F<:Number}
  vi = zeros(Float64,length(trees),length(stc));
  res = zeros(Float64,length(trees),length(stc));
  tmp = [Vector{Vector{Any}}(undef,2) for _ in 1:Threads.nthreads()];
  Threads.@threads for i in eachindex(trees)
    @inbounds begin
      if duplication_score(trees[i]) >= dupth
        t2 = [Vector{Vector{Any}}(undef,2) for _ in 1:2]
        for k in 1:2
          treeVI!(trees[i][k],stc,t2[k];fv=fv)
        end
        tmp[Threads.threadid()] = t2[1][1]>t2[1][2] ? t2[1] : t2[2];
      else
        treeVI!(trees[i],stc,tmp[Threads.threadid()];fv=fv)
      end
      for j in 1:length(tmp[Threads.threadid()][1])
        if isfinite(tmp[Threads.threadid()][1][j])
          vi[i,j] = tmp[Threads.threadid()][1][j]
          res[i,j] = tmp[Threads.threadid()][2][j]
        end
      end
    end
  end
  return (treevi = vi, resolution = res)
end

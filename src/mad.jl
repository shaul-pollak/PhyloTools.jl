function make_distmat(tr::Node{T,I}) where {T,I}
    pw = postwalk(tr);
    nn = length(pw);
    o = zeros(Float64, nn, nn);
    id2i = Dict(id(pw[i]) => i for i in eachindex(pw));
    for n in filter(!isleaf,pw)
        ni = id2i[id(n)]
        nc = children(n)
        for c in nc
            ci = id2i[id(c)]
            o[ni, ci] = o[ci, ni] = distance(c)
            if !isleaf(c)
                for cc in postwalk(c)
                    cc==c && continue
                    cci = id2i[id(cc)]
                    o[ni, cci] = o[cci, ni] = distance(c) + o[ci, cci]
                end
            end
        end
        for i in eachindex(nc)
            for j in i+1:length(nc)
                Threads.@threads for nc1 in postwalk(nc[i])
                    i1 = id2i[id(nc1)]
                    for nc2 in postwalk(nc[j])
                        i2 = id2i[id(nc2)]
                        o[i1, i2] = o[i2, i1] = o[i1, ni] + o[i2, ni]
                    end
                end
            end
        end
    end
    return o, id2i
end


function reldev(n::Node{T,I}, dists::Matrix{Float64}, id2i::Dict{T,Int}) where {T,I}
    nc = children(n)
    ni = haskey(id2i, id(n)) ? id2i[id(n)] : 0
    s = 0.0
    for i in 1:length(nc)
        for l1 in getleaves(nc[i])
            i1 = id2i[id(l1)]
            d_i_l1 = ni != 0 ? dists[ni, i1] : (dists[id2i[id(nc[i])], i1] + distance(nc[i]))
            for j in i+i:length(nc)
                for l2 in getleaves(nc[j])
                    i2 = id2i[id(l2)]
                    if dists[i1, i2] > 0
                        nom = 2 * d_i_l1
                        denom = dists[i1, i2]
                        v = (nom / denom - 1)^2
                        s += v
                    end
                end
            end
        end
    end
    return s
end

function walk2root(tr::Node{T,I}, nid::T) where {T,I}
    n = filter(x -> id(x) == nid, prewalk(tr))[]
    o = Vector{T}()
    p = parent(n)
    while !isnothing(p)
        push!(o, id(p))
        p = parent(p)
    end
    return o
end



function node2graph(tr::Node{T,I}, id2i::Dict{T,Int64}) where {T,I}
    g = SimpleGraph(length(id2i))
    for n in postwalk(tr)
        if !isroot(n) | length(children(n))==3
            nid = id2i[id(n)]
            for c in children(n)
                add_edge!(g, nid, id2i[id(c)])
            end
        elseif length(children(n))==2
            ids = [id2i[x] for x in id.(children(n))]
            add_edge!(g, ids...)
        end
    end
    return g
end


graph_children(g::SimpleGraph, c, p) = filter(x -> x != p, Graphs.neighbors(g, c))

function findtips!(g, c::T, p::T, o::Vector{T}=Vector{T}()) where {T}
    ccv = graph_children(g, c, p)
    if isempty(ccv)
        push!(o, c)
        return nothing
    end
    for cc in ccv
        findtips!(g, cc, c, o)
    end
end


function mad_graph(g::SimpleGraph, dists::Matrix{Float64}, leafind::Vector{Int64})
    elist = edges(g) |> collect
    dss = zeros(Float64, length(elist))
    l2i = Dict(x[2] => x[1] for x in enumerate(leafind))
    nterms = length(l2i) * (length(l2i)-1)
    lk = ReentrantLock()
    prg = Progress(length(elist); desc="going over all tree edges", dt=.1)
    rsd = [Dict{Int,Float64}() for _ in 1:Graphs.nv(g)]
    Threads.@threads for tid in 1:length(elist)
        e = elist[tid]
        cv = [(src(e), dst(e)); (dst(e), src(e))]
        dss[tid] = 0.0
        while !isempty(cv)
            (p, c) = pop!(cv)
            if Graphs.degree(g, c) != 1
                if !haskey(rsd[c],p)
                    lock(lk) do
                        rsd[c][p] = anc_norm_div(dists, c, straddling_node_tips(g, c, p)...)
                    end
                end
                dss[tid] += rsd[c][p]
            end
            for c2 in graph_children(g, c, p)
                Graphs.degree(g, c2) != 1 && push!(cv, (c, c2))
            end
        end
        ρ = min_rho(g, e, dists)
        d_io = dists[dst(e), src(e)] * ρ
        bs, cs = straddling_edge_tips(g, e)
        i = src(e)
        for b in bs
            for c in cs
                tpsh = (2 * (dists[i, b] + d_io) / dists[b, c] - 1)^2
                dss[tid] += isnan(tpsh) ? 0.0 : tpsh
            end
        end
        dss[tid] = sqrt(dss[tid] / nterms)
        next!(prg)
    end
    return elist[sortperm(dss)[1]]
end

function min_rho(g::SimpleGraph, e::Edge, dists::Matrix{Float64})
    tid = Threads.threadid()
    ii, jj = straddling_edge_tips(g, e)
    nom = 0.0
    denom = 0.0
    for x in ii
        for y in jj
            d_bc_sq = dists[x, y] == 0.0 ? 0.0 : dists[x, y]^(-2)
            nom += (dists[x, y] - 2 * dists[src(e), x]) * d_bc_sq
            denom += d_bc_sq
        end
    end
    denom = 2 * dists[dst(e), src(e)] * denom
    o = nom / denom
    return min(max(0, o), 1)
end

function straddling_edge_tips(g::SimpleGraph, e::Edge)
    i = src(e)
    j = dst(e)
    ii = Int[]
    jj = Int[]
    findtips!(g, i, j, ii)
    findtips!(g, j, i, jj)
    return ii, jj
end

function straddling_node_tips(g, i, p)
    ics = graph_children(g, i, p)
    length(ics) == 0 && return nothing
    ii = Int[]
    jj = Int[]
    findtips!(g, ics[1], i, ii)
    findtips!(g, ics[2], i, jj)
    return ii, jj
end

function anc_norm_div(dists, i, bs, cs)
    s = 0.0
    for b in bs
        for c in cs
            tpsh = (2 * dists[i, b] / dists[b, c] - 1)^2
            s += isnan(tpsh) ? 0.0 : tpsh
        end
    end
    return s
end


function mad(tr::Node{T,I}) where {T,I}
    wtr = deepcopy(tr)
    dists, id2i = make_distmat(wtr)
    i2id = id.(postwalk(wtr))
    leafind = [id2i[k] for k in id.(getleaves(wtr))]
    g = node2graph(wtr, id2i)
    oe = mad_graph(g, dists, leafind)
    pw = postwalk(wtr)
    n1 = pw[src(oe)]
    n1id = id(n1)
    n2 = pw[dst(oe)]
    tr2 = (parent(n1) == n2 ? reroot!(n1) : reroot!(n2)) |> ladderize!
    # fix root distances
    ρ = min_rho(g, oe, dists)
    d_io = dists[dst(oe), src(oe)] * ρ
    d_jo = dists[dst(oe), src(oe)] * (1-ρ)
    if id(children(tr2)[1])==n1id
        children(tr2)[1].data.distance = d_io
        children(tr2)[2].data.distance = d_jo
    else
        children(tr2)[1].data.distance = d_jo
        children(tr2)[2].data.distance = d_io
    end
    fix_node_ids!(tr2)
    return tr2
end

function fix_node_ids!(tr::Node{T,I}) where {T,I}
    # leaves get ids in 1:ntip range
    leaf_range = T.(collect(1:ntip(tr)))
    for (i,l) in enumerate(getleaves(tr))
        l.id = leaf_range[i]
    end
    # nodes get ids in range ntip+1:nnode from the root to the tips (preorder)
    i = T(ntip(tr)+1)
    for x in filter(!isleaf,prewalk(tr))
        x.id = i
        i+=1
    end
end
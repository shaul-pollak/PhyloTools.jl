leafnames(x::Node) = [name(a) for a in getleaves(x)];

function Ntip(x::Node)
    length(getleaves(x))
end

dist2root(n) = sum([isfinite(distance(x)) ? distance(x) : 0 for x in getpath(n)]);

function node_children_dist!(n, nd, ld, o)
    @inbounds begin
        for c in children(n)
            cci = id(c)
            ocs = filter(x -> id(x) != cci, children(n))
            for l1 in leafnames(c)
                i = ld[l1]
                d1 = isleaf(c) ? distance(c) : nd[cci][l1] + distance(c)
                for oc in ocs
                    for l2 in leafnames(oc)
                        j = ld[l2]
                        d2 = isleaf(oc) ? distance(oc) : nd[id(oc)][l2] + distance(oc)
                        o[i, j] = o[j, i] = d1 + d2
                    end
                end
            end
        end
    end
end

function _outdist(x::Node{I,T}, nd::Dict{I,Dict{String,Float64}}) where {I,T}
    nl = Ntip(x)
    o = zeros(Float64, nl, nl)
    lns = leafnames(x)
    ld = Dict(lns[i] => i for i in eachindex(lns))
    # calculate pairwise distance between leaves under node
    for n in filter(!isleaf, prewalk(x))
        node_children_dist!(n, nd, ld, o)
    end
    return o
end

function cophenetic(x::Node{T,I}) where {T,I}
    leaves = leafnames(x)
    # leaf name to row/col index
    ld = Dict(leaves[i] => i for i in eachindex(leaves))
    # node dict
    nd = Dict{T,Dict{String,Float64}}()
    # populate node dict using postwalk traversal
    o = postwalk(x)
    @inbounds begin
        for k in 1:length(o)
            n = o[k]
            if !isleaf(n)
                lv = Vector{String}()
                vv = Vector{Float64}()
                for c in children(n)
                    if isleaf(c)
                        push!(lv, name(c))
                        push!(vv, distance(c))
                    else
                        ci = id(c)
                        append!(lv, keys(nd[ci]))
                        append!(vv, values(nd[ci]) .+ distance(c))
                    end
                end
                # create entry in nd
                nd[id(n)] = Dict(zip(lv, vv))
            end
        end
        ## populate output
        out = _outdist(x, nd)
    end
    return out
end

cart2lin(x::CartesianIndex, dims::Tuple) = dims[1] * (x[2] - 1) + x[1];

make_independent_tree(x::Node) = readnw(nwstr(x));

function reroot!(node::Node{T,I}) where {T,I}
    curr_root = getroot(node)
    p = parent(node)
    # return original tree if node is already outgroup
    p == curr_root && length(children(curr_root)) == 2 && return p
    i = maximum(id, prewalk(curr_root)) + 1 # assign new root maximum(id) + 1
    newroot = Node(T(i), NewickData()) # create new root node
    alongtheway = [node]
    while !isnothing(p)
        push!(alongtheway, p)
        p = parent(p)
    end
    while length(alongtheway) > 1
        p = pop!(alongtheway)
        # this will become the parent of p
        a = filter(x -> id(x) == id(alongtheway[end]), children(p))[1]
        d = distance(a)
        delete!(p, id(a))
        if node != a
            push!(a, p)
            setdistance!(p, d)
        else
            push!(newroot, p)
            push!(newroot, a)
            setdistance!(p, d / 2)
            setdistance!(a, d / 2)
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
    o = postwalk(x)
    out = Dict(id(n) => 0.0 for n in o)
    for n in o
        if !isleaf(n)
            out[id(n)] = maximum([out[id(c)] + distance(c) for c in children(n)])
        end
    end
    return out
end

function cluster_tree(tree::Node{I,T}, f::Number, tp="color") where {I,T}
    lns = leafnames(tree)
    out = Dict(lns[i] => i for i in eachindex(lns))
    ds = getdivs(tree)
    tot = ds[id(tree)]
    th = tot / f
    po = prewalk(tree)
    i = 1
    while length(po) > 0
        n = po[1]
        if ds[id(n)] <= th
            lvs = leafnames(n)
            for k in lvs
                out[k] = i
            end
            setdiff!(po, postwalk(n))
            i += 1
        else
            popfirst!(po)
        end
    end
    cv = Colors.distinguishable_colors(min(50, maximum(values(out))))
    out2 = Dict(l => cv[1] for l in keys(out))
    for k in keys(out)
        i = mod(out[k], length(cv)) + 1
        out2[k] = cv[i]
    end
    if tp == "color"
        return out2
    elseif tp == "int"
        return out
    else
        return out, out2
    end
end

function extract_with_tips(N::Node{T,I}, x::Vector{F}, nparents::Int=0) where {T,I,F<:AbstractString}
    checkall(p, x) = all([ln ∈ x for ln in leafnames(p)])
    ls = getleaves(N)
    c = ls[name.(ls).==x[1]][1]
    p = parent(c)
    flg = true
    while flg
        if checkall(p, x)
            c = p
            p = parent(c)
        else
            flg = false
        end
    end
    if nparents > 0
        for i in 1:nparents
            c = parent(c)
        end
    end
    return c
end

function extract_with_tips(N::Node{T,I}, x::AbstractDict{K,V}, nparents::Int=0) where {T,I,K<:AbstractString,V}
    x = collect(keys(x))
    extract_with_tips(N, x, nparents)
end

# create dictionary where keys are internal nodes and values are dictionaries
# where keys are tips and values are distance to that tip
function max_dist_to_descending_tip(tr::Node{I,T}) where {I,T}
    pw = postwalk(tr)
    d = Dict{I,Float64}()
    for n in pw # leaves to root (find farthest descending node)
        i = id(n)
        if !isleaf(n)
            m = 0
            for c in children(n)
                j = id(c)
                tmp = d[j] + distance(c)
                if tmp > m
                    m = tmp
                end
            end
            d[i] = m
        else
            d[i] = 0.0
        end
    end
    return d
end

function cross_parent_maxdist(p, n, d, m)
    for c in children(p)
        if id(c) != id(n)
            cdist = haskey(d, id(c)) ? d[id(c)] + distance(c) : distance(c)
            if cdist > m
                m = cdist
            end
        end
    end
    return distance(n) + m
end

function walk_up(n, d, m)
    totdist = 0.0
    p = n
    while !isnothing(parent(p))
        tmp = cross_parent_maxdist(parent(p), p, d, m) + totdist
        if tmp > m
            m = tmp
        end
        totdist += distance(p)
        p = parent(p)
    end
    return m
end

#function add_incoming_tips(tr::Node{I,T}, d::Dict{I, Float64}) where {I,T}
#  d2 = deepcopy(d);
#  pw = filter(!isroot,prewalk(tr));
#  for n in pw
#    m = 0.
#    if isroot(parent(n))
#      m = cross_parent_maxdist(tr,n,d,m);
#    else
#      m = walk_up(n,d,m) 
#    end 
#    if (haskey(d,id(n)) && m>d[id(n)]) | !haskey(d,id(n))
#      d2[id(n)] = m
#    end
#  end
#  return d2
#end


function add_incoming_tips(tr::Node{I,T}, d::Dict{I,Float64}) where {I,T}
    d2 = deepcopy(d)
    pw = filter(!isroot, prewalk(tr))
    for n in pw
        ni = id(n)
        m = d[ni]
        m2p = 0.0
        cn = n
        while !isnothing(parent(cn))
            m2p += distance(cn)
            cnc = filter(x -> id(x) != id(cn), children(parent(cn))) # finds sisters
            for c in cnc
                v = d[id(c)] + distance(c) + m2p
                if v > m
                    m = v
                end
            end
            cn = parent(cn)
        end
        d2[ni] = m
    end
    return d2
end


function midpoint_root(tr::Node{I,T}) where {I,T}
    all([isleaf(x) for x in children(tr)]) && return tr
    d = max_dist_to_descending_tip(tr)
    d = add_incoming_tips(tr, d)
    m = Inf
    i = 0
    for (k, v) in d
        if v < m
            i = k
            m = v
        end
    end
    n = filter(x -> id(x) == i, prewalk(tr))[1]
    if !isroot(n)
        tr = reroot!(n)
    end
    d = max_dist_to_descending_tip(tr)
    cs = children(tr)
    md = [isleaf(x) ? 0.0 : maximum(values(d[id(x)])) for x in cs]
    mdiff = diff(sort(md))[1]
    mdiff == 0 && return tr
    i = md[1] > md[2] ? 1 : 2
    j = md[1] > md[2] ? 2 : 1
    th = distance(cs[1])
    newdist = mdiff / 2 <= th ? mdiff / 2 : th / 2
    setdistance!(cs[i], max(eps(), newdist))
    setdistance!(cs[j], distance(cs[j]) + newdist / 2)
    return tr
end

function whichmax(x::Matrix{T}) where {T<:Number}
    o = zeros(Int8, size(x, 1))
    m = maximum(x, dims=2)
    for i in 1:size(x, 1)
        o[i] = findfirst(x[i, :] .== m[i])
    end
    return o
end

function whichmin(x::Matrix{T}) where {T<:Number}
    o = zeros(Int8, size(x, 1))
    m = minimum(x, dims=2)
    for i in 1:size(x, 1)
        o[i] = findfirst(x[i, :] .== m[i])
    end
    return o
end

function MRCA(tree::Node{I,T}, tips::Vector{K}) where {I,T,K<:AbstractString}
    findleaf(tree, tip) = filter(x -> name(x) == tip, getleaves(tree))[1]
    p1 = getpath(findleaf(tree, tips[1]))
    paths = [id.(getpath(findleaf(tree, x))) for x in tips]
    i = reduce(intersect, paths)[1]
    return p1[@. id(p1) == i][1]
end

function ladderize!(phy::Node{I,T}, reverse_ord::Bool=true) where {I,T}
    # make size dict
    sdict = make_node_size_dict(phy)
    # walk from root to tips
    no = filter(x -> !isleaf(x), prewalk(phy))
    for n in no
        tmp = children(n)[sortperm(.-[sdict[k] for k in id.(children(n))])]
        if reverse_ord
            tmp = reverse(tmp)
        end
        n.children = tmp
    end
    return phy
end

function make_node_size_dict(x::Node{I,T}) where {I,T}
    o = postwalk(x)
    out = Dict(id(n) => Int(0) for n in o)
    for n in o
        if !isleaf(n)
            for c in children(n)
                out[id(n)] += out[id(c)]
            end
        else
            out[id(n)] = 1
        end
    end
    return out
end

function readfasta(p::String)
    o = Vector{Tuple{String,String}}()
    hb = IOBuffer()
    sb = IOBuffer()
    io = open(p, "r")
    flg = false
    for l in eachline(io)
        if l[1] == '>'
            if flg
                push!(o, (String(take!(hb)), String(take!(sb))))
            end
            write(hb, l[2:end])
        else
            if !flg
                flg = true
            end
            write(sb, l)
        end
    end
    close(io)
    push!(o, (String(take!(hb)), String(take!(sb))))
    return o
end

"""
    read_gff(p::String, tp::String="")

reads a gff file into a vector of strings.
the string tp filters the resulting vector by the type of entry
examples of tp strings are "CDS" or "gene"
"""
function read_gff(p::String, tp::String="")
    r = readlines(p)
    filter!(!startswith("#"), r)
    rs = split.(r, '\t')
    filter!(x -> length(x) >= 3, rs)
    if tp != ""
        filter!(x -> x[3] == tp, rs)
    end
    return rs
end

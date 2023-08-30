leafnames(x::Node) = [name(a) for a in getleaves(x)]

using Base: readuntil_vector!
ntip(x::Node)::Int = length(getleaves(x))
nnode(x::Node)::Int = length(postwalk(x))

dist2root(n) = sum([isfinite(distance(x)) ? distance(x) : 0 for x in getpath(n)])

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

function _outdist(x::Node{T,I}, nd::Dict{T,Dict{String,Float64}}) where {T,I}
    nl = ntip(x)
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
    while length(alongtheway) > 0
        p = pop!(alongtheway)
        # this will become the parent of p
        a = filter(x -> id(x) == id(alongtheway[end]), children(p))[]
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
        c = children(curr_root)[]
        d = distance(curr_root)
        p = parent(curr_root)
        delete!(p, id(curr_root))
        push!(p, curr_root[1])
        !isnan(d) && setdistance!(c, distance(c) + d)
    end
    return newroot
end

function getdivs(x::Node{T,I}) where {T,I}
    o = postwalk(x)
    out = Dict(id(n) => 0.0 for n in o)
    for n in o
        if !isleaf(n)
            out[id(n)] = maximum([out[id(c)] + distance(c) for c in children(n)])
        end
    end
    return out
end

function cluster_tree(tree::Node{T,I}, f::Number, tp="color") where {T,I}
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
function max_dist_to_descending_tip(tr::Node{T,I}) where {T,I}
    pw = postwalk(tr)
    d = Dict{T,Float64}()
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

function add_incoming_tips(tr::Node{T,I}, d::Dict{T,Float64}) where {T,I}
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


function midpoint_root(tr::Node{T,I}) where {T,I}
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

function MRCA(tree::Node{T,I}, tips::Vector{K}) where {T,I,K<:AbstractString}
    findleaf(tree, tip) = filter(x -> name(x) == tip, getleaves(tree))[1]
    p1 = getpath(findleaf(tree, tips[1]))
    paths = [id.(getpath(findleaf(tree, x))) for x in tips]
    i = reduce(intersect, paths)[1]
    return p1[@. id(p1) == i][1]
end

function ladderize!(phy::Node{T,I}, reverse_ord::Bool=true) where {T,I}
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

function make_node_size_dict(x::Node{T,I}) where {T,I}
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
    fx(x::String)::Vector{String} = split(x, "\n>")
    fx2(x::String)::Vector{String} = split(x, '\n')
    r = read(p, String)
    rs = fx(r)
    rs[1] = rs[1][2:end] # first record still has '>'
    o = Dict{String,String}()
    sizehint!(o, length(rs))
    lk = ReentrantLock()
    Threads.@threads for r in rs
        x = fx2(r)
        if length(x) > 1
            lock(lk) do
                o[x[1]] = join(x[2:end])
            end
        end
    end
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


"""
    readclu(p::T, opt::Int=0) where {T<:AbstractString}

read protein clustering generated by diamond deepclust or mmseqs.
in these files, the first column is the cluster name and the second is the protein id
columns are separated by tabs.

opt=0 means both plot2clu and clu2prot dicts
opt=1 means clu2prot only
opt=2 means prot2clu only
opt=3 means from new format of clu
"""
function readclu(p::T; opt::Int=0, id::Int=50) where {T<:AbstractString}
    !(opt ∈ [0, 1, 2, 3]) && error("opt must be in {0, 1, 2, 3}")
    if opt==3
        return clu3(p, id)
    end
    open(p, "r") do io
        # functions
        fx(x)::Vector{String} = split(x, '\t')
        function safepush!(d, k, v)
            if haskey(d, k)
                push!(d[k], v)
            else
                d[k] = [v]
            end
        end
        function create_chunks(data::Vector{UInt8})
            # create chunks that end on newline or eof
            chunk_size = (length(data) / Threads.nthreads()) |> ceil |> Int
            chunks = [((i-1)*chunk_size+1):min(i * chunk_size, length(data)) for i in 1:Threads.nthreads()]
            nlis = zeros(Int, length(chunks))
            Threads.@threads for i in 1:Threads.nthreads()
                ci = chunks[i][end-1000:end]
                @views chunk = data[ci]
                nlis[i] = ci[findlast(chunk .== UInt('\n'))]
            end
            for i in eachindex(chunks)
                if i == 1
                    nc = chunks[i].start:nlis[i]
                elseif i == length(chunks)
                    nc = nlis[i-1]+1:chunks[i].stop
                else
                    nc = nlis[i-1]+1:nlis[i]
                end
                chunks[i] = nc
            end
            return chunks
        end
        function findnl(x, j, c)
            while (x[j]) != c
                j += 1
            end
            return j
        end
        findtb(x, c) = findfirst(x .== c)
        nl = UInt8('\n')
        tb = UInt8('\t')
        # memory map data
        data = mmap(io)
        chunks = create_chunks(data)
        # prepare output
        lk = ReentrantLock()
        o1 = Dict{String,Vector{String}}()
        o2 = Dict{String,String}()
        Threads.@threads for i in 1:Threads.nthreads()
            @views d = data[chunks[i]]
            j1 = 1
            while j1 < length(d)
                j2 = findnl(d, j1, nl)
                l = d[j1:j2]
                j1 = j2 + 1
                tbi = findtb(l, tb)
                lc = String(l[1:tbi-1])
                lp = String(l[tbi+1:end]) |> strip
                lock(lk) do
                    if opt == 0
                        safepush!(o1, lc, lp)
                        o2[lp] = lc
                    elseif opt == 1
                        safepush!(o1, lc, lp)
                    elseif opt == 2
                        o2[lp] = lc
                    end
                end
            end
        end
        if opt == 0
            return o1, o2
        elseif opt == 1
            return o1
        elseif opt == 2
            return o2
        end
    end
end

function clu3(p, id)
    f(x::String) = split(x, '\t')[2:end]
    function fx(x::Vector{SubString{String}}, T::DataType)
        o = zeros(T, length(x))
        @inbounds begin
            Threads.@threads for i in eachindex(x)
                o[i] = @views parse(T, x[i])
            end
        end
        return o
    end
    function read_zst(p)
        cio = open(p,"r")
        lns = readlines(ZstdDecompressorStream(io))
        close(cio)
        return lns
    end
    hdrs = read_idx0(p);
    printstyled("\nreading clu into memory\n"; color=:green)
    lns = endswith(".zst")(p) ? read_zst(p) : readlines(p)
    i = startswith("min_id=$id").(lns) |> findfirst
    lns = lns[i]
    d = f(lns) |> x -> fx(x, UInt32)
    printstyled("initializing clu2prot dict\n"; color=:green)
    o = Dict(hdrs[i] => String[] for i in IterTools.distinct(d) if i>0)
    prg = Progress(length(d); desc="populating clu2prot dict")
    for i in eachindex(d)
        if d[i]>0
            k = hdrs[d[i]]
            v = @views o[k]
            push!(v, hdrs[i])
        end
        next!(prg)
    end
    return o
end

function read_idx0(p)
    s(x)::String = split(x,'\t')[1]
    io = open(p,"r")
    lns = s.(readlines(ZstdDecompressorStream(io)))
    close(io)
    return lns
end
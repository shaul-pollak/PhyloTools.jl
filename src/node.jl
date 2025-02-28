# defaults
name(x) = string(x)
distance(x) = NaN
support(x) = NaN

"""
    NewickData{T,S<:AbstractString}
A simple container for the allowed fields in a newick tree. Those are the
`distance` (expected number of substitutions, time, what have you), `support`
(e.g. bootstrap support value, posterior clade probability) and `name` (leaf
names). Note that internal node labels are not allowed.
"""
mutable struct NewickData{T,S}
    distance::T
    support ::T
    name    ::S
end

NewickData(; d=NaN, s=NaN, n="") = NewickData(promote(d, s)..., n)

name(n::NewickData) = n.name
support(n::NewickData) = n.support
distance(n::NewickData) = n.distance
setname!(n::NewickData, s) = n.name = string(s)
setsupport!(n::NewickData, x) = n.support = x
setdistance!(n::NewickData, x) = n.distance = x

"""
    Node{I,T}
A node of a (phylogenetic) tree.
"""
mutable struct Node{I,T}
    id      ::I
    data    ::T
    parent  ::Node{I,T}
    children::Vector{Node{I,T}}

    Node(id::I, data::T) where {I,T} = new{I,T}(id, data)
    Node(id::I, data::T, p::Nothing) where {I,T} = new{I,T}(id, data)
    function Node(id::I, data::T, p) where {I,T}
        n = new{I,T}(id, data, p)
        push!(p, n)
        return n
    end
end

Node(i; kwargs...) = Node(i, NewickData(; kwargs...))
Node(i, p::Node; kwargs...) = Node(i, NewickData(; kwargs...), p)

id(n::Node) = n.id
data(n::Node) = n.data
name(n::Node) = name(n.data)
support(n::Node) = support(n.data)
distance(n::Node) = distance(n.data)

function reindex!(n::Node{I}) where I
    for (i,n) in enumerate(prewalk(n))
        n.id = I(i)
    end
end

setname!(n::Node, s) = setname!(n.data, s)
setsupport!(n::Node, s) = setsupport(n.data, s)
setdistance!(n::Node, d) = setdistance!(n.data, d)

isroot(n::Node) = !isdefined(n, :parent)
isleaf(n::Node) = !isdefined(n, :children) || length(n) == 0
degree(n::Node) = isleaf(n) ? 0 : length(n.children)
nv(n::Node) = length(prewalk(n))

function getnode(f, n::Node) 
    for x in prewalk(n)
        f(x) && return x
    end
end

Base.pop!(n::Node) = pop!(n.children)
Base.last(n::Node) = last(children(n))
Base.first(n::Node) = first(children(n))
Base.length(n::Node) = length(n.children)
Base.eltype(::Type{Node{T}}) where T = Node{T}

Base.getindex(n::Node, i) = n.children[i]
Base.setindex!(n::Node, x, i) = n.children[i] = x

Base.parent(n::Node) = isdefined(n, :parent) ? n.parent : nothing
Base.parent(root, n::Node) = isdefined(n, :parent) ? n.parent : nothing

function Base.delete!(n::Node, i::Integer)
    deleteat!(n.children, findfirst(c->id(c) == i, children(n)))
end

function Base.delete!(n::Node, m::Node) 
    deleteat!(n.children, findfirst(c->c == m, children(n)))
end

Base.push!(n::Node, args...) = for x in args; push!(n, x); end

function Base.push!(n::Node, m::Node)
    isdefined(n, :children) ? 
        push!(n.children, m) : 
        n.children = [m]
    m.parent = n
end

function getroot(n::Node)
    while !isroot(n) 
        n = parent(n) 
    end
    return n
end

function findleaf(n::Node, s::String)
    for x in getleaves(n)
        name(x) == s && return x
    end
end

function sister(n::Node)
    xs = children(parent(n))
    xs[findfirst(x->x != n, xs)]
end

Base.show(io::IO, n::Node{I,<:NewickData}) where I = write(io, "$(nwstr(n; internal=true))")
Base.show(io::IO, n::Node) = write(io,"Node($(id(n)), $(n.data))")

# AbstractTrees interface
AbstractTrees.children(n::Node) = isdefined(n, :children) ? n.children : typeof(n)[]
Base.eltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Node{I,T}
AbstractTrees.nodetype(::Type{Node{I,T}}) where {I,T} = Node{I,T}

# AbstractTrees traits
# Base.IteratorEltype(::Type{<:TreeIterator{Node{I,T}}}) where {I,T} = Base.HasEltype()
# AbstractTrees.parentlinks(::Type{Node{I,T}}) where {I,T} = nothing #AbstractTrees.StoredParents()

isbifurcating(n::Node) = all(map(x->isleaf(x) || degree(x) == 2, prewalk(n)))

"""
    postwalk(n)
Recursive postorder traversal.
"""
function postwalk(t)
    ns = typeof(t)[]
    function walk!(n)
        !isleaf(n) && for c in n.children walk!(c) end
        push!(ns, n)
    end
    walk!(t)
    ns
end

"""
    prewalk(n)
Recursive preorder traversal.
"""
function prewalk(t)
    ns = typeof(t)[]
    function walk!(n)
        push!(ns, n)
        !isleaf(n) && for c in n.children walk!(c) end
    end
    walk!(t)
    ns
end

"""
    getpath(n)
    getpath(n, m)
Get the path from n to the root, or the path connecting nodes n and m.
"""
function getpath(n::Node, m::Node)
    p1 = getpath(n)
    p2 = getpath(m)
    mrca = intersect(p1, p2)[1]
    i1 = findfirst(x->x==mrca, p1)
    i2 = findfirst(x->x==mrca, p2)
    p1[1:i1], p2[1:i2]
end

function getpath(n::Node)
    path = typeof(n)[]
    while !isnothing(n)
        push!(path, n)
        n = parent(n)
    end
    path
end

"""
    getdistance(n, m)
Get the distance from node n to m.
"""
function getdistance(n::Node, m::Node)
    p1, p2 = getpath(n, m)
    sum(distance.(p1[1:end-1])) + sum(distance.(p2[1:end-1]))
end

"""
    heights(n::Node)
Get node heights for the subtree rooted in `n`. The root has height 0!.
"""
function getheights(n::Node)
    d = Dict(id(n) => !isfinite(distance(n)) ? 0. : distance(n))
    function walk(n)
        if !haskey(d, id(n)) 
            x = isnan(distance(n)) ? 1. : distance(n)
            d[id(n)] = d[id(parent(n))] + x
        end
        for c in children(n)
            walk(c)
        end
    end
    walk(n)
    return d
end

height(n::Node) = distance.(getpath(n)) |> x->filter(!isnan, x) |> sum

# """
#     extract(n::Node, l::AbstractVector{String})
# Extract the tree with leaves in `l` from a given tree, preserving
# distances if relevant.
# """
# function extract(n::Node{I,T}, l::AbstractVector) where {I,T}
#     function walk(n)
#         if isleaf(n)
#             return name(n) ∈ l ? deepcopy(n) : nothing
#         else
#             below = Node{I,T}[]
#             for c in children(n)
#                 m = walk(c)
#                 !isnothing(m) && push!(below, m)
#             end
#             if length(below) == 0
#                 return nothing
#             elseif length(below) == 1
#                 setdistance!(below[1].data,
#                     distance(below[1]) + distance(n))
#                 return below[1]
#             else
#                 m = deepcopy(n)
#                 m.children = below
#                 for c in below c.parent = m end
#                 return m
#             end
#         end
#     end
#     walk(n)
# end


# """
#     extract(n::Node, l::AbstractVector{String})
# Extract the tree with leaves in `l` from a given tree, preserving
# distances if relevant.
# """
# function extract(n::Node{T,I}, l::Vector{F}) where {T,I,F<:AbstractString}
#     otr = deepcopy(n);
#     ls = Set(l);
#     pw = getleaves(otr);
#     lk = ReentrantLock();
#     Threads.@threads for xi in 1:length(pw)
#         x = pw[xi]
#         if isleaf(x) && !(name(x) ∈ ls)
#             px = parent(x)
#             lock(lk) do
#                 filter!(y -> y!=x, px.children)
#                 while length(px)==0
#                     filter!(y -> y!=px, parent(px).children)
#                     px = parent(px)
#                 end
#             end
#         end
#     end
#     pw = getleaves(otr)
#     Threads.@threads for x in pw
#         if length(parent(x))<2
#             p = parent(x)
#             s = distance(x)
#             torem = x
#             while length(p)<2
#                 s += distance(p)
#                 torem = p
#                 p = parent(p)
#             end
#             lock(lk) do
#                 filter!(y -> y!=torem, children(p))
#                 push!(p,x)
#                 x.data.distance = s
#             end
#         end
#     end
#     l1 = filter(x->!isleaf(x) && length(x)==1,postwalk(otr));
#     for x in l1
#         r = children(x)
#         while length(r)==1
#             r = children(r)[1]
#         end
#         p = parent(r)
#         s = distance(r)
#         torem = r
#         while length(p)<2
#             s += distance(p)
#             torem = p
#             p = parent(p)
#         end
#         filter!(y -> y!=torem, children(p))
#         push!(p,r)
#         r.data.distance = s
#     end
#     return otr
# end


"""
    insertnode!(n, m)
Insert node `m` between `n` and its parent at a distance `distance(n) -
distance(m)` from node `n`.
"""
function insertnode!(n::Node{I,T}, m::Node{I,T}) where {I,T}
    if !isroot(n)
        a = parent(n)
        delete!(a, n)
        push!(a, m)
    end
    push!(m, n)
    setdistance!(n.data, distance(n) - distance(m))
end

"""
    insertnode!(n, dist=NaN, name="")
Insert a new node with name `name` above `n` at a distance `dist` from `n`.
"""
function insertnode!(n::Node{I,<:NewickData}; dist=NaN, name="") where I
    i = maximum(id, postwalk(getroot(n))) + 1
    dist = isnan(dist) ? distance(n) / 2 : dist
    m = Node(I(i), NewickData(d=dist, n=name))
    insertnode!(n, m)
    return m
end

"""
    getlca(n, taxona, taxonb)
Get the last common ancestor of taxona and taxonb in the tree rooted in node n.
"""
getlca(n::Node, a) = getlca(n, a, a)
getlca(n::Node, a, b) = getlca(findleaf(n, a), findleaf(n, b))

function getlca(a::Node, b::Node)
    while !(b ∈ getleaves(a))
        a = parent(a)
    end
    return a
end

"""
    getleaves(n)
Get the leaves of the tree rooted in node n.
"""
function getleaves(n::N) where N<:Node  # mostly faster than Leaves...
    xs = N[]
    for node in postwalk(n)
        isleaf(node) && push!(xs, node)
    end
    xs
end


"""
    set_outgroup!(node::Node)
Set the node `node` as an outgroup to the tree to which it belongs. This will
introduce a root node between `node` and its parent. When the tree has a
trifurcating root, this will resolve the trifurcation (i.e. root the tree).
When the tree is already bifuracting, this will reroot the tree. (perhaps it
would be wise to split these two functionalities, since they do different
things... i.e. rooting/rerooting).
```julia
julia> tr = nw"(A,(B,(C,D)),E);"
(A,(B,(C,D)),E);
julia> NewickTree.set_outgroup!(tr[2][1])
(((C,D),(A,E)),B);
```
"""
function set_outgroup!(node::Node{I}) where I
    curr_root = getroot(node)
    # someone tried to set the root as outgroup, return just the root
    # (or should we error? or should we proceed, returning `(node);`)
    curr_root == node && return node
    p = parent(node)
    # node is already an outgroup
    p == curr_root && length(children(curr_root)) == 2 && return p
    i = maximum(id, prewalk(curr_root)) + 1
    newroot = Node(I(i), NewickData())    
    # nodes above where the root will be inserted, this is the part of the
    # tree that has to be reoriented
    alongtheway = [node]
    while !isnothing(p)
        push!(alongtheway, p) 
        p = parent(p)
    end
    # reorient the part of the tree above the new root node
    while length(alongtheway) > 1
        p = pop!(alongtheway)
        # this will become the parent of p
        a = filter(x->id(x) == id(alongtheway[end]), children(p))[1]  
        d = distance(a)
        delete!(p, id(a))
        if a != node
            push!(a, p)
            setdistance!(p, d)
        else
            push!(newroot, p)
            setdistance!(p, d/2)
        end
    end 
    push!(newroot, node)
    setdistance!(node, distance(node)/2)
    if length(children(curr_root)) == 1  
        # when the previous root was not multifurcating, we end up with a
        # superfluous node, so we delete it
        p = parent(curr_root)
        delete!(p, id(curr_root))
        push!(p, curr_root[1])
    end
    return newroot
end

set_outgroup(node) = set_outgroup!(deepcopy(node))
set_outgroup(tree, leafname)  = set_outgroup( getlca(tree, leafname)) 
set_outgroup!(tree, leafname) = set_outgroup!(getlca(tree, leafname)) 

"""
    prune!(tree, leaves)
Prune `leaves` from the tree, i.e. get the tree structure for the
complement of `leaves`. This is actually the dual of `extract`.
"""
function prune!(n, leaves)
    for x in leaves
        x ∈ getleaves(n) && recursive_prune!(x)
    end
    return n
end

function recursive_prune!(node)
    degree(node) == 0 && isroot(node) && return node
    if degree(node) == 0
        delete!(parent(node), node)
        pnode = parent(node)
        node.parent = node
        recursive_prune!(pnode)
    elseif degree(node) == 1 && !isroot(node)
        p = parent(node)
        c = node[1]
        setdistance!(c, distance(c) + distance(node))
        node.parent = node
        delete!(p, node)
        push!(p, c)
        c.parent = p
    elseif degree(node) == 1 && isroot(node)
        setdistance!(node, distance(node) + distance(node[1]))
        node.children = children(node[1])
        for c in children(node)
            c.parent = node
        end
    end
end

relabel(tree, m) = relabel!(deepcopy(tree), m)
function relabel!(tree, m::Dict{T,T}) where T<:AbstractString
    for n in getleaves(tree)
        n.data.name = m[name(n)]
    end
    return tree
end


"""
    fixdist!(tr::Node{I,T}) where {I,T}

fixes negative distances for fastme generated trees
"""
function fixdist!(tr::Node{I,T}) where {I,T}
    for n in postwalk(tr) # fix negative lengths
        dn = distance(n)
        if !isnan(dn)
            n.data.distance = abs(dn)
        end
    end
end

"""
    isbinarytree(pw::Vector{Node{I,T}}) where {I,T}

checks that the tree is a binary bifurcating tree
"""
function isbinarytree(pw::Vector{Node{I,T}}) where {I,T}
     all(unique([isleaf(x) ? 0 : length(x) for x in pw]) .∈ Ref(Set([0,2])))
end

"""
    rmleaf(n::Node, x::Node)

removes the leaf node x from tree n, returning a new tree with the removed leaf
"""
function rmleaf(n::Node, x::Node)
    p = parent(x)
    s = sister(x)
    if !isroot(p)
        pp = parent(p)
        # connect pp to s
        push!(pp.children, s)
        deleteat!(pp.children, findfirst(c->c==p,children(pp)))
        s.parent = pp
        s.data.distance = distance(s) + distance(p)
        return n
    else
        # while length(parent(s))>0 pop!(s.parent) end
        # @assert isroot(s) "s is not the root"
        return make_independent_tree(s)
    end
end


"""
    extract(n::Node, l::Vector{T}) where T<:AbstractString

extracts a new tree from tree n that only contains leaf names in l
"""
function extract(n::Node, l::Vector{T}) where T<:AbstractString
    pow = postwalk(n);
    @assert isbinarytree(pow) "not a binary tree. fix."
    torem = setdiff(leafnames(n), l)
    for k in torem
        x = filter(x->name(x)==k, getleaves(n))
        if length(x)>0
            n = rmleaf(n, x[1])
        else
            @error "couldnt find leaf $k"
        end
    end
    return n
end

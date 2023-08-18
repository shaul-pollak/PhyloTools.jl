"""
Ancestral reconstruction of discrete (integer) traits based on Narushima and Hanazawa (1997).

The vector x takes integer values and should be ordered according to the leaves on the tree.

This is a copy of the MPR function from the R package "ape"
"""
function mpr(tr::Node, x::Vector{T})::Matrix{T} where {T<:Integer}

    if length(children(tr)) != 2
        error("tree must be rooted")
    end

    n = ntip(tr)
    m = length(postwalk(tr)) - n # number of internal nodes
    edges = edgemat(tr; postorder=true)

    I = zeros(Union{T,Missing}, n + m)
    I[1:n] .= x
    I[n+1:n+m] .= missing
    I = [I'; I']

    function med(x)
        i = length(x) / 2 |> Int
        return sort(x[:])[[i, i + 1]]
    end

    ## 1st pass
    anc = edges[1:2:end, 1]
    des = reshape(edges[:, 2], 2, :)
    for i in eachindex(anc)
        I[:, anc[i]] = med(I[:, des[:, i]])
    end

    ## 2nd pass
    out = zeros(Union{T,Missing}, 2, m)
    out[:] .= missing
    out = [[x';x'] out]
    # we do the most basal node before lopping
    Iw = I[:, des[:, m]][:]
    out[:, anc[m]] .= extrema(med(Iw))
    for i in m-1:-1:1
        j = anc[i]
        Iw = I[:, des[:,i]][:]
        k = findfirst(edges[:,2] .== j)
        tmp = out[:, edges[k, 1]]
        out[1, j] = minimum(med([tmp[1]; tmp[1]; Iw]))
        out[2, j] = maximum(med([tmp[2]; tmp[2]; Iw]))
    end

    return out'

end

"""
This function returns the edge matrix as in the APE package in R.
The edges are preorder ordered (first edge is the root)
The entries of the matrix represent node ids.
"""
function edgemat(n::Node{T,I}; postorder::Bool=false) where {T,I}
    pw = prewalk(n)
    ne = [length(children(x)) for x in pw] |> sum
    edges = zeros(T, ne, 2)
    i = 0
    for x in pw
        v1 = id(x)
        for c in children(x)
            i += 1
            v2 = id(c)
            @inbounds edges[i, 1] = v1
            @inbounds edges[i, 2] = v2
        end
    end
    return postorder ? reverse(edges, dims=1) : edges
end
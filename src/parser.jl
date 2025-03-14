# Parsing and I/O
"""
A string macro to read a newick string to a tree. Just to have the
syntactic sugar enabling `nw"((A,B),C);"`
"""
macro nw_str(s)
    readnw(s)
end

"""
    readnw(s::AbstractString, I::Type)
Read a newick string to a tree. Supports the original Newick standard
(http://evolution.genetics.washington.edu/phylip/newicktree.html). One can
have either support values for internal nodes or a node label, but not both.
"""
readnw(s::AbstractString, I::Type=UInt32) =
    try
        if (length(s) <= 100) && ispath(s) && isfile(s)
            l = readlines(s) |> join
            l = replace(l, " " => "")
            if l[end] != ';'
                l = "$l;"
            end
            tr = readnw(IOBuffer(l), I)
            fix_node_ids!(tr)
            return tr
        else
            tr = readnw(IOBuffer(s), I)
            fix_node_ids!(tr)
            return tr
        end
    catch EOFError
        more = s[end] != ";" ? "(no trailing semicolon?)" : ""
        more = !ispath(s) ? more : "(path-like arg instead of Newick string?)"
        throw("Malformed Newick string '$s' $more")
    end

"""
    nwstr(n::Node{I,N}; internal=false)
Generate a newick tree string for the tree rooted in `n`. To make this work,
`N` (the type of `n.data`) should implement `name()` and `distance()` methods,
and optionally `support()`. If `support` is implemented for the data type and
it is not NaN, it will supersede `name` for the labeling of internal nodes. See
for instance the `NewickData` type.
"""
function nwstr(n; internal=false, dist=true)
    function walk(n)
        d = dist ? stringify(':', distance(n)) : ""
        isleaf(n) && return "$(name(n))$d"
        s = join([walk(c) for c in children(n)], ",")
        sv = hasmethod(support, Tuple{typeof(n)}) ?
             stringify(support(n)) : ""
        sv = sv == "" && internal ? name(n) : sv
        d = dist ? stringify(':', distance(n)) : ""
        return "($s)$sv$d"
    end
    s = walk(n)
    s * ";"
end

stringify(x) = isnan(x) ? "" : string(x)
stringify(x, y) = isnan(y) ? "" : string(x, y)

"""
    writenw(io::IO, n)
    writenw(fname::AbstractString, n)
Write a newick representation of the tree rooted in `n`. Note that the tree
data type should allow `nwstr(n)` to work. See the `nwstr` docstring.
"""
writenw(io::IO, n) = write(io, nwstr(n))
writenw(fname::AbstractString, n) = write(fname, nwstr(n) * "\n")


# function readnw(io::IOBuffer, I::Type=UInt32)

#     i = I(1)
#     c = read(io, Char)
#     stack = Node[]
#     currdata = NewickData()
#     nodedata = NewickData()
#     nums = reduce(vcat, string.(0:9) .|> collect)
#     buff = IOBuffer()

#     while (!eof(io)) || (c != ';')

#         if c == '('

#             push!(stack, Node(i, NewickData()))
#             i += one(i)
#             c = read(io, Char)

#         elseif c == ')' || c == ','

#             target = pop!(stack)
#             source = last(stack)
#             push!(source, target)

#             if name(nodedata) != ""

#                 target.data = deepcopy(nodedata)
#                 nodedata.name = ""

#             else

#                 target.data = deepcopy(currdata)
#                 # currdata = NewickData()

#             end

#             if c == ')'

#                 c = read(io, Char)

#                 if eof(io) || c == ';'

#                     break

#                 else

#                     nodename = ""
#                     sv = nothing

#                     if (c == '\'') || !(c ∈ nums)

#                         nodename, c = get_nodename(io, buff, c)
#                         s1 = split(nodename, ':')

#                         if length(s1)>1
#                             sv = string(s1[1])
#                             nodename = string(s1[2])
#                         else
#                             nodename = string(s1[1])
#                         end

#                     end

#                     if eof(io) || (c==';')
#                         last(stack).data = NewickData(0.0, 100.0, nodename)
#                         break
#                     end

#                     nodedata, c = get_nodedata(io, c, nodename; support=sv, buff=buff)

#                 end

#             else

#                 c = read(io, Char)

#             end

#         elseif isspace(c)

#             c = read(io, Char)

#         else

#             push!(stack, Node(i, NewickData()))
#             i += one(i)
#             leafname, c = get_leafname(io, c, buff)
#             currdata, c = get_nodedata(io, c, leafname; buff=buff)

#         end

#     end

#     last(stack)

# end

function get_leafname(io::IOBuffer, c, buff::IOBuffer)

    while !_isnwdelim(c)

        write(buff, c)
        c = read(io, Char)

    end

    strip(take!(buff) |> String) |> String, c

end




function get_nodedata(io::IOBuffer, c, name=""; support::Union{String,Nothing}=nothing, buff::IOBuffer)

    # get everything up to the next comma or )
    if isnothing(support)

        support, c = _readwhile!(io, c, buff)

    end

    distance = ""

    if c == ':'

        c = read(io, Char)
        distance, c = _readwhile!(io, c, buff)

    end

    sv = nanparse(support)

    if typeof(sv) == String
        name = String(strip(sv))
        sv = NaN
    end

    dv = nanparse(distance)
    NewickData(dv, sv, name), c

end






_isquote(c::Char) = c == '\'' || c == '"';




function _readwhile!(io::IOBuffer, c, buff)
    if _isquote(c)
        write(buff, c)
        c = read(io, Char)
        while !_isquote(c)
            write(buff, c)
            c = read(io, Char)
        end
    end
    while !_isnwdelim(c)
        write(buff, c)
        c = read(io, Char)
    end
    String(take!(buff)), c
end



_isnwdelim(c::Char) = c == ',' || c == ')' || c == ':' || c == ';'

function nanparse(x)
    y = tryparse(Float64, x)
    isnothing(y) ? (x == "" ? NaN : x) : parse(Float64, x)
end

function quotednanparse(x)
    x = strip(x, ['\'', '"'])
    support = Float64(0)
    name = ""
    if contains(x, ':')
        xs = split(x, ':')
        support = nanparse(xs[1])
        name = string(xs[2])
    end
    return support, name
end



function get_nodename(io, buff, c)

    c != '\'' && write(buff,c)
    copyuntil(buff, io, '\'')

    cret = eof(io) ? '\0' : read(io, Char)

    return String(take!(buff)), cret

end


function readnw(io::IOBuffer, I::Type=UInt32)

    i = I(1)
    c = read(io, Char)
    stack = []
    currdata = NewickData()
    buff = IOBuffer()

    while c != ';'

        if c == '('

            push!(stack, Node(i, NewickData())); i += one(i)
            c = read(io, Char)

        elseif c == ')' || c == ','

            target = pop!(stack)
            source = last(stack)
            push!(source, target)
            target.data = currdata

            if c == ')'
                c = read(io, Char)
                eof(io) || c == ';' ? (break) :
                    (currdata, c) = get_nodedata(io, c; buff=buff)

            else

                c = read(io, Char)
            end

        elseif isspace(c)

            c = read(io, Char)

        else

            push!(stack, Node(i, NewickData())); i += one(i)
            leafname, c = get_leafname(io, c, buff)
            currdata, c = get_nodedata(io, c, leafname; buff=buff)

        end

    end

    last(stack)
end

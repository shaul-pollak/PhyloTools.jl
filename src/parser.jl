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
        if (length(s)<=100) && ispath(s) && isfile(s)
            l = readlines(s) |> join
            l = replace(l, " " => "")
            if l[end]!=';'
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
    s*";"
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


function readnw(io::IOBuffer, I::Type=UInt32)
    i = I(1)
    c = read(io, Char)
    stack = []
    currdata = NewickData()
    nodedata = NewickData()
    while c != ';'
        if c == '('
            push!(stack, Node(i, NewickData())); i += one(i)
            c = read(io, Char)
        elseif c == ')' || c == ','
            target = pop!(stack)
            source = last(stack)
            push!(source, target)
            if nodedata.name != ""
                target.data = nodedata
                # target.id = parse(I, target.data.name)
                nodedata = NewickData()
            else
                target.data = currdata
            end
            if c == ')'
                c = read(io, Char)
                if eof(io) || c == ';'
                    break
                else
                    if c=='\''
                        nodename, c = PhyloTools.get_leafname(io, c)
                        sv, nodename = split(nodename,':') .|> string
                    else
                        nodename = ""
                        sv = nothing
                    end
                    nodedata, c = PhyloTools.get_nodedata(io, c, nodename, support=sv)
                end
            else
                c = read(io, Char)
            end
        elseif isspace(c)
            c = read(io, Char)
        else
            push!(stack, Node(i, NewickData())); i += one(i)
            leafname, c = PhyloTools.get_leafname(io, c)
            currdata, c = PhyloTools.get_nodedata(io, c, leafname)
        end
    end
    last(stack)
end

function get_leafname(io::IOBuffer, c)
    leafname = ""
    if c!='\''
        while !_isnwdelim(c)
            leafname *= c
            c = read(io, Char)
        end
    else
        c = read(io, Char)
        while c!='\''
            leafname *= c
            c = read(io, Char)
        end
        c = read(io, Char)
    end
    String(strip(leafname)), c
end

function get_nodedata(io::IOBuffer, c, name=""; support::Union{String,Nothing}=nothing)
    # get everything up to the next comma or )
    if isnothing(support)
        support, c = _readwhile!(io, c)
    end
    distance = ""
    if c == ':'
        c = read(io, Char)
        distance, c = _readwhile!(io, c)
    end
    if length(support)>0 && support[1]=='\''
        sv,name = quotednanparse(support)
    else
        sv = nanparse(support)
        if typeof(sv) == String
            name = String(strip(sv))
            sv = NaN
        end
    end
    NewickData(nanparse(distance), sv, name), c
end

_isquote(c::Char) = c == '\'' || c == '"';

function _readwhile!(io::IOBuffer, c)
    out = ""
    if _isquote(c)
        out *= c;
        c = read(io,Char);
        while !_isquote(c)
            out *= c;
            c = read(io,Char)
        end
    end
    while !_isnwdelim(c)
        out *= c
        c = read(io, Char)
    end
    out, c
end



_isnwdelim(c::Char) = c == ',' || c == ')' || c == ':' || c == ';'

function nanparse(x)
    y = tryparse(Float64, x)
    isnothing(y) ? (x == "" ? NaN : x) : parse(Float64, x)
end

function quotednanparse(x)
    x = strip(x,['\'','"']);
    support = Float64(0);
    name = "";
    if contains(x,':')
        xs = split(x,':');
        support = nanparse(xs[1]);
        name = string(xs[2]);
    end
    return support, name
end

module PhyloTools

using AbstractTrees
using Colors
using MakieCore
using GeometryBasics
using ProgressMeter
using Graphs
using Mmap
using IterTools
using CodecZstd, CodecZlib
using StringViews

# export Node, NewickData
export isroot, isleaf, postwalk, prewalk, children, sister, getleaves
export readnw, writenw
export distance, name, id, getheights
export leafnames, ntip, cophenetic, make_independent_tree, reroot!, midpoint_root, mad
export MRCA, ladderize!
export readfasta, read_gff, readclu, treepositions
export mpr, prune

include("node.jl")
include("parser.jl")
include("functions.jl")
include("treemi.jl")
include("makie.jl")
include("mad.jl")
include("ancrec.jl")


end

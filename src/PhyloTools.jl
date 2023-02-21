module PhyloTools

using AbstractTrees
using RecipesBase
using Colors
using MakieCore
using GeometryBasics

export Node, NewickData
export isroot, isleaf, postwalk, prewalk, children, sister
export getroot, getlca, getleaves
export insertnode!, print_tree, readnw, writenw, @nw_str
export distance, name, id, nwstr, degree, getheights
export set_outgroup!, set_outgroup
export leafnames, Ntip, cophenetic, make_independent_tree, reroot!, midpoint_root
export tippoint!, nodepoint!, cluster_tree, extract_with_tips, gene
export treeVI!, treeVI, duplication_score, whichmin, whichmax
export MRCA, ladderize!
export readfasta, read_gff

include("node.jl")
include("parser.jl")
include("functions.jl")
include("plotting.jl")
include("treemi.jl")
include("makie.jl")

end

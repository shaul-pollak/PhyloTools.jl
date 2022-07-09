# Phy
A julia package for phylogenetics that includes some visualization functions. \
The basic tree data structure and the plot recipes are based off of NewickTree.jl

## Tree manipulation
To read a tree:
```Julia
tree = readnw(readline(path_to_tree))
```

We can also root the tree automatically based on a midpoint rooting criterion
(at the node which has the minimal distance to the tip that is farthest away
from it):\
tree_rooted = root_midpoint(tree)

Trees can be clustered based on a divergence threshold:\
clu = cluster_tree(tree_rooted,2)

## Plotting
Plotting is based on Plots.jl. \
To plot a tree (without showing its tip labels): \
plot(tree_rooted)

or with tip labels: \
plot(tree_rooted; labels=true)

To plot the tree in an upward direction:\
plot(tree_rooted; upwards=true)

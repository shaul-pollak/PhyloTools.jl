# PhyloTools
A julia package for phylogenetics that includes some visualization functions. \
The basic tree data structure and the plot recipes are based off of NewickTree.jl

## Tree manipulation
To read a tree:
```Julia
tree = readnw(path_to_tree)
```

We can also root the tree automatically based on a midpoint rooting criterion
(at the node which has the minimal distance to the tip that is farthest away
from it):
```Julia
tree_rooted = root_midpoint(tree)
```

Trees can be clustered based on the fraction of the total divergence allowed in
each subtree:
```Julia
clu = cluster_tree(tree_rooted, 2)
```

## Plotting
Plotting is based on Makie.jl. \
To plot a tree (without showing its tip labels): 
```Julia
plot(tree_rooted)
```

or with tip labels: 
```Julia
plot(tree_rooted; labels=true)
```

To plot the tree in an upward direction:
```Julia
plot(tree_rooted; upwards=true)
```

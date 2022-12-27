# local_node_connectivity_with_GPGPU

about local_node_connectivity is [here](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.approximation.connectivity.local_node_connectivity.html#networkx.algorithms.approximation.connectivity.local_node_connectivity).

## Usage

### generate graph
graph generator is based on the [GNP random graph](https://networkx.org/documentation/stable/reference/generated/networkx.generators.random_graphs.gnp_random_graph.html#networkx.generators.random_graphs.gnp_random_graph).

```bash
python gengraph (# of node) (rate of edge) (seed)
```

### build

```bash
make
```

### test

```bash
timecheck.sh
```

## Credits
https://junstar92.tistory.com/272
https://networkx.org/documentation/stable/reference/generated/networkx.generators.random_graphs.gnp_random_graph.html#networkx.generators.random_graphs.gnp_random_graph

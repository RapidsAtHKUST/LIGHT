# UnlabeldSubgraphEnumeration
## Licenses
The code package can be ONLY used for your personal research purpose. Please DO NOT distribute the code to anyone else, as our research on set intersections is under going.

## Compile
Under the src directory, execute the following commands to compile the source code. In our experiments, we use icpc 16.0.0 to compile the source code. If you do not install icpc, then replace line 3 with "cmake ..". After compiling, You will get two binaries "LIGHT" and "Converter".

```zsh
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=path_to_icpc ..
make -j
```

## Preprocess Dataset
Use "Converter" to preprocess the dataset. "Converter" takes a graph file stored in the edge list as input and converts it to the required format. The vertex id is rangled from 0 to |V| - 1 where V is the vertex set. You can find a sample input file under the sample_graph folder.

Example (triangle):

```zsh
0 1
0 2
1 2
```

Following TwinTwig and SEED, we preprocess the graph in two aspects.

* Decompose the graph and only keep 2-core of the graph.
* Reorder the vertices in the graph with the ascending order of the degree.

Execute the binary with the following command ./Converter -o Convert -ift EdgeList -ief input_edge_list -oft CSR -odf output_degree_file -oaf output_adjacent_list_file in which -ief specifies the input edge list, -odf specifies the degree output file and -oaf specificies the adjacent output file. The -o, -ift and -oft configure the operation, the input file format and the output file format respectively. Please do not change them.

Example:

```zsh
./Converter -o Convert -ift EdgeList -ief ../../sample_graph/triangle.edge_list -oft CSR -odf ../../sample_graph/triangle_degree.bin -oaf ../../sample_graph/triangle_adj.bin
```

## Execute

Execute the binary with the following command ./LIGHT -df input_degree_file -af input_adjacent_list_file -p pattern_name -n num_threads. Currently, our program supports 9 patterns, which are
named as p0-8. Specifically, p1-7 are listed in the paper. p0 is triangle, and p8 is a combination of three triangles. In our experiment, we execute the binary with 64 threads on a machine equipped with 20 physical cores that enable hyper-threading.

Example:

```zsh
./LIGHT -df ../../sample_graph/triangle_degree.bin -af ../../sample_graph/triangle_adj.bin -p p0 -n 1
```

The python script "generate_query_plan.py" is used to generate the query plan.

## Optimization

Following TwinTwig, SEED, CRYSTAL and BigJoin, we adopt two optimizations.

* Use symmetry breaking to reduce the vertices involved in the computation. Specifically, if the symmetry breaking requires that v1 < v2, then we can only consider the neighbors of v1 the id of which is greater than v1 instead of the entire neighbor set.
* The candidate sets of the vertices that have no constraint on each other are kept as lists instead of performing the Cartesian products.

Except DUALSIM that we only have the binary, all competing algorithms including both algorithms working in a single machine (e.g. SE) and distributed algorithms adopt the two optimization to make a fair comparison in our experiment. And they take the graph after the preprocessing as the input.

If you want to examine the performance of enumerating each result by performing Cartesian products, then you can add "#define ENUMERATE_RESULTS" in Config.h to disable the second optimization. Under this circumstance, please make sure each competing algorithm disable the optimization to make a fair comparison. Note that performing the Cartesian products can take much longer time, because there can be a huge number of results (e.g., 9.56e+12 subgraphs identical to p6 in livejournal). 
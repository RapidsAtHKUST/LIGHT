//
// Created by ssunah on 11/28/16.
//

#include "Config.h"
#include "Graph.h"
#include "Utility.h"
#include <chrono>

void Graph::DestroyGraph(Graph *graph) {
    delete graph;
}

void Graph::ComputeNeighborPlus() {
    edges_plus_offset_ = new unsigned int[vertices_num_];

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        unsigned int temp_count;
        const unsigned int* neighbors = Neighbors(i, temp_count);
        int begin = 0;
        int end = (int)temp_count - 1;

        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] > i) {
                end = mid - 1;
            }
            else {
                begin = mid + 1;
            }
        }

        edges_plus_offset_[i] = (unsigned int)begin;
    }
}

void Graph::ConvertEdgeListToCSR(vector<pair<unsigned int, unsigned int>> &edge_list) {
    // Compute the degree of vertices.
    map<unsigned int, unsigned int> vertex_degree_dict;
    for (auto edge : edge_list) {
        if (vertex_degree_dict.find(edge.first) == vertex_degree_dict.end()) {
            vertex_degree_dict.insert(make_pair(edge.first, 0));
        }
        if (vertex_degree_dict.find(edge.second) == vertex_degree_dict.end()) {
            vertex_degree_dict.insert(make_pair(edge.second, 0));
        }

        vertex_degree_dict[edge.first] += 1;
        vertex_degree_dict[edge.second] += 1;
    }

    vertices_num_ = vertex_degree_dict.size();
    edges_offset_ = new unsigned int[vertices_num_ + 1];

    // Compute the offset of vertex and the max degree.

    max_degree_ = 0;

    edges_offset_[0] = 0;
    for (auto vertex_degree : vertex_degree_dict) {
        edges_offset_[vertex_degree.first + 1] = edges_offset_[vertex_degree.first] + vertex_degree.second;

        if (vertex_degree.second > max_degree_)
            max_degree_ = vertex_degree.second;
    }

    // PatternGraphRelease the memory of vertex_degree_dict.
    map<unsigned int, unsigned int>().swap(vertex_degree_dict);

    edges_num_ = edges_offset_[vertices_num_];
    edges_ = new unsigned int[edges_num_];

    // Store edges.
    for (auto edge : edge_list) {
        edges_[edges_offset_[edge.first]++] = edge.second;
        edges_[edges_offset_[edge.second]++] = edge.first;
    }

    // Recompute the offset of vertex.
    for (unsigned int i = vertices_num_; i > 0; --i) {
        edges_offset_[i] = edges_offset_[i - 1];
    }
    edges_offset_[0] = 0;

    // Sort the neighbors.
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        sort(edges_ + edges_offset_[i], edges_ + edges_offset_[i + 1],
             [](const unsigned int left, const unsigned int right) -> bool {
                return left < right; }
        );
    }

    // Compute the neighbor plus.
    ComputeNeighborPlus();
    ComputeMetaData();
}

void Graph::RemoveNonCoreVertex() {
    cout << PRINT_SEPARATOR << endl;
    cout << "Remove non core vertices start..." << endl;
    auto start = chrono::high_resolution_clock::now();

    unsigned int* core_table = new unsigned int[vertices_num_];
    Utility::GetKCore(this, core_table);

    unsigned int* updated_id = new unsigned int[vertices_num_];
    bool* valid_id = new bool[vertices_num_];
    memset(valid_id, 0, sizeof(bool) * vertices_num_);

    unsigned int core_vertex_count = 0;
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        if (core_table[i] >= 2) {
            valid_id[i] = true;
            updated_id[i] = core_vertex_count;
            core_vertex_count += 1;
        }
    }

    delete[] core_table;

    unsigned int* core_edge_offset = new unsigned int[core_vertex_count + 1];

    unsigned int core_edge_count = 0;
    core_edge_offset[0] = 0;

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        if (!valid_id[i])
            continue;

        unsigned int src_new_id = updated_id[i];
        unsigned int neighbors_count;
        const unsigned int *neighbors = Neighbors(i, neighbors_count);

        for (unsigned int j = 0; j < neighbors_count; ++j) {
            if (!valid_id[neighbors[j]])
                continue;

            core_edge_count += 1;
            core_edge_offset[src_new_id + 1] = core_edge_count;
        }
    }

    unsigned int* core_edges = new unsigned int[core_edge_count];
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        if (!valid_id[i])
            continue;

        unsigned int src_new_id = updated_id[i];
        unsigned int neighbors_count;
        const unsigned int* neighbors = Neighbors(i, neighbors_count);

        for (unsigned int j = 0; j < neighbors_count; ++j) {
            if (!valid_id[neighbors[j]])
                continue;

            unsigned int dst_new_id = updated_id[neighbors[j]];
            core_edges[core_edge_offset[src_new_id]++] = dst_new_id;

        }

    }

    for (unsigned int i = core_vertex_count; i > 0; --i) {
        core_edge_offset[i] = core_edge_offset[i - 1];
    }
    core_edge_offset[0] = 0;


    delete[] updated_id;
    delete[] valid_id;
    delete[] edges_offset_;
    delete[] edges_;

    vertices_num_ = core_vertex_count;
    edges_num_ = core_edge_count;

    edges_offset_ = core_edge_offset;
    edges_ = core_edges;

    max_degree_ = 0;

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        if (edges_offset_[i + 1] - edges_offset_[i] > max_degree_)
            max_degree_ = edges_offset_[i + 1] - edges_offset_[i];
    }

    auto end = chrono::high_resolution_clock::now();
    cout << "Remove non core vertices time:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Remove non core vertices end." << endl;
}

void Graph::Relabel() {
    cout << PRINT_SEPARATOR << endl;
    cout << "Relabel start..." << endl;
    auto start = chrono::high_resolution_clock::now();

    unsigned int* labels = new unsigned int[vertices_num_];
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        labels[i] = i;
    }

    unsigned int* offset = edges_offset_;
    sort(labels, labels + vertices_num_, [offset](const unsigned int left, const unsigned int right) -> bool {
        unsigned int left_degree = offset[left + 1] - offset[left];
        unsigned int right_degree = offset[right + 1] - offset[right];

        if (left_degree < right_degree) {
            return true;
        }
        else if (left_degree == right_degree) {
            return left < right;
        }
        else {
            return false;
        }
    });

    unsigned int* map_id_label = new unsigned int[vertices_num_];
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        map_id_label[labels[i]] = i;
    }

    delete[] labels;

    unsigned int* label_offset = new unsigned int[vertices_num_ + 1];
    unsigned int* label_edges = new unsigned int[edges_num_];

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        unsigned int degree = Degree(i);
        unsigned int label = map_id_label[i];
        label_offset[label + 1] = degree;
    }

    label_offset[0] = 0;
    for (unsigned int i = 1; i <= vertices_num_; ++i) {
        label_offset[i] += label_offset[i - 1];
    }

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        unsigned int neighbors_count;
        const unsigned int* neighbors = Neighbors(i, neighbors_count);

        unsigned int src_label = map_id_label[i];
        for (unsigned int j = 0; j < neighbors_count; ++j) {
            unsigned int dst_label = map_id_label[neighbors[j]];
            label_edges[label_offset[src_label]++] = dst_label;
        }
    }

    for (unsigned int i = vertices_num_; i > 0; --i) {
        label_offset[i] = label_offset[i - 1];
    }
    label_offset[0] = 0;

    delete[] edges_offset_;
    delete[] edges_;
    delete[] map_id_label;

    edges_offset_ = label_offset;
    edges_ = label_edges;

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        sort(edges_ + edges_offset_[i], edges_ + edges_offset_[i + 1],
             [](const unsigned int left, const unsigned int right) -> bool {
                 return left < right; }
        );
    }

    auto end = chrono::high_resolution_clock::now();
    cout << "Relabel time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Relabel end." << endl;
}

Graph *Graph::LoadEdgeList(const string &file_path) {
    cout << PRINT_SEPARATOR << endl;
    cout << "Load edge list start..." << endl;

    auto start = chrono::high_resolution_clock::now();

    std::ifstream infile(file_path);

    if (infile.is_open()) {
        cout << "Open edge list file " << file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot open edge list file " << file_path << " ." << endl;
        exit(-1);
    }

    vector<pair<unsigned int, unsigned int> > edge_list;

    // Reserve memory.
    edge_list.reserve(256 * 1024 * 1024);

    string line;
    unsigned int line_count = 0;

    while (getline(infile, line))
    {
        line_count += 1;

        // Skip the line starting with '#'.
        if (line[0] == '#') {
            continue;
        }

        istringstream iss(line);
        int begin_vertex, end_vertex;
        if (!(iss >> begin_vertex >> end_vertex)) {
            cerr << "Line " << line_count << " error: cannot convert the vertex id to integer." << endl;
            exit(-1);
        }

        if (begin_vertex >= end_vertex) {
            cerr << "Line " << line_count << " error: the id of the end vertex is greater than or equal to the id of the begin vertex." << endl;
            exit(-1);
        }

        edge_list.push_back(make_pair(begin_vertex, end_vertex));
    }

    // Sort the edge list based on the begin vertex.
    sort(edge_list.begin(), edge_list.end(),  [](const pair<unsigned int, unsigned int> left, const pair<unsigned int, unsigned int> right) -> bool {
        return left.first < right.first;
    });

    Graph* graph = new Graph();
    graph->ConvertEdgeListToCSR(edge_list);
    auto end = chrono::high_resolution_clock::now();

    cout << "Load edge list time:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Load edge list end." << endl;
    return graph;
}

Graph *Graph::LoadCSR(const string &degree_file_path, const string &edge_file_path) {
    cout << PRINT_SEPARATOR << endl;
    cout << "Load CSR start..." << endl;
    auto start = chrono::high_resolution_clock::now();

    ifstream deg_file(degree_file_path, ios::binary);

    if (deg_file.is_open()) {
        cout << "Open degree file " << degree_file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot open degree file " << degree_file_path << " ." << endl;
        exit(-1);
    }

    int int_size;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);

    Graph* graph = new Graph();

    deg_file.read(reinterpret_cast<char *>(&graph->vertices_num_), 4);
    deg_file.read(reinterpret_cast<char *>(&graph->edges_num_), 4);

    graph->edges_offset_ = new unsigned int[graph->vertices_num_ + 1];
    unsigned int* degrees = new unsigned int[graph->vertices_num_];

    deg_file.read(reinterpret_cast<char *>(degrees), sizeof(int) * graph->vertices_num_);


    deg_file.close();
    deg_file.clear();

    ifstream adj_file(edge_file_path, ios::binary);

    if (adj_file.is_open()) {
        cout << "Open edge file " << edge_file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot open edge file " << edge_file_path << " ." << endl;
        exit(-1);
    }

    graph->edges_ = new unsigned int[graph->edges_num_];

    graph->edges_offset_[0] = 0;
    for (unsigned int i = 1; i <= graph->vertices_num_; ++i) {
        graph->edges_offset_[i] = graph->edges_offset_[i - 1] + degrees[i - 1];
    }

    graph->max_degree_ = 0;

    for (unsigned int i = 0; i < graph->vertices_num_; ++i) {
        if (degrees[i] > 0) {
            if (degrees[i] > graph->max_degree_)
                graph->max_degree_ = degrees[i];

            adj_file.read(reinterpret_cast<char *>(graph->edges_ + graph->edges_offset_[i]), degrees[i] * sizeof(int));
        }
    }

    adj_file.close();
    adj_file.clear();

    graph->ComputeNeighborPlus();
    graph->ComputeMetaData();
    delete[] degrees;
    auto end = chrono::high_resolution_clock::now();

    cout << "Load CSR time:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Load CSR end." << endl;

    return graph;
}

void Graph::PreProcessing(Graph *graph) {
    graph->RemoveNonCoreVertex();
    graph->Relabel();
    graph->ComputeNeighborPlus();
    graph->ComputeMetaData();
}

void Graph::PrintGraphStatistics() {
    cout << PRINT_SEPARATOR << endl;
    cout << "|V| = " << vertices_num_ << " , " << "|E| = " << edges_num_ / 2 << " , max degree = " << max_degree_ << " , max degree sum = " << max_degree_sum_ << endl;
}

void Graph::StoreEdgeList(const Graph *graph, const string &file_path) {
    cout << PRINT_SEPARATOR << endl;
    cout << "Store edge list start..." << endl;
    auto start = chrono::high_resolution_clock::now();

    ofstream outfile(file_path);

    if (outfile.is_open()) {
        cout << "Open degree file " << file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot open degree file " << file_path << " ." << endl;
        exit(-1);
    }

    for (unsigned int i = 0; i < graph->vertices_num_; ++i) {
        unsigned int begin = i;
        unsigned int neighbors_count;
        const unsigned int* neighbors = graph->Neighbors(begin, neighbors_count);

        for (unsigned int j = 0; j < neighbors_count; ++j) {
            unsigned int end = neighbors[j];

            if (begin < end) {
                outfile << begin << ' ' << end << '\n';
            }
        }
    }

    outfile.close();
    outfile.clear();

    auto end = chrono::high_resolution_clock::now();
    cout << "Store edge list time:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Store edge list end." << endl;
}

void Graph::StoreCSR(const Graph *graph, const string &degree_file_path, const string &edge_file_path) {
    cout << PRINT_SEPARATOR << endl;
    cout << "Store CSR start..." << endl;
    auto start = chrono::high_resolution_clock::now();

    unsigned int* degrees = new unsigned int[graph->vertices_num_];
    for (unsigned int i = 0; i < graph->vertices_num_; ++i) {
        degrees[i] = graph->edges_offset_[i + 1] - graph->edges_offset_[i];
    }

    ofstream deg_outputfile(degree_file_path, ios::binary);

    if (deg_outputfile.is_open()) {
        cout << "Open degree file " << degree_file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot degree edge file " << degree_file_path << " ." << endl;
        exit(-1);
    }

    int int_size = sizeof(int);
    size_t vertex_array_bytes = ((size_t)graph->vertices_num_) * 4;
    deg_outputfile.write(reinterpret_cast<const char *>(&int_size), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&graph->vertices_num_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&graph->edges_num_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(degrees), vertex_array_bytes);

    deg_outputfile.close();
    deg_outputfile.clear();

    delete[] degrees;

    ofstream edge_outputfile(edge_file_path, ios::binary);

    if (edge_outputfile.is_open()) {
        cout << "Open edge file " << edge_file_path << " successfully." << endl;
    }
    else {
        cerr << "Cannot edge file " << edge_file_path << " ." << endl;
        exit(-1);
    }

    size_t edge_array_bytes = ((size_t)graph->edges_num_) * 4;
    edge_outputfile.write(reinterpret_cast<const char *>(graph->edges_), edge_array_bytes);

    edge_outputfile.close();
    edge_outputfile.clear();

    auto end = chrono::high_resolution_clock::now();
    cout << "Store CSR time:" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms." << endl;
    cout << "Store CSR end." << endl;
}

void Graph::ComputeMetaData() {
    // Generate max degree sum.
    max_degree_sum_ = 0;

    for (unsigned int i = 0; i < vertices_num_; ++i) {
        unsigned int neighbors_count;
        const unsigned int *neighbors = Neighbors(i, neighbors_count);

        unsigned int temp_sum = 0;

        for (unsigned int j = 0; j < neighbors_count; ++j) {
            unsigned int v = neighbors[j];

            temp_sum += Degree(v);
        }

        if (temp_sum > max_degree_sum_)
            max_degree_sum_ = temp_sum;
    }

    // Compute neighbor plus sum.

    neighbor_plus_sum_ = 0;
    for (unsigned int i = 0; i < vertices_num_; ++i) {
        neighbor_plus_sum_ += DegreePlus(i);
    }
}
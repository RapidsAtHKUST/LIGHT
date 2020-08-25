//
// Created by ssunah on 11/28/16.
//

#ifndef FSM_GRAPH_H
#define FSM_GRAPH_H

#include <cstring>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
using namespace std;

/*
 * The graph is stored as CSR format.
 * The id of vertices starts from 0, and the id of vertices is required to be consecutive.
 */

class Graph {
private:
    // Meta info of graph.
    unsigned int vertices_num_;                      // The number of vertices.
    unsigned int edges_num_;                         // The number of edges.
    unsigned int max_degree_;                        // The maximum degree of vertices.
    unsigned int max_degree_sum_;                    // The maximum degree sum of vertices.
    unsigned int neighbor_plus_sum_;                 // The sum of neighbor plus.
    // Graph data.
    unsigned int* edges_plus_offset_;                // The offset of neighbors plus.
    unsigned int* edges_offset_;                     // The offset of neighbors.
    unsigned int* edges_;                            // The neighbors of vertices.

private:
    void ConvertEdgeListToCSR(vector<pair<unsigned int, unsigned int>>& edge_list);
    void RemoveNonCoreVertex();
    void Relabel();
    void ComputeNeighborPlus();
    void ComputeMetaData();
public:
    // Static function.
    static void DestroyGraph(Graph* graph);
    static Graph* LoadEdgeList(const string& file_path);
    static Graph* LoadCSR(const string& degree_file_path, const string& edge_file_path);
    static void StoreEdgeList(const Graph* graph, const string& file_path);
    static void StoreCSR(const Graph* graph, const string& degree_file_path, const string& edge_file_path);
    static void PreProcessing(Graph* graph);
public:
    void PrintGraphStatistics();

public:
    // Inline function.
    Graph() {
        vertices_num_ = 0;
        edges_num_ = 0;
        max_degree_ = 0;
        edges_offset_ = NULL;
        edges_ = NULL;
        edges_plus_offset_ = NULL;
    }

    ~Graph() {
        delete[] edges_offset_;
        edges_offset_ = NULL;
        delete[] edges_;
        edges_ = NULL;
        delete[] edges_plus_offset_;
        edges_plus_offset_ = NULL;
    }

    inline const unsigned int Degree(const unsigned int id) const {
        return edges_offset_[id + 1] - edges_offset_[id];
    }

    inline const unsigned int DegreePlus(const unsigned int id) const {
        return edges_offset_[id + 1] - edges_offset_[id] - edges_plus_offset_[id];
    }

    inline const unsigned int MaxDegreeSum() const {
        return max_degree_sum_;
    }

    inline const unsigned int MaxNeighborPlusSum() const {
        return neighbor_plus_sum_;
    }

    inline const unsigned int* Neighbors(const unsigned int id, unsigned int& count) const {
        count = edges_offset_[id + 1] - edges_offset_[id];
        return edges_ + edges_offset_[id];
    }

    inline const unsigned int * NeighborsPlus(const unsigned int id, unsigned int & count) const {
        count = edges_offset_[id + 1] - edges_offset_[id] - edges_plus_offset_[id];
        return edges_ + edges_offset_[id] + edges_plus_offset_[id];
    }

    inline const unsigned int VerticesCount() const {
        return vertices_num_;
    }

    inline const unsigned int EdgesCount() const {
        return edges_num_;
    }

    inline const unsigned int MaxDegree() const {
        return max_degree_;
    }

    inline const bool IsEdge(const unsigned int u, const unsigned int v) const {
        unsigned int count = 0;
        const unsigned int* neighbors =  Neighbors(v, count);
        int begin = 0;
        int end = (int)count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
};


#endif //FSE_GRAPH_H

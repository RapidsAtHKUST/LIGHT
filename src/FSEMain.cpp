#include <iostream>
#include <iomanip>
#include "Config.h"
#include "Utility.h"
#include "FSECommand.h"
#include "FSE.h"

using namespace std;

int main(int argc, char** argv) {
    // Parse command.
    FSECommand command(argc, argv);

    const string data_graph_adj_file_path = command.GetDataGraphAdjFilePath();
    const string data_graph_degree_file_path = command.GetDataGraphDegreeFilePath();
    const string pattern_name = command.GetPatternName();
    const unsigned int thread_num = command.GetThreadNum();

    cout << "Query start..." << endl;

    Graph* data_graph = Graph::LoadCSR(data_graph_degree_file_path, data_graph_adj_file_path);
    data_graph->PrintGraphStatistics();

    cout << PRINT_SEPARATOR << endl;
    cout << "Enumeration start..." << endl;

    FSE fse;
    size_t result = fse.Enumerate(data_graph, pattern_name, thread_num);

    cout << "#Embeddings: " << result << " ." << endl;
    cout << "Search time: " << setprecision(4) << fse.search_time_ << " seconds." << endl;
    cout << "Total time: "  << setprecision(4) << fse.total_time_ << " seconds." << endl;
    cout << "Enumeration end." << endl;
    cout << PRINT_SEPARATOR << endl;
    cout << "Query end." << endl;

    Graph::DestroyGraph(data_graph);
    return 0;
}
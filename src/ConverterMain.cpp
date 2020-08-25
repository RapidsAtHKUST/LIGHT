//
// Created by ssunah on 11/30/17.
//

#include "ConverterCommand.h"
#include "Graph.h"
#include "Config.h"

int main(int argc, char** argv) {
    ConverterCommand command(argc, argv);

    string input_file_format = command.GetCommandValue(OptionKeyword::InputFileFormat);

    Graph* graph = NULL;

    if (input_file_format == "CSR") {
        string input_degree_file = command.GetCommandValue(OptionKeyword::InputDegreeFile);
        string input_adj_file = command.GetCommandValue(OptionKeyword::InputAdjFile);

        graph = Graph::LoadCSR(input_degree_file, input_adj_file);
    }
    else if (input_file_format == "EdgeList") {
        string input_edgelist_file = command.GetCommandValue(OptionKeyword::InputEdgeListFile);

        graph = Graph::LoadEdgeList(input_edgelist_file);
    }
    else {
        cout << "The input file format is not correct." << endl;
        exit(-1);
    }

    graph->PrintGraphStatistics();

    string operation = command.GetCommandValue(OptionKeyword::Command);

    if (operation == "Convert") {
        Graph::PreProcessing(graph);
    }
    else if (operation == "Analyze") {

    }
    else {
        cout << "The specified operation is not correct." << endl;
        exit(-1);
    }

    graph->PrintGraphStatistics();

    string output_file_format = command.GetCommandValue(OptionKeyword::OutputFileFormat);

    if (output_file_format == "CSR") {
        string output_degree_file = command.GetCommandValue(OptionKeyword::OutputDegreeFile);
        string output_adj_file = command.GetCommandValue(OptionKeyword::OutputAdjFile);

        Graph::StoreCSR(graph, output_degree_file, output_adj_file);
    }
    else if (output_file_format == "EdgeList") {
        string ouput_edgelist_file = command.GetCommandValue(OptionKeyword::OutputEdgeListFile);

        Graph::StoreEdgeList(graph, ouput_edgelist_file);
    }
    else {
        cout << PRINT_SEPARATOR << endl;
        cout << "No output file." << endl;
    }

    Graph::DestroyGraph(graph);

    return 0;
}
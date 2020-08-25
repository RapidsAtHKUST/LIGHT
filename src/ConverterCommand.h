//
// Created by ssunah on 11/30/17.
//

#ifndef FSE_CONVERTERCOMMAND_H_H
#define FSE_CONVERTERCOMMAND_H_H

#include "CommandParser.h"
#include <map>
using namespace std;

enum OptionKeyword {
    Command = 0,            // -o, convert or analysis
    InputFileFormat = 1,    // -ift, The input file format, CSR or EdgeList
    InputEdgeListFile = 2,  // -ief, The input edge list file path.
    InputDegreeFile = 3,    // -idf, The input degree file path.
    InputAdjFile = 4,       // -iaf, The input adjacent list file path.
    OutputFileFormat = 5,   // -oft, The output file format, CSR or EdgeList.
    OutputEdgeListFile = 6, // -oef, The output edge list file path.
    OutputDegreeFile = 7,   // -odf, The input degree file path.
    OutputAdjFile = 8       // -oaf, The input adjacent list file path.
};

class ConverterCommand : public CommandParser {
private:
    map<OptionKeyword, string> options_key;
    map<OptionKeyword, string> options_value;

private:
    void ProcessOptions() {
        for (auto option : options_key) {
            options_value[option.first] = GetCommandOption(option.second);
        }
    }

public:
    ConverterCommand(const int &argc, char **argv) : CommandParser(argc, argv) {
        // Initialize options value
        options_key[OptionKeyword::Command] = "-o";
        options_key[OptionKeyword::InputFileFormat] = "-ift";
        options_key[OptionKeyword::InputEdgeListFile] = "-ief";
        options_key[OptionKeyword::InputDegreeFile] = "-idf";
        options_key[OptionKeyword::InputAdjFile] = "-iaf";
        options_key[OptionKeyword::OutputFileFormat] = "-oft";
        options_key[OptionKeyword::OutputEdgeListFile] = "-oef";
        options_key[OptionKeyword::OutputDegreeFile] = "-odf";
        options_key[OptionKeyword::OutputAdjFile] = "-oaf";
        ProcessOptions();
    }

    const string& GetCommandValue(OptionKeyword option) {
        return options_value[option];
    }
};

#endif //FSE_CONVERTERCOMMAND_H_H

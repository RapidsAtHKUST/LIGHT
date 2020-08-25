//
// Created by ssunah on 11/28/16.
//

#ifndef FSM_FSMCOMMAND_H
#define FSM_FSMCOMMAND_H

#include "CommandParser.h"
#include <map>
using namespace std;

enum FSEOptionKeyword {
    InputDataGraphDegreeFilePath = 1,        // -df, The query graph file path, compulsive parameter.
    InputDataGraphAdjFilePath = 2,           // -af, The data graph file path, compulsive parameter.
    PatternName = 4,                         // -p,  The name of pattern.
    ThreadNum = 5                            // -n,  The number of threads.
};

class FSECommand : public CommandParser {
private:
    map<FSEOptionKeyword, string> options_key;
    map<FSEOptionKeyword, string> options_value;

private:
    void ProcessOptions() {

        // Degree file path
        string value = GetCommandOption(options_key[FSEOptionKeyword::InputDataGraphDegreeFilePath]);
        if (value == "") {
            printf("The degree file path cannot be empty.\n");
            exit(-1);
        }
        options_value[FSEOptionKeyword::InputDataGraphDegreeFilePath] = value;

        // Adj file path
        value = GetCommandOption(options_key[FSEOptionKeyword::InputDataGraphAdjFilePath]);
        if (value == "") {
            printf("The adjacent file path cannot be empty.\n");
            exit(-1);
        }
        options_value[FSEOptionKeyword::InputDataGraphAdjFilePath] = value;

        // Pattern name
        value = GetCommandOption(options_key[FSEOptionKeyword::PatternName]);
        if (value == "") {
            printf("The pattern name cannot be empty.\n");
            exit(-1);
        }
        options_value[FSEOptionKeyword::PatternName] = value;

        // Thread num
        value = GetCommandOption(options_key[FSEOptionKeyword::ThreadNum]);
        if (value == "") {
            printf("The thread number cannot be empty.\n");
            exit(-1);
        }
        options_value[FSEOptionKeyword::ThreadNum] = value;
    }

public:
    FSECommand(const int &argc, char **argv) : CommandParser(argc, argv) {
        // Initialize options value
        options_key[FSEOptionKeyword::InputDataGraphDegreeFilePath] = "-df";
        options_key[FSEOptionKeyword::InputDataGraphAdjFilePath] = "-af";
        options_key[FSEOptionKeyword::PatternName] = "-p";
        options_key[FSEOptionKeyword::ThreadNum] = "-n";
        ProcessOptions();
    }

    const string& GetDataGraphAdjFilePath() {
        return options_value[FSEOptionKeyword::InputDataGraphAdjFilePath];
    }

    const string& GetDataGraphDegreeFilePath() {
        return options_value[FSEOptionKeyword::InputDataGraphDegreeFilePath];
    }

    const string& GetPatternName() {
        return options_value[FSEOptionKeyword::PatternName];
    }

    const unsigned int GetThreadNum() {
        return (unsigned int)stoul(options_value[FSEOptionKeyword::ThreadNum]);
    }
};


#endif //FSM_FSMCOMMAND_H

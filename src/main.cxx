
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cvodesolver.h"
#if EXTERNAL_VARIABLES
extern "C"
{
#include "external_variables.h"
}
#endif
#include "json.hpp"
#include "model_header.h"


using DataStore = std::vector<std::vector<double>>;

static const auto StartingPoint = "StartingPoint";
static const auto EndingPoint = "EndingPoint";
static const auto PointInterval = "PointInterval";

void writeState(std::ofstream &file, double voi, size_t statesCount, double *states, size_t variablesCount, double *variables)
{
    file << voi << ",";
    for(size_t i = 0; i < statesCount; ++i) {
        file << states[i] << ",";
    }
    for(size_t i = 0; i < variablesCount; ++i) {
        file << variables[i] << ",";
    }
    file << std::endl;
}

void writeDataStore(std::ofstream &file, const DataStore& dataStore)
{
    for(const auto &entry: dataStore) {
        for(const auto& e: entry) {
            file << e << ",";
        }
        file << std::endl;
    }
}

void storeState(DataStore& dataStore, double voi, size_t statesCount, double *states, size_t variablesCount, double *variables)
{
    std::vector<double> entry;
    entry.push_back(voi);
    entry.insert(entry.end(), states, &states[statesCount]);
    entry.insert(entry.end(), variables, &variables[variablesCount]);

    dataStore.push_back(entry);
}

void initialiseArray(size_t size, double *array)
{
    for(size_t i = 0; i < size; ++i) {
        array[i]  = 0.0;
    }
}

Properties processConfig(const std::string &configFile)
{
    // Read a JSON config file.
    std::ifstream i(configFile);
    nlohmann::json j;
    i >> j;

    // Write out prettified JSON.
    Properties properties;
    // std::cout << std::setw(4) << j << std::endl;
    for (auto& el : j.items()) {
        if(j[el.key()].is_number_integer()) {
            properties[el.key()] = j[el.key()].get<int>();
        } else if (j[el.key()].is_boolean()) {
            properties[el.key()] = j[el.key()].get<bool>();
        } else if (j[el.key()].is_string()) {
            properties[el.key()] = j[el.key()].get<std::string>();
        } else if (j[el.key()].is_number_float()) {
            properties[el.key()] = j[el.key()].get<double>();
        }
    }

    return properties;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        return 1;
    }

    std::vector<std::string> params(argv, argv+argc);

    std::string solverConfigFile = params[1];
    std::string simulationConfigFile = params[2];
    std::string outputFile = "output.csv";
    if (argc > 3) {
      outputFile = params[3];
    }

    Properties solverProperties = processConfig(solverConfigFile);
    Properties simulationProperties = processConfig(simulationConfigFile);

    OdeSolver *solver = new CvodeSolver();
    solver->setProperties(solverProperties);

    double startingPoint = std::any_cast<double>(simulationProperties.at(StartingPoint));
    double endingPoint = std::any_cast<double>(simulationProperties.at(EndingPoint));
    double pointInterval = std::any_cast<double>(simulationProperties.at(PointInterval));

    size_t pointCounter = 0;

    double voi = startingPoint;
\
    double *states = createStatesArray();
    double *rates = createStatesArray();
    double *variables = createVariablesArray();

    initialiseArray(STATE_COUNT, states);
    initialiseArray(STATE_COUNT, rates);
    initialiseArray(VARIABLE_COUNT, variables);

    initialiseVariables(
#if EXTERNAL_VARIABLES
                  voi,
#endif
    states, variables
#if EXTERNAL_VARIABLES
                  , computeExternalVariable
#endif
);
    computeComputedConstants(variables);

//    std::ofstream dataMeta("output.meta");

//    if (dataMeta.is_open()) {
//        nlohmann::json j = {{"state_count", STATE_COUNT}, {"variable_count", VARIABLE_COUNT}};
//        dataMeta << j;
//        dataMeta.close();
//    }


    DataStore dataStore;
#if STORE_FILE
    std::ofstream dataFile(outputFile);
    if (dataFile.is_open()) {
#endif

    computeVariables(voi, states, rates, variables
#if EXTERNAL_VARIABLES
                  , computeExternalVariable
#endif
);
#if STORE_FILE
        writeState(dataFile, voi, STATE_COUNT, states, VARIABLE_COUNT, variables);
#else
    storeState(dataStore, voi, STATE_COUNT, states, VARIABLE_COUNT, variables);
#endif
    solver->initialize(voi, STATE_COUNT, states, rates, variables, computeRates);

    while (std::abs(voi - endingPoint) > 1e-08) {

        solver->solve(voi, std::min(endingPoint, startingPoint + double(++pointCounter) * pointInterval));

        computeVariables(voi, states, rates, variables
#if EXTERNAL_VARIABLES
                  , computeExternalVariable
#endif
);
#if STORE_FILE
        writeState(dataFile, voi, STATE_COUNT, states, VARIABLE_COUNT, variables);
#else
        storeState(dataStore, voi, STATE_COUNT, states, VARIABLE_COUNT, variables);
#endif
    }


//        writeDataStore(dataFile, dataStore);
#if STORE_FILE
        dataFile.close();
    }
#endif

    deleteArray(states);
    deleteArray(rates);
    deleteArray(variables);

    delete solver;

    return 0;
}

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "mlpcppwrapper.hpp"
#include "MLPCpp/include/CLookUp_ANN.hpp"

MLPCppEvaluator::MLPCppEvaluator(){};
MLPCppEvaluator::~MLPCppEvaluator() {
    delete lookup_mlp;
}
void MLPCppEvaluator::AddMLP(std::string mlp_file_name) {
    mlp_filenames_vector.push_back(mlp_file_name);
}

void MLPCppEvaluator::DisplayMLPs() {
    for (auto f : mlp_filenames_vector) std::cout << f << std::endl;
}

void MLPCppEvaluator::GenerateMLP() {
    size_t n_mlps = mlp_filenames_vector.size();
    std::string mlp_filenames[n_mlps];
    for (size_t i_mlp=0; i_mlp<n_mlps; i_mlp++) {
        mlp_filenames[i_mlp] = mlp_filenames_vector[i_mlp];
    }
    lookup_mlp = new MLPToolbox::CLookUp_ANN(n_mlps, mlp_filenames);
}

void MLPCppEvaluator::SetQueryInputs(std::vector<std::string>input_names){
    query_vars_in.clear();
    for (auto inp : input_names) query_vars_in.push_back(inp);
}
void MLPCppEvaluator::SetQueryOutputs(std::vector<std::string>output_names){
    query_vars_out.clear();
    for (auto outp : output_names) query_vars_out.push_back(outp);
    query_out.resize(query_vars_out.size());
    query_refs_out.resize(query_vars_out.size());
    for (size_t iq=0; iq<query_out.size(); iq++){
        query_refs_out[iq] = &query_out[iq];
    }
}

std::vector<std::vector<double>> MLPCppEvaluator::EvaluateMLP(std::vector<std::vector<double>> network_inputs) {
    std::vector<std::vector<double>> mlp_output;
    mlp_output.resize(network_inputs.size());
    for (size_t iOut=0; iOut<network_inputs.size(); iOut++) {
        mlp_output[iOut].resize(query_out.size());
    }
    iomap = new MLPToolbox::CIOMap(query_vars_in, query_vars_out);
    lookup_mlp->PairVariableswithMLPs(*iomap);
    for (size_t iq=0; iq<network_inputs.size(); iq++){
        lookup_mlp->Predict(*iomap, network_inputs[iq],query_refs_out);
        std::copy(query_out.begin(), query_out.end(), mlp_output[iq].begin());
    }
    delete iomap;
    return mlp_output;
}

double MLPCppEvaluator::GetMLPOuput(int i_out){
    return query_out[i_out];
}
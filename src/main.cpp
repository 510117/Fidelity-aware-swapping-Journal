#include <iostream>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <ctime>
#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/MyAlgo1/MyAlgo1.h"
#include "Algorithm/MyAlgo2/MyAlgo2.h"
#include "Algorithm/MyAlgo3/MyAlgo3.h"
#include "Algorithm/MyAlgo4/MyAlgo4.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
#include "Network/PathMethod/PathMethodBase/PathMethod.h"
#include "Network/PathMethod/Greedy/Greedy.h"
#include "Network/PathMethod/QCAST/QCAST.h"

using namespace std;

pair<int, int> generate_new_request(int num_of_node){
    //亂數引擎 
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node-1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);
    
    return make_pair(node1, node2);
}


int main(){
    string file_path = "../data/";

    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 50;
    default_setting["area_alpha"] = 0.0005;
    default_setting["time_limit"] = 14;
    default_setting["avg_memory"] = 7;
    default_setting["tao"] = 0.2;
    default_setting["path_length"] = -1;
    default_setting["min_fidelity"] = 0.7;
    default_setting["max_fidelity"] = 0.95;
    default_setting["swap_prob"] = 0.9;
    default_setting["fidelity_threshold"] = 0.5;

    map<string, vector<double>> change_parameter;
    change_parameter["request_cnt"] = {30, 40, 50, 60, 70};
    change_parameter["num_nodes"] = {40, 70, 100, 130, 160};
    change_parameter["min_fidelity"] = {0.5, 0.7, 0.75, 0.85, 0.95};
    change_parameter["avg_memory"] = {3, 5, 7, 9, 11};
    change_parameter["tao"] = {0.2, 0.4, 0.6, 0.8, 1};
    change_parameter["path_length"] = {11, 13, 15, 17};
    change_parameter["swap_prob"] = {0.5, 0.6, 0.7, 0.8, 0.9}; 
    change_parameter["fidelity_threshold"] = {0.5, 0.6, 0.7, 0.8, 0.9};

    // vector<string> X_names = {"time_limit", "request_cnt", "num_nodes", "avg_memory", "tao"};
    vector<string> X_names = {"request_cnt", "path_length", "avg_memory", "min_fidelity"};
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt", "utilization"};
    vector<string> algo_names = {"MyAlgo1", "MyAlgo2", "MyAlgo3", "Merge", "Linear"};
    // init result

    


    int round = 5;
    vector<PathMethod*> path_methods;
    path_methods.emplace_back(new Greedy());
    path_methods.emplace_back(new QCAST());
    for(PathMethod *path_method : path_methods) {

        for(string X_name : X_names) {
            for(string Y_name : Y_names){
                string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                fstream file( file_path + filename, ios::out );
            }
        }

        for(string X_name : X_names) {
            map<string, double> input_parameter = default_setting;

            for(double change_value : change_parameter[X_name]) {
                vector<map<string, map<string, double>>> result(round);
                input_parameter[X_name] = change_value;
                
                int num_nodes = input_parameter["num_nodes"];
                int avg_memory = input_parameter["avg_memory"];
                int memory_up = avg_memory + 1;
                int memory_lb = avg_memory - 1;
                int request_cnt = input_parameter["request_cnt"];
                int time_limit = input_parameter["time_limit"];
                double min_fidelity = input_parameter["min_fidelity"];
                double max_fidelity = input_parameter["max_fidelity"];
                double area_alpha = input_parameter["area_alpha"];
                double swap_prob = input_parameter["swap_prob"];
                double fidelity_threshold = input_parameter["fidelity_threshold"];
                int length_upper, length_lower;
                if(input_parameter["path_length"] == -1) {
                    length_upper = num_nodes;
                    length_lower = 5;
                } else {
                    length_upper = input_parameter["path_length"] + 1;
                    length_lower = input_parameter["path_length"] - 1;
                }

                int sum_has_path = 0;
                #pragma omp parallel for
                for(int r = 0; r < round; r++){
                    string round_str = to_string(r);
                    ofstream ofs;
                    ofs.open(file_path + "log/" + path_method->get_name() + "_" + X_name + "_in_" + to_string(change_value) + "_Round_" + round_str + ".log");

                    time_t now = time(0);
                    char* dt = ctime(&now);
                    cerr  << "時間 " << dt << endl << endl; 
                    ofs << "時間 " << dt << endl << endl; 

                    string filename = file_path + "input/round_" + round_str + ".input";
                    string command = "python3 graph_generator.py ";
                    string parameter = to_string(num_nodes) + " " + to_string(memory_lb) + " " + to_string(memory_up) + " " + " " + to_string(min_fidelity) + " " + to_string(max_fidelity) + " " + to_string(area_alpha);
                    if(system((command + filename + " " + parameter).c_str()) != 0){
                        cerr<<"error:\tsystem proccess python error"<<endl;
                        exit(1);
                    }

                    double A = 0.25, B = 0.75, tao = input_parameter["tao"], T = 10, n = 2;

                    Graph graph(filename, time_limit, swap_prob, fidelity_threshold, A, B, n, T, tao);

                    ofs << "--------------- in round " << r << " -------------" <<endl;
                    vector<pair<int, int>> requests;
                    for(int i = 0; i < request_cnt; i++) {
                        pair<int, int> new_request = generate_new_request(num_nodes);
                        int len = graph.distance(new_request.first, new_request.second);
                        int cnt = 1000;
                        while(len < length_lower || len > length_upper) {
                            new_request = generate_new_request(num_nodes);
                            len = graph.distance(new_request.first, new_request.second);
                            if(cnt == 0) break;
                            cnt--;
                        }
                        requests.push_back(new_request);
                    }

                    Graph path_graph = graph;
                    path_graph.increase_resources(10);
                    PathMethod *new_path_method;
                    if(path_method->get_name() == "Greedy") new_path_method = new Greedy();
                    else if(path_method->get_name() == "QCAST") new_path_method = new QCAST();
                    else {
                        cerr << "unknown path method" << endl;
                        assert(false);
                    }

                    new_path_method->build_paths(path_graph, requests);
                    map<SDpair, vector<Path>> paths = new_path_method->get_paths();

                    int path_len = 0, path_cnt = 0, mx_path_len = 0;

                    int has_path = 0;
                    for(auto P : paths) {
                        int mi_path_len = INF;
                        has_path += !P.second.empty();
                        for(Path path : P.second) {
                            mi_path_len = min(mi_path_len, (int)path.size());
                            for(int i = 1; i < (int)path.size(); i++) {
                                assert(graph.adj_set[path[i]].count(path[i - 1]));
                            }
                        }
                        if(mi_path_len != INF) {
                            mx_path_len = max(mx_path_len, mi_path_len);
                            path_cnt++;
                            path_len += mi_path_len;
                        }
                    }
                    
                    sum_has_path += has_path;
                    cerr << "Path method: " << path_method->get_name() << "\n";
                    cerr << "Request cnt: " << request_cnt << "\n";
                    cerr << "Has Path cnt: " << has_path << "\n";
                    cerr << "Avg path length = " << path_len / (double)path_cnt << "\n";
                    cerr << "Max path length = " << mx_path_len << "\n";
                    vector<AlgorithmBase*> algorithms;
                    algorithms.emplace_back(new MyAlgo1(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo2(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo3(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo4(graph, requests, paths));
                    algorithms.emplace_back(new MyAlgo5(graph, requests, paths));


                    #pragma omp parallel for
                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        algorithms[i]->run();
                    }


                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        for(string Y_name : Y_names) {
                            result[r][algorithms[i]->get_name()][Y_name] = algorithms[i]->get_res(Y_name);
                        }
                    }

                    now = time(0);
                    dt = ctime(&now);
                    cerr << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl; 
                    ofs.close();
                
                    for(auto &algo : algorithms){
                        delete algo;
                    }
                    algorithms.clear();
                
                }
                
                map<string, map<string, double>> sum_res;
                // for(string algo_name : algo_names){
                //     for(int r = 0; r < round; r++){
                //         result[r][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
                //         result[r][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
                //         result[r][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
                //         result[r][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
                //         result[r][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                //     }
                // }

                for(string Y_name : Y_names) {
                    string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                    ofstream ofs;
                    ofs.open(file_path + filename, ios::app);
                    ofs << change_value << ' ';
                    
                    for(string algo_name : algo_names){
                        for(int r = 0; r < round; r++){
                            sum_res[algo_name][Y_name] += result[r][algo_name][Y_name];
                        }
                        ofs << sum_res[algo_name][Y_name] / round << ' ';
                    }
                    ofs << sum_has_path / (double)round;
                    ofs << endl;
                    ofs.close();
                }
            }
        }
    }
    return 0;
}
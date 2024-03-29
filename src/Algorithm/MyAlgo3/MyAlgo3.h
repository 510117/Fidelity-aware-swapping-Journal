#ifndef __MYALGO3_H
#define __MYALGO3_H

#include    "../AlgorithmBase/AlgorithmBase.h"
#include    "../../Network/Graph/Graph.h"
#include    "../../config.h"

using namespace std;

class MyAlgo3 : public AlgorithmBase {
    vector<vector<vector<vector<double>>>> dp, dp2;
    vector<vector<vector<vector<pair<int, int>>>>> par, par2;
    vector<vector<vector<vector<bool>>>> caled, caled2;
    pair<Shape, double> calculate_best_shape(int src, int dst);
    pair<Shape, double> calculate_best_shape2(int src, int dst);
    double solve_fidelity(int left, int right, int t, int state, vector<int> &path);
    double solve2(int left, int right, int t, int state, vector<int> &path);
    double cp_value(Shape shape);
    Shape_vector backtracing_shape(int left, int right, int t, int state, vector<int> &path);
    Shape_vector backtracing_shape2(int left, int right, int t, int state, vector<int> &path);
public:
    MyAlgo3(Graph graph, vector<pair<int, int>> requests, map<SDpair, vector<Path>> paths);
    void run();
};

#endif
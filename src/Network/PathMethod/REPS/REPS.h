#ifndef __REPS_H
#define __REPS_H

#include "../PathMethodBase/PathMethod.h"

using namespace std;

class REPS : public PathMethod {
    void PFT_LP(vector<double> &t_plum, vector<map<pair<int, int>, double>> &f_plum);
    pair<vector<int>, double> REPS::dijkstra(int src, int dst, map<pair<int, int>, double>&f_plum_i);
    void update_f_plum(Path path, double width, map<pair<int, int>, double>&f_plum_i);
public:
    REPS();
    ~REPS();
    void build_paths(Graph graph, vector<SDpair> requests);
};

#endif
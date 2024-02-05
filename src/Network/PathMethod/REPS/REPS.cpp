#include "REPS.h"

using namespace std;

REPS::REPS() {
    method_name = "REPS";
}

REPS::~REPS() {}

void REPS::PFT_LP(vector<double> &t_plum, vector<map<pair<int, int>, double>> &f_plum) {
    
    //return value is store in t_plum and f_plum
    t_plum.clear();
    f_plum.clear();
    
    //do LP
try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
        vector<map<pair<int, int>, GRBVar>> f(1000);    //fi(u, v)
    GRBModel model = GRBModel(env);
    for(int i = 0; i < 100; i++){
        for(int u = 0; u < 100; u++){
            for(int v = 0; v < 100; v++){
                                f[i][make_pair(u, v)] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "f" + to_string(i) + "("+to_string(u)+", "+to_string(v) + ")");
                        }
        }
    }
        GRBLinExpr expr = 0;

    // Create variables
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    // Set objective: maximize x + y + 2 z
    model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

    // Add constraint: x + 2 y + 3 z <= 4
    model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

    // Add constraint: x + y >= 1
    model.addConstr(x + y >= 1, "c1");

    // Optimize model
    model.optimize();

    cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
    cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
    cout << z.get(GRB_StringAttr_VarName) << " "
         << z.get(GRB_DoubleAttr_X) << endl;

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
}

void REPS::build_paths(Graph _graph, vector<SDpair> _requests) {
    paths.clear();
    graph = _graph;
    requests = _requests;

    //PFT Using Progressive Rounding
    vector<double> t_plum;
    vector<map<pair<int, int>, double>> f_plum;
    PFT_LP(t_plum, f_plum);
    vector<int> neighbor;
    bool flag = true;
    double width;
    vector<int> path;

    while(flag){
        paths.clear();
        flag = false;

        vector<pair<double, Path>> paths;
        for(int i = 0; i < (int)requests.size(); i++){
            int src = requests[i].first, dst = requests[i].second;
            while(true) {
                tie(path, width) = dijkstra(src, dst, f_plum[i]);
                if(path.empty() || width >= EPS) break;
                int width_floor = floor(width);
                if(width_floor >= 1) {
                    reserve_path(path, width_floor);
                    flag = true;
                }
                update_f_plum(path, width, f_plum[i]);
                paths.emplace_back(width - width_floor, path);
            }
        }

        sort(paths.rbegin(), paths.rend());
        for(auto P : paths) {
            Path path = P.second;
            if(graph.check_path_resource(path, 1)) {
                flag = true;
                reserve_path(path, 1);
            }
        }
        // cout << "call PFT_LP in REPS::path_assignment()" << endl;
        PFT_LP(t_plum, f_plum);
        // cout << "call PFT_LP in REPS::path_assignment()--end" << endl;
    }
}

pair<Path, double> REPS::dijkstra(int src, int dst, map<pair<int, int>, double> &f_plum_i) {
    vector<bool> vis(graph.get_num_nodes(), false);
    vector<double> flow(graph.get_num_nodes(), 0);
    vector<int> par(graph.get_num_nodes(), -1);
    priority_queue<pair<double, int>> pq;
    pq.push({INF, src});
    flow[src] = INF;
    while(!pq.empty()) {
        int cur = pq.top().second;
        if(vis[cur]) continue;
        vis[cur] = true;
        for(int v : graph.adj_list[cur]) {
            if(f_plum_i[{cur, v}] <= EPS) continue;
            double new_flow = min(flow[cur], f_plum_i[{cur, v}]);
            if(flow[v] < new_flow) {
                flow[v] = new_flow;
                pq.push({flow[v], v});
                par[v] = cur;
            }
        }
    }

    if(!vis[dst]) return {{}, -1};

    Path path;
    int cur = dst;

    while(cur != -1) {
        path.push_back(cur);
        cur = par[cur];
    }

    reverse(path.begin(), path.end());
    return {path, flow[dst]};
}

void REPS::update_f_plum(Path path, double width, map<pair<int, int>, double>&f_plum_i) {
    for(int i = 1; i < (int)path.size(); i++) {
        int u = path[i - 1], v = path[i];
        if(f_plum_i[{u, v}] + EPS < width) {
            cerr << "error: f_plum_i is not enough" << endl;
            cerr << "f_plum = " << f_plum_i[{u, v}] << " width = " << width << endl;
            exit(1);
        }
        f_plum_i[{u, v}] -= width;
    }
}
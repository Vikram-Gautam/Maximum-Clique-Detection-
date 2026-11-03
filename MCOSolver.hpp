#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <chrono>
using namespace std;

struct Graph {
    int n = 0;       //vetices                     
    int m = 0;       //edges                       
    vector<vector<int>> adj_d;   //for dense
    vector<vector<int>> adj_l;     //for sparse      
    vector<int> degree;                     

    Graph() = default;

    void init(int vertices) {
        n = vertices;
        m = 0;
        adj_d.assign(n, vector<int>(n, 0));
        adj_l.assign(n, vector<int>());
        degree.assign(n, 0);
    }

    void add_edge(int u,int v) {
        if (u<0 || u>=n || v<0 || v>=n || u==v) return;
        if (adj_d[u][v]==0){
            adj_d[u][v] = 1;
            adj_d[v][u] = 1;
            adj_l[u].push_back(v);
            adj_l[v].push_back(u);
            degree[u]++;
            degree[v]++;
            m++;
        }
    }

    bool loadFromFile(const string& filename) {
        ifstream infile(filename);
        if (!infile.is_open()) {
            cerr<<"couldn't open"<<filename<<endl;
            return false;
        }

        string line;
        bool flag = false;

        while (getline(infile, line)) {
            if (line.empty()||line[0]=='c') continue; 

            if (line[0]=='p') {
                stringstream ss(line);
                string p,edge;
                int verticesNum,edgeNum;
                ss>>p>>edge>>verticesNum>>edgeNum;
                init(verticesNum);
                flag = true;
            } 
            else if (line[0]=='e') {
                stringstream ss(line);
                char e;
                int u,v;
                ss>>e>>u>> v;
                add_edge(u-1,v-1);
            } 
            else if(!flag){
                stringstream ss(line);
                int verticesNum, edgeNum;
                if (ss>>verticesNum>>edgeNum) {
                    init(verticesNum);
                    flag=true;
                    int u,v;
                    while(infile>>u>>v) {
                        add_edge(u-1,v-1);
                    }
                }
            }
        }
        infile.close();
        return true;
    }
};

class MCOSolver {
private:
    const Graph& G;
    int maxSize;
    vector<int> bestSol;
    long long nodes;
    double executionTime;
    
    bool timeOut;
    chrono::time_point<chrono::high_resolution_clock> startTime;
    double timeLimit;

    void saveSol(const vector<int>& C) {
        bestSol= C;
        maxSize =C.size();
    }

    void expand(vector<int>& C, vector<int>& P) {
        if (timeOut) return;

        if ((nodes&4095)==0) {
            auto current = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = current - startTime;
            if (elapsed.count()>=timeLimit) {
                timeOut=true;
                return;
            }
        }

        nodes++;

        for (int i=static_cast<int>(P.size()) -1;i>=0;i--) {
            if (timeOut) return; 
            if (C.size()+P.size()<=static_cast<size_t>(maxSize)) {
                return;
            }

            int v = P[i];
            C.push_back(v);

            vector<int> newP;
            newP.reserve(i); 

            for (int j=0; j<=i;j++) {
                int w = P[j];
                if (G.adj_d[v][w] == 1) {
                    newP.push_back(w);
                }
            }

            if (newP.empty() && C.size() > static_cast<size_t>(maxSize)) {
                saveSol(C);
            }

            if (!newP.empty()) {
                expand(C,newP);
            }

            C.pop_back();
            P.pop_back();
        }
    }

public:
    explicit MCOSolver(const Graph& g,double limitSec =3600.0) 
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false), timeLimit(limitSec){}

    vector<int> solve() {
        nodes =0;
        maxSize =0;
        bestSol.clear();
        timeOut= false;

        vector<int> C;
        vector<int> P;

        P.reserve(G.n);
        for (int i=0; i< G.n;i++) {
            P.push_back(i);
        }

        startTime =chrono::high_resolution_clock::now();

        expand(C,P);

        auto end =chrono::high_resolution_clock::now();
        chrono::duration<double,milli> elapsed =end -startTime;
        executionTime =elapsed.count();

        return bestSol;
    }

    bool hastimeOut() const{return timeOut;}
    long long getNodesExplored() const{return nodes;}
    double getexecutionTime() const{return executionTime;}
    int getMaxSize() const{return maxSize;}
};
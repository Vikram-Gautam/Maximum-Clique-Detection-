#pragma once
#include "MCOSolver.hpp"   // for Graph

class MCRSolver {
private:
    const Graph& G;
    int maxSize;
    vector<int> bestSol;
    long long nodes;
    double executionTime;

    bool timeOut;
    chrono::time_point<chrono::high_resolution_clock> startTime;
    double timeLimit;

    vector<int> origId;            // new index -> original vertex
    vector<vector<char>> adj;      // adjacency matrix in new indexing

    void saveSol(const vector<int>& C) {
        bestSol.clear();
        for (int v : C) bestSol.push_back(origId[v]);
        maxSize = C.size();
    }

    // minimum-width ordering: repeatedly remove the min-degree vertex and put it last.
    // ties: smaller sum of neighbour degrees
    void buildOrdering() {
        int n = G.n;
        vector<int> curDeg(G.degree);
        vector<char> removed(n, 0);
        origId.assign(n, 0);

        for (int pos = n - 1; pos >= 0; pos--) {
            int best = -1;
            long long bestSum = 0;
            for (int v = 0; v < n; v++) {
                if (removed[v]) continue;
                if (best == -1 || curDeg[v] < curDeg[best]) {
                    best = v;
                    bestSum = -1;
                } else if (curDeg[v] > curDeg[best]) {
                    continue;
                }
                // equal degree: compare neighbour degree sums
                if (bestSum == -1) {
                    bestSum = 0;
                    for (int u : G.adj_l[best]) if (!removed[u]) bestSum += curDeg[u];
                }
                if (v != best) {
                    long long s = 0;
                    for (int u : G.adj_l[v]) if (!removed[u]) s += curDeg[u];
                    if (s < bestSum) { best = v; bestSum = s; }
                }
            }
            origId[pos] = best;
            removed[best] = 1;
            for (int u : G.adj_l[best]) if (!removed[u]) curDeg[u]--;
        }

        adj.assign(n, vector<char>(n, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                adj[i][j] = G.adj_d[origId[i]][origId[j]];
    }

    // greedy sequential colouring; output vertices class by class with their colour number
    void colorSort(const vector<int>& P, vector<int>& order, vector<int>& color) {
        vector<vector<int>> classes;
        for (int v : P) {
            size_t k = 0;
            for (; k < classes.size(); k++) {
                bool ok = true;
                for (int u : classes[k]) {
                    if (adj[v][u]) { ok = false; break; }
                }
                if (ok) break;
            }
            if (k == classes.size()) classes.emplace_back();
            classes[k].push_back(v);
        }
        order.clear();
        color.clear();
        for (size_t k = 0; k < classes.size(); k++) {
            for (int v : classes[k]) {
                order.push_back(v);
                color.push_back(static_cast<int>(k) + 1);
            }
        }
    }

    void expand(vector<int>& C, const vector<int>& P) {
        if (timeOut) return;

        if ((nodes & 4095) == 0) {
            auto current = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = current - startTime;
            if (elapsed.count() >= timeLimit) {
                timeOut = true;
                return;
            }
        }

        nodes++;

        vector<int> order, color;
        colorSort(P, order, color);

        for (int i = static_cast<int>(order.size()) - 1; i >= 0; i--) {
            if (timeOut) return;
            if (static_cast<int>(C.size()) + color[i] <= maxSize) return;

            int v = order[i];
            C.push_back(v);

            vector<int> newP;
            newP.reserve(i);
            for (int j = 0; j < i; j++) {
                int w = order[j];
                if (adj[v][w]) newP.push_back(w);
            }

            if (newP.empty()) {
                if (static_cast<int>(C.size()) > maxSize) saveSol(C);
            } else {
                expand(C, newP);
            }

            C.pop_back();
        }
    }

public:
    explicit MCRSolver(const Graph& g, double limitSec = 3600.0)
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false), timeLimit(limitSec) {}

    vector<int> solve() {
        nodes = 0;
        maxSize = 0;
        bestSol.clear();
        timeOut = false;

        buildOrdering();

        vector<int> C;
        vector<int> P(G.n);
        for (int i = 0; i < G.n; i++) P[i] = i;

        startTime = chrono::high_resolution_clock::now();

        expand(C, P);

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> elapsed = end - startTime;
        executionTime = elapsed.count();

        return bestSol;
    }

    bool hastimeOut() const { return timeOut; }
    long long getNodesExplored() const { return nodes; }
    double getexecutionTime() const { return executionTime; }
    int getMaxSize() const { return maxSize; }
};
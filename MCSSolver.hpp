#pragma once
#include "MCOSolver.hpp"   // for Graph

class MCSSolver {
private:
    const Graph& G;
    int maxSize;
    vector<int> bestSol;
    long long nodes;
    double executionTime;

    bool timeOut;
    chrono::time_point<chrono::high_resolution_clock> startTime;
    double timeLimit;

    int reorderDepth;              // subproblems with |C| < reorderDepth get re-sorted (0 = off)

    vector<int> origId;
    vector<vector<char>> adj;

    void saveSol(const vector<int>& C) {
        bestSol.clear();
        for (int v : C) bestSol.push_back(origId[v]);
        maxSize = C.size();
    }

    // same root ordering as MCR
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

    // minimum-width sort of a subproblem on its induced subgraph
    void sortByMinWidth(vector<int>& P) {
        int k = P.size();
        vector<int> deg(k, 0);
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                if (a != b && adj[P[a]][P[b]]) deg[a]++;

        vector<char> used(k, 0);
        vector<int> res(k);
        for (int pos = k - 1; pos >= 0; pos--) {
            int best = -1;
            for (int a = 0; a < k; a++) {
                if (used[a]) continue;
                if (best == -1 || deg[a] < deg[best]) best = a;
            }
            res[pos] = P[best];
            used[best] = 1;
            for (int b = 0; b < k; b++)
                if (!used[b] && adj[P[best]][P[b]]) deg[b]--;
        }
        P = res;
    }

    int countNeighbors(int v, const vector<int>& cls, int& lastNeighbor) {
        int cnt = 0;
        for (int u : cls) {
            if (adj[v][u]) {
                cnt++;
                lastNeighbor = u;
                if (cnt > 1) return cnt;
            }
        }
        return cnt;
    }

    // Re-NUMBER: try to place p in a class with index < thr by moving one conflicting vertex.
    // classes are 0-indexed, class k has colour k+1; thr = maxSize - |C|
    bool reNumber(int p, vector<vector<int>>& classes, int thr) {
        for (int k1 = 0; k1 <= thr - 2; k1++) {
            int q = -1;
            int cnt = countNeighbors(p, classes[k1], q);
            if (cnt != 1) continue;
            for (int k2 = k1 + 1; k2 <= thr - 1; k2++) {
                int dummy = -1;
                if (countNeighbors(q, classes[k2], dummy) == 0) {
                    // move q from k1 to k2, put p into k1
                    vector<int>& c1 = classes[k1];
                    c1.erase(find(c1.begin(), c1.end(), q));
                    classes[k2].push_back(q);
                    c1.push_back(p);
                    return true;
                }
            }
        }
        return false;
    }

    void colorSort(const vector<int>& P, int thr, vector<int>& order, vector<int>& color) {
        vector<vector<int>> classes;
        for (int v : P) {
            size_t k = 0;
            for (; k < classes.size(); k++) {
                int dummy = -1;
                if (countNeighbors(v, classes[k], dummy) == 0) break;
            }
            // v would become a branching vertex -> try to recolour
            if (static_cast<int>(k) >= thr && thr >= 2) {
                if (reNumber(v, classes, thr)) continue;
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

    void expand(vector<int>& C, vector<int>& P) {
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

        if (static_cast<int>(C.size()) < reorderDepth && P.size() > 2) {
            sortByMinWidth(P);
        }

        int thr = maxSize - static_cast<int>(C.size());
        if (thr < 0) thr = 0;

        vector<int> order, color;
        colorSort(P, thr, order, color);

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
    explicit MCSSolver(const Graph& g, double limitSec = 3600.0, int reorderDepth_ = 2)
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false),
          timeLimit(limitSec), reorderDepth(reorderDepth_) {}

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
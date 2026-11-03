#pragma once
#include <vector>
#include <chrono>
#include <algorithm>
#include <cstdint>
#include "MCOSolver.hpp"   // for Graph

class MCFSolver {
private:
    const Graph& G;
    int maxSize;
    std::vector<int> bestSol;
    long long nodes;
    double executionTime;

    bool timeOut;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    double timeLimit;

    bool enableES;                 // Enable Essential Set preprocessing
    std::vector<int> origId;       // Maps current active vertices -> original vertex id
    std::vector<std::vector<int>> adj;

    void saveSol(const std::vector<int>& C, const std::vector<int>& mapping) {
        bestSol.clear();
        for (int v : C) {
            bestSol.push_back(mapping[v]);
        }
        maxSize = static_cast<int>(bestSol.size());
    }

    // Checks timeLimit periodically
    inline void checkTimeout() {
        if ((nodes & 4095) == 0) {
            auto current = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = current - startTime;
            if (elapsed.count() >= timeLimit) {
                timeOut = true;
            }
        }
    }

    // Fast greedy coloring upper bound for candidate set P[cite: 4]
    int greedyColoringBound(const std::vector<int>& P, const std::vector<std::vector<int>>& A) {
        int pSize = static_cast<int>(P.size());
        if (pSize == 0) return 0;

        std::vector<int> colors(pSize, -1);
        int maxColor = 0;

        for (int i = 0; i < pSize; i++) {
            int u = P[i];
            std::vector<char> used(pSize + 1, 0);
            for (int j = 0; j < i; j++) {
                if (A[u][P[j]] && colors[j] != -1) {
                    used[colors[j]] = 1;
                }
            }
            int cr = 0;
            while (cr < pSize && used[cr]) cr++;
            colors[i] = cr;
            if (cr + 1 > maxColor) maxColor = cr + 1;
        }
        return maxColor;
    }

    // Greedy lower-bound heuristic for MCF (Eblen et al.)[cite: 4]
    std::vector<int> findGreedyClique(int n, const std::vector<std::vector<int>>& A) {
        std::vector<int> deg(n, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (A[i][j]) deg[i]++;
            }
        }

        // Pick highest degree vertex[cite: 4]
        int bestV = 0;
        for (int i = 1; i < n; i++) {
            if (deg[i] > deg[bestV]) bestV = i;
        }

        std::vector<int> C = { bestV };
        std::vector<int> cand;
        for (int i = 0; i < n; i++) {
            if (A[bestV][i]) cand.push_back(i);
        }

        while (!cand.empty()) {
            int nextV = cand[0];
            for (int u : cand) {
                if (deg[u] > deg[nextV]) nextV = u;
            }
            C.push_back(nextV);

            std::vector<int> nextCand;
            for (int u : cand) {
                if (u != nextV && A[nextV][u]) {
                    nextCand.push_back(u);
                }
            }
            cand = std::move(nextCand);
        }
        return C;
    }

    // MCF core branch-and-bound recursion[cite: 4]
    void mcfBranch(std::vector<int>& C, std::vector<int>& P, 
                   const std::vector<std::vector<int>>& A, 
                   const std::vector<int>& mapping, int& currentMax, std::vector<int>& currentBest) {
        if (timeOut) return;
        checkTimeout();
        nodes++;

        if (static_cast<int>(C.size() + P.size()) <= currentMax) return;

        // Bounding by vertex coloring[cite: 4]
        int colBound = greedyColoringBound(P, A);
        if (static_cast<int>(C.size()) + colBound <= currentMax) return;

        for (int i = static_cast<int>(P.size()) - 1; i >= 0; i--) {
            if (timeOut) return;
            if (static_cast<int>(C.size()) + i + 1 <= currentMax) return;

            int v = P[i];
            C.push_back(v);

            std::vector<int> newP;
            newP.reserve(i);
            for (int j = 0; j < i; j++) {
                int u = P[j];
                if (A[v][u]) newP.push_back(u);
            }

            if (newP.empty()) {
                if (static_cast<int>(C.size()) > currentMax) {
                    currentMax = static_cast<int>(C.size());
                    currentBest.clear();
                    for (int x : C) currentBest.push_back(mapping[x]);
                }
            } else {
                mcfBranch(C, newP, A, mapping, currentMax, currentBest);
            }

            C.pop_back();
        }
    }

    // Subroutine to run standalone MCF on an arbitrary adjacency matrix[cite: 4]
    std::vector<int> runMCF(int n, const std::vector<std::vector<int>>& A, const std::vector<int>& mapping) {
        if (n == 0) return {};

        std::vector<int> greedyC = findGreedyClique(n, A);
        int lowerBound = static_cast<int>(greedyC.size());
        std::vector<int> currentBest;
        for (int v : greedyC) currentBest.push_back(mapping[v]);
        int currentMax = lowerBound;

        // Degree and color preprocessing[cite: 4]
        std::vector<char> active(n, 1);
        bool changed = true;
        while (changed) {
            changed = false;
            for (int i = 0; i < n; i++) {
                if (!active[i]) continue;
                int d = 0;
                for (int j = 0; j < n; j++) {
                    if (active[j] && A[i][j]) d++;
                }
                if (d < currentMax - 1) {
                    active[i] = 0;
                    changed = true;
                }
            }
        }

        // Universal vertices extraction (deg == active_vertices - 1)[cite: 4]
        std::vector<int> universal;
        int activeCount = 0;
        for (int i = 0; i < n; i++) if (active[i]) activeCount++;

        for (int i = 0; i < n; i++) {
            if (!active[i]) continue;
            int d = 0;
            for (int j = 0; j < n; j++) {
                if (active[j] && A[i][j]) d++;
            }
            if (d == activeCount - 1) {
                universal.push_back(i);
                active[i] = 0;
            }
        }

        std::vector<int> subMap;
        for (int i = 0; i < n; i++) {
            if (active[i]) subMap.push_back(i);
        }

        int subN = static_cast<int>(subMap.size());
        if (subN == 0) {
            if (static_cast<int>(universal.size()) > currentMax) {
                currentBest.clear();
                for (int v : universal) currentBest.push_back(mapping[v]);
            }
            return currentBest;
        }

        std::vector<std::vector<int>> subA(subN, std::vector<int>(subN, 0));
        for (int i = 0; i < subN; i++) {
            for (int j = 0; j < subN; j++) {
                subA[i][j] = A[subMap[i]][subMap[j]];
            }
        }

        std::vector<int> fullSubMapping(subN);
        for (int i = 0; i < subN; i++) fullSubMapping[i] = mapping[subMap[i]];

        std::vector<int> C;
        std::vector<int> P(subN);
        for (int i = 0; i < subN; i++) P[i] = i;

        int adjustedMax = currentMax - static_cast<int>(universal.size());
        std::vector<int> branchBest;
        mcfBranch(C, P, subA, fullSubMapping, adjustedMax, branchBest);

        if (static_cast<int>(branchBest.size() + universal.size()) > currentMax) {
            currentBest = branchBest;
            for (int u : universal) currentBest.push_back(mapping[u]);
        }
        return currentBest;
    }

    // Essential Set (ES) Reduction (Eblen et al., Fig. 3)[cite: 6]
    void applyEssentialSetReduction() {
        int n = G.n;
        origId.resize(n);
        for (int i = 0; i < n; i++) origId[i] = i;

        adj = G.adj_d;

        // Step 1: Initial maximum clique M <- MCF(G)[cite: 6]
        std::vector<int> M = runMCF(n, adj, origId);
        if (!M.empty() && static_cast<int>(M.size()) > maxSize) {
       bestSol = M;
       maxSize = static_cast<int>(bestSol.size());
        }


        if (M.empty() || timeOut) return;

        bestSol = M;
        maxSize = static_cast<int>(bestSol.size());

        // Step 2: For each vertex v in M, test if v is essential: omega(G \ {v}) < |M|[cite: 5, 6]
        std::vector<char> isEssential(n, 0);
        for (int vOrig : M) {
            if (timeOut) break;

            std::vector<int> gPrimeMap;
            for (int i = 0; i < n; i++) {
                if (i != vOrig) gPrimeMap.push_back(i);
            }

            int nPrime = static_cast<int>(gPrimeMap.size());
            std::vector<std::vector<int>> gPrimeAdj(nPrime, std::vector<int>(nPrime, 0));
            for (int i = 0; i < nPrime; i++) {
                for (int j = 0; j < nPrime; j++) {
                    gPrimeAdj[i][j] = adj[gPrimeMap[i]][gPrimeMap[j]];
                }
            }

            std::vector<int> mPrime = runMCF(nPrime, gPrimeAdj, gPrimeMap);
            if (static_cast<int>(mPrime.size()) < maxSize) {
                // v is essential: it belongs to every maximum clique[cite: 5, 6]
                isEssential[vOrig] = 1;
            }
        }

        // Step 3: Remove all vertices not adjacent to essential vertices[cite: 5, 6]
        std::vector<char> keep(n, 1);
        bool filtered = false;
        for (int i = 0; i < n; i++) {
            if (!isEssential[i]) continue;
            for (int j = 0; j < n; j++) {
                if (i != j && !adj[i][j]) {
                    keep[j] = 0;
                    filtered = true;
                }
            }
        }

        if (filtered) {
            std::vector<int> newOrigId;
            for (int i = 0; i < n; i++) {
                if (keep[i]) newOrigId.push_back(origId[i]);
            }
            int newN = static_cast<int>(newOrigId.size());
            std::vector<std::vector<int>> newAdj(newN, std::vector<int>(newN, 0));
            for (int i = 0; i < newN; i++) {
                for (int j = 0; j < newN; j++) {
                    newAdj[i][j] = G.adj_d[newOrigId[i]][newOrigId[j]];
                }
            }
            origId = std::move(newOrigId);
            adj = std::move(newAdj);
        }
    }

public:
    explicit MCFSolver(const Graph& g, double limitSec = 3600.0, bool useES = true)
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false),
          timeLimit(limitSec), enableES(useES) {}

    std::vector<int> solve() {
        nodes = 0;
        maxSize = 0;
        bestSol.clear();
        timeOut = false;

        startTime = std::chrono::high_resolution_clock::now();

        if (enableES) {
            applyEssentialSetReduction();
        } else {
            origId.resize(G.n);
            for (int i = 0; i < G.n; i++) origId[i] = i;
            adj = G.adj_d;
        }

      if (!timeOut) {
    std::vector<int> finalSol = runMCF(static_cast<int>(origId.size()), adj, origId);
    if (static_cast<int>(finalSol.size()) > maxSize) {
        bestSol = std::move(finalSol);
        maxSize = static_cast<int>(bestSol.size());
    }
}

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - startTime;
        executionTime = elapsed.count();

        return bestSol;
    }

    bool hastimeOut() const { return timeOut; }
    long long getNodesExplored() const { return nodes; }
    double getexecutionTime() const { return executionTime; }
    int getMaxSize() const { return maxSize; }
};
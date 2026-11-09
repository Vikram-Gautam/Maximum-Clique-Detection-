#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>
#include <iostream>
#include "MCOSolver.hpp"   // for Graph

class RMCSolver {
public:
    struct Seed {
        std::vector<int> s;           // Current growing clique s
        std::vector<int> cand;        // Candidate set Lambda(s)
        int ub;                       // Upper bound of omega(s U Lambda(s))

        Seed() : ub(0) {}
        Seed(const std::vector<int>& s_, const std::vector<int>& cand_, int ub_ = 0)
            : s(s_), cand(cand_), ub(ub_) {}
    };

private:
    const Graph& G;
    int maxSize;
    std::vector<int> bestSol;
    long long nodes;
    double executionTime;

    bool timeOut;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    double timeLimit;

    double cParam;
    int theta;
    std::mt19937_64 rng;
    long long opsCounter;

    inline void checkTimeout() {
        if ((++opsCounter & 511) == 0) {
            auto current = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = current - startTime;
            if (elapsed.count() >= timeLimit) {
                timeOut = true;
            }
        }
    }

    void updateChampion(const std::vector<int>& clique) {
        if (static_cast<int>(clique.size()) > maxSize) {
            maxSize = static_cast<int>(clique.size());
            bestSol = clique;
        }
    }

    // O(m) Core Decomposition on a candidate subset
    std::vector<int> computeSubCoreNumbers(int m, const std::vector<int>& cand) {
        if (m == 0) return {};
        std::vector<int> deg(m, 0);
        int maxDeg = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i != j && G.adj_d[cand[i]][cand[j]]) deg[i]++;
            }
            if (deg[i] > maxDeg) maxDeg = deg[i];
        }

        std::vector<int> bin(maxDeg + 1, 0);
        for (int i = 0; i < m; i++) bin[deg[i]]++;
        int start = 0;
        for (int d = 0; d <= maxDeg; d++) {
            int num = bin[d];
            bin[d] = start;
            start += num;
        }

        std::vector<int> pos(m, 0), vert(m, 0);
        for (int i = 0; i < m; i++) {
            pos[i] = bin[deg[i]];
            vert[pos[i]] = i;
            bin[deg[i]]++;
        }
        for (int d = maxDeg; d >= 1; d--) bin[d] = bin[d - 1];
        bin[0] = 0;

        std::vector<int> core(m, 0);
        for (int i = 0; i < m; i++) {
            int v = vert[i];
            core[v] = deg[v];
            for (int u = 0; u < m; u++) {
                if (G.adj_d[cand[v]][cand[u]] && deg[u] > deg[v]) {
                    int du = deg[u];
                    int pu = pos[u];
                    int pw = bin[du];
                    int w = vert[pw];
                    if (u != w) {
                        pos[u] = pw; vert[pu] = w;
                        pos[w] = pu; vert[pw] = u;
                    }
                    bin[du]++;
                    deg[u]--;
                }
            }
        }
        return core;
    }

    // Greedy Degeneracy Coloring on candidate subset
    int computeSubColoring(int m, const std::vector<int>& cand, std::vector<int>& colors) {
        if (m == 0) { colors.clear(); return 0; }
        std::vector<int> deg(m, 0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i != j && G.adj_d[cand[i]][cand[j]]) deg[i]++;
            }
        }

        std::vector<char> used(m, 0);
        std::vector<int> order(m);
        for (int pos = m - 1; pos >= 0; pos--) {
            int best = -1;
            for (int i = 0; i < m; i++) {
                if (!used[i] && (best == -1 || deg[i] < deg[best])) best = i;
            }
            order[pos] = best;
            used[best] = 1;
            for (int j = 0; j < m; j++) {
                if (!used[j] && G.adj_d[cand[best]][cand[j]]) deg[j]--;
            }
        }

        colors.assign(m, -1);
        int maxCol = 0;
        std::vector<char> colUsed(m + 2, 0);
        for (int u : order) {
            std::fill(colUsed.begin(), colUsed.end(), 0);
            for (int j = 0; j < m; j++) {
                if (G.adj_d[cand[u]][cand[j]] && colors[j] != -1) {
                    if (colors[j] <= m + 1) colUsed[colors[j]] = 1;
                }
            }
            int c = 1;
            while (colUsed[c]) c++;
            colors[u] = c;
            if (c > maxCol) maxCol = c;
        }
        return maxCol;
    }

    // Candidate degree truss filter
    std::vector<int> fastTrussFilter(const std::vector<int>& cand, int k_truss) {
        int m = static_cast<int>(cand.size());
        if (m == 0 || k_truss <= 1) return cand;

        int reqDeg = k_truss - 1;
        std::vector<int> current = cand;
        bool changed = true;

        while (changed) {
            changed = false;
            int sz = static_cast<int>(current.size());
            if (sz < k_truss) return {};

            std::vector<int> nextCand;
            nextCand.reserve(sz);

            for (int i = 0; i < sz; i++) {
                int u = current[i];
                int d = 0;
                for (int j = 0; j < sz; j++) {
                    if (i != j && G.adj_d[u][current[j]]) d++;
                }
                if (d >= reqDeg) {
                    nextCand.push_back(u);
                } else {
                    changed = true;
                }
            }
            current = std::move(nextCand);
        }
        return current;
    }

    // Algorithm 3: Init(G)
    bool initBounds(int& lb, int& ub) {
        int n = G.n;
        if (n == 0) { lb = 0; ub = 0; return true; }

        std::vector<int> core(n, 0);
        std::vector<int> deg(n, 0);
        int maxDeg = 0;
        for (int i = 0; i < n; i++) {
            deg[i] = G.degree[i];
            if (deg[i] > maxDeg) maxDeg = deg[i];
        }

        std::vector<int> bin(maxDeg + 1, 0);
        for (int i = 0; i < n; i++) bin[deg[i]]++;
        int start = 0;
        for (int d = 0; d <= maxDeg; d++) {
            int num = bin[d];
            bin[d] = start;
            start += num;
        }

        std::vector<int> pos(n, 0), vert(n, 0);
        for (int i = 0; i < n; i++) {
            pos[i] = bin[deg[i]];
            vert[pos[i]] = i;
            bin[deg[i]]++;
        }
        for (int d = maxDeg; d >= 1; d--) bin[d] = bin[d - 1];
        bin[0] = 0;

        for (int i = 0; i < n; i++) {
            int v = vert[i];
            core[v] = deg[v];
            for (int u : G.adj_l[v]) {
                if (deg[u] > deg[v]) {
                    int du = deg[u];
                    int pu = pos[u];
                    int pw = bin[du];
                    int w = vert[pw];
                    if (u != w) {
                        pos[u] = pw; vert[pu] = w;
                        pos[w] = pu; vert[pw] = u;
                    }
                    bin[du]++;
                    deg[u]--;
                }
            }
        }

        int cm = 0;
        for (int v : core) if (v > cm) cm = v;

        ub = cm + 1;
        lb = 0;

        // Check if max-core is an exact complete clique
        std::vector<int> maxCoreVertices;
        for (int i = 0; i < n; i++) {
            if (core[i] == cm) maxCoreVertices.push_back(i);
        }
        if (static_cast<int>(maxCoreVertices.size()) == cm + 1) {
            bool isClique = true;
            for (size_t i = 0; i < maxCoreVertices.size() && isClique; i++) {
                for (size_t j = i + 1; j < maxCoreVertices.size(); j++) {
                    if (!G.adj_d[maxCoreVertices[i]][maxCoreVertices[j]]) {
                        isClique = false;
                        break;
                    }
                }
            }
            if (isClique) {
                lb = ub;
                updateChampion(maxCoreVertices);
                return true;
            }
        }

        // Greedy maximal clique initialization with candidate degree heuristic
        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return core[a] > core[b];
        });

        for (int v : order) {
            if (core[v] < lb) break;
            nodes++;
            std::vector<int> C = { v };
            std::vector<int> cand;
            for (int u : G.adj_l[v]) {
                if (core[u] >= lb) cand.push_back(u);
            }
            while (!cand.empty()) {
                int nextV = cand[0];
                int maxCandDeg = -1;
                for (int u : cand) {
                    int d = 0;
                    for (int w : cand) if (u != w && G.adj_d[u][w]) d++;
                    if (d > maxCandDeg) {
                        maxCandDeg = d;
                        nextV = u;
                    }
                }
                C.push_back(nextV);
                std::vector<int> nextCand;
                for (int u : cand) {
                    if (u != nextV && G.adj_d[nextV][u]) {
                        nextCand.push_back(u);
                    }
                }
                cand = std::move(nextCand);
            }
            if (static_cast<int>(C.size()) > lb) {
                lb = static_cast<int>(C.size());
                updateChampion(C);
            }
        }

        std::vector<int> colors;
        int cn = computeSubColoring(n, order, colors);
        if (ub > cn) ub = cn;

        if (lb >= ub) return true;
        return false;
    }

    // Algorithm 6: compIS(g = (V, E))
    void compIS(int m, const std::vector<int>& cand, std::vector<int>& P, std::vector<int>& I) {
        P.clear();
        I.clear();
        if (m == 0) return;

        std::vector<int> deg(m, 0);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i != j && !G.adj_d[cand[i]][cand[j]]) deg[i]++;
            }
        }

        std::vector<char> state(m, 0);
        bool t = true;

        while (true) {
            int bestV = -1;
            int minDeg = m + 1;
            for (int i = 0; i < m; i++) {
                if (state[i] == 0 && deg[i] < minDeg) {
                    minDeg = deg[i];
                    bestV = i;
                }
            }
            if (bestV == -1) break;

            state[bestV] = 1;
            I.push_back(bestV);

            if (minDeg > 1) t = false;
            if (t) {
                P.push_back(bestV);
            }

            for (int u = 0; u < m; u++) {
                if (u != bestV && !G.adj_d[cand[bestV]][cand[u]] && state[u] == 0) {
                    state[u] = 1;
                    for (int w = 0; w < m; w++) {
                        if (w != u && !G.adj_d[cand[u]][cand[w]] && state[w] == 0) {
                            deg[w]--;
                        }
                    }
                }
            }
        }
    }

    // Algorithm 4: scSeed(G, k)
    std::vector<Seed> scSeed(int k, int& lb, int& ub) {
        std::vector<Seed> S;
        int n = G.n;
        if (k <= 2) return S;

        // Sampling probability p per Theorem 7.1
        double denom = (double)k * (k - 1) * (k - 2);
        double p = 1.0;
        if (denom > 0.0) {
            p = std::sqrt((2.0 * cParam * std::log((double)n)) / denom);
        }
        if (p > 1.0) p = 1.0;

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // Fast Edge-list Triangle Sampling
        std::vector<std::vector<int>> sampled_adj(n);
        for (int u = 0; u < n; u++) {
            if (G.degree[u] + 1 < k) continue;
            for (int v : G.adj_l[u]) {
                if (u < v && G.degree[v] + 1 >= k) {
                    if (p >= 1.0 || dist(rng) <= p) {
                        sampled_adj[u].push_back(v);
                        sampled_adj[v].push_back(u);
                    }
                }
            }
        }

        for (int u = 0; u < n; u++) {
            checkTimeout();
            if (timeOut) return S;
            if (sampled_adj[u].size() < 2) continue;

            for (size_t i = 0; i < sampled_adj[u].size(); i++) {
                int v = sampled_adj[u][i];
                for (size_t j = i + 1; j < sampled_adj[u].size(); j++) {
                    int w = sampled_adj[u][j];
                    if (G.adj_d[v][w]) { // (v,w) is an edge in G
                        int sUb = std::min({G.degree[u], G.degree[v], G.degree[w]}) + 1;
                        if (sUb >= k) {
                            std::vector<int> cand;
                            for (int x : G.adj_l[u]) {
                                if (x != v && x != w && G.adj_d[v][x] && G.adj_d[w][x] && G.degree[x] + 1 >= k) {
                                    cand.push_back(x);
                                }
                            }
                            if (static_cast<int>(cand.size()) + 3 >= k) {
                                S.emplace_back(std::vector<int>{u, v, w}, cand, sUb);
                            }
                        }
                    }
                }
            }
        }

        std::sort(S.begin(), S.end(), [](const Seed& a, const Seed& b) {
            return a.ub > b.ub;
        });

        std::vector<Seed> filteredS;
        for (auto& seed : S) {
            checkTimeout();
            if (timeOut) return filteredS;

            if (seed.ub < k || static_cast<int>(seed.cand.size() + seed.s.size()) < k) continue;

            int m = static_cast<int>(seed.cand.size());
            std::vector<int> subCore = computeSubCoreNumbers(m, seed.cand);
            std::vector<int> filteredCand;
            int reqCore = k - static_cast<int>(seed.s.size()) - 1;

            for (int i = 0; i < m; i++) {
                if (subCore[i] >= reqCore) {
                    filteredCand.push_back(seed.cand[i]);
                }
            }
            seed.cand = std::move(filteredCand);

            int newM = static_cast<int>(seed.cand.size());
            std::vector<int> F;
            std::vector<int> remCand;
            for (int i = 0; i < newM; i++) {
                int u = seed.cand[i];
                int d = 0;
                for (int j = 0; j < newM; j++) {
                    if (i != j && G.adj_d[u][seed.cand[j]]) d++;
                }
                if (d == newM - 1) {
                    F.push_back(u);
                } else {
                    remCand.push_back(u);
                }
            }

            seed.s.insert(seed.s.end(), F.begin(), F.end());
            seed.cand = std::move(remCand);

            if (static_cast<int>(seed.s.size()) > lb) {
                lb = static_cast<int>(seed.s.size());
                updateChampion(seed.s);
            }

            if (lb >= k) k = lb + 1;

            if (!seed.cand.empty() && static_cast<int>(seed.s.size() + seed.cand.size()) >= k) {
                seed.ub = std::min(seed.ub, static_cast<int>(seed.s.size() + seed.cand.size()));
                filteredS.push_back(seed);
            }
        }

        if (!filteredS.empty()) {
            int maxSeedUb = 0;
            for (const auto& seed : filteredS) {
                if (seed.ub > maxSeedUb) maxSeedUb = seed.ub;
            }
            ub = std::min(ub, maxSeedUb);
        }

        return filteredS;
    }

    // Algorithm 5: tciSeed(G, S, k)
    std::vector<Seed> tciSeed(std::vector<Seed>& S, int k, int& lb, int& ub) {
        std::vector<Seed> resS;
        for (auto& seed : S) {
            nodes++;
            checkTimeout();
            if (timeOut) break;

            if (seed.ub < k || static_cast<int>(seed.s.size() + seed.cand.size()) < k) continue;

            int k_truss = k - static_cast<int>(seed.s.size());
            seed.cand = fastTrussFilter(seed.cand, k_truss);
            if (seed.cand.empty()) continue;

            int m = static_cast<int>(seed.cand.size());
            std::vector<int> colors;
            int numColors = computeSubColoring(m, seed.cand, colors);
            if (static_cast<int>(seed.s.size()) + numColors < k) continue;

            int reqColors = k - static_cast<int>(seed.s.size()) - 1;
            std::vector<int> filteredCand;
            for (int i = 0; i < m; i++) {
                int u = seed.cand[i];
                int maxNeighborColor = 0;
                for (int j = 0; j < m; j++) {
                    if (i != j && G.adj_d[u][seed.cand[j]]) {
                        if (colors[j] > maxNeighborColor) maxNeighborColor = colors[j];
                    }
                }
                if (maxNeighborColor >= reqColors) {
                    filteredCand.push_back(u);
                }
            }
            seed.cand = std::move(filteredCand);
            m = static_cast<int>(seed.cand.size());
            if (m == 0) continue;

            numColors = computeSubColoring(m, seed.cand, colors);
            if (static_cast<int>(seed.s.size()) + numColors < k) continue;

            std::vector<int> P_local, I_local;
            compIS(m, seed.cand, P_local, I_local);

            // Independent set I in complement graph is a CLIQUE in G
            std::vector<int> cliqueI = seed.s;
            for (int idx : I_local) cliqueI.push_back(seed.cand[idx]);
            if (static_cast<int>(cliqueI.size()) > lb) {
                lb = static_cast<int>(cliqueI.size());
                updateChampion(cliqueI);
            }

            for (int idx : P_local) {
                seed.s.push_back(seed.cand[idx]);
            }

            std::vector<int> newCand;
            for (int u : seed.cand) {
                bool inP = false;
                for (int pIdx : P_local) {
                    if (u == seed.cand[pIdx]) { inP = true; break; }
                }
                if (inP) continue;

                bool adjAllP = true;
                for (int pIdx : P_local) {
                    if (!G.adj_d[u][seed.cand[pIdx]]) {
                        adjAllP = false;
                        break;
                    }
                }
                if (adjAllP) newCand.push_back(u);
            }
            seed.cand = std::move(newCand);

            if (static_cast<int>(seed.s.size()) > lb) {
                lb = static_cast<int>(seed.s.size());
                updateChampion(seed.s);
            }
            if (lb >= k) k = lb + 1;

            if (!seed.cand.empty() && static_cast<int>(seed.s.size() + seed.cand.size()) >= k) {
                seed.ub = static_cast<int>(seed.s.size()) + numColors;
                resS.push_back(seed);
            }
        }

        if (!resS.empty()) {
            int maxSeedUb = 0;
            for (const auto& seed : resS) {
                if (seed.ub > maxSeedUb) maxSeedUb = seed.ub;
            }
            ub = std::min(ub, maxSeedUb);
        }

        return resS;
    }

    // Algorithm 7: divSeed(G, S, k)
    std::vector<Seed> divSeed(std::vector<Seed>& S, int k, int& lb, int& ub) {
        std::vector<Seed> S_prime;

        for (const auto& seed : S) {
            nodes++;
            checkTimeout();
            if (timeOut) break;

            int m = static_cast<int>(seed.cand.size());
            if (m == 0 || static_cast<int>(seed.s.size() + m) < k) continue;

            std::vector<int> colors;
            int numColors = computeSubColoring(m, seed.cand, colors);
            if (static_cast<int>(seed.s.size()) + numColors < k) continue;

            int reqColor = k - static_cast<int>(seed.s.size());
            std::vector<int> V1;
            for (int i = 0; i < m; i++) {
                if (colors[i] >= reqColor) {
                    V1.push_back(seed.cand[i]);
                }
            }

            for (int v : V1) {
                checkTimeout();
                if (timeOut) break;

                std::vector<int> s_prime = seed.s;
                s_prime.push_back(v);

                std::vector<int> cand_prime;
                for (int u : seed.cand) {
                    if (u != v && G.adj_d[v][u]) {
                        cand_prime.push_back(u);
                    }
                }

                if (static_cast<int>(s_prime.size() + cand_prime.size()) < k) continue;

                int m_prime = static_cast<int>(cand_prime.size());
                std::vector<int> primeCore = computeSubCoreNumbers(m_prime, cand_prime);
                int reqCorePrime = k - static_cast<int>(s_prime.size()) - 1;
                std::vector<int> filteredCandPrime;
                for (int i = 0; i < m_prime; i++) {
                    if (primeCore[i] >= reqCorePrime) {
                        filteredCandPrime.push_back(cand_prime[i]);
                    }
                }

                int szPrime = static_cast<int>(filteredCandPrime.size());
                std::vector<int> F;
                std::vector<int> remPrime;
                for (int i = 0; i < szPrime; i++) {
                    int u = filteredCandPrime[i];
                    int d = 0;
                    for (int j = 0; j < szPrime; j++) {
                        if (i != j && G.adj_d[u][filteredCandPrime[j]]) d++;
                    }
                    if (d == szPrime - 1) {
                        F.push_back(u);
                    } else {
                        remPrime.push_back(u);
                    }
                }

                s_prime.insert(s_prime.end(), F.begin(), F.end());
                cand_prime = std::move(remPrime);

                if (static_cast<int>(s_prime.size()) > lb) {
                    lb = static_cast<int>(s_prime.size());
                    updateChampion(s_prime);
                }
                if (lb >= k) k = lb + 1;

                if (!cand_prime.empty() && static_cast<int>(s_prime.size() + cand_prime.size()) >= k) {
                    S_prime.emplace_back(s_prime, cand_prime, seed.ub);
                }
            }
        }

        return tciSeed(S_prime, k, lb, ub);
    }

public:
    explicit RMCSolver(const Graph& g, double limitSec = 3600.0, double c = 0.2, int th = 3)
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false),
          timeLimit(limitSec), cParam(c), theta(th),
          rng(std::random_device{}()), opsCounter(0) {}

    // Algorithm 2: RMC(G)
    std::vector<int> solve() {
        nodes = 0;
        maxSize = 0;
        bestSol.clear();
        timeOut = false;
        opsCounter = 0;

        startTime = std::chrono::high_resolution_clock::now();

        int lb = 0, ub = 0;
        if (initBounds(lb, ub)) {
            auto end = std::chrono::high_resolution_clock::now();
            executionTime = std::chrono::duration<double, std::milli>(end - startTime).count();
            return bestSol;
        }

        std::vector<Seed> S;

        // Binary Search Phase (Algorithm 2, lines 3-17)
        while (ub - lb >= theta && !timeOut) {
            checkTimeout();
            if (timeOut) break;

            int wt = (lb + ub) / 2;
            int prevLb = lb, prevUb = ub;

            S = scSeed(wt, lb, ub);

            if (S.empty()) {
                if (lb >= wt) break;
                ub = std::min(ub, wt - 1);
            } else {
                S = tciSeed(S, wt, lb, ub);
                if (S.empty()) {
                    if (lb >= wt) break;
                    ub = std::min(ub, wt - 1);
                } else if (lb == prevLb && ub == prevUb) {
                    break;
                }
            }
        }

        // Terminal divSeed Phase (Algorithm 2, lines 18-23)
        if (!timeOut && lb < ub) {
            int wt = lb + 1;
            S = scSeed(wt, lb, ub);

            if (!S.empty()) {
                S = tciSeed(S, wt, lb, ub);

                while (!S.empty() && !timeOut && lb < ub) {
                    checkTimeout();
                    if (timeOut) break;
                    S = divSeed(S, lb + 1, lb, ub);
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        executionTime = std::chrono::duration<double, std::milli>(end - startTime).count();
        return bestSol;
    }

    bool hastimeOut() const { return timeOut; }
    long long getNodesExplored() const { return nodes; }
    double getexecutionTime() const { return executionTime; }
    int getMaxSize() const { return maxSize; }
};
#pragma once
#include <vector>
#include <chrono>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include "MCOSolver.hpp"   // for Graph

class BBMCSolver {
public:
    // Fixed-capacity or dynamically sized 64-bit word BitSet
    struct BitSet {
        int nWords;
        std::vector<uint64_t> words;

        BitSet() : nWords(0) {}
        explicit BitSet(int n) {
            init(n);
        }

        void init(int n) {
            nWords = (n + 63) / 64;
            words.assign(nWords, 0ULL);
        }

        inline void setBit(int i) {
            words[i >> 6] |= (1ULL << (i & 63));
        }

        inline void clearBit(int i) {
            words[i >> 6] &= ~(1ULL << (i & 63));
        }

        inline bool getBit(int i) const {
            return (words[i >> 6] & (1ULL << (i & 63))) != 0ULL;
        }

        inline void setAll(int n) {
            for (int i = 0; i < nWords; i++) words[i] = ~0ULL;
            if (n % 64 != 0) {
                words[nWords - 1] = (1ULL << (n % 64)) - 1ULL;
            }
        }

        inline void clearAll() {
            std::fill(words.begin(), words.end(), 0ULL);
        }

        inline bool isEmpty() const {
            for (int i = 0; i < nWords; i++) {
                if (words[i] != 0ULL) return false;
            }
            return true;
        }

        inline int popcount() const {
            int cnt = 0;
            for (int i = 0; i < nWords; i++) {
                cnt += __builtin_popcountll(words[i]);
            }
            return cnt;
        }

        inline void andInto(const BitSet& a, const BitSet& b) {
            for (int i = 0; i < nWords; i++) {
                words[i] = a.words[i] & b.words[i];
            }
        }

        inline void andWith(const BitSet& a) {
            for (int i = 0; i < nWords; i++) {
                words[i] &= a.words[i];
            }
        }

        inline int nextSetBit(int startBit = 0) const {
            int w = startBit >> 6;
            if (w >= nWords) return -1;
            uint64_t word = words[w] & (~0ULL << (startBit & 63));
            while (true) {
                if (word != 0ULL) {
                    return (w << 6) + __builtin_ctzll(word);
                }
                w++;
                if (w >= nWords) break;
                word = words[w];
            }
            return -1;
        }
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

    // Bit-board representations
    std::vector<int> origId;           // maps internal bit index -> original vertex id
    std::vector<BitSet> N;             // neighbourhood bitsets
    std::vector<BitSet> invN;          // inverse neighbourhood bitsets

    void saveSol(const BitSet& C) {
        bestSol.clear();
        for (int i = 0; i < G.n; i++) {
            if (C.getBit(i)) {
                bestSol.push_back(origId[i]);
            }
        }
        maxSize = bestSol.size();
    }

    // Minimum-width / Smallest-Last degeneracy ordering (Style 2 in Prosser's study)
    void orderVertices() {
        int n = G.n;
        std::vector<int> curDeg(G.degree);
        std::vector<char> removed(n, 0);
        origId.assign(n, 0);

        // Compute smallest-last ordering
        for (int pos = n - 1; pos >= 0; pos--) {
            int best = -1;
            for (int v = 0; v < n; v++) {
                if (removed[v]) continue;
                if (best == -1 || curDeg[v] < curDeg[best]) {
                    best = v;
                }
            }
            origId[pos] = best;
            removed[best] = 1;
            for (int u : G.adj_l[best]) {
                if (!removed[u]) curDeg[u]--;
            }
        }

        // Initialize BitBoard representations
        N.resize(n);
        invN.resize(n);
        for (int i = 0; i < n; i++) {
            N[i].init(n);
            invN[i].init(n);
        }

        for (int i = 0; i < n; i++) {
            int u = origId[i];
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                int v = origId[j];
                if (G.adj_d[u][v]) {
                    N[i].setBit(j);
                } else {
                    invN[i].setBit(j);
                }
            }
        }
    }

    // Class-based bit-parallel coloring (BBColour from San Segundo et al.)
    void BBColour(const BitSet& P, std::vector<int>& U, std::vector<int>& colour) {
        BitSet copyP = P;
        int colourClass = 0;
        int idx = 0;

        BitSet Q(G.n);

        while (!copyP.isEmpty()) {
            colourClass++;
            Q = copyP;
            while (!Q.isEmpty()) {
                int v = Q.nextSetBit(0);
                copyP.clearBit(v);
                Q.clearBit(v);
                Q.andWith(invN[v]);

                U[idx] = v;
                colour[idx] = colourClass;
                idx++;
            }
        }
    }

    // Recursive search procedure
    void BBMaxClique(BitSet& C, const BitSet& P, int cSize) {
        if (timeOut) return;

        if ((nodes & 4095) == 0) {
            auto current = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = current - startTime;
            if (elapsed.count() >= timeLimit) {
                timeOut = true;
                return;
            }
        }

        nodes++;

        int m = P.popcount();
        if (m == 0) return;

        std::vector<int> U(m);
        std::vector<int> colour(m);
        BBColour(P, U, colour);

        BitSet candP = P;
        BitSet newP(G.n);

        for (int i = m - 1; i >= 0; i--) {
            if (timeOut) return;
            // Prune branches whose color bound cannot unseat current champion
            if (colour[i] + cSize <= maxSize) {
                return;
            }

            int v = U[i];
            C.setBit(v);
            newP.andInto(candP, N[v]);

            int newCard = newP.popcount();
            if (newCard == 0) {
                if (cSize + 1 > maxSize) {
                    saveSol(C);
                }
            } else {
                BBMaxClique(C, newP, cSize + 1);
            }

            candP.clearBit(v);
            C.clearBit(v);
        }
    }

public:
    explicit BBMCSolver(const Graph& g, double limitSec = 3600.0)
        : G(g), maxSize(0), nodes(0), executionTime(0.0), timeOut(false),
          timeLimit(limitSec) {}

    std::vector<int> solve() {
        nodes = 0;
        maxSize = 0;
        bestSol.clear();
        timeOut = false;

        orderVertices();

        BitSet C(G.n);
        BitSet P(G.n);
        P.setAll(G.n);

        startTime = std::chrono::high_resolution_clock::now();

        BBMaxClique(C, P, 0);

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
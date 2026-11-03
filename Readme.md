# Exact Maximum Clique Solvers in C++

A comprehensive suite of exact algorithms and benchmarking drivers for solving the **Maximum Clique Problem (MCP)** on undirected graphs[cite: 1, 2, 3, 4, 5, 6]. The repository implements classic branch-and-bound strategies, degeneracy/minimum-width ordering, vertex coloring bounds, bitset parallelism, essential set reductions, and randomized seed-based pruning[cite: 1, 2, 4, 5, 6].

---

## Algorithms Implemented

| Algorithm | Engine File | Driver File | Core Techniques & Features |
| :--- | :--- | :--- | :--- |
| **MCO** | `MCOSolver.hpp`[cite: 3] | `main_MCO.cpp`[cite: 9] | Baseline recursive branch-and-bound search[cite: 3]. |
| **MCR** | `MCRSolver.hpp`[cite: 4] | `main_MCR.cpp`[cite: 10] | Degeneracy / minimum-width ordering with sequential greedy coloring upper bound[cite: 4]. |
| **MCS** | `MCSSolver.hpp`[cite: 5] | `main_MCS.cpp`[cite: 11] | Minimum-width ordering with exact vertex recoloring (`Re-NUMBER`) to tighten bounds[cite: 5]. |
| **BBMC** | `BBMCSolver.hpp`[cite: 1] | `main_BBMC.cpp`[cite: 7] | Bit-parallel branch-and-bound using a custom 64-bit `BitSet` bitboard representation[cite: 1]. |
| **MCF / MCF+ES** | `MCFSolver.hpp`[cite: 2] | `main_MCF.cpp`[cite: 8] | Max Clique Fast with degree/color bounds, universal vertex contraction, and optional Essential Set (ES) reduction[cite: 2, 4, 6]. |
| **RMC** | `RMCSolver.hpp`[cite: 6] | `main_RMC.cpp`[cite: 12] | Randomized framework using $k$-core decomposition, triangle sampling, truss filtering, and seed diversification[cite: 6]. |

---

## Project Structure

```text
.
├── BBMCSolver.hpp       # Bit-parallel branch-and-bound solver[cite: 1]
├── MCFSolver.hpp        # MCF engine with Essential Set reduction[cite: 2]
├── MCOSolver.hpp        # Baseline graph loader and solver[cite: 3]
├── MCRSolver.hpp        # Minimum-width ordering + greedy coloring solver[cite: 4]
├── MCSSolver.hpp        # Recoloring branch-and-bound solver[cite: 5]
├── RMCSolver.hpp        # Randomized seed-based reduction solver[cite: 6]
├── main_BBMC.cpp        # Driver for BBMC[cite: 7]
├── main_MCF.cpp         # Driver for MCF / MCF+ES[cite: 8]
├── main_MCO.cpp         # Driver for MCO[cite: 9]
├── main_MCR.cpp         # Driver for MCR[cite: 10]
├── main_MCS.cpp         # Driver for MCS[cite: 11]
├── main_RMC.cpp         # Driver for RMC[cite: 12]
├── brock200-1.mtx       # DIMACS / Benchmark graph sample[cite: 13]
├── frb40-19-1.mtx       # MatrixMarket benchmark graph sample[cite: 14]
├── frb59-26-5.mtx       # MatrixMarket benchmark graph sample[cite: 15]
└── README.md
```

---

## Input File Formats

The parser in `MCOSolver.hpp` automatically recognizes two formats[cite: 3]:

1. **DIMACS Graph Format (`.clq`, `.dimacs`)**[cite: 3]:
   ```text
   p edge <num_vertices> <num_edges>
   e <u_1> <v_1>
   e <u_2> <v_2>
   ```
2. **MatrixMarket Coordinate Format (`.mtx`)**[cite: 3, 14]:
   ```text
   %%MatrixMarket matrix coordinate pattern symmetric
   <num_vertices> <num_vertices> <num_edges>
   <u_1> <v_1>
   <u_2> <v_2>
   ```

*Note: Graph vertex numbering in files is 1-based and internally converted to 0-based indexing[cite: 3].*

---

## Compilation

The drivers use Windows API functions (`<windows.h>`, `<psapi.h>`) to monitor peak memory consumption, which requires linking `psapi`[cite: 7, 8, 9, 10, 11, 12].

### Windows (GCC / MinGW)

```bash
g++ -O3 -std=c++17 main_BBMC.cpp -o bbmc.exe -lpsapi
g++ -O3 -std=c++17 main_MCF.cpp -o mcf.exe -lpsapi
g++ -O3 -std=c++17 main_MCR.cpp -o mcr.exe -lpsapi
g++ -O3 -std=c++17 main_MCS.cpp -o mcs.exe -lpsapi
g++ -O3 -std=c++17 main_RMC.cpp -o rmc.exe -lpsapi
g++ -O3 -std=c++17 main_MCO.cpp -o mco.exe -lpsapi
```

### Linux / POSIX Compatibility

To compile on Linux, uncomment the `getrusage` blocks and comment out the `<psapi.h>` logic in the `main_*.cpp` files[cite: 7, 8, 9, 10, 11, 12], then compile with:

```bash
g++ -O3 -std=c++17 main_BBMC.cpp -o bbmc
```

---

## Execution & CLI Parameters

### BBMC
```bash
./bbmc.exe <graph_file> [timeout_seconds]
```
* Example: `./bbmc.exe brock200-1.mtx 1800`[cite: 7]

### MCF / MCF+ES
```bash
./mcf.exe <graph_file> [timeout_seconds] [enable_es (0/1)]
```
* `enable_es`: `1` runs MCF with Essential Set preprocessing (default); `0` runs basic MCF[cite: 8].
* Example: `./mcf.exe frb40-19-1.mtx 1800 1`

### MCR
```bash
./mcr.exe <graph_file> [timeout_seconds]
```
* Example: `./mcr.exe brock200-1.mtx 1800`[cite: 10]

### MCS
```bash
./mcs.exe <graph_file> [timeout_seconds] [reorder_depth]
```
* `reorder_depth`: Depth cutoff to re-sort subproblems by minimum width (default: `2`, `0` turns it off)[cite: 5, 11].
* Example: `./mcs.exe brock200-1.mtx 1800 2`[cite: 11]

### RMC
```bash
./rmc.exe <graph_file> [timeout_seconds] [c_param] [theta]
```
* `c_param`: Sampling parameter $c$ (default: `0.2`)[cite: 6, 12].
* `theta`: Stopping gap threshold $\theta$ between upper and lower bounds (default: `3`)[cite: 6, 12].
* Example: `./rmc.exe frb59-26-5.mtx 1800 0.2 3`

### MCO
```bash
./mco.exe <graph_file> [timeout_seconds]
```
* Example: `./mco.exe brock200-1.mtx 1800`[cite: 9]

---

## Output Metrics

Each execution prints standard benchmarking statistics[cite: 7, 8, 9, 10, 11, 12]:
- **Status**: `Optimal` or `Time out`[cite: 7, 8, 9, 10, 11, 12]
- **MaxClique Size**: The maximum clique cardinality $\omega(G)$ found[cite: 7, 8, 9, 10, 11, 12]
- **Nodes Explored**: Total search tree states visited[cite: 7, 8, 9, 10, 11, 12]
- **Execution Time**: High-resolution wall-clock duration in ms and seconds[cite: 7, 8, 9, 10, 11, 12]
- **Peak Memory Used**: Maximum resident set size in KB[cite: 7, 8, 9, 10, 11, 12]
- **Clique Vertices**: Space-separated list of 1-based vertex labels in the optimal clique[cite: 7, 8, 9, 10, 11, 12]
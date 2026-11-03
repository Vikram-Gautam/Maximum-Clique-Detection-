#include <iostream>
#include <iomanip>
//#include <sys/resource.h>
#include "MCSSolver.hpp"
#include <windows.h>
#include <psapi.h>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <graph file> [timeout seconds] [reorder depth]" << endl;
        return 1;
    }

    string filePath = argv[1];
    double timeoutSec = 1800.0;
    int reorderDepth = 2; // default from MCSSolver; 0 turns re-sorting off
    if (argc >= 3)
    {
        timeoutSec = stod(argv[2]);
    }
    if (argc >= 4)
    {
        reorderDepth = stoi(argv[3]);
    }

    Graph g;
    cout << "Loading " << filePath << " " << flush;
    if (!g.loadFromFile(filePath))
    {
        return 1;
    }
    cout << "vertices:" << g.n << " Edges:" << g.m << endl;

    MCSSolver solver(g, timeoutSec, reorderDepth);
    vector<int> maxClique = solver.solve();

    // struct rusage usage;
    // getrusage(RUSAGE_SELF, &usage);

    cout << "MCS" << endl;
    cout << "Status:" << (solver.hastimeOut() ? "Time out" : "Optimal") << endl;
    cout << "MaxClique Size:" << solver.getMaxSize() << endl;
    cout << "Nodes Explored:" << solver.getNodesExplored() << endl;
    cout << "Execution Time:" << fixed << setprecision(3)
         << solver.getexecutionTime() << " ms (" << solver.getexecutionTime() / 1000.0 << " s)" << endl;
    // cout << "Peak Memory Used:" << usage.ru_maxrss << " KB" << endl;
       PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
                cout << "Peak Memory Used:" << pmc.PeakWorkingSetSize / 1024 << " KB" << endl;
        }

    cout << "Clique Vertices:";
    for (int v : maxClique)
    {
        cout << " " << (v + 1);
    }
    cout << endl;

    return 0;
}
#include <iostream>
#include <iomanip>
//#include <sys/resource.h>

#include <windows.h>
#include <psapi.h>
#include "MCOSolver.hpp"
using namespace std;

int main(int argc, char *argv[])
{

    string filePath = argv[1];
    double timeoutSec = 1800.0;
    if (argc >= 3)
    {
        timeoutSec = stod(argv[2]);
    }

    Graph g;
    cout << "Loading" << filePath << " " << flush;
    if (!g.loadFromFile(filePath))
    {
        return 1;
    }
    cout << "vertices:" << g.n << "Edges:" << g.m << endl;

    MCOSolver solver(g, timeoutSec);
    vector<int> maxClique = solver.solve();

    // struct rusage usage;
    // getrusage(RUSAGE_SELF, &usage);

    cout << "MCO" << endl;
    cout << "Status:" << (solver.hastimeOut() ? "Time out" : "Optimal") << endl;
    cout << "MaxClique Size:" << solver.getMaxSize() << endl;
    cout << "Nodes Explored:" << solver.getNodesExplored() << endl;
    cout << "Execution Time:" << fixed << setprecision(3) << solver.getexecutionTime() << "ms" << solver.getexecutionTime() / 1000.0 << " s" << endl;
    // cout << "Peak Memory Used:" << usage.ru_maxrss << " KB" << endl;
       PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
                cout << "Peak Memory Used:" << pmc.PeakWorkingSetSize / 1024 << " KB" << endl;
        }

    cout << "Clique Vertices:";
    for (int v : maxClique)
    {
        cout << (v + 1) << " ";
    }
    cout << endl;

    return 0;
}
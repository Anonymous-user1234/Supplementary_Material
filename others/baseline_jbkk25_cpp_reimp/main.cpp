#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

#include "New_BP.cpp"
#include "extract_circuit_formal.cpp"
#include "BPD_opti.cpp"


#ifndef HH
#define HH (23)
#endif

#ifndef MC
#define MC (0)
#endif


#ifndef CACUB
#define CACUB (999)
#endif

void print_vecnum(const vector<int>& data) {
    if (data.empty()) {
        cout << "Empty vector" << endl;
        return;
    }
    
    map<int, int> freq;
    for (int value : data) {
        freq[value]++;
    }
    
    cout << "Total items: " << data.size() << ", Unique values:" << freq.size() << endl;
    cout << "-----------------" << endl;
    for (const auto& pair : freq) {
        cout << setw(5) << pair.first << " : " << setw(5) << pair.second 
                  << " (" << fixed << setprecision(1) 
                  << (static_cast<double>(pair.second) / data.size() * 100) << "%)" << endl;
    }
    cout << "==================" << endl;
}


int main(int argc, char**argv) {

    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <matrix_file>" << endl;
        return 1;
    }
    string filename = argv[1];
    int Hdepth = HH;


    int best_xor = INT_MAX;
    string best_log = "";
    vector<int> xors_record;

    
    //step1: read and init circuit
    PRINTF_STAMP("BPD Optimize filename: %s, depth: %d\n\n",filename.c_str(), Hdepth);
    
    
    CirParser circuit(filename);
    
    circuit.cir_formal(filename);
    circuit.extract_info(filename);
    
    // if use raandom circuit
    if(MC) circuit.modify_cir();
    
    circuit.init_bit();
    
    circuit.print_cir();

    
    while(true){
        //step2: BPD optimize
        BPDOptimizer BPDopti(CACUB, MC, filename, circuit.getn(), circuit.getm(), circuit.getk(), circuit.getS(), circuit.getD(), circuit.getY(), circuit.getSname(), circuit.getNLs(), circuit.getNOTs());
    
        auto result = BPDopti.optimize(Hdepth);
        if(result.first < best_xor){
            best_xor = result.first;
            best_log = result.second;
        }

        xors_record.push_back(result.first);
    
        PRINTF_STAMP("Caculate finish~~~ (best xors: %d -- %s) \n current xors: %d, log filename: %s\n\n", best_xor, best_log.c_str(), result.first, result.second.c_str());
        print_vecnum(xors_record);

    }


    return 0;
}
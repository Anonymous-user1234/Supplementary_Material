#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <utility>
#include <thread>
#include <iostream>
#include <atomic>
#include <unistd.h>
#include <sys/syscall.h>
#include <map>
#include <set>

void print_unordered_map(const unordered_map<uint64_t, int>& S_map) {
    cout << "\nmap contents:" << endl;
    for (const auto& pair : S_map) {
        cout << bitset<40>(pair.first) << ": " << pair.second << endl;
    }
}



void print_vecunordered_map(const vector<unordered_map<uint64_t, int>>& map_collection, 
                       const string& name = "Map_Collection") {
    cout << "=== " << name << " ===" << endl;
    cout << "Total levels: " << map_collection.size() << endl;
    
    for (size_t i = 0; i < map_collection.size(); i++) {
        cout << "Level " << i << ": " << map_collection[i].size() << " entries" << endl;
    }
    cout << "==================" << endl << endl;
}


void print_map_fre(const unordered_map<uint64_t, int>& input_map, 
                        const string& map_name = "Map") {
    
    map<int, int> frequency;
    for (const auto& pair : input_map) {
        frequency[pair.second]++;
    }
    for (const auto& freq_pair : frequency) {
        cout << freq_pair.first << ": " << freq_pair.second << endl;
    }
}

void print_vec(const vector<int>& x){
    for(auto &ele: x){
        cout << ele << " ";  
    } 
    cout << endl;
}


pid_t get_threadid() {
    return syscall(SYS_gettid);
}



class BPDOptimizer {
private:
    string Filename;
    string LOG_filename;
    int n, m, k;
    vector<uint64_t> S;
    uint64_t Smask;
    vector<uint64_t> D;     
    vector<uint64_t> Y;
    // vector<string> Sname;
    vector<string> NLs;    
    vector<string> NOTs;
    
    vector<int> Dist;
    vector<int> prev_Dist;
    vector<int> HY;  
    int depth;
    
    vector<unordered_map<uint64_t, int>> S_XORs_2PD;
    vector<unordered_map<uint64_t, int>> S_XORs_2PD_for_new_pair;


    vector<op_t> ops;
    unordered_map<uint64_t, int> S_map;

    
    mt19937 rng;
    
    ofstream log_file;

    bool modify_f;

    int cacu_b;

    
    
public:
    BPDOptimizer(int caculate_bound, bool modify_flag, string filename, int input_size, int output_size, int nonlinear_gates, 
                 const vector<uint64_t>& initial_S,
                 const vector<uint64_t>& initial_D,
                 const vector<uint64_t>& targets_Y,
                 const vector<string>& signal_names,
                 const vector<string>& nonlinear_ops,
                 const vector<string>& not_vars)
        : cacu_b(caculate_bound), modify_f(modify_flag), n(input_size), m(output_size), k(nonlinear_gates),
          S(initial_S), D(initial_D), Y(targets_Y),
          NLs(nonlinear_ops), NOTs(not_vars) , rng(random_device{}()){
        
        Filename = filename; 
        
        S_XORs_2PD.resize(5);
        S_XORs_2PD_for_new_pair.resize(5);

        S_map.clear();
        Smask = 0;
        for (size_t i = 0; i < S.size(); i++) {
            S_map[S[i]] = i;
            Smask |= S[i];
        }

    }


    
    pair<int, string> optimize(int H = 999, string log_filename = "log") {

        if(modify_f) log_filename += "_modify_";
        init_log(log_filename + "_" + Filename + "_core" + to_string(get_threadid()));
        init_opti(H);

        
        int t = 0; 
        uint64_t best_choose = 0;
        
        while (any_of(Dist.begin(), Dist.end(), [](int d) { return d > 0; })) {
            
            string add_info = "";
            
            if (auto [g_index, i1, i2] = find_NLGate(); g_index >= 0) {// Step 1: Easy case II 
                best_choose = add_NLGate(g_index, i1, i2);

                add_info = "NO1. Easy case II";

            }else if (min_d999()) {// Step 2: process dist999
                int index = get_999index();
                best_choose = add_d999T(index);
                
                add_info = "NO2. Dist999";
            }
            else if (!has_d1T()) { // Step 3: BPD alg
                best_choose = add_BPD_Opt(H);

                add_info = "NO3. BPD opti";
            }
            // Step 4: BP easy case
            else {
                best_choose = add_BPD_EC();

                add_info = "NO4. BPD easy case";
            }
            
            t++;
            update_state(best_choose);
            log_optstep(t, add_info);


            if(t > cacu_b + k){
                PRINTF_STAMP("%d reached bound, cut caculate!\n\n", t);
                return {999, LOG_filename};
            }
            
        }
        
        
        
        log_ops();
        cleanup();

        write_result();
        assert(*max_element(D.begin(), D.end()) <= depth);
        assert(depth <= H);

        return {S.size() - n - k, LOG_filename};
        
    }

private:

    // init HY, depth, Dist, XOR table
    void init_opti(int H) {
        
        // 1. caculate HY 
        auto result =limit_depths(H);
        depth = result.first;
        HY = result.second;
        log_file << "Max depth: " << depth << endl;
        log_file << "target depth bound: ";
        for (size_t i = 0; i < HY.size(); i++) {
            log_file << HY[i] << ", ";
        }
        log_file << "\n" << endl;
        
        // 2. init dist
        Dist = init_dist();
        log_file << "Dist: ";
        for (size_t i = 0; i < Dist.size(); i++) {
            log_file << Dist[i] << ", ";
        }
        log_file << "\n" << endl;

        prev_Dist = Dist;

        
        // 3. init XOR table
        init_SXORs();
        
        log_file << "S_XORs_2PD size: " << S_XORs_2PD.size() << endl;
        for (size_t i = 0; i < S_XORs_2PD.size(); i++) {
            log_file << "  level " << i << " size: " << S_XORs_2PD[i].size() << endl;
        }
        log_file << endl;
        
        PRINTF_STAMP("init HY, dist. All init finish!\n\n");
    }

    bool is_inS(uint64_t signal) {
        return S_map.find(signal) != S_map.end();
    }

    bool is_inY(uint64_t signal) const {
        return find(Y.begin(), Y.end(), signal) != Y.end();
    }

    tuple<int, int, int> find_NLGate() {

        for (int i = 0; i < k; i++) {
            if (is_inS(1ULL << (n + i))) {
                continue; // g[i] already exists in S, next find
            }
            
            bool input1 = false; bool input2 = false;
            int i1 = -1; int i2 = -1;

            for (size_t j = 0; j < S.size(); j++) {
                // cout << S[j] << " " << D[j] << endl;
                // cout << Y[2 * i] << " " << HY[2 * i] << endl;
                if (S[j] == Y[2 * i] && D[j] <= HY[2 * i]) {
                    input1 = true;
                    i1 = j;                   
                }
                if (S[j] == Y[2 * i + 1] && D[j] <= HY[2 * i + 1]) {
                    input2 = true;
                    i2 = j;
                }
            }
            
            if (input1 && input2) {
                return make_tuple(i, i1, i2);
            }
        }
        return make_tuple(-1, -1, -1);
    }

    uint64_t add_NLGate(int g_index, int i1, int i2) {
        uint64_t g = 1ULL << (n + g_index);

        int new_index = S.size();

        S.push_back(g);
        Smask |= g;
        D.push_back(max(D[i1], D[i2]) + 1);
        // Sname.push_back("t["+ to_string(new_index - n) + "]");
        S_map[g] = new_index;
        ops.emplace_back(op_t{new_index, {i1, i2}, "g", g_index});

        PRINTF_STAMP("%dth Add: %s\n\n", int(S.size() - n), signalToString(S.back()).c_str());
        return g;
    }

    void init_log(const string& log_filename) {
        int i = 0;
        string final_filename;
        
        while (true) {
            final_filename = "log/" + log_filename + "_" + 
                               (i < 10 ? "00" : (i < 100 ? "0" : "")) + 
                               to_string(i) + ".txt";
            ifstream test_file(final_filename);
            if (!test_file.good()) {
                break;
            }
            i++;
        }

        LOG_filename = final_filename;

        log_file.open(final_filename);
        if (!log_file.is_open()) {
            cerr << "warning: can't create file " << final_filename << endl;
            return;
        }
        cout << "\n\n";
        PRINTF_STAMP("Log file: %s \n\n", final_filename.c_str() );
        
        log_file << "===== BPD optimize =====" << endl;
        log_file << "Begin time: " << getCurrentSystemTime() << endl;
        log_file << "Log file: " << final_filename << endl;
        log_file << "========================\n" << endl;
    }

    vector<int> init_dist() {
        vector<int> dist;
        for (size_t i = 0; i < Y.size(); i++) {
            if (Y[i] >= (1ULL << n)) {
                dist.push_back(999);
            } else {
                dist.push_back(__builtin_popcountll(Y[i]) - 1);
            }
        }
        return dist;
    }

    void init_SXORs() {
        
        // level0: 0
        S_XORs_2PD[0][0] = 0;
        
        // level1 
        for (size_t i = 0; i < S.size(); i++) {
            S_XORs_2PD[1][S[i]] = (1 << D[i]);  
        }
        
        // level2
        for (size_t i = 0; i < S.size(); i++) {
            for (size_t j = i + 1; j < S.size(); j++) {
                uint64_t xor_val = S[i] ^ S[j];
                int cost = (1 << D[i]) + (1 << D[j]);
                S_XORs_2PD[2][xor_val] = cost;// D[xor_val] = cost(sigma 2^D)
            }
        }
        
        // level3
        for (size_t i = 0; i < S.size(); i++) {
            for (size_t j = i + 1; j < S.size(); j++) {
                for (size_t k = j + 1; k < S.size(); k++) {
                    uint64_t xor_val = S[i] ^ S[j] ^ S[k];
                    int cost = (1 << D[i]) + (1 << D[j]) + (1 << D[k]);
                    S_XORs_2PD[3][xor_val] = cost;
                }
            }
        }
        
        // level4
        for (size_t i = 0; i < S.size(); i++) {
            for (size_t j = i + 1; j < S.size(); j++) {
                for (size_t k = j + 1; k < S.size(); k++) {
                    for (size_t l = k + 1; l < S.size(); l++) {
                        uint64_t xor_val = S[i] ^ S[j] ^ S[k] ^ S[l];
                        int cost = (1 << D[i]) + (1 << D[j]) + (1 << D[k]) + (1 << D[l]);
                        S_XORs_2PD[4][xor_val] = cost;
                    }
                }
            }
        }
    }

    pair<int, vector<int>> limit_depths(int H) {
        // X： all Y[i]'s impletation 
        vector<vector<vector<uint64_t>>> X(Y.size());
        for (size_t i = 0; i < Y.size(); i++) {
            uint64_t y = Y[i];
            vector<uint64_t> xors;
            // change r,y to vars
            for (int j = 0; j < 64; j++) {
                if (y & (1ULL << j)) {
                    xors.push_back(1ULL << j);
                }
            }
            X[i].push_back(xors);//xors: x,g
        }
        
        // init depth
        unordered_map<uint64_t, int> D;
        // D: all vars depth
        for (int i = 0; i < n; i++) {
            D[1ULL << i] = 0;  //init input var depth = 0
        }
        // init r var 
        for (int i = 0; i < k; i++) {
            //process r[2*i]
            uint64_t r0 = Y[2 * i];
            int sum0 = 0;
            for (uint64_t x : X[2 * i][0]) {
                sum0 += (1 << D[x]);
            }
            D[r0] = ceil(log2(sum0));
            
            // process r[2*i+1]
            uint64_t r1 = Y[2 * i + 1];
            int sum1 = 0;
            for (uint64_t x : X[2 * i + 1][0]) {
                sum1 += (1 << D[x]);
            }
            D[r1] = ceil(log2(sum1));
            
            // process g[i]
            uint64_t g = 1ULL << (n + i);
            D[g] = max(D[Y[2 * i]], D[Y[2 * i + 1]]) + 1;
        }
        // init y var
        for (int i = 0; i < m; i++) {
            uint64_t y = Y[2 * k + i];
            int sum = 0;
            for (uint64_t x : X[2 * k + i][0]) {
                sum += (1 << D[x]);
            }
            D[y] = ceil(log2(sum));
        }
        // print_unordered_map(D);
        
        // let DY = D
        vector<int> DY(Y.size());
        for (size_t i = 0; i < Y.size(); i++) {
            DY[i] = D[Y[i]];
        }
        // first check
        auto [d, _] = check_max_depth(D, X);
        
        D.clear();
        for (int i = 0; i < n; i++) {
            D[1ULL << i] = 0;
        }
        for (int i = 0; i < 2 * k; i++) {// y = H ,r by add
            if (DY[i] <= 1) {
                D[Y[i]] = DY[i]; 
                continue;
            }
            if (D.find(Y[i]) != D.end()) {
                continue;
            }
            while (true) {
                D[Y[i]] = DY[i];
                auto [current_d, new_DY] = check_max_depth(D, X);
                if (current_d > H) {
                    D[Y[i]] = new_DY[i] - 1;  
                    break;
                }
                DY[i]++;
            }
        }
        
        return check_max_depth(D, X);
    }

    pair<int, vector<int>> check_max_depth(const unordered_map<uint64_t, int>& D_based, const vector<vector<vector<uint64_t>>>& X) {
        
        unordered_map<uint64_t, int> D = D_based;
        
        for (int i = 0; i < n; i++) {
            uint64_t x_mask = 1ULL << i;
            if (D.find(x_mask) == D.end()) {
                D[x_mask] = 0;
            }
        }
        
        for (int i = 0; i < k; i++) {
            uint64_t r0 = Y[2 * i];
            if (D.find(r0) == D.end()) {
                int min_sum = INT_MAX;
                for (const auto& XOR : X[2 * i]) {
                    int sum = 0;
                    for (uint64_t x : XOR) {
                        if (D.find(x) != D.end()) {
                            sum += (1 << D[x]);
                        }
                    }
                    min_sum = min(min_sum, sum);
                }
                D[r0] = ceil(log2(min_sum));
            }
            
            uint64_t r1 = Y[2 * i + 1];
            if (D.find(r1) == D.end()) {
                int min_sum = INT_MAX;
                for (const auto& XOR : X[2 * i + 1]) {
                    int sum = 0;
                    for (uint64_t x : XOR) {
                        if (D.find(x) != D.end()) {
                            sum += (1 << D[x]);
                        }
                    }
                    min_sum = min(min_sum, sum);
                }
                D[r1] = ceil(log2(min_sum));
            }
            
            uint64_t g = 1ULL << (n + i);
            if (D.find(g) == D.end()) {
                D[g] = max(D[Y[2 * i]], D[Y[2 * i + 1]]) + 1;
            }
        }
        
        for (int i = 0; i < m; i++) {
            uint64_t y = Y[2 * k + i];
            if (D.find(y) == D.end()) {
                int min_sum = INT_MAX;
                for (const auto& XOR : X[2 * k + i]) {
                    int sum = 0;
                    for (uint64_t x : XOR) {
                        if (D.find(x) != D.end()) {
                            sum += (1 << D[x]);
                        }
                    }
                    min_sum = min(min_sum, sum);
                }
                D[y] = ceil(log2(min_sum));
            }
        }
        
        int max_depth = 0;
        vector<int> DY(Y.size());
        for (size_t i = 0; i < Y.size(); i++) {
            DY[i] = D[Y[i]];
            if (DY[i] > max_depth) {
                max_depth = DY[i];
            }
        }

        return {max_depth, DY};

    }

    bool min_d999() {
        int min_dist = INT_MAX;
        for (size_t i = 0; i < Dist.size(); i++) {
            if (Dist[i] > 0 && Dist[i] < 999) {
                return false;
            }
        }
        return true;
    }

    int get_999index() {
        for (size_t i = 0; i < Dist.size(); i++) {
            if (Dist[i] == 999) return i;
        }
        return -1;
    }

    void add_toS(const uint64_t& elem, const int& elem_d, const pair_t& sources, int y_index = -1) {
        int new_index = S.size();
        S.push_back(elem);
        Smask |= elem;
        assert(elem_d < 999); 
        D.push_back(elem_d);
        S_map[elem] = new_index;

        if(y_index == -1){
            // Sname.push_back("t[" + to_string(new_index - n) + "]");
            ops.emplace_back(op_t{new_index, sources});    
        }else{
            assert(y_index >= 0);
            // Sname.push_back("t[" + to_string(new_index - n) + "]");s
            ops.emplace_back(op_t{new_index, sources, "y", y_index});    

        }
        
        PRINTF_STAMP("%dth Add: %s\n\n", int(S.size() - n), signalToString(S.back()).c_str());
        
    }

    uint64_t add_d999T(int& t_index) {
     
        assert(t_index != -1);
        
        uint64_t y = Y[t_index];
        
        // Extract all vars from y
        unordered_map<uint64_t, int> XD;
        for (int j = 0; j < 64; j++) {
            if (y & (1ULL << j)) {
                uint64_t var = 1ULL << j;
                for (size_t k = 0; k < S.size(); k++) {
                    if (S[k] == var) {
                        XD[var] = D[k];
                        break;
                    }
                }
            }
        }
        
        vector<uint64_t> Xs;
        for (const auto& pair : XD) {
            Xs.push_back(pair.first);
        }
        
        uint64_t w;
        int w_depth;
        uint64_t v0, v1;
        
        while (true) {
            // Sort depth
            sort(Xs.begin(), Xs.end(), [&](uint64_t a, uint64_t b) {
                return XD[a] < XD[b];
            });
            
            v0 = Xs[0];
            v1 = Xs[1];
            w = v0 ^ v1; // the same as bound caculate method
            
            // Check if w is already in S
            bool w_in_S = false;
            for (size_t i = 0; i < S.size(); i++) {
                if (S[i] == w) {
                    w_in_S = true;
                    XD[w] = D[i];
                    break;
                }
            }
            
            if (!w_in_S) {
                w_depth = max(XD[v0], XD[v1]) + 1;
                break;
            } else {
                // Remove v0, v1 and add w
                Xs.erase(Xs.begin(), Xs.begin() + 2);
                Xs.push_back(w);
                
                // Check depth
                int sum_d = 0;
                for (uint64_t x : Xs) {
                    sum_d += (1 << XD[x]);
                }
                if (ceil(log2(sum_d)) > HY[t_index]) {
                    w_depth = max(XD[v0], XD[v1]) + 1;// should not happen
                    assert(0);
                    PRINTF_STAMP("Warning: The depth has exceed, the combination not finish, HY index --  %d\n\n", t_index);
                    break;
                }
            }
        }
        
        add_toS(w, w_depth, {S_map[v0], S_map[v1]});
        return w;

    }

    bool has_d1T() {
        for (int dist : Dist) {
            if (dist == 1) return true;
        }
        return false;
    }

    uint64_t add_BPD_Opt(int H) {

        // 1: Generate candidate and depth bound
        unordered_map<uint64_t, int> WD_depth; 
        unordered_map<uint64_t, pair<uint64_t, uint64_t>> WD_pairs;
        make_WD(H, WD_depth, WD_pairs);


        // 2: Filter candidates that can reduce distance
        std::vector<uint64_t> valid_candi;
        std::unordered_map<uint64_t, std::vector<int>> WD_Dist;

        std::vector<int> original_dist = Dist;

        int        best_999 = INT_MAX;
        long long  best_sum = LLONG_MAX;
        long long  best_norm = LLONG_MIN;

        for (const auto& kv : WD_depth) {
            const uint64_t w       = kv.first;
            const int      w_depth = kv.second;
            const uint64_t w_mask  = (Smask | w);

            std::vector<int> test_dist = Dist;

            int       cnt999 = 0;
            long long sum    = 0;
            long long norm   = 0;
            bool      reduces_dist = false;
            bool      pruned_early = false;

            for (size_t u = 0; u < Y.size(); ++u) {
                if (Dist[u] == 0) {
                    test_dist[u] = 0;
                } else {
                    if ( (w_mask & Y[u]) != Y[u] ) {
                        test_dist[u] = 999;
                        ++cnt999;
                    } else {
                        const int nd = update_distance_for_new_pair(w, w_depth, Y[u], HY[u], Dist[u]);
                        test_dist[u] = nd;
                        if (nd != original_dist[u]) reduces_dist = true;

                        if (nd == 999) {
                            ++cnt999;
                        } else {
                            sum  += nd;
                            norm += 1LL * nd * nd;
                        }
                    }
                }

                if (best_999 != INT_MAX) {
                    if (cnt999 > best_999) {
                        pruned_early = true;
                        break; 
                    }
                    if (cnt999 == best_999 && sum > best_sum) {
                        pruned_early = true;
                        break; 
                    }
                }
            }

            if (pruned_early || !reduces_dist) continue;

            valid_candi.push_back(w);
            WD_Dist[w] = std::move(test_dist);

            if (cnt999 < best_999 ||
                (cnt999 == best_999 && (sum < best_sum ||
                (sum == best_sum && norm > best_norm)))) {
                best_999 = cnt999;
                best_sum = sum;
                best_norm = norm;
            }
        }
                
        //3: Minimize 999 dist, Minimize sum of distances, Maximize Euclidean norm
        int min_sum = INT_MAX;
        int max_norm = 0;
        int min_999_count = INT_MAX;
        vector<uint64_t> best_candi;


        map<uint64_t, tuple<int, int, int>> candidates_with_sum;
        
        for (uint64_t w : valid_candi) {
            int sum = 0;
            double norm = 0;
            int count_999 = 0;

            for (int d : WD_Dist[w]) {
                if (d != 999){
                    norm += d * d;
                    sum += d;
                } 
                else count_999++;
            }
            
            candidates_with_sum[w] = {count_999, sum, norm};

            if(count_999 < min_999_count){
                min_sum = sum;
                max_norm = norm;
                min_999_count = count_999;
                best_candi = {w};
            }else if(count_999 == min_999_count){
                if(sum < min_sum || (sum == min_sum && norm > max_norm)){
                    min_sum = sum;
                    max_norm = norm;
                    best_candi = {w};
                }else if(sum == min_sum && norm == max_norm){
                    best_candi.push_back(w);
                }
            }
        }
        ana_candi(candidates_with_sum);
        
        // 4: Random selection
        uniform_int_distribution<int> dist(0, best_candi.size() - 1);
        int select_index = dist(rng);
        uint64_t selected_w = best_candi[select_index];
        log_file << "Candidate size: " << to_string(best_candi.size()) << ", Select index:" << to_string(select_index) << endl;
        add_toS(selected_w, WD_depth[selected_w], WD_pairs[selected_w]);
        
        return selected_w;
    }

    uint64_t add_BPD_EC() {
        // Find the first target with distance 1
        int target_index = -1;
        for (size_t i = 0; i < Dist.size(); i++) {
            if (Dist[i] == 1) {
                target_index = i;
                break;
            }
        }
        
        assert(target_index != -1);
        
        uint64_t y = Y[target_index];

        uint64_t src0 = 0, src1 = 0;
        int y_depth = -1;
        // use min depth
        vector<pair<set<uint64_t>, int>> ec_candi;
        for (size_t i = 0; i < S.size(); i++) {
            // uint64_t needed = S[i] ^ y;
            src0 = S[i];
            src1 = S[i] ^ y;
            if(is_inS(src1)){
                y_depth = max(D[S_map[src0]], D[S_map[src1]]) + 1;
                if(y_depth <= HY[target_index]){
                    pair<set<uint64_t>, int> src = {{src0, src1}, y_depth};
                    if(find(ec_candi.begin(),ec_candi.end(),src) == ec_candi.end()) ec_candi.push_back(src);
                }
            }
        }

        sort(ec_candi.begin(), ec_candi.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

        if(ec_candi.empty()) assert(0);

        print_eccandi(ec_candi);
        pair<set<uint64_t>, int> best_candi = ec_candi[0];
        src0 = *best_candi.first.begin();
        src1 = *best_candi.first.rbegin();
        y_depth = best_candi.second;
        log_file << "Easy case best choose: " << src0 << ", " << src1 << " -- " << y_depth << endl;


        assert(src0 != 0 && src1 != 0);
        y_depth = max(D[S_map[src0]], D[S_map[src1]]) + 1;
        assert(y_depth <= HY[target_index]);

        int y_index = -1;
        if(target_index >= 2*k) y_index = target_index - 2*k;
        assert(y_index < m);
        add_toS(y, y_depth, {S_map[src0], S_map[src1]}, y_index);

        return y;
        
    }

    // Update xor table and dist
    void update_state(const uint64_t& candidate) {

        // Update dist for all targets
        for (size_t i = 0; i < Y.size(); i++) {
            if (Dist[i] == 0) continue;
            Dist[i] = update_dist(i, candidate);
        }
        
        // Update S_XORs structures
        update_S_XORs_for_new_pair(candidate, D[S_map[candidate]]);
        update_SXORs();
            
        
    }



    
    pair<int, int> get_diststa() const {
        int count_999 = 0;
        int count_non_999 = 0;
        
        for (int d : Dist) {
            if (d == 999) count_999++;
            else count_non_999++;
        }
        
        return {count_non_999, count_999};
    }

    string signalToString(uint64_t signal) const {
        string binary_str;
        
        // only display n+k bit
        for (int i = (n + k - 1); i >= 0; i--) {
            if (signal & (1ULL << i)) {
                binary_str += '1';
            } else {
                binary_str += '0';
            }
        }
        
        return binary_str;
    }

    void log_optstep(int step, string add_info) {

        log_file << "\n"<< getCurrentSystemTime() << endl;
        log_file << step << "th Basis (S, D) : " << ", " 
                 << S.back() << "( " << signalToString(S.back()) <<" ), " << D.back() << endl;
        
        log_file  << add_info << endl;
        log_file << "New gate " << ops.back().type << "[" << ops.back().type_index << "] index: " << ops.back().dst << " = " << ops.back().src.first << " , " << ops.back().src.second << endl;
        auto [n999, y999] = get_diststa();
        log_file << "#999: " << n999 << ", #<999: " << y999 << ", Mask: " << signalToString(Smask) << endl; 
        log_file << "Dist : ";
        for (size_t i = 0; i < Dist.size(); i++) {
            log_file << Dist[i];
            if (i < Dist.size() - 1) log_file << ", ";
        }
        log_file << endl;

        log_file << "       ";
        for (size_t i = 0; i < Dist.size(); i++) {
            if (Dist[i] < prev_Dist[i]) {
                assert(prev_Dist[i] == 999 || Dist[i] + 1 == prev_Dist[i]);
                log_file << "_"; 
            } else {
                if(Dist[i] == 999) log_file << "   ";
                else  log_file << " ";
            }
            if (i < Dist.size() - 1) log_file << "  ";
        }
        log_file << endl;

        prev_Dist = Dist;

    
        log_file << endl;
    }


    void log_ops() {
        int xor_count = S.size() - n - k;
        log_file << "XOR count: " << xor_count << endl;
        log_file << "Circuit:" << endl;
        for(int i = 0; i < n; i++){
            log_file << "t[" << i << "]=x[" << i << "]" << endl;
        }
        for (const op_t& gate : ops) {
            if(gate.type == "g") log_file << "t[" << gate.dst << "]=t[" << gate.src.first << "]" << NLs[gate.type_index] <<"t["<< gate.src.second << "]" << endl;
            else if(gate.type == "y"){
                log_file << "t[" << gate.dst << "]=t[" << gate.src.first << "]^t["<< gate.src.second << "]" << endl;
                
                string not_add = "y[" + to_string(gate.type_index) + "]";
                if(find(NOTs.begin(), NOTs.end(), not_add) != NOTs.end()) log_file << "y[" << gate.type_index << "]=t[" << gate.dst << "]^1" << endl;
                else log_file << "y[" << gate.type_index << "]=t[" << gate.dst << "]" << endl;
            }
            else log_file << "t[" << gate.dst << "]=t[" << gate.src.first << "]^t["<< gate.src.second << "]" << endl;
        }
        log_file << "\n\nOptimization completed successfully!" << endl;
    }

    void cleanup() {
        // Clear the S_XORs structures
        S_XORs_2PD.clear();
        S_XORs_2PD.resize(5);
        S_XORs_2PD_for_new_pair.clear();
        S_XORs_2PD_for_new_pair.resize(5);
        
        log_file.close();
    }

    void make_WD(int H, unordered_map<uint64_t, int>& WD_depth, unordered_map<uint64_t, pair<uint64_t, uint64_t>>& WD_pairs) {
        WD_depth.clear();
        WD_pairs.clear();
        
        for (size_t i = 0; i < S.size(); i++) {
            for (size_t j = i + 1; j < S.size(); j++) {
                uint64_t w = S[i] ^ S[j];
                int d = max(D[i], D[j]) + 1;
                
                if(d >= H) continue;

                bool flag = false;
                for(int i = 0; i < Y.size(); i++){
                    if(Y[i] == w && d > HY[i]) flag = true;
                }
                if(flag) continue;
                
                // Check if w is already in S with better depth
                if(is_inS(w) && d >= D[S_map[w]]) continue;
                
                // be a candidate if w is new  or w with better depth
                if (WD_depth.find(w) == WD_depth.end() || d < WD_depth[w]) {
                    WD_depth[w] = d;
                    WD_pairs[w] = {i,j};
                }
            }
        }
        
    }

    void update_S_XORs_for_new_pair(uint64_t w, int d) {
        S_XORs_2PD_for_new_pair[0].clear();
        S_XORs_2PD_for_new_pair[1][w] = (1 << d);
        
        // Level 2: w XOR each element in level 1
        for (const auto& pair : S_XORs_2PD[1]) {
            uint64_t xor_val = w ^ pair.first;
            int cost = (1 << d) + pair.second;
            S_XORs_2PD_for_new_pair[2][xor_val] = cost;
        }
        
        // Level 3: w XOR each element in level 2
        for (const auto& pair : S_XORs_2PD[2]) {
            uint64_t xor_val = w ^ pair.first;
            int cost = (1 << d) + pair.second;
            S_XORs_2PD_for_new_pair[3][xor_val] = cost;
        }
        
        // Level 4: w XOR each element in level 3
        for (const auto& pair : S_XORs_2PD[3]) {
            uint64_t xor_val = w ^ pair.first;
            int cost = (1 << d) + pair.second;
            S_XORs_2PD_for_new_pair[4][xor_val] = cost;
        }
    }

    void update_SXORs() {
        for (int i = 1; i <= 4; i++) {
            for (const auto& pair : S_XORs_2PD_for_new_pair[i]) {
                uint64_t w = pair.first;
                int cost = pair.second;
                
                if (S_XORs_2PD[i].find(w) != S_XORs_2PD[i].end()) {
                    S_XORs_2PD[i][w] = min(S_XORs_2PD[i][w], cost);
                } else {
                    S_XORs_2PD[i][w] = cost;
                }
            }
        }
    }

    int update_dist(int u, const uint64_t& candidate) {
        
        if (is_inS(Y[u]) || (candidate == Y[u])) return 0;

        uint64_t mask = Smask | candidate;
        if ((mask & Y[u]) != Y[u]) return 999;
        
        if (Dist[u] != 999) {
            if (reachable(Y[u] ^ candidate, Dist[u]-1, 0, (1 << HY[u]) - (1 << D[S_map[candidate]]))) {
                return Dist[u] - 1;
            } else {
                return Dist[u];
            }
        } else {
            int test_size = __builtin_popcountll(Y[u]) - 1;
            test_size = 8;
            for (int dist = 1; dist <= test_size; dist++) {
                if (reachable(Y[u] ^ candidate, dist - 1, 0, (1 << HY[u]) - (1 << D[S_map[candidate]]))) {
                    return dist - 1;
                }
            }
            return 999;
        }
    }

    int update_distance_for_new_pair(uint64_t w, int w_depth, uint64_t y, int Hy, int current_dist) {

        if (is_inS(y) || (w == y)) return 0;
        
        uint64_t mask = Smask | w;
        if ((mask & y) != y) return 999;
        
        if (current_dist != 999) {
            if (reachable(y ^ w, current_dist-1, 0, (1 << Hy) - (1 << w_depth))) {
                return current_dist - 1;
            } else {
                return current_dist;
            }
        } else {
            int test_size = __builtin_popcountll(y) - 1;
            test_size = 8;
            for (int dist = 1; dist <= test_size; dist++) {
                if (reachable(y ^ w, dist - 1, 0, (1 << Hy) - (1 << w_depth))) {
                    return dist - 1;
                }
            }
            return 999;

        }
    }

    bool reachable(uint64_t T, int K, int index, int d2) const {
        if ((S.size() - index) < K) return false;
       
        if (d2 < 0) return false;
        
        if (K <= 4) {
            auto it = S_XORs_2PD[K].find(T);
            if (it != S_XORs_2PD[K].end()) {
                int cost = it->second; 
                if (d2 - cost < 0) return false;
                else return true;
            } else {
                return false;
            }
        }
              
        return reachable(T ^ S[index], K - 1, index + 1, d2 - (1 << D[index])) || reachable(T, K, index + 1, d2);

    }

    string getAfterSlash(const string& str) {
        size_t pos = str.find_last_of('/');
        string pp = (pos == string::npos) ? "" : str.substr(pos + 1);
        size_t pos1 = pp.find_last_of('.');

        return (pos1 == string::npos) ? "" : pp.substr(0, pos1);
    }

    bool write_result() {
        ifstream infile("code_target_imps/" + Filename + ".py");
        if (!infile.is_open()) {
            return false;
        }
        
        stringstream buffer;
        buffer << infile.rdbuf();
        string content = buffer.str();
        infile.close();
        
        string delimiter = "################### Here is your code !! ###################";
        size_t pos1 = content.find(delimiter);
        size_t pos2 = content.find(delimiter, pos1 + delimiter.length());
        
        if (pos1 == string::npos || pos2 == string::npos) {
            cerr << "not found -- Here is your code !! " << endl;
            return false;
        }
        
        string head = content.substr(0, pos1);
        string body = content.substr(pos1 + delimiter.length(), 
                                    pos2 - pos1 - delimiter.length());
        string tail = content.substr(pos2 + delimiter.length());

        string out_filename = Filename + "_H" + to_string(depth) + "_XORS" + to_string(S.size() - n - k) + "_" +  getAfterSlash(LOG_filename);
        ofstream outfile("Results/" + out_filename + ".py");
        if (!outfile.is_open()) {
            cerr << "output file error!" << endl;
            return false;
        }
        
        outfile << head;
        outfile << "        ################### Here is your code !! ###################\n";

        string space_str = "        ";
        outfile << space_str << "t=[0]*" << to_string(ops.size() + n) << endl;
        for(int i = 0; i < n; i++){
            outfile << space_str << "t[" << i << "]=x[" << i << "]" << endl;
        }
        for (const op_t& gate : ops) {
            if(gate.type == "g") outfile << space_str << "t[" << gate.dst << "]=t[" << gate.src.first << "]" << NLs[gate.type_index] <<"t["<< gate.src.second << "]" << endl;
            else if(gate.type == "y"){
                outfile << space_str << "t[" << gate.dst << "]=t[" << gate.src.first << "]^t["<< gate.src.second << "]" << endl;
                
                string not_add = "y[" + to_string(gate.type_index) + "]";
                if(find(NOTs.begin(), NOTs.end(), not_add) != NOTs.end()) outfile <<  space_str << "y[" << gate.type_index << "]=t[" << gate.dst << "]^1" << endl;
                else outfile <<  space_str << "y[" << gate.type_index << "]=t[" << gate.dst << "]" << endl;
            }
            else outfile << space_str << "t[" << gate.dst << "]=t[" << gate.src.first << "]^t["<< gate.src.second << "]" << endl;
        }

        outfile << "        ################### Here is your code !! ###################\n";
        outfile << tail;
        outfile.close();
        
        PRINTF_STAMP("Generate vartify code: Results/%s\n\n", out_filename.c_str() );
        return true;
    }

    
    void ana_candi(const map<uint64_t, tuple<int, int, int>>& candidates_with_sum) {
        if (candidates_with_sum.empty()) {
            log_file << "Map is empty!" << endl;
            return;
        }

        int min_first = numeric_limits<int>::max();
        int max_first = numeric_limits<int>::min();
        int min_second = numeric_limits<int>::max();
        int max_second = numeric_limits<int>::min();
        int min_third = numeric_limits<int>::max();
        int max_third = numeric_limits<int>::lowest();
    
        for (const auto& entry : candidates_with_sum) {
            const auto& tuple_val = entry.second;
            min_first = min(min_first, get<0>(tuple_val));
            max_first = max(max_first, get<0>(tuple_val));
            min_second = min(min_second, get<1>(tuple_val));
            max_second = max(max_second, get<1>(tuple_val));
            min_third = min(min_third, get<2>(tuple_val));
            max_third = max(max_third, get<2>(tuple_val));
        }
    
        log_file << "\n============" << endl;
        log_file << "Count999: " << min_first << " --- " << max_first << endl;
        log_file << "Sum: " << min_second << " --- " << max_second << endl;
        log_file << "Norm: " << min_third << " --- " << max_third << endl;
        log_file << endl;

        vector<pair<uint64_t, tuple<int, int, int>>> sorted_candidates(
            candidates_with_sum.begin(), candidates_with_sum.end()
        );
    
        sort(sorted_candidates.begin(), sorted_candidates.end(),
        [](const auto& a, const auto& b) {
            const auto& tuple_a = a.second;
            const auto& tuple_b = b.second;

            if (get<0>(tuple_a) != get<0>(tuple_b)) {
                return get<0>(tuple_a) < get<0>(tuple_b);
            }

            if (get<1>(tuple_a) != get<1>(tuple_b)) {
                return get<1>(tuple_a) < get<1>(tuple_b);
            }
  
            return get<2>(tuple_a) > get<2>(tuple_b);
        });
    


        log_file << " ---- " << endl;
        int count = 0;
        for (const auto& candidate : sorted_candidates) {
            if (count >= 50) break;

            if (count % 3 == 0) {
                if (count > 0) log_file << endl; 
            } else {
                log_file << "   "; 
            }
    
            log_file << candidate.first << ": " 
                      << get<0>(candidate.second) << ","
                      << get<1>(candidate.second) << ","
                      << fixed << setprecision(2) << get<2>(candidate.second);
    
            count++;
        }
        log_file << "\n============" << endl;
    }
    
    void print_eccandi(const vector<pair<set<uint64_t>, int>>& ec_candi) {
        if (ec_candi.empty()) {
            log_file << "Map is empty!" << endl;
            return;
        }

        int count = 0;
        for (const auto& entry : ec_candi) {

            if (count % 3 == 0) {
                if (count > 0) log_file << endl;
            } else {
                log_file << "   "; 
            }
    
            stringstream ss;
            bool first = true;
            for (uint64_t num : entry.first) {
                if (!first) {
                    ss << ",";
                }
                ss << num;
                first = false;
            }

            log_file << ss.str() << ": " << entry.second;
            
            count++;
        }
        log_file << endl;
    }

};
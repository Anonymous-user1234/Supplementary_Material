#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdint>
#include <bitset>
#include <fstream>
#include <sstream>
#include <chrono> 
#include <algorithm>
#include <utility>
#include <assert.h>
using namespace std;
using namespace std::chrono;

#ifndef OUTFLAG
#define OUTFLAG true
// #define OUTFLAG false
#endif

#ifndef INDIM
#define INDIM (64)
#endif
#define COL_SIZE (INDIM)

#ifndef OUTDIM
#define OUTDIM (INDIM)
#endif
#define ROW_SIZE (OUTDIM)

typedef uint64_t bscol;
typedef uint64_t bsrow;

typedef std::bitset<COL_SIZE> bscol_display;
typedef std::pair<int, int> pair_t;

struct op_t {
    int dst = -1; 
    pair_t src = {-1,-1};
    // string sname = "";
    string type = "";
    int type_index = -1;
};

pair<int, int> getBitTwoPos(const bscol& bs) {
    int pos1 = -1, pos2 = -1;
    
    bscol temp = bs;
    if (temp != 0) {
        pos1 = __builtin_ctzll(temp); 
        temp &= temp - 1; // eliminate rightest 1
        if (temp != 0) {
            pos2 = __builtin_ctzll(temp);
        }
    }
    
    return {pos1, pos2};
}

class matrix {
private:
    vector<bscol> rows_;
    vector<bscol> basis_vec_;           
    unordered_map<bscol, int> basis_map_;
    vector<int> dist_; 
    vector<op_t> ops_;

    vector<int> temp_dist_;

public:
    matrix(const string& filename) : rows_(ROW_SIZE), dist_(ROW_SIZE), temp_dist_(ROW_SIZE) {
        // init rows
        ifstream file(filename);
        if (!file.is_open()) {
            exit(1);
        }
        string line;
        int rowIndex = 0;
        while (getline(file, line) && rowIndex < ROW_SIZE) {
            stringstream ss(line);
            string token;
            bscol rowValue = 0;
            int bitIndex = 0;
    
            while (getline(ss, token, ',')) {
                if (token == "1") {
                    // rowValue |= (1ULL << (COL_SIZE - 1 - bitIndex)); 
                    rowValue |= (1ULL << bitIndex); 
                }
                bitIndex++;
            }
    
            rows_[rowIndex] = rowValue;
            ++rowIndex;
        }
    
        file.close();

        if(OUTFLAG){    
            cout << "Loaded matrix: " << filename << "\nROW: " << ROW_SIZE << ", COL:" << COL_SIZE << endl;
            for (int i = 0; i < ROW_SIZE; ++i) {
                cout << bscol_display(rows_[i]) << endl;
            }
        }

        // init basis_vec_ and basis_map_
        basis_vec_.clear();
        basis_map_.clear();
        for (int i = 0; i < COL_SIZE; ++i) {
            bscol base_elem = (1ULL << i);
            basis_vec_.push_back(base_elem);
            basis_map_[base_elem] = i;
        }
    
        for (int i = 0; i < ROW_SIZE; ++i) {
            dist_[i] = __builtin_popcountll(rows_[i]) - 1; 
        }

        if(OUTFLAG){
            cout << "Initial basis set:" << endl;
            for (auto b : basis_vec_) {
                cout << bscol_display(b) << " " << basis_map_[b] << endl;
            }
            cout << "Initial dist array: ";
            for (int i = 0; i < ROW_SIZE; ++i) {
                cout << dist_[i] << " ";
            }
            cout << endl;
        }    
    }

    bool isInBasis(const bscol& elem) const {
        return basis_map_.find(elem) != basis_map_.end();
    }

    void addToBasis(const bscol& elem, const pair_t& sources) {
        int new_index = basis_vec_.size();
        basis_vec_.push_back(elem);
        basis_map_[elem] = new_index;
        ops_.emplace_back(op_t{new_index, sources});
    }

    bool reachable(bscol T, int K, int S) const {
        if ((basis_vec_.size() - S) < K) return false;
        if (K == 0) return false;
        if (K == 1) {
            for (int i = S; i < basis_vec_.size(); ++i) {
                if (T == basis_vec_[i]) return true;
            }
            return false;
        }
        return reachable(T ^ basis_vec_[S], K - 1, S + 1) || reachable(T, K, S + 1);
    }

    int NewDistance(int u, const bscol& candidate) const {
        if (isInBasis(rows_[u]) || (candidate == rows_[u])) return 0;
        if (reachable(rows_[u] ^ candidate, dist_[u]-1, 0)) return dist_[u]-1;
        else return dist_[u];
    } 

    void updateDist(const bscol& candidate) {
        for (int i = 0; i < ROW_SIZE; ++i) {
            dist_[i] = NewDistance(i, candidate);
        }
        if(OUTFLAG){
            cout << "dist: " << bscol_display(candidate) << endl;
            cout << "[ ";
            for (int i = 0; i < ROW_SIZE; ++i) {
                cout << dist_[i] << " ";
            }
            cout << " ]" << endl;
        }
    }

    pair<int, int> evaluateCandidate(const bscol& candidate) {
        int sum = 0;
        int norm = 0;
        
        for (int i = 0; i < ROW_SIZE; ++i) {
            temp_dist_[i] = NewDistance(i, candidate);
            sum += temp_dist_[i];
            norm += temp_dist_[i] * temp_dist_[i];
        }
        // norm = sqrt(norm);
        return {sum, norm};
    }

    void find_imp() {
        bscol best_choose = 0;
        
        while (true) {
            // finish check
            bool done = true;
            for (int d : dist_) {
                if (d != 0) {
                    done = false;
                    break;
                }
            }
            if (done) break;
    
            // step1: easy case
            bool found_easy = false;
            for (int i = 0; i < ROW_SIZE; i++) {
                if (dist_[i] == 1) {
                    if(best_choose != 0){
                        bscol src0 = best_choose;
                        bscol src1 = best_choose ^ rows_[i];
                        
                        assert(isInBasis(src0) && isInBasis(src1));
                        addToBasis(rows_[i], {basis_map_[src0], basis_map_[src1]});
                    } else {
                        assert(__builtin_popcountll(rows_[i]) == 2);
                        auto [src0_i, src1_i] = getBitTwoPos(rows_[i]);
                        addToBasis(rows_[i], {src0_i, src1_i});
                    }
                    updateDist(rows_[i]);
                    found_easy = true;
                    break; 
                }
            }
            if (found_easy) continue;
    
            // step2: generate candidate
            int min_sum = INT_MAX;
            int max_norm = -1;
            bscol best_candidate = 0;
            pair_t best_sources;
            
            int current_size = basis_vec_.size();
            
            for (int i = 0; i < current_size; ++i) {
                for (int j = i + 1; j < current_size; ++j) {
                    bscol candidate = basis_vec_[i] ^ basis_vec_[j];
                    
                    if (isInBasis(candidate)) continue;
                    
                    auto [current_sum, current_norm] = evaluateCandidate(candidate);
                    
                    if (current_sum < min_sum || 
                       (current_sum == min_sum && current_norm > max_norm)) {
                        min_sum = current_sum;
                        max_norm = current_norm;
                        best_candidate = candidate;
                        best_sources = {basis_map_[basis_vec_[i]], basis_map_[basis_vec_[j]]};
                    }
                }
            }
            
            if (best_candidate == 0) {
                cout << "no candidate, break" << endl;
                break;
            }
            if(OUTFLAG) {
                cout << "\n\nSelected best candidate: " << bscol_display(best_candidate) 
                     << " (sum: " << min_sum << ", norm: " << max_norm << ")" << endl;
            }
    
            //step3: add best chooes to basis
            addToBasis(best_candidate, best_sources);
            updateDist(best_candidate);
            best_choose = best_candidate;
        }
    
        cout << "Optimal circuit size: " << basis_vec_.size() - COL_SIZE << endl;
        //output circuit
        if(OUTFLAG) {
            cout << "\nConstruction sequence:" << endl;
            for (const auto& op : ops_) {
                cout << "x" << op.dst << " = x" << op.src.first << " ^ x" << op.src.second << endl;
            }
        }
    }
};


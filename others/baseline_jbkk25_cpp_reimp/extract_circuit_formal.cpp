#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdint>
#include <bitset>
#include <chrono> 
#include <algorithm>
#include <utility>
#include <assert.h>
#include <regex>
#include "time.inl"
#include <random>

using namespace std;

int calculate_num(unordered_map<string, vector<string>>& XORs) {
    int count = 0;
    for (const auto& pair : XORs) {
        count += pair.second.size() - 1;
    }
    return count;
}


int caculate_depth(int n, int m, int k, const unordered_map<string, vector<string>>& XORs, const vector<string>& NLs, const vector<string>& NOTs) {
    unordered_map<string, int> D;
    
    for (int i = 0; i < n; ++i) {
        D["x[" + to_string(i) + "]"] = 0;
    }
    
    for (int i = 0; i < k; ++i) {
        string r0 = "r[" + to_string(2 * i) + "]";
        string r1 = "r[" + to_string(2 * i + 1) + "]";
        
        if (XORs.find(r0) != XORs.end()) {
            uint64_t sum = 0;
            for (const auto& input : XORs.at(r0)) {
                sum += (1ULL << D[input]);
            }
            D[r0] = static_cast<int>(ceil(log2(sum)));
        }
        
        if (XORs.find(r1) != XORs.end()) {
            uint64_t sum = 0;
            for (const auto& input : XORs.at(r1)) {
                sum += (1ULL << D[input]);
            }
            D[r1] = static_cast<int>(ceil(log2(sum)));
        }
        
        D["g[" + to_string(i) + "]"] = max(D[r0], D[r1]) + 1;
    }
    
    int max_depth = 0;
    for (int i = 0; i < m; ++i) {
        string y_var = "y[" + to_string(i) + "]";
        if (XORs.find(y_var) != XORs.end()) {
            uint64_t sum = 0;
            for (const auto& input : XORs.at(y_var)) {
                sum += (1ULL << D[input]);
            }
            D[y_var] = static_cast<int>(ceil(log2(sum)));
            max_depth = max(max_depth, D[y_var]);
        }
    }
    
    return max_depth;
}

void AND_AND_ver2(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    string g_var = "g[" + to_string(i) + "]";
    
    // handle XORs: r1 -> r1^r0
    vector<string> r0_inputs = XORs[r0];
    for (const auto& x : r0_inputs) {
        auto it = find(XORs[r1].begin(), XORs[r1].end(), x);
        if (it != XORs[r1].end()) {
            XORs[r1].erase(it);
        } else {
            XORs[r1].push_back(x);
        }
    }
    
    // handle NOTs: r1 -> r1^r0
    auto r0_not_it = find(NOTs.begin(), NOTs.end(), r0);
    if (r0_not_it != NOTs.end()) {
        auto r1_not_it = find(NOTs.begin(), NOTs.end(), r1);
        if (r1_not_it != NOTs.end()) {
            NOTs.erase(r1_not_it);
        } else {
            NOTs.push_back(r1);
        }
    }
    
    // handle g: g -> g^r0
    for (auto& pair : XORs) {
        if (find(pair.second.begin(), pair.second.end(), g_var) != pair.second.end()) {
            // XORs: g -> g^r0
            for (const auto& x : XORs[r0]) {
                auto it = find(pair.second.begin(), pair.second.end(), x);
                if (it != pair.second.end()) {
                    pair.second.erase(it);
                } else {
                    pair.second.push_back(x);
                }
            }
            
            // NOTs: g -> g^r0
            if (find(NOTs.begin(), NOTs.end(), r0) != NOTs.end()) {
                auto not_it = find(NOTs.begin(), NOTs.end(), pair.first);
                if (not_it != NOTs.end()) {
                    NOTs.erase(not_it);
                } else {
                    NOTs.push_back(pair.first);
                }
            }
        }
    }
}

void AND_AND_ver3(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    string g_var = "g[" + to_string(i) + "]";
    
    // handle XORs: r0 -> r0^r1
    vector<string> r1_inputs = XORs[r1];
    for (const auto& x : r1_inputs) {
        auto it = find(XORs[r0].begin(), XORs[r0].end(), x);
        if (it != XORs[r0].end()) {
            XORs[r0].erase(it);
        } else {
            XORs[r0].push_back(x);
        }
    }
    
    // handle NOTs: r0 -> r0^r1
    auto r1_not_it = find(NOTs.begin(), NOTs.end(), r1);
    if (r1_not_it != NOTs.end()) {
        auto r0_not_it = find(NOTs.begin(), NOTs.end(), r0);
        if (r0_not_it != NOTs.end()) {
            NOTs.erase(r0_not_it);
        } else {
            NOTs.push_back(r0);
        }
    }
    
    // handle g: g -> g^r1
    for (auto& pair : XORs) {
        if (find(pair.second.begin(), pair.second.end(), g_var) != pair.second.end()) {
            // XORs: g -> g^r1
            for (const auto& x : XORs[r1]) {
                auto it = find(pair.second.begin(), pair.second.end(), x);
                if (it != pair.second.end()) {
                    pair.second.erase(it);
                } else {
                    pair.second.push_back(x);
                }
            }
            
            // NOTs: g -> g^r0
            if (find(NOTs.begin(), NOTs.end(), r1) != NOTs.end()) {
                auto not_it = find(NOTs.begin(), NOTs.end(), pair.first);
                if (not_it != NOTs.end()) {
                    NOTs.erase(not_it);
                } else {
                    NOTs.push_back(pair.first);
                }
            }
        }
    }
}


void OR_OR_ver2(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    
    // r1 -> r1^r0
    vector<string> r0_inputs = XORs[r0];
    for (const auto& x : r0_inputs) {
        auto it = find(XORs[r1].begin(), XORs[r1].end(), x);
        if (it != XORs[r1].end()) {
            XORs[r1].erase(it);
        } else {
            XORs[r1].push_back(x);
        }
    }

    // handle NOTs: r1 -> r1^r0
    auto r0_not_it = find(NOTs.begin(), NOTs.end(), r0);
    if (r0_not_it != NOTs.end()) {
        auto r1_not_it = find(NOTs.begin(), NOTs.end(), r1);
        if (r1_not_it != NOTs.end()) {
            NOTs.erase(r1_not_it);
        } else {
            NOTs.push_back(r1);
        }
    }
}

void OR_OR_ver3(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    
    // handle XORs: r0 -> r0^r1
    vector<string> r1_inputs = XORs[r1];
    for (const auto& x : r1_inputs) {
        auto it = find(XORs[r0].begin(), XORs[r0].end(), x);
        if (it != XORs[r0].end()) {
            XORs[r0].erase(it);
        } else {
            XORs[r0].push_back(x);
        }
    }
    // handle NOTs: r0 -> r0^r1
    auto r1_not_it = find(NOTs.begin(), NOTs.end(), r1);
    if (r1_not_it != NOTs.end()) {
        auto r0_not_it = find(NOTs.begin(), NOTs.end(), r0);
        if (r0_not_it != NOTs.end()) {
            NOTs.erase(r0_not_it);
        } else {
            NOTs.push_back(r0);
        }
    }
}


void AND_to_OR(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    string g_var = "g[" + to_string(i) + "]";
    
    for (auto& pair : XORs) {
        auto& outputs = pair.second;
        auto it = find(outputs.begin(), outputs.end(), g_var);
        if (it != outputs.end()) {
            // outputs.erase(it);
            
            // add r[2i]
            for (const auto& input : XORs[r0]) {
                auto input_it = find(outputs.begin(), outputs.end(), input);
                if (input_it != outputs.end()) {
                    outputs.erase(input_it);
                } else {
                    outputs.push_back(input);
                }
            }
            
            // add r[2i+1]
            for (const auto& input : XORs[r1]) {
                auto input_it = find(outputs.begin(), outputs.end(), input);
                if (input_it != outputs.end()) {
                    outputs.erase(input_it);
                } else {
                    outputs.push_back(input);
                }
            }
            

            bool r0_not = (find(NOTs.begin(), NOTs.end(), r0) != NOTs.end());
            bool r1_not = (find(NOTs.begin(), NOTs.end(), r1) != NOTs.end());
            
            if (r0_not != r1_not) {
                auto not_it = find(NOTs.begin(), NOTs.end(), pair.first);
                if (not_it != NOTs.end()) {
                    NOTs.erase(not_it);
                } else {
                    NOTs.push_back(pair.first);
                }
            }
        }
    }
    
    // change gate type OR
    NLs[i] = '|';
}


void OR_to_AND(unordered_map<string, vector<string>>& XORs, vector<string>& NLs, vector<string>& NOTs, int i) {
    string r0 = "r[" + to_string(2 * i) + "]";
    string r1 = "r[" + to_string(2 * i + 1) + "]";
    string g_var = "g[" + to_string(i) + "]";
    
    for (auto& pair : XORs) {
        auto& outputs = pair.second;
        auto it = find(outputs.begin(), outputs.end(), g_var);
        if (it != outputs.end()) {
            // outputs.erase(it);
            
            // add r[2i] 
            for (const auto& input : XORs[r0]) {
                auto input_it = find(outputs.begin(), outputs.end(), input);
                if (input_it != outputs.end()) {
                    outputs.erase(input_it);
                } else {
                    outputs.push_back(input);
                }
            }
            
            // add r[2i+1]
            for (const auto& input : XORs[r1]) {
                auto input_it = find(outputs.begin(), outputs.end(), input);
                if (input_it != outputs.end()) {
                    outputs.erase(input_it);
                } else {
                    outputs.push_back(input);
                }
            }
            
            bool r0_not = (find(NOTs.begin(), NOTs.end(), r0) != NOTs.end());
            bool r1_not = (find(NOTs.begin(), NOTs.end(), r1) != NOTs.end());

            if (r0_not != r1_not) {
                auto not_it = find(NOTs.begin(), NOTs.end(), pair.first);
                if (not_it != NOTs.end()) {
                    NOTs.erase(not_it);
                } else {
                    NOTs.push_back(pair.first);
                }
            }
        }
    }
    
    // change gate type AND
    NLs[i] = '&';
}

class CirParser {
private:
    int n;  // input var num
    int m;  // output var num
    int k; // #NLgate

    //circuit information
    unordered_map<string, vector<string>> XORs;
    vector<string> NLs;
    vector<string> NOTs;
    
    //circuit bit 
    vector<uint64_t> S; 
    vector<uint64_t> D; 
    vector<uint64_t> Y; 
    vector<string> Sname;

    string remove_char(const string& line) {
        string result;
        for (char c : line) {
            if (c != ' ' && c != ';') {
                result += c;
            }
        }
        return result;
    }
    vector<string> split_assign(const string& line) {
        vector<string> parts;
        size_t pos = line.find('=');
        if (pos != string::npos) {
            parts.push_back(line.substr(0, pos));
            parts.push_back(line.substr(pos + 1));
        }
        return parts;
    }
    bool is_outvar(const string& var) {
        return var.find("y[") == 0 && var.back() == ']';
    }
    
    bool is_invar(const string& var) {
        return var.find("x[") == 0 && var.back() == ']';
    }
    
    
public:
    CirParser(const string& filename) {
        ifstream file("code_target_imps/" + filename + ".py");
        if (!file.is_open()) {
            return;
        }
        
        // read file
        stringstream buffer;
        buffer << file.rdbuf();
        string content = buffer.str();
        file.close();
        
        string delimiter = "################### Here is your code !! ###################";
        size_t pos = content.find(delimiter);
        if (pos == string::npos) {
            cerr << "read file error!" << endl;
            return;
        }
        
        string body = content.substr(pos + delimiter.length());
        pos = body.find(delimiter);
        if (pos != string::npos) {
            body = body.substr(0, pos);
        }
        
        n = 0;
        m = 0;
        
        for (int i = 100; i >= 0; --i) {
            string x_pattern = "x[" + to_string(i) + "]";
            if (body.find(x_pattern) != string::npos) {
                n = i + 1;
                break;
            }
        }
        
        for (int i = 100; i >= 0; --i) {
            string y_pattern = "y[" + to_string(i) + "]";
            if (body.find(y_pattern) != string::npos) {
                m = i + 1;
                break;
            }
        }
        
        PRINTF_STAMP("n = %d, m = %d\n\n", n, m);
    }
    
    //step1： change circuit to standard formal 
    bool cir_formal(const string& filename) {
        ifstream infile("code_target_imps/" + filename + ".py");
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
        
        // handle (^=, &=, |=)
        vector<string> lines;
        istringstream bodyStream(body);
        string line;
        
        while (getline(bodyStream, line)) {
            string cleanLine;
            for (char c : line) {
                if (c != ' ' && c != ';') {
                    cleanLine += c;
                }
            }
            if (!cleanLine.empty()) {
                lines.push_back(cleanLine);
            }
        }
        
        vector<string> new_lines;
        int cnt = 0;
        int len_lines = lines.size();
        
        for (int i = 0; i < len_lines; i++) {
            if (lines[i].find('#') != string::npos) continue;
            if (lines[i].find('=') == string::npos) continue;
            
            if (lines[i].find("^=") != string::npos || lines[i].find("&=") != string::npos || lines[i].find("|=") != string::npos) {                
                size_t equal_pos = lines[i].find('=');
                string Y_op = lines[i].substr(0, equal_pos);
                string X = lines[i].substr(equal_pos + 1);
                
                char op = Y_op.back();
                string Y = Y_op.substr(0, Y_op.length() - 1);
                
                for (int j = i + 1; j < len_lines; j++) {
                    size_t pos = 0;
                    while ((pos = lines[j].find(Y, pos)) != string::npos) {
                        if ((pos == 0 || !isalnum(lines[j][pos-1])) && 
                            (pos + Y.length() >= lines[j].length() || !isalnum(lines[j][pos + Y.length()]))) {
                            lines[j].replace(pos, Y.length(), "DJ_master[" + to_string(cnt + 1) + "]");
                            pos += 12 + to_string(cnt + 1).length();
                        } else {
                            pos += Y.length();
                        }
                    }
                }
                
                new_lines.push_back("DJ_master[" + to_string(cnt) + "]=" + X);
                new_lines.push_back("DJ_master[" + to_string(cnt + 1) + "]=" + Y + op + "DJ_master[" + to_string(cnt) + "]");
                cnt += 2;
            } else {
                new_lines.push_back(lines[i]);
            }
        }
        
        vector<string> new2_lines;
        len_lines = new_lines.size();
        for (int i = 0; i < len_lines; i++) {
            if (new_lines[i].find('#') != string::npos) continue;
            if (new_lines[i].find('=') == string::npos) continue;
            
            size_t equal_pos = new_lines[i].find('=');
            string Y = new_lines[i].substr(0, equal_pos);
            string X = new_lines[i].substr(equal_pos + 1);
            
            bool is_output = false;
            for (int idx = 0; idx < m; idx++) {
                if (Y == "y[" + to_string(idx) + "]") {
                    is_output = true;
                    break;
                }
            }
            
            if (!is_output) {
                for (int j = i + 1; j < len_lines; j++) {
                    size_t pos = 0;
                    while ((pos = new_lines[j].find(Y, pos)) != string::npos) {
                        if ((pos == 0 || !isalnum(new_lines[j][pos-1])) && 
                            (pos + Y.length() >= new_lines[j].length() || !isalnum(new_lines[j][pos + Y.length()]))) {
                            new_lines[j].replace(pos, Y.length(), "DJ_master[" + to_string(cnt) + "]");
                            pos += 12 + to_string(cnt).length();
                        } else {
                            pos += Y.length();
                        }
                    }
                }
                new2_lines.push_back("DJ_master[" + to_string(cnt) + "]=" + X);
                cnt++;
            } else {
                new2_lines.push_back(new_lines[i]);
            }
        }
        
        vector<string> new3_lines;
        len_lines = new2_lines.size();
        for (int i = 0; i < len_lines; i++) {
            if (new2_lines[i].find('#') != string::npos) continue;
            if (new2_lines[i].find('=') == string::npos) continue;
            
            size_t equal_pos = new2_lines[i].find('=');
            string Y = new2_lines[i].substr(0, equal_pos);
            string X = new2_lines[i].substr(equal_pos + 1);
            
            bool is_output = false;
            for (int idx = 0; idx < m; idx++) {
                if (Y == "y[" + to_string(idx) + "]") {
                    is_output = true;
                    break;
                }
            }
            
            if (is_output) {
                for (int j = i + 1; j < len_lines; j++) {
                    size_t pos = 0;
                    while ((pos = new2_lines[j].find(Y, pos)) != string::npos) {
                        if ((pos == 0 || !isalnum(new2_lines[j][pos-1])) && 
                            (pos + Y.length() >= new2_lines[j].length() || !isalnum(new2_lines[j][pos + Y.length()]))) {
                            new2_lines[j].replace(pos, Y.length(), "DJ_master[" + to_string(cnt) + "]");
                            pos += 12 + to_string(cnt).length();
                        } else {
                            pos += Y.length();
                        }
                    }
                }
                new3_lines.push_back("DJ_master[" + to_string(cnt) + "]=" + X);
                new3_lines.push_back(Y + "=DJ_master[" + to_string(cnt) + "]");
                cnt++;
            } else {
                new3_lines.push_back(new2_lines[i]);
            }
        }
        
        ofstream outfile("code_formal/" + filename + ".py");
        if (!outfile.is_open()) {
            cerr << "output file error!" << endl;
            return false;
        }
        
        outfile << head;
        outfile << "        ################### Here is your code !! ###################\n";
        outfile << "        DJ_master=[0]*" << cnt << "\n";
        
        for (const auto& line : new3_lines) {
            outfile << "        " << line << "\n";
        }
        
        outfile << "        ################### Here is your code !! ###################\n";
        outfile << tail;
        outfile.close();
        
        PRINTF_STAMP("Generate standard formal, #DJ_master: %d\n\n", cnt );
        return true;
    }
    
    //step2: init XORs, NLs, NOTs
    bool extract_info(const string& filename) {
        ifstream infile("code_formal/" + filename + ".py");
        if (!infile.is_open()) {
            return false;
        }
        stringstream buffer;
        buffer << infile.rdbuf();
        string content = buffer.str();
        infile.close();
        
        string delimiter = "################### Here is your code !! ###################";
        size_t pos = content.find(delimiter);
        if (pos == string::npos) {
            cerr << "not found -- Here is your code !!" << endl;
            return false;
        }
        
        string body = content.substr(pos + delimiter.length());
        pos = body.find(delimiter);
        if (pos != string::npos) {
            body = body.substr(0, pos);
        }

        vector<string> lines;
        istringstream bodyStream(body);
        string line;
        
        while (getline(bodyStream, line)) {
            string cleanLine;
            for (char c : line) {
                if (c != ' ' && c != ';' && !isspace(c)) {
                    cleanLine += c;
                }
            }
            
            if (!cleanLine.empty() && cleanLine.find("DJ_master=[0]*") == string::npos) {
                lines.push_back(cleanLine);
            }
        }
        
        unordered_map<string, vector<string>> pureXORs;
        k = 0; 
        
        for (const auto& line : lines) {
            if (line.find('=') == string::npos) continue;
            
            auto parts = split_assign(line);
            if (parts.size() != 2) continue;
            
            string Y = parts[0], X = parts[1];
            auto cleanStr = [](string& s) {
                s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
            };
            cleanStr(Y);
            cleanStr(X);

            if (is_outvar(Y)) {
                XORs[Y] = {X};
            } else if (X.find('&') != string::npos || X.find('|') != string::npos) {
                string op = (X.find("&") != string::npos) ? "&" : "|";
                NLs.push_back(op);
                
                size_t opPos = X.find(op);
                string r0 = X.substr(0, opPos);
                string r1 = X.substr(opPos + 1);
                cleanStr(r0);
                cleanStr(r1);
                string g_var = "g[" + to_string(k) + "]";
                // g replace y
                for (size_t j = 0; j < lines.size(); j++) {
                    if (lines[j] != line) { 
                        size_t pos = 0;
                        while ((pos = lines[j].find(Y, pos)) != string::npos) {
                            if ((pos == 0 || !isalnum(lines[j][pos-1])) && 
                                (pos + Y.length() >= lines[j].length() || !isalnum(lines[j][pos + Y.length()]))) {
                                lines[j].replace(pos, Y.length(), g_var);
                                pos += g_var.length();
                            } else {
                                pos += Y.length();
                            }
                        }
                    }
                }

                XORs["r[" + to_string(2 * k) + "]"] = {r0};
                XORs["r[" + to_string(2 * k + 1) + "]"] = {r1};
                k++;
            } else if (X.find('^') != string::npos) {
                vector<string> xorInputs;
                size_t start = 0, end = 0;
                while ((end = X.find('^', start)) != string::npos) {
                    string token = X.substr(start, end - start);
                    cleanStr(token);
                    if (!token.empty()) {
                        xorInputs.push_back(token);
                    }
                    start = end + 1;
                }
                string lastToken = X.substr(start);
                cleanStr(lastToken);
                if (!lastToken.empty()) {
                    xorInputs.push_back(lastToken);
                }
                pureXORs[Y] = xorInputs;
            } else {
                pureXORs[Y] = {X};
            }
        }


        // expand XOR
        for (auto& pair : XORs) {
            bool changed;
            do {
                changed = false;
                vector<string> newInputs;
                for (const auto& input : pair.second) {
                    if (is_invar(input) || input.find("g[") == 0 || input == "1") {
                        //directly add x,g,1 
                        auto it = find(newInputs.begin(), newInputs.end(), input);
                        if (it != newInputs.end()) {
                            newInputs.erase(it);
                        } else {
                            newInputs.push_back(input);
                        }
                    } else if (pureXORs.find(input) != pureXORs.end()) {
                        changed = true;
                        for (const auto& expanded : pureXORs[input]) {
                            // remove if it already exists, otherwise add
                            auto it = find(newInputs.begin(), newInputs.end(), expanded);
                            if (it != newInputs.end()) {
                                newInputs.erase(it);//1+1=0
                            } else {
                                newInputs.push_back(expanded);
                            }
                        }
                    } else {
                        newInputs.push_back(input); //not happen
                    }
                }
                pair.second = newInputs;
            } while (changed);// all vars in {x,g,1}
        }
        
        // init NOTs
        for (auto& pair : XORs) {
            auto it = find(pair.second.begin(), pair.second.end(), "1");
            if (it != pair.second.end()) {
                pair.second.erase(it); // delete 1
                NOTs.push_back(pair.first);
            }
        }
        
        PRINTF_STAMP("Init XORs,NLs,NOTs finish. k= %d\n\n", k );
        return true;
    }
    
    bool modify_cir() {
        auto origin_XORs = XORs;
        auto origin_NLs = NLs;
        auto origin_NOTs = NOTs;

        int origin_depth = caculate_depth(n, m, k, origin_XORs, origin_NLs, origin_NOTs);
        int origin_xornum = calculate_num(origin_XORs);

        PRINTF_STAMP("origin_depth: %d\n\n", origin_depth);
        PRINTF_STAMP("origin_XORnum: %d\n\n", origin_xornum);

        random_device rd;
        mt19937 gen(rd());

        int best_depth = 0;
        int best_xornum = 0;
        
        while (true) {

            auto new_XORs = origin_XORs;
            auto new_NLs = origin_NLs;
            auto new_NOTs = origin_NOTs;
            
            vector<string> gate_types = {"&", "|"};
            uniform_int_distribution<> type_dist(0, gate_types.size() - 1);
            
            vector<string> new_NL;
            for (size_t i = 0; i < new_NLs.size(); ++i) {
                new_NL.push_back(gate_types[type_dist(gen)]);
            }
            
            vector<string> versions = {"0", "1", "2", "3", "4", "5"};
            uniform_int_distribution<> ver_dist(0, versions.size() - 1);
            vector<string> new_ver;

            for (size_t i = 0; i < new_NLs.size(); ++i) {
                new_ver.push_back(versions[ver_dist(gen)]);

                string current_gate = new_NLs[i];
                string target_gate = new_NL[i];
                
                // if (current_gate == "&" && target_gate == "|") {
                //     AND_to_OR(new_XORs, new_NLs, new_NOTs, i);
                // } else if (current_gate == "|" && target_gate == "&") {
                //     OR_to_AND(new_XORs, new_NLs, new_NOTs, i);
                // }

                if (new_ver[i] == "2") {
                    if (new_NLs[i] == "&") {
                        AND_AND_ver2(new_XORs, new_NLs, new_NOTs, i);
                    } 
                    // else if (new_NLs[i] == "|") {
                    //     OR_OR_ver2(new_XORs, new_NLs, new_NOTs, i);
                    // }
                } else if (new_ver[i] == "3") {
                    if (new_NLs[i] == "&") {
                        AND_AND_ver3(new_XORs, new_NLs, new_NOTs, i);
                    } 
                    // else if (new_NLs[i] == "|" ) {
                    //     OR_OR_ver3(new_XORs, new_NLs, new_NOTs, i);
                    // }
                }

            }
            
            
            int depth = caculate_depth(n, m, k, new_XORs, new_NLs, new_NOTs);
            int xornum = calculate_num(new_XORs);
            
            if (depth <= origin_depth && xornum < static_cast<int>(origin_xornum * 0.9)) {
                XORs = new_XORs;
                NLs = new_NLs;
                NOTs = new_NOTs;
                best_depth = depth;
                best_xornum = xornum;
                break;
            }
        }

        PRINTF_STAMP("modified_depth: %d\n\n", best_depth);
        PRINTF_STAMP("modified_XORnum: %d\n\n", best_xornum);

        PRINTF_STAMP("Circuit modification completed. \n\n");
        // exit(1);
        return true;
    }



    //step3：init S, D, Y, Sname
    void init_bit() {
        //init S, Sname
        for (int i = 0; i < n; ++i) {
            S.push_back(1ULL << i);
            Sname.push_back("x[" + to_string(i) + "]");
        }
        
        D = vector<uint64_t>(n, 0);//init deep = 0 
        
        for (int i = 0; i < 2 * k; ++i) {
            string r_var = "r[" + to_string(i) + "]";
            uint64_t v = 0;
            
            if (XORs.find(r_var) != XORs.end()) {
                for (const auto& input : XORs[r_var]) {
                    if (is_invar(input)) {
                        int index = stoi(input.substr(2, input.length() - 3));
                        v |= 1ULL << index;
                    } else if (input.find("g[") == 0) {
                        int index = stoi(input.substr(2, input.length() - 3));
                        v |= 1ULL << (n + index);
                    }
                }
            }
            Y.push_back(v);
        }
        
        for (int i = 0; i < m; ++i) {
            string y_var = "y[" + to_string(i) + "]";
            uint64_t v = 0;
            
            if (XORs.find(y_var) != XORs.end()) {
                for (const auto& input : XORs[y_var]) {
                    if (is_invar(input)) {
                        int index = stoi(input.substr(2, input.length() - 3));
                        v |= 1ULL << index;
                    } else if (input.find("g[") == 0) {
                        int index = stoi(input.substr(2, input.length() - 3));
                        v |= 1ULL << (n + index);
                    }
                }
            }
            Y.push_back(v);
        }
        
        PRINTF_STAMP("init S, D, Y, Sname finish");
    }

    const unordered_map<string, vector<string>>& getXORs() const { return XORs; }
    const vector<string>& getNLs() const { return NLs; }
    const vector<string>& getNOTs() const { return NOTs; }
    const vector<uint64_t>& getS() const { return S; }
    const vector<uint64_t>& getD() const { return D; }
    const vector<uint64_t>& getY() const { return Y; }
    const vector<string>& getSname() const { return Sname; }
    int getn() const { return n; }
    int getm() const { return m; }
    int getk() const { return k; }

    void print_cir() const {
        cout << "\n\n=== Circuit information ===" << endl;
        cout << "Input var number: n = " << n << endl;
        cout << "Output var number: m = " << m << endl;
        cout << "Nonliner gate number: k = " << k << endl;
        
        auto extractIndex = [](const string& var) -> int {
            if (var.length() < 4 || var[1] != '[') return -1;
            return stoi(var.substr(2, var.length() - 3));
        };

        cout << "\n=== XORs ===" << endl;

        vector<pair<string, vector<string>>> sortedXORs(XORs.begin(), XORs.end());
        sort(sortedXORs.begin(), sortedXORs.end(), 
            [&](const auto& a, const auto& b) {
                // sort：r < y < others
                char type_a = a.first[0];
                char type_b = b.first[0];
                if (type_a != type_b) return type_a < type_b;
                
                // same type: index
                return extractIndex(a.first) < extractIndex(b.first);
            });
        
        for (const auto& pair : sortedXORs) {
            cout << pair.first << ": ";
            for (const auto& input : pair.second) {
                cout << input << " ";
            }
            cout << endl;
        }
        
        cout << "\n=== NLs ===" << endl;
        for (string op : NLs) {
            cout << op << " ";
        }
        cout << endl;
        
        cout << "\n=== NOTs ===" << endl;
        for (const auto& notVar : NOTs) {
            cout << notVar << " ";
        }
        cout << std::endl;
    }


    void print_cirbit() const {
        
        cout << "S:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "  " << Sname[i] << ": " << S[i];
            cout << " (" << i << ")" << endl;
        }
        
        cout << "Y:" << endl;
        for (int i = 0; i < 2 * k; i++) {
            cout << signalToString(Y[i]) << endl;
            // cout << "  r[" << i << "]: " << Y[i];
            // cout << " [";
            // bool first = true;
            // for (int j = 0; j < n; j++) {
            //     if (Y[i] & (1ULL << j)) {
            //         if (!first) cout << ", ";
            //         cout << "x[" << j << "]";
            //         first = false;
            //     }
            // }
            // for (int j = 0; j < k; j++) {
            //     if (Y[i] & (1ULL << (n + j))) {
            //         if (!first) cout << ", ";
            //         cout << "g[" << j << "]";
            //         first = false;
            //     }
            // }
            // cout << "]" << endl;
        }
        
        for (int i = 0; i < m; i++) {
            int index = 2 * k + i;
            cout << signalToString(Y[index]) << endl;
            // cout << "  y[" << i << "]: " << Y[index];
            // cout << " [";
            // bool first = true;
            // for (int j = 0; j < n; j++) {
            //     if (Y[index] & (1ULL << j)) {
            //         if (!first) cout << ", ";
            //         cout << "x[" << j << "]";
            //         first = false;
            //     }
            // }
            // for (int j = 0; j < k; j++) {
            //     if (Y[index] & (1ULL << (n + j))) {
            //         if (!first) cout << ", ";
            //         cout << "g[" << j << "]";
            //         first = false;
            //     }
            // }
            // cout << "]" << endl;
        }
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

    
};






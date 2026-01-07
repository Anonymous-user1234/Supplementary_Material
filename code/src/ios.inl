	


	///////////////////////////////////////// basic functions /////////////////////////////////////////
	inline std::string get_outs(Bvector const & v) { 
		std::string st;
		if (del_vec_infos.find(v) != del_vec_infos.end()) {
			if (del_vec_infos.at(v).key_len == 1)
				st = "1";
			else
				st = "h";
		} else {
			st = "*";
		}
		return st;
	}

	inline wind_t get_y_index(vec_t const & v) const { 
		auto const it = y_map.find(v);
		if (it != y_map.end()) return it->second;
		return -1;
	}


	inline wind_t get_target_index(vec_t const & v) const { 
		auto const it = target_map.find(v);
		if (it != target_map.end()) return it->second;
		return -1;
	}



	void table_print() const { 
		std::cout << "sbox table: " << std::endl;
		wind_t ctr = 0;
		for (auto const & f : table) {		
			if (ctr && (ctr % (word64_t(1) << (DIMENSION / 2)) == 0))
				std::cout << std::endl;
			printf(" 0x%02lX(%3lu)", f, f);
			++ctr;
		}
		std::cout << std::endl << std::endl;

	}


	void basis_print(std::vector<Bvector> const & _basis) const { 
		std::cout << "basis number: " << _basis.size() << std::endl;

		wind_t n_xor = _basis.size() - DIMENSION - NL_NUM;
		if (n_xor < 0) n_xor = 0;
		std::cout << "Xor number: " << n_xor << std::endl;

		wind_t ctr = 0;
		for(auto const & b : _basis) {
			std::cout << "[" << ctr++ << "]\t";
			b.print();

			auto ti = get_target_index(b.vector);
			if (ti != -1) {
				std::cout << "\t target " << ti;

				if (targets[ti].used_by_gs.size() > 0) std::cout << "\t by";
				for (auto const & gi : targets[ti].used_by_gs)
					std::cout << " g" << gi;
			} else {				
				auto idx = get_y_index(b.vector);
				if (idx != -1)
					std::cout << "\t y[" << idx << "]";				
			}

		
			std::cout << std::endl;
		}
		std::cout << "max ad depth: " << max_ad_depth << std::endl;		
		std::cout << "max xor depth: " << max_imp_depth << std::endl;
		std::cout << "Xor number: " << n_xor << std::endl;
		std::cout << "All number: " << _basis.size() - DIMENSION << std::endl << std::endl;
	}


	void basis_print() const {
		basis_print(basis);
	}



	inline void ops_print(std::vector<op_t> const & _ops) const {
		std::cout << "\tTotal Count = " << _ops.size() << std::endl;

		wind_t n_xor = 0;
		wind_t n_nxor = 0;
		wind_t n_and = 0;
		wind_t n_or = 0;

		for(auto const & op : _ops) {

			std::string st = "t[" + std::to_string(op.dst) + "] = ";
			if (op.src_0 < DIMENSION) st += "x";
			else st += "t";

			st += "[" + std::to_string(op.src_0) + "]";
			switch (op.type) {
				case TYPE_XOR:  { st += " ^ "; ++n_xor;  break; }
				case TYPE_NXOR: { st += " ^ "; ++n_nxor; break; }
				case TYPE_AND:  { st += " & "; ++n_and;  break; }
				case TYPE_OR:   { st += " | "; ++n_or;   break; }
				default:
					std::cout << "this should not happen in ops" << std::endl;
					assert(false);
			}

			if (op.src_1 < DIMENSION) st += "x";
			else st += "t";

			st += "[" + std::to_string(op.src_1) + "]";

			if (op.type == TYPE_NXOR) st += " ^ 1";
			else st += "    ";

			if (op.tar_idx >= 0) st += "\t ===> y[" + std::to_string(op.tar_idx) + "]";
			st += "\n";
			std::cout << st;
		}
		if (n_xor) std::cout << "Xor Count = " << n_xor << std::endl;
		if (n_nxor) std::cout << "NXOR Count = " << n_nxor << std::endl;
		if (n_and) std::cout << "AND Count = " << n_and << std::endl;
		if (n_or) std::cout << "OR Count = " << n_or << std::endl;
		std::cout << "total = " << _ops.size() << std::endl;
		std::cout << "max_ad_depth = " << max_ad_depth << std::endl;
		std::cout << "max_imp_depth = " << max_imp_depth << std::endl << std::endl;
	}

	void ops_print() const { 
		ops_print(ops);
	}



	bool ops_check(std::vector<op_t> const & _ops) const { 

		bool pass = true;

		for (wind_t ctr = 0; pass && ctr < table.size(); ++ctr) {
			std::vector<bool> t(DIMENSION + _ops.size());
			word64_t y_value = 0;

			for (wind_t i = 0; i < DIMENSION; ++i)
				if ((ctr >> i) & 1) t[i] = 1;
				else t[i] = 0;

			for (auto const & op : _ops) {

				if (op.dst > DIMENSION + _ops.size()) {
					std::cout << "!(op.dst > DIMENSION + _ops.size()), op.dst = " << op.dst << ", ops.size() = " << _ops.size() << std::endl;
					pass = false;
					break;
				}
				if (op.src_0 >= op.dst) {
					std::cout << "!(op.src_0 >= op.dst), op.src_0 = " << op.src_0 << ", op.dst = " << op.dst << std::endl;
					pass = false;
					break;
				}
				if (op.src_1 >= op.dst) {
					std::cout << "!(op.src_1 >= op.dst), op.src_1 = " << op.src_1 << ", op.dst = " << op.dst << std::endl;
					pass = false;
					break;
				}
				if (op.src_0 >= op.src_1) {
					std::cout << "!(op.src_0 >= op.src_1), op.src_1 = " << op.src_0 << ", op.src_1 = " << op.src_1 << std::endl;
					pass = false;
					break;
				}

				auto const & vs0 = t[op.src_0];
				auto const & vs1 = t[op.src_1];

				switch (op.type) {
					case TYPE_XOR:  { t[op.dst] = vs0 ^ vs1; break; }
					case TYPE_NXOR: { t[op.dst] = vs0 ^ vs1 ^ 1; break; }
					case TYPE_AND:  { t[op.dst] = vs0 & vs1; break; }
					case TYPE_OR:   { t[op.dst] = vs0 | vs1; break; }
					default: assert(false);
				}

				if (op.tar_idx >= 0 && t[op.dst])
					y_value ^= (word64_t(1) << op.tar_idx);

			}


			if (pass && y_value != table[ctr]) {
				std::cout << "y_value != table[x] error "
				          << "y_value = " << y_value
				          << ", x = " << ctr 
				          << ", table[x] = " << table[ctr] 
				          << std::endl;
				pass = false;
			}
		}


		if (pass) {
			if (SHOW_DETAILS_MATRIX) PRINTF_STAMP("operators CHECK passed ~~~~\n");
		} else {
			PRINTF_STAMP("CHECK not pass !!!!!!!!!!!!!!!!!!\n");
			exit(1);
		}
		return pass;
	}

	bool ops_check() const { 
		return ops_check(ops);
	}



	void ops_write(std::vector<op_t> _ops, bool use_hash){

		std::string str = std::string("results/") + EXPAND_STRINGIFY(FILENAME) + "_" + TIE_BREAKING_STRA + "_" + TAR_STRATEGY_STRA + "_d" + std::to_string(DEPTH_LIMIT);
		mkdir("results/", S_IRWXU);
    	mkdir(str.c_str(), S_IRWXU);

    	std::string out_file;

		out_file = str + "/ops_" + std::to_string(DEPTH_LIMIT) + "_" + std::to_string(_ops.size());

    	if (use_hash) {	

    		std::hash<vec_t> vechash;

    		size_t hash_val = 0;
    		for (auto const & b : _ops){
    			hash_val ^= vechash(b.src_0) ^ vechash(b.src_1) ^ vechash(b.type);	
			}
    		
    		out_file += "_hash" + std::to_string(hash_val);
    	}
    	out_file += ".txt";
		
    	std::cout << "\twrite file: " << out_file << std::endl << std::endl;


		wind_t n_xor = 0;
		wind_t n_nxor = 0;
		wind_t n_and = 0;
		wind_t n_or = 0;

		std::fstream fout;
     	fout.open(out_file.c_str(), std::ios::out);

		for(auto const & op : _ops) {

			std::string st = "t[" + std::to_string(op.dst) + "] = ";
			if (op.src_0 < DIMENSION) st += "x";
			else st += "t";

			st += "[" + std::to_string(op.src_0) + "]";
			switch (op.type) {
				case TYPE_XOR:  { st += " ^ "; ++n_xor;  break; }
				case TYPE_NXOR: { st += " ^ "; ++n_nxor; break; }
				case TYPE_AND:  { st += " & "; ++n_and;  break; }
				case TYPE_OR:   { st += " | "; ++n_or;   break; }
				default:
					fout << "this should not happen in ops" << std::endl;
					assert(false);
			}

			if (op.src_1 < DIMENSION) st += "x";
			else st += "t";

			st += "[" + std::to_string(op.src_1) + "]";

			if (op.type == TYPE_NXOR) st += " ^ 1";
			else st += "    ";

			if (op.tar_idx >= 0) st += "\t ===> y[" + std::to_string(op.tar_idx) + "]";
			st += "\n";
			fout << st;
		}

		if (n_xor) fout << "Xor Count = " << n_xor << std::endl;
		if (n_nxor) fout << "NXOR Count = " << n_nxor << std::endl;
		if (n_and) fout << "AND Count = " << n_and << std::endl;
		if (n_or) fout << "OR Count = " << n_or << std::endl;
		fout << "total = " << _ops.size() << std::endl;
		fout << "max_ad_depth = " << max_ad_depth << std::endl;
		fout << "max_imp_depth = " << max_imp_depth << std::endl << std::endl;
		fout.close();
	}


	// write ops to disk
	void ops_write() { 
		ops_write(ops, true);
	}



	// load ops and check
	bool ops_read(const std::string & in_file) { 

		std::ifstream inputFile(in_file);
		if (!inputFile.is_open()) {
		    std::cout << "error in open file!!! " << in_file << std::endl;
		    return false;
		}

		PRINTF_STAMP("loading solutions from file: %s\n", in_file.c_str());		    

		std::string line;
		int n_line = 0;
		ops.clear();
		while (std::getline(inputFile, line)) {
		  	if (line.size() == 0) continue;

		  	std::vector<wind_t> row;
		  	std::stringstream lineStream(line);
		  	std::string cell;

		  	while (std::getline(lineStream, cell, ' ')) 
		  		row.emplace_back(std::stoi(cell));
		  	assert(row.size() == 5);

		  	assert(row[1] != row[2]);
		  	if (row[1] < row[2])
		  		ops.emplace_back(op_t({row[0], row[1], row[2], row[3], row[4]}));
		  	else
		  		ops.emplace_back(op_t({row[0], row[2], row[1], row[3], row[4]}));

		  	++n_line;
		}
		inputFile.close();
		PRINTF_STAMP("loaded %d ops in from %s\n\n", n_line, in_file.c_str());

		if (SHOW_DETAILS_MATRIX) ops_print();

		ops_check();

		return true;
	}	



	///////////////////////////////////////// initialize functions /////////////////////////////////////////

	// b[i < DIMENSION] refers to ei.
	// b[DIMENSION] means 1.
	void set_init_basis() { 
		basis.clear();
		basis_map.clear();
		vector_mask.reset();

		for (wind_t c = 0; c < DIMENSION; ++c) {
			Bvector b(c, 0, 0);
			basis_map[b] = basis.size();
			basis.emplace_back(b);

			vector_mask[c] = 1;
		}
		max_imp_depth = 0;
		max_ad_depth = 0;

	}


	inline wind_t extend_target(vec_t const & v, wind_t const & gi) { 
		wind_t ti = get_target_index(v);
		
		assert(ti != -2);
		if (ti == -1) {
			ti = targets.size();
			Target T(v, ti, get_y_index(v));
			T.used_by_gs.emplace_back(gi);
			target_map[v] = ti;
			targets.emplace_back(T);
		} else {
			targets[ti].used_by_gs.emplace_back(gi);
		}
		assert(ti >= 0);
		return ti;
	}



	void set_targets() { 

		assert(ops.size() > 0);
		set_init_basis();


		max_imp_depth = 0;
		max_ad_depth = 0;

		targets.clear();
		std::vector<Target> ys;
		target_map.clear();
		y_map.clear();

		gys.clear();
		std::vector<GY> tmp_gy;


		for (auto const & op : ops) {
			assert(basis.size() > op.src_0);
			assert(basis.size() > op.src_1);
			assert(basis.size() == op.dst);

			Bvector b;
			assert(op.src_0 < op.src_1);
			auto const & vs0 = basis[op.src_0];
			auto const & vs1 = basis[op.src_1];


			if (op.type == TYPE_XOR || op.type == TYPE_NXOR) {
				b.linear_op(vs0, vs1);
				if (op.type == TYPE_NXOR)
					assert(op.tar_idx >= 0); // only in targets. the input ops can be normalized.

			} else if (op.type == TYPE_AND || op.type == TYPE_OR) {
				b.vector[DIMENSION + gys.size()] = 1;
				b.nonlinear_depth(vs0, vs1);

				GY g;
				g.ts.emplace_back(extend_target(vs0.vector, gys.size()));
				g.ts.emplace_back(extend_target(vs1.vector, gys.size()));
				g.t_vec.emplace_back(vs0.vector);
				g.t_vec.emplace_back(vs1.vector);
				g.type = op.type;
				g.ad_depth = b.ad_depth;



				gys.emplace_back(g);

			} else {
				assert(false);
			}


			// update basis
			basis.emplace_back(b);
			
			if (max_imp_depth < b.imp_depth)
				max_imp_depth = b.imp_depth;

			if (max_ad_depth < b.ad_depth)
				max_ad_depth = b.ad_depth;


			if (op.tar_idx >= 0) {
				Target T(b.vector, -2, op.tar_idx);
				T.has_one = (op.type == TYPE_NXOR);
				target_map[b.vector] = -2;
				y_map[b.vector] = op.tar_idx;
				
				ys.emplace_back(T);

				GY y;

				y.set_bound = DEPTH_LIMIT;

				y.type = (T.has_one) ? TYPE_NXOR : TYPE_XOR;
				y.ad_depth = b.ad_depth;

				tmp_gy.emplace_back(y);
			}
		}

		assert(gys.size() == NL_NUM);
		assert(ys.size() == tmp_gy.size());

		// deal with ys
		for (wind_t i = 0; i < ys.size(); ++i) {
			auto & T = ys[i];
			T.tar_index = targets.size();
			target_map[T.tar] = targets.size();

			tmp_gy[i].ts.emplace_back(targets.size());

			targets.emplace_back(T);
			gys.emplace_back(tmp_gy[i]);
		}
		assert(gys.size() == NL_NUM + DIMENSION);
		assert(target_map.size() <= targets.size());

		basis_print();

		std::cout << "total targets: " << targets.size()
				  << ", max_ad_depth: " << max_ad_depth
				  << ", max_imp_depth: " << max_imp_depth
				  << std::endl;
	}



	void gys_print() const { 
		assert(gys.size() == NL_NUM + DIMENSION);

		for (wind_t idx = 0; idx < NL_NUM; ++idx)
			gys[idx].print(idx);

		std::cout << "----------------" << std::endl;

		for (wind_t idx = 0; idx < DIMENSION; ++idx) {
			assert(gys[NL_NUM + idx].ts.size() == 1);
			auto const ti = gys[NL_NUM + idx].ts[0];
			assert(targets[ti].y_index >= 0);
			gys[NL_NUM + idx].print(NL_NUM + idx, targets[ti].y_index);
		}
		std::cout << "----------------" << std::endl;

	}


	inline void bp_inputs_check() { 
		assert(basis.size() >= DIMENSION);
		assert(basis.size() == basis_map.size());
		assert(targets.size() > 0);
		if (ops.size() > 0)
			assert(ops.size() + DIMENSION == basis.size());
		wind_t n_targets = 0;
		for (auto const & T : targets)
			n_targets += (T.dist > 0);
		assert(n_targets > 0);
	}



	// ///////////////////////////////////////// other basic functions /////////////////////////////////////////









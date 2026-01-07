
class LinComb { 

public:

    idx_t lin;
    word64_t value = 0;
    wind_t depth = -1;

    LinComb(idx_t const _lin) { 
    	lin = _lin;
    	value = lin.size();
    	depth = upper_log2(value);
    }


    LinComb(idx_t const _lin,
    		std::vector<Bvector> const & basis
    		) { 
    	lin = _lin;
    	value = 0;
    	for (auto const & idx : lin)
    		value += dvalue(basis[idx].imp_depth);
    	depth = upper_log2(value);
    }


    inline bool reduce_by_pair(idx_t & res, 
   							   pair_t const & pair, 
   							   wind_t const bidx
   							   ) const { 
    	res.clear();
    	auto it_0 = std::lower_bound(lin.begin(), lin.end(), pair.first);
    	if (it_0 == lin.end() || *it_0 != pair.first) return false;

    	auto it_1 = std::lower_bound(it_0 + 1, lin.end(), pair.second);
    	if (it_1 == lin.end() || *it_1 != pair.second) return false;

    	res.insert(res.end(), lin.begin(), it_0);
    	res.insert(res.end(), it_0 + 1, it_1);
    	res.insert(res.end(), it_1 + 1, lin.end());
    	res.emplace_back(bidx);
    	return true;
    }



    inline void print(bool enter = true) const { 
    	for (auto & idx : lin) std::cout << " " << idx;
    	std::cout << " [val " << value << ", d " << depth << "]";
    	if (enter) std::cout << std::endl;
	}
};







class Target {

public:

	vec_t tar; 
	vec_t masked_tar;

	wind_t tar_index = -1;
	wind_t y_index = -1; 

	std::vector<wind_t> used_by_gs; 



	wind_t dist = -1; 

	word64_t minimal_value = 0; 
	word64_t value_bound = 0;

	wind_t depth = -1; // only make sense when dist == 0
	wind_t depth_bound = -1; // could be dimension + nl_num, but use -1 to detect uninitialized cases.

	std::vector<LinComb> lincombs;
	std::set<pair_t> supports; 

	std::set<wind_t> unavai_gs; 


	bool has_one = false; 




	Target(vec_t const & _tar, wind_t const & _tar_index, wind_t const & _y_index) { 
		tar = _tar;
		masked_tar.reset();

		tar_index = _tar_index;
		y_index = _y_index;

		has_one = false;

		lincombs.clear();
		idx_t init_lin;
		for (wind_t j = 0; j < DIMENSION; ++j)
			if (tar.test(j)) {
				init_lin.emplace_back(j);
				masked_tar[j] = 1;
			}
		LinComb lc(init_lin);
		lincombs.emplace_back(lc);
		assert(lincombs.size() > 0);
		supports.clear();

		minimal_value = lc.value;
		depth = -1;
		value_bound = 0;
		depth_bound = -1;

		unavai_gs.clear();
		for (wind_t j = DIMENSION; j < DIMENSION + NL_NUM; ++j)
			if (tar.test(j))
				unavai_gs.emplace(j - DIMENSION);

		dist = init_lin.size() + unavai_gs.size() - 1;		

		if (dist == 0 && (unavai_gs.size() == 0)) {
			assert(unavai_gs.size() == 0);
			depth = 0;
			lincombs.clear();
		}

	}




	inline void append_supports(std::vector<LinComb> const & newlcs,
								std::vector<Bvector> const & basis
	) { 
		for (auto const & lc : newlcs) {
			assert(lc.value <= value_bound);
			assert(lc.depth <= depth_bound);
			auto const & lin = lc.lin;
			if (lin.size() == 0) {
				assert(newlcs.size() == 1);
				return;
			}
			for (wind_t i = 0; i < lin.size() - 1; ++i) {
				auto const & idx0 = lin[i];
				auto const & d0 = basis[idx0].imp_depth;

				for (wind_t j = i + 1; j < lin.size(); ++j) {
					auto const & idx1 = lin[j];
					auto const & d1 = basis[idx1].imp_depth;

					word64_t offset = dvalue(upper_log2(dvalue(d0) + dvalue(d1))) - dvalue(d0) - dvalue(d1);
					if (offset + lc.value <= value_bound) 
						supports.emplace(std::make_pair(idx0, idx1));

				}
			}
		}
	}



	void print() const { 
		std::cout << "target vector " << tar_index;
		if (y_index >= 0) std::cout << ", y_index " << y_index;
		std::cout << ": " << std::endl;
		vector_print(tar);

		if (used_by_gs.size() > 0) {
			std::cout << "used by:";
			for (auto const & gi : used_by_gs)
				std::cout << " g" << gi;
			std::cout << std::endl;
		}


		std::cout << "depth: " << depth
		          << ", minimal_value: " << minimal_value
		          << ", depth_bound: " << depth_bound
		          << ", value_bound: " << value_bound;
		if (has_one) std::cout << ", has one: " << has_one;
		std::cout << std::endl;

		std::cout << "gs:";
		for (auto const & gi : unavai_gs) std::cout << " " << gi;
		std::cout << std::endl;

		std::cout << "optimal dist: " << dist << std::endl;
		for (wind_t i = 0; i < lincombs.size(); ++i) {
			auto & lc = lincombs[i];
			std::cout << "[" << i << "]:";
			lc.print();
		}
		std::cout << "supports (" << supports.size() << "): ";
		for (auto & pa : supports) {
			std::cout << " ";
			pair_print(pa);
		}
		std::cout << std::endl << std::endl;
	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////


	void solve_by_enum(std::vector<idx_t> & lins, 
					   wind_t const tar_n_basis, 
					   std::vector<vec_t> const & basis, 
					   wind_t const not_fixed_basis, 
					   vec_t const & target) { //ok5.2

		std::unordered_map<vec_t, wind_t> bumap;
		for (wind_t i = 0; i < not_fixed_basis; ++i)
			bumap[basis[i]] = i;

		if (tar_n_basis == 1) {
			if (bumap.find(target) != bumap.end()) {
				idx_t lin = {bumap.at(target)};
				lins.emplace_back(lin);
			}
			return;
		}


		Enumerator enumer(tar_n_basis - 1, not_fixed_basis);

		while (enumer.step()) {
			auto sel_idxs = enumer.current();

			auto vec = target;
			for (auto const & idx : sel_idxs)
				vec ^= basis[idx];

			auto it = bumap.find(vec);
			if (it != bumap.end() && it->second > *sel_idxs.rbegin()) {
				auto lin = sel_idxs;
				lin.emplace_back(it->second);
				lins.emplace_back(lin);
			}
		}	

	}



	wind_t solve_by_gurobi(std::vector<idx_t> & lins, 
						   wind_t const tar_n_basis, 
						   std::vector<Bvector> const & basis, 
						   wind_t const not_fixed_basis,
						   vec_t const & target
	) { 

		assert(target.count() > 0);
		assert(tar_n_basis > 0);
		assert(not_fixed_basis <= basis.size());

		// Create an environment
		GRBEnv env = GRBEnv(true);
		// env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();
		// Create an empty model
		GRBModel model = GRBModel(env);


		std::vector<GRBVar> sel;
		GRBLinExpr n_sel = 0;
		GRBLinExpr n_fix = 0;
		for (wind_t i = 0; i < basis.size(); ++i) {
			auto v = model.addVar(0., 1., 0., GRB_BINARY);
			n_sel += v;
			if (i >= not_fixed_basis)
				n_fix += v;
			sel.emplace_back(v);
		}
		model.addConstr(n_sel <= tar_n_basis);
		model.addConstr(n_fix >= 1);


		// depth constraint
		GRBLinExpr sum_depth = 0;
		for (wind_t i = 0; i < basis.size(); ++i)
			sum_depth += sel[i] * dvalue(basis[i].imp_depth);
		model.addConstr(sum_depth <= value_bound);



		// linear combination
		for (wind_t j = 0; j < DIMENSION + NL_NUM; ++j) {

			// deal with target[j] bit
			GRBLinExpr l = 0;
			for (wind_t i = 0; i < sel.size(); ++i)
				if (basis[i].vector.test(j))
					l += sel[i];

			if (target.test(j)) l += 1;
			
			GRBVar m = model.addVar(0., (tar_n_basis + 1) / 2 + 1, 0., GRB_INTEGER);
			model.addConstr(l == 2 * m);
		}


		// paramters of gurobi

		model.set(GRB_IntParam_RINS, 0);
		model.set(GRB_IntParam_VarBranch, 2);
		model.set(GRB_IntParam_MIPFocus, 3);

		model.set(GRB_IntParam_Threads, GUROBI_THREADS);

		// model.set(GRB_DoubleParam_MIPGap, GRB_INFINITY);
		model.set(GRB_IntParam_PoolSearchMode, 2);
		model.set(GRB_IntParam_PoolSolutions, 1000);

		// model.set(GRB_IntParam_Presolve, PRESOLVE_LEVEL);

		model.setObjective(n_sel, GRB_MINIMIZE);

		model.optimize();



		int solCount = model.get(GRB_IntAttr_SolCount);

		if (solCount == 0)
			return tar_n_basis;


		wind_t new_lincomb_num = tar_n_basis;
		for (int si = 0; si < solCount; ++si) {
		    model.set(GRB_IntParam_SolutionNumber, si);

		    idx_t lin;
		    for (int i = 0; i < sel.size(); ++i)
		        if (sel[i].get(GRB_DoubleAttr_Xn) > 0.5)
		        	lin.emplace_back(i);

		    if (new_lincomb_num > lin.size()) {
		    	new_lincomb_num = lin.size();
		    	lins.clear();
		    }

		    if (new_lincomb_num == lin.size())
		    	lins.emplace_back(lin);
		}
		assert(new_lincomb_num > 0);

		return new_lincomb_num;
	}




  

	void update(std::vector<Bvector> const & basis, 
				pair_t const & sel_pair,
				bool const show_flag
	) { 

		assert(dist > 0);

		std::string st = "Target " + std::to_string(tar_index) + ":";
		if (tar_index < 10) st += " ";
		if (tar_index < 100) st += " ";
		st += " \tdist = " + std::to_string(dist) + ",";
		if (dist < 10) st += " ";
		st += " \tbasis = " + std::to_string(basis.size()) + ", \tlins = " + std::to_string(lincombs.size());
		if (lincombs.size() < 10) st += " ";
		if (lincombs.size() < 100) st += " ";



		bool display = show_flag;
		const wind_t bidx = basis.size() - 1;
		std::vector<LinComb> newlcs;

		if (supports.find(sel_pair) != supports.end()) {
		
			--dist;

			if (dist == 0) {
				assert(unavai_gs.size() == 0);
				assert(tar == basis.rbegin()->vector);

				depth = basis[bidx].imp_depth;
				assert(depth <= depth_bound);
				minimal_value = dvalue(depth);
				value_bound = minimal_value;

				st += " \tdecrease " + std::to_string(lincombs.size()) + " => 0";
				st += " (dist: " +  std::to_string(dist) + ")";
				
				lincombs.clear();
				supports.clear();

			} else {

				word64_t min_ss = word64_t(-1);
				for (auto const & lc : lincombs) {
					idx_t lin;
					if (!lc.reduce_by_pair(lin, sel_pair, bidx)) continue;

					assert(lin.size() > 0);
					LinComb nlc(lin, basis);
					
					if (min_ss > nlc.value)
						min_ss = nlc.value;

					if (nlc.value <= value_bound)
						newlcs.emplace_back(nlc);
				}
				assert(newlcs.size() > 0);
				assert(min_ss != word64_t(-1));
				assert(min_ss >= minimal_value);

				minimal_value = min_ss;

				if (lincombs.size() - newlcs.size() > 0) {
					st += " \tdecrease " + std::to_string(lincombs.size() - newlcs.size());
					st += " => " + std::to_string(newlcs.size()) + " (dist: " +  std::to_string(dist) + ")";
				}
				
				lincombs = newlcs;
				supports.clear();
				append_supports(lincombs, basis);
			}			

		} else {

			if (masked_tar.count() > 0) {
				Timer tt;
				tt.start();

				std::vector<idx_t> newlins;
				wind_t new_lincomb_num = -1;
				if (dist + 1 - 1 <= ENUM_BOUND && 0) {
					assert(false);

				} else {
					assert(dist + 1 > unavai_gs.size());
					new_lincomb_num = solve_by_gurobi(newlins, dist + 1 - unavai_gs.size(), basis, basis.size() - 1, masked_tar);
					assert(new_lincomb_num + unavai_gs.size() <= dist + 1);
				}

				std::string timestr = std::to_string(tt.used_time());
				timestr.resize(timestr.size() - 3);
				st += " \ttime: " + timestr;

				if (newlins.size() > 0) {

					if (new_lincomb_num + unavai_gs.size() == dist + 1) {
						st += ", found " + std::to_string(newlins.size()) + " ~~~~~~~~ ";

						word64_t min_ss = minimal_value;
						for (auto const & lin : newlins) {
							LinComb nlc(lin, basis);

							if (min_ss > nlc.value)
								min_ss = nlc.value;

							assert(nlc.value <= value_bound);
							newlcs.emplace_back(nlc);
						}

						minimal_value = min_ss;						

						lincombs.insert(lincombs.end(), newlcs.begin(), newlcs.end());
						append_supports(newlcs, basis);

					} else {

						dist = new_lincomb_num + unavai_gs.size() - 1;
						st += ", found " + std::to_string(newlins.size()) + " dist = " + std::to_string(dist) + " ~~~~~~~~ ";					

						word64_t min_ss = word64_t(-1);
						for (auto const & lin : newlins) {
							LinComb nlc(lin, basis);

							if (min_ss > nlc.value)
								min_ss = nlc.value;

							assert(nlc.value <= value_bound);
							newlcs.emplace_back(nlc);
						}
						assert(newlcs.size() > 0);
						assert(min_ss != word64_t(-1));	

						minimal_value = min_ss;

						lincombs = newlcs;
						supports.clear();
						append_supports(lincombs, basis);
					}
				
				} else {
					if (tt.used_time() < 1.)
						display = false;	
				}
			} else {
				assert(unavai_gs.size() > 0);
				assert(minimal_value == 0);
				assert(depth == -1);
				st += " masked_tar.count() == 0";
			}
		}


		if (display)
			std::cout << st << std::endl;

	}


};

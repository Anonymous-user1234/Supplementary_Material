

	#define SHOW_DETAILS_BP (0)


	wind_t predict_supports(pair_t pa, std::vector<Target> targets){
		std::vector<Bvector> basis_set;
		for(auto const & b: basis){
			basis_set.emplace_back(b);
		}

		std::vector<Target> tmp_targets; 

		for(auto  &t: targets){
			tmp_targets.emplace_back(t);
		}
		Bvector newb;

		newb.linear_op(basis_set[pa.first], basis_set[pa.second]);
		basis_set.emplace_back(newb);

		for (wind_t ti = 0; ti < tmp_targets.size(); ++ti) {
            if (tmp_targets[ti].dist <= 0) continue;

            threadpool.push([this, ti, &tmp_targets, &basis_set, &pa]() {
                tmp_targets[ti].update(basis_set, pa, 0);
            });
        }
        threadpool.wait_work();
		        

		
		int dist_max = 0;
		int dist_min = 1000;
		int unfinish = -1;
		
		int support_min = 1000;
		for (auto & T:tmp_targets){
			if(T.dist > 0 || T.dist == -1){
				unfinish ++;
			}

			if (T.unavai_gs.size() > 0) continue;
			if(dist_max < T.dist){
				dist_max = T.dist;
			}
			if(dist_min > T.dist && T.dist > 0){
				dist_min = T.dist;
			}
			if(T.supports.size() < support_min && T.supports.size()>0){
				support_min = T.supports.size();
			}
			
		}
		int dist_bound = DIMENSION +NL_NUM;

		#if TARGET_COND == 2
			dist_bound = dist_min + 3;
		#endif


		int res = 0;

		#if STRATEGY == 3
			int support_num = 0;

			std::map<pair_t, wind_t> votes;
			for (auto & T : tmp_targets) {
				#if TARGET_COND == 2
					if (T.dist > dist_bound) continue;
				#endif
				if (T.unavai_gs.size() > 0) continue;
				support_num += T.supports.size();
				for (auto const & pa : T.supports)
					votes[pa] ++;
			}

			res=support_num-votes.size();
		#endif

		return res;
	}
	

	pair_t bp_strategy(bool const show_flag, bool const greedy) { 
		int display = 0;
		// voting best pairs.
		std::map<pair_t, wind_t> votes;
		int dist_max = 0;
		int dist_min = 1000;
		int unfinish = -1;
		if(display) printf("ops size:%ld\n", ops.size());
		int support_min = 1000;
		for (auto const & T:targets){
			if(T.dist > 0 || T.dist == -1){
				unfinish ++;
			}

			if (T.unavai_gs.size() > 0) continue;
			if(dist_max < T.dist){
				dist_max = T.dist;
			}
			if(dist_min > T.dist && T.dist > 0){
				dist_min = T.dist;
			}
			if(T.supports.size() < support_min && T.supports.size()>0){
				support_min = T.supports.size();
			}
			if(display){
				if(T.dist != 0){
					printf("support size:%ld, dist:%d\n", T.supports.size(), T.dist);
					T.print();
				}
			}
		}
		if(display) printf("\n");
		
		int dist_bound = DIMENSION +NL_NUM;

		#if TARGET_COND == 2
			dist_bound = dist_min + 3;
		#endif

		for (auto const & T : targets) {
#if TARGET_COND == 2
			if (T.dist > dist_bound) continue;
#endif
			if (T.unavai_gs.size() > 0) continue;
			wind_t w = 1;
			if (T.y_index == -1)
				w = T.used_by_gs.size();

			for (auto const & pa : T.supports)
				votes[pa] += w;
		}
	

		wind_t highest_val = 0;
		std::vector<pair_t> max_dist_pairs;
		for (auto const & it : votes) {

			if (SHOW_DETAILS) {
				pair_print(it.first, false);
				std::cout << ": " << it.second << ", ";
			}

			if (highest_val < it.second) {
				highest_val = it.second;
				max_dist_pairs.clear();
			}

			if (highest_val == it.second)
				max_dist_pairs.emplace_back(it.first);
		}
    	if (SHOW_DETAILS)
    		std::cout << std::endl << std::endl;



		wind_t best_score = -1;
		std::vector<pair_t> bp_greedy_pairs;

		for (auto const & pa : max_dist_pairs) {
			wind_t sc = 0;
#if STRATEGY==1 //MN
			for (auto const & T : targets) {
				if (T.dist <= 0) continue;
				if (T.unavai_gs.size() > 0) continue;
#if TARGET_COND == 2
			if (T.dist > dist_bound) continue;
#endif
				wind_t w = 1;
				if (T.y_index == -1)
					w = T.used_by_gs.size();

				assert(T.dist > 0);
				if (T.supports.find(pa) != T.supports.end()){
					sc += (T.dist * T.dist - (T.dist - 1) * (T.dist - 1)) * w;
				}
			}
			if (best_score == -1 || best_score > sc) {
				best_score = sc;
				bp_greedy_pairs.clear();
			}
			if (sc == best_score)
				bp_greedy_pairs.emplace_back(pa);
#endif


#if STRATEGY== 2 //ASM
			for (auto const & T : targets) {
				if (T.dist <= 0) continue;
				if (T.unavai_gs.size() > 0) continue;
#if TARGET_COND == 2
			if (T.dist > dist_bound) continue;
#endif
				wind_t w = 1;
				if (T.y_index == -1)
					w = T.used_by_gs.size();

				assert(T.dist > 0);
				if (T.supports.find(pa) != T.supports.end()){
					sc += T.supports.size()/T.lincombs.size(); 
				}
			}
			if (best_score == -1 || best_score > sc) {
				best_score = sc;
				bp_greedy_pairs.clear();
			}
			if (sc == best_score)
				bp_greedy_pairs.emplace_back(pa);

#endif



#if STRATEGY==3 //SOM
			sc = predict_supports(pa, targets)/10;
			if (best_score == -1 || best_score < sc) {
				best_score = sc;
				bp_greedy_pairs.clear();
			}
			if (sc == best_score)
				bp_greedy_pairs.emplace_back(pa);

#endif
			
		}


		wind_t greedy_num = bp_greedy_pairs.size();

		if (show_flag) {
			std::cout << "highest = " << highest_val;
			std::cout << ", candidates " << max_dist_pairs.size();
			std::cout << " ---> " << bp_greedy_pairs.size() << std::endl;
		}
	

		if (max_dist_pairs.size() == 0)
			return std::make_pair(-1, -1);

		if (greedy)
			return bp_greedy_pairs[rndint(0, bp_greedy_pairs.size() - 1)];
		else
			return max_dist_pairs[rndint(0, max_dist_pairs.size() - 1)];
	}

	





	bool append_newvector(select_t const & sel, bool const show_flag) { 

		// step 1: add new basis
		ops.emplace_back(op_t({wind_t(basis.size()), sel.sel_pair.first, sel.sel_pair.second, sel.type, sel.tar_idx}));

		if (basis_map.find(sel.newb) != basis_map.end()) {
			std::cout << "this should not happen (append_newvector) !!!!" << std::endl;
			return false;
		}
		basis_map[sel.newb] = basis.size();
		basis.emplace_back(sel.newb);
		
		if (max_imp_depth < sel.newb.imp_depth)
			max_imp_depth = sel.newb.imp_depth;

		if (max_ad_depth < sel.newb.ad_depth)
			max_ad_depth = sel.newb.ad_depth;


		if (sel.type == TYPE_AND || sel.type == TYPE_OR) {
			assert(sel.newb.vector.count() == 1);
			vector_mask |= sel.newb.vector;
		}

		
		if (show_flag && SHOW_DETAILS_BP)
			std::cout << "add basis " << vector_to_uint64(sel.newb.vector) << std::endl;


		if (show_flag && SHOW_DETAILS_BP) {
			basis_print();
			PRINT_ENTER()
		}

		// step 2: build linear representations
		if (sel.type != TYPE_XOR) {

			assert(sel.newb.vector.count() == 1);
			assert(sel.newb == *basis.rbegin());

			wind_t gi = 0;
			while (gi < NL_NUM && !sel.newb.vector.test(gi + DIMENSION))
				++gi;
			assert(gi < NL_NUM);



			for (wind_t ti = 0; ti < targets.size(); ++ti) {
				if (targets[ti].dist == 0 && targets[ti].unavai_gs.size() == 0) continue;

				auto & T = targets[ti];

				if (T.unavai_gs.find(gi) == T.unavai_gs.end()) continue;

				T.masked_tar = T.tar & vector_mask;
				T.unavai_gs.erase(gi);

				T.value_bound += dvalue(gys[gi].set_bound);
				assert(T.value_bound <= dvalue(T.depth_bound));
				if (T.unavai_gs.size() == 0)
					assert(T.value_bound == dvalue(T.depth_bound));

				for (auto & lc : T.lincombs) {

					lc.lin.emplace_back(basis.size() - 1);
					lc.value += dvalue(sel.newb.imp_depth);
					lc.depth = upper_log2(lc.value);

					assert(lc.value <= T.value_bound);
					assert(lc.depth <= T.depth_bound);
				}
				T.minimal_value += dvalue(sel.newb.imp_depth);

				T.supports.clear();
				T.append_supports(T.lincombs, basis);

				if (T.unavai_gs.size() > 0) continue;

				threadpool.push([this, ti]() {

					auto & T = targets[ti];

					std::string st = "Target " + std::to_string(T.tar_index) + ":";
					if (T.tar_index < 10) st += " ";
					if (T.tar_index < 100) st += " ";
					st += " \tdist = " + std::to_string(T.dist) + ",";
					if (T.dist < 10) st += " ";
					st += " \tbasis = " + std::to_string(basis.size()) + ", \tlins = " + std::to_string(T.lincombs.size());
					if (T.lincombs.size() < 10) st += " ";
					if (T.lincombs.size() < 100) st += " ";

					Timer tt;
					tt.start();
					
					std::vector<idx_t> newlins;
					assert(T.tar == T.masked_tar);
					wind_t new_lincomb_num = T.solve_by_gurobi(newlins, T.dist + 1, basis, 0, T.masked_tar);
					assert(new_lincomb_num <= T.dist + 1);
				

					std::string timestr = std::to_string(tt.used_time());
					timestr.resize(timestr.size() - 3);
					st += " \ttime: " + timestr;

					
					assert(newlins.size() > 0);

					T.dist = new_lincomb_num - 1;
					st += ", found " + std::to_string(newlins.size()) + " dist = " + std::to_string(T.dist) + " ~~~~~~~~ ";					

					word64_t min_ss = word64_t(-1);
					T.lincombs.clear();
					for (auto const & lin : newlins) {
						LinComb nlc(lin, basis);

						if (min_ss > nlc.value)
							min_ss = nlc.value;

						assert(nlc.value <= T.value_bound);
						T.lincombs.emplace_back(nlc);
					}
					assert(T.lincombs.size() > 0);
					assert(min_ss != word64_t(-1));	

					T.minimal_value = min_ss;

					T.supports.clear();
					T.append_supports(T.lincombs, basis);							
					
				});		
			}
			threadpool.wait_work();

		} else {

			for (wind_t ti = 0; ti < targets.size(); ++ti) {
				if (targets[ti].dist == 0) continue;

				threadpool.push([this, ti, sel, show_flag]() {
					targets[ti].update(basis, sel.sel_pair, show_flag);
				});		
			}
			threadpool.wait_work();
		}

		return true;
	}




	bool update_result(wind_t & best_gxor) { 

		bool updated = false;

		my_local_reconstruct();


		ops_check();
		printf("Finish local optimization checking! best_gxor:%d, ops size:%ld\n", best_gxor, ops.size());
		if (best_gxor > ops.size()) {

			updated = true;
			ops_write();

			ops_print();
			best_gxor = ops.size();

		}
		printf("lr done!");
		return updated;
	}



	///////////////////////////////// universal bp //////////////////////////////////////////


	bool universal_bp(bp_param_t const & param, wind_t & best_gxor) { //

		bool show_flag = param.display;	

		if (show_flag) {
			PRINT_ENTER()
			PRINTF_STAMP("start universal BP algorithm\n");
			param_print(param);
			PRINTF_STAMP("status: targets %ld, basis %ld, ops %ld\n", targets.size(), basis.size(), ops.size());
		}
		// check inputs
		bp_inputs_check();




		// start
		Timer tt;
		tt.start();


		wind_t n_targets = 0;
		for (auto const & T : targets) n_targets += (T.dist != 0);
		assert(n_targets > 0);


		wind_t last_dist = -1;
		wind_t iter = 0;
		ops.clear();


		// main loops
		while (true) {	

			// update statistics
			wind_t n_targets = 0;
			wind_t new_dist = 0;
			wind_t lbs_num = 0;
			wind_t supp_num = 0;
			for (auto const & T : targets) {
				n_targets += (T.dist != 0);
				new_dist += T.dist;
				lbs_num += T.lincombs.size();
				supp_num += T.supports.size();	
			}


			if (n_targets == 0) break;



			// display information
			if (show_flag) {
				std::string st = "";
				for (wind_t ti = 0; ti < targets.size(); ++ti) {
					auto const & T = targets[ti];
					if (ti % 32 == 0) {
						if (ti > 0) st += "]";
						st += "\n\t\t[";
					}
					st += " ";
					if (T.dist == 0) st += pad_to_len(DIM_LEN);
					else if (T.dist == -1) st += "*" + pad_to_len(DIM_LEN - 1);
					else st += pad_to_len(T.dist, DIM_LEN);
				}

				PRINT_ENTER()
				PRINTF_STAMP("iter %d, targets %d, dist %d ===> %d [%d], xor upper bound %d, #lincombs %d, #supp %d", iter, n_targets, last_dist, new_dist, last_dist - new_dist, iter + new_dist, lbs_num, supp_num);
				std::cout << st << "]" << std::endl;
				last_dist = new_dist;
			}
			++iter;



			if (SHOW_DETAILS_BP) {
				for (auto const & T : targets) 
					if (T.dist != 0) T.print();
				PRINT_ENTER()
			}
	

			for (auto const & T : targets) {
				for (auto const & lc : T.lincombs)
					assert(lc.value <= T.value_bound);
				assert(T.value_bound <= dvalue(T.depth_bound));
			}



			select_t sel;

			// step 1: easy case, some dist is 1, or nonlinear operation is activated!
			{
				// type I
				for (auto const & T : targets) {
					if (T.dist != 1 || T.unavai_gs.size() > 0) continue;

					sel.tar_idx = T.tar_index;
					auto const & lin = T.lincombs[0].lin;
					assert(lin.size() == 2);
					sel.sel_pair = std::make_pair(lin[0], lin[1]);
					assert(T.supports.find(sel.sel_pair) != T.supports.end());
					
					sel.newb.linear_op(basis[lin[0]], basis[lin[1]]);					

					if (show_flag) 
						std::cout << "easy case type I !" << std::endl;

					sel.easy_case = true;
					sel.type = TYPE_XOR;

					break;
				}



				// type II
				for (wind_t gi = 0; gi < NL_NUM; ++gi) {

					if (sel.easy_case) break;

					Bvector g(DIMENSION + gi);

					bool found = false;
					for (auto const & b : basis)
						if (b.vector == g.vector) {
							found = true;
							break;
						}
					if (found) continue;

					assert(gys[gi].ts.size() == 2);
					auto sub0_ti = gys[gi].ts[0];
					auto sub1_ti = gys[gi].ts[1];

					if (targets[sub0_ti].dist != 0 || targets[sub0_ti].unavai_gs.size() != 0) continue;
					if (targets[sub1_ti].dist != 0 || targets[sub1_ti].unavai_gs.size() != 0) continue;


					// hit
					wind_t s0 = -1;
					for (wind_t i = 0; i < basis.size(); ++i)
						if (basis[i].vector == targets[sub0_ti].tar && basis[i].imp_depth < gys[gi].set_bound) {
							s0 = i;
							break;
						}
					assert(s0 >= 0);

					wind_t s1 = -1;
					for (wind_t i = 0; i < basis.size(); ++i)
						if (basis[i].vector == targets[sub1_ti].tar && basis[i].imp_depth < gys[gi].set_bound) {
							s1 = i;
							break;
						}
					assert(s1 >= 0);
					assert(s0 != s1);
					if (s0 < s1) sel.sel_pair = std::make_pair(s0, s1);
					else sel.sel_pair = std::make_pair(s1, s0);

					assert(basis[s0].imp_depth < gys[gi].set_bound);
					assert(basis[s1].imp_depth < gys[gi].set_bound);

					sel.newb = g;
					sel.newb.nonlinear_depth(basis[s0], basis[s1]);
					assert(sel.newb.imp_depth <= gys[gi].set_bound);

					sel.tar_idx = get_target_index(sel.newb.vector);

					if (show_flag) 
						std::cout << "easy case type II !" << std::endl;

					sel.easy_case = true;
					sel.type = gys[gi].type;

					break;
				}
			}



			// step 2: voting best pairs.
			if (!sel.easy_case) {

				sel.sel_pair = bp_strategy(show_flag, param.greedy);
				if (sel.sel_pair.first == -1) {
					std::cout << "no pair left!!! break " << std::endl;
					return false;
				}

				sel.newb.linear_op(basis[sel.sel_pair.first], basis[sel.sel_pair.second]);

				sel.tar_idx = -1;
				sel.type = TYPE_XOR;

			}


			if (show_flag)
				sel.print();


			// step 3: append new vector to basis
			if (!append_newvector(sel, show_flag))
				return false;

		}


		if (show_flag) std::cout << std::endl;



		// convert xor to nxor
		for (wind_t oi = 0; oi < ops.size(); ++oi) {
			auto & op = ops[oi];
			if (op.tar_idx < 0) continue;

			auto ti = op.tar_idx;
			op.tar_idx = get_y_index(basis[op.dst].vector);
			if (op.tar_idx < 0) continue;

			if (!targets[ti].has_one) continue;


			if (op.type != TYPE_XOR && op.type != TYPE_NXOR) {
				std::cout << "y[" << op.tar_idx << "] cannot be generated by xor/nxor" << std::endl;
				return false;
			}

			if (op.type == TYPE_XOR) op.type = TYPE_NXOR;
			else if (op.type == TYPE_NXOR) op.type = TYPE_XOR;

			for (wind_t oj = oi + 1; oj < ops.size(); ++oj) {
				if (ops[oj].src_0 == op.dst || ops[oj].src_1 == op.dst) {
					if (ops[oj].type != TYPE_XOR && ops[oj].type != TYPE_NXOR) {
						std::cout << "y[" << op.tar_idx << "] cannot be used latter" << std::endl;
						return false;
					}
					if (ops[oj].type == TYPE_XOR) ops[oj].type = TYPE_NXOR;
					else if (ops[oj].type == TYPE_NXOR) ops[oj].type = TYPE_XOR;
				}
			}


		}


		ops_check();


		if (ops.size() < best_gxor) {
			best_gxor = ops.size();
			ops_print();
			ops_write();
		}


		if (ops.size() < best_gxor + 7) {
			printf("best_gxor:%d, ops size:%ld\n", best_gxor, ops.size());

			update_result(best_gxor);
			printf("In local opti: best_gxor:%d, ops size:%ld\n", best_gxor, ops.size());
		}


		PRINTF_STAMP("xor count %ld, max_ad %d, max_imp_depth %d, time used: %.3f\n", ops.size() - NL_NUM, max_ad_depth, max_imp_depth, tt.used_time());

		return true;
		
	}




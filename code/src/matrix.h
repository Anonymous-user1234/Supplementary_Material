

#ifndef BP_MATRIX_H
#define BP_MATRIX_H



#define VEC_ERASE(vec, ele) vec.erase(std::remove(vec.begin(), vec.end(), ele), vec.end());


#define SHOW_DETAILS_MATRIX (1)





class GY {

public:

	std::vector<wind_t> ts;
	std::vector<vec_t> t_vec;
    wind_t set_bound = -1;
    wind_t ad_depth = -1;

    wind_t type = TYPE_NULL;

    
    inline void print(wind_t const & idx, wind_t const & y_index = -1) const { 
        if (idx < NL_NUM) {
        	assert(ts.size() == 2);
            std::cout << "g" << idx;
            if (idx < 10) std::cout << " ";
            std::cout << " from t" << ts[0] << ",";
            if (ts[0] < 10) std::cout << " ";
            std::cout << " t" << ts[1] << ",";
            if (ts[1] < 10) std::cout << " ";
        }
        else {
        	assert(y_index != -1);
            std::cout << "y" << y_index << " from t" << idx;
        }

        std::cout << "\t type " << type
                  << ",\t ad_depth " << ad_depth
                  << ", \tset to " << set_bound
                  << std::endl;
    }

};






class Bvector {

public:

	vec_t vector;

	wind_t imp_depth;
	wind_t ad_depth;

	Bvector() {}

	Bvector(wind_t const & idx, wind_t const _imp = 0, wind_t const _ad = 0) { 
		vector.reset();
		vector[idx] = 1;
		imp_depth = _imp;
		ad_depth = _ad;
	}

	

	inline void linear_op(Bvector const & a, Bvector const & b) { 
		vector = a.vector ^ b.vector;
		imp_depth = std::max(a.imp_depth, b.imp_depth) + 1;
		ad_depth = std::max(a.ad_depth, b.ad_depth);
	}

	inline void nonlinear_depth(Bvector const & a, Bvector const & b) {
		imp_depth = std::max(a.imp_depth, b.imp_depth) + 1;
		ad_depth = std::max(a.ad_depth, b.ad_depth) + 1;
	}

	bool operator==(Bvector const & other) const {
		return (vector == other.vector && imp_depth == other.imp_depth && ad_depth == other.ad_depth);
	}

	bool operator<(const Bvector& other) const {

        if (vector != other.vector) {

            for (int i = DIMENSION + NL_NUM - 1; i >= 0; --i) {
                if (vector[i] != other.vector[i]) {
                    return vector[i] < other.vector[i];
                }
            }
        }
        if (imp_depth != other.imp_depth) {
            return imp_depth < other.imp_depth;
        }
        return ad_depth < other.ad_depth;
    }
	void print(bool const & enter = false) const {
		vector_print(vector, false);
		std::cout << " depth [" << imp_depth << ", ad " << ad_depth << "]";
		if (enter) std::cout << std::endl;
	}

};


namespace std {
	template<>
	struct hash<Bvector> {
		size_t operator()(Bvector const & b) const {
			return (hash<vec_t>{}(b.vector) * 929) ^ (hash<wind_t>{}(b.imp_depth) * 733) ^ (hash<wind_t>{}(b.ad_depth) * 331);
		}
	};
}


struct BvectorComparator {
    bool operator()(const Bvector & a, const Bvector & b) const {
        if (a.vector.count() != b.vector.count()) {
            return a.vector.count() < b.vector.count();
        }
        
        for (wind_t i = 0; i < DIMENSION; ++i) {
            if (a.vector[i] != b.vector[i]) {
                return a.vector[i] < b.vector[i];  
            }
        }
        
        if (a.imp_depth != b.imp_depth) {
            return a.imp_depth < b.imp_depth;
        }
        return a.ad_depth < b.ad_depth;
    }
};



typedef std::set<Bvector, BvectorComparator> vecset_t;

typedef std::set<vec_t, VectorComparator> my_vecset_t;


typedef std::vector<Bvector> Bvvec_t;



//return true if the target vector can be generated.
inline bool extendible(std::unordered_map<Bvector, wind_t> const & gen_map, 
					   Bvector  & tar_vec) { 
	assert(tar_vec.vector.count() > 0);
	for (auto const & it : gen_map) {
		Bvector bnew;
		bnew.linear_op(it.first, tar_vec);
		for (auto const & res: gen_map){
			if(res.first.vector == bnew.vector && tar_vec.imp_depth > res.first.imp_depth && tar_vec.imp_depth > it.first.imp_depth){
				tar_vec.imp_depth = std::max(res.first.imp_depth, it.first.imp_depth)+1;

				return true;
			}
		}
			
	}
	return false;
}

inline bool activate_g(std::unordered_map<Bvector, wind_t> const & gen_map, Bvector  & tar_vec, std::vector<GY> & gys){
	for(int gpos = DIMENSION; gpos < DIMENSION + NL_NUM; gpos++){

		if(tar_vec.vector.test(gpos)){
			vec_t t0 = gys[gpos - DIMENSION].t_vec[0];
			vec_t t1 = gys[gpos - DIMENSION].t_vec[1];
			int app = 0;
			int depth = -1;
			for (auto const & it: gen_map){
				if(it.first.vector == t0 || it.first.vector == t1){
					app ++;
					if(it.first.imp_depth > depth){
						depth = it.first.imp_depth;
					}
				}
			}
			if(app == 2 && tar_vec.imp_depth > depth){
				tar_vec.imp_depth = depth + 1;
				return true;
			}
		}
	}
	
	return false;
}


bool extend_vectors(std::vector<Bvector> & useful_set, 
					Bvvec_t & candidates, std::vector<GY> & gys) { 

	bool display = 0;

	std::unordered_map<Bvector, wind_t> bmap;
	wind_t idx = 0;
	for (auto const & v : useful_set) bmap[v] = idx++;

	wind_t level = 0;

	while (true) {
		Bvvec_t added;

		for (auto  & test_b : candidates) {
			assert(test_b.vector.count() > 0);
			if(test_b.imp_depth == 0){
				added.emplace_back(test_b);
				continue;				
			}

			if (extendible(bmap, test_b))
				added.emplace_back(test_b);

			if (test_b.imp_depth != 0 && test_b.vector.count() == 1){
				if(activate_g(bmap, test_b, gys)){
					added.emplace_back(test_b);
				}
			}	
		}

		if (display)
			std::cout << "level " << level++ << ": " << added.size()
			          << ", \t current useful_set: " << useful_set.size() 
			          << ", \t todo: " << candidates.size() 
			          << std::endl;

		if (added.size() == 0)
			break;

		for (auto const & b : added) {
			bmap[b] = useful_set.size();
			useful_set.emplace_back(b);				
			candidates.erase(std::remove(candidates.begin(), candidates.end(), b), candidates.end());

		}
	}

	if (display)
		std::cout << "\t useful_set: " << useful_set.size() 
				  << ", \t remain: " << candidates.size() 
				  << std::endl << std::endl;

	return candidates.size() == 0;
}





struct select_t { 

    Bvector newb;
    pair_t sel_pair;

    wind_t tar_idx = -1;
    wind_t type = TYPE_NULL;

    bool easy_case = false;


	void print() const { 
	    std::cout << "selection: sel pair = ";
	    pair_print(sel_pair, false);
	    std::cout << ", tar_index = " << tar_idx
	              << ", type = ";

	    switch (type) {
	        case TYPE_XOR : { std::cout << "xor"; break; }
	        case TYPE_NXOR : { std::cout << "nxor"; break; }
	        case TYPE_AND : { std::cout << "and"; break; }
	        case TYPE_OR : { std::cout << "or"; break; }
	        default: {
	            assert(false);
	            break;
	        }
	    }

	    std::cout << ", easy case = " << easy_case
	              << std::endl;
	    newb.print(true);
	}

};





bool calc_key_len(std::vector<Bvector> & useful_set, Bvvec_t & affected, std::vector<GY> & gys) { 
		wind_t key_len = -1;

		assert(affected.size() > 0);
		assert(useful_set.size() > 0);

		if (affected.size() <= 1) {
			key_len = affected.size();
			return key_len == 1;
		}
		key_len = 0;
		for (auto const & to_test : affected) {

			auto tmp_useful_set = useful_set;
			tmp_useful_set.emplace_back(to_test);

			auto candidates = affected;
            candidates.erase(std::remove(candidates.begin(), candidates.end(), to_test), candidates.end());

			if (extend_vectors(tmp_useful_set, candidates, gys)) {
				key_len = 1;
				break;
			}
		}
		return key_len == 1;
	}




#include "enum.inl"
#include "target.inl"

class Sub_target {

public:

	Bvector target;

	wind_t min_len = -1;

	Bvvec_t useful_set;
	Bvvec_t to_delete;

	std::vector<LinComb> minimal_lincombs;
	std::set<pair_t> supports;




	Bvvec_t possible_keys(Bvvec_t const & basis_set, std::vector<GY> & gys) {

		Bvvec_t affected = basis_set;
		for (auto const & b : to_delete)
			VEC_ERASE(affected, b)

		useful_set.clear();
		extend_vectors(useful_set, affected, gys);

		assert(affected.size() > 0);
		assert(useful_set.size() > 0);

		if (affected.size() <= 1)
			return affected;

		Bvvec_t keys;
		for (auto const & to_test : affected) {

			auto tmp_useful_set = useful_set;
			tmp_useful_set.emplace_back(to_test);

			auto candidates = affected;
            candidates.erase(std::remove(candidates.begin(), candidates.end(), to_test), candidates.end());

			if (extend_vectors(tmp_useful_set, candidates, gys))
				keys.emplace_back(to_test);
		}

		return keys;
	}


	inline void append_supports(std::vector<LinComb> const & newlcs, std::vector<Bvector> const & useful_set,
								std::vector<Bvector> const & basis
	) { 
		for (auto const & lc : newlcs) {
			assert(lc.value <= dvalue(target.imp_depth));
			assert(lc.depth <= target.imp_depth);
			auto const & lin = lc.lin;
			if (lin.size() == 0) {
				assert(newlcs.size() == 1);
				return;
			}
			for (wind_t i = 0; i < lin.size() - 1; ++i) {
				auto const & idx0_in_usefulset = lin[i];
				wind_t idx0 = -1;
				auto it = std::find(basis.begin(), basis.end(), useful_set[idx0_in_usefulset]);
				
				if (it != basis.end()) {
					idx0 = std::distance(basis.begin(), it);

				}
				assert(idx0 > -1);
				auto const & d0 = basis[idx0].imp_depth;

				for (wind_t j = i + 1; j < lin.size(); ++j) {
					auto const & idx1_in_usefulset = lin[j];
					wind_t idx1 = -1;
					auto it = std::find(basis.begin(), basis.end(), useful_set[idx1_in_usefulset]);
					
					if (it != basis.end()) {
						idx1 = std::distance(basis.begin(), it);

					}
					assert(idx1 > -1);


					auto const & d1 = basis[idx1].imp_depth;

					word64_t offset = dvalue(upper_log2(dvalue(d0) + dvalue(d1))) - dvalue(d0) - dvalue(d1);
					if (offset + lc.value <= dvalue(target.imp_depth)) 
						supports.emplace(std::make_pair(idx0, idx1));
				}
			}
		}
	}

};



// store information when deleting this vector
class Delete_vector {

public:

	Bvvec_t useful_set;
	Bvvec_t affected;
	wind_t key_len = -1;

	
};


class Matrix {

public:

	std::vector<word64_t> table;
	
	std::vector<Target> targets;  
	std::unordered_map<vec_t, wind_t> target_map; 
	std::unordered_map<vec_t, wind_t> y_map; 

	std::vector<GY> gys;

	std::vector<Bvector> basis;
	std::unordered_map<Bvector, wind_t> basis_map;

	wind_t max_imp_depth = -1; 
	wind_t max_ad_depth = -1; 
	vec_t vector_mask;

	std::vector<op_t> ops;

	std::unordered_map<Bvector, Delete_vector> del_vec_infos; 

	wind_t threads;
	thread_pool::thread_pool threadpool;

	Matrix & operator=(const Matrix & mat) {

		table = mat.table;

		targets = mat.targets;
		target_map = mat.target_map;
		y_map = mat.y_map;

		gys = mat.gys;

		basis = mat.basis;
		basis_map = mat.basis_map;

		max_imp_depth = mat.max_imp_depth;
		max_ad_depth = mat.max_ad_depth;
		vector_mask = mat.vector_mask;

		ops = mat.ops;

		
		threads = WORKERS;
		threadpool.resize(threads);	

		return *this;
	}

	Matrix() {		
		threads = WORKERS;
		threadpool.resize(threads);	
	}

	Matrix(const Matrix & mat) {
		*this = mat;
	}





	Matrix( std::string const & sbox_path ) { 
		table.clear();

		std::ifstream input(sbox_path);
		std::string element;

		while (getline(input, element))
			table.emplace_back(std::stoi(element, nullptr, 16));
		assert(table.size() == (word64_t(1) << DIMENSION));

	    if (SHOW_DETAILS_MATRIX) table_print();

	    threads = WORKERS;
	    threadpool.resize(threads);
	}





	#include "ios.inl"



	inline bool check_increase_ti(wind_t const & rise_ti) const { 

		bool display = 0;

		assert(targets.size() - DIMENSION >= 0);

		std::vector<wind_t> ts_depth(NL_NUM * 2 + DIMENSION);
		std::vector<wind_t> gys_depth;

		for (wind_t j = 0; j < ts_depth.size(); ++j)
			ts_depth[j] = -1;

		for (wind_t gi = 0; gi < NL_NUM; ++gi) {
			for (wind_t j = 0; j < gys[gi].ts.size(); ++j) {
				auto ti = gys[gi].ts[j];
				if (targets[ti].depth_bound >= 0)
					ts_depth[gi * 2 + j] = targets[ti].depth_bound + (ti == rise_ti);
			}
		}


		for (wind_t j = 0; j < ts_depth.size(); ++j) {

			if (display && (j == NL_NUM * 2)) std::cout << "-------------" << std::endl;

			if (ts_depth[j] == -1) {

				wind_t ti;
				if (j < NL_NUM * 2) ti = gys[j / 2].ts[j % 2];
				else ti = gys[j - NL_NUM * 2 + NL_NUM].ts[0];


				auto const & T = targets[ti];

				word64_t ss = T.minimal_value;
				for (auto const & gi : T.unavai_gs) {
					assert(gi < gys_depth.size());
					ss += dvalue(gys_depth[gi]);
				}
				ts_depth[j] = upper_log2(ss);

			}


			// set gys lower bound
			if (ts_depth[j] > DEPTH_LIMIT)
				return false;


			// merge to generat gi
			if (j < NL_NUM * 2 && j % 2 == 1) {
				assert(ts_depth[j - 1] >= 0);
				assert(ts_depth[j] >= 0);

				wind_t est_d = std::max(ts_depth[j - 1], ts_depth[j]) + 1;
				if (est_d > DEPTH_LIMIT)
					return false;
				gys_depth.emplace_back(est_d);
			}
		}


		return true;
	}




	// check whether ti's set_bound can be increased.
	inline bool check_target_bounds() const { 

		assert(targets.size() - DIMENSION >= 0);

		std::vector<wind_t> ts_depth(NL_NUM * 2 + DIMENSION);
		std::vector<wind_t> gys_depth;

		for (wind_t j = 0; j < ts_depth.size(); ++j)
			ts_depth[j] = -1;

		for (wind_t gi = 0; gi < NL_NUM; ++gi) {
			for (wind_t j = 0; j < gys[gi].ts.size(); ++j) {
				auto ti = gys[gi].ts[j];
				assert(targets[ti].depth_bound >= 0);
				ts_depth[gi * 2 + j] = targets[ti].depth_bound;
			}
		}
		for (wind_t gi = NL_NUM; gi < NL_NUM + DIMENSION; ++gi) {
			auto ti = gys[gi].ts[0];
			assert(targets[ti].depth_bound >= 0);
			ts_depth[NL_NUM + gi] = targets[ti].depth_bound;
		}

		for (auto const & g : gys)
			gys_depth.emplace_back(g.set_bound);

		for (wind_t j = 0; j < ts_depth.size(); ++j) {

			assert(ts_depth[j] >= 0);

			wind_t ti;
			if (j < NL_NUM * 2) ti = gys[j / 2].ts[j % 2];
			else ti = gys[j - NL_NUM * 2 + NL_NUM].ts[0];

				
			auto const & T = targets[ti];

			word64_t ss = T.minimal_value;
			for (auto const & gi : T.unavai_gs) {
				assert(gi < gys_depth.size());
				ss += dvalue(gys_depth[gi]);
			}
			assert(upper_log2(ss) <= ts_depth[j]);


			// merge to generat gi
			if (j < NL_NUM * 2 && j % 2 == 1) {
				assert(ts_depth[j - 1] >= 0);
				assert(ts_depth[j] >= 0);

				assert(std::max(ts_depth[j - 1], ts_depth[j]) + 1 == gys_depth[j / 2]);
			}
		}

		std::cout << "bounds check pass ~~~~" << std::endl;

		return true;
	}




	// update arget.depth_bounds and then gys's set bound
	void set_target_bounds() { 

		// step 1: initialize
		for (auto & T : targets) {
			if (T.dist == 1 && T.unavai_gs.size() == 0)
				T.depth = 1;

			if (T.depth == 0 || T.depth == 1)
				T.depth_bound = T.depth;
			else
				T.depth_bound = -1;

			if (T.y_index >= 0)
				T.depth_bound = DEPTH_LIMIT;
		}

		// step 2: increase the smallest one, and check		
		for (auto & T : targets) {
			std::cout << "test t" << T.tar_index << std::endl;
			{
				word64_t ss = T.minimal_value;
				for (auto const & gi : T.unavai_gs) {
					assert(gys[gi].set_bound >= 0);
					ss += dvalue(gys[gi].set_bound);
				}
				T.depth_bound = upper_log2(ss);
			}

			while (T.depth_bound < DEPTH_LIMIT && check_increase_ti(T.tar_index))
				++T.depth_bound;

			std::cout << "t" << T.tar_index << ": depth_bound = " << T.depth_bound << std::endl;

			assert(T.depth_bound <= DEPTH_LIMIT);

			for (auto const & gi : T.used_by_gs) {
				auto & g = gys[gi];
				if (targets[g.ts[0]].depth_bound >= 0 && targets[g.ts[1]].depth_bound >= 0)
					g.set_bound = std::max(targets[g.ts[0]].depth_bound, targets[g.ts[1]].depth_bound) + 1;
			}


		}


		gys_print();
			
		for (wind_t gi = 0; gi < NL_NUM; ++gi) {
			std::cout << "g" << gi << ": ";
			for (auto const & ti : gys[gi].ts)
				std::cout << "\t " << targets[ti].depth_bound << "(" << ti << ")";
			std::cout << std::endl;
		}

		for (wind_t gi = NL_NUM; gi < NL_NUM + DIMENSION; ++gi) {
			auto ti = gys[gi].ts[0];
			std::cout << "y" << targets[ti].y_index << ": ";			
			std::cout << "\t " << targets[ti].depth_bound << "(" << ti << ")";
			std::cout << std::endl;	
		}

		check_target_bounds();


		// set targets value bounds
		for (auto & T : targets) {
			T.value_bound = dvalue(T.depth_bound);
			for (auto const & gi : T.unavai_gs) {
				assert(T.value_bound >= dvalue(gys[gi].set_bound));
				T.value_bound -= dvalue(gys[gi].set_bound);
			}

			assert(T.supports.size() == 0);
			T.append_supports(T.lincombs, basis);
		}


	}






	/////////////////////////////////// algorithms ///////////////////////////////////


	#include "super_bp.inl"


	#include "my_refine.inl"



	//////////////////////////////////////////////////////////////////////////////////







};


#endif

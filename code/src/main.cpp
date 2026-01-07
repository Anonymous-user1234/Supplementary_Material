


/////////////////////////////////////////////////////////////////


#include "linearlayer.hpp"


/////////////////////////////////////////////////////////////////




int main(int argc, char**argv) {


	// parameters
	assert(argc >= 3);
	std::string sbox_file = argv[1];
	std::string ops_file = argv[2];

	std::random_device rd;
	int seed = rd();

	std::cout << "seed = " << seed << std::endl;
	e.seed(seed);


	print_paramters();


	// generate initial mat

	// setting basic values
	wind_t best_gxor = (DIMENSION + NL_NUM) * (DIMENSION + NL_NUM);
	bp_param_t param;

	param.display = false; // show details


	// read matrix from file. update table only.
	Matrix init_mat(sbox_file);

	// load existing ops
	init_mat.ops_read(ops_file);



	// analyze nonlinear inputs.
	init_mat.set_targets();


	
	// setting inital basis.
	init_mat.set_init_basis();
	init_mat.ops.clear();


	init_mat.basis_print();

	for (auto const & T : init_mat.targets) 
		T.print();

	init_mat.gys_print();
	
	
	init_mat.set_target_bounds();


	PRINTF_STAMP("starts rounds\n");


	wind_t i = 0;
	int new_it = 0;
	int last_best= best_gxor;
	while (true) {	

		Matrix mat;
		mat = init_mat;

		if (!mat.universal_bp(param, best_gxor)){
			continue;
		}
		
		if(best_gxor < last_best){
			new_it = i;
			last_best = best_gxor;
		}
		PRINTF_STAMP("round %d done,\t\t\t\t\t\t best Xors: %d(iter:%d), current Xors: %ld\n\n", i++, best_gxor - NL_NUM, new_it, mat.ops.size() - NL_NUM);
		if(i > 30000 && i > new_it*2){
			printf("stable\n");
			break;
		}
	}
	
	if (CHECK_CORRECT) PRINTF_STAMP("All done correctly\n");
	
    return 0;
}















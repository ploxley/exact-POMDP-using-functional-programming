/* LP for findRegionPoint */

/* Three changes (1-3) below are required for switching between concave and convex 
peicewise linear functions (pwlf). Max for convex pwlf, Min for concave pwlf. */

#include "export.h"
#include <glpk.h>
#include "math.h"

int LP_solve(Pointer mat_in, Int32 mat_length, Pointer obj_out, Pointer x_out, Int32 x_length, Pointer stats)
{
	// Set up LP
	glp_prob *lp;
	double* mat = (double*) mat_in;
	double* obj = (double*) obj_out;
	double* x = (double*) x_out;
	int* stat = (int*) stats;
	int num_vars = x_length;  // all arrays padded with start zero to index from 1
	int num_rows = (mat_length-1)/num_vars;
	int ia[mat_length];		
	int ja[mat_length];
	lp = glp_create_prob();  
	glp_set_obj_dir(lp, GLP_MIN); // 1. maximize: GLP_MAX, minimize: GLP_MIN
	
	// Set objective variables
	glp_add_cols(lp, num_vars);	// cols are variables
	glp_add_rows(lp, num_rows);	// rows are constraint eqns
	for (int i = 1; i != num_vars; ++i) {
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
		glp_set_obj_coef(lp, i, 0.0);	// variable coeffs in objective
	}
	glp_set_col_bnds(lp, num_vars, GLP_UP, 0.0, 0.0);  // 2. GLP_LO for max, GLP_UP for min
	glp_set_obj_coef(lp, num_vars, 1.0);	// final variable coeff in objective
	for (int i = 1; i != num_rows; ++i) {
		glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); // 3. GLP_LO for max, GLP_UP for min
	}
	glp_set_row_bnds(lp, num_rows, GLP_FX, 1.0, 1.0);
	
	// Fill out constraint matrix indices (111 222 333, 123 123 123 pattern) 
	for (int i = 1; i != mat_length; ++i) {
		ia[i] = 1 + (i-1) / num_vars;
		ja[i] = num_vars - (num_vars-1)*i % num_vars;
	}
	glp_load_matrix(lp, mat_length-1, ia, ja, mat);
	
	
	// Solve LP and return result
	glp_simplex(lp, NULL); 			// NULL -- solver uses default settings
	*stat = glp_get_status(lp);
	*obj = glp_get_obj_val(lp);
        for (int i = 1; i != num_vars; ++i) {
	        x[i] = glp_get_col_prim(lp, i);
	}
	
	glp_delete_prob(lp);
	
	return 0;
}

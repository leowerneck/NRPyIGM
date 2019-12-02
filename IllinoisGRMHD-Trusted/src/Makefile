CXX=g++

#Default compiler flags:
CXXFLAGS=-O2 -Wall -fopenmp
DEFINES=-DENABLE_STANDALONE_IGM_C2P_SOLVER

default: driver_conserv_to_prims

driver_conserv_to_prims: apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C driver_conserv_to_prims.C font_fix_gamma_law.C  \
			 harm_primitives_lowlevel.C  harm_utoprim_2d.c inlined_functions.C standalone_conserv_to_prims_main_function.h IllinoisGRMHD_Noble2D_con2prim_onepoint.C \
			 convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij.C eigen.C harm_primitives_headers.h harm_u2p_util.c IllinoisGRMHD_headers.h Makefile \
			 NRPy+ADMgammaUU.h  NRPy+g4DD.h  NRPy+g4UU.h
	$(CXX)  $(CXXFLAGS) $(DEFINES) driver_conserv_to_prims.C -o driver_conserv_to_prims

clean: 
	rm -f driver_conserv_to_prims

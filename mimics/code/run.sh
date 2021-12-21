#!/bin/bash
set -ex

gfortran   -o mimics1.5  ../data/trunk_hght_f.f ./src/branch_phase.f ./src/crown_layer.f ./src/initialize.f ./src/matrix_ops.f ./src/pdf_setup.f  ./src/solve_canopy.f \
./src/branch_phase2.f  ./src/cyl_size_int.f ./src/leaf_phase.f ./src/pdf_setup2.f ./src/branch_phase3.f ./src/cylinf_sub.f ./src/leaf_po_smat.f ./src/ndl_rlgh_smat.f \
./src/pdf_sz_setup.f ./src/stokes_sub.f ./src/branch_phase4.f ./src/dielectric.f ./src/leaf_rlgh_smt.f ./src/needle_phase.f ./src/trunk_layer.f ./src/branch_phase5.f \
./src/format_output.f ./src/lngthncyl_smt.f ./src/pd_funcs.f ./src/res_cyl_smat.f ./src/write_data.f ./src/branch_phase6.f ./src/ground_layer.f \
./src/pd_funcs2.f ./src/snow_layer.f ./src/write_error.f ./src/mimics_1_5a.f ./src/start_loop.f ./src/read_input.f ./src/kappa_mat_sub.f ./src/dbl_cmx_fnct.f \
./src/dbl_functions.f ./src/functions.f ./src/naas_subs.f ./src/vector_ops.f

./mimics1.5

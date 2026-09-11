Vendored copy of scripts/mpc/bsolver/{core,model/generated}.
DO NOT EDIT — regenerate with scripts/mpc/bsolver/sync_firmware.sh.
Byte-identical to the source; the include path in core/bs_solver.h resolves
to the sibling model/generated/ in both trees, so no rewriting is needed.

SIX generated headers are here, four configurations:
  bs_model_data.h        + bs_reference_data.h          record,        N = 40
  bs_model_corner_data.h + bs_reference_corner_data.h   corner-online, N = 30
  bs_model_runtime_const.h     runtime constants, model v1 (NU 3, 50 rows)
  bs_model_runtime_const_v2.h  runtime constants, model v2 (NU 2, 48 rows,
                               exact-ZOH e_l; BS_MODEL_V2 defined)
Each pair / runtime header defines the same symbols, so exactly one may be
included in a translation unit.  The selection is made in ../wscript through
-DBS_DATA_HEADER (consumed by core/bs_solver.h) and -DBS_REF_HEADER (consumed
by ../AP_BSolver.cpp), never by editing a #include.  Under the v2 header the
fixed-clock mission builder (core/bs_mission_build.c) compiles to stubs; the
arc-cell clock runtime (core/bs_arc.c) is dimension-generic and serves both.

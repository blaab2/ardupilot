Vendored copy of scripts/mpc/bsolver/{core,model/generated}.
DO NOT EDIT — regenerate with scripts/mpc/bsolver/sync_firmware.sh.
Byte-identical to the source; the include path in core/bs_solver.h resolves
to the sibling model/generated/ in both trees, so no rewriting is needed.

FOUR generated headers are here, two configurations:
  bs_model_data.h        + bs_reference_data.h          record,        N = 40
  bs_model_corner_data.h + bs_reference_corner_data.h   corner-online, N = 30
Each pair defines the same symbols, so exactly one pair may be included in a
translation unit.  The selection is made in ../wscript through
-DBS_DATA_HEADER (consumed by core/bs_solver.h) and -DBS_REF_HEADER (consumed
by ../AP_BSolver.cpp), never by editing a #include.

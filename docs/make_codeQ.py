#!/usr/bin/env python3
# Derive Code_Q_CapeNested.m from Code_P_Nested.m.
# The CAPE-vs-no-CAPE study (PROP figures) still used the ECA and built its
# geometries with array_geometry(..., N_ref, ...), which gives a 64-element
# ULA against a 23-element nested array - an unequal comparison. It now uses
# the NESTED array and builds both arrays at 32 physical antennas, matching
# Part 2.
import io

SRC = "/home/user/csc/Code_P_Nested.m"
DST = "/home/user/csc/Code_Q_CapeNested.m"
with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ---- nested instead of ECA in the CAPE study --------------------------------
rep("d_arr       = {'ula', 'eca'};",
    "d_arr       = {'ula', 'nested'};   % [CODE Q] nested, not ECA\n"
    "d_Nphys     = 32;                  % [CODE Q] equal PHYSICAL antennas, as in Part 2")
rep("d_names     = {'ULA + local est.', 'ECA + local est.', 'ULA + CAPE', 'ECA + CAPE (proposed)'};",
    "d_names     = {'Dense ULA + local est.', 'Nested + local est.', ...\n"
    "               'Dense ULA + CAPE', 'Nested + CAPE (proposed)'};")

# ---- equal-N geometry in both places the study builds arrays ----------------
rep("    D_geo{a} = array_geometry(d_arr{a}, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa, nest_N1, nest_N2);",
    "    % [CODE Q] both arrays at the SAME physical antenna count, otherwise the\n"
    "    % CAPE comparison is confounded by array gain (array_geometry would give\n"
    "    % a 64-element ULA against a 23-element nested array).\n"
    "    D_geo{a} = geom_targetN(d_arr{a}, d_Nphys, lambda);")
rep("    Ga = array_geometry(d_arr{a}, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa, nest_N1, nest_N2);",
    "    Ga = geom_targetN(d_arr{a}, d_Nphys, lambda);   % [CODE Q] equal-N geometry")

rep("save('codeP_nested_results.mat'", "save('codeQ_capenested_results.mat'")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code Q written:", t.count("\n") + 1, "lines")

#!/usr/bin/env python3
# Derive Code_O_Final.m from Code_I_Comparison.m.
#
# ONLY the three changes requested - nothing else is added:
#   1. Fix the per-AP vs fused likelihood figure (PROP-1). It passed XM(:)
#      alone to geo_from_pos, so only one axis of the grid was used and the
#      map rendered as stripes. The PROP study is switched back on so the
#      figure is actually produced.
#   2. PART 1 (dense array) unchanged except the array: DENSE ULA with 64
#      antennas per AP.
#   3. PART 2 (sparse-array comparison) unchanged except the antenna count:
#      32 physical antennas per AP, same two cases as before (one lightly
#      loaded, one overloaded).
import io

SRC = "/home/user/csc/Code_I_Comparison.m"
DST = "/home/user/csc/Code_O_Final.m"
with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ---- 1. PART 1: DENSE array, 64 antennas per AP -----------------------------
rep("array_type = 'eca';",
    "array_type = 'ula';   % [CODE O] PART 1 uses the DENSE array")
rep("N_ref      = 64;     % elements of the dense ULA and of the USA",
    "N_ref      = 64;     % [CODE O] 64 antennas per AP in PART 1 (dense)")

# ---- 2. the likelihood-map figure --------------------------------------------
rep("d_study     = false;   % [CODE I] PROP figures off (superseded)",
    "d_study     = true;    % [CODE O] on, so the fixed PROP-1 map is produced")
rep("    [d3g, sg] = geo_from_pos(APw_map, XM(:), hDiff);",
    "    % [CODE O] FIX: the grid must carry BOTH coordinates. Passing XM(:)\n"
    "    % alone dropped the y-axis, so the map was evaluated along one axis\n"
    "    % only and came out as vertical stripes instead of a 2-D likelihood.\n"
    "    [d3g, sg] = geo_from_pos(APw_map, XM(:) + 1j * YM(:), hDiff);")

# ---- 3. PART 2: 32 antennas per AP, same two cases ---------------------------
rep("en_N0       = [24 8];       % equal PHYSICAL antennas within each case",
    "en_N0       = [32 32];      % [CODE O] 32 PHYSICAL antennas per AP, both cases")
rep("en_K        = [8 16];       % users per case",
    "en_K        = [8 64];       % [CODE O] same two cases: K<N then K>N")

# ---- required so the K>N case does not overflow the per-user buffer ----------
rep("EN_seD  = nan(4, 16, nSchE, nCaseE, nSetE);         % per-user SE at ref SNR (4 detectors)",
    "EN_seD  = nan(4, max(en_K), nSchE, nCaseE, nSetE);  % [CODE O] sized by the largest K")
rep("    s_seD  = nan(4, 16, nSchE, nCaseE);",
    "    s_seD  = nan(4, max(en_K), nSchE, nCaseE);")

rep("save('codeI_comparison_results.mat'", "save('codeO_final_results.mat'")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code O written:", t.count("\n") + 1, "lines")

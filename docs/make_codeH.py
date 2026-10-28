#!/usr/bin/env python3
# Derive Code_H_CoArray_Proposed.m from Code_G_CoArray_EqualN.m.
# Keeps ONE idea only: symmetric nested array + difference co-array + CAPE
# (the near-field-valid, best-justified choice). The three other ideas are
# removed because they rested on the same DOF gain and performed identically.
# The equal-N two-regime comparison stays (ECA+CAPE vs the proposed idea);
# the old unequal-N co_study (CODE F) is disabled.
import io

SRC = "/home/user/csc/Code_G_CoArray_EqualN.m"
DST = "/home/user/csc/Code_H_CoArray_Proposed.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:160], c)
    t = t.replace(old, new)

# 1. disable the old unequal-N CODE F co_study (its 63-vs-23 figures are unfair)
rep("co_study    = true;", "co_study    = false;   % [CODE H] old unequal-N study off")

# 2. EN study comment block -> two schemes
rep(
"%  geometry builds is left to whatever it yields. Five schemes:\n"
"%    1. ECA + CAPE (original)   eca geometry, NO co-array        (co_idea 0)\n"
"%    2. Idea 1                  symmetric nested + diff co-array (co_idea 1)\n"
"%    3. Idea 2                  nested + cross-AP co-array        (co_idea 2)\n"
"%    4. Idea 3                  nested + sum+diff co-array        (co_idea 3)\n"
"%    5. Idea 4                  nested + hybrid-field             (co_idea 4)",
"%  geometry builds is left to whatever it yields. Two schemes only:\n"
"%    1. ECA + CAPE (original)              eca geometry, NO co-array  (co_idea 0)\n"
"%    2. Proposed (symmetric nested array) difference co-array + CAPE (co_idea 1)\n"
"%  The three other ideas (cross-AP, sum+diff, hybrid-field) all rested on the\n"
"%  same difference-co-array DOF gain and performed identically, so the single\n"
"%  best-justified one is kept: a SYMMETRIC nested array keeps the difference\n"
"%  co-array valid in the near field, and CAPE is used unchanged (only the\n"
"%  physical array is replaced; the co-array just feeds CAPE extra candidates).")

# 3. EN scheme definition -> two schemes
rep(
"en_names    = {'ECA + CAPE (original)', 'Idea 1 (sym-nested)', 'Idea 2 (cross-AP)', ...\n"
"               'Idea 3 (sum+diff)', 'Idea 4 (hybrid)'};\n"
"en_geom     = {'eca', 'symnested', 'nested', 'nested', 'nested'};\n"
"en_idea     = [0 1 2 3 4];",
"en_names    = {'ECA + CAPE (original)', 'Proposed: sym-nested co-array + CAPE'};\n"
"en_geom     = {'eca', 'symnested'};\n"
"en_idea     = [0 1];")

# 4. save tag
rep("save('codeG_equalN_results.mat'", "save('codeH_proposed_results.mat'")

# 5. top-of-file note
rep(
"%   * eqN_study (CODE G): FAIR equal-physical-antenna comparison of all five\n"
"%     schemes in two regimes (antenna-rich K<N and user-rich K>N), EN-1..EN-5.",
"%   * eqN_study (CODE H): FAIR equal-physical-antenna comparison of the ONE\n"
"%     proposed scheme (symmetric nested + co-array + CAPE) against the original\n"
"%     ECA + CAPE, in two regimes (antenna-rich K<N and user-rich K>N), EN-1..EN-5.")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code H written:", t.count("\n") + 1, "lines")

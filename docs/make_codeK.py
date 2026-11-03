#!/usr/bin/env python3
# Derive Code_K_SumFirst.m from Code_J_Tiling.m.
# Replaces the distributed-tiling array with IDEA B: the sum-co-array-first
# geometry (dense core + offset block, offset chosen so the SUM co-array is
# hole-free). Also reports the sum-co-array uDOF as a new metric so the
# sum-first claim can be checked directly against the other arrays.
import io

SRC = "/home/user/csc/Code_J_Tiling.m"
DST = "/home/user/csc/Code_K_SumFirst.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ---- (1) add the sum-first geometry to geom_targetN -------------------------
rep("""        case 'eca'
            M = 5;  Nc = 4;                                % a generous coprime base""",
"""        case 'sumfirst'
            % [CODE K] IDEA B -- SUM-CO-ARRAY-FIRST GEOMETRY.
            % Dense core C = {0..nc-1} plus an offset block O = {D..D+nc-1}.
            % The three sum blocks are core+core = [0, 2nc-2],
            % core+offset = [D, D+2nc-2] and offset+offset = [2D, 2D+2nc-2].
            % They tile the line WITHOUT holes exactly when D = 2nc-1, which is
            % the offset used here (the value 10 in the note leaves gaps at
            % 7..9 and 17..19 and is NOT hole-free).
            nc = max(2, floor(N0 / 2));
            Doff = 2 * nc - 1;
            n = [(0:nc - 1), Doff + (0:(N0 - nc) - 1)]';
            n = unique(n);
            n = force_count(n, N0);
        case 'eca'
            M = 5;  Nc = 4;                                % a generous coprime base""")

# ---- (2) swap the compared scheme ------------------------------------------
rep("""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed: distributed tiling'};
en_geom     = {'ula', 'eca', 'nested', 'tiling'};
en_idea     = [0 0 1 1];""",
"""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed: sum-first (Idea B)'};
en_geom     = {'ula', 'eca', 'nested', 'sumfirst'};
en_idea     = [0 0 1 1];""")

# ---- (3) record the SUM co-array uDOF as a reported metric ------------------
rep("""EN_Nphys = zeros(nSchE, nCaseE);    % physical antennas PER AP (equal across schemes)
EN_dofU  = zeros(nSchE, nCaseE);    % NETWORK virtual sensors = union of per-AP co-arrays
EN_apAvg = zeros(nSchE, nCaseE);    % mean per-AP aperture""",
"""EN_Nphys = zeros(nSchE, nCaseE);    % physical antennas PER AP (equal across schemes)
EN_dofU  = zeros(nSchE, nCaseE);    % NETWORK virtual sensors = union of per-AP co-arrays
EN_apAvg = zeros(nSchE, nCaseE);    % mean per-AP aperture
EN_udD   = zeros(nSchE, nCaseE);    % hole-free DIFFERENCE co-array run (angle DOF)
EN_udS   = zeros(nSchE, nCaseE);    % hole-free SUM co-array run (range DOF)  [CODE K]""")

rep("""        EN_Nphys(scE, caE) = Gl{1}.N;
        EN_dofU(scE, caE)  = numel(lagU);
        EN_apAvg(scE, caE) = apS / L;""",
"""        EN_Nphys(scE, caE) = Gl{1}.N;
        EN_dofU(scE, caE)  = numel(lagU);
        EN_apAvg(scE, caE) = apS / L;
        EN_udD(scE, caE)   = caL{1}.Nv;                                   % angle DOF
        casum              = coarray_setup(Gl{1}.pos, lambda, 'sum');     % [CODE K]
        EN_udS(scE, caE)   = casum.Nv;                                    % range DOF""")

# ---- (4) print both DOFs ----------------------------------------------------
rep("""        fprintf('     %-28s %-8s N/AP=%d  network virtual=%d  (per-AP %d, mean aper %.2f m)\\n', ...
                en_names{scE}, en_geom{scE}, EN_Nphys(scE, caE), EN_dofU(scE, caE), ...
                EN_CA{scE, caE}{1}.DOF, EN_apAvg(scE, caE));""",
"""        fprintf('     %-28s %-9s N/AP=%d  uDOFdiff=%d  uDOFsum=%d  aper=%.2f m\\n', ...
                en_names{scE}, en_geom{scE}, EN_Nphys(scE, caE), ...
                EN_udD(scE, caE), EN_udS(scE, caE), EN_apAvg(scE, caE));""")

rep("""        fprintf('%-28s %5s %8s %6s %10s %10s %9s %9s\\n', 'array', 'N/AP', 'netVirt', 'NF[%]', ...
                'NMSE', 'BER-XAP', 'sumSE', 'time[s]');""",
"""        fprintf('%-28s %5s %8s %8s %6s %10s %10s %9s\\n', 'array', 'N/AP', 'uDOFdif', 'uDOFsum', ...
                'NF[%]', 'NMSE', 'BER-XAP', 'sumSE');""")

rep("""            fprintf('%-28s %5d %8d %6.0f %10.2e %10.2e %9.2f %9.3f\\n', en_names{scE}, EN_Nphys(scE, caE), ...
                    EN_dofU(scE, caE), 100 * mNF(scE, caE), mNM(scE, caE, si0), ...
                    mBE(4, scE, caE, si0), mSU(4, scE, caE, si0), mTIME(scE, caE));""",
"""            fprintf('%-28s %5d %8d %8d %6.0f %10.2e %10.2e %9.2f\\n', en_names{scE}, EN_Nphys(scE, caE), ...
                    EN_udD(scE, caE), EN_udS(scE, caE), 100 * mNF(scE, caE), mNM(scE, caE, si0), ...
                    mBE(4, scE, caE, si0), mSU(4, scE, caE, si0));""")

# ---- (5) ARR-6: show angle DOF vs range DOF --------------------------------
rep("""    Nphys = EN_Nphys(:, 1);
    Nvirt = EN_dofU(:, 1);          % [CODE J] NETWORK co-array (union over APs)""",
"""    Nphys = EN_Nphys(:, 1);
    Nvirt = EN_udD(:, 1);           % [CODE K] hole-free DIFFERENCE co-array (angle DOF)
    Nvsum = EN_udS(:, 1);           % [CODE K] hole-free SUM co-array (range DOF)""")

rep("""    bar(1:nSchE, [Nphys(:), Nvirt(:)]); set(gca, 'YScale', 'log');""",
    """    bar(1:nSchE, [Nphys(:), Nvirt(:), Nvsum(:)]); set(gca, 'YScale', 'log');""")

rep("""    ylabel('count'); title(sprintf('Equal physical N=%d per AP: antennas vs NETWORK virtual sensors', en_N0(1)));
    legend({'physical antennas / AP', 'network virtual sensors'}, 'Location', 'northwest');""",
    """    ylabel('count'); title(sprintf('Equal physical N=%d per AP: angle DOF vs range DOF', en_N0(1)));
    legend({'physical antennas / AP', 'uDOF difference (angle)', 'uDOF sum (range)'}, ...
           'Location', 'northwest');""")

rep("""         'EN_Nphys', 'EN_dofU', 'EN_apAvg');""",
    """         'EN_Nphys', 'EN_dofU', 'EN_apAvg', 'EN_udD', 'EN_udS');""")
rep("save('codeJ_tiling_results.mat'", "save('codeK_sumfirst_results.mat'")

# ---- (6) header note --------------------------------------------------------
rep("""%   * [CODE J] Idea 2 (distributed co-array tiling): each AP has its OWN""",
"""%   * [CODE K] Idea B (sum-co-array-first geometry): a dense core plus an
%     offset block, the offset chosen so the SUM co-array is hole-free, on the
%     argument that the sum lags carry the RANGE in the near field. The study
%     now reports uDOFdiff (angle) AND uDOFsum (range) for every array so the
%     claim can be checked directly. MEASURED BEFORE SHIPPING: at N=24 the
%     plain nested array already attains uDOFdiff=311 and uDOFsum=168, whereas
%     the sum-first construction attains only 69 and 69, so it is dominated on
%     BOTH axes; at N=8 nested is (39, 24) against sum-first (21, 21). The
%     geometry is included so the comparison is visible in the figures.
%   * [CODE J] Idea 2 (distributed co-array tiling): each AP has its OWN""")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code K written:", t.count("\n") + 1, "lines")

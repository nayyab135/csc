#!/usr/bin/env python3
"""Build Code_U_v2.m from Code_U_MovableSparse.m.
 * many more drops (smooth curves), Code T run size per drop (2 x 400 symbols)
 * faster movable-array estimator (one dictionary start, 6 refinement rounds)
 * CAPE: final full-window refinement from the chosen lobe
 * no-CAPE localisation baseline (AoA triangulation, same APs / lobes as CAPE)
 * fraction of users located within 1 m (EN_loc1)
 * figures V1-V9 (old U4-U11, U13 behind fig_old)
Every edit is anchored and must match exactly once."""
import pathlib
import sys

REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent
src = (REPO / 'Code_U_MovableSparse.m').read_text()
figu = (HERE / 'figblock_u.m').read_text()
figv = (HERE / 'figs_v2.m').read_text()
codeT = (REPO / 'Code_T_ArrayComparison.m').read_text()


def rep(a, b):
    global src
    n = src.count(a)
    if n != 1:
        sys.exit(f'anchor found {n} times:\n{a[:300]}')
    src = src.replace(a, b)


ind = lambda x: ''.join(('    ' + ln if ln.strip() else ln) for ln in x.splitlines(True))

# ------------------------------------------------------------------ header
rep("%  CODE U : MOVABLE (SLIDING-MIRROR) SPARSE ARRAY vs STATIC SPARSE ARRAYS",
    "%  CODE U v2 : MOVABLE (SLIDING-MIRROR) SPARSE ARRAY vs STATIC SPARSE ARRAYS")
rep("""%  fast_mode = true  -> quick run;  fast_mode = false -> full paper run.
%  Figures U1..U13. Results saved to codeU_movable_results.mat BEFORE plotting;
%  plot_codeU_figures.m redraws the figures without rerunning.""",
"""%  v2: many more random drops (smooth curves; the BER noise comes from the
%  drops, not from the symbols), the run size per drop of Code T; a faster
%  movable-array estimator (its stacked lambda/4 aperture has no grating
%  lobes, so one dictionary start + the two mirror-pair angles suffice, and 6
%  refinement rounds are 200x finer than the hit tolerance); CAPE ends with a
%  full-window refinement from its chosen lobe; a no-CAPE localisation
%  baseline (angle-of-arrival triangulation from the same APs and lobes); the
%  share of users located within 1 m. Figures V1-V9 (U1-U3, U12 geometry;
%  the older U4-U11, U13 with fig_old = true).
%  fast_mode = true  -> quick check;  fast_mode = false -> full run.
%  Results saved to codeU_v2_results.mat BEFORE plotting;
%  plot_codeU_v2_figures.m redraws the figures without rerunning.""")

# ---------------------------------------------------------------- settings
rep("""fast_mode = true;     % true: quick run;  false: full paper run (see the run-time print)
if fast_mode
    nSetE = 6;   nRealE = 2;  nSymE = 300;  Ksw = [8 16 24];
else
    nSetE = 10;  nRealE = 4;  nSymE = 800;  Ksw = [8 12 16 20 24];   % = Code T
end""",
"""fast_mode = false;    % true: quick check;  false: full run (the script prints its run-time estimate)
if fast_mode
    nSetE = 28;   nSetK = 10;
else
    nSetE = 140;  nSetK = 50;     % random drops: SNR study / user-count sweep (about the cost of Code T)
end
nRealE = 2;  nSymE = 400;         % per drop, as Code T (the BER noise comes from the drops)
Ksw = [8 12 16 20 24];            % user-count sweep (schemes with CAPE)""")

# ------------------------------------------------------ movable-array speed-up
rep("""        PDe.smax = pd_smax;  PDe.stop = pd_stop;  PDe.nref = pd_nref;  PDe.nstart = pd_nstart;
        ES.PDe = PDe;  ES.Ne = numel(G.frames);  ES.mir = G.mir;""",
"""        PDe.smax = pd_smax;  PDe.stop = pd_stop;  PDe.nref = 6;  PDe.nstart = pd_nstart;   % [v2] 6 rounds
        ES.PDe = PDe;  ES.Ne = numel(G.frames);  ES.mir = G.mir;
        ES.J = 1;          % [v2] stacked lambda/4 aperture: no grating lobes, one dictionary start""")
rep("    ES = struct('kind', bs_kind{b}, 'F', nPilotBlk, 'PDe', PDb, 'Ne', G.N, 'mir', [], 'CA', []);",
    "    ES = struct('kind', bs_kind{b}, 'F', nPilotBlk, 'PDe', PDb, 'Ne', G.N, 'mir', [], 'CA', [], 'J', cape_J);")
rep("""    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    if strcmp(ES.kind, 'sma')
        ES.he""",
"""    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    FG.J = ES.J;                                   % [v2] dictionary starts per link
    if strcmp(ES.kind, 'sma')
        ES.he""")

# --------------------------------------------------------------- main study
rep("""fprintf(['  estimated run time with %d worker(s): main study about %.0f min, user sweep about %.0f min ' ...
         '(fast_mode = %d)\\n'], nW, 2 * sum(tB) * nSNR * nSetE / nW / 60, ...
        sum(tB) * nSetE * numel(Ksw) / nW / 60, fast_mode);""",
"""fprintf(['  estimated run time with %d worker(s): main study about %.0f min, user sweep about %.0f min ' ...
         '(fast_mode = %d)\\n'], nW, 2 * sum(tB) * nSNR * nSetE / nW / 60, ...
        sum(tB) * nSetK * numel(Ksw) / nW / 60, fast_mode);""")
rep("T_nf   = zeros(1, nTask);  T_seD = nan(4, Ku, nTask);",
    "T_nf   = zeros(1, nTask);  T_seD = nan(4, Ku, nTask);  T_loc1 = nan(1, nTask);")
rep("    if ~isempty(O.loc), T_loc(tk) = median(O.loc); end\n    T_hit(tk) = O.hit;",
    "    if ~isempty(O.loc), T_loc(tk) = median(O.loc);  T_loc1(tk) = mean(O.loc < 1); end\n    T_hit(tk) = O.hit;")
rep("EN_rmse = reshape(permute(reshape(T_rm, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);",
    "EN_rmse = reshape(permute(reshape(T_rm, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);\n"
    "EN_loc1 = reshape(permute(reshape(T_loc1, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);   % users within 1 m")

# ---------------------------------------------------------------- K sweep
a = "% ---- [CODE U] user-count sweep at the reference SNR (schemes with CAPE) ----"
b = "if eqN_study\n    % save BEFORE plotting"
i, j = src.index(a), src.index(b)
src = src[:i] + src[i:j].replace('nSetE', 'nSetK') + src[j:]

# ------------------------------------------------------------ save + figures
rep("""    save('codeU_movable_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', 'Dtgt', 'tolS', ...
         'nPilotBlk', 'sma_sig_d', 'bs_name', 'bs_geom', 'bs_kind', 'en_base', 'en_cape', 'en_names', ...
         'en_caseName', 'en_N0', 'en_K', 'fast_mode', 'nSetE', 'nRealE', 'nSymE', 'modOrder', 'clu_div', ...
         'EN_geo', 'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', ...
         'EN_nf', 'EN_seD', 'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\\nResults saved to codeU_movable_results.mat\\n');""",
"""    save('codeU_v2_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', 'Dtgt', 'tolS', ...
         'nPilotBlk', 'sma_sig_d', 'bs_name', 'bs_geom', 'bs_kind', 'en_base', 'en_cape', 'en_names', ...
         'en_caseName', 'en_N0', 'en_K', 'fast_mode', 'nSetE', 'nSetK', 'nRealE', 'nSymE', 'modOrder', 'clu_div', ...
         'EN_geo', 'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_loc1', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', ...
         'EN_nf', 'EN_seD', 'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\\nResults saved to codeU_v2_results.mat\\n');""")

def seg(t, a, b=None):
    i = t.index(a); j = t.index(b) if b else len(t); return t[i:j]
pre  = figu[:figu.index("% ---- U1:")].replace(
    "%% ===================== [CODE U] FIGURES U1-U13 =====================",
    "%% ===================== [CODE U v2] FIGURES =====================\n"
    "%  V1-V9: smooth comparison figures (see below). U1-U3, U12: geometry.\n"
    "%  The older figures U4-U11, U13: set fig_old = true.")
u1_3 = seg(figu, "% ---- U1:", "% ---- U4-U6:")
u4_11 = seg(figu, "% ---- U4-U6:", "% ---- U12:")
u12 = seg(figu, "% ---- U12:", "% ---- U13:")
u13 = seg(figu, "% ---- U13:")
fig = (pre + "if ~exist('fig_old', 'var'), fig_old = false; end   % true: also draw U4-U11, U13\n"
       + u1_3 + u12 + "if fig_old\n" + ind(u4_11 + u13) + "end\n\n" + figv)
i = src.index("    %% ===================== [CODE U] FIGURES U1-U13 =====================")
j = src.index("function [Hhat, C, Hpre, Cpre] = csi_acquire(")
k = src.rindex("end\n", i, j)
src = src[:i] + ind(fig) + src[k:]

# --------------------------------------------------- estimator: CAPE + lobes
rep("function [Hhat, C, xhat, nfus, aErr] = csi_nf_estimate(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG, EST, useCape)",
    "function [Hhat, C, xhat, nfus, aErr, cSa, cQa] = csi_nf_estimate(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG, EST, useCape)")
rep("""    gamL = log(size(PDe.W, 2)) + 0.58;
    nc = FG.J + EST.nExtra;""",
"""    gamL = log(size(PDe.W, 2)) + 0.58;
    nc = FG.J + EST.nExtra;
    cSa = nan(L, K, nReal, nc);  cQa = nan(L, K, nReal, nc);   % [v2] local candidate lobes (triangulation baseline)""")
rep("""            if useCape
                cand = [];""",
"""            cSa(:, k, m, :) = reshape(cS, L, 1, 1, nc);  cQa(:, k, m, :) = reshape(cQ, L, 1, 1, nc);
            if useCape
                cand = [];""")
rep("""                    if ~((sameLobe && mLoc(l) > mfu) || clearWin)
                        b = bf;  s1 = sf;  q1 = qf;
                    end
                end""",
"""                    if ~((sameLobe && mLoc(l) > mfu) || clearWin)
                        b = bf;  s1 = sf;  q1 = qf;
                    end
                    [b2, m2, s2, q2] = refine_atom(zl, s1, q1, PDe);   % [v2] full window from the chosen lobe
                    if m2 > abs(b' * zl)
                        b = b2;  s1 = s2;  q1 = q2;
                    end
                end""")

# ---------------------------------------------- eval_scheme_u: local positions
rep("""    [Hh, Cs, xh, nfus, aErr] = csi_nf_estimate(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG, EST, useCape);
    O.time = toc(tE);
    O.nfus = mean(nfus);
    O.loc = [];
    if useCape
        dx""",
"""    [Hh, Cs, xh, nfus, aErr, cSa, cQa] = csi_nf_estimate(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG, EST, useCape);
    O.time = toc(tE);
    O.nfus = mean(nfus);
    O.loc = [];
    if ~useCape                                       % [v2] no-CAPE localisation baseline: at most 4
        nU = min(size(cSa, 4), 4);                    % lobes per AP (dictionary lobes + strongest extra),
        xh = local_triangulate(cSa(:, :, :, 1:nU), cQa(:, :, :, 1:nU), ST.NF, bet, FG);   % as in Code T
    end
    if true
        dx""")

# ------------------------------------------ triangulation (from Code T v3)
i = codeT.index("function xh = local_triangulate(sE, qE, NFm, beta, FG)")
src = src.rstrip() + "\n\n" + codeT[i:].rstrip() + "\n"
(REPO / 'Code_U_v2.m').write_text(src)

plot = (r"""%% PLOT_CODEU_V2_FIGURES  Redraw the Code U v2 figures without rerunning the simulation.
%  V1-V9 (smooth comparison figures), U1-U3 and U12 (geometry). Set
%  fig_old = true before running this script to also draw U4-U11 and U13.
%  Uses the workspace variables after Code_U_v2.m, or loads codeU_v2_results.mat
%  if they are not there. It never clears anything. No local functions.
if ~exist('EN_be', 'var') || ~exist('EN_loc1', 'var') || ~exist('bs_kind', 'var')
    load('codeU_v2_results.mat');
end
""" + fig)
(REPO / 'plot_codeU_v2_figures.m').write_text(plot)
print('Code_U_v2 lines:', src.count('\n'), ' plot lines:', plot.count('\n'))

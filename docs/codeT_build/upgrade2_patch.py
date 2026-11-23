#!/usr/bin/env python3
"""Second upgrade of Code_T_ArrayComparison.m (in place), on top of upgrade_patch.py:
 * no-CAPE localisation baseline: hard-decision triangulation of the local
   per-AP (angle, range) estimates over the same fused APs
 * CAPE final step: full local refinement (polar OMP) started from the CAPE
   lobe; the estimate with the lower data-estimated error is kept
 * fraction of users located within 1 m per task (EN_loc1)
 * new figure groups A (arrays), B (co-array), C (CAPE); old T3-T17 behind fig_old
Every edit is anchored and must match exactly once."""
import pathlib
import sys

REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent
f = REPO / 'Code_T_ArrayComparison.m'
src = f.read_text()
fig = (HERE / 'figblock_t3.m').read_text()


def rep(a, b):
    global src
    n = src.count(a)
    if n != 1:
        sys.exit(f'anchor found {n} times:\n{a[:300]}')
    src = src.replace(a, b)


def cut(start, end, new):
    global src
    for a in (start, end):
        if src.count(a) != 1:
            sys.exit(f'cut anchor found {src.count(a)} times:\n{a[:300]}')
    i, j = src.index(start), src.index(end)
    assert j > i
    src = src[:i] + new + src[j:]


# ------------------------------------------------------------------ header
rep("""%  CO-ARRAY METRICS: (a) true-lobe angle hit rate / RMSE of the channel
%  estimator; (b) a classic DOA test of the same antennas: physical MUSIC vs
%  co-array MUSIC, sources found vs number of sources, RMSE vs SNR.""",
"""%  CO-ARRAY METRICS: (a) true-lobe angle hit rate / RMSE of the channel
%  estimator; (b) a classic DOA test of the same antennas: physical MUSIC vs
%  co-array MUSIC, sources found vs number of sources, RMSE vs SNR.
%  LOCALISATION BASELINE WITHOUT CAPE: every AP's own (angle, range) estimate
%  gives candidate positions; the CPU picks the most consistent ones over the
%  same fused APs (hard-decision triangulation). CAPE fuses the pilot
%  LIKELIHOODS instead, then refines every link locally from the fused lobe.
%  FIGURES: A1-A4 array comparison, B1-B3 role of the difference co-array,
%  C1-C3 benefit of CAPE (T1, T2, T9 geometry; T3-T17 with fig_old = true).""")
rep("%  Figures T1..T17. Results are saved to codeT_arraycomp_results.mat BEFORE",
    "%  Results are saved to codeT_arraycomp_results.mat BEFORE")

# --------------------------------------------------------- main study I/O
rep("T_seD  = nan(4, max(en_K), nTask);  T_hit = nan(1, nTask);  T_rm = nan(1, nTask);",
    "T_seD  = nan(4, max(en_K), nTask);  T_hit = nan(1, nTask);  T_rm = nan(1, nTask);  T_loc1 = nan(1, nTask);")
rep("    if ~isempty(O.loc), T_loc(tk) = median(O.loc); end",
    "    if ~isempty(O.loc), T_loc(tk) = median(O.loc);  T_loc1(tk) = mean(O.loc < 1); end")
rep("EN_rmse = permute(reshape(T_rm, sz), [2 3 1 4]);",
    "EN_rmse = permute(reshape(T_rm, sz), [2 3 1 4]);\n"
    "EN_loc1 = permute(reshape(T_loc1, sz), [2 3 1 4]);   % fraction of users located within 1 m")
rep("         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...",
    "         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_loc1', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...")

# ----------------------------------------------------------- figure block
i = src.index("    %% ===================== [CODE T] FIGURES T1-T17 =====================")
j = src.index("function [Hhat, C, Hpre, Cpre, aErr] = csi_acquire(")
k = src.rindex("end\n", i, j)              # the 'end' that closes 'if eqN_study'
fig_ind = ''.join(('    ' + ln if ln.strip() else ln) for ln in fig.splitlines(True))
src = src[:i] + fig_ind + src[k:]

# ------------------------------------------------- csi_acquire: (s, q) out
rep("function [Hhat, C, Hpre, Cpre, aErr] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...",
    "function [Hhat, C, Hpre, Cpre, aErr, sEst, qEst] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...")
rep("""    aErr = nan(L, K, nReal);                   % [CODE T] angle error of the estimate (scoring only)
    hasFG = nargin >= 17 && ~isempty(FG);""",
"""    aErr = nan(L, K, nReal);                   % [CODE T] angle error of the estimate (scoring only)
    sEst = nan(L, K, nReal);  qEst = nan(L, K, nReal);   % local (angle, inverse range) of the first atom
    hasFG = nargin >= 17 && ~isempty(FG);""")
rep("""                    [hh, se2, sA] = pdomp_estimate(z, PD, beta(l, k), sw2, caAng{l});
                    Hhat(idx, :, k) = hh;
                    C(:, :, l, k) = mean(se2) * eye(N);""",
"""                    [hh, se2, sA, qA] = pdomp_estimate(z, PD, beta(l, k), sw2, caAng{l});
                    Hhat(idx, :, k) = hh;
                    C(:, :, l, k) = mean(se2) * eye(N);
                    sEst(l, k, :) = reshape(sA, 1, 1, []);  qEst(l, k, :) = reshape(qA, 1, 1, []);""")

# ------------------------------------------------- pdomp: q out + forced seed
rep("function [hh, se2, sAng] = pdomp_estimate(z, PD, beta_lk, sw2, sx)",
    "function [hh, se2, sAng, qAng] = pdomp_estimate(z, PD, beta_lk, sw2, sx, sq0)")
rep("""% [CODE T] sx (optional): extra start angles for the FIRST atom (co-array
% MUSIC angles, range by a 1-D search). sAng: angle of the first atom.
    if nargin < 5, sx = []; end""",
"""% [CODE T] sx (optional): extra start angles for the FIRST atom (co-array
% MUSIC angles, range by a 1-D search). sq0 (optional, [s q]): the FIRST atom
% is refined from this point ONLY (used by CAPE to keep its lobe). sAng/qAng:
% angle and inverse range of the first atom.
    if nargin < 5, sx = []; end
    if nargin < 6, sq0 = []; end""")
rep("""    sAng = nan(1, M);
    gam = log(size(PD.W, 2)) + 0.58 + 3;""",
"""    sAng = nan(1, M);  qAng = nan(1, M);
    gam = log(size(PD.W, 2)) + 0.58 + 3;""")
rep("""            bBest = [];
            mBest = -inf;
            sBest = NaN;
            for st = 1:min(PD.nstart, numel(ordc))
                j = ordc(st);
                [b, mg, sb] = refine_atom(r, PD.s(j), PD.q(j), PD);
                if mg > mBest
                    mBest = mg;
                    bBest = b;
                    sBest = sb;
                end
            end
            if it == 1
                for jj = 1:numel(sx)                       % [CODE T] co-array seeds
                    q0 = coarray_range_search(r, sx(jj), PD, PD.nRng);
                    [b, mg, sb] = refine_atom(r, sx(jj), q0, PD);
                    if mg > mBest
                        mBest = mg;
                        bBest = b;
                        sBest = sb;
                    end
                end
                sAng(m) = sBest;
            end""",
"""            bBest = [];
            mBest = -inf;
            sBest = NaN;  qBest = NaN;
            if it == 1 && ~isempty(sq0)                    % [CODE T] forced start (CAPE lobe)
                [bBest, mBest, sBest, qBest] = refine_atom(r, sq0(1), sq0(2), PD);
            else
                for st = 1:min(PD.nstart, numel(ordc))
                    j = ordc(st);
                    [b, mg, sb, qb] = refine_atom(r, PD.s(j), PD.q(j), PD);
                    if mg > mBest
                        mBest = mg;
                        bBest = b;
                        sBest = sb;  qBest = qb;
                    end
                end
            end
            if it == 1
                for jj = 1:numel(sx)                       % [CODE T] co-array seeds
                    q0 = coarray_range_search(r, sx(jj), PD, PD.nRng);
                    [b, mg, sb, qb] = refine_atom(r, sx(jj), q0, PD);
                    if mg > mBest
                        mBest = mg;
                        bBest = b;
                        sBest = sb;  qBest = qb;
                    end
                end
                sAng(m) = sBest;  qAng(m) = qBest;
            end""")

# ------------------------------------------- CAPE final step: local refinement
rep("""                [b, mf, sSel] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);  % fused lobe, local precision
                sameLobe = abs(bLoc(:, l)' * b) >= FG.consist;
                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level
                if (sameLobe && mLoc(l) > mf) || clearWin                  % consistency / likelihood gate
                    b = bLoc(:, l);
                    sSel = sLoc(l, jLoc(l));
                end
                if hasTrue
                    dS = abs(sSel - FG.sTrue(l, k));
                    aErr(l, k, m) = min(dS, 2 - dS);
                end
                c = b' * zl;
                g = N * beta(l, k) / (N * beta(l, k) + sw2);
                Hhat(idx, m, k) = b * (g * c);
                res = zl - b * c;
                off = max(real(res' * res) - (N - 1) * sw2, 0);
                se2(l, m) = (N * beta(l, k) * sw2 / (N * beta(l, k) + sw2) + off) / N;""",
"""                [b, mf, sSel, qSel] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);   % fused lobe, local precision
                sameLobe = abs(bLoc(:, l)' * b) >= FG.consist;
                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level
                if (sameLobe && mLoc(l) > mf) || clearWin                  % consistency / likelihood gate
                    b = bLoc(:, l);
                    sSel = sLoc(l, jLoc(l));  qSel = qLoc(l, jLoc(l));
                end
                c = b' * zl;
                g = N * beta(l, k) / (N * beta(l, k) + sw2);
                hC = b * (g * c);
                res = zl - b * c;
                off = max(real(res' * res) - (N - 1) * sw2, 0);
                seC = (N * beta(l, k) * sw2 / (N * beta(l, k) + sw2) + off) / N;
                % [CODE T] full local refinement STARTED FROM THE CAPE LOBE (polar
                % OMP: wide window, a second path if significant, data-driven
                % gains). The lobe chosen by CAPE is kept; of the two estimates the
                % one with the lower data-estimated error is used, so CAPE never
                % discards a better local fit (high SNR) and keeps its fused atom
                % on weak links (low SNR).
                [hP, seP, sP] = pdomp_estimate(zl, PD, beta(l, k), sw2, [], [sSel, qSel]);
                if seP < seC
                    Hhat(idx, m, k) = hP;  se2(l, m) = seP;  sSel = sP;
                else
                    Hhat(idx, m, k) = hC;  se2(l, m) = seC;
                end
                if hasTrue
                    dS = abs(sSel - FG.sTrue(l, k));
                    aErr(l, k, m) = min(dS, 2 - dS);
                end""")

# ----------------------------------------------- eval_scheme: local positions
rep("""        FGl = FG;  FGl.on = false;              % [CODE T] local estimation (+ co-array seeds)
        [Hh, Cs, ~, ~, aErr] = csi_acquire(Ht, ST.R2, ST.NF, ST.hD, bet, (1:K)', K, p, ...
                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD, FGl);""",
"""        FGl = FG;  FGl.on = false;              % [CODE T] local estimation (+ co-array seeds)
        [Hh, Cs, ~, ~, aErr, sEl, qEl] = csi_acquire(Ht, ST.R2, ST.NF, ST.hD, bet, (1:K)', K, p, ...
                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD, FGl);
        xh = local_triangulate(sEl, qEl, ST.NF, bet, FG);   % [CODE T] no-CAPE localisation baseline
        dx = abs(real(xh) - real(UEm(:)));  dx = min(dx, FG.sL - dx);
        dy = abs(imag(xh) - imag(UEm(:)));  dy = min(dy, FG.sL - dy);
        e = sqrt(dx.^2 + dy.^2);
        O.loc = e(isfinite(e));""")

# -------------------------------------------------------------- new function
NEWFUN = r"""

function xh = local_triangulate(sE, qE, NFm, beta, FG)
% [CODE T] Localisation WITHOUT CAPE (baseline): every one of the FG.nAP
% strongest near-field APs of a user makes its own hard decision (its local
% angle and inverse range), which gives candidate positions (two mirror
% solutions; a ray of ranges for a far-field atom). The CPU picks, among all
% candidates, the one most consistent with the other APs (smallest sum of
% distances to their nearest candidates) and averages it with those nearest
% candidates. Same APs and pilots as CAPE; only decisions are fused, not the
% pilot likelihoods.
    [~, K, M] = size(sE);
    xh = nan(K, M);
    sL = FG.sL;
    td = @(a, b) sqrt(min(abs(real(a) - real(b)), sL - abs(real(a) - real(b))) .^ 2 + ...
                      min(abs(imag(a) - imag(b)), sL - abs(imag(a) - imag(b))) .^ 2);
    for k = 1:K
        nf = find(NFm(:, k));
        if isempty(nf), continue; end
        [~, ob] = sort(beta(nf, k), 'descend');
        A = nf(ob(1:min(FG.nAP, numel(nf))));
        for m = 1:M
            P = cell(numel(A), 1);  ok = true(numel(A), 1);
            for i = 1:numel(A)
                if ~isfinite(sE(A(i), k, m)), ok(i) = false; continue; end
                P{i} = cand_positions(FG.AP(A(i)), sE(A(i), k, m), qE(A(i), k, m), sL, FG.hD);
            end
            P = P(ok);
            if isempty(P), continue; end
            if numel(P) == 1, xh(k, m) = P{1}(1); continue; end
            allC = vertcat(P{:});  best = inf;  xb = allC(1);
            for c = allC(:).'
                sc = 0;  nb = c;
                for i = 1:numel(P)
                    [dmin, jm] = min(td(c, P{i}));
                    sc = sc + dmin;  nb = [nb; P{i}(jm)]; %#ok<AGROW>
                end
                if sc < best, best = sc;  xb = nb; end
            end
            off = xb - xb(1);                              % average on the torus around xb(1)
            ox = mod(real(off) + sL / 2, sL) - sL / 2;  oy = mod(imag(off) + sL / 2, sL) - sL / 2;
            x = xb(1) + mean(ox) + 1j * mean(oy);
            xh(k, m) = mod(real(x), sL) + 1j * mod(imag(x), sL);
        end
    end
end
"""
src = src.rstrip() + '\n' + NEWFUN
f.write_text(src)

plot = (r"""%% PLOT_CODET_FIGURES  Redraw the Code T figures without rerunning the simulation.
%  A1-A4: which sparse array performs better; B1-B3: role of the difference
%  co-array; C1-C3: benefit of CAPE; T1, T2, T9: geometry. Set fig_old = true
%  before running this script to also draw the older figures T3-T17.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_loc1', 'var') || ~exist('DOA_fC', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
""" + fig)
(REPO / 'plot_codeT_figures.m').write_text(plot)
print('Code T lines:', src.count('\n'), ' plot lines:', plot.count('\n'))

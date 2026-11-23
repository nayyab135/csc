#!/usr/bin/env python3
"""Upgrade Code_T_ArrayComparison.m (in place):
 * real fast_mode switch and many more drops in the Code T study
 * six schemes WITHOUT CAPE next to the six with CAPE (paired)
 * co-array angles also seed the local (no-CAPE) estimator
 * in-system angle hit rate / RMSE of the channel estimator
 * classic DOA test (physical vs co-array MUSIC, resolvable sources)
 * figures T11-T17, list-detector memory fix, stale-pool close
Every edit is anchored and must match exactly once."""
import pathlib
import sys

REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent
f = REPO / 'Code_T_ArrayComparison.m'
src = f.read_text()
fig = (HERE / 'figblock_t2.m').read_text()


def rep(a, b, text=None):
    global src
    t = src if text is None else text
    n = t.count(a)
    if n != 1:
        sys.exit(f'anchor found {n} times:\n{a[:300]}')
    t = t.replace(a, b)
    if text is None:
        src = t
    return t


def cut(start, end, new):
    global src
    for a in (start, end):
        if src.count(a) != 1:
            sys.exit(f'cut anchor found {src.count(a)} times:\n{a[:300]}')
    i, j = src.index(start), src.index(end)
    assert j > i
    src = src[:i] + new + src[j:]


# ---------------------------------------------------------------- header
rep("clear; close all; clc;\n%% ====================================================================\n%  CODE T : ARRAY COMPARISON",
    "clear; close all; clc;\ntry, delete(gcp('nocreate')); catch, end   % close a parallel pool left open by a previous run (frees memory)\n"
    "%% ====================================================================\n%  CODE T : ARRAY COMPARISON")
rep("""%  SPEED: vectorised detector (same decisions as the symbol-by-symbol one),
%  no unused perfect-CSI pass, one parallel task per (drop, scheme, SNR).
%  The script prints an estimated run time after the first task.
%  Figures T1..T10. Results are saved to codeT_arraycomp_results.mat BEFORE
%  plotting; plot_codeT_figures.m redraws the figures without rerunning.""",
"""%  WITH vs WITHOUT CAPE: every scheme also runs without CAPE (local per-AP
%  polar-domain OMP; with co-array processing the co-array angles seed it),
%  on the same drops, channels, pilots and symbols, so the CAPE benefit for
%  every detector (incl. the proposed Cross-AP List) is a paired comparison.
%  CO-ARRAY METRICS: (a) true-lobe angle hit rate / RMSE of the channel
%  estimator; (b) a classic DOA test of the same antennas: physical MUSIC vs
%  co-array MUSIC, sources found vs number of sources, RMSE vs SNR.
%  SETTINGS OF THE STUDY: fast_mode, nSetE, nSetK, nRealE, nSymE in the
%  [CODE T] ARRAY SET block (the older settings at the top are NOT used here).
%  SPEED: vectorised detector (same decisions as the symbol-by-symbol one),
%  no unused perfect-CSI pass, one parallel task per (drop, scheme, SNR).
%  The script prints an estimated run time after the first task.
%  Figures T1..T17. Results are saved to codeT_arraycomp_results.mat BEFORE
%  plotting; plot_codeT_figures.m redraws the figures without rerunning.""")

# ------------------------------------------- legacy runtime block (not used)
rep("""fast_mode = false;
if fast_mode
    nSetups = 3;""",
"""legacy_fast_mode = false;   % [CODE T] NOT used by the Code T study: see fast_mode in the ARRAY SET block
if legacy_fast_mode
    nSetups = 3;""")

# ------------------------------------------------------ settings + schemes
NEWSET = r"""% ---- [CODE T] run size. The BER noise comes from the random DROPS (a few bad
%      drops carry most errors), not from the symbols, so the drops are many
%      and the symbols per drop fewer. ----
fast_mode = false;    % true: quick check;  false: paper run
if fast_mode
    nSetE = 30;   nSetK = 20;
else
    nSetE = 200;  nSetK = 100;   % random AP/user drops: SNR study / user-count sweep
end
nRealE = 2;                   % channel realisations per drop
nSymE  = 400;                 % symbols per realisation: 200 x 2 x 400 x 16 users x 4 bits = 1e7 bits per point

eqN_study   = true;
tB_names    = {'ULA lambda/2', sprintf('ULA %g lambda', ulaw_d_lambda), 'Nested', ...
               'Nested + co-array', 'Coprime', 'Coprime + co-array'};
tB_leg      = {'ULA (\lambda/2)', sprintf('ULA (%g\\lambda)', ulaw_d_lambda), 'Nested', ...
               'Nested + co-array', 'Coprime', 'Coprime + co-array'};
if sp_equal_aperture
    tB_leg(3:6) = {'Nested (stretched)', 'Nested (stretched) + co-array', ...
                   'Coprime (stretched)', 'Coprime (stretched) + co-array'};
end
% schemes 1-6: with CAPE (as before); schemes 7-12: the same arrays WITHOUT CAPE
en_names    = [cellfun(@(s) [s ' + CAPE'], tB_names, 'UniformOutput', false), ...
               cellfun(@(s) [s ', no CAPE'], tB_names, 'UniformOutput', false)];
en_leg      = [tB_leg, cellfun(@(s) [s ', no CAPE'], tB_leg, 'UniformOutput', false)];
en_geom     = repmat({'ula', 'ulaw', 'nested', 'nested', 'coprime', 'coprime'}, 1, 2);
en_idea     = repmat([0 0 0 1 0 1], 1, 2);   % 1 = add difference co-array processing (same antennas)
en_cape     = [ones(1, 6), zeros(1, 6)];     % 1 = CAPE, 0 = local per-AP estimation
en_arr      = repmat([1 2 3 3 4 4], 1, 2);   % physical array (1 ULA, 2 wide ULA, 3 nested, 4 coprime)
nSchE       = numel(en_names);
en_caseName = {'Crowded (N = 8, K = 16)'};
en_N0       = 8;               % 8 PHYSICAL antennas per AP
en_K        = 16;              % 16 users (K > N: loaded cell-free system)
nCaseE      = numel(en_N0);
Ksw         = [8 12 16 20 24];          % user-count sweep (all K >= N = 8: crowded)
Ksw_sch     = [1:6, 8, 9, 11];          % CAPE schemes + no-CAPE wide ULA, nested, coprime
                                        % (dense ULA: no near-field links, no CAPE effect)
if sp_equal_aperture           % target aperture = aperture of the wide ULA
    sp_Dtgt = (en_N0(1) - 1) * ulaw_eta * lambda / 2;
else
    sp_Dtgt = [];
end
tolS = 0.5 * lambda * ((en_N0(1) - 1) * ulaw_eta * lambda / 2)^(-1);   % angle hit tolerance: half a beamwidth of 2.8 m

"""
cut("eqN_study   = true;\nen_names    = {'ULA lambda/2 + CAPE'",
    "% ---- prebuild geometry, polar dictionary and co-arrays per (scheme, case) ----", NEWSET)

rep("EN_seD  = nan(4, max(en_K), nSchE, nCaseE, nSetE);  % [CODE O] sized by the largest K",
    "EN_seD  = nan(4, max(en_K), nSchE, nCaseE, nSetE);  % [CODE O] sized by the largest K\n"
    "EN_hit  = nan(nSchE, nCaseE, nSNR, nSetE);          % true-lobe angle hit rate of the estimator\n"
    "EN_rmse = nan(nSchE, nCaseE, nSNR, nSetE);          % angle RMSE of the hits")
rep("            'mod', modOrder, 'tau_c', tau_c);",
    "            'mod', modOrder, 'tau_c', tau_c, 'tolS', tolS);")
rep("T_seD  = nan(4, max(en_K), nTask);",
    "T_seD  = nan(4, max(en_K), nTask);  T_hit = nan(1, nTask);  T_rm = nan(1, nTask);")
rep("""fprintf('  one task takes about %.0f s -> main study at most about %.0f min with %d worker(s)\\n', ...
        tP, tP * nTask / nW / 60, nW);""",
"""fprintf(['  one task takes about %.0f s -> main study at most about %.0f min, user sweep about %.0f min, ' ...
         'with %d worker(s)\\n'], tP, tP * nTask / nW / 60, ...
        tP * numel(Ksw) * numel(Ksw_sch) * nSetK / nW / 60, nW);""")
rep("    T_time(tk) = O.time;  T_fh(tk) = O.nfus;  T_nf(tk) = nfv;",
    "    T_time(tk) = O.time;  T_fh(tk) = O.nfus;  T_nf(tk) = nfv;  T_hit(tk) = O.hit;  T_rm(tk) = O.rmse;")
rep("EN_loc  = permute(reshape(T_loc, sz), [2 3 1 4]);",
    "EN_loc  = permute(reshape(T_loc, sz), [2 3 1 4]);\n"
    "EN_hit  = permute(reshape(T_hit, sz), [2 3 1 4]);\n"
    "EN_rmse = permute(reshape(T_rm, sz), [2 3 1 4]);")

# ------------------------------------------------------- user-count sweep
a = "% ---- [CODE T] user-count sweep at the reference SNR, all six schemes ----"
b = "if eqN_study\n    % [CODE T] save BEFORE plotting"
i, j = src.index(a), src.index(b)
blk = src[i:j]
blk = rep("Ksw     = [8 12 16 20 24];       % all K >= N = 8: crowded\nKsw_sch = 1:nSchE;\n",
          "% Ksw and Ksw_sch are set in the ARRAY SET block; nSetK drops per K.\n", blk)
blk = blk.replace("% ---- [CODE T] user-count sweep at the reference SNR, all six schemes ----",
                  "% ---- [CODE T] user-count sweep at the reference SNR ----")
blk = blk.replace('nSetE', 'nSetK')
DOA = r"""
% ---- [CODE T] classic DOA test of the four physical arrays: physical MUSIC vs
%      co-array MUSIC on the SAME antennas. D far-field sources at random
%      angles inside the array's unambiguous sector (minimum separation),
%      DOA_T snapshots. Shows the co-array's degrees-of-freedom advantage
%      (more sources than antennas), which BER cannot show. ----
DOA_D   = 2:2:20;             % numbers of sources
DOA_D0  = 6;                  % sources for the RMSE-vs-SNR curve
DOA_T   = 200;                % snapshots (one coherence block)
DOA_snr = SNR_dB;             % SNR per antenna per source [dB]
[~, iD] = min(abs(DOA_snr - 10));  DOA_snrD = DOA_snr(iD);
if fast_mode, DOA_nTr = 60; else, DOA_nTr = 300; end
DOA_fP = zeros(4, numel(DOA_D));  DOA_fC = DOA_fP;  DOA_pP = DOA_fP;  DOA_pC = DOA_fP;
DOA_rP = zeros(4, numel(DOA_snr));  DOA_rC = DOA_rP;  DOA_lim = zeros(4, 2);
if eqN_study
fprintf('\n=== [CODE T] DOA test: physical vs co-array MUSIC, %d trials per point ===\n', DOA_nTr);
tDoa = tic;
DOA_arr = [1 2 3 5];
for a = 1:4
    sc = DOA_arr(a);
    rng(4e7 + 100 * a);
    [DOA_fP(a, :), DOA_fC(a, :), DOA_pP(a, :), DOA_pC(a, :), ~, ~, DOA_lim(a, :)] = ...
        doa_test(EN_geo{sc, 1}, EN_CA{sc, 1}, DOA_D, 10^(DOA_snrD / 10), DOA_T, DOA_nTr, ca_Nsub_max, lambda);
    [~, ~, ~, ~, DOA_rP(a, :), DOA_rC(a, :)] = ...
        doa_test(EN_geo{sc, 1}, EN_CA{sc, 1}, DOA_D0, 10.^(DOA_snr / 10), DOA_T, DOA_nTr, ca_Nsub_max, lambda);
end
fprintf('  DOA test done in %.1f min\n', toc(tDoa) / 60);
end

"""
src = src[:i] + blk + DOA + src[j:]

# ------------------------------------------------------ save + figures
fig_ind = ''.join(('    ' + ln if ln.strip() else ln) for ln in fig.splitlines(True))
SAVE = r"""if eqN_study
    % [CODE T] save BEFORE plotting, so a plotting problem can never lose the run
    save('codeT_arraycomp_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', ...
         'en_names', 'en_leg', 'en_geom', 'en_idea', 'en_cape', 'en_arr', 'en_caseName', ...
         'en_N0', 'en_K', 'ulaw_d_lambda', 'sp_equal_aperture', 'modOrder', 'clu_div', ...
         'fast_mode', 'nSetE', 'nSetK', 'nRealE', 'nSymE', 'tolS', 'EN_geo', ...
         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...
         'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum', ...
         'DOA_D', 'DOA_D0', 'DOA_T', 'DOA_snr', 'DOA_snrD', 'DOA_nTr', 'DOA_fP', 'DOA_fC', ...
         'DOA_pP', 'DOA_pC', 'DOA_rP', 'DOA_rC', 'DOA_lim');
    fprintf('\nResults saved to codeT_arraycomp_results.mat\n');

""" + fig_ind + "end\n\n"
cut("if eqN_study\n    % [CODE T] save BEFORE plotting",
    "function [Hhat, C, Hpre, Cpre] = csi_acquire(", SAVE)

# ------------------------------------------- csi_acquire: local + co-array
rep("function [Hhat, C, Hpre, Cpre] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...",
    "function [Hhat, C, Hpre, Cpre, aErr] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...")
rep("""    LN = N * L;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);
    C = zeros(N, N, L, K);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        for t = 1:tau""",
"""    LN = N * L;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);
    C = zeros(N, N, L, K);
    aErr = nan(L, K, nReal);                   % [CODE T] angle error of the estimate (scoring only)
    hasFG = nargin >= 17 && ~isempty(FG);
    hasTrue = hasFG && isfield(FG, 'sTrue');
    caAng = cell(L, 1);                        % [CODE T] co-array angles seed the local estimator
    if hasFG && ~genie && isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)
        caAng = coarray_pilot_angles(H, p, N, L, K, tau, FG);
    end
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        for t = 1:tau""")
rep("""                    [hh, se2] = pdomp_estimate(z, PD, beta(l, k), sw2);
                    Hhat(idx, :, k) = hh;
                    C(:, :, l, k) = mean(se2) * eye(N);""",
"""                    [hh, se2, sA] = pdomp_estimate(z, PD, beta(l, k), sw2, caAng{l});
                    Hhat(idx, :, k) = hh;
                    C(:, :, l, k) = mean(se2) * eye(N);
                    if hasTrue
                        dS = abs(sA - FG.sTrue(l, k));
                        aErr(l, k, :) = reshape(min(dS, 2 - dS), 1, 1, []);
                    end""")

# ---------------------------------------------- pdomp: seeds + angle output
rep("function [hh, se2] = pdomp_estimate(z, PD, beta_lk, sw2)",
    "function [hh, se2, sAng] = pdomp_estimate(z, PD, beta_lk, sw2, sx)")
rep("""    [N, M] = size(z);
    hh  = zeros(N, M);
    se2 = zeros(1, M);
    gam = log(size(PD.W, 2)) + 0.58 + 3;""",
"""% [CODE T] sx (optional): extra start angles for the FIRST atom (co-array
% MUSIC angles, range by a 1-D search). sAng: angle of the first atom.
    if nargin < 5, sx = []; end
    [N, M] = size(z);
    hh  = zeros(N, M);
    se2 = zeros(1, M);
    sAng = nan(1, M);
    gam = log(size(PD.W, 2)) + 0.58 + 3;""")
rep("""            bBest = [];
            mBest = -inf;
            for st = 1:min(PD.nstart, numel(ordc))
                j = ordc(st);
                [b, mg] = refine_atom(r, PD.s(j), PD.q(j), PD);
                if mg > mBest
                    mBest = mg;
                    bBest = b;
                end
            end""",
"""            bBest = [];
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
            end""")

# ----------------------------------------------- CAPE: angle of the estimate
rep("function [Hhat, C, xhat, nfus] = csi_acquire_cape(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG)",
    "function [Hhat, C, xhat, nfus, aErr] = csi_acquire_cape(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG)")
rep("""    nfus = zeros(K, 1);                        % APs whose pilots were fused (pilot fronthaul)
    sw2 = 1 / (p * tau);""",
"""    nfus = zeros(K, 1);                        % APs whose pilots were fused (pilot fronthaul)
    aErr = nan(L, K, nReal);                   % [CODE T] angle error of the estimate (scoring only)
    hasTrue = isfield(FG, 'sTrue');
    sw2 = 1 / (p * tau);""")
rep("""            bLoc = zeros(N, L);  mLoc = zeros(L, 1);  sLoc = zeros(L, FG.J);  qLoc = zeros(L, FG.J);""",
    """            bLoc = zeros(N, L);  mLoc = zeros(L, 1);  sLoc = zeros(L, FG.J);  qLoc = zeros(L, FG.J);
            jLoc = ones(L, 1);""")
rep("""                    if mj > mLoc(l)
                        mLoc(l) = mj;
                        bLoc(:, l) = bj;
                    end""",
"""                    if mj > mLoc(l)
                        mLoc(l) = mj;
                        bLoc(:, l) = bj;
                        jLoc(l) = j;
                    end""")
rep("""                [b, mf] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);        % fused lobe, local precision
                sameLobe = abs(bLoc(:, l)' * b) >= FG.consist;
                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level
                if (sameLobe && mLoc(l) > mf) || clearWin                  % consistency / likelihood gate
                    b = bLoc(:, l);
                end""",
"""                [b, mf, sSel] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);  % fused lobe, local precision
                sameLobe = abs(bLoc(:, l)' * b) >= FG.consist;
                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level
                if (sameLobe && mLoc(l) > mf) || clearWin                  % consistency / likelihood gate
                    b = bLoc(:, l);
                    sSel = sLoc(l, jLoc(l));
                end
                if hasTrue
                    dS = abs(sSel - FG.sTrue(l, k));
                    aErr(l, k, m) = min(dS, 2 - dS);
                end""")

# ------------------------------------------------------------- eval_scheme
rep("""        [Hh, Cs, xh, nfus] = csi_acquire_cape(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG);""",
    """        [Hh, Cs, xh, nfus, aErr] = csi_acquire_cape(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG);""")
rep("""        [Hh, Cs] = csi_acquire(Ht, ST.R2, ST.NF, ST.hD, bet, (1:K)', K, p, ...
                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD);""",
"""        FGl = FG;  FGl.on = false;              % [CODE T] local estimation (+ co-array seeds)
        [Hh, Cs, ~, ~, aErr] = csi_acquire(Ht, ST.R2, ST.NF, ST.hD, bet, (1:K)', K, p, ...
                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD, FGl);""")
rep("""    O.time = toc(tE);
    O.nmse = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);""",
"""    O.time = toc(tE);
    a = aErr(isfinite(aErr));                  % [CODE T] true-lobe hit rate of the estimator
    O.hit = NaN;  O.rmse = NaN;
    if ~isempty(a)
        hit = a < EP.tolS;
        O.hit = mean(hit);
        if any(hit), O.rmse = sqrt(mean(a(hit) .^ 2)); end
    end
    O.nmse = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);""")

# ---------------------------------------------------------------- run_task
rep("    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;",
    "    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;\n"
    "    FG.sTrue = D.sE;                          % [CODE T] true angles: used ONLY to score the angle hit rate\n"
    "    PD.nRng = CF.nRng;")
rep("                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF);",
    "                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF, 'tolS', CF.tolS);")

# --------------------------------------------- list detector memory (chunk)
rep("    chunk = 256;", "    chunk = 64;    % symbol block of the list detector (memory only; results do not change)")

# ------------------------------------------------------------ new functions
NEWFUN = r"""

%% ====================================================================
%  [CODE T] CO-ARRAY ANGLES FOR THE LOCAL ESTIMATOR AND THE DOA TEST
%% ====================================================================
function caAng = coarray_pilot_angles(H, p, N, L, K, tau, FG)
% Per-AP difference-co-array MUSIC angles from the raw pilot snapshots, with
% the SAME side random stream and the same steps as csi_acquire_cape, so the
% local and the CAPE variant of a co-array scheme see identical angles.
    caAng = cell(L, 1);
    Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);
    useSide = isfield(FG, 'caSeed');
    if useSide, rsMain = rng; rng(FG.caSeed); end
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        Hl = reshape(H(idx, 1, :), [N, K]);
        Yl = sqrt(p) * Hl * Phi(1:K, :) + sqrt(0.5) * (randn(N, tau) + 1j * randn(N, tau));
        Rl = (Yl * Yl') / tau;
        caAng{l} = coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
    end
    if useSide, rng(rsMain); end
end

function [fP, fC, pP, pC, rP, rC, lim] = doa_test(G, CA, Dlist, snrList, T, nTr, nSubMax, lambda)
% Classic DOA test of one array: D far-field equal-power sources at random
% angles in sin(theta) inside the unambiguous sector of the array (width
% 2/u for a grid of u half-wavelengths, at most 1.9), with a minimum
% separation of half the mean spacing; T snapshots. Physical MUSIC (N x N covariance, at most N - 1 sources) vs
% difference-co-array MUSIC (spatially smoothed, at most Ns - 1 sources) on
% the SAME snapshots. A source is found if a peak lies within half the
% minimum separation (one-to-one matching).
%   fP/fC : fraction of sources found, pP/pC : probability that ALL are found
%   rP/rC : angle RMSE of the found sources [sin(theta)]
%   lim   : [physical, co-array] maximum number of sources
% Dlist or snrList may be a vector (the other one a scalar).
    pos = G.pos(:);  N = numel(pos);
    per = min(2, lambda / G.unit);  w = min(per, 1.9);
    sg = linspace(-w / 2, w / 2, 4001);
    Ap = exp(1j * 2 * pi * pos * sg / lambda);
    Ns = min(CA.Nsub, nSubMax);
    Av = exp(1j * 2 * pi * (0:Ns - 1)' * CA.d0 * sg / lambda);
    lim = [N - 1, Ns - 1];
    nP = max(numel(Dlist), numel(snrList));
    fP = zeros(1, nP);  fC = fP;  pP = fP;  pC = fP;  rP = fP;  rC = fP;
    for ip = 1:nP
        D = Dlist(min(ip, numel(Dlist)));  snr = snrList(min(ip, numel(snrList)));
        Lw = 0.95 * w;  dmin = 0.5 * Lw / D;  tol = dmin / 2;
        nFP = 0; nFC = 0; aP = 0; aC = 0; eP = 0; eC = 0;
        for tr = 1:nTr
            s = -Lw / 2 + dmin * (0:D - 1) + sort(rand(1, D)) * (Lw - (D - 1) * dmin);
            A = exp(1j * 2 * pi * pos * s / lambda);
            Y = sqrt(snr / 2) * A * (randn(D, T) + 1j * randn(D, T)) + sqrt(0.5) * (randn(N, T) + 1j * randn(N, T));
            R = (Y * Y') / T;
            % physical MUSIC
            [U, Dg] = eig((R + R') / 2);  [~, o] = sort(real(diag(Dg)), 'descend');  U = U(:, o);
            En = U(:, min(D, N - 1) + 1:end);
            Pm = 1 ./ max(sum(abs(En' * Ap) .^ 2, 1), eps);
            [k1, e1] = match_peaks(s, coarray_peaks(sg, Pm, D), tol, per);
            % co-array MUSIC (same snapshots)
            z = R(:);  zu = accumarray(CA.ic, z) ./ max(accumarray(CA.ic, ones(size(z))), 1);
            zv = zu(CA.loc);  nWin = CA.Nv - Ns + 1;  Rss = zeros(Ns);
            for i = 1:nWin
                zi = zv(i:i + Ns - 1);  Rss = Rss + zi * zi';
            end
            [U, Dg] = eig((Rss + Rss') / 2);  [~, o] = sort(real(diag(Dg)), 'descend');  U = U(:, o);
            En = U(:, min(D, Ns - 1) + 1:end);
            Pc = 1 ./ max(sum(abs(En' * Av) .^ 2, 1), eps);
            [k2, e2] = match_peaks(s, coarray_peaks(sg, Pc, D), tol, per);
            nFP = nFP + k1;  nFC = nFC + k2;  aP = aP + (k1 == D);  aC = aC + (k2 == D);
            eP = eP + e1;  eC = eC + e2;
        end
        fP(ip) = nFP / (D * nTr);  fC(ip) = nFC / (D * nTr);
        pP(ip) = aP / nTr;  pC(ip) = aC / nTr;
        rP(ip) = sqrt(eP / max(nFP, 1));  rC(ip) = sqrt(eC / max(nFC, 1));
        if nFP == 0, rP(ip) = NaN; end
        if nFC == 0, rC(ip) = NaN; end
    end
end

function [nF, e2] = match_peaks(s, pk, tol, per)
% Greedy one-to-one matching of true angles s and peaks pk (periodic distance
% with period per). Returns the number matched within tol and the sum of
% their squared errors.
    nF = 0;  e2 = 0;
    if isempty(pk), return; end
    Dm = abs(mod(s(:) - pk(:).' + per / 2, per) - per / 2);
    while true
        [m, ix] = min(Dm(:));
        if isempty(m) || ~isfinite(m) || m >= tol, break; end
        [r, c] = ind2sub(size(Dm), ix);
        nF = nF + 1;  e2 = e2 + m ^ 2;
        Dm(r, :) = inf;  Dm(:, c) = inf;
    end
end
"""
src = src.rstrip() + '\n' + NEWFUN
f.write_text(src)

plot = (r"""%% PLOT_CODET_FIGURES  Redraw Code T figures T1-T17 without rerunning the simulation.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_hit', 'var') || ~exist('DOA_fC', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
""" + fig)
(REPO / 'plot_codeT_figures.m').write_text(plot)
print('Code T lines:', src.count('\n'), ' plot lines:', plot.count('\n'))

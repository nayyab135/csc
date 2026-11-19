"""Speed patch for Code T: vectorised detector, no unused perfect-CSI pass, no meshgrid in the
CAPE refinement, and task-level parallelism. Every edit is anchored and must match exactly once."""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


def _once(src, old):
    n = src.count(old)
    if n != 1:
        sys.exit(f'speed_patch: anchor found {n} times:\n{old[:160]}')


def _rep(src, old, new):
    _once(src, old)
    return src.replace(old, new)


def _cut(src, start, end, new):
    _once(src, start)
    _once(src, end)
    i, j = src.index(start), src.index(end)
    assert j > i
    return src[:i] + new + src[j:]


MAIN = r"""r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);
% ---- [CODE T] speed: the random drops are drawn once (reproducible) and the
%      work is split into one parallel task per (drop, case, scheme, SNR), so
%      every worker stays busy. Each task uses the same seed as before. ----
CF = struct('lambda', lambda, 'eps_NF', eps_NF, 'asd_phi', ASD_varphi, 'asd_th', ASD_theta, ...
            'clu_div', clu_div, 'sL', squareLen, 'hD', hDiff, 'nAP', cape_nAP, 'nCand', cape_nCand, ...
            'J', cape_J, 'step0', cape_step0, 'nIter', cape_nIter, 'consist', cape_consist, ...
            'nAng', ca_nAngGrid, 'nRng', ca_nRangeGrid, 'nPeak', ca_nPeak, 'nSubMax', ca_Nsub_max, ...
            'nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, 'maxBr', maxBr, ...
            'mod', modOrder, 'tau_c', tau_c);
rng(20261119);
DRP = cell(nSetE, nCaseE);
for ns = 1:nSetE
    for caE = 1:nCaseE
        DRP{ns, caE} = make_drop(en_K(caE), L, n_collinear_pairs, r_min_NF, r_max_en, hs_r, squareLen, ...
                                 wrap_fix, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    end
end
sz    = [nSNR, nSchE, nCaseE, nSetE];
nTask = prod(sz);
T_nmse = zeros(1, nTask);  T_be = zeros(4, nTask);  T_sum = zeros(4, nTask);  T_loc = nan(1, nTask);
T_time = zeros(1, nTask);  T_fh = zeros(1, nTask);  T_nf  = zeros(1, nTask);
T_seD  = nan(4, max(en_K), nTask);
nW = n_workers();
[~, scP] = max(en_idea);                         % heaviest scheme (co-array) for the time estimate
tP = tic;
run_task(DRP{1, 1}, EN_geo{scP, 1}, EN_PD{scP, 1}, EN_CA{scP, 1}, EN_CAs{scP, 1}, en_idea(scP), ...
         en_cape(scP), en_K(1), SNR_lin(snr_ref_idx), 1, CF);
tP = toc(tP);
fprintf('  one task takes about %.0f s -> main study at most about %.0f min with %d worker(s)\n', ...
        tP, tP * nTask / nW / 60, nW);
tMain = tic;
parfor tk = 1:nTask
    [si, scE, caE, ns] = ind2sub(sz, tk);
    seed = 9e6 + 1e4 * ns + 1000 * caE + si;     % same seed for every scheme (paired)
    [O, nfv] = run_task(DRP{ns, caE}, EN_geo{scE, caE}, EN_PD{scE, caE}, EN_CA{scE, caE}, EN_CAs{scE, caE}, ...
                        en_idea(scE), en_cape(scE), en_K(caE), SNR_lin(si), seed, CF);
    T_nmse(tk) = O.nmse;  T_be(:, tk) = O.be;  T_sum(:, tk) = O.sum;
    if ~isempty(O.loc), T_loc(tk) = median(O.loc); end
    T_time(tk) = O.time;  T_fh(tk) = O.nfus;  T_nf(tk) = nfv;
    sd = nan(4, max(en_K));
    if si == snr_ref_idx, sd(:, 1:size(O.se, 2)) = O.se; end
    T_seD(:, :, tk) = sd;
end
fprintf('  main study: %d tasks done in %.1f min\n', nTask, toc(tMain) / 60);
EN_nmse = permute(reshape(T_nmse, sz), [2 3 1 4]);
EN_be   = permute(reshape(T_be, [4 sz]), [1 3 4 2 5]);
EN_sum  = permute(reshape(T_sum, [4 sz]), [1 3 4 2 5]);
EN_loc  = permute(reshape(T_loc, sz), [2 3 1 4]);
EN_time = reshape(mean(reshape(T_time, sz), 1), nSchE, nCaseE, nSetE);
EN_fh   = reshape(mean(reshape(T_fh, sz), 1), nSchE, nCaseE, nSetE);
EN_nf   = reshape(mean(reshape(T_nf, sz), 1), nSchE, nCaseE, nSetE);
seAll   = reshape(T_seD, [4, max(en_K), sz]);
EN_seD  = reshape(seAll(:, :, snr_ref_idx, :, :, :), [4, max(en_K), nSchE, nCaseE, nSetE]);
end

"""

KSW = r"""if eqN_study
fprintf('\n=== [CODE T] user-count sweep at %d dB: K = %s ===\n', SNR_dB(snr_ref_idx), mat2str(Ksw));
DRK = cell(nSetE, nKsw);
for ns = 1:nSetE
    for ik = 1:nKsw
        DRK{ns, ik} = make_drop(Ksw(ik), L, n_collinear_pairs, r_min_NF, r_max_en, hs_r, squareLen, ...
                                wrap_fix, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    end
end
nTk  = nKsch * nKsw * nSetE;
K_be = zeros(1, nTk);  K_sum = zeros(1, nTk);
tSw  = tic;
parfor tk = 1:nTk
    [jj, ik, ns] = ind2sub([nKsch, nKsw, nSetE], tk);
    scE = Ksw_sch(jj);
    O = run_task(DRK{ns, ik}, EN_geo{scE, 1}, EN_PD{scE, 1}, EN_CA{scE, 1}, [], en_idea(scE), ...
                 en_cape(scE), Ksw(ik), SNR_lin(snr_ref_idx), 3e7 + 1e4 * ns + 100 * ik, CF);
    K_be(tk) = O.be(4);  K_sum(tk) = O.sum(4);
end
EN_Kbe  = reshape(K_be, nKsch, nKsw, nSetE);
EN_Ksum = reshape(K_sum, nKsch, nKsw, nSetE);
fprintf('  user-count sweep: %d tasks done in %.1f min\n', nTk, toc(tSw) / 60);
end

"""

HELPERS = r"""
function D = make_drop(Kc, L, nCol, rmin, rmax, hs_r, sL, wfix, hD, fc_GHz, hBS, hUT, c0, noisePow_dBm)
% [CODE T] one random AP/user drop with uplink power control (received SNR at
% each user's strongest AP equals the SNR axis value).
    D.APp = (rand(L, 1) + 1j * rand(L, 1)) * sL;
    UEp = drop_scenario(D.APp, Kc, L, min(nCol, floor(Kc / 2)), rmin, rmax, 'hotspot', hs_r, sL);
    D.UEm = mod(real(UEp), sL) + 1j * mod(imag(UEp), sL);
    D.APw = wrap_images(D.APp, sL, wfix);
    [D.d3, D.sE, D.phL, D.thL, bet] = link_geometry(D.APw, UEp, hD, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    D.bet = bet ./ max(bet, [], 1);
end

function [O, nf] = run_task(D, G, PD, CA, CAs, idea, useCape, Kc, p, seed, CF)
% [CODE T] one simulation task: one drop, one array/processing scheme, one SNR.
    ST = array_stats(G, D.d3, D.sE, D.phL, D.thL, D.bet, CF.lambda, CF.eps_NF, CF.asd_phi, CF.asd_th, CF.clu_div);
    nf = mean(ST.NF(:));
    FG = cape_cfg(true, D.APp, D.APw, CF.sL, CF.hD, CF.nAP, CF.nCand, CF.J, CF.step0, CF.nIter, CF.consist);
    FG.useCA = (idea > 0);  FG.co_idea = idea;  FG.CA = CA;  FG.CAs = CAs;
    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    EP = struct('nReal', CF.nReal, 'nSym', CF.nSym, 'M', CF.M, 'dth', CF.dth, 'maxBr', CF.maxBr, ...
                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF);
    O = eval_scheme(ST, D.bet, useCape, PD, FG, p, EP, false, D.UEm, seed);
end

function nW = n_workers()
% number of parallel workers (1 without a pool or without the Parallel Toolbox)
    nW = 1;
    try
        pl = gcp('nocreate');
        if isempty(pl), pl = gcp; end
        nW = max(1, pl.NumWorkers);
    catch
    end
end
"""


def patch(src):
    # 1. task-level parallel main study
    src = _cut(src, "r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);\nparfor ns = 1:nSetE\n",
               "% ---- [CODE T] user-count sweep at the reference SNR", MAIN)
    # 2. task-level parallel user sweep
    src = _cut(src, "if eqN_study\nfprintf('\\n=== [CODE T] user-count sweep at %d dB: K = %s ===\\n'",
               "if eqN_study\n    % [CODE T] save BEFORE plotting", KSW)
    # 3. vectorised detector (bit-identical decisions, verified against the original)
    fast = (HERE / 'prof' / 'metrics_case_fast.m').read_text()
    fast = fast.replace('function [berU, nmseU, bits_pu, eng_pu, etaXap] = metrics_case_fast(',
                        'function [berU, nmseU, bits_pu, eng_pu, etaXap] = metrics_case(', 1)
    src = _cut(src, "function [berU, nmseU, bits_pu, eng_pu, etaXap] = metrics_case(",
               "function [beL, beH, bits] = ber_linsic(", fast.rstrip() + "\n\n")
    # 4. eval_scheme: drop the perfect-CSI detection pass (never used in Code T)
    src = _rep(src, "    Cs0 = zeros(size(Cs));\n    bE = zeros(4, K);  nE = zeros(4, K);  bP = zeros(4, K);  nb = 0;  ne = 0;  eX = 0;\n",
               "    bE = zeros(4, K);  nE = zeros(4, K);  nb = 0;  ne = 0;  eX = 0;\n")
    src = _rep(src, "        bu = metrics_case(Hm, Hm, Cs0, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);\n        bP = bP + bu;\n", "")
    src = _rep(src, "    O.bp  = sum(bP, 2) / max(K * nb, 1);\n",
               "    O.bp  = nan(4, 1);   % [CODE T] perfect-CSI pass removed: never used here, it doubled the run time\n")
    # 5. CAPE refinement window without meshgrid (same 25 points, same order)
    src = _rep(src, "    k0 = (1 - s0 ^ 2) * q0;\n    for rnd = 1:PD.nref\n",
               "    k0 = (1 - s0 ^ 2) * q0;\n    iS = ceil((1:25) / 5);  iQ = mod(0:24, 5) + 1;   % [CODE T] 5 x 5 window, meshgrid order\n"
               "    for rnd = 1:PD.nref\n")
    src = _rep(src, "            [SS, KK] = meshgrid(sc, kc);\n            SS = SS(:).';\n            QQ = min(KK(:).' ./ max(1 - SS .^ 2, 1e-6), PD.q_max);\n",
               "            SS = sc(iS);  KK = kc(iQ);\n            QQ = min(KK ./ max(1 - SS .^ 2, 1e-6), PD.q_max);\n")
    src = _rep(src, "            [SS, QQ] = meshgrid(sc, qc);\n            SS = SS(:).';\n            QQ = QQ(:).';\n",
               "            SS = sc(iS);  QQ = qc(iQ);\n")
    # 6. helpers
    src = src.rstrip() + "\n" + HELPERS
    # 7. header note
    src = _rep(src, "%  Figures T1..T10. Results are saved to codeT_arraycomp_results.mat BEFORE\n",
               "%  SPEED: vectorised detector (same decisions as the symbol-by-symbol one),\n"
               "%  no unused perfect-CSI pass, one parallel task per (drop, scheme, SNR).\n"
               "%  The script prints an estimated run time after the first task.\n"
               "%  Figures T1..T10. Results are saved to codeT_arraycomp_results.mat BEFORE\n")
    return src

#!/usr/bin/env python3
"""Build Code_U_MovableSparse.m from Code_T_ArrayComparison.m (the sped-up Code T).
Every edit is anchored and must match exactly once."""
import pathlib
import sys

REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent   # docs/codeU_build
src = (REPO / 'Code_T_ArrayComparison.m').read_text()
fig = (HERE / 'figblock_u.m').read_text()


def once(text, a):
    n = text.count(a)
    if n != 1:
        sys.exit(f'anchor found {n} times:\n{a[:200]}')


def cut(text, start, end, new):
    once(text, start)
    once(text, end)
    i, j = text.index(start), text.index(end)
    assert j > i
    return text[:i] + new + text[j:]


HEADER = r"""clear; close all; clc;
delete(gcp('nocreate'));   % close any parallel pool left open by a previous run (frees memory)
%% ====================================================================
%  CODE U : MOVABLE (SLIDING-MIRROR) SPARSE ARRAY vs STATIC SPARSE ARRAYS
%           crowded near-field / hybrid-field cell-free uplink (from Code T)
%  --------------------------------------------------------------------
%  Same setting as Code T: L = 8 APs, N = 8 antennas per AP, K = 16 users,
%  16-QAM, uplink power control (SNR axis = received SNR per antenna at the
%  user's strongest AP), clusters within 10 dB, hotspot users. The sparse
%  arrays all have the same ~2.8 m aperture, so most AP-user links are near
%  field; the dense lambda/2 ULA (0.35 m) is far field by physics.
%
%  Arrays (each run WITHOUT CAPE = local per-AP estimation, and WITH CAPE):
%    dense ULA | nested | nested + difference co-array | coprime |
%    coprime + difference co-array | MISC | IMISC | EMISC |
%    movable array, uniform half | movable array, sparse (increasing-gap) half
%  MOVABLE ARRAY (Idea 1): two mirror halves at +u and -u slide outward by
%  lambda/4 per frame (nPilotBlk frames). Every mirror pair gives the lag 2u
%  with the near-field range term cancelled; over all frames the lags fill a
%  hole-free virtual array (angle candidates). The pilots of all frames are
%  also stacked into a synthetic aperture for the range and gain fit. Data
%  are detected with the array parked at frame sma_sig_d.
%  FAIRNESS: every scheme gets the SAME total pilot energy, that of one Code T
%  pilot block. Static arrays use it as in Code T; the movable array splits it
%  over its nPilotBlk frames (each frame 1/nPilotBlk of the energy), so it gets
%  no extra pilot energy. Static arrays therefore run exactly as in Code T.
%  MISC / IMISC / EMISC need their inter-element spacing sets (ies_misc,
%  ies_imisc, ies_emisc below); an empty set skips that array.
%  fast_mode = true  -> quick run;  fast_mode = false -> full paper run.
%  Figures U1..U13. Results saved to codeU_movable_results.mat BEFORE plotting;
%  plot_codeU_figures.m redraws the figures without rerunning.
%% ====================================================================
"""
src = cut(src, 'clear; close all; clc;\n%% ====', '\nclear;\nclose all;\nclc;\n', HEADER)

UMAIN = r"""%% ====================================================================
%  [CODE U] SETTINGS (the Code T crowded setting, unchanged) AND SCHEMES
%% ====================================================================
fast_mode = true;     % true: quick run;  false: full paper run (see the run-time print)
if fast_mode
    nSetE = 6;   nRealE = 2;  nSymE = 300;  Ksw = [8 16 24];
else
    nSetE = 10;  nRealE = 4;  nSymE = 800;  Ksw = [8 12 16 20 24];   % = Code T
end
ulaw_d_lambda = 4;                    % spacing that sets the common aperture (as Code T)
ulaw_eta      = 2 * ulaw_d_lambda;
ca_Nsub_max   = 1000;                 % co-array MUSIC: full smoothed virtual aperture
modOrder = 16;                        % 16-QAM (as Code T)
clu_div  = 10;                        % clusters: APs within 10 dB of the strongest AP
SNR_dB   = 0:4:28;                    % received SNR per antenna at the strongest AP [dB]
SNR_lin  = 10.^(SNR_dB / 10);
nSNR     = numel(SNR_dB);
[~, snr_ref_idx] = min(abs(SNR_dB - 10));
snr_ref_dB = SNR_dB(snr_ref_idx);
N0u  = 8;                             % antennas per AP
Ku   = 16;                            % users
Dtgt = (N0u - 1) * ulaw_eta * lambda / 2;   % common aperture of the sparse arrays (2.8 m)
nPilotBlk = 24;                       % frames of the movable array (they share ONE Code T pilot block's energy)
sma_sig_d = 7;                        % data frame of the movable array (frame 7: +-2 ... +-14 lambda)
tolS = 0.5 * lambda / Dtgt;           % angle hit tolerance in sin(theta): half a beamwidth of 2.8 m
% ---- MISC family: inter-element spacing (IES) sets for N0u = 8 antennas, in units of the grid step.
%      Copy them from the papers (N0u - 1 = 7 spacings each). An empty set skips that array.
ies_misc  = [];   % MISC : Zheng, Wang, Kong, Zhang, IEEE TSP 2019
ies_imisc = [];   % IMISC: arXiv 2208.09574
ies_emisc = [];   % EMISC: Sheng, Lu, Li, de Lamare, arXiv 2309.09044

eqN_study = true;
bs_name = {'Dense ULA', 'Nested', 'Nested + co-array', 'Coprime', 'Coprime + co-array', ...
           'MISC', 'IMISC', 'EMISC', 'Movable (uniform half)', 'Movable (sparse half)'};
bs_geom = {'ula', 'nested', 'nested', 'coprime', 'coprime', 'misc', 'imisc', 'emisc', 'sma_u', 'sma_s'};
bs_kind = {'plain', 'plain', 'coarray', 'plain', 'coarray', 'plain', 'plain', 'plain', 'sma', 'sma'};
bs_ies  = {[], [], [], [], [], ies_misc, ies_imisc, ies_emisc, [], []};
keepB = true(1, numel(bs_name));
for b = 6:8
    if isempty(bs_ies{b})
        keepB(b) = false;
        fprintf('NOTE: %s skipped: its spacing set ies_%s is empty (fill it in from the paper).\n', ...
                bs_name{b}, bs_geom{b});
    end
end
bs_name = bs_name(keepB);  bs_geom = bs_geom(keepB);  bs_kind = bs_kind(keepB);  bs_ies = bs_ies(keepB);
nB = numel(bs_name);
en_base = kron(1:nB, [1 1]);          % scheme -> base scheme
en_cape = repmat([0 1], 1, nB);       % 0 = local estimation, 1 = CAPE
nSchE   = numel(en_base);
en_names = cell(1, nSchE);
for sc = 1:nSchE
    en_names{sc} = sprintf('%s, %s', bs_name{en_base(sc)}, ternary_str(en_cape(sc), 'CAPE', 'local'));
end
en_caseName = {'Crowded (N = 8, K = 16)'};
en_N0 = N0u;  en_K = Ku;  nCaseE = 1;

% ---- geometry, dictionaries and estimation set-up of every base scheme ----
EN_geo = cell(1, nB);  BPD = cell(1, nB);  BES = cell(1, nB);
fprintf('\n=== [CODE U] arrays (N = %d antennas, common aperture %.2f m, movable array %d frames) ===\n', N0u, Dtgt, nPilotBlk);
for b = 1:nB
    G = geom_u(bs_geom{b}, N0u, lambda, Dtgt, nPilotBlk, sma_sig_d, bs_ies{b});
    PDb = build_polar_dict(G.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
    PDb.smax = pd_smax;  PDb.stop = pd_stop;  PDb.nref = pd_nref;  PDb.nstart = pd_nstart;
    ES = struct('kind', bs_kind{b}, 'F', nPilotBlk, 'PDe', PDb, 'Ne', G.N, 'mir', [], 'CA', []);
    if strcmp(bs_kind{b}, 'sma')
        PDe = build_polar_dict(G.frames(:), lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
        PDe.smax = pd_smax;  PDe.stop = pd_stop;  PDe.nref = pd_nref;  PDe.nstart = pd_nstart;
        ES.PDe = PDe;  ES.Ne = numel(G.frames);  ES.mir = G.mir;
    elseif strcmp(bs_kind{b}, 'coarray')
        ES.CA = coarray_setup(G.pos, lambda, 'diff', G.unit);
    end
    EN_geo{b} = G;  BPD{b} = PDb;  BES{b} = ES;
    fprintf('  %-26s aperture %5.2f m  Rayleigh %6.0f m  smallest gap %4.2f lambda\n', bs_name{b}, ...
            G.aperture, 2 * G.aperture^2 / lambda, G.min_gap / lambda);
end

% ---- main study: one parallel task per (SNR, scheme, drop) ----
r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);
CF = struct('lambda', lambda, 'eps_NF', eps_NF, 'asd_phi', ASD_varphi, 'asd_th', ASD_theta, ...
            'clu_div', clu_div, 'sL', squareLen, 'hD', hDiff, 'nAP', cape_nAP, 'nCand', cape_nCand, ...
            'J', cape_J, 'step0', cape_step0, 'nIter', cape_nIter, 'consist', cape_consist, ...
            'nAng', ca_nAngGrid, 'nRng', ca_nRangeGrid, 'nPeak', ca_nPeak, 'nSubMax', ca_Nsub_max, ...
            'nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, 'maxBr', maxBr, ...
            'mod', modOrder, 'tau_c', tau_c, 'tolS', tolS);
rng(20261119);
DRP = cell(nSetE, 1);
for ns = 1:nSetE
    DRP{ns} = make_drop(Ku, L, n_collinear_pairs, r_min_NF, r_max_en, hs_r, squareLen, ...
                        wrap_fix, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
end
sz    = [nSNR, nSchE, nSetE];
nTask = prod(sz);
T_nmse = zeros(1, nTask);  T_be = zeros(4, nTask);  T_sum = zeros(4, nTask);  T_loc = nan(1, nTask);
T_hit  = nan(1, nTask);    T_rm = nan(1, nTask);    T_time = zeros(1, nTask);  T_fh = zeros(1, nTask);
T_nf   = zeros(1, nTask);  T_seD = nan(4, Ku, nTask);
nW = n_workers();
% run-time estimate: one CAPE task of every array at the reference SNR
tB = zeros(1, nB);
parfor b = 1:nB
    t0 = tic;
    run_task_u(DRP{1}, EN_geo{b}, BPD{b}, BES{b}, 1, Ku, SNR_lin(snr_ref_idx), 1, CF);
    tB(b) = toc(t0);
end
fprintf('\n  task time per array [s]: %s\n', sprintf('%.1f ', tB));
fprintf(['  estimated run time with %d worker(s): main study about %.0f min, user sweep about %.0f min ' ...
         '(fast_mode = %d)\n'], nW, 2 * sum(tB) * nSNR * nSetE / nW / 60, ...
        sum(tB) * nSetE * numel(Ksw) / nW / 60, fast_mode);
fprintf('\n=== [CODE U] main study: %d schemes x %d SNR x %d drops = %d tasks ===\n', nSchE, nSNR, nSetE, nTask);
tMain = tic;
parfor tk = 1:nTask
    [si, scE, ns] = ind2sub(sz, tk);
    b = en_base(scE);
    seed = 9e6 + 1e4 * ns + 1000 + si;           % same seed for every scheme (paired; as Code T)
    [O, nfv] = run_task_u(DRP{ns}, EN_geo{b}, BPD{b}, BES{b}, en_cape(scE), Ku, SNR_lin(si), seed, CF);
    T_nmse(tk) = O.nmse;  T_be(:, tk) = O.be;  T_sum(:, tk) = O.sum;
    if ~isempty(O.loc), T_loc(tk) = median(O.loc); end
    T_hit(tk) = O.hit;  T_rm(tk) = O.rmse;  T_time(tk) = O.time;  T_fh(tk) = O.nfus;  T_nf(tk) = nfv;
    sd = nan(4, Ku);
    if si == snr_ref_idx, sd(:, 1:size(O.se, 2)) = O.se; end
    T_seD(:, :, tk) = sd;
end
fprintf('  main study done in %.1f min\n', toc(tMain) / 60);
EN_nmse = reshape(permute(reshape(T_nmse, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_be   = reshape(permute(reshape(T_be, [4 sz]), [1 3 2 4]), 4, nSchE, 1, nSNR, nSetE);
EN_sum  = reshape(permute(reshape(T_sum, [4 sz]), [1 3 2 4]), 4, nSchE, 1, nSNR, nSetE);
EN_loc  = reshape(permute(reshape(T_loc, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_hit  = reshape(permute(reshape(T_hit, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_rmse = reshape(permute(reshape(T_rm, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_time = reshape(mean(reshape(T_time, sz), 1), nSchE, 1, nSetE);
EN_fh   = reshape(mean(reshape(T_fh, sz), 1), nSchE, 1, nSetE);
EN_nf   = reshape(mean(reshape(T_nf, sz), 1), nSchE, 1, nSetE);
seAll   = reshape(T_seD, [4, Ku, sz]);
EN_seD  = reshape(seAll(:, :, snr_ref_idx, :, :), 4, Ku, nSchE, 1, nSetE);

"""
src = cut(src, '%% ====================================================================\n%  [CODE T] ARRAY SET',
          '% ---- [CODE T] user-count sweep at the reference SNR', UMAIN)

UKSW = r"""% ---- [CODE U] user-count sweep at the reference SNR (schemes with CAPE) ----
Ksw_sch = find(en_cape == 1);
nKsw = numel(Ksw);  nKsch = numel(Ksw_sch);
fprintf('\n=== [CODE U] user-count sweep at %d dB: K = %s ===\n', SNR_dB(snr_ref_idx), mat2str(Ksw));
rng(20261120);
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
    b = en_base(Ksw_sch(jj));
    O = run_task_u(DRK{ns, ik}, EN_geo{b}, BPD{b}, BES{b}, 1, Ksw(ik), SNR_lin(snr_ref_idx), ...
                   3e7 + 1e4 * ns + 100 * ik, CF);
    K_be(tk) = O.be(4);  K_sum(tk) = O.sum(4);
end
EN_Kbe  = reshape(K_be, nKsch, nKsw, nSetE);
EN_Ksum = reshape(K_sum, nKsch, nKsw, nSetE);
fprintf('  user-count sweep done in %.1f min\n', toc(tSw) / 60);

"""
src = cut(src, '% ---- [CODE T] user-count sweep at the reference SNR',
          'if eqN_study\n    % [CODE T] save BEFORE plotting', UKSW)

fig_ind = ''.join(('    ' + ln if ln.strip() else ln) for ln in fig.splitlines(True))
USAVE = r"""if eqN_study
    % save BEFORE plotting, so a plotting problem can never lose the run
    save('codeU_movable_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', 'Dtgt', 'tolS', ...
         'nPilotBlk', 'sma_sig_d', 'bs_name', 'bs_geom', 'bs_kind', 'en_base', 'en_cape', 'en_names', ...
         'en_caseName', 'en_N0', 'en_K', 'fast_mode', 'nSetE', 'nRealE', 'nSymE', 'modOrder', 'clu_div', ...
         'EN_geo', 'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', ...
         'EN_nf', 'EN_seD', 'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\nResults saved to codeU_movable_results.mat\n');

""" + fig_ind + "end\n\n"
src = cut(src, 'if eqN_study\n    % [CODE T] save BEFORE plotting',
          'function [Hhat, C, Hpre, Cpre] = csi_acquire(', USAVE)

UFUN = r"""
%% ====================================================================
%  [CODE U] MOVABLE ARRAY, GEOMETRIES AND THE UNIFIED ESTIMATOR
%% ====================================================================
function G = geom_u(type, N0, lambda, Dtgt, F, sig_d, ies)
% Geometry of one array with N0 physical antennas. Sparse arrays are put on a
% grid of u half-wavelengths so that their aperture is about Dtgt.
    d0 = lambda / 2;
    switch lower(type)
        case 'ula'
            G = geom_targetN('ula', N0, lambda, 8);
        case {'nested', 'coprime'}
            G = geom_targetN(type, N0, lambda, 8, Dtgt);
        case {'misc', 'imisc', 'emisc'}
            n = [0; cumsum(ies(:))];
            assert(numel(n) == N0, 'geom_u: the %s spacing set needs N0 - 1 = %d entries.', type, N0 - 1);
            u = max(1, round(Dtgt / ((max(n) - min(n)) * d0)));
            G = finish_geom(n * u * d0, u * d0, type, d0);
        case {'sma_u', 'sma_s'}
            % MOVABLE MIRROR ARRAY: positive half at u = (b_i + sigma) lambda/4,
            % negative half at -u, frames sigma = 0 .. F-1.
            nh = N0 / 2;  q4 = lambda / 4;
            if strcmpi(type, 'sma_u')
                b = 1 + 16 * (0:nh - 1);                   % uniform half, pitch 4 lambda
            else
                b = 1 + 4 * (0:nh - 1) .* (1:nh);          % increasing gaps 2, 4, 6 ... lambda
            end
            U = (b(:) + (0:F - 1)) * q4;                   % nh x F
            G = finish_geom([U(:, sig_d + 1); -U(:, sig_d + 1)], q4, type, d0);
            G.frames = [U; -U];                            % rows 1..nh: +u_i, rows nh+1..N0: -u_i
            G.b = b;  G.F = F;  G.sig_d = sig_d;
            % mirror pairs of the stacked frames: lag 2u = (b_i + sigma) lambda/2
            [ii, ss] = ndgrid(1:nh, 1:F);
            MI.iPos = (ss(:) - 1) * N0 + ii(:);
            MI.iNeg = (ss(:) - 1) * N0 + nh + ii(:);
            MI.q    = b(ii(:)).' + ss(:) - 1;
            MI.Q    = max(MI.q);
            MI.cnt  = accumarray(MI.q, 1, [MI.Q 1]);
            MI.nfft = 2 ^ nextpow2(8 * MI.Q);
            G.mir = MI;
            if any(MI.cnt == 0)
                warning('geom_u: %s mirror lags have %d holes (increase F).', type, sum(MI.cnt == 0));
            end
        otherwise
            error('geom_u: unknown array %s', type);
    end
end

function G = finish_geom(p, unit, type, d0)
    p = p(:);
    G.pos = p - (max(p) + min(p)) / 2;
    G.unit = unit;
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    G.sparsity = G.aperture / ((G.N - 1) * d0);
    G.min_gap = min(diff(sort(G.pos)));
end

function h = nf_channel(pos, d3, sE, bet, lambda)
% NUSW near-field channel (same model as array_stats) at arbitrary positions.
    [L, K] = size(d3);  pos = pos(:);  h = zeros(numel(pos), L, K);
    for l = 1:L
        for k = 1:K
            rn = sqrt(d3(l, k)^2 + pos.^2 - 2 * d3(l, k) * pos * sE(l, k));
            h(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
        end
    end
end

function sx = sma_mirror_angles(z, MI, nPk)
% Idea 1: z(+u) conj(z(-u)) ~ exp(j pi q sin(theta)) on the lag q = 2u/(lambda/2);
% the near-field range term is the same on both twins and cancels. Average the
% pairs per lag, take the spectrum of the hole-free virtual array, keep peaks.
    w = z(MI.iPos) .* conj(z(MI.iNeg));
    v = accumarray(MI.q, w, [MI.Q 1]) ./ max(MI.cnt, 1);
    sp = abs(fft([0; v], MI.nfft));                   % index k <-> sin(theta) = 2k/nfft
    pk = find(sp > circshift(sp, 1) & sp >= circshift(sp, -1));
    [~, o] = sort(sp(pk), 'descend');
    pk = pk(o(1:min(nPk, numel(o))));
    sx = mod(2 * (pk(:).' - 1) / MI.nfft + 1, 2) - 1;
end

function [Hhat, C, xhat, nfus, aErr] = csi_nf_estimate(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG, EST, useCape)
% Unified estimator. Far-field links: statistical LMMSE (data-block pilot).
% Near-field links: single-atom polar estimation on the estimation observation
% EST.Ze (static arrays: one pilot block, as Code T; movable array: the
% frames stacked into a synthetic aperture, each frame with 1/F of the energy). Candidates per link = FG.J refined
% dictionary atoms + extra angles (co-array MUSIC or mirror-pair angles) with a
% 1-D range search. Without CAPE each AP keeps its best-fitting candidate; with
% CAPE the candidates become positions and the cross-AP fused likelihood picks
% one (as csi_acquire_cape, same consistency gate). The channel at the data
% positions is rebuilt from the chosen (angle, range).
    LN = N * L;
    PDe = EST.PDe;  Ne = EST.Ne;  sw2e = EST.sw2e;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);  C = zeros(N, N, L, K);
    xhat = nan(K, nReal);  nfus = zeros(K, 1);  aErr = nan(L, K, nReal);
    PDw = PDe;  PDw.ds = PDe.ds / 2;  PDw.dq = PDe.dq / 2;
    gamL = log(size(PDe.W, 2)) + 0.58;
    nc = FG.J + EST.nExtra;
    for k = 1:K
        for l = 1:L
            if ~NFm(l, k)                               % far-field link: statistical LMMSE
                idx = (l - 1) * N + 1:l * N;
                yp = sqrt(p) * tau * H(idx, :, k) + sqrt(tau) * Np(:, :, l, k);
                Psi = p * tau * R(:, :, l, k) + eye(N);
                RPsi = R(:, :, l, k) / Psi;
                Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                C(:, :, l, k) = R(:, :, l, k) - p * tau * RPsi * R(:, :, l, k);
            end
        end
        nf = find(NFm(:, k));
        if isempty(nf), continue; end
        [~, ob] = sort(beta(nf, k), 'descend');
        candAP = nf(ob(1:min(FG.nCand, numel(nf))));
        fusAP  = nf(ob(1:min(FG.nAP, numel(nf))));
        if useCape, nfus(k) = numel(fusAP); end
        se2 = zeros(L, nReal);
        for m = 1:nReal
            Z = reshape(EST.Ze(:, m, :, k), Ne, L);
            cS = nan(L, nc);  cQ = nan(L, nc);
            bLoc = zeros(Ne, L);  mLoc = zeros(L, 1);  sLoc = zeros(L, 1);  qLoc = zeros(L, 1);
            for l = nf(:)'
                zl = Z(:, l);
                [~, oc] = sort(abs(PDe.W' * zl), 'descend');
                for j = 1:FG.J
                    [bj, mj, cS(l, j), cQ(l, j)] = refine_atom(zl, PDe.s(oc(j)), PDe.q(oc(j)), PDe);
                    if mj > mLoc(l), mLoc(l) = mj;  bLoc(:, l) = bj;  sLoc(l) = cS(l, j);  qLoc(l) = cQ(l, j); end
                end
                if strcmp(EST.kind, 'coarray')
                    sx = EST.angCA{l};
                elseif strcmp(EST.kind, 'sma')
                    sx = sma_mirror_angles(zl, EST.mir, EST.nPkSMA);
                else
                    sx = [];
                end
                for j = 1:min(numel(sx), EST.nExtra)
                    q0 = coarray_range_search(zl, sx(j), PDe, FG.nRng);
                    [bj, mj, cS(l, FG.J + j), cQ(l, FG.J + j)] = refine_atom(zl, sx(j), q0, PDe);
                    if mj > mLoc(l)
                        mLoc(l) = mj;  bLoc(:, l) = bj;  sLoc(l) = cS(l, FG.J + j);  qLoc(l) = cQ(l, FG.J + j);
                    end
                end
            end
            if useCape
                cand = [];
                for l = candAP(:)'
                    for j = 1:nc
                        if isfinite(cS(l, j))
                            cand = [cand; cand_positions(FG.AP(l), cS(l, j), cQ(l, j), FG.sL, FG.hD)];  %#ok<AGROW>
                        end
                    end
                end
                if isempty(cand), cand = FG.AP(candAP(1)) + 1; end
                Zf = Z(:, fusAP);
                sc = fused_score(cand, Zf, fusAP, FG, PDe);
                [best, jb] = max(sc);
                x = cand(jb);  st = FG.step0;  nbr = exp(1j * pi / 4 * (0:7)).';
                for it = 1:FG.nIter
                    xn = x + st * nbr;
                    sn = fused_score(xn, Zf, fusAP, FG, PDe);
                    [mx, jn] = max(sn);
                    if mx > best
                        best = mx;  x = xn(jn);
                    else
                        st = st / 2;
                    end
                end
                xhat(k, m) = mod(real(x), FG.sL) + 1j * mod(imag(x), FG.sL);
                [d3x, sxp] = geo_from_pos(FG.APw, x, FG.hD);
            end
            for l = nf(:)'
                zl = Z(:, l);
                b = bLoc(:, l);  s1 = sLoc(l);  q1 = qLoc(l);
                if useCape
                    [bf, mfu, sf, qf] = refine_atom(zl, sxp(l), 1 / d3x(l), PDw);   % fused lobe
                    sameLobe = abs(bLoc(:, l)' * bf) >= FG.consist;
                    clearWin = mLoc(l) ^ 2 - mfu ^ 2 > gamL * sw2e;
                    if ~((sameLobe && mLoc(l) > mfu) || clearWin)
                        b = bf;  s1 = sf;  q1 = qf;
                    end
                end
                c = b' * zl;
                g = Ne * beta(l, k) / (Ne * beta(l, k) + sw2e);
                idx = (l - 1) * N + 1:l * N;
                Hhat(idx, m, k) = polar_atoms(PD.pos, s1, q1, PD.lambda) * (g * c * sqrt(N / Ne));
                res = zl - b * c;
                off = max(real(res' * res) - (Ne - 1) * sw2e, 0);
                se2(l, m) = (Ne * beta(l, k) * sw2e / (Ne * beta(l, k) + sw2e) + off) / Ne;
                dS = abs(s1 - EST.sTrue(l, k));
                aErr(l, k, m) = min(dS, 2 - dS);
            end
        end
        for l = nf(:)'
            C(:, :, l, k) = mean(se2(l, :)) * eye(N);
        end
    end
end

function O = eval_scheme_u(ST, bet, useCape, PD, FG, p, EP, UEm, seed, ES)
% One scheme at one SNR: channels, pilots of all blocks / frames, estimation,
% clustering and the four detectors (estimated CSI).
    [L, K] = size(bet);  Na = ST.N;
    rng(seed);
    Ht = draw_channel(ST.hD, ST.Rs2, ST.NF, Na, L, K, EP.nReal);
    tE = tic;
    sw2 = 1 / (p * K);                                   % de-spread pilot noise of one block (Code T)
    if strcmp(ES.kind, 'sma')
        sw2e = sw2 * ES.F;                               % F frames, each with 1/F of the pilot energy
    else
        sw2e = sw2;                                      % static: one pilot block, as Code T
    end
    Ze = repmat(reshape(ES.he, ES.Ne, 1, L, K), [1, EP.nReal, 1, 1]) + ...
         sqrt(sw2e / 2) * (randn(ES.Ne, EP.nReal, L, K) + 1j * randn(ES.Ne, EP.nReal, L, K));
    EST = struct('kind', ES.kind, 'Ze', Ze, 'Ne', ES.Ne, 'sw2e', sw2e, 'PDe', ES.PDe, ...
                 'sTrue', ES.sTrue, 'nPkSMA', 2, 'nExtra', 0);
    EST.mir = ES.mir;  EST.angCA = cell(L, 1);
    if strcmp(ES.kind, 'coarray')
        rsMain = rng;  rng(mod(seed + 5e8, 2^32 - 1));     % covariance noise from a side stream (paired)
        Phi = exp(-1j * 2 * pi * (0:K - 1)' * (0:K - 1) / K) / sqrt(K);
        for l = 1:L
            idx = (l - 1) * Na + 1:l * Na;
            Hl = reshape(Ht(idx, 1, :), [Na, K]);
            Yl = sqrt(p) * Hl * Phi + sqrt(0.5) * (randn(Na, K) + 1j * randn(Na, K));   % one block, as Code T
            EST.angCA{l} = coarray_music_angles((Yl * Yl') / K, ES.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
        end
        rng(rsMain);
        EST.nExtra = FG.nPk;
    elseif strcmp(ES.kind, 'sma')
        EST.nExtra = EST.nPkSMA;
    end
    [Hh, Cs, xh, nfus, aErr] = csi_nf_estimate(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG, EST, useCape);
    O.time = toc(tE);
    O.nfus = mean(nfus);
    O.loc = [];
    if useCape
        dx = abs(real(xh) - real(UEm(:)));  dx = min(dx, FG.sL - dx);
        dy = abs(imag(xh) - imag(UEm(:)));  dy = min(dy, FG.sL - dy);
        e = sqrt(dx.^2 + dy.^2);
        O.loc = e(isfinite(e));
    end
    a = aErr(isfinite(aErr));
    O.hit = NaN;  O.rmse = NaN;
    if ~isempty(a)
        hit = a < EP.tolS;
        O.hit = mean(hit);
        if any(hit), O.rmse = sqrt(mean(a(hit) .^ 2)); end
    end
    O.nmse = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);
    Drt = cluster_rate_rule(ST.R2, Cs, p, ST.Dcn);
    bE = zeros(4, K);  nE = zeros(4, K);  nb = 0;  ne = 0;
    for mc = 1:EP.nReal
        Hm = reshape(Ht(:, mc, :), [L * Na, K]);
        He = reshape(Hh(:, mc, :), [L * Na, K]);
        [bu, nu, bpu, epu] = metrics_case(He, Hm, Cs, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);
        bE = bE + bu;  nE = nE + nu;  nb = nb + bpu;  ne = ne + epu;
    end
    O.be  = sum(bE, 2) / max(K * nb, 1);
    seU   = EP.prelog * log2(1 + 1 ./ max(nE / max(ne, 1), 1 / (EP.nReal * EP.nSym)));
    O.se  = seU;
    O.sum = sum(seU, 2);
end

function [O, nf] = run_task_u(D, G, PD, ES, useCape, Kc, p, seed, CF)
% One task: one drop, one scheme, one SNR.
    ST = array_stats(G, D.d3, D.sE, D.phL, D.thL, D.bet, CF.lambda, CF.eps_NF, CF.asd_phi, CF.asd_th, CF.clu_div);
    nf = mean(ST.NF(:));
    FG = cape_cfg(true, D.APp, D.APw, CF.sL, CF.hD, CF.nAP, CF.nCand, CF.J, CF.step0, CF.nIter, CF.consist);
    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    if strcmp(ES.kind, 'sma')
        ES.he = nf_channel(G.frames(:), D.d3, D.sE, D.bet, CF.lambda) .* reshape(ST.NF, [1, size(ST.NF)]);
    else
        ES.he = ST.hD;
    end
    ES.sTrue = D.sE;
    EP = struct('nReal', CF.nReal, 'nSym', CF.nSym, 'M', CF.M, 'dth', CF.dth, 'maxBr', CF.maxBr, ...
                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF, 'tolS', CF.tolS);
    O = eval_scheme_u(ST, D.bet, useCape, PD, FG, p, EP, D.UEm, seed, ES);
end
"""
src = src.rstrip() + '\n' + UFUN
src = src.replace('    chunk = 256;', '    chunk = 64;   % symbol block of the list detector (memory only; results do not change)')
for bad in ['codeT_arraycomp_results', '[CODE T] save BEFORE']:
    if bad in src:
        sys.exit('leftover: ' + bad)
(REPO / 'Code_U_MovableSparse.m').write_text(src)
plot = (r"""%% PLOT_CODEU_FIGURES  Redraw Code U figures U1-U13 without rerunning the simulation.
%  Uses the workspace variables after Code_U_MovableSparse.m, or loads
%  codeU_movable_results.mat if they are not there. It never clears anything.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('bs_kind', 'var') || ~exist('EN_hit', 'var')
    load('codeU_movable_results.mat');
end
""" + fig)
(REPO / 'plot_codeU_figures.m').write_text(plot)
print('Code_U lines:', src.count('\n'), ' plot lines:', plot.count('\n'))

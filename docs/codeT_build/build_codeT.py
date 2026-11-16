#!/usr/bin/env python3
"""Assemble Code_T_ArrayComparison.m from Code_S_IdeaUniform.m (which ran cleanly).

Every edit is anchored on exact text that must occur exactly once."""
import sys, pathlib

REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent
src = (REPO / 'Code_S_IdeaUniform.m').read_text()
fig = (HERE / 'figblock.m').read_text()


def once(text, anchor):
    n = text.count(anchor)
    if n != 1:
        sys.exit(f'anchor found {n} times (need 1):\n{anchor[:200]}')


def replace_once(text, old, new):
    once(text, old)
    return text.replace(old, new)


def cut_between(text, start, end_anchor, new, keep_end=True):
    """Replace text[start : end_anchor] (end anchor kept unless keep_end=False)."""
    once(text, start)
    once(text, end_anchor)
    i = text.index(start)
    j = text.index(end_anchor)
    if j <= i:
        sys.exit('end anchor before start anchor')
    if not keep_end:
        j += len(end_anchor)
    return text[:i] + new + text[j:]


# ---------------------------------------------------------------- 1. header
HEADER = r"""clear; close all; clc;
%% ====================================================================
%  CODE T : ARRAY COMPARISON IN A CROWDED CELL-FREE SYSTEM (from Code R)
%  --------------------------------------------------------------------
%  Crowded (loaded) cell-free uplink: N = 8 antennas per AP, K = 16 users,
%  L = 8 APs, CAPE channel estimation, 16-QAM. Compares:
%    1. ULA, d = lambda/2            dense reference (aperture 0.35 m)
%    2. ULA, d = 4 lambda            wide-spaced ULA, aperture 2.8 m
%    3. Nested array                 stretched to the same aperture (~2.8 m)
%    4. Nested array + co-array      SAME antennas as 3, co-array processing
%    5. Coprime array                stretched to the same aperture (~3.0 m)
%    6. Coprime array + co-array     SAME antennas as 5, co-array processing
%  SNR axis = received SNR per antenna at each user's strongest AP (uplink
%  power control: every user scales its power so that this SNR equals the
%  axis value). User clusters keep APs within 10 dB of the strongest AP.
%  "Difference co-array" is a PROCESSING option of a sparse array, not a
%  separate array: schemes 3/4 and 5/6 use identical antennas. Every scheme
%  sees the same random drops, channels, pilots and noise in a setup
%  (paired comparison), so a 3-vs-4 or 5-vs-6 difference is caused only by
%  the co-array processing. The co-array MUSIC subarray is NOT capped here,
%  so the co-array gets its full smoothed virtual aperture.
%  Figures T1..T10. Results are saved to codeT_arraycomp_results.mat BEFORE
%  plotting; plot_codeT_figures.m redraws the figures without rerunning.
%% ====================================================================
"""
src = cut_between(src, 'clear; close all; clc;\n%% ====', '\nclear;\nclose all;\nclc;\n', HEADER)

# ---------------------------------------------------------- 2. array set
ARRAYSET = r"""%% ====================================================================
%  [CODE T] ARRAY SET (equal hardware: N0 = 32 physical antennas per AP)
%  --------------------------------------------------------------------
%  'ula'     dense ULA, spacing d = lambda/2 (reference)
%  'ulaw'    wide-spaced ULA, spacing d = ulaw_d_lambda * lambda. Uniform
%            spacing, so it is a ULA; wide spacing enlarges the aperture
%            with the same N but creates grating lobes (d > lambda/2).
%  'nested'  two-level nested array (Pal & Vaidyanathan, IEEE TSP 2010)
%  'coprime' extended coprime array (Pal & Vaidyanathan, DSP/SPE 2011):
%            N1 elements at spacing M*d0 plus 2M elements at spacing N1*d0,
%            M < N1 coprime, 2M + N1 - 1 = N0 (M = 8, N1 = 17 for N0 = 32)
%  Co-array processing (en_idea = 1) is only meaningful for nested/coprime:
%  a wide ULA has no hole-free co-array (all its lags are multiples of eta).
%% ====================================================================
ulaw_d_lambda = 4;                    % wide-ULA spacing in wavelengths
ulaw_eta      = 2 * ulaw_d_lambda;    % same spacing in units of d0 = lambda/2
ca_Nsub_max   = 1000;                 % no cap: full smoothed virtual aperture for the co-array
sp_equal_aperture = true;             % stretch nested/coprime to the wide-ULA aperture
                                      % (false -> classic lambda/2 grid, ~1 m, almost no near field)

% ---- [CODE T] crowded cell-free setting (few antennas per AP, many users) ----
modOrder = 16;                        % 16-QAM (QPSK reaches zero errors by 4-8 dB: no visible waterfall)
clu_div  = 10;                        % user cluster: APs within 10 dB of the strongest AP (was 3 dB)
SNR_dB   = 0:4:28;                    % received SNR per antenna at the strongest AP [dB]
SNR_lin  = 10.^(SNR_dB / 10);
nSNR     = numel(SNR_dB);
[~, snr_ref_idx] = min(abs(SNR_dB - 10));
snr_ref_dB = SNR_dB(snr_ref_idx);

"""
src = cut_between(src, '%% ====================================================================\n%  [CODE S] UNIFORM SPARSE ARRAY DESIGN',
                  'eqN_study   = true;\nen_names', ARRAYSET)

# ------------------------------------------------------ 3. scheme set
SCHEMES = r"""eqN_study   = true;
en_names    = {'ULA lambda/2 + CAPE', sprintf('ULA %g lambda + CAPE', ulaw_d_lambda), 'Nested + CAPE', ...
               'Nested + co-array + CAPE', 'Coprime + CAPE', 'Coprime + co-array + CAPE'};
en_leg      = {'ULA (\lambda/2)', sprintf('ULA (%g\\lambda)', ulaw_d_lambda), 'Nested', ...
               'Nested + co-array', 'Coprime', 'Coprime + co-array'};
if sp_equal_aperture
    en_leg(3:6) = {'Nested (stretched)', 'Nested (stretched) + co-array', ...
                   'Coprime (stretched)', 'Coprime (stretched) + co-array'};
end
en_geom     = {'ula', 'ulaw', 'nested', 'nested', 'coprime', 'coprime'};
en_idea     = [0 0 0 1 0 1];   % 1 = add difference co-array processing (same antennas)
en_cape     = [1 1 1 1 1 1];   % CAPE for every scheme: the array is the only variable
en_arr      = [1 2 3 3 4 4];   % physical array of each scheme (1 ULA, 2 wide ULA, 3 nested, 4 coprime)
nSchE       = numel(en_names);
en_caseName = {'Crowded (N = 8, K = 16)'};
en_N0       = 8;               % 8 PHYSICAL antennas per AP
en_K        = 16;              % 16 users (K > N: loaded cell-free system)
nCaseE      = numel(en_N0);
nSetE       = 10;              % random AP/user drops (>= 10 for the paper; 4 for a quick look)
nRealE      = 4;               % channel realisations per drop
nSymE       = 800;             % symbols per realisation: 10 x 4 x 800 x 16 users x 4 bits
                               % = 2.0e6 bits per point, enough to show BER ~1e-5
if sp_equal_aperture           % target aperture = aperture of the wide ULA
    sp_Dtgt = (en_N0(1) - 1) * ulaw_eta * lambda / 2;
else
    sp_Dtgt = [];
end
"""
src = cut_between(src, 'eqN_study   = true;\nen_names', 'nSymE       = 100;\n', SCHEMES, keep_end=False)

# ------------------------------------------- 4. geometry call + paired seed
src = replace_once(src, 'Gtmp = geom_targetN(en_geom{scE}, en_N0(caE), lambda, uspa_eta);',
                   'Gtmp = geom_targetN(en_geom{scE}, en_N0(caE), lambda, ulaw_eta, sp_Dtgt);')
src = replace_once(src, "        EN_CA{scE, caE}  = coarray_setup(Gtmp.pos, lambda, 'diff');",
                   "        EN_CA{scE, caE}  = coarray_setup(Gtmp.pos, lambda, 'diff', Gtmp.unit);   % co-array on the array's own grid")
src = replace_once(src, "            EN_CAs{scE, caE} = coarray_setup(Gtmp.pos, lambda, 'sum');",
                   "            EN_CAs{scE, caE} = coarray_setup(Gtmp.pos, lambda, 'sum', Gtmp.unit);")
src = replace_once(src, 'seed = 9e6 + 1e4 * ns + 1000 * caE + 10 * scE + si;',
                   'seed = 9e6 + 1e4 * ns + 1000 * caE + si;   % [CODE T] same seed for every scheme (paired)')
src = replace_once(src, "fprintf('\\n=== [CODE I] Array + detector comparison: %d schemes x %d cases ===\\n', nSchE, nCaseE);",
                   "fprintf('\\n=== [CODE T] array comparison: %d schemes x %d cases x %d setups ===\\n', nSchE, nCaseE, nSetE);")

# ----------------------------------------------------------- 5. K-sweep
KSWEEP = r"""% ---- [CODE T] user-count sweep at the reference SNR, all six schemes ----
% N is fixed at en_N0(1) = 32, so the prebuilt geometry, polar dictionary and
% co-array of case 1 are reused (identical to rebuilding them).
Ksw     = [8 12 16 20 24];       % all K >= N = 8: crowded
Ksw_sch = 1:nSchE;
nKsw    = numel(Ksw);  nKsch = numel(Ksw_sch);
EN_Kbe  = zeros(nKsch, nKsw, nSetE);    % BER, Cross-AP detector
EN_Ksum = zeros(nKsch, nKsw, nSetE);    % sum SE, Cross-AP detector
if eqN_study
fprintf('\n=== [CODE T] user-count sweep at %d dB: K = %s ===\n', SNR_dB(snr_ref_idx), mat2str(Ksw));
parfor ns = 1:nSetE
    kbe = zeros(nKsch, nKsw);  ksum = zeros(nKsch, nKsw);
    for ik = 1:nKsw
        Kc  = Ksw(ik);
        APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
        UEp = drop_scenario(APp, Kc, L, min(n_collinear_pairs, floor(Kc / 2)), ...
                            r_min_NF, r_max_en, 'hotspot', hs_r, squareLen);
        UEm = mod(real(UEp), squareLen) + 1j * mod(imag(UEp), squareLen);
        APw = wrap_images(APp, squareLen, wrap_fix);
        [d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
        EPc = struct('nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, ...
                     'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);
        for jj = 1:nKsch
            scE = Ksw_sch(jj);
            Gc  = EN_geo{scE, 1};
            PDc = EN_PD{scE, 1};
            ST  = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);
            FG  = cape_cfg(true, APp, APw, squareLen, hDiff, cape_nAP, cape_nCand, ...
                           cape_J, cape_step0, cape_nIter, cape_consist);
            FG.useCA = (en_idea(scE) > 0);  FG.co_idea = en_idea(scE);
            FG.CA = EN_CA{scE, 1};  FG.CAs = [];
            FG.nAng = ca_nAngGrid;  FG.nRng = ca_nRangeGrid;
            FG.nPk = max(ca_nPeak, ceil(Kc / 2));  FG.nSubMax = ca_Nsub_max;
            seed = 3e7 + 1e4 * ns + 100 * ik;   % same seed for every scheme (paired)
            O = eval_scheme(ST, bet, en_cape(scE), PDc, FG, SNR_lin(snr_ref_idx), EPc, false, UEm, seed);
            kbe(jj, ik)  = O.be(4);
            ksum(jj, ik) = O.sum(4);
        end
    end
    EN_Kbe(:, :, ns) = kbe;  EN_Ksum(:, :, ns) = ksum;
    fprintf('  user-count sweep setup %d/%d done\n', ns, nSetE);
end
end

"""
src = cut_between(src, '% ---- [CODE R] K-sweep at the reference SNR',
                  'if eqN_study\n    %% ============ [CODE R] FOCUSED SPARSE-ARRAY COMPARISON FIGURES', KSWEEP)

# ----------------- 5b. power control and per-K pilot overhead (both loops)
LG = "[d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);"
assert src.count(LG) == 2, src.count(LG)
src = src.replace(LG, LG + "\n        bet = bet ./ max(bet, [], 1);   % [CODE T] uplink power control: received SNR at each\n"
                  "                                         % user's strongest AP = the SNR axis value")
PL = "'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);"
assert src.count(PL) == 2, src.count(PL)
src = src.replace(PL, "'maxBr', maxBr, 'mod', modOrder, 'prelog', 1 - Kc / tau_c, 'eps_NF', eps_NF);   % pilots: tau_p = Kc")

# --------------------------------------------- 6. save first, then figures
fig_ind = ''.join(('    ' + ln if ln.strip() else ln) for ln in fig.splitlines(True))
FIGS = r"""if eqN_study
    % [CODE T] save BEFORE plotting, so a plotting problem can never lose the run
    save('codeT_arraycomp_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', ...
         'en_names', 'en_leg', 'en_geom', 'en_idea', 'en_cape', 'en_arr', 'en_caseName', ...
         'en_N0', 'en_K', 'ulaw_d_lambda', 'sp_equal_aperture', 'modOrder', 'clu_div', ...
         'nSetE', 'nRealE', 'nSymE', 'EN_geo', ...
         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...
         'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\nResults saved to codeT_arraycomp_results.mat\n');

""" + fig_ind + "end\n\n"
src = cut_between(src, 'if eqN_study\n    %% ============ [CODE R] FOCUSED SPARSE-ARRAY COMPARISON FIGURES',
                  'function [Hhat, C, Hpre, Cpre] = csi_acquire(', FIGS)

# ------------------------- 7. co-array covariance noise on its own stream
src = replace_once(src,
    "        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);\n        for l = 1:L\n",
    "        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);\n"
    "        useSide = isfield(FG, 'caSeed');      % [CODE T] covariance noise from a separate\n"
    "        if useSide, rsMain = rng; rng(FG.caSeed); end   % stream, so the rest of the run\n"
    "        for l = 1:L                           % stays paired with the no-co-array scheme\n")
src = replace_once(src,
    "            end\n        end\n    end\n    PDw = PD;",
    "            end\n        end\n        if useSide, rng(rsMain); end          % [CODE T] back to the main stream\n    end\n    PDw = PD;")
# eval_scheme hands the side-stream seed to the co-array branch
src = replace_once(src, "    rng(seed);\n    Ht = draw_channel(",
                   "    rng(seed);\n    FG.caSeed = mod(seed + 5e8, 2^32 - 1);   % [CODE T] side stream for co-array covariance noise\n"
                   "    Ht = draw_channel(")

# ------------------------------------------------ 8. geometry: ulaw + coprime
src = replace_once(src, "% 'usparse' is a uniform wide-spacing array (new idea, scheme 4).",
                   "% [CODE T] 'ulaw' is a wide-spaced ULA (spacing eta*d0); 'coprime' is the\n"
                   "% extended coprime array of Pal & Vaidyanathan.")
src = replace_once(src, "function G = geom_targetN(type, N0, lambda, eta)\n",
                   "function G = geom_targetN(type, N0, lambda, eta, Dtgt)\n")
src = replace_once(src, "    if nargin < 4 || isempty(eta), eta = 9; end\n    d0 = lambda / 2;\n",
                   "    if nargin < 4 || isempty(eta), eta = 8; end\n    d0 = lambda / 2;\n"
                   "    u  = 1;                       % [CODE T] grid unit in half-wavelengths\n")
_gi = src.index("function G = geom_targetN(")
_ge = src.index("\nend\n", _gi) + 5
src = src[:_gi] + replace_once(src[_gi:_ge], "    end\n    p = n * d0;\n    G.pos = p - (max(p) + min(p)) / 2;\n",
                   "    end\n"
                   "    % [CODE T] stretch nested/coprime onto a coarser grid so that their aperture\n"
                   "    % matches Dtgt (the wide-ULA aperture). All positions then lie on a grid of\n"
                   "    % u half-wavelengths, so for u > 1 the array also has grating lobes.\n"
                   "    if nargin >= 5 && ~isempty(Dtgt) && Dtgt > 0 && any(strcmpi(type, {'nested', 'coprime'}))\n"
                   "        u = max(1, round(Dtgt / ((max(n) - min(n)) * d0)));\n"
                   "    end\n"
                   "    p = n * u * d0;\n    G.pos = p - (max(p) + min(p)) / 2;\n"
                   "    G.unit = u * d0;              % grid unit [m] (co-array lags are multiples of it)\n") + src[_ge:]
src = replace_once(src, "function CA = coarray_setup(pos, lambda, mode)\n",
                   "function CA = coarray_setup(pos, lambda, mode, unit)\n")
src = replace_once(src, "    if nargin < 3, mode = 'diff'; end\n    d0 = lambda / 2;\n",
                   "    if nargin < 3, mode = 'diff'; end\n"
                   "    if nargin < 4 || isempty(unit), unit = lambda / 2; end   % [CODE T] grid of the array\n"
                   "    d0 = unit;\n")
GEOM = r"""        case 'ulaw'
            % [CODE T] wide-spaced ULA: uniform spacing eta*d0 (eta = 8 ->
            % d = 4 lambda). Still a uniform linear array; the wide spacing
            % enlarges the aperture at the same N and creates grating lobes.
            n = (0:N0 - 1)';  u = eta;
        case 'coprime'
            % [CODE T] extended coprime array (Pal & Vaidyanathan 2011): Nc
            % elements at spacing Mc*d0 and 2*Mc elements at spacing Nc*d0,
            % Mc < Nc coprime, sharing the origin -> 2*Mc + Nc - 1 = N0. Among
            % the valid pairs the one with the longest hole-free co-array
            % (Mc*Nc + Mc - 1 lags each side) is used: (2, 5) for N0 = 8,
            % (8, 17) for N0 = 32.
            best = -1;  Mc = 0;  Nc = 0;
            for Mt = 2:floor(N0 / 2)
                Nt = N0 + 1 - 2 * Mt;
                if Nt > Mt && gcd(Mt, Nt) == 1 && Mt * Nt + Mt - 1 > best
                    best = Mt * Nt + Mt - 1;  Mc = Mt;  Nc = Nt;
                end
            end
            if Mc == 0, error('geom_targetN: no coprime pair for N0 = %d', N0); end
            n = unique([(0:Nc - 1) * Mc, (0:2 * Mc - 1) * Nc])';
"""
src = cut_between(src, "        case 'usparse'\n", "        otherwise\n            error('geom_targetN: unknown type %s', type);", GEOM)

# ----------------------------------- 9. drop the Code S design functions
tail = 'function [eta_max_snr, Pf] = ambiguity_budget('
once(src, tail)
src = src[:src.index(tail)].rstrip() + '\n'

# ------------------------------------------------------------ sanity checks
for bad in ['uspa_eta', 'usd_', 'eta_max', 'Pf_eta', 'ambiguity_budget', 'optimise_positions',
            'steer_nf', 'psl_of', "'usparse'", 'codeS_ideaUniform']:
    if bad in src:
        sys.exit(f'leftover reference: {bad}')

(REPO / 'Code_T_ArrayComparison.m').write_text(src)
plot = (r"""%% PLOT_CODET_FIGURES  Redraw Code T figures T1-T10 without rerunning the simulation.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_geo', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
""" + fig)
(REPO / 'plot_codeT_figures.m').write_text(plot)
print('Code_T lines:', src.count('\n'), ' plot lines:', plot.count('\n'))

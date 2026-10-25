#!/usr/bin/env python3
"""Build Code E from Code D: add the two-level nested array, the difference
co-array + vectorization of Pinto & de Lamare (arXiv:2312.06033), and a
comparison study with five new figures (CO-1..CO-5). Code D's own figures and
pipeline are left untouched."""
import re, os

src = open('/home/user/csc/Code_D_Sparse_CAPE.m').read()


def rep(old, new, count=1):
    global src
    n = src.count(old)
    assert n == count, (n, old[:80])
    src = src.replace(old, new)


# ---------------------------------------------------------------- 1. header
rep("""%  CODE D : SPARSE APs MADE USEFUL BY CROSS-AP POSITION-FUSED ESTIMATION""",
    """%  CODE E : CODE D + DIFFERENCE CO-ARRAY AND VECTORIZATION (VIRTUAL ARRAY)
%  ------------------------------------------------------------------
%  Everything in Code D is kept unchanged. What is added here:
%   * 'nested' array geometry: the two-level nested array of Pal & Vaidyanathan
%     (IEEE TSP 2010), the array used by Pinto & de Lamare (arXiv:2312.06033)
%     for super-resolution-array communications.
%   * Difference co-array + vectorization: the sample covariance of the raw
%     pilot snapshots is vectorized, duplicate lags are averaged, the
%     contiguous virtual ULA is extracted and spatially smoothed, and MUSIC on
%     that VIRTUAL array supplies high-resolution angle candidates. The range
%     is then found by a 1-D search at that angle on the physical observation,
%     which is the procedure used in the near-field co-array literature.
%   * co_study: compares ECA + CAPE (proposed) against the nested array with
%     and without co-array processing, with five new figures CO-1 ... CO-5.
%
%  HONEST NOTES (please read before using the numbers)
%   1. A virtual sensor adds DEGREES OF FREEDOM and RESOLUTION, not received
%      energy. It cannot raise the detector's array gain by itself.
%   2. The standard difference co-array assumes PLANE waves. On near-field
%      links the curved wavefront breaks the clean virtual-ULA structure, so
%      the co-array angles here are approximate candidates that CAPE and the
%      local refinement then correct.
%   3. The covariance must be estimated from only tau = K pilot snapshots, so
%      it is noisy; this limits the achievable co-array gain.
%   4. With K = 8 users and N = 63 antennas per AP the O(N^2) degree-of-freedom
%      advantage is NOT the binding constraint here; the benefit to look for is
%      resolution (and antenna count), not source count.
%
%  ORIGINAL CODE D HEADER FOLLOWS.
%  ==================================================================
%  CODE D : SPARSE APs MADE USEFUL BY CROSS-AP POSITION-FUSED ESTIMATION""")

# ---------------------------------------------------------------- 2. params
rep("""wrap_fix   = true;   % full 3 x 3 wrap-around (see header)
""",
    """wrap_fix   = true;   % full 3 x 3 wrap-around (see header)

%% ====================================================================
%  [CODE E] NESTED ARRAY AND DIFFERENCE CO-ARRAY (VIRTUAL ARRAY)
%  ------------------------------------------------------------------
%  Two-level nested array (Pal & Vaidyanathan, IEEE TSP 2010): an inner ULA of
%  nest_N1 elements at spacing d0 = lambda/2 and an outer ULA of nest_N2
%  elements at spacing (nest_N1 + 1) d0. Its difference co-array is a HOLE-FREE
%  virtual ULA with 2*nest_N2*(nest_N1+1) - 1 virtual sensors, i.e. about
%  N^2/2 virtual sensors from N = nest_N1 + nest_N2 physical ones.
%  The default (11, 12) gives 23 antennas over 7.2 m: about the SAME aperture
%  as the 63-element ECA with roughly one third of the antennas, which is the
%  antenna-saving comparison made by Pinto & de Lamare (arXiv:2312.06033).
%% ====================================================================
nest_N1      = 11;       % inner (dense) elements, spacing d0
nest_N2      = 12;       % outer (sparse) elements, spacing (nest_N1+1) d0
ca_use       = true;     % use the difference co-array to propose candidates
ca_nAngGrid  = 721;      % angle grid of the co-array MUSIC spectrum
ca_nRangeGrid = 60;      % inverse-range grid of the 1-D range search
ca_nPeak     = 3;        % co-array angle candidates kept per proposing AP
ca_Nsub_max  = 80;       % cap on the spatial-smoothing subarray length (speed)
""")

# ---------------------------------------------------------------- 3. geometry
rep("""        case {'eca', 'eca_sa'}""",
    """        case 'nested'
            % [CODE E] two-level nested array (Pal & Vaidyanathan, TSP 2010):
            % inner ULA at spacing d0 plus outer ULA at spacing (N1+1) d0.
            assert(nargin >= 10 && ~isempty(nN1) && ~isempty(nN2), ...
                   'array_geometry: nested needs nest_N1 and nest_N2.');
            n = unique([(1:nN1), (1:nN2) * (nN1 + 1)])';
        case {'eca', 'eca_sa'}""")
rep("function G = array_geometry(type, N_ref, lambda, eta, M, Nc, Lc, Lc_sa)",
    "function G = array_geometry(type, N_ref, lambda, eta, M, Nc, Lc, Lc_sa, nN1, nN2)")
rep("""%  'eca_sa' the same ECA with Lc_sa basic arrays (about the dense aperture)""",
    """%  'eca_sa' the same ECA with Lc_sa basic arrays (about the dense aperture)
%  'nested' [CODE E] two-level nested array: inner nN1 elements at spacing d0
%           plus outer nN2 elements at spacing (nN1+1) d0 (Pal & Vaidyanathan,
%           IEEE TSP 2010). Its difference co-array is a hole-free virtual ULA.""")

# ---------------------------------------------------------------- 4. CAPE: co-array candidates
rep("""            % ---- 1. candidate positions from the strongest APs ----
            cand = [];
            for l = candAP(:)'
                for j = 1:FG.J
                    cand = [cand; cand_positions(FG.AP(l), sLoc(l, j), qLoc(l, j), FG.sL, FG.hD)];  %#ok<AGROW>
                end
            end""",
    """            % ---- 1. candidate positions from the strongest APs ----
            cand = [];
            for l = candAP(:)'
                for j = 1:FG.J
                    cand = [cand; cand_positions(FG.AP(l), sLoc(l, j), qLoc(l, j), FG.sL, FG.hD)];  %#ok<AGROW>
                end
            end
            % ---- 1b. [CODE E] extra candidates from the DIFFERENCE CO-ARRAY ----
            % The virtual array (vectorized covariance) gives high-resolution
            % angles; the range is then found by a 1-D search at that angle on
            % the physical observation, as in the near-field co-array
            % literature. These only ADD candidates: the fused likelihood of
            % step 2 still decides which one is the true position.
            if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)
                for l = candAP(:)'
                    sMu = FG.caAng{l};
                    for jj = 1:numel(sMu)
                        qj = coarray_range_search(Z(:, m, l), sMu(jj), PD, FG.nRng);
                        cand = [cand; cand_positions(FG.AP(l), sMu(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>
                    end
                end
            end""")

# covariance of the raw pilot snapshots, built once per call
rep("""    sw2 = 1 / (p * tau);
    PDw = PD;  PDw.ds = PD.ds / 2;  PDw.dq = PD.dq / 2;   % local window around the fused atom""",
    """    sw2 = 1 / (p * tau);
    % [CODE E] sample covariance of the RAW pilot snapshots at every AP (all K
    % users present, orthogonal pilots). Only used when the co-array is on.
    if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)
        FG.caAng = cell(L, 1);
        % orthogonal pilot matrix (explicit DFT: no toolbox dependency)
        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            Hl = reshape(H(idx, 1, :), [N, K]);
            Yl = sqrt(p) * Hl * Phi(1:K, :) ...
                 + sqrt(0.5) * (randn(N, tau) + 1j * randn(N, tau));
            Rl = (Yl * Yl') / tau;                     % sample covariance, all K users
            % The MUSIC spectrum of the virtual array is a PER-AP quantity: it
            % already contains the peaks of every user, so it is computed ONCE
            % per AP here and reused for all users (a large speed saving).
            FG.caAng{l} = coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
        end
    end
    PDw = PD;  PDw.ds = PD.ds / 2;  PDw.dq = PD.dq / 2;   % local window around the fused atom""")

# ---------------------------------------------------------------- 5. new functions
CA_FUNCS = r"""

%% ====================================================================
%  [CODE E] DIFFERENCE CO-ARRAY AND VECTORIZATION (VIRTUAL ARRAY)
%  Pal & Vaidyanathan (IEEE TSP 2010); used for communications by
%  Pinto & de Lamare (arXiv:2312.06033).
%% ====================================================================
function CA = coarray_setup(pos, lambda)
% Difference co-array of a linear array whose elements sit on the lambda/2
% grid. Every pair (m, n) of physical elements contributes one VIRTUAL sensor
% at the lag x_m - x_n. The contiguous run of lags around zero is the virtual
% ULA that can be used by subspace methods.
    d0 = lambda / 2;
    pos = pos(:);
    % Lags MUST come from position DIFFERENCES: the array is centred, so
    % pos/d0 can be half-integer and rounding it would corrupt the lag map,
    % while differences of centred positions are exact multiples of d0.
    Dm = round((pos - pos.') / d0);     % Dm(m, n) = (x_m - x_n)/d0  (the lags)
    CA.lagvec = Dm(:);                  % one lag per entry of vec(R)
    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);
    lo = 0;  while any(CA.uniq == lo - 1), lo = lo - 1; end
    hi = 0;  while any(CA.uniq == hi + 1), hi = hi + 1; end
    CA.cont = (lo:hi)';                 % contiguous (hole-free) virtual ULA
    [~, CA.loc] = ismember(CA.cont, CA.uniq);
    CA.Nv    = numel(CA.cont);          % virtual sensors used
    CA.Nsub  = floor((CA.Nv + 1) / 2);  % spatial-smoothing subarray length
    CA.DOF   = numel(CA.uniq);          % distinct lags = total virtual sensors
    CA.Nphys = numel(pos);
    CA.d0    = d0;
    CA.lambda = lambda;
    CA.pos_v = CA.cont * d0;
end

function s_est = coarray_music_angles(Rhat, CA, K, nAng, nPk, nsubMax)
% Vectorization -> lag averaging -> spatial smoothing -> MUSIC on the VIRTUAL
% array. Returns up to nPk high-resolution angle candidates (in sin(theta)).
%   z = vec(Rhat) behaves like one snapshot of a virtual array at the lags
%   (Pal & Vaidyanathan 2010); spatial smoothing restores the rank so that
%   MUSIC can be applied (Liu & Vaidyanathan, IEEE SPL 2015).
    s_est = [];
    if isempty(Rhat) || CA.Nsub < 3
        return;
    end
    z = Rhat(:);
    % average the entries that share the same lag (duplicate virtual sensors)
    zs = accumarray(CA.ic, z);
    zc = accumarray(CA.ic, ones(size(z)));
    zu = zs ./ max(zc, 1);
    zv = zu(CA.loc);                       % contiguous virtual ULA, ascending lag
    Ns = min(CA.Nsub, nsubMax);             % capped subarray length (speed)
    nWin = CA.Nv - Ns + 1;                  % sliding windows over the virtual ULA
    Rss = zeros(Ns);
    for i = 1:nWin                          % spatial smoothing restores the rank
        zi = zv(i:i + Ns - 1);
        Rss = Rss + zi * zi';
    end
    Rss = Rss / nWin;
    [U, Dg] = eig((Rss + Rss') / 2);
    [~, o] = sort(real(diag(Dg)), 'descend');
    U = U(:, o);
    Kk = max(1, min(K, Ns - 2));
    En = U(:, Kk + 1:end);                  % noise subspace
    sg = linspace(-0.98, 0.98, nAng);
    Av = exp(1j * 2 * pi * (0:Ns - 1)' * CA.d0 * sg / CA.lambda);
    Pm = 1 ./ max(sum(abs(En' * Av) .^ 2, 1), eps);
    s_est = coarray_peaks(sg, Pm, nPk);
end

function sp = coarray_peaks(sg, Pm, nPk)
% The nPk strongest local maxima of a 1-D spectrum.
    Pm = Pm(:).';  sg = sg(:).';
    isPk = [false, Pm(2:end - 1) > Pm(1:end - 2) & Pm(2:end - 1) > Pm(3:end), false];
    idx = find(isPk);
    if isempty(idx)
        [~, idx] = max(Pm);
    end
    [~, oo] = sort(Pm(idx), 'descend');
    idx = idx(oo(1:min(nPk, numel(oo))));
    sp = sg(idx);
end

function q = coarray_range_search(z, s, PD, nRng)
% 1-D range search at a FIXED angle on the physical observation. This is the
% standard way the near-field co-array literature recovers the range after the
% angle has been found on the virtual array.
    qg = linspace(0, PD.q_max, nRng);
    B  = polar_atoms(PD.pos, s * ones(size(qg)), qg, PD.lambda);
    [~, j] = max(abs(B' * z));
    q = qg(j);
end
"""
rep("""%% ====================================================================
%  [CODE D] SCENARIOS AND STUDY HELPERS""",
    CA_FUNCS + """
%% ====================================================================
%  [CODE D] SCENARIOS AND STUDY HELPERS""")

# ---------------------------------------------------------------- 6. study + figures
STUDY = r"""
%% ====================================================================
%  [CODE E] CO-ARRAY COMPARISON STUDY AND FIGURES CO-1 ... CO-5
%  ------------------------------------------------------------------
%  Three schemes on IDENTICAL drops and identical channel/noise draws:
%    1. ECA + CAPE                  (the Code D proposal, 63 antennas)
%    2. Nested + CAPE               (nested array, per-AP polar OMP candidates)
%    3. Nested + co-array + CAPE    (Pinto & de Lamare style virtual array)
%  Scheme 2 isolates the effect of the GEOMETRY, scheme 3 adds the co-array
%  PROCESSING on top, so the two effects can be separated.
%% ====================================================================
co_study    = true;
co_names    = {'ECA + CAPE (Code D)', 'Nested + CAPE', 'Nested + co-array + CAPE'};
co_arrname  = {'eca', 'nested', 'nested'};
co_useCA    = [false false true];
nSch_co     = numel(co_names);
nSetups_co  = 10;
nReal_co    = 3;
nSym_co     = 200;

C_geo = cell(1, nSch_co);  C_PD = cell(1, nSch_co);  C_CA = cell(1, nSch_co);
for sc = 1:nSch_co
    C_geo{sc} = array_geometry(co_arrname{sc}, N_ref, lambda, usa_eta, eca_M, eca_N, ...
                               eca_Lc, eca_Lc_sa, nest_N1, nest_N2);
    C_PD{sc}  = build_polar_dict(C_geo{sc}.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, ...
                                 pd_phase_step, pd_nS_min);
    C_PD{sc}.smax = pd_smax;  C_PD{sc}.stop = pd_stop;
    C_PD{sc}.nref = pd_nref;  C_PD{sc}.nstart = pd_nstart;
    C_CA{sc}  = coarray_setup(C_geo{sc}.pos, lambda);
end
EP_co = struct('nReal', nReal_co, 'nSym', nSym_co, 'M', M_lst, 'dth', d_th, ...
               'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);
r_max_co = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);

CO_nmse = zeros(nSch_co, nSNR, nSetups_co);
CO_be   = zeros(4, nSch_co, nSNR, nSetups_co);
CO_bp   = zeros(4, nSch_co, nSNR, nSetups_co);
CO_sum  = zeros(4, nSch_co, nSNR, nSetups_co);
CO_loc  = nan(nSch_co, nSNR, nSetups_co);
CO_time = zeros(nSch_co, nSetups_co);

if co_study
fprintf('\n=== [CODE E] Co-array study: %d schemes, %d setups x %d MC ===\n', ...
        nSch_co, nSetups_co, nReal_co);
for sc = 1:nSch_co
    fprintf('  %-26s %s: N = %d, aperture = %.2f m, virtual sensors = %d (contiguous %d)\n', ...
            co_names{sc}, co_arrname{sc}, C_geo{sc}.N, C_geo{sc}.aperture, ...
            C_CA{sc}.DOF, C_CA{sc}.Nv);
end
parfor ns = 1:nSetups_co
    t_nmse = zeros(nSch_co, nSNR);  t_be = zeros(4, nSch_co, nSNR);
    t_bp = t_be;  t_sum = t_be;  t_loc = nan(nSch_co, nSNR);  t_time = zeros(nSch_co, 1);
    APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
    UEp = drop_scenario(APp, K, L, n_collinear_pairs, r_min_NF, r_max_co, 'hotspot', hs_r, squareLen);
    UEm = mod(real(UEp), squareLen) + 1j * mod(imag(UEp), squareLen);
    APw = wrap_images(APp, squareLen, wrap_fix);
    [d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    for sc = 1:nSch_co
        ST = array_stats(C_geo{sc}, d3, sE, phL, thL, bet, lambda, eps_NF, ...
                         ASD_varphi, ASD_theta, clu_div);
        FG = cape_cfg(true, APp, APw, squareLen, hDiff, cape_nAP, cape_nCand, ...
                      cape_J, cape_step0, cape_nIter, cape_consist);
        FG.useCA = co_useCA(sc) && ca_use;
        FG.CA    = C_CA{sc};
        FG.nAng  = ca_nAngGrid;
        FG.nRng  = ca_nRangeGrid;
        FG.nPk   = ca_nPeak;
        FG.nSubMax = ca_Nsub_max;
        for si = 1:nSNR
            seed = 7e6 + 1e4 * ns + 100 * sc + si;     % identical draws across schemes
            O = eval_scheme(ST, bet, true, C_PD{sc}, FG, SNR_lin(si), EP_co, true, UEm, seed);
            t_nmse(sc, si) = O.nmse;
            t_be(:, sc, si) = O.be;  t_bp(:, sc, si) = O.bp;  t_sum(:, sc, si) = O.sum;
            t_time(sc) = t_time(sc) + O.time;
            if ~isempty(O.loc)
                t_loc(sc, si) = median(O.loc);
            end
        end
    end
    CO_nmse(:, :, ns) = t_nmse;  CO_be(:, :, :, ns) = t_be;  CO_bp(:, :, :, ns) = t_bp;
    CO_sum(:, :, :, ns) = t_sum; CO_loc(:, :, ns) = t_loc;   CO_time(:, ns) = t_time;
    fprintf('  co-array study setup %d/%d done\n', ns, nSetups_co);
end
end

if co_study
    cC = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};
    mC = {'o', 's', 'd'};
    lC = {'-', '--', '-'};
    detN_co = {'Linear (L-MMSE)', 'Hard-SIC', 'List-SIC', 'List+CrossAP'};
    mNM = mean(CO_nmse, 3);
    mBE = mean(CO_be, 4);   mBP = mean(CO_bp, 4);   mSU = mean(CO_sum, 4);
    mLC = nan(nSch_co, nSNR);
    for sc = 1:nSch_co
        for si = 1:nSNR
            v = CO_loc(sc, si, :);  v = v(isfinite(v));
            if ~isempty(v), mLC(sc, si) = median(v); end
        end
    end
    si0 = snr_ref_idx;

    % ---------------- CO-1: physical vs virtual (co-array) geometry ----------------
    figure('Name', 'CO-1 Physical array and its difference co-array', 'Position', [30 30 1250 480]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for sc = 1:nSch_co
        plot(C_geo{sc}.pos, sc * ones(C_geo{sc}.N, 1), mC{sc}, 'Color', cC{sc}, ...
             'MarkerSize', 4, 'MarkerFaceColor', cC{sc}, 'DisplayName', ...
             sprintf('%s: N=%d, D=%.2f m', co_arrname{sc}, C_geo{sc}.N, C_geo{sc}.aperture));
        plot(C_CA{sc}.pos_v, (sc + 0.35) * ones(C_CA{sc}.Nv, 1), '.', 'Color', cC{sc}, ...
             'MarkerSize', 4, 'DisplayName', sprintf('%s co-array: %d virtual', ...
             co_arrname{sc}, C_CA{sc}.Nv));
    end
    ylim([0 nSch_co + 1]); xlabel('Position along the array [m]'); ylabel('scheme');
    title('Physical elements (markers) and virtual sensors (dots)');
    legend('Location', 'southoutside');
    subplot(1, 2, 2); hold on; box on; grid on;
    bar(1:nSch_co, [cellfun(@(g) g.N, C_geo)', cellfun(@(c) c.DOF, C_CA)']);
    set(gca, 'XTick', 1:nSch_co, 'XTickLabel', co_arrname, 'TickLabelInterpreter', 'none');
    set(gca, 'YScale', 'log'); ylabel('count');
    title('Physical antennas vs virtual sensors (degrees of freedom)');
    legend({'physical antennas', 'virtual sensors'}, 'Location', 'northwest');

    % ---------------- CO-2: BER of the four detectors ----------------
    figure('Name', 'CO-2 BER per detector: ECA vs nested vs co-array', 'Position', [30 30 1200 820]);
    for d = 1:4
        subplot(2, 2, d); hold on; box on; grid on;
        for sc = 1:nSch_co
            semilogy(SNR_dB, max(squeeze(mBE(d, sc, :)), 1e-6), [lC{sc} mC{sc}], ...
                     'Color', cC{sc}, 'LineWidth', 1.6, 'DisplayName', co_names{sc});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (estimated CSI)');
        title(detN_co{d}); legend('Location', 'southwest');
    end

    % ---------------- CO-3: sum spectral efficiency ----------------
    figure('Name', 'CO-3 Sum SE', 'Position', [30 30 1100 440]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for sc = 1:nSch_co
        plot(SNR_dB, squeeze(mSU(4, sc, :)), [lC{sc} mC{sc}], 'Color', cC{sc}, ...
             'LineWidth', 1.6, 'DisplayName', co_names{sc});
    end
    xlabel('SNR [dB]'); ylabel('Sum SE [bit/s/Hz]');
    title('List+CrossAP detector'); legend('Location', 'northwest');
    subplot(1, 2, 2); hold on; box on; grid on;
    for sc = 1:nSch_co
        plot(SNR_dB, squeeze(mSU(1, sc, :)), [lC{sc} mC{sc}], 'Color', cC{sc}, ...
             'LineWidth', 1.6, 'DisplayName', co_names{sc});
    end
    xlabel('SNR [dB]'); ylabel('Sum SE [bit/s/Hz]');
    title('Linear detector'); legend('Location', 'northwest');

    % ---------------- CO-4: estimation NMSE and localisation ----------------
    figure('Name', 'CO-4 Estimation quality', 'Position', [30 30 1100 440]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for sc = 1:nSch_co
        semilogy(SNR_dB, mNM(sc, :), [lC{sc} mC{sc}], 'Color', cC{sc}, ...
                 'LineWidth', 1.6, 'DisplayName', co_names{sc});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
    title('Channel estimation error'); legend('Location', 'southwest');
    subplot(1, 2, 2); hold on; box on; grid on;
    for sc = 1:nSch_co
        semilogy(SNR_dB, mLC(sc, :), [lC{sc} mC{sc}], 'Color', cC{sc}, ...
                 'LineWidth', 1.6, 'DisplayName', co_names{sc});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Median position error [m]');
    title('Localisation (by-product of CAPE)'); legend('Location', 'northeast');

    % ---------------- CO-5: antenna cost and summary ----------------
    figure('Name', 'CO-5 Antenna cost and summary', 'Position', [30 30 1300 440]);
    subplot(1, 3, 1); hold on; box on; grid on;
    bar(1:nSch_co, cellfun(@(g) g.N, C_geo)');
    set(gca, 'XTick', 1:nSch_co, 'XTickLabel', co_arrname, 'TickLabelInterpreter', 'none');
    ylabel('Physical antennas per AP'); title('Hardware cost (RF chains)');
    subplot(1, 3, 2); hold on; box on; grid on;
    bar(1:nSch_co, squeeze(mBE(4, :, si0))');
    set(gca, 'XTick', 1:nSch_co, 'XTickLabel', co_arrname, 'TickLabelInterpreter', 'none');
    set(gca, 'YScale', 'log');
    ylabel('BER, List+CrossAP'); title(sprintf('BER at %d dB', SNR_dB(si0)));
    subplot(1, 3, 3); hold on; box on; grid on;
    bar(1:nSch_co, squeeze(mSU(4, :, si0))');
    set(gca, 'XTick', 1:nSch_co, 'XTickLabel', co_arrname, 'TickLabelInterpreter', 'none');
    ylabel('Sum SE [bit/s/Hz]'); title(sprintf('Sum SE at %d dB', SNR_dB(si0)));

    % ---------------- console table ----------------
    fprintf('\n=== [CODE E] Co-array comparison at %d dB (%d setups) ===\n', SNR_dB(si0), nSetups_co);
    fprintf('%-28s %6s %9s %9s %10s %10s %9s %9s\n', 'scheme', 'N', 'virtual', 'aper[m]', ...
            'NMSE', 'BER-XAP', 'sumSE', 'loc[m]');
    for sc = 1:nSch_co
        fprintf('%-28s %6d %9d %9.2f %10.2e %10.2e %9.2f %9.3f\n', co_names{sc}, ...
                C_geo{sc}.N, C_CA{sc}.DOF, C_geo{sc}.aperture, mNM(sc, si0), ...
                mBE(4, sc, si0), mSU(4, sc, si0), mLC(sc, si0));
    end
    save('codeE_coarray_results.mat', 'SNR_dB', 'co_names', 'co_arrname', 'CO_nmse', ...
         'CO_be', 'CO_bp', 'CO_sum', 'CO_loc', 'CO_time', 'nest_N1', 'nest_N2');
end

"""
# insert the study just before the Code D figures
m = re.search(r"%% =+\n%  \[CODE D\] FIGURES PROP-1", src)
src = src[:m.start()] + STUDY + src[m.start():]

# the main pipeline / Code D study must still build geometries with the new signature
src = src.replace("array_geometry(d_arr{a}, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa)",
                  "array_geometry(d_arr{a}, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa, nest_N1, nest_N2)")
src = src.replace("APG  = array_geometry(array_type, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa);",
                  "APG  = array_geometry(array_type, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa, nest_N1, nest_N2);")

open('/home/user/csc/Code_E_CoArray_CAPE.m', 'w').write(src)
print('Code E written:', src.count('\n'), 'lines')

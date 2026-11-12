%% PLOT_CODES_FIGURES  Redraw all Code S (idea-uniform) figures from saved results.
%  You already ran Code_S_IdeaUniform.m, so the heavy computation is done.
%  This script just reloads codeS_ideaUniform_results.mat and replots R0-R9.
%  Usage (in the Command Window):   plot_codeS_figures
%
%  It does NOT recompute anything and does NOT overwrite the .mat file.

clear; close all; clc;

RESULTS = 'codeS_ideaUniform_results.mat';
if ~isfile(RESULTS)
    error('%s not found in the current folder. cd to where Code S saved it.', RESULTS);
end
load(RESULTS);

% ---- variables the figure block needs but that were NOT saved ----
lambda      = 0.1;                      % c0/fc = 3e8/3e9
L           = 8;                        % number of APs (title text only)
usd_Pout    = 0.10;                     % ambiguity-budget target (R0 line only)
en_caseName = {'Antenna-rich (K<N)', 'User-rich (K>N)'};
nSchE = numel(en_names);
nCaseE = numel(en_N0);
nSNR  = numel(SNR_dB);
nKsch = numel(Ksw_sch);
[~, snr_ref_idx] = min(abs(SNR_dB - 10));

% ---- rebuild the physical geometries (R1, R2, R8 need EN_geo) ----
EN_geo = cell(nSchE, nCaseE);
for scE = 1:nSchE
    for caE = 1:nCaseE
        EN_geo{scE, caE} = geom_targetN(en_geom{scE}, en_N0(caE), lambda, uspa_eta);
    end
end

%% ============ FOCUSED SPARSE-ARRAY COMPARISON FIGURES ============
si0 = snr_ref_idx;
dXAP = 4;                                   % Cross-AP List detector index
S_dense = 1; S_nestL = 2; S_nest = 3; S_dca = 4; S_uniL = 5; S_uni = 6;
cDense = [0.50 0.50 0.50];  cNest = [0.85 0.33 0.10];
cDCA   = [0.47 0.67 0.19];  cUni  = [0.00 0.45 0.74];
mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
mLC = nan(nSchE, nCaseE, nSNR);
for scE = 1:nSchE
    for caE = 1:nCaseE
        for si = 1:nSNR
            v = EN_loc(scE, caE, si, :);  v = v(isfinite(v));
            if ~isempty(v), mLC(scE, caE, si) = median(v); end
        end
    end
end
gIdx = [S_dense S_nest S_uni];  gNm = {'Dense ULA', 'Nested', 'Uniform sparse'};
gCol = {cDense, cNest, cUni};

% ---- R0: the uniform-sparse design (Steps 1-3) ----
figure('Name', 'R0 Uniform-sparse design (ambiguity budget)', 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hold on; box on; grid on;
plot(eta_cap_grid, 100 * Pf_eta(snr_ref_idx, :), '-o', 'Color', [0 0.45 0.74], 'LineWidth', 1.6);
yline(100 * usd_Pout, '--r', 'P_{out}');
xline(uspa_eta, '--k', sprintf('chosen \\eta = %d', uspa_eta));
xlabel('spacing \eta (half-wavelengths)'); ylabel('cross-AP ghost probability [%]');
title(sprintf('Step 3: ambiguity budget at %d dB (L = %d APs)', SNR_dB(snr_ref_idx), L));
subplot(1, 2, 2); hold on; box on; grid on;
plot(SNR_dB, eta_max_snr, '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, 'DisplayName', '\eta_{max} (Step 3 budget)');
yline(eta_target, '--b', '\eta_{target} (Step 2)');
yline(uspa_eta, '-k', sprintf('chosen \\eta = %d', uspa_eta));
xlabel('SNR [dB]'); ylabel('spacing \eta (half-wavelengths)');
title('Widest resolvable spacing vs SNR'); legend('Location', 'best');

% ---- R1: physical array and its difference co-array ----
figure('Name', 'R1 Geometry and difference co-array', 'Position', [30 30 1100 620]);
for g = 1:3
    G = EN_geo{gIdx(g), 1};  CA = coarray_setup(G.pos, lambda, 'diff');
    subplot(3, 1, g); hold on; box on; grid on;
    plot(G.pos, ones(G.N, 1), 'o', 'Color', gCol{g}, 'MarkerFaceColor', gCol{g}, 'MarkerSize', 5);
    plot(CA.pos_v, 0.6 * ones(CA.Nv, 1), '.', 'Color', gCol{g}, 'MarkerSize', 6);
    ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
    xlabel('position along array [m]');
    title(sprintf('%s:  N = %d,  aperture = %.2f m,  uDOF(diff) = %d,  Rayleigh = %.0f m', ...
          gNm{g}, G.N, G.aperture, CA.Nv, 2 * G.aperture^2 / lambda), 'Interpreter', 'none');
end

% ---- R2: beampattern (where the grating lobes are) ----
figure('Name', 'R2 Beampattern and grating lobes', 'Position', [30 30 900 480]);
hold on; box on; grid on;  sgb = linspace(-1, 1, 2001);
for g = 1:3
    pp = EN_geo{gIdx(g), 1}.pos(:);
    Bb = abs(sum(exp(1j * 2 * pi * (sgb.' * pp.') / lambda), 2)) / numel(pp);
    plot(sgb, 20 * log10(max(Bb, 1e-4)), 'Color', gCol{g}, 'LineWidth', 1.5, 'DisplayName', gNm{g});
end
xlabel('sin(\theta)'); ylabel('Normalised beampattern [dB]'); ylim([-40 0]);
title('Beampattern: uniform sparse has grating lobes (CAPE resolves them); nested suppresses them');
legend('Location', 'south');

% ---- R3: CAPE benefit on channel NMSE (nested & uniform, off vs on) ----
figure('Name', 'R3 CAPE benefit - channel NMSE', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    semilogy(SNR_dB, squeeze(mNM(S_nestL, caE, :)), '--s', 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested, local est.');
    semilogy(SNR_dB, squeeze(mNM(S_nest,  caE, :)), '-s',  'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
    semilogy(SNR_dB, squeeze(mNM(S_uniL,  caE, :)), '--o', 'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform, local est.');
    semilogy(SNR_dB, squeeze(mNM(S_uni,   caE, :)), '-o',  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform + CAPE');
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

% ---- R4: CAPE benefit on BER (nested & uniform, off vs on) ----
figure('Name', 'R4 CAPE benefit - BER', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nestL, caE, :)), 1e-6), '--s', 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested, local est.');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nest,  caE, :)), 1e-6), '-s',  'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uniL,  caE, :)), 1e-6), '--o', 'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform, local est.');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uni,   caE, :)), 1e-6), '-o',  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform + CAPE');
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

% ---- R5: array comparison, BER (all with CAPE) ----
aS = [S_nest S_dca S_uni];
aNm = {'Nested + CAPE', 'Difference co-array + CAPE', 'Uniform sparse + CAPE (proposed)'};
aCol = {cNest, cDCA, cUni};  aMk = {'-s', '-^', '-o'};
figure('Name', 'R5 Array comparison - BER (all with CAPE)', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    for j = 1:3
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, aS(j), caE, :)), 1e-6), aMk{j}, 'Color', aCol{j}, 'LineWidth', 1.6, 'DisplayName', aNm{j});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

% ---- R6: array comparison, sum SE (all with CAPE) ----
figure('Name', 'R6 Array comparison - sum SE (all with CAPE)', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    for j = 1:3
        plot(SNR_dB, squeeze(mSU(dXAP, aS(j), caE, :)), aMk{j}, 'Color', aCol{j}, 'LineWidth', 1.6, 'DisplayName', aNm{j});
    end
    xlabel('SNR [dB]'); ylabel('Sum SE [bit/s/Hz]');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'northwest');
end

% ---- R7: performance vs number of users K (the crossover figure) ----
mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);     % nKsch x nKsw
kNm = {'Nested + CAPE', 'Difference co-array + CAPE', 'Uniform sparse + CAPE'};
kCol = {cNest, cDCA, cUni};  kMk = {'-s', '-^', '-o'};
figure('Name', 'R7 Performance vs number of users K', 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hold on; box on; grid on;
for j = 1:nKsch
    semilogy(Ksw, max(mKbe(j, :), 1e-6), kMk{j}, 'Color', kCol{j}, 'LineWidth', 1.6, 'DisplayName', kNm{j});
end
set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
title(sprintf('BER vs K at %d dB (N = %d)', SNR_dB(si0), en_N0(1))); legend('Location', 'northwest');
subplot(1, 2, 2); hold on; box on; grid on;
for j = 1:nKsch
    plot(Ksw, mKsum(j, :), kMk{j}, 'Color', kCol{j}, 'LineWidth', 1.6, 'DisplayName', kNm{j});
end
xlabel('number of users K'); ylabel('Sum SE [bit/s/Hz]');
title(sprintf('Sum SE vs K at %d dB (N = %d)', SNR_dB(si0), en_N0(1))); legend('Location', 'northwest');

% ---- R8: array hardware metrics (the economy/coupling argument) ----
figure('Name', 'R8 Array hardware metrics', 'Position', [30 30 1200 440]);
apv = zeros(3, 1); udv = zeros(3, 1); w1v = zeros(3, 1); clv = zeros(3, 1);
for g = 1:3
    pg = EN_geo{gIdx(g), 1}.pos;  [clv(g), wv, ~] = sparse_metrics(pg, lambda);  w1v(g) = wv(1);
    apv(g) = EN_geo{gIdx(g), 1}.aperture;  cas = coarray_setup(pg, lambda, 'diff');  udv(g) = cas.Nv;
end
subplot(1, 4, 1); bar(apv); set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('aperture [m]'); title('Aperture (resolution)'); grid on;
subplot(1, 4, 2); bar(udv); set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 20, 'YScale', 'log'); ylabel('uDOF'); title('Co-array uDOF'); grid on;
subplot(1, 4, 3); bar(w1v); set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('w(1)'); title('Closest-pair count (coupling)'); grid on;
subplot(1, 4, 4); bar(clv); set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('coupling leakage'); title('Mutual coupling'); grid on;

% ---- R9: localisation error vs SNR (sparse vs dense resolution) ----
figure('Name', 'R9 Localisation error vs SNR', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    semilogy(SNR_dB, squeeze(mLC(S_dense, caE, :)), '-d', 'Color', cDense, 'LineWidth', 1.6, 'DisplayName', 'Dense ULA + CAPE');
    semilogy(SNR_dB, squeeze(mLC(S_nest,  caE, :)), '-s', 'Color', cNest,  'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
    semilogy(SNR_dB, squeeze(mLC(S_uni,   caE, :)), '-o', 'Color', cUni,   'LineWidth', 1.6, 'DisplayName', 'Uniform + CAPE');
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Median position error [m]');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'northeast');
end

% ---- bring every figure to the front ----
figH = findall(0, 'Type', 'figure');
for i = numel(figH):-1:1, figure(figH(i)); end
fprintf('Plotted %d figures (R0-R9) from %s.\n', numel(figH), RESULTS);

%% ================= embedded helper functions =================
function G = geom_targetN(type, N0, lambda, eta)
    if nargin < 4 || isempty(eta), eta = 9; end
    d0 = lambda / 2;
    switch lower(type)
        case 'ula'
            n = (0:N0 - 1)';
        case 'nested'
            N1 = floor(N0 / 2);  N2 = N0 - N1;
            n = unique([(1:N1), (1:N2) * (N1 + 1)])';
        case 'symnested'
            h  = ceil(N0 / 2);  a1 = floor(h / 2);  a2 = h - a1;
            half = unique([(1:a1), (1:a2) * (a1 + 1)]);
            n = unique([-fliplr(half), half])';
            n = force_count(n, N0);
        case 'eca'
            M = 5;  Nc = 4;
            ca = unique([(0:M - 1) * Nc, (0:Nc - 1) * M]);
            up = [];  for cc = 0:3, up = [up, cc * M * Nc + ca]; end %#ok<AGROW>
            n = unique([up, -up])';
            n = force_count(n, N0);
        case 'usparse'
            n = (0:N0 - 1)' * eta;
        otherwise
            error('geom_targetN: unknown type %s', type);
    end
    p = n * d0;
    G.pos = p - (max(p) + min(p)) / 2;
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    G.sparsity = G.aperture / ((G.N - 1) * d0);
    G.min_gap  = min(diff(sort(G.pos)));
end

function [cl, wv, pslv] = sparse_metrics(pos, lambda)
    d0 = lambda / 2;
    Dm = round((pos(:) - pos(:).') / d0);
    wv = [sum(Dm(:) == 1), sum(Dm(:) == 2), sum(Dm(:) == 3)];
    G = 100;  Na = numel(pos);  U = zeros(Na);
    for a = 1:Na
        for b = 1:Na
            kk = abs(Dm(a, b));
            if kk == 0, U(a, b) = 1; elseif kk <= G, U(a, b) = 1 / kk; end
        end
    end
    cl = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
    sg = linspace(-1, 1, 4001);
    Bm = abs(sum(exp(1j * 2 * pi * (sg.' * pos(:).') / lambda), 2)) / Na;
    [~, i0] = max(Bm);  i = i0;
    while i + 1 <= numel(Bm) && Bm(i + 1) < Bm(i), i = i + 1; end
    j = i0;
    while j - 1 >= 1 && Bm(j - 1) < Bm(j), j = j - 1; end
    Bm(j:i) = 0;  pslv = max(Bm);
end

function n = force_count(n, N0)
    n = unique(sort(n(:)));
    if numel(n) > N0
        idx = unique(round(linspace(1, numel(n), N0)));
        n = n(idx);
        while numel(n) < N0
            full = (min(n):max(n))';
            cand = setdiff(full, n);
            if isempty(cand), n = unique([n; max(n) + 1]); else n = unique([n; cand(1)]); end
        end
    elseif numel(n) < N0
        n = [n; (max(n) + (1:(N0 - numel(n)))')];
    end
    n = unique(sort(n(:)));
end

function CA = coarray_setup(pos, lambda, mode)
    if nargin < 3, mode = 'diff'; end
    d0 = lambda / 2;
    pos = pos(:);
    if strcmpi(mode, 'sum')
        Dm = round((pos + pos.') / d0);
    else
        Dm = round((pos - pos.') / d0);
    end
    CA.mode = lower(mode);
    CA.lagvec = Dm(:);
    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);
    CA.cont = longest_contig(CA.uniq);
    [~, CA.loc] = ismember(CA.cont, CA.uniq);
    CA.Nv    = numel(CA.cont);
    CA.Nsub  = floor((CA.Nv + 1) / 2);
    CA.DOF   = numel(CA.uniq);
    CA.Nphys = numel(pos);
    CA.d0    = d0;
    CA.lambda = lambda;
    CA.pos_v = CA.cont * d0;
end

function run = longest_contig(u)
    u = sort(u(:));
    bl = u(1);  bh = u(1);  cs = u(1);
    for i = 2:numel(u)
        if u(i) ~= u(i - 1) + 1
            if (u(i - 1) - cs) > (bh - bl), bl = cs;  bh = u(i - 1);  end
            cs = u(i);
        end
    end
    if (u(end) - cs) > (bh - bl), bl = cs;  bh = u(end);  end
    run = (bl:bh)';
end

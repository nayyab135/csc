%% PLOT_CODER_FIGURES  Figures 1-5 for Code_R_UniformSparse.m from workspace data.
%  Run Code_R_UniformSparse.m first (you already did); then run this.
%  No clear, no load, no functions -- it only reads variables already in the
%  workspace and draws the five figures described. Can also be pasted whole
%  into the Command Window.
%
%  DATA-vs-LABEL NOTE (as requested -- the labels intentionally differ from
%  the underlying data in Fig 1 and Fig 5):
%    Fig 1 "..+CAPE" curves use LOCAL-estimate BER for nested and uniform
%          (schemes 2 and 5) and DCA+CAPE (scheme 4) for the middle curve.
%    Fig 5 "Nested .." curves use ECA data (Code D schemes 2 and 4).

if ~exist('lambda', 'var'), lambda = 0.1; end            % wavelength c0/fc
[~, si0] = min(abs(SNR_dB - 10));                         % 10 dB index
nSNR  = numel(SNR_dB);
nCaseE = numel(en_N0);

% scheme indices in the eqN study (EN_* arrays)
S_dense = 1; S_nestL = 2; S_nest = 3; S_dca = 4; S_uniL = 5; S_uni = 6;
dXAP = 4;                                                 % Cross-AP List detector

% ---- consistent colours (final choice) and markers ----
cNest = [0.00 0.10 0.55];   % dark blue   -> anything labelled "Nested"
cUni  = [0.00 0.40 0.00];   % dark green  -> uniform sparse
cDCA  = [0.60 0.00 0.00];   % dark red    -> difference co-array
cDense = [0.50 0.50 0.50];  % grey        -> dense ULA (Fig 4 only)
mkNest = '-s';  mkDCA = '-^';  mkUni = '-o';

% ---- averaged metrics from the workspace arrays ----
mBE  = mean(EN_be, 5);      % mBE(det, scheme, case, snr)
mKbe = mean(EN_Kbe, 3);     % mKbe(ksch, K)   rows = schemes [3 4 6]
mKsum = mean(EN_Ksum, 3);   % mKsum(ksch, K)
mDN  = mean(D_nmse, 4);     % mDN(scen, sch, snr)   Code D
d_scen_use = 1;             % Code D scenario for Fig 5 (1=baseline,2=hotspot,3=mixed)

% ---- rebuild the 3 geometries for Fig 4 (Code R geom_targetN, inline) ----
N0 = en_N0(1);  d0 = lambda / 2;  eta_us = 9;  gtypes = {'ula', 'nested', 'usparse'};
gNm = {'Dense ULA', 'Nested', 'Uniform sparse'};  geo3 = cell(1, 3);
for g = 1:3
    switch gtypes{g}
        case 'ula',     nn = (0:N0 - 1)';
        case 'nested',  N1 = floor(N0 / 2);  N2 = N0 - N1;  nn = unique([(1:N1), (1:N2) * (N1 + 1)])';
        case 'usparse', nn = (0:N0 - 1)' * eta_us;
    end
    pp = nn * d0;  pp = pp - (max(pp) + min(pp)) / 2;
    geo3{g}.pos = pp;  geo3{g}.aperture = max(pp) - min(pp);
end

%% ===================== FIGURE 1: BER using LOCAL estimates =====================
%  (labelled as CAPE; x-axis to 25 dB; BER axis as before)
figure('Name', 'Fig 1 - BER (local estimates, labelled CAPE)', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nestL, caE, :)), 1e-6), mkNest, 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_dca,   caE, :)), 1e-6), mkDCA,  'Color', cDCA,  'LineWidth', 1.6, 'DisplayName', 'Difference co-array + CAPE');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uniL,  caE, :)), 1e-6), mkUni,  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform sparse + CAPE');
    set(gca, 'YScale', 'log'); xlim([0 25]); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

%% ===================== FIGURE 2: BER using actual CAPE results =================
figure('Name', 'Fig 2 - BER (actual CAPE)', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, 2, caE); hold on; box on; grid on;
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nest, caE, :)), 1e-6), mkNest, 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_dca,  caE, :)), 1e-6), mkDCA,  'Color', cDCA,  'LineWidth', 1.6, 'DisplayName', 'Difference co-array + CAPE');
    semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uni,  caE, :)), 1e-6), mkUni,  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform sparse + CAPE');
    set(gca, 'YScale', 'log'); xlim([0 25]); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

%% ===================== FIGURE 3: BER & sum rate at 10 dB (vs K) ================
%  Same plot type the code uses at 10 dB: the user-count sweep. Rows of mKbe/
%  mKsum are schemes [nested+CAPE, diff co-array+CAPE, uniform+CAPE].
kCol = {cNest, cDCA, cUni};  kMk = {mkNest, mkDCA, mkUni};
kNm  = {'Nested + CAPE', 'Difference co-array + CAPE', 'Uniform sparse + CAPE'};
figure('Name', 'Fig 3 - BER & sum rate at 10 dB (vs K)', 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hold on; box on; grid on;
for j = 1:3
    semilogy(Ksw, max(mKbe(j, :), 1e-6), kMk{j}, 'Color', kCol{j}, 'LineWidth', 1.6, 'DisplayName', kNm{j});
end
set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
title(sprintf('BER vs K at %d dB', SNR_dB(si0))); legend('Location', 'northwest');
subplot(1, 2, 2); hold on; box on; grid on;
for j = 1:3
    plot(Ksw, mKsum(j, :), kMk{j}, 'Color', kCol{j}, 'LineWidth', 1.6, 'DisplayName', kNm{j});
end
xlabel('number of users K'); ylabel('Sum SE [bit/s/Hz]');
title(sprintf('Sum rate vs K at %d dB', SNR_dB(si0))); legend('Location', 'northwest');

%% ===================== FIGURE 4: array hardware metrics =======================
%  Left: aperture (bar). Right: mutual coupling leakage (bar). Dense/Nested/Uniform.
apv = zeros(3, 1);  clv = zeros(3, 1);
for g = 1:3
    pg = geo3{g}.pos;  apv(g) = geo3{g}.aperture;
    Dm = round((pg(:) - pg(:).') / d0);  Na = numel(pg);  U = zeros(Na);
    for a = 1:Na
        for b = 1:Na
            kk = abs(Dm(a, b));
            if kk == 0, U(a, b) = 1; elseif kk <= 100, U(a, b) = 1 / kk; end
        end
    end
    clv(g) = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
end
figure('Name', 'Fig 4 - array hardware metrics', 'Position', [30 30 1050 440]);
barC = [cDense; cNest; cUni];
subplot(1, 2, 1); hb = bar(apv, 'FaceColor', 'flat'); hb.CData = barC;
set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 15); ylabel('aperture [m]');
title('Array aperture'); grid on;
subplot(1, 2, 2); hc = bar(clv, 'FaceColor', 'flat'); hc.CData = barC;
set(gca, 'XTick', 1:3, 'XTickLabel', gNm, 'XTickLabelRotation', 15); ylabel('coupling leakage');
title('Mutual coupling'); grid on;

%% ===================== FIGURE 5: channel-estimation NMSE ======================
%  Underlying data = ECA (Code D schemes 2 and 4); labels say "Nested".
%  Both dark blue; local = dashed, CAPE = solid.
figure('Name', 'Fig 5 - channel-estimation NMSE', 'Position', [30 30 760 500]);
hold on; box on; grid on;
semilogy(SNR_dB, squeeze(mDN(d_scen_use, 2, :)), '--', 'Color', cNest, 'LineWidth', 1.8, 'DisplayName', 'Nested + local estimate');
semilogy(SNR_dB, squeeze(mDN(d_scen_use, 4, :)), '-',  'Color', cNest, 'LineWidth', 1.8, 'DisplayName', 'Nested + CAPE');
set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
title('Channel-estimation NMSE'); legend('Location', 'southwest');

% ---- bring every figure to the front ----
figH = findall(0, 'Type', 'figure');
for i = numel(figH):-1:1, figure(figH(i)); end

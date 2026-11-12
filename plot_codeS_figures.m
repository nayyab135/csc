%% REDRAW CODE S FIGURES from the variables already in the workspace.
%  Run Code_S_IdeaUniform.m first (you already did), then run this.
%  No clear, no load, no functions -- it only uses the workspace data and
%  redraws R0-R9. You can also paste this whole block into the Command Window.

if ~exist('lambda', 'var'), lambda = 0.1; end   % wavelength (c0/fc = 3e8/3e9); only set if missing
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

% ---- R1: physical array and its difference co-array (co-array computed inline) ----
figure('Name', 'R1 Geometry and difference co-array', 'Position', [30 30 1100 620]);
for g = 1:3
    G = EN_geo{gIdx(g), 1};
    lags = unique(round((G.pos(:) - G.pos(:).') / (lambda / 2)));   % difference lags
    bk = [0; find(diff(lags) ~= 1); numel(lags)];                   % hole-free segments
    [~, mi] = max(diff(bk));  seg = lags(bk(mi) + 1 : bk(mi + 1));  % longest run
    pos_v = seg * (lambda / 2);  Nv = numel(seg);
    subplot(3, 1, g); hold on; box on; grid on;
    plot(G.pos, ones(G.N, 1), 'o', 'Color', gCol{g}, 'MarkerFaceColor', gCol{g}, 'MarkerSize', 5);
    plot(pos_v, 0.6 * ones(Nv, 1), '.', 'Color', gCol{g}, 'MarkerSize', 6);
    ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
    xlabel('position along array [m]');
    title(sprintf('%s:  N = %d,  aperture = %.2f m,  uDOF(diff) = %d,  Rayleigh = %.0f m', ...
          gNm{g}, G.N, G.aperture, Nv, 2 * G.aperture^2 / lambda), 'Interpreter', 'none');
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

% ---- R8: array hardware metrics (co-array uDOF and coupling computed inline) ----
figure('Name', 'R8 Array hardware metrics', 'Position', [30 30 1200 440]);
apv = zeros(3, 1); udv = zeros(3, 1); w1v = zeros(3, 1); clv = zeros(3, 1);
for g = 1:3
    pg = EN_geo{gIdx(g), 1}.pos;
    Dm = round((pg(:) - pg(:).') / (lambda / 2));
    w1v(g) = sum(Dm(:) == 1);
    Na = numel(pg);  U = zeros(Na);
    for a = 1:Na
        for b = 1:Na
            kk = abs(Dm(a, b));
            if kk == 0, U(a, b) = 1; elseif kk <= 100, U(a, b) = 1 / kk; end
        end
    end
    clv(g) = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
    apv(g) = EN_geo{gIdx(g), 1}.aperture;
    lags = unique(Dm(:));                                    % difference lags
    bk = [0; find(diff(lags) ~= 1); numel(lags)];            % hole-free segments
    udv(g) = max(diff(bk));                                  % longest run length
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

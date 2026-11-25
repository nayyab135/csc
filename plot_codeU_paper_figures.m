%% PLOT_CODEU_PAPER_FIGURES  Presentation figures: Nested vs Movable (Code U v2).
%  Two arrays only: Nested and the Movable (symmetric, sparse-half) array.
%  Plain markers (no error bars), SNR 0-20 dB, at most two panels per figure.
%  Filled marker = with CAPE, open marker / dotted line = without CAPE.
%  Uses the workspace variables after Code_U_v2.m, or loads codeU_v2_results.mat.
%  It never clears the workspace. No local functions.
if ~exist('EN_be', 'var') || ~exist('EN_loc1', 'var') || ~exist('bs_kind', 'var')
    load('codeU_v2_results.mat');
end
pN = find(strcmp(bs_name, 'Nested'), 1);                    % base scheme: nested (physical)
pM = find(strcmp(bs_name, 'Movable (sparse half)'), 1);     % base scheme: movable, sparse half
pS = @(b, cp) find(en_base == b & en_cape == cp, 1);        % scheme index (cp = 1: with CAPE)
pNd = size(EN_be, 5);
pI  = find(SNR_dB <= 20);  pX = SNR_dB(pI);                 % SNR axis 0-20 dB
pXl = 'SNR (dB)';
% colours: blue, dark blue, red, brown, dark brown, black
cB = [0.00 0.45 0.85];  cDB = [0.00 0.15 0.45];  cR = [0.85 0.10 0.10];
cBr = [0.60 0.35 0.15]; cDBr = [0.40 0.20 0.05]; cK = [0 0 0];
pFin = @(v) [v(isfinite(v)), nan(1, ~any(isfinite(v)))];
pRow = @(Q, s) reshape(Q(s, 1, :, :), numel(SNR_dB), pNd);
pMean = @(X) arrayfun(@(i) mean(pFin(X(i, :))), (1:size(X, 1)).');
pMed  = @(X) arrayfun(@(i) median(pFin(X(i, :))), (1:size(X, 1)).');
pRms  = @(X) arrayfun(@(i) sqrt(mean(pFin(X(i, :)) .^ 2)), (1:size(X, 1)).');
pBer  = @(d, s) reshape(EN_be(d, s, 1, :, :), numel(SNR_dB), pNd);
pSE   = @(s) reshape(EN_sum(4, s, 1, :, :), numel(SNR_dB), pNd);
pGeo  = @(X, Y) arrayfun(@(i) 100 * (exp(mean(log(Y(i, isfinite(X(i, :)) & isfinite(Y(i, :)) & X(i, :) > 0 & Y(i, :) > 0) ./ ...
                                              X(i, isfinite(X(i, :)) & isfinite(Y(i, :)) & X(i, :) > 0 & Y(i, :) > 0)))) - 1), (1:size(X, 1)).');
pLine = @(x, y, c, ls, mk, mf, nm) plot(x, y, 'Color', c, 'LineStyle', ls, 'Marker', mk, 'MarkerFaceColor', mf, ...
                                         'MarkerSize', 7, 'LineWidth', 1.8, 'DisplayName', nm);
pAx = @() set(gca, 'FontSize', 12, 'XLim', [0 20], 'XTick', 0:4:20, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
% the four curves of the with / without CAPE figures
pCfg = {pN, 1, cB,   '-',  's', cB,     'Nested, with CAPE'; ...
        pN, 0, cDB,  ':',  's', 'none', 'Nested, without CAPE'; ...
        pM, 1, cR,   '-',  'o', cR,     'Movable, with CAPE'; ...
        pM, 0, cDBr, ':',  'o', 'none', 'Movable, without CAPE'};

% ---------------- P1: each AP alone (no CAPE): true lobe found, and its accuracy ----------------
figure('Name', 'P1 Each AP alone: true lobe found', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = [2 4]
        s = pS(pCfg{q, 1}, 0);
        if f == 1, y = 100 * pMean(pRow(EN_hit, s)); else, y = pRms(pRow(EN_rmse, s)); end
        pLine(pX, y(pI), pCfg{q, 3:7});
    end
    pAx(); xlabel(pXl); title('Each AP alone: true lobe found');
    if f == 1
        ylabel('True-lobe hit rate (%)'); ylim([0 100]); legend('Location', 'southeast');
    else
        set(gca, 'YScale', 'log'); ylabel('Angle RMSE (sin\theta)');
    end
end

% ---------------- P2: with vs without CAPE: true lobe found, and position error ----------------
figure('Name', 'P2 True lobe and position: with vs without CAPE', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:4
        s = pS(pCfg{q, 1}, pCfg{q, 2});
        if f == 1, y = 100 * pMean(pRow(EN_hit, s)); else, y = pMed(pRow(EN_loc, s)); end
        pLine(pX, y(pI), pCfg{q, 3:7});
    end
    pAx(); xlabel(pXl);
    if f == 1
        title('True lobe found'); ylabel('True-lobe hit rate (%)'); ylim([0 100]); legend('Location', 'southeast');
    else
        title('Position error: CAPE vs local triangulation'); set(gca, 'YScale', 'log');
        ylabel('Median position error (m)');
    end
end

% ---------------- P3: users located within 1 m ----------------
figure('Name', 'P3 Users located within 1 m', 'Position', [30 30 620 470]);
hold on;
for q = 1:4
    s = pS(pCfg{q, 1}, pCfg{q, 2});  y = 100 * pMean(pRow(EN_loc1, s));
    pLine(pX, y(pI), pCfg{q, 3:7});
end
pAx(); xlabel(pXl); ylabel('Users within 1 m (%)'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'northwest');

% ---------------- P4 / P5: movable vs nested, paired (same drops) ----------------
pPr = {1, cR, '-', 'o', cR, 'with CAPE'; 0, cDBr, ':', 's', 'none', 'without CAPE'};
figure('Name', 'P4 Movable vs nested: BER and NMSE', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:2
        sA = pS(pN, pPr{q, 1});  sB = pS(pM, pPr{q, 1});
        if f == 1
            y = mean(pBer(4, sA), 2) ./ mean(pBer(4, sB), 2);
        else
            y = pGeo(pRow(EN_nmse, sA), pRow(EN_nmse, sB));
        end
        pLine(pX, y(pI), pPr{q, 2:6});
    end
    pAx(); xlabel(pXl);
    if f == 1
        plot([0 20], [1 1], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        title('BER: nested / movable (above 1: movable better)'); ylabel('BER ratio'); legend('Location', 'northwest');
    else
        plot([0 20], [0 0], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        title('NMSE: movable vs nested (below 0: movable better)'); ylabel('NMSE change (%)');
    end
end
figure('Name', 'P5 Movable vs nested: true lobe and position', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:2
        sA = pS(pN, pPr{q, 1});  sB = pS(pM, pPr{q, 1});
        if f == 1
            y = 100 * (pMean(pRow(EN_hit, sB)) - pMean(pRow(EN_hit, sA)));
        else
            y = pGeo(pRow(EN_loc, sA), pRow(EN_loc, sB));
        end
        pLine(pX, y(pI), pPr{q, 2:6});
    end
    pAx(); xlabel(pXl);
    if f == 1
        title('True-lobe hit rate: movable minus nested'); ylabel('Hit-rate gain (points)'); ylim([0 60]);
        legend('Location', 'southeast');
    else
        plot([0 20], [0 0], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        title('Position error: movable vs nested (below 0: movable better)'); ylabel('Position error change (%)');
        ylim([-100 0]);
    end
end

% ---------------- P6 / P7: detection (Cross-AP List, with CAPE) ----------------
pA2 = {pN, cB, '-', 's', cB, 'Nested'; pM, cR, '-', 'o', cR, 'Movable'};
figure('Name', 'P6 Detection: BER and reliability (with CAPE)', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:2
        X = pBer(4, pS(pA2{q, 1}, 1));
        if f == 1, y = mean(X, 2); else, y = 100 * mean(X <= 1e-3, 2); end
        pLine(pX, y(pI), pA2{q, 2:6});
    end
    pAx(); xlabel(pXl);
    if f == 1
        set(gca, 'YScale', 'log'); title('BER (Cross-AP List, with CAPE)'); ylabel('BER'); legend('Location', 'southwest');
    else
        title('Reliability: drops with BER \leq 10^{-3}'); ylabel('Drops (%)'); ylim([0 100]);
    end
end
[~, pI12] = min(abs(SNR_dB - 12));
figure('Name', 'P7 Detection: sum SE and per-drop BER (with CAPE)', 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hold on;
for q = 1:2
    y = mean(pSE(pS(pA2{q, 1}, 1)), 2);  pLine(pX, y(pI), pA2{q, 2:6});
end
pAx(); xlabel(pXl); ylabel('Sum SE (bit/s/Hz)'); title('Sum spectral efficiency (with CAPE)'); legend('Location', 'northwest');
subplot(1, 2, 2); hold on;
for q = 1:2
    X = pBer(4, pS(pA2{q, 1}, 1));  v = sort(max(X(pI12, :), 1e-6));
    stairs(v, (1:pNd) / pNd, 'Color', pA2{q, 2}, 'LineWidth', 1.8, 'DisplayName', pA2{q, 6});
end
set(gca, 'FontSize', 12, 'XScale', 'log', 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on'); xlim([1e-6 1]); ylim([0 1]);
xlabel('Per-drop BER (error-free drops at 10^{-6})'); ylabel('CDF');
title(sprintf('Per-drop BER at %d dB (CDF over %d drops)', SNR_dB(pI12), pNd));

% ---------------- P8: number of users (with CAPE) ----------------
figure('Name', 'P8 Number of users (with CAPE)', 'Position', [30 30 1150 470]);
nSK = size(EN_Kbe, 3);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:2
        j = find(Ksw_sch == pS(pA2{q, 1}, 1), 1);
        if f == 1, Q = EN_Kbe; else, Q = EN_Ksum; end
        y = mean(reshape(Q(j, :, :), numel(Ksw), nSK), 2);
        pLine(Ksw, y, pA2{q, 2:6});
    end
    set(gca, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'XTick', Ksw); xlabel('Number of users K');
    if f == 1
        set(gca, 'YScale', 'log'); ylabel('BER'); title(sprintf('BER vs number of users (%d dB)', SNR_dB(snr_ref_idx)));
        legend('Location', 'southeast');
    else
        ylabel('Sum SE (bit/s/Hz)'); title(sprintf('Sum SE vs number of users (%d dB)', SNR_dB(snr_ref_idx)));
    end
end

% ---------------- P9: CAPE gain for List SIC and for the proposed Cross-AP List ----------------
figure('Name', 'P9 CAPE gain: List SIC vs Cross-AP List', 'Position', [30 30 1150 470]);
pTt = {'Nested: BER without CAPE / with CAPE', 'Movable: BER without CAPE / with CAPE'};
pBb = [pN pM];
for f = 1:2
    subplot(1, 2, f); hold on;
    plot([0 20], [1 1], 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    y3 = mean(pBer(3, pS(pBb(f), 0)), 2) ./ mean(pBer(3, pS(pBb(f), 1)), 2);
    y4 = mean(pBer(4, pS(pBb(f), 0)), 2) ./ mean(pBer(4, pS(pBb(f), 1)), 2);
    pLine(pX, y3(pI), cK, ':', 's', 'none', 'List SIC');
    pLine(pX, y4(pI), cR, '-', 'o', cR, 'Cross-AP List (proposed)');
    pAx(); xlabel(pXl); ylabel('BER ratio'); ylim([0.95 1.2]);
    title(pTt{f});  if f == 1, legend('Location', 'northwest'); end
end

% ---------------- P10-P12: physical reasons and hardware cost (bars) ----------------
rs0 = rng;  rng(11);
pV = zeros(2, 6);  pB2 = [pN pM];
for ia = 1:2
    b = pB2(ia);  G = EN_geo{b};  pa = G.pos(:);
    if strcmp(bs_kind{b}, 'sma'), pe = G.frames(:); pV(ia, 5) = 100 * (G.F - 1) * lambda / 4; pV(ia, 6) = G.F;
    else, pe = pa; pV(ia, 5) = 0; pV(ia, 6) = 1; end
    pV(ia, 1) = max(pe) - min(pe);
    gg = round(abs(diff(sort(unique(round(pe / (lambda / 4)))))));  gu = gg(1);
    for q = gg(:).', gu = gcd(gu, q); end
    pV(ia, 2) = min(2, 4 / gu);
    pV(ia, 3) = 100 * mean(EN_nf(pS(b, 1), 1, :));
    mu = 0;
    for t = 1:300
        c = (40 + 110 * rand) * exp(1j * 2 * pi * rand);
        U = c + 25 * sqrt(rand(en_K(1), 1)) .* exp(1j * 2 * pi * rand(en_K(1), 1));
        dH = abs(U);  d3 = sqrt(dH .^ 2 + 100);  s = sin(angle(U)) .* dH ./ d3;
        rn = sqrt(d3.' .^ 2 + pa .^ 2 - 2 * pa * (d3.' .* s.'));
        Am = (d3.' ./ rn) .* exp(-1j * 2 * pi * rn / lambda);  Am = Am ./ sqrt(sum(abs(Am) .^ 2, 1));
        G2 = abs(Am' * Am) .^ 2;  mu = mu + mean(G2(~eye(en_K(1))));
    end
    pV(ia, 4) = mu / 300;
end
rng(rs0);
pBt = {'Estimation aperture', 'Grating-lobe period (2 = no grating lobes)', 'Near-field links', ...
       'User-channel correlation (lower = easier)', 'Element travel', 'Pilot frames per estimate'};
pBy = {'Aperture (m)', 'Period (sin\theta)', 'Links (%)', 'Mean |a_i^H a_j|^2', 'Travel (cm)', 'Frames'};
pBn = {'P10 Benefit: aperture and angle ambiguity', 'P11 Near-field links and user correlation', ...
       'P12 Cost: travel and frames'};
for g = 1:3
    figure('Name', pBn{g}, 'Position', [30 30 1000 420]);
    for f = 1:2
        m = 2 * (g - 1) + f;
        subplot(1, 2, f); hold on;
        patch(1 + [-0.3 0.3 0.3 -0.3], [0 0 pV(1, m) pV(1, m)], cB, 'EdgeColor', 'k');
        patch(2 + [-0.3 0.3 0.3 -0.3], [0 0 pV(2, m) pV(2, m)], cR, 'EdgeColor', 'k');
        set(gca, 'FontSize', 12, 'XTick', 1:2, 'XTickLabel', {'Nested', 'Movable'}, 'Box', 'on', 'YGrid', 'on');
        xlim([0.4 2.6]);  ylim([0 1.15 * max(max(pV(:, m)), eps)]);
        ylabel(pBy{m});  title(pBt{m});
    end
end

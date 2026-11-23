%% PLOT_CODEU_V2_FIGURES  Redraw the Code U v2 figures without rerunning the simulation.
%  V1-V9 (smooth comparison figures), U1-U3 and U12 (geometry). Set
%  fig_old = true before running this script to also draw U4-U11 and U13.
%  Uses the workspace variables after Code_U_v2.m, or loads codeU_v2_results.mat
%  if they are not there. It never clears anything. No local functions.
if ~exist('EN_be', 'var') || ~exist('EN_loc1', 'var') || ~exist('bs_kind', 'var')
    load('codeU_v2_results.mat');
end
%% ===================== [CODE U v2] FIGURES =====================
%  V1-V9: smooth comparison figures (see below). U1-U3, U12: geometry.
%  The older figures U4-U11, U13: set fig_old = true.
%  Every base scheme is run twice: WITHOUT CAPE (each AP estimates its own
%  links) and WITH CAPE (cross-AP position fusion). Colour/marker = array.
%  Dashed line + open marker = the same antennas with difference co-array
%  processing. SNR axis = received SNR per antenna at each user's strongest AP.
uNB = numel(bs_name);  uNS = numel(en_names);  uNsnr = numel(SNR_dB);
uSet = size(EN_be, 5);
[~, uSi] = min(abs(SNR_dB - 10));
uDet = 4;                                      % Cross-AP List detector
uXl  = 'SNR [dB] (received, strongest AP)';
uCapeT = {'without CAPE (local estimation)', 'with CAPE (cross-AP fusion)'};
uCol = zeros(uNB, 3);  uMk = cell(1, uNB);  uLs = cell(1, uNB);  uMfc = cell(1, uNB);
for b = 1:uNB
    switch bs_geom{b}
        case 'ula',     c = [0.35 0.35 0.35];  mk = 'd';
        case 'nested',  c = [0.00 0.10 0.55];  mk = 's';
        case 'coprime', c = [0.50 0.10 0.60];  mk = '^';
        case 'misc',    c = [0.85 0.45 0.00];  mk = 'v';
        case 'imisc',   c = [0.55 0.30 0.10];  mk = '>';
        case 'emisc',   c = [0.00 0.55 0.55];  mk = '<';
        case 'sma_u',   c = [0.00 0.50 0.00];  mk = 'o';
        case 'sma_s',   c = [0.80 0.00 0.00];  mk = 'p';
        otherwise,      c = [0.00 0.00 0.00];  mk = 'x';
    end
    uCol(b, :) = c;  uMk{b} = mk;  uLs{b} = '-';  uMfc{b} = c;
    if strcmp(bs_kind{b}, 'coarray'), uLs{b} = '--';  uMfc{b} = 'none'; end
end
uArr = find(~strcmp(bs_kind, 'coarray'));      % one entry per physical array
uSc  = zeros(uNB, 2);                          % scheme index of (base, without/with CAPE)
for b = 1:uNB
    for cp = 0:1
        uSc(b, cp + 1) = find(en_base == b & en_cape == cp, 1);
    end
end
mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
uTmp = EN_hit;   uCnt = sum(isfinite(uTmp), 4);  uTmp(~isfinite(uTmp)) = 0;
mHT = 100 * sum(uTmp, 4) ./ max(uCnt, 1);  mHT(uCnt == 0) = NaN;
uTmp = EN_rmse;  uCnt = sum(isfinite(uTmp), 4);  uTmp(~isfinite(uTmp)) = 0;
mRM = sum(uTmp, 4) ./ max(uCnt, 1);  mRM(uCnt == 0) = NaN;
uTmp = EN_loc;   mLC = nan(uNS, 1, uNsnr);
for sc = 1:uNS
    for si = 1:uNsnr
        v = uTmp(sc, 1, si, :);  v = v(isfinite(v));
        if ~isempty(v), mLC(sc, 1, si) = median(v); end
    end
end

if ~exist('fig_old', 'var'), fig_old = false; end   % true: also draw U4-U11, U13
% ---- U1: physical arrays (positions at the data frame) ----
figure('Name', 'U1 Physical arrays', 'Position', [30 30 1100 120 + 95 * numel(uArr)]);
for a = 1:numel(uArr)
    b = uArr(a);  G = EN_geo{b};  pa = G.pos(:);  ap = max(pa) - min(pa);
    subplot(numel(uArr), 1, a); hold on; box on; grid on;
    plot(pa / lambda, zeros(size(pa)), uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', uCol(b, :), 'MarkerSize', 6);
    plot([0 0], [-1 1], 'k:');
    ylim([-1 1]);  set(gca, 'YTick', []);
    title(sprintf('%s:  aperture %.2f m,  Rayleigh %.0f m,  smallest gap %.2f\\lambda', bs_name{b}, ap, ...
          2 * ap^2 / lambda, min(diff(sort(pa))) / lambda));
end
xlabel('position along the array [\lambda]');

% ---- U2: the movable array (frames and mirror-pair virtual lags) ----
uSma = find(strcmp(bs_kind, 'sma'));
if ~isempty(uSma)
    figure('Name', 'U2 Movable array: frames and virtual lags', 'Position', [30 30 520 * numel(uSma) + 80 760]);
    for a = 1:numel(uSma)
        b = uSma(a);  G = EN_geo{b};  nh = size(G.frames, 1) / 2;  F = size(G.frames, 2);
        subplot(2, numel(uSma), a); hold on; box on; grid on;
        for sg = 1:F
            plot(G.frames(1:nh, sg) / lambda, (sg - 1) * ones(nh, 1), 'b.', 'MarkerSize', 9);
            plot(G.frames(nh + 1:end, sg) / lambda, (sg - 1) * ones(nh, 1), 'r.', 'MarkerSize', 9);
        end
        plot(G.pos / lambda, G.sig_d * ones(size(G.pos)), 'ko', 'MarkerSize', 7);
        plot([0 0], [-1 F], 'k--');
        xlabel('position [\lambda]  (blue +u, red -u, circles = data frame)');  ylabel('frame \sigma');
        title(sprintf('%s: %d frames, each half moves %.2f cm', bs_name{b}, F, 100 * (F - 1) * lambda / 4));
        subplot(2, numel(uSma), numel(uSma) + a); hold on; box on; grid on;
        bar(1:numel(G.mir.cnt), G.mir.cnt, 1, 'FaceColor', uCol(b, :));
        xlabel('mirror-pair lag 2u  [\lambda/2]');  ylabel('pairs per lag');
        title(sprintf('virtual array: lags 1..%d, %d holes', numel(G.mir.cnt), sum(G.mir.cnt == 0)));
    end
end

% ---- U3: far-field beampattern at the data frame ----
figure('Name', 'U3 Beampattern at the data frame', 'Position', [30 30 1100 170 * ceil(numel(uArr) / 2) + 80]);
sg = linspace(-1, 1, 4001);
for a = 1:numel(uArr)
    b = uArr(a);  pa = EN_geo{b}.pos(:).';
    B = abs(exp(1j * 2 * pi * (sg.') * pa / lambda) * ones(numel(pa), 1)) / numel(pa);
    subplot(ceil(numel(uArr) / 2), 2, a);
    plot(sg, 20 * log10(max(B, 1e-4)), 'Color', uCol(b, :), 'LineWidth', 1.1);
    grid on; box on; ylim([-40 0]); xlabel('sin(\theta)'); ylabel('dB'); title(bs_name{b});
end

% ---- U12: geometry and hardware metrics ----
uM = zeros(numel(uArr), 6);
for a = 1:numel(uArr)
    b = uArr(a);  G = EN_geo{b};  pa = G.pos(:);
    Dm = round((pa - pa.') / (lambda / 2));
    U = zeros(numel(pa));  nz = (Dm ~= 0) & (abs(Dm) <= 100);
    U(nz) = 1 ./ abs(Dm(nz));  U(Dm == 0) = 1;
    if strcmp(bs_kind{b}, 'sma')
        vl = sum(G.mir.cnt > 0) * lambda / 2;               % mirror-pair virtual aperture
    else
        lg = unique(round((pa - pa.') / G.unit));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
        vl = (max(diff(bk)) - 1) * G.unit;                  % hole-free difference co-array span
    end
    ap = max(pa) - min(pa);
    uM(a, :) = [ap, 2 * ap^2 / lambda, 100 * mean(EN_nf(uSc(b, 2), 1, :)), sum(Dm(:) == 1), ...
                norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps), vl];
end
figure('Name', 'U12 Geometry and hardware metrics', 'Position', [30 30 1250 720]);
uT12 = {'aperture [m]', 'Rayleigh distance [m]', 'near-field links [%]', 'w(1): pairs \lambda/2 apart', ...
        'coupling leakage', 'virtual (hole-free) aperture [m]'};
for m = 1:6
    subplot(2, 3, m); hold on; box on; grid on;
    for a = 1:numel(uArr)
        bar(a, uM(a, m), 0.6, 'FaceColor', uCol(uArr(a), :));
    end
    set(gca, 'XTick', 1:numel(uArr), 'XTickLabel', bs_name(uArr));
    if exist('xtickangle', 'file') || exist('xtickangle', 'builtin'), xtickangle(35); end
    title(uT12{m});
end

if fig_old
    % ---- U4-U6: BER, sum SE and NMSE vs SNR (without / with CAPE) ----
    uFig = {'U4 BER vs SNR', 'U5 Sum SE vs SNR', 'U6 Channel NMSE vs SNR'};
    uYl  = {'BER (Cross-AP List)', 'sum SE [bit/s/Hz]', 'channel NMSE'};
    uLoc = {'southwest', 'northwest', 'southwest'};
    for f = 1:3
        figure('Name', uFig{f}, 'Position', [30 30 1200 500]);
        for cp = 0:1
            subplot(1, 2, cp + 1); hold on; box on; grid on;
            for b = 1:uNB
                sc = uSc(b, cp + 1);
                if f == 1
                    y = max(squeeze(mBE(uDet, sc, 1, :)), 1e-6);
                elseif f == 2
                    y = squeeze(mSU(uDet, sc, 1, :));
                else
                    y = squeeze(mNM(sc, 1, :));
                end
                plot(SNR_dB, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                     'LineWidth', 1.4, 'DisplayName', bs_name{b});
            end
            if f ~= 2, set(gca, 'YScale', 'log'); end
            if f == 1, ylim([1e-6 1]); end
            xlabel(uXl); ylabel(uYl{f}); title(sprintf('%s, %s', en_caseName{1}, uCapeT{cp + 1}));
            legend('Location', uLoc{f});
        end
    end

    % ---- U7: angle acquisition of the near-field links (hit rate and RMSE) ----
    figure('Name', 'U7 Angle acquisition of near-field links', 'Position', [30 30 1200 760]);
    for row = 1:2
        for cp = 0:1
            subplot(2, 2, (row - 1) * 2 + cp + 1); hold on; box on; grid on;
            for b = 1:uNB
                sc = uSc(b, cp + 1);
                if row == 1, y = squeeze(mHT(sc, 1, :)); else, y = 1e4 * squeeze(mRM(sc, 1, :)); end
                plot(SNR_dB, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                     'LineWidth', 1.4, 'DisplayName', bs_name{b});
            end
            xlabel(uXl);
            if row == 1
                ylabel(sprintf('hit rate [%%]  (|\\Delta sin\\theta| < %.3f)', tolS));  ylim([0 100]);
            else
                ylabel('RMSE of sin\theta on hits [\times 10^{-4}]');  set(gca, 'YScale', 'log');
            end
            title(uCapeT{cp + 1});  legend('Location', 'best');
        end
    end

    % ---- U8: localisation error (CAPE only: positions are estimated by the fusion) ----
    figure('Name', 'U8 Localisation error (CAPE)', 'Position', [30 30 760 500]);
    hold on; box on; grid on;
    for b = 1:uNB
        sc = uSc(b, 2);
        plot(SNR_dB, squeeze(mLC(sc, 1, :)), [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
             'LineWidth', 1.4, 'DisplayName', bs_name{b});
    end
    set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
    title(sprintf('%s, with CAPE', en_caseName{1})); legend('Location', 'northeast');

    % ---- U9: what CAPE adds, at 10 dB ----
    figure('Name', 'U9 CAPE gain at 10 dB', 'Position', [30 30 1200 480]);
    uV = {zeros(uNB, 2), zeros(uNB, 2)};
    for b = 1:uNB
        for cp = 0:1
            sc = uSc(b, cp + 1);
            uV{1}(b, cp + 1) = max(mBE(uDet, sc, 1, uSi), 1e-6);
            uV{2}(b, cp + 1) = mNM(sc, 1, uSi);
        end
    end
    uYl9 = {'BER (Cross-AP List)', 'channel NMSE'};
    for m = 1:2
        subplot(1, 2, m); hold on; box on; grid on;
        bar(1:uNB, uV{m}, 'grouped');
        set(gca, 'YScale', 'log', 'XTick', 1:uNB, 'XTickLabel', bs_name);
        if exist('xtickangle', 'file') || exist('xtickangle', 'builtin'), xtickangle(35); end
        ylabel(uYl9{m});  title(sprintf('%s at %d dB', uYl9{m}, SNR_dB(uSi)));
        legend({'without CAPE', 'with CAPE'}, 'Location', 'best');
    end

    % ---- U10: BER and sum SE vs number of users K at 10 dB (CAPE schemes) ----
    mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);
    figure('Name', 'U10 BER and sum SE vs number of users', 'Position', [30 30 1200 480]);
    for f = 1:2
        subplot(1, 2, f); hold on; box on; grid on;
        for j = 1:numel(Ksw_sch)
            b = en_base(Ksw_sch(j));
            if f == 1, y = max(mKbe(j, :), 1e-6); else, y = mKsum(j, :); end
            plot(Ksw, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                 'LineWidth', 1.4, 'DisplayName', bs_name{b});
        end
        xlabel('number of users K');
        if f == 1, set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)'); else, ylabel('sum SE [bit/s/Hz]'); end
        title(sprintf('with CAPE, %d dB, N = %d per AP', SNR_dB(uSi), en_N0(1)));  legend('Location', 'northwest');
    end

    % ---- U11: does co-array processing help? (paired: same antennas, same channels) ----
    uPairs = [];  uPairNm = {};
    for g = {'nested', 'coprime'}
        bp = find(strcmp(bs_geom, g{1}) & strcmp(bs_kind, 'plain'));
        bc = find(strcmp(bs_geom, g{1}) & strcmp(bs_kind, 'coarray'));
        if ~isempty(bp) && ~isempty(bc), uPairs = [uPairs; bp bc]; uPairNm{end + 1} = bs_name{bp}; end %#ok<AGROW>
    end
    if ~isempty(uPairs)
        figure('Name', 'U11 Effect of co-array processing (paired)', 'Position', [30 30 1200 760]);
        uYl11 = {'BER change with co-array [%]  (<0 = better)', 'sum-SE change with co-array [%]  (>0 = better)'};
        for row = 1:2
            for cp = 0:1
                subplot(2, 2, (row - 1) * 2 + cp + 1); hold on; box on; grid on;
                plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'HandleVisibility', 'off');
                for pr = 1:size(uPairs, 1)
                    a = uSc(uPairs(pr, 1), cp + 1);  c2 = uSc(uPairs(pr, 2), cp + 1);
                    if row == 1, X = EN_be(uDet, a, 1, :, :);  Y = EN_be(uDet, c2, 1, :, :);
                    else,        X = EN_sum(uDet, a, 1, :, :); Y = EN_sum(uDet, c2, 1, :, :); end
                    X = reshape(X, uNsnr, uSet);  Y = reshape(Y, uNsnr, uSet);  Dd = Y - X;
                    ref = max(abs(mean(X, 2)), eps);
                    hEb = errorbar(SNR_dB, 100 * mean(Dd, 2) ./ ref, 100 * std(Dd, 0, 2) / sqrt(uSet) ./ ref, ...
                                   ['-' uMk{uPairs(pr, 1)}]);
                    set(hEb, 'Color', uCol(uPairs(pr, 1), :), 'LineWidth', 1.3, 'DisplayName', uPairNm{pr});
                end
                xlabel(uXl); ylabel(uYl11{row}); title(uCapeT{cp + 1}); legend('Location', 'best');
            end
        end
    end

    % ---- U13: the four detectors, with CAPE, for every array ----
    uDn = {'Linear MMSE', 'Hard SIC', 'List SIC', 'Cross-AP List'};  uDs = {':', '-.', '--', '-'};
    figure('Name', 'U13 BER of the four detectors (with CAPE)', 'Position', [30 30 1100 200 * ceil(numel(uArr) / 2) + 80]);
    for a = 1:numel(uArr)
        b = uArr(a);  sc = uSc(b, 2);
        subplot(ceil(numel(uArr) / 2), 2, a); hold on; box on; grid on;
        for d = 1:4
            plot(SNR_dB, max(squeeze(mBE(d, sc, 1, :)), 1e-6), [uDs{d} uMk{b}], 'Color', uCol(b, :), ...
                 'MarkerFaceColor', uCol(b, :), 'LineWidth', 1.2, 'DisplayName', uDn{d});
        end
        set(gca, 'YScale', 'log'); ylim([1e-6 1]); xlabel(uXl); ylabel('BER'); title(bs_name{b});
        legend('Location', 'southwest');
    end

    % ---- console summary at 10 dB ----
    fprintf('\n=== [CODE U] %s at %d dB, mean +- std. error over %d drops ===\n', en_caseName{1}, SNR_dB(uSi), uSet);
    fprintf('%-38s %6s %9s %21s %15s %7s %9s %8s\n', 'scheme', 'NF[%]', 'NMSE', 'BER (Cross-AP List)', 'sum SE', 'hit[%]', 'RMSE e-4', 'loc[m]');
    for sc = 1:uNS
        be = reshape(EN_be(uDet, sc, 1, uSi, :), [], 1);  su = reshape(EN_sum(uDet, sc, 1, uSi, :), [], 1);
        fprintf('%-38s %6.0f %9.2e %10.3e +- %7.1e %7.2f +- %4.2f %7.1f %9.1f %8.3f\n', en_names{sc}, ...
                100 * mean(EN_nf(sc, 1, :)), mNM(sc, 1, uSi), mean(be), std(be) / sqrt(uSet), ...
                mean(su), std(su) / sqrt(uSet), mHT(sc, 1, uSi), 1e4 * mRM(sc, 1, uSi), mLC(sc, 1, uSi));
    end
end

% ======================================================================
%  [CODE U v2] FIGURES V1-V9. Every curve is the mean over the drops with
%  +- one standard error; every scheme sees the same drops, channels,
%  pilots, noise and symbols (paired). With CAPE: solid, filled marker.
%  Without CAPE: dotted, open marker. Co-array processing: dashed.
% ======================================================================
vS12  = find(SNR_dB >= 12, 1);  if isempty(vS12), vS12 = uSi; end    % summary SNR (12 dB)
vSEm  = SNR_dB <= 16;                           % sum SE shown up to 16 dB (measured SE saturates above)
vFl   = 1 / (max(en_K) * nRealE * nSymE * log2(modOrder) * uSet);   % one bit error over all drops
vFin  = @(v) [v(isfinite(v)), nan(1, ~any(isfinite(v)))];
vMS   = @(X) deal(arrayfun(@(i) mean(vFin(X(i, :))), (1:size(X, 1)).'), ...
                  arrayfun(@(i) std(vFin(X(i, :))) / sqrt(max(nnz(isfinite(X(i, :))), 1)), (1:size(X, 1)).'));
vMed  = @(X) arrayfun(@(i) median(vFin(X(i, :))), (1:size(X, 1)).');
vRow  = @(Q, sc) reshape(Q(sc, 1, :, :), uNsnr, uSet);
vRowB = @(d, sc) reshape(EN_be(d, sc, 1, :, :), uNsnr, uSet);
vRowS = @(d, sc) reshape(EN_sum(d, sc, 1, :, :), uNsnr, uSet);
vPhys = find(~strcmp(bs_kind, 'coarray'));     % physical processing: one per array
vSpar = vPhys(~strcmp(bs_geom(vPhys), 'ula')); % sparse arrays (incl. movable)
vSma  = find(strcmp(bs_kind, 'sma'));
bCop  = find(strcmp(bs_geom, 'coprime') & strcmp(bs_kind, 'plain'), 1);   % best static array of Code T
vEB   = @(h, b, ls, mf, nm) set(h, 'LineStyle', ls, 'Marker', uMk{b}, 'Color', uCol(b, :), ...
                                'MarkerFaceColor', mf, 'LineWidth', 1.4, 'DisplayName', nm);

% ---------------- V1: detection of all arrays (with CAPE) ----------------
figure('Name', 'V1 Which array: detection (Cross-AP List, with CAPE)', 'Position', [30 30 1150 820]);
for b = vPhys
    sc = uSc(b, 2);  X = vRowB(uDet, sc);  Y = vRowS(uDet, sc);
    m = max(mean(X, 2), vFl);  s = std(X, 0, 2) / sqrt(uSet);
    subplot(2, 2, 1); hold on;  vEB(errorbar(SNR_dB, m, min(s, 0.9 * m), s), b, '-', uCol(b, :), bs_name{b});
    r = mean(X <= 1e-3, 2);  rs = sqrt(r .* (1 - r) / uSet);
    subplot(2, 2, 2); hold on;  vEB(errorbar(SNR_dB, 100 * r, 100 * rs), b, '-', uCol(b, :), bs_name{b});
    subplot(2, 2, 3); hold on;
    vEB(errorbar(SNR_dB(vSEm), mean(Y(vSEm, :), 2), std(Y(vSEm, :), 0, 2) / sqrt(uSet)), b, '-', uCol(b, :), bs_name{b});
    v = sort(max(X(vS12, :), 1e-6));
    subplot(2, 2, 4); hold on;  stairs(v, (1:uSet) / uSet, 'Color', uCol(b, :), 'LineWidth', 1.6, 'DisplayName', bs_name{b});
end
subplot(2, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('BER (Cross-AP List)');
title('Mean BER \pm s.e.'); legend('Location', 'southwest');
subplot(2, 2, 2); box on; grid on; xlabel(uXl); ylabel('drops with BER \leq 10^{-3} [%]'); ylim([0 100]);
title('Reliability: share of drops with BER \leq 10^{-3}'); legend('Location', 'northwest');
subplot(2, 2, 3); box on; grid on; xlabel(uXl); ylabel('sum SE [bit/s/Hz]');
title('Sum SE \pm s.e. (up to 16 dB: above it the measured SE saturates)'); legend('Location', 'northwest');
subplot(2, 2, 4); box on; grid on; set(gca, 'XScale', 'log'); xlim([1e-6 1]); ylim([0 1]);
xlabel('per-drop BER (error-free drops shown at 10^{-6})'); ylabel('CDF over the drops');
title(sprintf('Distribution of the per-drop BER at %d dB', SNR_dB(vS12))); legend('Location', 'northwest');

% ---------------- V2: angle and localisation (with CAPE) ----------------
figure('Name', 'V2 Which array: angle and localisation (with CAPE)', 'Position', [30 30 1500 470]);
for b = vSpar
    sc = uSc(b, 2);
    [m, s] = vMS(vRow(EN_hit, sc));
    subplot(1, 3, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, '-', uCol(b, :), bs_name{b});
    subplot(1, 3, 2); hold on;
    plot(SNR_dB, vMed(vRow(EN_loc, sc)), '-', 'Marker', uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', uCol(b, :), ...
         'LineWidth', 1.4, 'DisplayName', bs_name{b});
    [m, s] = vMS(vRow(EN_loc1, sc));
    subplot(1, 3, 3); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, '-', uCol(b, :), bs_name{b});
end
subplot(1, 3, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
title(sprintf('Angle: true lobe found (|\\Delta sin\\theta| < %.3f)', tolS)); legend('Location', 'southeast');
subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
title('Position error (median over the drops)'); legend('Location', 'northeast');
subplot(1, 3, 3); box on; grid on; xlabel(uXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'southeast');

% ---------------- V3: the movable array resolves the angle WITHOUT cross-AP help ----------------
figure('Name', 'V3 Movable array: angle from ONE AP (no CAPE)', 'Position', [30 30 1150 470]);
for b = vSpar
    sc = uSc(b, 1);
    [m, s] = vMS(vRow(EN_hit, sc));
    subplot(1, 2, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ':', 'none', bs_name{b});
    R = vRow(EN_rmse, sc);  rm = arrayfun(@(i) sqrt(mean(vFin(R(i, :)) .^ 2)), (1:uNsnr).');
    subplot(1, 2, 2); hold on;
    plot(SNR_dB, rm, ':', 'Marker', uMk{b}, 'Color', uCol(b, :), 'LineWidth', 1.4, 'DisplayName', bs_name{b});
end
subplot(1, 2, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
title('Each AP alone (no cross-AP fusion): true lobe found'); legend('Location', 'southeast');
subplot(1, 2, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('angle RMSE of the hits [sin\theta]');
title('Angle accuracy when the true lobe is found (no CAPE)'); legend('Location', 'northeast');

% ---------------- V4: load (number of users), with CAPE ----------------
nSK = size(EN_Kbe, 3);
figure('Name', 'V4 Which array vs number of users (with CAPE)', 'Position', [30 30 1150 470]);
for j = 1:numel(Ksw_sch)
    sc = Ksw_sch(j);  b = en_base(sc);
    if strcmp(bs_kind{b}, 'coarray'), continue; end
    Xb = reshape(EN_Kbe(j, :, :), numel(Ksw), nSK);  Xs = reshape(EN_Ksum(j, :, :), numel(Ksw), nSK);
    m = max(mean(Xb, 2), vFl);  s = std(Xb, 0, 2) / sqrt(nSK);
    subplot(1, 2, 1); hold on;  vEB(errorbar(Ksw, m, min(s, 0.9 * m), s), b, '-', uCol(b, :), bs_name{b});
    subplot(1, 2, 2); hold on;  vEB(errorbar(Ksw, mean(Xs, 2), std(Xs, 0, 2) / sqrt(nSK)), b, '-', uCol(b, :), bs_name{b});
end
subplot(1, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
title(sprintf('BER vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southeast');
subplot(1, 2, 2); box on; grid on; xlabel('number of users K'); ylabel('sum SE [bit/s/Hz]');
title(sprintf('Sum SE vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southeast');

% ---------------- V5: physical and hardware metrics (benefit vs cost) ----------------
rs0 = rng;  rng(11);
nA = numel(vPhys);  vM = zeros(nA, 6);
for ia = 1:nA
    b = vPhys(ia);  G = EN_geo{b};  pa = G.pos(:);
    if strcmp(bs_kind{b}, 'sma')
        pe = G.frames(:);  vM(ia, 5) = 100 * (G.F - 1) * lambda / 4;  vM(ia, 6) = G.F;
    else
        pe = pa;  vM(ia, 5) = 0;  vM(ia, 6) = 1;
    end
    vM(ia, 1) = max(pe) - min(pe);                                 % estimation aperture
    vM(ia, 2) = 100 * mean(EN_nf(uSc(b, 2), 1, :));                % near-field links
    gg = round(abs(diff(sort(unique(round(pe / (lambda / 4))))))); % grid in lambda/4 units
    gu = gg(1);  for q = gg(:).', gu = gcd(gu, q); end
    vM(ia, 3) = min(2, 2 / (gu / 2));                              % grating-lobe period of the estimation aperture
    mu = 0;
    for t = 1:300                                                  % 16 users in a 25 m hotspot, 40-150 m away
        c = (40 + 110 * rand) * exp(1j * 2 * pi * rand);
        U = c + 25 * sqrt(rand(en_K(1), 1)) .* exp(1j * 2 * pi * rand(en_K(1), 1));
        dH = abs(U);  d3 = sqrt(dH .^ 2 + 100);  s = sin(angle(U)) .* dH ./ d3;
        rn = sqrt(d3.' .^ 2 + pa .^ 2 - 2 * pa * (d3.' .* s.'));
        Am = (d3.' ./ rn) .* exp(-1j * 2 * pi * rn / lambda);  Am = Am ./ sqrt(sum(abs(Am) .^ 2, 1));
        G2 = abs(Am' * Am) .^ 2;  mu = mu + mean(G2(~eye(en_K(1))));
    end
    vM(ia, 4) = mu / 300;
end
rng(rs0);
vShort = strrep(strrep(strrep(bs_name(vPhys), 'Movable (uniform half)', 'Mov. unif.'), 'Movable (sparse half)', 'Mov. sparse'), 'Dense ULA', 'ULA');
figure('Name', 'V5 Why the arrays differ: benefit and hardware cost', 'Position', [30 30 1500 760]);
vY = {'estimation aperture [m]', 'near-field links [%]', 'grating-lobe period in sin\theta (2 = none)', ...
      'mean |a_i^H a_j|^2 between hotspot users', 'travel of the outer element [cm]', 'pilot frames per estimate'};
vT = {'Aperture used for estimation', 'Links in the near field', 'Angle ambiguity of the estimation aperture', ...
      'User-channel correlation at the data positions', 'Mechanical cost: travel', 'Latency cost: frames (user must stay still)'};
for mI = 1:6
    subplot(2, 3, mI); hold on; box on; grid on;
    for ia = 1:nA, bar(ia, vM(ia, mI), 0.6, 'FaceColor', uCol(vPhys(ia), :)); end
    set(gca, 'XTick', 1:nA, 'XTickLabel', vShort); xlim([0.4 nA + 0.6]); ylabel(vY{mI}); title(vT{mI});
end

% ---------------- V6: movable array vs the best static array (paired) ----------------
figure('Name', 'V6 Movable array vs coprime (paired, same drops)', 'Position', [30 30 1150 820]);
vT6 = {'BER(coprime) / BER(movable)  (>1: movable better)', 'true-lobe hit rate gain [points]', ...
       'position error change [%]  (<0: movable better)', 'NMSE change [%]  (<0: movable better)'};
for f = 1:4
    subplot(2, 2, f); hold on; box on; grid on;
    if f == 1, yr = ones(size(SNR_dB)); else, yr = zeros(size(SNR_dB)); end
    plot(SNR_dB, yr, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for b = vSma(:).'
        for cp = [1 0]
            sA = uSc(bCop, cp + 1);  sB = uSc(b, cp + 1);
            if f == 1
                X = vRowB(uDet, sA);  Y = vRowB(uDet, sB);  mx = mean(X, 2);  my = mean(Y, 2);
                R = mx ./ max(my, vFl);
                cxy = arrayfun(@(i) sum((X(i, :) - mx(i)) .* (Y(i, :) - my(i))) / (uSet - 1), (1:uNsnr).');
                se = R .* sqrt(max(var(X, 0, 2) ./ (uSet * max(mx, vFl) .^ 2) + var(Y, 0, 2) ./ (uSet * max(my, vFl) .^ 2) ...
                     - 2 * cxy ./ (uSet * max(mx, vFl) .* max(my, vFl)), 0));
                h = errorbar(SNR_dB, R, se);
            elseif f == 2
                [m, s] = vMS(100 * (vRow(EN_hit, sB) - vRow(EN_hit, sA)));  h = errorbar(SNR_dB, m, s);
            else
                if f == 3, X = vRow(EN_loc, sA);  Y = vRow(EN_loc, sB); else, X = vRow(EN_nmse, sA);  Y = vRow(EN_nmse, sB); end
                mu = nan(uNsnr, 1);  lo = mu;  hi = mu;
                for si = 1:uNsnr
                    ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
                    if nnz(ok) < 2, continue; end
                    rr = log(Y(si, ok) ./ X(si, ok));  mm = mean(rr);  ss = std(rr) / sqrt(nnz(ok));
                    mu(si) = 100 * (exp(mm) - 1);  lo(si) = 100 * (exp(mm - ss) - 1);  hi(si) = 100 * (exp(mm + ss) - 1);
                end
                h = errorbar(SNR_dB, mu, mu - lo, hi - mu);
            end
            if cp == 1, vEB(h, b, '-', uCol(b, :), [bs_name{b} ', CAPE']);
            else, vEB(h, b, ':', 'none', [bs_name{b} ', no CAPE']); end
        end
    end
    xlabel(uXl); ylabel(vT6{f}); legend('Location', 'best');
end

% ---------------- V7: co-array, paired gains (same antennas) ----------------
vCoP = [];  vCoN = {};
for b = find(strcmp(bs_kind, 'coarray'))
    bp = find(strcmp(bs_geom, bs_geom{b}) & strcmp(bs_kind, 'plain'), 1);
    if ~isempty(bp), vCoP = [vCoP; bp b]; vCoN{end + 1} = bs_name{bp}; end %#ok<AGROW>
end
figure('Name', 'V7 Co-array: paired gains vs SNR (same antennas)', 'Position', [30 30 1150 820]);
vT7 = {'position error change [%]  (<0 = co-array better)', 'true-lobe hit rate change [points]', ...
       'NMSE change [%]  (<0 = better)', 'BER change [%] (Cross-AP List)  (no loss = around 0)'};
for f = 1:4
    subplot(2, 2, f); hold on; box on; grid on;
    plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for pr = 1:size(vCoP, 1)
        for cp = [1 0]
            sA = uSc(vCoP(pr, 1), cp + 1);  sB = uSc(vCoP(pr, 2), cp + 1);  b = vCoP(pr, 1);
            if f == 1 || f == 3
                if f == 1, X = vRow(EN_loc, sA);  Y = vRow(EN_loc, sB); else, X = vRow(EN_nmse, sA);  Y = vRow(EN_nmse, sB); end
                mu = nan(uNsnr, 1);  lo = mu;  hi = mu;
                for si = 1:uNsnr
                    ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
                    if nnz(ok) < 2, continue; end
                    rr = log(Y(si, ok) ./ X(si, ok));  mm = mean(rr);  ss = std(rr) / sqrt(nnz(ok));
                    mu(si) = 100 * (exp(mm) - 1);  lo(si) = 100 * (exp(mm - ss) - 1);  hi(si) = 100 * (exp(mm + ss) - 1);
                end
                h = errorbar(SNR_dB, mu, mu - lo, hi - mu);
            elseif f == 2
                [m, s] = vMS(100 * (vRow(EN_hit, sB) - vRow(EN_hit, sA)));  h = errorbar(SNR_dB, m, s);
            else
                X = vRowB(uDet, sA);  Y = vRowB(uDet, sB);  ref = max(mean(X, 2), eps);
                h = errorbar(SNR_dB, 100 * mean(Y - X, 2) ./ ref, 100 * std(Y - X, 0, 2) / sqrt(uSet) ./ ref);
            end
            if cp == 1, vEB(h, b, '-', uCol(b, :), [vCoN{pr} ', CAPE']);
            else, vEB(h, b, ':', 'none', [vCoN{pr} ', no CAPE']); end
        end
    end
    xlabel(uXl); ylabel(vT7{f}); legend('Location', 'best');
end

% ---------------- V8: CAPE vs no CAPE (angle and localisation) ----------------
figure('Name', 'V8 CAPE: angle and localisation (with vs without CAPE)', 'Position', [30 30 1500 470]);
for b = vSpar
    for cp = [1 0]
        sc = uSc(b, cp + 1);
        if cp == 1, ls = '-';  mf = uCol(b, :);  nmc = ', CAPE'; else, ls = ':';  mf = 'none';  nmc = ', no CAPE'; end
        [m, s] = vMS(vRow(EN_hit, sc));
        subplot(1, 3, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ls, mf, [bs_name{b} nmc]);
        subplot(1, 3, 2); hold on;
        plot(SNR_dB, vMed(vRow(EN_loc, sc)), ls, 'Marker', uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', mf, ...
             'LineWidth', 1.4, 'DisplayName', [bs_name{b} nmc]);
        [m, s] = vMS(vRow(EN_loc1, sc));
        subplot(1, 3, 3); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ls, mf, [bs_name{b} nmc]);
    end
end
subplot(1, 3, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
title('Angle: true lobe found'); legend('Location', 'southeast');
subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
title('Position: CAPE vs local triangulation'); legend('Location', 'northeast');
subplot(1, 3, 3); box on; grid on; xlabel(uXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'northwest');

% ---------------- V9: CAPE gain per detector (List SIC vs proposed Cross-AP List) ----------------
vG = [vSpar(:).'];
figure('Name', 'V9 CAPE gain: List SIC vs the proposed Cross-AP List', 'Position', [30 30 1500 420 * ceil(numel(vG) / 3)]);
vDs = {'--', '-'};  vDl = {'List SIC', 'Cross-AP List (proposed)'};
for ig = 1:numel(vG)
    b = vG(ig);  sc = uSc(b, 2);  sn = uSc(b, 1);
    subplot(ceil(numel(vG) / 3), 3, ig); hold on; box on; grid on;
    plot(SNR_dB, ones(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for id = 1:2
        X = vRowB(id + 2, sn);  Y = vRowB(id + 2, sc);  mx = mean(X, 2);  my = mean(Y, 2);  ok = mx > 0 & my > 0;
        R = nan(uNsnr, 1);  se = R;
        cxy = arrayfun(@(i) sum((X(i, :) - mx(i)) .* (Y(i, :) - my(i))) / (uSet - 1), (1:uNsnr).');
        R(ok) = mx(ok) ./ my(ok);
        se(ok) = R(ok) .* sqrt(max(var(X(ok, :), 0, 2) ./ (uSet * mx(ok) .^ 2) + var(Y(ok, :), 0, 2) ./ (uSet * my(ok) .^ 2) ...
                 - 2 * cxy(ok) ./ (uSet * mx(ok) .* my(ok)), 0));
        mf = uCol(b, :);  if id == 1, mf = 'none'; end
        vEB(errorbar(SNR_dB, R, se), b, vDs{id}, mf, vDl{id});
    end
    xlabel(uXl); ylabel('BER(no CAPE) / BER(CAPE)  (>1: CAPE helps)'); title(bs_name{b}); legend('Location', 'best');
end

% ---------------- console summary (V figures) ----------------
fprintf('\n=== [CODE U v2] SUMMARY AT %d dB (mean over %d drops) ===\n', SNR_dB(vS12), uSet);
fprintf('%-26s %10s %10s %8s %10s %9s %8s %10s\n', 'array (CAPE)', 'BER', 'BER<=1e-3', 'sum SE', 'pos. [m]', 'within1m', 'hit', 'hit noCAPE');
for b = vPhys
    sc = uSc(b, 2);  X = vRowB(uDet, sc);  Y = vRowS(uDet, sc);
    md = vMed(vRow(EN_loc, sc));  [w1, ~] = vMS(vRow(EN_loc1, sc));  [ht, ~] = vMS(vRow(EN_hit, sc));
    [h0, ~] = vMS(vRow(EN_hit, uSc(b, 1)));
    fprintf('%-26s %10.2e %9.0f%% %8.1f %10.2f %8.0f%% %7.0f%% %9.0f%%\n', bs_name{b}, mean(X(vS12, :)), ...
            100 * mean(X(vS12, :) <= 1e-3), mean(Y(vS12, :)), md(vS12), 100 * w1(vS12), 100 * ht(vS12), 100 * h0(vS12));
end

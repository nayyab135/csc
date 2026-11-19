%% PLOT_CODEU_FIGURES  Redraw Code U figures U1-U13 without rerunning the simulation.
%  Uses the workspace variables after Code_U_MovableSparse.m, or loads
%  codeU_movable_results.mat if they are not there. It never clears anything.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('bs_kind', 'var') || ~exist('EN_hit', 'var')
    load('codeU_movable_results.mat');
end
%% ===================== [CODE U] FIGURES U1-U13 =====================
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

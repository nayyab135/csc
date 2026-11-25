% ======================================================================
%  plot_codeU_coarray_figures.m : two figures on the DIFFERENCE CO-ARRAY
%  (Nested physical, Nested + co-array, Movable), from codeU_v2_results.mat.
%    Fig 6  Why: physical antennas and virtual (co-array) lags | angle spectrum
%    Fig 7  In the system, with CAPE: users within 1 m | median position error
%  It never clears the workspace. No local functions.
% ======================================================================
if ~exist('EN_loc1', 'var') || ~exist('EN_geo', 'var') || ~exist('bs_kind', 'var')
    load('codeU_v2_results.mat');
end
cN = find(strcmp(bs_name, 'Nested'), 1);
cC = find(strcmp(bs_name, 'Nested + co-array'), 1);
cM = find(strcmp(bs_name, 'Movable (sparse half)'), 1);
cB = [0.00 0.45 0.85];  cDB = [0.00 0.15 0.45];  cR = [0.85 0.10 0.10];
lamC = Dtgt / 28;                                   % lambda (Dtgt = 28 lambda = 2.8 m)

% ---------------- geometry: physical positions and virtual lags ----------------
gN = EN_geo{cN};  if iscell(gN), gN = gN{1}; end
gM = EN_geo{cM};  if iscell(gM), gM = gM{1}; end
xN = sort(gN.pos(:)).';  xN = xN - xN(1);           % nested, physical (m)
lN = unique(round((xN - xN.') / gN.unit));  lN = lN(lN >= 0).' * gN.unit;   % nested co-array lags (m)
xM = sort(gM.pos(:)).';  xM = xM - xM(1);           % movable, data frame (m)
lM = (1:gM.mir.Q) * lamC / 2;                       % movable mirror-pair lags 2u (m)

% ---------------- Fig 6: physical antennas, virtual lags, angle spectrum ----------------
figure('Name', 'Fig 6 Difference co-array: nested vs movable', 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hold on;
plot(xN, 4 * ones(size(xN)), 's', 'Color', cB, 'MarkerFaceColor', cB, 'MarkerSize', 8);
plot(lN, 3 * ones(size(lN)), 's', 'Color', cDB, 'MarkerFaceColor', 'none', 'MarkerSize', 6);
plot(xM, 2 * ones(size(xM)), 'o', 'Color', cR, 'MarkerFaceColor', cR, 'MarkerSize', 8);
plot(lM, 1 * ones(size(lM)), 'o', 'Color', cR, 'MarkerFaceColor', 'none', 'MarkerSize', 3);
text(0, 4.35, sprintf('Nested: %d antennas', numel(xN)), 'FontSize', 11);
text(0, 3.35, sprintf('Nested co-array: lags 0 to %d \\times %.1f\\lambda (no holes)', numel(lN) - 1, gN.unit / lamC), 'FontSize', 11);
text(0, 2.35, sprintf('Movable (one frame): %d antennas', numel(xM)), 'FontSize', 11);
text(0, 1.35, sprintf('Movable mirror lags: 1 to %d \\times \\lambda/2 (no holes)', numel(lM)), 'FontSize', 11);
set(gca, 'FontSize', 12, 'YTick', [], 'YLim', [0.5 4.8], 'XLim', [-0.1 3.8], 'Box', 'on', 'XGrid', 'on');
xlabel('Position or lag (m)');
title('Physical antennas and virtual (co-array) lags');
pAxL6 = gca;

subplot(1, 2, 2); hold on;
s0 = 0.2;  sg = linspace(-1, 1, 4001);              % one user at sin(theta) = 0.2
af = @(x) abs(exp(1j * 2 * pi / lamC * (sg(:) - s0) * x(:).') * ones(numel(x), 1)) / numel(x);
lNs = unique(round((xN - xN.') / gN.unit)).' * gN.unit;   % full (two-sided) nested co-array
plot(sg, af(xN), '-', 'Color', cB, 'LineWidth', 1.2, 'DisplayName', 'Nested, physical');
plot(sg, af(lNs), '--', 'Color', cDB, 'LineWidth', 1.8, 'DisplayName', 'Nested, co-array');
plot(sg, af(lM), '-', 'Color', cR, 'LineWidth', 1.8, 'DisplayName', 'Movable, mirror lags');
set(gca, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'XLim', [-1 1], 'YLim', [0 1.05]);
xlabel('sin\theta'); ylabel('Normalised spectrum');
title(sprintf('Angle spectrum, one user at sin\\theta = %.1f', s0));
legend('Location', 'southoutside', 'NumColumns', 3);
pP6 = get(gca, 'Position');  pL6 = get(pAxL6, 'Position');   % same height for both panels
set(pAxL6, 'Position', [pL6(1), pP6(2), pL6(3), pP6(4)]);

% ---------------- Fig 7: in the system, with CAPE ----------------
cI = find(SNR_dB <= 20);  cX = SNR_dB(cI);
cRow = @(X, s) reshape(X(s, 1, :, :), numel(SNR_dB), []);
cCfg = {cN, cB, '-', 's', cB, 'Nested'; cC, cDB, '--', 'd', 'none', 'Nested + co-array'; ...
        cM, cR, '-', 'o', cR, 'Movable'};
figure('Name', 'Fig 7 Difference co-array in the system (with CAPE)', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on;
    for q = 1:3
        s = find(en_base == cCfg{q, 1} & en_cape == 1, 1);
        if f == 1
            y = 100 * mean(cRow(EN_loc1, s), 2);
        else
            X = cRow(EN_loc, s);  y = nan(size(X, 1), 1);
            for i = 1:size(X, 1), v = X(i, ~isnan(X(i, :)));  if ~isempty(v), y(i) = median(v); end, end
        end
        plot(cX, y(cI), 'Color', cCfg{q, 2}, 'LineStyle', cCfg{q, 3}, 'Marker', cCfg{q, 4}, ...
             'MarkerFaceColor', cCfg{q, 5}, 'MarkerSize', 7, 'LineWidth', 1.8, 'DisplayName', cCfg{q, 6});
    end
    set(gca, 'FontSize', 12, 'XLim', [0 20], 'XTick', 0:4:20, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on');
    xlabel('SNR (dB)');
    if f == 1
        ylabel('Users within 1 m (%)'); ylim([0 100]);
        title('Users located within 1 m (with CAPE)'); legend('Location', 'northwest');
    else
        set(gca, 'YScale', 'log'); ylabel('Median position error (m)');
        title('Position error (with CAPE)');
    end
end

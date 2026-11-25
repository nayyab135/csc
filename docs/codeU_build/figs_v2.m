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
title('Position error (median over the drops)'); legend('Location', 'southwest');
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
figure('Name', 'V4 Which array vs number of users (with CAPE)', 'Position', [30 30 1150 560]);
for j = 1:numel(Ksw_sch)
    sc = Ksw_sch(j);  b = en_base(sc);
    if strcmp(bs_kind{b}, 'coarray'), continue; end
    Xb = reshape(EN_Kbe(j, :, :), numel(Ksw), nSK);  Xs = reshape(EN_Ksum(j, :, :), numel(Ksw), nSK);
    m = max(mean(Xb, 2), vFl);  s = std(Xb, 0, 2) / sqrt(nSK);
    subplot(1, 2, 1); hold on;  vEB(errorbar(Ksw, m, min(s, 0.9 * m), s), b, '-', uCol(b, :), bs_name{b});
    subplot(1, 2, 2); hold on;  vEB(errorbar(Ksw, mean(Xs, 2), std(Xs, 0, 2) / sqrt(nSK)), b, '-', uCol(b, :), bs_name{b});
end
subplot(1, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
title(sprintf('BER vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southoutside', 'NumColumns', 3);
subplot(1, 2, 2); box on; grid on; xlabel('number of users K'); ylabel('sum SE [bit/s/Hz]');
title(sprintf('Sum SE vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southoutside', 'NumColumns', 3);

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
vShort = strrep(strrep(strrep(strrep(strrep(bs_name(vPhys), 'Movable (uniform half)', 'Mov-U'), 'Movable (sparse half)', 'Mov-S'), 'Dense ULA', 'ULA'), 'Coprime', 'Copr.'), 'Nested', 'Nest.');
figure('Name', 'V5 Why the arrays differ: benefit and hardware cost', 'Position', [30 30 1500 760]);
vY = {'estimation aperture [m]', 'near-field links [%]', 'grating-lobe period in sin\theta (2 = none)', ...
      'mean |a_i^H a_j|^2 between hotspot users', 'travel of the outer element [cm]', 'pilot frames per estimate'};
vT = {'Aperture used for estimation', 'Links in the near field', 'Angle ambiguity of the estimation aperture', ...
      'User-channel correlation at the data positions', 'Mechanical cost: travel', 'Latency cost: frames (user must stay still)'};
for mI = 1:6
    subplot(2, 3, mI); hold on; box on; grid on;
    for ia = 1:nA
        patch(ia + [-0.3 0.3 0.3 -0.3], [0 0 vM(ia, mI) vM(ia, mI)], uCol(vPhys(ia), :), 'EdgeColor', 'k');
    end
    set(gca, 'XTick', 1:nA, 'XTickLabel', vShort); xlim([0.4 nA + 0.6]); ylabel(vY{mI}); title(vT{mI});
end

% ---------------- V6: movable array vs the best static array (paired) ----------------
figure('Name', 'V6 Movable array vs coprime (paired, same drops)', 'Position', [30 30 1150 950]);
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
    xlabel(uXl); ylabel(vT6{f}); legend('Location', 'southoutside', 'NumColumns', 2);
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
figure('Name', 'V8 CAPE: angle and localisation (with vs without CAPE)', 'Position', [30 30 1500 640]);
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
title('Angle: true lobe found'); legend('Location', 'southoutside', 'NumColumns', 2);
subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
title('Position: CAPE vs local triangulation'); legend('Location', 'southoutside', 'NumColumns', 2);
subplot(1, 3, 3); box on; grid on; xlabel(uXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'southoutside', 'NumColumns', 2);

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

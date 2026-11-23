% ======================================================================
%  NEW FIGURES. Every curve is the mean over the drops with +- one standard
%  error (paired drops: every scheme sees the same drops, channels, pilots,
%  noise and symbols). Colour = physical array. With CAPE: solid line,
%  filled marker. Without CAPE: dotted line, open marker.
%  A: which sparse array performs better, and why
%  B: the role of the difference co-array (same antennas, processing only)
%  C: the benefit of CAPE (cross-AP position-fused estimation)
% ======================================================================
tS12  = find(SNR_dB >= 12, 1);  if isempty(tS12), tS12 = tSi; end   % summary SNR (12 dB)
tSEmx = 16;  tSEm = SNR_dB <= tSEmx;           % sum SE shown up to 16 dB (measured SE saturates above)
tSp   = 2:4;                                   % sparse arrays (index into tArrSch): wide ULA, nested, coprime
tCoP  = [3 4; 5 6];                            % [physical, co-array] scheme pairs, with CAPE
tCoL  = tPairL;                                % the same pairs without CAPE
tRow  = @(Q, sc) reshape(Q(sc, 1, :, :), tNsnr, tSet);           % SNR x drops of a scheme
tRowB = @(d, sc) reshape(EN_be(d, sc, 1, :, :), tNsnr, tSet);    % BER of detector d
tRowS = @(d, sc) reshape(EN_sum(d, sc, 1, :, :), tNsnr, tSet);   % sum SE of detector d
tFin  = @(v) [v(isfinite(v)), nan(1, ~any(isfinite(v)))];       % finite entries (NaN if none)
tNan  = @(X) deal(arrayfun(@(i) mean(tFin(X(i, :))), (1:size(X, 1)).'), ...
                  arrayfun(@(i) std(tFin(X(i, :))) / sqrt(max(nnz(isfinite(X(i, :))), 1)), (1:size(X, 1)).'));
tMedN = @(X) arrayfun(@(i) median(tFin(X(i, :))), (1:size(X, 1)).');
tLegC = {'ULA (\lambda/2)', 'ULA (4\lambda)', 'Nested', 'Coprime'};

% ---------------- A1: detection performance of the four arrays ----------------
figure('Name', 'A1 Which array: detection (Cross-AP List, with CAPE)', 'Position', [30 30 1150 820]);
for a = 1:4
    sc = tArrSch(a);  X = tRowB(tDet, sc);  Y = tRowS(tDet, sc);
    m = mean(X, 2);  s = std(X, 0, 2) / sqrt(tSet);
    subplot(2, 2, 1); hold on;
    h = errorbar(SNR_dB, max(m, tFl), min(s, 0.9 * max(m, tFl)), s, ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
    r = mean(X <= 1e-3, 2);  rs = sqrt(r .* (1 - r) / tSet);
    subplot(2, 2, 2); hold on;
    h = errorbar(SNR_dB, 100 * r, 100 * rs, ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
    subplot(2, 2, 3); hold on;
    h = errorbar(SNR_dB(tSEm), mean(Y(tSEm, :), 2), std(Y(tSEm, :), 0, 2) / sqrt(tSet), ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
    v = sort(max(X(tS12, :), 1e-6));
    subplot(2, 2, 4); hold on;
    stairs(v, (1:tSet) / tSet, 'Color', tArrCol{a}, 'LineWidth', 1.6, 'DisplayName', tLegC{a});
end
subplot(2, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel(tXl); ylabel('BER (Cross-AP List)');
title('Mean BER \pm s.e.'); legend('Location', 'southwest');
subplot(2, 2, 2); box on; grid on; xlabel(tXl); ylabel('drops with BER \leq 10^{-3} [%]'); ylim([0 100]);
title('Reliability: share of drops with BER \leq 10^{-3}'); legend('Location', 'northwest');
subplot(2, 2, 3); box on; grid on; xlabel(tXl); ylabel('sum SE [bit/s/Hz] (Cross-AP List)');
title(sprintf('Sum SE \\pm s.e. (up to %d dB: above it the measured SE saturates)', tSEmx)); legend('Location', 'northwest');
subplot(2, 2, 4); box on; grid on; set(gca, 'XScale', 'log'); xlim([1e-6 1]); ylim([0 1]);
xlabel('per-drop BER (error-free drops shown at 10^{-6})'); ylabel('CDF over the drops');
title(sprintf('Distribution of the per-drop BER at %d dB', SNR_dB(tS12))); legend('Location', 'northwest');

% ---------------- A2: localisation and angle of the three sparse arrays ----------------
figure('Name', 'A2 Which array: localisation and angle (with CAPE)', 'Position', [30 30 1500 470]);
for a = tSp
    sc = tArrSch(a);
    subplot(1, 3, 1); hold on;
    plot(SNR_dB, tMedN(tRow(EN_loc, sc)), ['-' tMk{sc}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, ...
         'LineWidth', 1.5, 'DisplayName', tLegC{a});
    [m, s] = tNan(tRow(EN_loc1, sc));
    subplot(1, 3, 2); hold on;
    h = errorbar(SNR_dB, 100 * m, 100 * s, ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
    [m, s] = tNan(tRow(EN_hit, sc));
    subplot(1, 3, 3); hold on;
    h = errorbar(SNR_dB, 100 * m, 100 * s, ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
end
subplot(1, 3, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel(tXl); ylabel('median position error [m]');
title('Position error (median over the drops)'); legend('Location', 'northeast');
subplot(1, 3, 2); box on; grid on; xlabel(tXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'northwest');
subplot(1, 3, 3); box on; grid on; xlabel(tXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
title(sprintf('Angle: true lobe found (|\\Delta sin\\theta| < %.3f)', tolS)); legend('Location', 'northwest');

% ---------------- A3: load (number of users) ----------------
nSK = size(EN_Kbe, 3);
figure('Name', 'A3 Which array: performance vs number of users', 'Position', [30 30 1150 470]);
for j = 1:numel(Ksw_sch)
    sc = Ksw_sch(j);
    if en_cape(sc) == 0 || en_idea(sc) > 0, continue; end
    a = en_arr(sc);  Xb = reshape(EN_Kbe(j, :, :), numel(Ksw), nSK);  Xs = reshape(EN_Ksum(j, :, :), numel(Ksw), nSK);
    subplot(1, 2, 1); hold on;
    m = mean(Xb, 2);  s = std(Xb, 0, 2) / sqrt(nSK);
    h = errorbar(Ksw, max(m, tFl), min(s, 0.9 * max(m, tFl)), s, ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
    subplot(1, 2, 2); hold on;
    h = errorbar(Ksw, mean(Xs, 2), std(Xs, 0, 2) / sqrt(nSK), ['-' tMk{sc}]);
    set(h, 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.5, 'DisplayName', tLegC{a});
end
subplot(1, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
title(sprintf('BER vs load, %d dB, N = %d antennas per AP', SNR_dB(tSi), en_N0(1))); legend('Location', 'southeast');
subplot(1, 2, 2); box on; grid on; xlabel('number of users K'); ylabel('sum SE [bit/s/Hz]');
title(sprintf('Sum SE vs load, %d dB', SNR_dB(tSi))); legend('Location', 'southeast');

% ---------------- A4: the physical reasons ----------------
rs0 = rng;  rng(11);                           % fixed seed: same hotspot geometries for every array
tMu = zeros(4, 1);  nTrMu = 300;
for a = 1:4
    pa = EN_geo{tArrSch(a), 1}.pos(:);  mu = 0;
    for t = 1:nTrMu                            % 16 users in a 25 m hotspot, 40-150 m from the AP, 10 m height difference
        c = (40 + 110 * rand) * exp(1j * 2 * pi * rand);
        U = c + 25 * sqrt(rand(en_K(1), 1)) .* exp(1j * 2 * pi * rand(en_K(1), 1));
        dH = abs(U);  d3 = sqrt(dH .^ 2 + 100);  s = sin(angle(U)) .* dH ./ d3;
        rn = sqrt(d3.' .^ 2 + pa .^ 2 - 2 * pa * (d3.' .* s.'));
        Am = (d3.' ./ rn) .* exp(-1j * 2 * pi * rn / lambda);  Am = Am ./ sqrt(sum(abs(Am) .^ 2, 1));
        G2 = abs(Am' * Am) .^ 2;  mu = mu + mean(G2(~eye(en_K(1))));
    end
    tMu(a) = mu / nTrMu;
end
rng(rs0);
tRes = zeros(4, 2);
for a = 1:4
    v = DOA_D(DOA_pP(a, :) >= 0.9);  if ~isempty(v), tRes(a, 1) = max(v); end
    v = DOA_D(DOA_pC(a, :) >= 0.9);  if ~isempty(v), tRes(a, 2) = max(v); end
end
tGl = zeros(4, 1);  tNfA = zeros(4, 1);
for a = 1:4
    Ga = EN_geo{tArrSch(a), 1};  tGl(a) = min(2, lambda / Ga.unit);
    tNfA(a) = 100 * mean(EN_nf(tArrSch(a), 1, :));
end
figure('Name', 'A4 Why the arrays differ: physical metrics', 'Position', [30 30 1500 760]);
tV = {tAp, tNfA, tGl, tMu, tUd};
tY = {'aperture [m]', 'near-field links [%]', 'grating-lobe period in sin\theta (2 = none)', ...
      'mean |a_i^H a_j|^2 between hotspot users', 'hole-free co-array length'};
tT = {'Aperture: beam width and range resolution', 'Links in the near field (range information)', ...
      'Grating lobes: how often a false angle repeats', 'User-channel correlation (lower = easier to separate)', ...
      'Difference co-array (virtual antennas)'};
for mI = 1:5
    subplot(2, 3, mI); hold on; box on; grid on;
    for a = 1:4, bar(a, tV{mI}(a), 0.6, 'FaceColor', tArrCol{a}); end
    set(gca, 'XTick', 1:4, 'XTickLabel', tLegC); xlim([0.4 4.6]); ylabel(tY{mI}); title(tT{mI});
end
subplot(2, 3, 6); hb = bar(tRes); box on; grid on;
set(gca, 'XTick', 1:4, 'XTickLabel', tLegC); ylabel('users all resolved in \geq 90% of trials');
legend(hb, {'physical MUSIC', 'co-array MUSIC'}, 'Location', 'northwest');
title(sprintf('Resolvable users (DOA test, %d dB)', DOA_snrD));

% ---------------- B1: co-array = more resolvable users (degrees of freedom) ----------------
figure('Name', 'B1 Co-array: resolvable users (degrees of freedom)', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on; box on; grid on;
    for a = 1:4
        if f == 1, yP = 100 * DOA_fP(a, :);  yC = 100 * DOA_fC(a, :);
        else, yP = 100 * DOA_pP(a, :);  yC = 100 * DOA_pC(a, :); end
        plot(DOA_D, yP, ['-' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s, physical (max %d)', tLegC{a}, DOA_lim(a, 1)));
        plot(DOA_D, yC, ['--' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', 'none', 'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s, co-array (max %d)', tLegC{a}, DOA_lim(a, 2)));
    end
    xlabel('number of users (sources) D'); ylim([0 105]);
    if f == 1
        ylabel('users found [%]'); title(sprintf('Users found, same 8 antennas, %d dB', DOA_snrD));
        legend('Location', 'southwest');
    else
        ylabel('ALL users resolved [%]'); title('Probability that all users are resolved');
    end
end

% ---------------- B2: co-array gains vs SNR (paired, same antennas) ----------------
figure('Name', 'B2 Co-array: paired gains vs SNR (same antennas)', 'Position', [30 30 1150 820]);
tB2t = {'position error change [%]  (<0 = co-array better)', 'true-lobe hit rate change [points]  (>0 = better)', ...
        'NMSE change [%]  (<0 = better)', 'BER change [%] (Cross-AP List)  (no loss = around 0)'};
for f = 1:4
    subplot(2, 2, f); hold on; box on; grid on;
    plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for cp = 1:2
        PP = tCoP;  if cp == 2, PP = tCoL; end
        for pr = 1:2
            a0 = PP(pr, 1);  b0 = PP(pr, 2);  aa = en_arr(a0);
            if f == 1 || f == 3                    % ratio of positive quantities: geometric mean
                if f == 1, X = tRow(EN_loc, a0);  Y = tRow(EN_loc, b0); else, X = tRow(EN_nmse, a0);  Y = tRow(EN_nmse, b0); end
                mu = nan(tNsnr, 1);  lo = mu;  hi = mu;
                for si = 1:tNsnr
                    ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
                    if nnz(ok) < 2, continue; end
                    rr = log(Y(si, ok) ./ X(si, ok));  mm = mean(rr);  ss = std(rr) / sqrt(nnz(ok));
                    mu(si) = 100 * (exp(mm) - 1);  lo(si) = 100 * (exp(mm - ss) - 1);  hi(si) = 100 * (exp(mm + ss) - 1);
                end
                h = errorbar(SNR_dB, mu, mu - lo, hi - mu);
            elseif f == 2
                Dd = 100 * (tRow(EN_hit, b0) - tRow(EN_hit, a0));  [mm, ss] = tNan(Dd);
                h = errorbar(SNR_dB, mm, ss);
            else
                X = tRowB(tDet, a0);  Y = tRowB(tDet, b0);  ref = max(mean(X, 2), eps);
                h = errorbar(SNR_dB, 100 * mean(Y - X, 2) ./ ref, 100 * std(Y - X, 0, 2) / sqrt(tSet) ./ ref);
            end
            ls = '-';  mf = tArrCol{aa};  nmc = ', CAPE';
            if cp == 2, ls = ':';  mf = 'none';  nmc = ', no CAPE'; end
            set(h, 'LineStyle', ls, 'Marker', tMk{a0}, 'Color', tArrCol{aa}, 'MarkerFaceColor', mf, ...
                'LineWidth', 1.4, 'DisplayName', [tLegC{aa} nmc]);
        end
    end
    xlabel(tXl); ylabel(tB2t{f}); legend('Location', 'best');
end

% ---------------- B3: where the co-array helps (summary) ----------------
tB3a = zeros(4, 2);  tB3b = zeros(2, 2);
for pr = 1:2
    a0 = tCoP(pr, 1);  b0 = tCoP(pr, 2);  aa = en_arr(a0);
    gm = @(X, Y) 100 * (exp(mean(log(Y(isfinite(X) & isfinite(Y) & X > 0 & Y > 0) ./ X(isfinite(X) & isfinite(Y) & X > 0 & Y > 0)))) - 1);
    X = tRowB(tDet, a0);  Y = tRowB(tDet, b0);  tB3a(1, pr) = 100 * (mean(Y(tS12, :)) / mean(X(tS12, :)) - 1);
    X = tRowS(tDet, a0);  Y = tRowS(tDet, b0);  tB3a(2, pr) = 100 * (mean(Y(tS12, :)) / mean(X(tS12, :)) - 1);
    X = tRow(EN_nmse, a0);  Y = tRow(EN_nmse, b0);  tB3a(3, pr) = gm(X(tS12, :), Y(tS12, :));
    X = tRow(EN_loc, a0);   Y = tRow(EN_loc, b0);   tB3a(4, pr) = gm(X(tS12, :), Y(tS12, :));
    Dd = 100 * (tRow(EN_hit, b0) - tRow(EN_hit, a0));  v = Dd(tS12, isfinite(Dd(tS12, :)));  tB3b(1, pr) = mean(v);
    tB3b(2, pr) = tRes(aa, 2) - tRes(aa, 1);
end
figure('Name', sprintf('B3 Co-array: where it helps (%d dB, with CAPE)', SNR_dB(tS12)), 'Position', [30 30 1150 470]);
subplot(1, 2, 1); hb = bar(tB3a); box on; grid on; hold on;
for pr = 1:2, set(hb(pr), 'FaceColor', tArrCol{en_arr(tCoP(pr, 1))}); end
set(gca, 'XTick', 1:4, 'XTickLabel', {'BER', 'sum SE', 'NMSE', 'position error'});
ylabel('change with co-array [%]'); legend(hb, {'Nested', 'Coprime'}, 'Location', 'southwest');
title('Communication metrics: small changes (BER, SE: no loss)');
subplot(1, 2, 2); hb = bar(tB3b); box on; grid on;
for pr = 1:2, set(hb(pr), 'FaceColor', tArrCol{en_arr(tCoP(pr, 1))}); end
set(gca, 'XTick', 1:2, 'XTickLabel', {'true-lobe hit rate [points]', 'resolvable users [extra users]'});
ylabel('gain with co-array'); legend(hb, {'Nested', 'Coprime'}, 'Location', 'northwest');
title('Angle-domain metrics: where the co-array helps');

% ---------------- C1: CAPE resolves the angle ambiguity and localises ----------------
figure('Name', 'C1 CAPE: angle ambiguity and localisation (with vs without CAPE)', 'Position', [30 30 1500 470]);
for a = tSp
    for cp = 1:2
        sc = tArrSch(a);  ls = '-';  mf = tArrCol{a};  nmc = ' + CAPE';
        if cp == 2, sc = tNoC(sc);  ls = ':';  mf = 'none';  nmc = ', no CAPE'; end
        [m, s] = tNan(tRow(EN_hit, sc));
        subplot(1, 3, 1); hold on;
        h = errorbar(SNR_dB, 100 * m, 100 * s);
        set(h, 'LineStyle', ls, 'Marker', tMk{sc}, 'Color', tArrCol{a}, 'MarkerFaceColor', mf, 'LineWidth', 1.4, 'DisplayName', [tLegC{a} nmc]);
        subplot(1, 3, 2); hold on;
        plot(SNR_dB, tMedN(tRow(EN_loc, sc)), 'LineStyle', ls, 'Marker', tMk{sc}, 'Color', tArrCol{a}, ...
             'MarkerFaceColor', mf, 'LineWidth', 1.4, 'DisplayName', [tLegC{a} nmc]);
        [m, s] = tNan(tRow(EN_loc1, sc));
        subplot(1, 3, 3); hold on;
        h = errorbar(SNR_dB, 100 * m, 100 * s);
        set(h, 'LineStyle', ls, 'Marker', tMk{sc}, 'Color', tArrCol{a}, 'MarkerFaceColor', mf, 'LineWidth', 1.4, 'DisplayName', [tLegC{a} nmc]);
    end
end
subplot(1, 3, 1); box on; grid on; xlabel(tXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
title('Angle ambiguity: true lobe selected'); legend('Location', 'northwest');
subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(tXl); ylabel('median position error [m]');
title('Position error: CAPE vs local triangulation'); legend('Location', 'northeast');
subplot(1, 3, 3); box on; grid on; xlabel(tXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
title('Users located within 1 m'); legend('Location', 'northwest');

% ---------------- C2: CAPE gain for List SIC vs the proposed Cross-AP List ----------------
figure('Name', 'C2 CAPE gain: List SIC vs the proposed Cross-AP List', 'Position', [30 30 1500 470]);
tDc = {'--', '-'};  tDl = {'List SIC', 'Cross-AP List (proposed)'};
for ia = 1:numel(tSp)
    a = tSp(ia);  sc = tArrSch(a);  sn = tNoC(sc);
    subplot(1, 3, ia); hold on; box on; grid on;
    plot(SNR_dB, ones(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for id = 1:2
        d = id + 2;  X = tRowB(d, sn);  Y = tRowB(d, sc);       % X: no CAPE, Y: CAPE
        mx = mean(X, 2);  my = mean(Y, 2);  ok = mx > 0 & my > 0;
        R = nan(tNsnr, 1);  se = R;
        cxy = arrayfun(@(i) sum((X(i, :) - mx(i)) .* (Y(i, :) - my(i))) / (tSet - 1), (1:tNsnr).');
        R(ok) = mx(ok) ./ my(ok);
        se(ok) = R(ok) .* sqrt(var(X(ok, :), 0, 2) ./ (tSet * mx(ok) .^ 2) + var(Y(ok, :), 0, 2) ./ (tSet * my(ok) .^ 2) ...
                 - 2 * cxy(ok) ./ (tSet * mx(ok) .* my(ok)));
        h = errorbar(SNR_dB, R, real(se));
        mf = tArrCol{a};  if id == 1, mf = 'none'; end
        set(h, 'LineStyle', tDc{id}, 'Marker', tMk{sc}, 'Color', tArrCol{a}, 'MarkerFaceColor', mf, ...
            'LineWidth', 1.4, 'DisplayName', tDl{id});
    end
    xlabel(tXl); ylabel('BER(no CAPE) / BER(CAPE)   (>1: CAPE helps)');
    title(sprintf('%s: CAPE gain per detector', tLegC{a})); legend('Location', 'best');
end

% ---------------- C3: building the proposed receiver ----------------
figure('Name', 'C3 Proposed receiver: List SIC -> Cross-AP List, without -> with CAPE', 'Position', [30 30 1500 470]);
tC3s = {':', '--', '-.', '-'};
tC3n = {'List SIC, no CAPE', 'List SIC + CAPE', 'Cross-AP List, no CAPE', 'Cross-AP List + CAPE (proposed)'};
for ia = 1:numel(tSp)
    a = tSp(ia);  sc = tArrSch(a);  sn = tNoC(sc);
    subplot(1, 3, ia); hold on; box on; grid on;
    cfg = [3 sn; 3 sc; 4 sn; 4 sc];
    for q = 1:4
        X = tRowB(cfg(q, 1), cfg(q, 2));  m = max(mean(X, 2), tFl);  s = std(X, 0, 2) / sqrt(tSet);
        h = errorbar(SNR_dB, m, min(s, 0.9 * m), s);
        mf = tArrCol{a};  if mod(q, 2) == 1, mf = 'none'; end
        set(h, 'LineStyle', tC3s{q}, 'Marker', tMk{sc}, 'Color', tArrCol{a}, 'MarkerFaceColor', mf, ...
            'LineWidth', 1.4, 'DisplayName', tC3n{q});
    end
    set(gca, 'YScale', 'log'); xlabel(tXl); ylabel('BER');
    title(tLegC{a}); legend('Location', 'southwest');
end

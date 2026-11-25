% ======================================================================
%  D1-D8: THE DIFFERENCE CO-ARRAY, ONE PARAMETER PER FIGURE, ONE PANEL PER
%  ARRAY. In every panel the SAME antennas are processed two ways:
%  physical array (solid, filled marker) and difference co-array (dashed,
%  open marker). D1-D4: array-processing test (DOA). D5-D8: inside the
%  cell-free system, with CAPE (thick) and without CAPE (thin).
% ======================================================================
dArr  = 1:4;                                   % ULA lambda/2, ULA 4 lambda, nested, coprime
dSys  = 2:4;                                   % arrays with near-field links (system figures)
dLimV = @(x, yl, c) plot([x x], yl, ':', 'Color', c, 'LineWidth', 1, 'HandleVisibility', 'off');
dFind = @(a, co, cp) find(en_arr == a & (en_idea > 0) == co & en_cape == cp, 1);
dVirt = zeros(4, 1);
for a = dArr
    Ga = EN_geo{tArrSch(a), 1};  pa = Ga.pos(:);
    lg = unique(round((pa - pa.') / Ga.unit));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
    dVirt(a) = max(diff(bk));                   % hole-free co-array length (virtual antennas)
end

% ---------------- D1 / D2: degrees of freedom ----------------
dNm = {'D1 Co-array vs physical: users found vs number of users', ...
       'D2 Co-array vs physical: ALL users resolved vs number of users'};
for f = 1:2
    figure('Name', dNm{f}, 'Position', [30 30 1500 430]);
    for a = dArr
        subplot(1, 4, a); hold on; box on; grid on;
        if f == 1, yP = 100 * DOA_fP(a, :);  yC = 100 * DOA_fC(a, :);
        else, yP = 100 * DOA_pP(a, :);  yC = 100 * DOA_pC(a, :); end
        plot(DOA_D, yP, ['-' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.6, ...
             'DisplayName', sprintf('physical (%d antennas)', numel(EN_geo{tArrSch(a), 1}.pos)));
        plot(DOA_D, yC, ['--' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', 'none', 'LineWidth', 1.6, ...
             'DisplayName', sprintf('co-array (%d virtual)', dVirt(a)));
        dLimV(DOA_lim(a, 1) + 0.5, [0 105], [0.3 0.3 0.3]);  dLimV(DOA_lim(a, 2) + 0.5, [0 105], tArrCol{a});
        ylim([0 105]); xlabel('number of users D');
        if f == 1, ylabel('users found [%]'); else, ylabel('all users resolved [%]'); end
        title(sprintf('%s: max %d physical, %d co-array', tLegC{a}, DOA_lim(a, 1), DOA_lim(a, 2)));
        legend('Location', 'southwest');
    end
end

% ---------------- D3: angular resolution of two close users ----------------
if exist('DOA_resP', 'var') && ~isempty(DOA_resP)
    figure('Name', 'D3 Co-array vs physical: resolving two close users', 'Position', [30 30 1500 430]);
    for a = dArr
        subplot(1, 4, a); hold on; box on; grid on;
        plot(DOA_sepB, 100 * DOA_resP(a, :), ['-' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, ...
             'LineWidth', 1.6, 'DisplayName', 'physical');
        plot(DOA_sepB, 100 * DOA_resC(a, :), ['--' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', 'none', ...
             'LineWidth', 1.6, 'DisplayName', 'co-array');
        ylim([0 105]); xlabel('separation of the two users [physical beamwidths]'); ylabel('both resolved [%]');
        title(sprintf('%s, %d dB', tLegC{a}, DOA_snrD)); legend('Location', 'southeast');
    end
end

% ---------------- D4: angle accuracy (RMSE) ----------------
figure('Name', sprintf('D4 Co-array vs physical: angle RMSE (%d users)', DOA_D0), 'Position', [30 30 1500 430]);
for a = dArr
    subplot(1, 4, a); hold on; box on; grid on;
    semilogy(DOA_snr, DOA_rP(a, :), ['-' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, ...
             'LineWidth', 1.6, 'DisplayName', 'physical');
    semilogy(DOA_snr, DOA_rC(a, :), ['--' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', 'none', ...
             'LineWidth', 1.6, 'DisplayName', 'co-array');
    set(gca, 'YScale', 'log'); xlabel('SNR per antenna per user [dB]'); ylabel('angle RMSE [sin\theta]');
    title(sprintf('%s, %d users', tLegC{a}, DOA_D0)); legend('Location', 'southwest');
end

% ---------------- D5-D8: inside the cell-free system ----------------
dQ  = {'hit', 'loc', 'nmse', 'ber'};
dNs = {'D5 Co-array vs physical: angle (true lobe) in the system', ...
       'D6 Co-array vs physical: position error in the system', ...
       'D7 Co-array vs physical: channel NMSE in the system', ...
       'D8 Co-array vs physical: BER in the system (no-loss check)'};
dYl = {'true-lobe hit rate [%]', 'median position error [m]', 'channel NMSE', 'BER (Cross-AP List)'};
dSt = {'-', '--', ':', '-.'};  dLw = [1.8 1.8 1.1 1.1];
dLg = {'physical, CAPE', 'co-array, CAPE', 'physical, no CAPE', 'co-array, no CAPE'};
for f = 1:4
    figure('Name', dNs{f}, 'Position', [30 30 1500 450]);
    for ia = 1:numel(dSys)
        a = dSys(ia);
        subplot(1, numel(dSys), ia); hold on; box on; grid on;
        cfg = [0 1; 1 1; 0 0; 1 0];               % [co-array, CAPE]
        nPl = 0;
        for q = 1:4
            sc = dFind(a, cfg(q, 1), cfg(q, 2));
            if isempty(sc), continue; end
            mf = tArrCol{a};  if cfg(q, 1) == 1, mf = 'none'; end
            switch dQ{f}
                case 'hit',  [m, s] = tNan(tRow(EN_hit, sc));  h = errorbar(SNR_dB, 100 * m, 100 * s);
                case 'loc',  h = plot(SNR_dB, tMedN(tRow(EN_loc, sc)));
                case 'nmse', X = tRow(EN_nmse, sc);  h = errorbar(SNR_dB, mean(X, 2), std(X, 0, 2) / sqrt(tSet));
                case 'ber'
                    X = tRowB(tDet, sc);  m = max(mean(X, 2), tFl);  s = std(X, 0, 2) / sqrt(tSet);
                    h = errorbar(SNR_dB, m, min(s, 0.9 * m), s);
            end
            set(h, 'LineStyle', dSt{q}, 'Marker', tMk{tArrSch(a)}, 'Color', tArrCol{a}, 'MarkerFaceColor', mf, ...
                'LineWidth', dLw(q), 'DisplayName', dLg{q});
            nPl = nPl + 1;
        end
        if f >= 2, set(gca, 'YScale', 'log'); end
        if f == 1, ylim([0 100]); end
        xlabel(tXl); ylabel(dYl{f});
        if isempty(dFind(a, 1, 1))
            title(sprintf('%s (no co-array scheme in this run)', tLegC{a}));
        else
            title(tLegC{a});
        end
        if nPl > 0, legend('Location', 'best'); end
    end
end

% ---------------- console: the co-array, array by array ----------------
fprintf('\n=== [CODE T] DIFFERENCE CO-ARRAY, ARRAY BY ARRAY ===\n');
fprintf('%-16s %9s %9s %16s %16s', 'array', 'antennas', 'virtual', 'max users P/C', 'resolved P/C');
if exist('DOA_resP', 'var') && ~isempty(DOA_resP), fprintf(' %20s', '50%-resolution P/C'); end
fprintf('\n');
for a = dArr
    fprintf('%-16s %9d %9d %10d / %-4d %10d / %-4d', strrep(tLegC{a}, '\lambda', 'l'), numel(EN_geo{tArrSch(a), 1}.pos), ...
            dVirt(a), DOA_lim(a, 1), DOA_lim(a, 2), tRes(a, 1), tRes(a, 2));
    if exist('DOA_resP', 'var') && ~isempty(DOA_resP)
        r1 = DOA_sepB(find(DOA_resP(a, :) >= 0.5, 1));  r2 = DOA_sepB(find(DOA_resC(a, :) >= 0.5, 1));
        if isempty(r1), r1 = NaN; end
        if isempty(r2), r2 = NaN; end
        fprintf('   %5.2f / %-5.2f bw', r1, r2);
    end
    fprintf('\n');
end
fprintf('  in the system at %d dB (physical -> co-array, with CAPE): hit [%%], position [m], NMSE, BER\n', SNR_dB(tS12));
for a = dSys
    s0 = dFind(a, 0, 1);  s1 = dFind(a, 1, 1);
    if isempty(s1), continue; end
    [h0, ~] = tNan(tRow(EN_hit, s0));  [h1, ~] = tNan(tRow(EN_hit, s1));
    l0 = tMedN(tRow(EN_loc, s0));  l1 = tMedN(tRow(EN_loc, s1));
    n0 = mean(tRow(EN_nmse, s0), 2);  n1 = mean(tRow(EN_nmse, s1), 2);
    b0 = mean(tRowB(tDet, s0), 2);  b1 = mean(tRowB(tDet, s1), 2);
    fprintf('    %-12s hit %3.0f -> %3.0f   pos %6.2f -> %6.2f   NMSE %.2e -> %.2e   BER %.2e -> %.2e\n', ...
            strrep(tLegC{a}, '\lambda', 'l'), 100 * h0(tS12), 100 * h1(tS12), l0(tS12), l1(tS12), ...
            n0(tS12), n1(tS12), b0(tS12), b1(tS12));
end

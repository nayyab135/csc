%% PLOT_CODET_FIGURES  Redraw Code T figures T1-T17 without rerunning the simulation.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_hit', 'var') || ~exist('DOA_fC', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
%% ===================== [CODE T] FIGURES T1-T17 =====================
%  T1-T10 : the six schemes WITH CAPE (as before). Colour = physical array.
%           Solid line + filled marker = processing on the PHYSICAL array;
%           dashed line + open marker = the SAME antennas with difference
%           co-array processing added.
%  T11-T14: CAPE versus NO CAPE (local per-AP polar-domain OMP), same drops,
%           channels, pilots and symbols. Dotted = no CAPE, physical processing;
%           dash-dot = no CAPE, co-array processing.
%  T15    : angle acquisition of the channel estimator (true-lobe hit rate).
%  T16    : paired effect of co-array processing on NMSE, localisation, angle.
%  T17    : classic DOA test of the four arrays: physical MUSIC vs co-array
%           MUSIC, sources found vs number of sources, RMSE vs SNR.
%  SNR axis = received SNR per antenna at each user's strongest AP.
tNs   = numel(en_names);  tNc = numel(en_N0);  tNsnr = numel(SNR_dB);
tSet  = size(EN_be, 5);
[~, tSi] = min(abs(SNR_dB - 10));             % 10 dB reference point
tDet  = 4;                                     % Cross-AP List detector
tXl   = 'SNR [dB] (received, strongest AP)';
tC    = find(en_cape == 1);                    % schemes with CAPE
tArrNm  = {'ULA (\lambda/2)', en_leg{2}, en_leg{3}, en_leg{5}};
tArrSch = [1 2 3 5];                           % one CAPE scheme per physical array
tArrCol = {[0.35 0.35 0.35], [0.00 0.45 0.00], [0.00 0.10 0.55], [0.50 0.10 0.60]};
tCol  = tArrCol(en_arr);
tLs   = repmat({'-'}, 1, tNs);   tLs(en_idea > 0) = {'--'};
tMk   = repmat({'d', 'o', 's', 's', '^', '^'}, 1, ceil(tNs / 6));
tMfc  = tCol;                    tMfc(en_idea > 0) = {'none'};
tLs2  = tLs;                                   % T11-T16 styles: no CAPE dotted / dash-dot
tLs2(en_cape == 0 & en_idea == 0) = {':'};
tLs2(en_cape == 0 & en_idea > 0)  = {'-.'};
tMfc2 = tMfc;                    tMfc2(en_cape == 0) = {'none'};
tPair = [3 4; 5 6];                            % [physical, co-array] pairs, CAPE
tPairNm = {'Nested', 'Coprime'};
tNoC  = zeros(1, tNs);                         % matching no-CAPE scheme of every CAPE scheme
for sc = tC
    jn = find(en_cape == 0 & strcmp(en_geom, en_geom{sc}) & en_idea == en_idea(sc), 1);
    if ~isempty(jn), tNoC(sc) = jn; end
end
tPairL = [tNoC(3) tNoC(4); tNoC(5) tNoC(6)];   % [physical, co-array] pairs, no CAPE
tW    = 80 + 560 * tNc;                        % figure width for tNc case panels
tFl   = 1 / (max(en_K) * nRealE * nSymE * log2(modOrder) * tSet);   % one bit error in all drops
mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
mLC = nan(tNs, tNc, tNsnr);  mHT = nan(tNs, tNc, tNsnr);  mRM = nan(tNs, tNc, tNsnr);
for sc = 1:tNs
    for ca = 1:tNc
        for si = 1:tNsnr
            v = EN_loc(sc, ca, si, :);   v = v(isfinite(v));
            if ~isempty(v), mLC(sc, ca, si) = median(v); end
            v = EN_hit(sc, ca, si, :);   v = v(isfinite(v));
            if ~isempty(v), mHT(sc, ca, si) = mean(v); end
            v = EN_rmse(sc, ca, si, :);  v = v(isfinite(v));
            if ~isempty(v), mRM(sc, ca, si) = sqrt(mean(v .^ 2)); end
        end
    end
end

% ---- T1: physical arrays and their hole-free difference co-arrays ----
figure('Name', 'T1 Physical arrays and difference co-arrays', 'Position', [30 30 1100 760]);
for a = 1:4
    Ga = EN_geo{tArrSch(a), 1};  pa = Ga.pos(:);
    uu = lambda / 2;  if isfield(Ga, 'unit'), uu = Ga.unit; end
    lg = unique(round((pa - pa.') / uu));
    bk = [0; find(diff(lg) ~= 1); numel(lg)];  [nv, mi] = max(diff(bk));
    cont = lg(bk(mi) + 1:bk(mi + 1));
    ap = max(pa) - min(pa);
    subplot(4, 1, a); hold on; box on; grid on;
    plot(pa, ones(size(pa)), 'o', 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'MarkerSize', 5);
    plot(cont * uu, 0.6 * ones(size(cont)), '.', 'Color', tArrCol{a}, 'MarkerSize', 6);
    ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
    xlabel('position along the array [m]');
    title(sprintf('%s:  N = %d,  aperture = %.2f m,  grid %.1f\\lambda,  hole-free co-array = %d,  Rayleigh = %.0f m', ...
          tArrNm{a}, numel(pa), ap, uu / lambda, nv, 2 * ap^2 / lambda));
end

% ---- T2: far-field beampattern (where the grating lobes are) ----
figure('Name', 'T2 Far-field beampattern (grating lobes)', 'Position', [30 30 1100 650]);
sg = linspace(-1, 1, 4001);
for a = 1:4
    pa = EN_geo{tArrSch(a), 1}.pos(:).';
    B = abs(exp(1j * 2 * pi * (sg.') * pa / lambda) * ones(numel(pa), 1)) / numel(pa);
    subplot(2, 2, a); plot(sg, 20 * log10(max(B, 1e-4)), 'Color', tArrCol{a}, 'LineWidth', 1.2);
    grid on; box on; ylim([-40 0]); xlabel('sin(\theta)'); ylabel('beampattern [dB]'); title(tArrNm{a});
end

% ---- T3-T5: NMSE, BER and sum SE vs SNR (the six schemes with CAPE) ----
tMet = {mNM, mBE, mSU};
tNm  = {'T3 Channel NMSE vs SNR', 'T4 BER vs SNR', 'T5 Sum SE vs SNR'};
tYl  = {'channel NMSE', 'BER (Cross-AP List)', 'sum SE [bit/s/Hz]'};
tLoc = {'southwest', 'southwest', 'northwest'};
for f = 1:3
    figure('Name', tNm{f}, 'Position', [30 30 tW 470]);
    for ca = 1:tNc
        subplot(1, tNc, ca); hold on; box on; grid on;
        for sc = tC
            if f == 1
                y = squeeze(tMet{f}(sc, ca, :));
            else
                y = squeeze(tMet{f}(tDet, sc, ca, :));
            end
            if f == 2, y = max(y, tFl); end
            plot(SNR_dB, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
                 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
        end
        if f < 3, set(gca, 'YScale', 'log'); end
        if f == 2, ylim([max(tFl, 1e-7) 1]); end
        xlabel(tXl); ylabel(tYl{f});
        title(sprintf('%s, all with CAPE', en_caseName{ca}));
        legend('Location', tLoc{f});
    end
end

% ---- T6: BER and sum SE vs number of users K at 10 dB (schemes with CAPE) ----
mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);
figure('Name', 'T6 BER and sum SE vs number of users', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on; box on; grid on;
    for j = 1:numel(Ksw_sch)
        sc = Ksw_sch(j);
        if en_cape(sc) == 0, continue; end
        if f == 1, y = max(mKbe(j, :), tFl); else, y = mKsum(j, :); end
        plot(Ksw, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
             'LineWidth', 1.5, 'DisplayName', en_leg{sc});
    end
    xlabel('number of users K');
    if f == 1
        set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)');
        title(sprintf('BER vs K at %d dB, N = %d per AP, with CAPE', SNR_dB(tSi), en_N0(1)));
    else
        ylabel('sum SE [bit/s/Hz]');
        title(sprintf('Sum SE vs K at %d dB, N = %d per AP, with CAPE', SNR_dB(tSi), en_N0(1)));
    end
    legend('Location', 'northwest');
end

% ---- T7: does co-array processing help? paired difference, same antennas ----
figure('Name', 'T7 Effect of co-array processing on BER and SE (paired)', 'Position', [30 30 tW 760]);
tYl7 = {'BER change with co-array [%]  (<0 = better)', 'sum-SE change with co-array [%]  (>0 = better)'};
for row = 1:2
    for ca = 1:tNc
        subplot(2, tNc, (row - 1) * tNc + ca); hold on; box on; grid on;
        plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        ttl = '';
        for pr = 1:2
            a = tPair(pr, 1);  b = tPair(pr, 2);
            if row == 1
                X = EN_be(tDet, a, ca, :, :);   Y = EN_be(tDet, b, ca, :, :);
            else
                X = EN_sum(tDet, a, ca, :, :);  Y = EN_sum(tDet, b, ca, :, :);
            end
            X = reshape(X, tNsnr, tSet);  Y = reshape(Y, tNsnr, tSet);  Dd = Y - X;
            ref = max(abs(mean(X, 2)), eps);
            mu  = 100 * mean(Dd, 2) ./ ref;
            se  = 100 * std(Dd, 0, 2) / sqrt(tSet) ./ ref;
            hEb = errorbar(SNR_dB, mu, se, ['-' tMk{a}]);
            set(hEb, 'Color', tCol{a}, 'MarkerFaceColor', tCol{a}, 'LineWidth', 1.4, 'DisplayName', tPairNm{pr});
            ttl = [ttl sprintf('  %s %d%%', tPairNm{pr}, round(100 * mean(Dd(:) == 0)))]; %#ok<AGROW>
        end
        xlabel(tXl); ylabel(tYl7{row});
        title({[en_caseName{ca} ', with CAPE'], ['identical:' ttl]});
        legend('Location', 'best');
    end
end

% ---- T8: localisation error vs SNR (CAPE only: local estimation gives no position) ----
figure('Name', 'T8 Localisation error vs SNR', 'Position', [30 30 tW 470]);
for ca = 1:tNc
    subplot(1, tNc, ca); hold on; box on; grid on;
    for sc = tC(2:end)
        plot(SNR_dB, squeeze(mLC(sc, ca, :)), [tLs{sc} tMk{sc}], 'Color', tCol{sc}, ...
             'MarkerFaceColor', tMfc{sc}, 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
    end
    set(gca, 'YScale', 'log'); xlabel(tXl); ylabel('median position error [m]');
    title([en_caseName{ca} ' (dense ULA: no near-field links, no position)']); legend('Location', 'northeast');
end

% ---- T9: array hardware metrics ----
tAp = zeros(4, 1);  tUd = zeros(4, 1);  tW1 = zeros(4, 1);  tCl = zeros(4, 1);
for a = 1:4
    Ga = EN_geo{tArrSch(a), 1};  pa = Ga.pos(:);
    uu = lambda / 2;  if isfield(Ga, 'unit'), uu = Ga.unit; end
    Dm = round((pa - pa.') / (lambda / 2));            % physical distance in half-wavelengths
    lg = unique(round((pa - pa.') / uu));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
    tAp(a) = max(pa) - min(pa);  tUd(a) = max(diff(bk));  tW1(a) = sum(Dm(:) == 1);
    U = zeros(numel(pa));  nz = (Dm ~= 0) & (abs(Dm) <= 100);
    U(nz) = 1 ./ abs(Dm(nz));  U(Dm == 0) = 1;
    tCl(a) = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
end
figure('Name', 'T9 Array hardware metrics', 'Position', [30 30 1000 720]);
tVal = {tAp, tUd, tW1, tCl};
tYl9 = {'aperture [m]', 'hole-free co-array length', 'w(1): pairs \lambda/2 apart', 'coupling leakage'};
tTt9 = {'Aperture (resolution)', 'Hole-free difference co-array', 'Closest pairs', 'Mutual coupling'};
for m = 1:4
    subplot(2, 2, m); hold on; box on; grid on;
    for a = 1:4
        bar(a, tVal{m}(a), 0.6, 'FaceColor', tArrCol{a});
    end
    set(gca, 'XTick', 1:4, 'XTickLabel', {'ULA \lambda/2', 'ULA wide', 'Nested', 'Coprime'});  xlim([0.4 4.6]);
    ylabel(tYl9{m});  title(tTt9{m});
end

% ---- T10: the four detectors for each physical array (with CAPE) ----
tDn = {'Linear MMSE', 'Hard SIC', 'List SIC', 'Cross-AP List'};
tDs = {':', '-.', '--', '-'};
figure('Name', 'T10 BER of the four detectors', 'Position', [30 30 1100 760]);
for a = 1:4
    sc = tArrSch(a);
    subplot(2, 2, a); hold on; box on; grid on;
    for d = 1:4
        plot(SNR_dB, max(squeeze(mBE(d, sc, 1, :)), tFl), [tDs{d} tMk{sc}], 'Color', tArrCol{a}, ...
             'MarkerFaceColor', tArrCol{a}, 'LineWidth', 1.3, 'DisplayName', tDn{d});
    end
    set(gca, 'YScale', 'log'); ylim([max(tFl, 1e-7) 1]); xlabel(tXl); ylabel('BER');
    title(sprintf('%s, %s, with CAPE', tArrNm{a}, en_caseName{1})); legend('Location', 'southwest');
end

% ---- T11: CAPE benefit on BER, every physical array (with vs without CAPE) ----
figure('Name', 'T11 CAPE benefit: BER with and without CAPE', 'Position', [30 30 1150 800]);
for a = 1:4
    subplot(2, 2, a); hold on; box on; grid on;
    sA = find(en_arr == a);
    for sc = sA
        plot(SNR_dB, max(squeeze(mBE(tDet, sc, 1, :)), tFl), [tLs2{sc} tMk{sc}], 'Color', tCol{sc}, ...
             'MarkerFaceColor', tMfc2{sc}, 'LineWidth', 1.4, 'DisplayName', en_leg{sc});
    end
    sc = tArrSch(a);  g = mBE(tDet, tNoC(sc), 1, tSi) / max(mBE(tDet, sc, 1, tSi), tFl);
    set(gca, 'YScale', 'log'); ylim([max(tFl, 1e-7) 1]); xlabel(tXl); ylabel('BER (Cross-AP List)');
    title(sprintf('%s: BER(no CAPE) / BER(CAPE) = %.2f at %d dB', tArrNm{a}, g, SNR_dB(tSi)));
    legend('Location', 'southwest');
end

% ---- T12: CAPE benefit on NMSE and sum SE ----
figure('Name', 'T12 CAPE benefit: NMSE and sum SE with and without CAPE', 'Position', [30 30 1500 760]);
for a = 1:4
    sA = find(en_arr == a);
    for row = 1:2
        subplot(2, 4, (row - 1) * 4 + a); hold on; box on; grid on;
        for sc = sA
            if row == 1, y = squeeze(mNM(sc, 1, :)); else, y = squeeze(mSU(tDet, sc, 1, :)); end
            plot(SNR_dB, y, [tLs2{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc2{sc}, ...
                 'LineWidth', 1.4, 'DisplayName', en_leg{sc});
        end
        xlabel(tXl);
        if row == 1
            set(gca, 'YScale', 'log'); ylabel('channel NMSE'); title(tArrNm{a});
            legend('Location', 'southwest');
        else
            ylabel('sum SE [bit/s/Hz] (Cross-AP List)'); title(tArrNm{a});
            legend('Location', 'northwest');
        end
    end
end

% ---- T13: CAPE gain per detector at 10 dB, and the gain of the proposed Cross-AP List ----
tGB = nan(4, 6);  tGS = nan(1, 6);  tXG = nan(2, 6);
for j = 1:6
    sc = tC(j);  sn = tNoC(sc);
    if sn == 0, continue; end
    for d = 1:4
        tGB(d, j) = max(mBE(d, sn, 1, tSi), tFl) / max(mBE(d, sc, 1, tSi), tFl);
    end
    tGS(j) = 100 * (mSU(tDet, sc, 1, tSi) / mSU(tDet, sn, 1, tSi) - 1);
    tXG(1, j) = max(mBE(3, sn, 1, tSi), tFl) / max(mBE(4, sn, 1, tSi), tFl);   % List SIC / Cross-AP, no CAPE
    tXG(2, j) = max(mBE(3, sc, 1, tSi), tFl) / max(mBE(4, sc, 1, tSi), tFl);   % List SIC / Cross-AP, CAPE
end
tXn = {'ULA \lambda/2', 'ULA wide', 'Nested', 'Nest.+co', 'Coprime', 'Copr.+co'};
figure('Name', 'T13 CAPE gain per detector and Cross-AP List gain (10 dB)', 'Position', [30 30 1500 470]);
subplot(1, 3, 1); hb = bar(tGB.'); grid on; box on; hold on;
plot([0.4 6.6], [1 1], 'k-', 'HandleVisibility', 'off');
set(gca, 'YScale', 'log', 'XTick', 1:6, 'XTickLabel', tXn);
ylabel('BER(no CAPE) / BER(CAPE)   (>1: CAPE helps)'); legend(hb, tDn, 'Location', 'northwest');
title(sprintf('CAPE gain on BER of each detector, %d dB', SNR_dB(tSi)));
subplot(1, 3, 2); bar(tGS); grid on; box on;
set(gca, 'XTick', 1:6, 'XTickLabel', tXn); ylabel('sum-SE gain of CAPE [%] (Cross-AP List)');
title(sprintf('CAPE gain on sum SE, %d dB', SNR_dB(tSi)));
subplot(1, 3, 3); hb = bar(tXG.'); grid on; box on; hold on;
plot([0.4 6.6], [1 1], 'k-', 'HandleVisibility', 'off');
set(gca, 'YScale', 'log', 'XTick', 1:6, 'XTickLabel', tXn);
ylabel('BER(List SIC) / BER(Cross-AP List)   (>1: Cross-AP helps)');
legend(hb, {'without CAPE', 'with CAPE'}, 'Location', 'northwest');
title(sprintf('Gain of the proposed Cross-AP List, %d dB', SNR_dB(tSi)));

% ---- T14: CAPE benefit vs number of users K at 10 dB ----
figure('Name', 'T14 CAPE benefit vs number of users', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on; box on; grid on;
    for j = 1:numel(Ksw_sch)
        sc = Ksw_sch(j);
        if en_idea(sc) > 0 || en_arr(sc) == 1, continue; end   % physical processing, sparse arrays
        if f == 1, y = max(mKbe(j, :), tFl); else, y = mKsum(j, :); end
        plot(Ksw, y, [tLs2{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc2{sc}, ...
             'LineWidth', 1.5, 'DisplayName', en_leg{sc});
    end
    xlabel('number of users K');
    if f == 1
        set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)');
    else
        ylabel('sum SE [bit/s/Hz] (Cross-AP List)');
    end
    title(sprintf('With CAPE (solid) and without CAPE (dotted), %d dB', SNR_dB(tSi)));
    legend('Location', 'northwest');
end

% ---- T15: angle acquisition of the channel estimator (near-field links) ----
figure('Name', 'T15 Angle acquisition of the channel estimator', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on; box on; grid on;
    for sc = find(en_arr > 1)
        if f == 1, y = 100 * squeeze(mHT(sc, 1, :)); else, y = squeeze(mRM(sc, 1, :)); end
        plot(SNR_dB, y, [tLs2{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc2{sc}, ...
             'LineWidth', 1.4, 'DisplayName', en_leg{sc});
    end
    xlabel(tXl);
    if f == 1
        ylabel('true-lobe hit rate [%]'); ylim([0 100]);
        title(sprintf('Estimated angle within %.4f of the true sin(\\theta)', tolS));
        legend('Location', 'southeast');
    else
        set(gca, 'YScale', 'log'); ylabel('angle RMSE of the hits [sin(\theta)]');
        title('Angle accuracy when the true lobe is found');
    end
end

% ---- T16: paired effect of co-array processing on NMSE, localisation and angle ----
%  Ratio Y/X per drop (co-array / physical, same antennas, same drop):
%  geometric mean over the drops +- one standard error. <0 % = co-array better.
figure('Name', 'T16 Effect of co-array processing on NMSE, localisation and angle (paired)', ...
       'Position', [30 30 1150 800]);
tQ  = {EN_nmse, EN_nmse, EN_loc};
tPp = {tPair, tPairL, tPair};
tT16 = {'NMSE change with co-array, with CAPE', 'NMSE change with co-array, without CAPE', ...
        'localisation-error change with co-array (CAPE)'};
for f = 1:3
    subplot(2, 2, f); hold on; box on; grid on;
    plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for pr = 1:2
        a = tPp{f}(pr, 1);  b = tPp{f}(pr, 2);
        X = reshape(tQ{f}(a, 1, :, :), tNsnr, tSet);  Y = reshape(tQ{f}(b, 1, :, :), tNsnr, tSet);
        mu = nan(tNsnr, 1);  lo = mu;  hi = mu;
        for si = 1:tNsnr
            ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
            if nnz(ok) < 2, continue; end
            r = log(Y(si, ok) ./ X(si, ok));  m = mean(r);  s = std(r) / sqrt(nnz(ok));
            mu(si) = 100 * (exp(m) - 1);  lo(si) = 100 * (exp(m - s) - 1);  hi(si) = 100 * (exp(m + s) - 1);
        end
        hEb = errorbar(SNR_dB, mu, mu - lo, hi - mu, ['-' tMk{a}]);
        set(hEb, 'Color', tCol{a}, 'MarkerFaceColor', tCol{a}, 'LineWidth', 1.4, 'DisplayName', tPairNm{pr});
    end
    xlabel(tXl); ylabel('change [%]  (<0 = co-array better)'); title(tT16{f}); legend('Location', 'best');
end
subplot(2, 2, 4); hold on; box on; grid on;
plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
tPP = [tPair; tPairL];  tPN = {'Nested, CAPE', 'Coprime, CAPE', 'Nested, no CAPE', 'Coprime, no CAPE'};
for pr = 1:4
    a = tPP(pr, 1);  b = tPP(pr, 2);
    X = reshape(EN_hit(a, 1, :, :), tNsnr, tSet);  Y = reshape(EN_hit(b, 1, :, :), tNsnr, tSet);
    Dd = 100 * (Y - X);  mu = nan(tNsnr, 1);  se = mu;
    for si = 1:tNsnr
        v = Dd(si, isfinite(Dd(si, :)));
        if numel(v) >= 2, mu(si) = mean(v);  se(si) = std(v) / sqrt(numel(v)); end
    end
    hEb = errorbar(SNR_dB, mu, se, [tLs2{a}(1:min(end, 2)) tMk{a}]);
    set(hEb, 'Color', tCol{a}, 'MarkerFaceColor', tMfc2{a}, 'LineWidth', 1.4, 'DisplayName', tPN{pr});
end
xlabel(tXl); ylabel('hit-rate change [percentage points]  (>0 = co-array better)');
title('true-lobe angle hit rate change with co-array'); legend('Location', 'best');

% ---- T17: classic DOA test: physical MUSIC vs co-array MUSIC on the same antennas ----
figure('Name', 'T17 Resolvable sources and angle accuracy: physical vs co-array MUSIC', ...
       'Position', [30 30 1500 470]);
for f = 1:3
    subplot(1, 3, f); hold on; box on; grid on;
    for a = 1:4
        if f == 1, yP = 100 * DOA_fP(a, :);  yC = 100 * DOA_fC(a, :);  x = DOA_D;
        elseif f == 2, yP = 100 * DOA_pP(a, :);  yC = 100 * DOA_pC(a, :);  x = DOA_D;
        else, yP = DOA_rP(a, :);  yC = DOA_rC(a, :);  x = DOA_snr;
        end
        plot(x, yP, ['-' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, ...
             'LineWidth', 1.4, 'DisplayName', [tArrNm{a} ', physical MUSIC']);
        plot(x, yC, ['--' tMk{tArrSch(a)}], 'Color', tArrCol{a}, 'MarkerFaceColor', 'none', ...
             'LineWidth', 1.4, 'DisplayName', sprintf('%s, co-array MUSIC (max %d)', tArrNm{a}, DOA_lim(a, 2)));
    end
    if f == 1
        xlabel('number of sources (users) D'); ylabel('sources found [%]'); ylim([0 105]);
        title(sprintf('Sources found, %d dB, %d snapshots', DOA_snrD, DOA_T));
        legend('Location', 'southwest');
    elseif f == 2
        xlabel('number of sources (users) D'); ylabel('all D sources resolved [%]'); ylim([0 105]);
        title(sprintf('Probability that ALL sources are resolved (physical max %d)', DOA_lim(1, 1)));
    else
        xlabel('SNR per antenna per source [dB]'); ylabel('angle RMSE [sin(\theta)]'); set(gca, 'YScale', 'log');
        title(sprintf('Angle RMSE, D = %d sources', DOA_D0));
    end
end

% ---- console summary at 10 dB ----
for ca = 1:tNc
    fprintf('\n=== [CODE T] %s, at %d dB (mean +- std. error over %d drops) ===\n', ...
            en_caseName{ca}, SNR_dB(tSi), tSet);
    fprintf('%-30s %10s %22s %18s %9s %10s\n', 'scheme', 'NMSE', 'BER (Cross-AP List)', 'sum SE', 'hit [%]', 'loc [m]');
    for sc = 1:tNs
        be = reshape(EN_be(tDet, sc, ca, tSi, :), [], 1);
        su = reshape(EN_sum(tDet, sc, ca, tSi, :), [], 1);
        fprintf('%-30s %10.2e %11.3e +- %8.1e %8.2f +- %5.2f %9.1f %10.3f\n', en_names{sc}, mNM(sc, ca, tSi), ...
                mean(be), std(be) / sqrt(tSet), mean(su), std(su) / sqrt(tSet), 100 * mHT(sc, ca, tSi), mLC(sc, ca, tSi));
    end
    fprintf('  CAPE gain at %d dB, BER(no CAPE)/BER(CAPE) for [Linear, Hard SIC, List SIC, Cross-AP List]:\n', SNR_dB(tSi));
    for j = 1:6
        fprintf('    %-26s %s\n', tXn{j}, sprintf('%7.2f ', tGB(:, j)));
    end
    for pr = 1:2
        a = tPair(pr, 1);  b = tPair(pr, 2);
        Dd = EN_be(tDet, b, ca, :, :) - EN_be(tDet, a, ca, :, :);
        fprintf('  %s: co-array vs physical processing (CAPE) identical in %d%% of (SNR, drop) points\n', ...
                tPairNm{pr}, round(100 * mean(Dd(:) == 0)));
    end
    fprintf('  Cross-AP List BER vs SNR (all schemes):\n');
    for sc = 1:tNs
        fprintf('    %-30s', en_names{sc});  fprintf(' %9.2e', squeeze(mBE(tDet, sc, ca, :)));  fprintf('\n');
    end
end
fprintf('\n  DOA test, sources found [%%] at D = %s:\n', mat2str(DOA_D));
for a = 1:4
    fprintf('    %-26s physical %s\n', tArrNm{a}, sprintf('%4.0f', 100 * DOA_fP(a, :)));
    fprintf('    %-26s co-array %s\n', '', sprintf('%4.0f', 100 * DOA_fC(a, :)));
end

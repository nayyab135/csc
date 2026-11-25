%% PLOT_CODET_FIGURES  Redraw the Code T figures without rerunning the simulation.
%  A1-A4: which sparse array performs better; B1-B3: role of the difference
%  co-array; C1-C3: benefit of CAPE; D1-D8: the co-array parameter by
%  parameter, one panel per array; T1, T2, T9: geometry. Set fig_old = true
%  before running this script to also draw the older figures T3-T17.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_loc1', 'var') || ~exist('DOA_fC', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
%% ===================== [CODE T] FIGURES =====================
%  A1-A4: which sparse array performs better, and why
%  B1-B3: role of the difference co-array (same antennas)
%  C1-C3: benefit of CAPE (with vs without, per detector)
%  T1, T2, T9: geometry. Older figures T3-T8, T10-T17: set fig_old = true.
%  ---- details of the older figures ----
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

if ~exist('fig_old', 'var'), fig_old = false; end   % true: also draw the older figures T3-T8, T10-T17
tGB = nan(4, 6);  tGS = nan(1, 6);  tXG = nan(2, 6);   % CAPE gain per detector (used by the summary)
for j = 1:6
    sc = tC(j);  sn = tNoC(sc);
    if sn == 0, continue; end
    for d = 1:4
        tGB(d, j) = max(mBE(d, sn, 1, tSi), tFl) / max(mBE(d, sc, 1, tSi), tFl);
    end
    tGS(j) = 100 * (mSU(tDet, sc, 1, tSi) / mSU(tDet, sn, 1, tSi) - 1);
    tXG(1, j) = max(mBE(3, sn, 1, tSi), tFl) / max(mBE(4, sn, 1, tSi), tFl);
    tXG(2, j) = max(mBE(3, sc, 1, tSi), tFl) / max(mBE(4, sc, 1, tSi), tFl);
end
tXn = {'ULA \lambda/2', 'ULA wide', 'Nested', 'Nest.+co', 'Coprime', 'Copr.+co'};

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

if fig_old
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
    figure('Name', sprintf('T13 CAPE gain per detector and Cross-AP List gain (%d dB)', SNR_dB(tSi)), 'Position', [30 30 1500 470]);
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

end

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

% ---- console summary of figures A, B, C ----
fprintf('\n=== [CODE T] SUMMARY AT %d dB (mean over %d drops) ===\n', SNR_dB(tS12), tSet);
fprintf('%-18s %10s %12s %9s %12s %10s %8s\n', 'array (CAPE)', 'BER', 'BER<=1e-3', 'sum SE', 'pos. err [m]', 'within 1 m', 'hit');
for a = 1:4
    sc = tArrSch(a);  X = tRowB(tDet, sc);  Y = tRowS(tDet, sc);
    md = tMedN(tRow(EN_loc, sc));  [w1, ~] = tNan(tRow(EN_loc1, sc));  [ht, ~] = tNan(tRow(EN_hit, sc));
    fprintf('%-18s %10.2e %11.0f%% %9.1f %12.2f %9.0f%% %7.0f%%\n', strrep(tLegC{a}, '\lambda', 'l'), mean(X(tS12, :)), ...
            100 * mean(X(tS12, :) <= 1e-3), mean(Y(tS12, :)), md(tS12), 100 * w1(tS12), 100 * ht(tS12));
end
fprintf('  users resolvable (all, >= 90%% of trials) physical / co-array: %s\n', mat2str(tRes));
fprintf('  co-array change [%%] (BER, sum SE, NMSE, position error)  nested %s  coprime %s\n', ...
        mat2str(round(tB3a(:, 1).')), mat2str(round(tB3a(:, 2).')));
fprintf('  co-array gain: hit rate [points] nested %.1f coprime %.1f; extra resolvable users nested %d coprime %d\n', ...
        tB3b(1, 1), tB3b(1, 2), tB3b(2, 1), tB3b(2, 2));
fprintf('  CAPE vs no CAPE (hit rate [%%], median position error [m], within 1 m [%%], BER ratio List SIC, BER ratio Cross-AP):\n');
for a = tSp
    sc = tArrSch(a);  sn = tNoC(sc);
    [h1, ~] = tNan(tRow(EN_hit, sc));  [h0, ~] = tNan(tRow(EN_hit, sn));
    l1 = tMedN(tRow(EN_loc, sc));  l0 = tMedN(tRow(EN_loc, sn));
    [w1, ~] = tNan(tRow(EN_loc1, sc));  [w0, ~] = tNan(tRow(EN_loc1, sn));
    B3n = tRowB(3, sn);  B3c = tRowB(3, sc);  B4n = tRowB(4, sn);  B4c = tRowB(4, sc);
    fprintf('    %-10s hit %3.0f -> %3.0f   pos %7.2f -> %7.2f   1 m %3.0f -> %3.0f   List %.2f   Cross-AP %.2f\n', ...
            strrep(tLegC{a}, '\lambda', 'l'), 100 * h0(tS12), 100 * h1(tS12), l0(tS12), l1(tS12), 100 * w0(tS12), 100 * w1(tS12), ...
            mean(B3n(tS12, :)) / max(mean(B3c(tS12, :)), eps), mean(B4n(tS12, :)) / max(mean(B4c(tS12, :)), eps));
end

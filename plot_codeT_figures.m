%% PLOT_CODET_FIGURES  Redraw Code T figures T1-T9 without rerunning the simulation.
%  Uses the variables already in the workspace after Code_T_ArrayComparison.m.
%  If they are not there (new MATLAB session), it loads them from
%  codeT_arraycomp_results.mat. It never clears the workspace.
%  No local functions, so it also works pasted into the Command Window.
if ~exist('EN_be', 'var') || ~exist('EN_geo', 'var') || ~exist('en_leg', 'var')
    load('codeT_arraycomp_results.mat');
end
%% ===================== [CODE T] FIGURES T1-T9 =====================
%  Every scheme uses CAPE. Colour = physical array. Solid line + filled
%  marker = processing on the PHYSICAL array. Dashed line + open marker =
%  the SAME antennas with difference co-array processing added.
tNs   = numel(en_names);  tNc = numel(en_N0);  tNsnr = numel(SNR_dB);
tSet  = size(EN_be, 5);
[~, tSi] = min(abs(SNR_dB - 10));             % 10 dB reference point
tDet  = 4;                                     % Cross-AP List detector
tArrNm  = {'ULA (\lambda/2)', en_leg{2}, 'Nested', 'Coprime'};
tArrSch = [1 2 3 5];                           % one scheme per physical array
tArrCol = {[0.35 0.35 0.35], [0.00 0.45 0.00], [0.00 0.10 0.55], [0.50 0.10 0.60]};
tCol  = tArrCol(en_arr);
tLs   = repmat({'-'}, 1, tNs);   tLs(en_idea > 0) = {'--'};
tMk   = {'d', 'o', 's', 's', '^', '^'};
tMfc  = tCol;                    tMfc(en_idea > 0) = {'none'};
tPair = [3 4; 5 6];                            % [physical, co-array] pairs
tPairNm = {'Nested', 'Coprime'};
mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
mLC = nan(tNs, tNc, tNsnr);
for sc = 1:tNs
    for ca = 1:tNc
        for si = 1:tNsnr
            v = EN_loc(sc, ca, si, :);  v = v(isfinite(v));
            if ~isempty(v), mLC(sc, ca, si) = median(v); end
        end
    end
end

% ---- T1: physical arrays and their hole-free difference co-arrays ----
figure('Name', 'T1 Physical arrays and difference co-arrays', 'Position', [30 30 1100 760]);
for a = 1:4
    pa = EN_geo{tArrSch(a), 1}.pos(:);
    lg = unique(round((pa - pa.') / (lambda / 2)));
    bk = [0; find(diff(lg) ~= 1); numel(lg)];  [nv, mi] = max(diff(bk));
    cont = lg(bk(mi) + 1:bk(mi + 1));
    ap = max(pa) - min(pa);
    subplot(4, 1, a); hold on; box on; grid on;
    plot(pa, ones(size(pa)), 'o', 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'MarkerSize', 5);
    plot(cont * lambda / 2, 0.6 * ones(size(cont)), '.', 'Color', tArrCol{a}, 'MarkerSize', 6);
    ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
    xlabel('position along the array [m]');
    title(sprintf('%s:  N = %d,  aperture = %.2f m,  hole-free co-array = %d lags,  Rayleigh = %.0f m', ...
          tArrNm{a}, numel(pa), ap, nv, 2 * ap^2 / lambda));
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

% ---- T3-T5: NMSE, BER and sum SE vs SNR (all six schemes, both cases) ----
tMet = {mNM, mBE, mSU};
tNm  = {'T3 Channel NMSE vs SNR', 'T4 BER vs SNR', 'T5 Sum SE vs SNR'};
tYl  = {'channel NMSE', 'BER (Cross-AP List)', 'sum SE [bit/s/Hz]'};
tLoc = {'southwest', 'southwest', 'northwest'};
for f = 1:3
    figure('Name', tNm{f}, 'Position', [30 30 1150 470]);
    for ca = 1:tNc
        subplot(1, 2, ca); hold on; box on; grid on;
        for sc = 1:tNs
            if f == 1
                y = squeeze(tMet{f}(sc, ca, :));
            else
                y = squeeze(tMet{f}(tDet, sc, ca, :));
            end
            if f == 2, y = max(y, 1e-6); end
            plot(SNR_dB, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
                 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
        end
        if f < 3, set(gca, 'YScale', 'log'); end
        xlabel('SNR [dB]'); ylabel(tYl{f});
        title(sprintf('%s (K = %d), all with CAPE', en_caseName{ca}, en_K(ca)));
        legend('Location', tLoc{f});
    end
end

% ---- T6: BER and sum SE vs number of users K at 10 dB ----
mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);
figure('Name', 'T6 BER and sum SE vs number of users', 'Position', [30 30 1150 470]);
for f = 1:2
    subplot(1, 2, f); hold on; box on; grid on;
    for j = 1:numel(Ksw_sch)
        sc = Ksw_sch(j);
        if f == 1, y = max(mKbe(j, :), 1e-6); else, y = mKsum(j, :); end
        plot(Ksw, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
             'LineWidth', 1.5, 'DisplayName', en_leg{sc});
    end
    xlabel('number of users K');
    if f == 1
        set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)');
        title(sprintf('BER vs K at %d dB, N = %d per AP', SNR_dB(tSi), en_N0(1)));
    else
        ylabel('sum SE [bit/s/Hz]');
        title(sprintf('Sum SE vs K at %d dB, N = %d per AP', SNR_dB(tSi), en_N0(1)));
    end
    legend('Location', 'northwest');
end

% ---- T7: does co-array processing help? paired difference, same antennas ----
figure('Name', 'T7 Effect of co-array processing (paired)', 'Position', [30 30 1150 760]);
tYl7 = {'BER change with co-array [%]  (<0 = better)', 'sum-SE change with co-array [%]  (>0 = better)'};
for row = 1:2
    for ca = 1:tNc
        subplot(2, 2, (row - 1) * 2 + ca); hold on; box on; grid on;
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
        xlabel('SNR [dB]'); ylabel(tYl7{row});
        title({sprintf('%s (K = %d)', en_caseName{ca}, en_K(ca)), ['identical:' ttl]});
        legend('Location', 'best');
    end
end

% ---- T8: localisation error vs SNR ----
figure('Name', 'T8 Localisation error vs SNR', 'Position', [30 30 1150 470]);
for ca = 1:tNc
    subplot(1, 2, ca); hold on; box on; grid on;
    for sc = 1:tNs
        plot(SNR_dB, squeeze(mLC(sc, ca, :)), [tLs{sc} tMk{sc}], 'Color', tCol{sc}, ...
             'MarkerFaceColor', tMfc{sc}, 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('median position error [m]');
    title(sprintf('%s (K = %d)', en_caseName{ca}, en_K(ca))); legend('Location', 'northeast');
end

% ---- T9: array hardware metrics ----
tAp = zeros(4, 1);  tUd = zeros(4, 1);  tW1 = zeros(4, 1);  tCl = zeros(4, 1);
for a = 1:4
    pa = EN_geo{tArrSch(a), 1}.pos(:);
    Dm = round((pa - pa.') / (lambda / 2));
    lg = unique(Dm(:));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
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
    set(gca, 'XTick', 1:4, 'XTickLabel', tArrNm);  xlim([0.4 4.6]);
    ylabel(tYl9{m});  title(tTt9{m});
end

% ---- console summary at 10 dB ----
for ca = 1:tNc
    fprintf('\n=== [CODE T] %s, K = %d, at %d dB (mean +- std. error over %d setups) ===\n', ...
            en_caseName{ca}, en_K(ca), SNR_dB(tSi), tSet);
    fprintf('%-28s %10s %22s %18s %10s\n', 'scheme', 'NMSE', 'BER (Cross-AP List)', 'sum SE', 'loc [m]');
    for sc = 1:tNs
        be = reshape(EN_be(tDet, sc, ca, tSi, :), [], 1);
        su = reshape(EN_sum(tDet, sc, ca, tSi, :), [], 1);
        fprintf('%-28s %10.2e %11.3e +- %8.1e %8.2f +- %5.2f %10.3f\n', en_names{sc}, mNM(sc, ca, tSi), ...
                mean(be), std(be) / sqrt(tSet), mean(su), std(su) / sqrt(tSet), mLC(sc, ca, tSi));
    end
    for pr = 1:2
        a = tPair(pr, 1);  b = tPair(pr, 2);
        Dd = EN_be(tDet, b, ca, :, :) - EN_be(tDet, a, ca, :, :);
        fprintf('  %s: co-array vs physical processing identical in %d%% of (SNR, setup) points\n', ...
                tPairNm{pr}, round(100 * mean(Dd(:) == 0)));
    end
end

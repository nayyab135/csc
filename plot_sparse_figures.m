%% ====================================================================
%  plot_sparse_figures.m  --  FIGURES ONLY (no simulation is run here)
%  --------------------------------------------------------------------
%  Loads a saved results .mat and draws the sparse-array comparison
%  figures R1..R9. It ADAPTS to what the file actually contains:
%    * Full set (R1..R9) when the file is the Code R 6-scheme run, i.e. it
%      has en_cape (CAPE on/off) and EN_Kbe (performance vs number of users).
%    * Otherwise it draws every figure the data supports and prints a
%      warning for the ones it cannot: R3/R4 (CAPE on/off) and R7 (vs-K)
%      are only possible from the Code R run, because those configurations
%      were never computed by the earlier runs.
%
%  USAGE: set RESULTS to your .mat, set lambda to the simulation wavelength,
%  then run. Re-run freely to rearrange plots -- no simulation is repeated.
%% ====================================================================
clear; close all; clc;

RESULTS = 'codeR_usparse_results.mat';   % <-- your results file
lambda  = 0.1;                           % <-- wavelength [m], must match the sim

assert(exist(RESULTS, 'file') == 2, 'Results file "%s" not found in this folder.', RESULTS);
S   = load(RESULTS);
has = @(n) isfield(S, n);

% ---- variables that every run saves ----
SNR_dB   = double(S.SNR_dB(:))';   nSNR  = numel(SNR_dB);
en_names = S.en_names;             nSchE = numel(en_names);
en_geom  = S.en_geom;
en_N0    = double(S.en_N0(:))';
en_K     = double(S.en_K(:))';     nCaseE = numel(en_N0);
en_caseName = {'Antenna-rich (K<N)', 'User-rich (K>N)'};
if nCaseE == 1, en_caseName = en_caseName(1); end
[~, si0] = min(abs(SNR_dB - 10));        % reference SNR ~ 10 dB
dXAP = 4;                                % Cross-AP List detector index

EN_nmse = S.EN_nmse;  EN_be = S.EN_be;  EN_sum = S.EN_sum;  EN_loc = S.EN_loc;
if has('EN_nf'), EN_nf = S.EN_nf; else, EN_nf = nan(nSchE, nCaseE, 1); end

is6 = (nSchE == 6) && has('en_cape');
if ~is6
    warning(['This .mat is NOT the Code R 6-scheme run (no en_cape, or nSchE ~= 6). ' ...
        'R3/R4 (CAPE on/off) and R7 (vs number of users) need that run and are skipped. ' ...
        'The array/geometry/localisation figures are drawn for the %d schemes present.'], nSchE);
end

% ---- averages over setups ----
mNM = mean(EN_nmse, 4);                  % nSchE x nCaseE x nSNR
mBE = mean(EN_be, 5);                    % 4 x nSchE x nCaseE x nSNR
mSU = mean(EN_sum, 5);
mLC = nan(nSchE, nCaseE, nSNR);
for scE = 1:nSchE
    for caE = 1:nCaseE
        for si = 1:nSNR
            v = EN_loc(scE, caE, si, :);  v = v(isfinite(v));
            if ~isempty(v), mLC(scE, caE, si) = median(v); end
        end
    end
end

% ---- rebuild geometries (positions are not saved in the .mat) ----
EN_geo = cell(nSchE, 1);
for scE = 1:nSchE, EN_geo{scE} = geom_targetN(en_geom{scE}, en_N0(1), lambda); end

cDense = [0.50 0.50 0.50];  cNest = [0.85 0.33 0.10];
cDCA   = [0.47 0.67 0.19];  cUni  = [0.00 0.45 0.74];
cPurp  = [0.49 0.18 0.56];

% ---- scheme role sets (adapt to the layout) ----
if is6
    S_dense = 1; S_nestL = 2; S_nest = 3; S_dca = 4; S_uniL = 5; S_uni = 6;
    arrSet = [S_nest S_dca S_uni];
    arrNm  = {'Nested + CAPE', 'Difference co-array + CAPE', 'Uniform sparse + CAPE (proposed)'};
    arrCol = {cNest, cDCA, cUni};
    locSet = [S_dense S_nest S_uni];
    locNm  = {'Dense ULA + CAPE', 'Nested + CAPE', 'Uniform + CAPE'};
    locCol = {cDense, cNest, cUni};
else
    arrSet = 1:nSchE;  arrNm = en_names;
    pal = {cDense, cNest, cDCA, cUni, cPurp, [0.30 0.30 0.30]};
    arrCol = pal(1:min(nSchE, 6));
    locSet = arrSet;  locNm = arrNm;  locCol = arrCol;
end
arrMk = {'-s', '-^', '-o', '-d', '-v', '-p'};

% ---- distinct physical geometries for R1/R2/R8 ----
[~, ia] = unique(en_geom, 'stable');  gIdx = ia(:)';  nG = numel(gIdx);
gNm = cell(1, nG);  gCol = cell(1, nG);  gpal = {cDense, cNest, cUni, cDCA, cPurp};
for g = 1:nG
    gNm{g}  = pretty_geom(en_geom{gIdx(g)});
    gCol{g} = gpal{min(g, numel(gpal))};
end

%% ---- R1: physical array and its difference co-array ----
figure('Name', 'R1 Geometry and difference co-array', 'Position', [30 30 1100 620]);
for g = 1:nG
    G = EN_geo{gIdx(g)};  CA = coarray_setup(G.pos, lambda, 'diff');
    subplot(nG, 1, g); hold on; box on; grid on;
    plot(G.pos, ones(G.N, 1), 'o', 'Color', gCol{g}, 'MarkerFaceColor', gCol{g}, 'MarkerSize', 5);
    plot(CA.pos_v, 0.6 * ones(CA.Nv, 1), '.', 'Color', gCol{g}, 'MarkerSize', 6);
    ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
    xlabel('position along array [m]');
    title(sprintf('%s:  N = %d,  aperture = %.2f m,  uDOF(diff) = %d,  Rayleigh = %.0f m', ...
        gNm{g}, G.N, G.aperture, CA.Nv, 2 * G.aperture^2 / lambda), 'Interpreter', 'none');
end

%% ---- R2: beampattern (grating lobes) ----
figure('Name', 'R2 Beampattern and grating lobes', 'Position', [30 30 900 480]);
hold on; box on; grid on;  sgb = linspace(-1, 1, 2001);
for g = 1:nG
    pp = EN_geo{gIdx(g)}.pos(:);
    Bb = abs(sum(exp(1j * 2 * pi * (sgb.' * pp.') / lambda), 2)) / numel(pp);
    plot(sgb, 20 * log10(max(Bb, 1e-4)), 'Color', gCol{g}, 'LineWidth', 1.5, 'DisplayName', gNm{g});
end
xlabel('sin(\theta)'); ylabel('Normalised beampattern [dB]'); ylim([-40 0]);
title('Beampattern: wide uniform spacing gives grating lobes (CAPE resolves them)');
legend('Location', 'south');

%% ---- R3 & R4: CAPE benefit (only from the 6-scheme run) ----
if is6
    figure('Name', 'R3 CAPE benefit - channel NMSE', 'Position', [30 30 1150 470]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        semilogy(SNR_dB, squeeze(mNM(S_nestL, caE, :)), '--s', 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested, local est.');
        semilogy(SNR_dB, squeeze(mNM(S_nest,  caE, :)), '-s',  'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
        semilogy(SNR_dB, squeeze(mNM(S_uniL,  caE, :)), '--o', 'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform, local est.');
        semilogy(SNR_dB, squeeze(mNM(S_uni,   caE, :)), '-o',  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform + CAPE');
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
        title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
    end
    figure('Name', 'R4 CAPE benefit - BER', 'Position', [30 30 1150 470]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nestL, caE, :)), 1e-6), '--s', 'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested, local est.');
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_nest,  caE, :)), 1e-6), '-s',  'Color', cNest, 'LineWidth', 1.6, 'DisplayName', 'Nested + CAPE');
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uniL,  caE, :)), 1e-6), '--o', 'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform, local est.');
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, S_uni,   caE, :)), 1e-6), '-o',  'Color', cUni,  'LineWidth', 1.6, 'DisplayName', 'Uniform + CAPE');
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
        title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
    end
else
    fprintf('[skip] R3/R4 (CAPE on/off) need the Code R 6-scheme run.\n');
end

%% ---- R5: array comparison, BER ----
figure('Name', 'R5 Array comparison - BER', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, nCaseE, caE); hold on; box on; grid on;
    for j = 1:numel(arrSet)
        semilogy(SNR_dB, max(squeeze(mBE(dXAP, arrSet(j), caE, :)), 1e-6), arrMk{mod(j-1,numel(arrMk))+1}, ...
            'Color', arrCol{j}, 'LineWidth', 1.6, 'DisplayName', arrNm{j});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (Cross-AP List)');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'southwest');
end

%% ---- R6: array comparison, sum SE ----
figure('Name', 'R6 Array comparison - sum SE', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, nCaseE, caE); hold on; box on; grid on;
    for j = 1:numel(arrSet)
        plot(SNR_dB, squeeze(mSU(dXAP, arrSet(j), caE, :)), arrMk{mod(j-1,numel(arrMk))+1}, ...
            'Color', arrCol{j}, 'LineWidth', 1.6, 'DisplayName', arrNm{j});
    end
    xlabel('SNR [dB]'); ylabel('Sum SE [bit/s/Hz]');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'northwest');
end

%% ---- R7: performance vs number of users K (only if the sweep was saved) ----
if has('EN_Kbe') && has('Ksw')
    Ksw = double(S.Ksw(:))';  mKbe = mean(S.EN_Kbe, 3);  mKsum = mean(S.EN_Ksum, 3);
    kNm = {'Nested + CAPE', 'Difference co-array + CAPE', 'Uniform sparse + CAPE'};
    kCol = {cNest, cDCA, cUni};  kMk = {'-s', '-^', '-o'};
    figure('Name', 'R7 Performance vs number of users K', 'Position', [30 30 1150 470]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for j = 1:size(mKbe, 1)
        semilogy(Ksw, max(mKbe(j, :), 1e-6), kMk{min(j,3)}, 'Color', kCol{min(j,3)}, 'LineWidth', 1.6, 'DisplayName', kNm{min(j,3)});
    end
    set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
    title(sprintf('BER vs K at %d dB', SNR_dB(si0))); legend('Location', 'northwest');
    subplot(1, 2, 2); hold on; box on; grid on;
    for j = 1:size(mKsum, 1)
        plot(Ksw, mKsum(j, :), kMk{min(j,3)}, 'Color', kCol{min(j,3)}, 'LineWidth', 1.6, 'DisplayName', kNm{min(j,3)});
    end
    xlabel('number of users K'); ylabel('Sum SE [bit/s/Hz]');
    title(sprintf('Sum SE vs K at %d dB', SNR_dB(si0))); legend('Location', 'northwest');
else
    fprintf('[skip] R7 (performance vs K) needs the Code R run that saves EN_Kbe/Ksw.\n');
end

%% ---- R8: array hardware metrics ----
figure('Name', 'R8 Array hardware metrics', 'Position', [30 30 1200 440]);
apv = zeros(nG, 1); udv = zeros(nG, 1); w1v = zeros(nG, 1); clv = zeros(nG, 1);
for g = 1:nG
    pg = EN_geo{gIdx(g)}.pos;  [clv(g), wv, ~] = sparse_metrics(pg, lambda);  w1v(g) = wv(1);
    apv(g) = EN_geo{gIdx(g)}.aperture;  cas = coarray_setup(pg, lambda, 'diff');  udv(g) = cas.Nv;
end
subplot(1, 4, 1); bar(apv); set(gca, 'XTick', 1:nG, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('aperture [m]'); title('Aperture (resolution)'); grid on;
subplot(1, 4, 2); bar(udv); set(gca, 'XTick', 1:nG, 'XTickLabel', gNm, 'XTickLabelRotation', 20, 'YScale', 'log'); ylabel('uDOF'); title('Co-array uDOF'); grid on;
subplot(1, 4, 3); bar(w1v); set(gca, 'XTick', 1:nG, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('w(1)'); title('Closest-pair count (coupling)'); grid on;
subplot(1, 4, 4); bar(clv); set(gca, 'XTick', 1:nG, 'XTickLabel', gNm, 'XTickLabelRotation', 20); ylabel('coupling leakage'); title('Mutual coupling'); grid on;

%% ---- R9: localisation error vs SNR ----
figure('Name', 'R9 Localisation error vs SNR', 'Position', [30 30 1150 470]);
for caE = 1:nCaseE
    subplot(1, nCaseE, caE); hold on; box on; grid on;
    for j = 1:numel(locSet)
        semilogy(SNR_dB, squeeze(mLC(locSet(j), caE, :)), arrMk{mod(j-1,numel(arrMk))+1}, ...
            'Color', locCol{j}, 'LineWidth', 1.6, 'DisplayName', locNm{j});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Median position error [m]');
    title(sprintf('%s (K = %d)', en_caseName{caE}, en_K(caE))); legend('Location', 'northeast');
end

fprintf('\nDone. Drew figures from "%s" (%d schemes, is6=%d).\n', RESULTS, nSchE, is6);

%% ==================== helper functions (geometry only) ====================
function s = pretty_geom(g)
    switch lower(g)
        case 'ula',       s = 'Dense ULA';
        case 'nested',    s = 'Nested';
        case 'usparse',   s = 'Uniform sparse';
        case 'symnested', s = 'Symmetric';
        case 'eca',       s = 'ECA';
        otherwise,        s = g;
    end
end
function G = geom_targetN(type, N0, lambda)
% [CODE G] build each array type with EXACTLY N0 PHYSICAL elements on the
% lambda/2 grid, so every scheme is compared at equal hardware. The virtual
% (co-)array that results is whatever each geometry yields; only the physical
% count is fixed. 'eca'/'symnested' are trimmed to N0 keeping both ends;
% 'usparse' is a uniform wide-spacing array (new idea, scheme 4).
    d0 = lambda / 2;
    switch lower(type)
        case 'ula'
            n = (0:N0 - 1)';
        case 'nested'
            N1 = floor(N0 / 2);  N2 = N0 - N1;            % exactly N0, no overlap
            n = unique([(1:N1), (1:N2) * (N1 + 1)])';
        case 'symnested'
            h  = ceil(N0 / 2);  a1 = floor(h / 2);  a2 = h - a1;
            half = unique([(1:a1), (1:a2) * (a1 + 1)]);   % a nested half (positive side)
            n = unique([-fliplr(half), half])';           % mirror -> symmetric
            n = force_count(n, N0);
        case 'eca'
            M = 5;  Nc = 4;                                % a generous coprime base
            ca = unique([(0:M - 1) * Nc, (0:Nc - 1) * M]);
            up = [];  for cc = 0:3, up = [up, cc * M * Nc + ca]; end  %#ok<AGROW>
            n = unique([up, -up])';
            n = force_count(n, N0);
        case 'usparse'
            % [CODE R] UNIFORM SPARSE ARRAY for near-field cell-free.
            %  Design (idea Steps 1-2): N0 elements at a wide UNIFORM spacing
            %  eta*d0. eta is set so the aperture D = (N0-1)*eta*d0 puts the
            %  whole coverage deep in the near field (Rayleigh 2D^2/lambda >>
            %  r_max) AND every gap >= lambda (eta >= 2), so w(1)=w(2)=0 ->
            %  minimum mutual coupling. eta = 9 gives D ~ 14 m and Rayleigh
            %  ~ 3.9 km at lambda = 0.1 m, N0 = 32 (same aperture as the
            %  nested array, for a fair equal-hardware comparison). The
            %  grating-lobe ambiguity this wide spacing creates is NOT fixed
            %  on the array (no dense sub-array, no hole-free co-array); it is
            %  resolved by non-coherent cross-AP CAPE fusion. That division of
            %  labour -- resolution on the array, ambiguity on the network --
            %  is the contribution, so this array carries en_idea = 0.
            eta_us = 9;
            n = (0:N0 - 1)' * eta_us;
        otherwise
            error('geom_targetN: unknown type %s', type);
    end
    p = n * d0;
    G.pos = p - (max(p) + min(p)) / 2;
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    G.sparsity = G.aperture / ((G.N - 1) * d0);
    G.min_gap  = min(diff(sort(G.pos)));
end

function [cl, wv, pslv] = sparse_metrics(pos, lambda)
% [CODE P] The design metrics reported in the sparse-array literature
% (nested / coprime / MISC / IMISC / EMISC papers):
%   cl   - coupling leakage, the standard mutual-coupling measure
%          ||U - diag(U)||_F / ||U||_F with the banded model u_k = 1/k
%   wv   - the weight function values w(1), w(2), w(3), i.e. how many antenna
%          pairs sit 1, 2 and 3 spacings apart. These dominate the coupling
%          matrix, which is why EMISC drives w(1) down to 1.
%   pslv - peak sidelobe of the beampattern (ambiguity / grating lobes)
    d0 = lambda / 2;
    Dm = round((pos(:) - pos(:).') / d0);
    wv = [sum(Dm(:) == 1), sum(Dm(:) == 2), sum(Dm(:) == 3)];
    G = 100;  Na = numel(pos);  U = zeros(Na);
    for a = 1:Na
        for b = 1:Na
            kk = abs(Dm(a, b));
            if kk == 0, U(a, b) = 1; elseif kk <= G, U(a, b) = 1 / kk; end
        end
    end
    cl = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
    sg = linspace(-1, 1, 4001);
    Bm = abs(sum(exp(1j * 2 * pi * (sg.' * pos(:).') / lambda), 2)) / Na;
    [~, i0] = max(Bm);  i = i0;
    while i + 1 <= numel(Bm) && Bm(i + 1) < Bm(i), i = i + 1; end
    j = i0;
    while j - 1 >= 1 && Bm(j - 1) < Bm(j), j = j - 1; end
    Bm(j:i) = 0;  pslv = max(Bm);
end

function n = force_count(n, N0)
% Force an integer-position set to EXACTLY N0 elements: drop interior points
% evenly if too many (keep both ends to preserve aperture), extend the grid if
% too few. Guaranteed to return exactly N0 distinct integers.
    n = unique(sort(n(:)));
    if numel(n) > N0
        idx = unique(round(linspace(1, numel(n), N0)));
        n = n(idx);
        while numel(n) < N0
            full = (min(n):max(n))';
            cand = setdiff(full, n);
            if isempty(cand), n = unique([n; max(n) + 1]); else n = unique([n; cand(1)]); end
        end
    elseif numel(n) < N0
        n = [n; (max(n) + (1:(N0 - numel(n)))')];
    end
    n = unique(sort(n(:)));
end

function CA = coarray_setup(pos, lambda, mode)
% Co-array of a linear array on the lambda/2 grid.
%   mode = 'diff' (default): lags x_m - x_n  -> DIFFERENCE co-array (the standard
%          Pinto building block; comes from the Hermitian covariance R = E[y y^H]).
%   mode = 'sum'          : lags x_m + x_n  -> SUM co-array (used by Idea 3; comes
%          from the complementary moment Rt = E[y y^T]). For near-field phase and
%          deterministic pilots the sum lags carry extra (range) structure.
    if nargin < 3, mode = 'diff'; end
    d0 = lambda / 2;
    pos = pos(:);
    % Lags come from position SUMS/DIFFERENCES: the array is centred, so pos/d0
    % can be half-integer and rounding it would corrupt the lag map, while sums
    % and differences of centred positions are exact multiples of d0.
    if strcmpi(mode, 'sum')
        Dm = round((pos + pos.') / d0);     % SUM lags
    else
        Dm = round((pos - pos.') / d0);     % DIFFERENCE lags
    end
    CA.mode = lower(mode);
    CA.lagvec = Dm(:);                  % one lag per entry of vec(R) / vec(Rt)
    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);
    CA.cont = longest_contig(CA.uniq);  % longest hole-free virtual ULA
    [~, CA.loc] = ismember(CA.cont, CA.uniq);
    CA.Nv    = numel(CA.cont);          % virtual sensors used
    CA.Nsub  = floor((CA.Nv + 1) / 2);  % spatial-smoothing subarray length
    CA.DOF   = numel(CA.uniq);          % distinct lags = total virtual sensors
    CA.Nphys = numel(pos);
    CA.d0    = d0;
    CA.lambda = lambda;
    CA.pos_v = CA.cont * d0;
end

function run = longest_contig(u)
% Longest run of consecutive integers in u (the hole-free virtual ULA).
    u = sort(u(:));
    bl = u(1);  bh = u(1);  cs = u(1);
    for i = 2:numel(u)
        if u(i) ~= u(i - 1) + 1
            if (u(i - 1) - cs) > (bh - bl), bl = cs;  bh = u(i - 1);  end
            cs = u(i);
        end
    end
    if (u(end) - cs) > (bh - bl), bl = cs;  bh = u(end);  end
    run = (bl:bh)';
end

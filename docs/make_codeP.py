#!/usr/bin/env python3
# Derive Code_P_Nested.m from Code_O_Final.m.
#   1. PART 2: nested array replaces the ECA for all sparse schemes.
#   2. PART 2: 32 physical antennas per AP; Case 1 K = 8, Case 2 K = 40.
#   3. PART 2: adds the performance figures the sparse-array literature uses
#      (physical/virtual geometry, design metrics, beampattern, uDOF vs N,
#      localisation RMSE vs SNR).
import io

SRC = "/home/user/csc/Code_O_Final.m"
DST = "/home/user/csc/Code_P_Nested.m"
with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ---- 1. nested instead of ECA -----------------------------------------------
rep("""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed (sym near-field)'};
en_geom     = {'ula', 'eca', 'nested', 'symnested'};""",
"""en_names    = {'Dense ULA', 'Nested sparse array', ...
               'Nested + difference co-array', 'Symmetric near-field co-array'};
en_geom     = {'ula', 'nested', 'nested', 'symnested'};""")

# ---- 2. users per case -------------------------------------------------------
rep("en_K        = [8 64];       % [CODE O] same two cases: K<N then K>N",
    "en_K        = [8 40];       % [CODE P] Case 1 K=8 (K<N), Case 2 K=40 (K>N)")

# ---- 3a. metrics helper ------------------------------------------------------
rep("""function n = force_count(n, N0)""",
"""function [cl, wv, pslv] = sparse_metrics(pos, lambda)
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

function n = force_count(n, N0)""")

# ---- 3b. the figures ---------------------------------------------------------
rep("""    %% ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\\n=== [CODE I] Array comparison, %s, at %d dB (N=%d, K=%d) ===\\n', ...""",
"""    %% ========== [CODE P] SPARSE-ARRAY FIGURES USED IN THE LITERATURE ==========
    % SPA-1: physical elements and their difference co-array (the opening
    % figure of essentially every sparse-array paper).
    figure('Name', 'SPA-1 Physical array and its difference co-array', 'Position', [30 30 1250 820]);
    for scE = 1:nSchE
        subplot(nSchE, 1, scE); hold on; box on; grid on;
        pp = EN_geo{scE, 1}.pos;
        plot(pp, ones(numel(pp), 1), 'o', 'Color', aC{scE}, ...
             'MarkerFaceColor', aC{scE}, 'MarkerSize', 5);
        plot(EN_CA{scE, 1}.pos_v, 0.6 * ones(EN_CA{scE, 1}.Nv, 1), '.', ...
             'Color', aC{scE}, 'MarkerSize', 6);
        ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
        xlabel('position along the array [m]');
        title(sprintf('%s:  N = %d physical,  %d virtual sensors (hole-free %d),  aperture %.2f m', ...
              en_names{scE}, EN_geo{scE, 1}.N, EN_CA{scE, 1}.DOF, EN_CA{scE, 1}.Nv, ...
              EN_geo{scE, 1}.aperture), 'Interpreter', 'none');
    end

    % SPA-2: the design metrics those papers tabulate
    figure('Name', 'SPA-2 Sparse-array design metrics', 'Position', [30 30 1300 820]);
    clv = zeros(nSchE, 1);  wvv = zeros(nSchE, 3);  pslw = zeros(nSchE, 1);
    apv = zeros(nSchE, 1);  rayv = zeros(nSchE, 1);
    udD = zeros(nSchE, 1);  udS = zeros(nSchE, 1);
    for scE = 1:nSchE
        pg = EN_geo{scE, 1}.pos;
        [clv(scE), wvv(scE, :), pslw(scE)] = sparse_metrics(pg, lambda);
        apv(scE)  = EN_geo{scE, 1}.aperture;
        rayv(scE) = 2 * apv(scE) ^ 2 / lambda;
        udD(scE)  = EN_CA{scE, 1}.Nv;
        cas_p     = coarray_setup(pg, lambda, 'sum');
        udS(scE)  = cas_p.Nv;
    end
    lblp = @() set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, ...
                   'XTickLabelRotation', 18, 'TickLabelInterpreter', 'none');
    subplot(2, 3, 1); bar(1:nSchE, [udD(:), udS(:)]); lblp(); grid on;
    ylabel('lags'); title('uDOF: difference (angle) and sum (range)');
    legend({'uDOF diff', 'uDOF sum'}, 'Location', 'northwest');
    subplot(2, 3, 2); bar(1:nSchE, clv); lblp(); grid on;
    ylabel('coupling leakage'); title('Mutual coupling (lower is better)');
    subplot(2, 3, 3); bar(1:nSchE, wvv); lblp(); grid on;
    ylabel('antenna pairs'); title('Weight function w(1), w(2), w(3)');
    legend({'w(1)', 'w(2)', 'w(3)'}, 'Location', 'northwest');
    subplot(2, 3, 4); bar(1:nSchE, pslw); lblp(); grid on;
    ylabel('peak sidelobe'); title('Ambiguity (peak sidelobe level)');
    subplot(2, 3, 5); bar(1:nSchE, apv); lblp(); grid on;
    ylabel('aperture [m]'); title('Physical aperture');
    subplot(2, 3, 6); bar(1:nSchE, rayv); lblp(); grid on; set(gca, 'YScale', 'log');
    ylabel('Rayleigh distance [m]'); title('Near-field reach  2D^2/\\lambda');

    % SPA-3: beampattern -- where the grating lobes are
    figure('Name', 'SPA-3 Beampattern', 'Position', [30 30 1100 520]);
    hold on; box on; grid on;
    sgb = linspace(-1, 1, 2001);
    for scE = 1:nSchE
        pp = EN_geo{scE, 1}.pos(:);
        Bb = abs(sum(exp(1j * 2 * pi * (sgb.' * pp.') / lambda), 2)) / numel(pp);
        plot(sgb, 20 * log10(max(Bb, 1e-4)), 'Color', aC{scE}, ...
             'LineWidth', 1.5, 'DisplayName', en_names{scE});
    end
    xlabel('sin(\\theta)'); ylabel('Normalised beampattern [dB]'); ylim([-40 0]);
    title('Beampattern: main lobe width and grating lobes'); legend('Location', 'south');

    % SPA-4: uDOF growth with the number of physical sensors
    figure('Name', 'SPA-4 uDOF vs number of physical antennas', 'Position', [30 30 1100 520]);
    hold on; box on; grid on;
    Nsw = 8:4:48;
    for scE = 1:nSchE
        udn = zeros(size(Nsw));
        for ii = 1:numel(Nsw)
            Gt = geom_targetN(en_geom{scE}, Nsw(ii), lambda);
            ct = coarray_setup(Gt.pos, lambda, 'diff');
            udn(ii) = ct.Nv;
        end
        plot(Nsw, udn, '-o', 'Color', aC{scE}, 'LineWidth', 1.6, 'DisplayName', en_names{scE});
    end
    xlabel('number of physical antennas N'); ylabel('uDOF (hole-free virtual sensors)');
    title('Degrees of freedom vs array size'); legend('Location', 'northwest');

    % SPA-5: localisation RMSE vs SNR (the DOA-paper form of accuracy)
    figure('Name', 'SPA-5 Localisation error vs SNR', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, squeeze(mLC(scE, caE, :)), '-o', 'Color', aC{scE}, ...
                     'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Median position error [m]');
        title(en_caseName{caE}); legend('Location', 'northeast');
    end

    %% ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\\n=== [CODE I] Array comparison, %s, at %d dB (N=%d, K=%d) ===\\n', ...""")

rep("save('codeO_final_results.mat'", "save('codeP_nested_results.mat'")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code P written:", t.count("\n") + 1, "lines")

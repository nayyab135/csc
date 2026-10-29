#!/usr/bin/env python3
# Derive Code_I_Comparison.m from Code_H_CoArray_Proposed.m.
# Consolidated comparison with two clean figure families:
#  (A) Detector comparison: Cross-AP List vs Linear MMSE / Hard-SIC / List-SIC
#      (BER, sum rate, per-user SE CDF) on the proposed array.
#  (B) Array comparison (equal physical N, curves vs SNR): Dense ULA, ECA+CAPE,
#      standard difference co-array + vectorization, and the proposed symmetric
#      near-field co-array -- BER, sum rate, NMSE, complexity, fronthaul,
#      per-user SE / localisation CDFs, geometry + near-field fraction.
# The old unequal-N co_study and the Code D PROP study are switched off.
import io

SRC = "/home/user/csc/Code_H_CoArray_Proposed.m"
DST = "/home/user/csc/Code_I_Comparison.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:160], c)
    t = t.replace(old, new)

# ---- turn off the Code D PROP study (superseded by the array comparison) ----
rep("d_study     = true;", "d_study     = false;   % [CODE I] PROP figures off (superseded)")

# ---- EN comment block -> four-array comparison ----
rep(
"%  geometry builds is left to whatever it yields. Two schemes only:\n"
"%    1. ECA + CAPE (original)              eca geometry, NO co-array  (co_idea 0)\n"
"%    2. Proposed (symmetric nested array) difference co-array + CAPE (co_idea 1)\n"
"%  The three other ideas (cross-AP, sum+diff, hybrid-field) all rested on the\n"
"%  same difference-co-array DOF gain and performed identically, so the single\n"
"%  best-justified one is kept: a SYMMETRIC nested array keeps the difference\n"
"%  co-array valid in the near field, and CAPE is used unchanged (only the\n"
"%  physical array is replaced; the co-array just feeds CAPE extra candidates).",
"%  geometry builds is left to whatever it yields. FOUR arrays, all at equal\n"
"%  physical N, all with CAPE channel estimation:\n"
"%    1. Dense ULA + CAPE          compact uniform array, NO co-array (co_idea 0)\n"
"%    2. ECA + CAPE                extended coprime sparse array      (co_idea 0)\n"
"%    3. Std diff co-array + vec   two-level nested + vectorization    (co_idea 1)\n"
"%    4. Proposed (sym near-field) symmetric nested co-array          (co_idea 1)\n"
"%  Detector comparison (Cross-AP List vs Linear MMSE / Hard-SIC / List-SIC) is\n"
"%  plotted on the proposed array. CAPE is used unchanged for every array; the\n"
"%  co-array only feeds CAPE extra candidates.")

# ---- four-scheme definition ----
rep(
"en_names    = {'ECA + CAPE (original)', 'Proposed: sym-nested co-array + CAPE'};\n"
"en_geom     = {'eca', 'symnested'};\n"
"en_idea     = [0 1];",
"en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...\n"
"               'Std diff co-array + vec', 'Proposed (sym near-field)'};\n"
"en_geom     = {'ula', 'eca', 'nested', 'symnested'};\n"
"en_idea     = [0 0 1 1];")

# ---- data-array declarations: add time / fronthaul / NF fraction / per-user SE ----
rep(
"EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);\n"
"EN_be   = zeros(4, nSchE, nCaseE, nSNR, nSetE);\n"
"EN_sum  = zeros(4, nSchE, nCaseE, nSNR, nSetE);\n"
"EN_loc  = nan(nSchE, nCaseE, nSNR, nSetE);",
"EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);\n"
"EN_be   = zeros(4, nSchE, nCaseE, nSNR, nSetE);\n"
"EN_sum  = zeros(4, nSchE, nCaseE, nSNR, nSetE);\n"
"EN_loc  = nan(nSchE, nCaseE, nSNR, nSetE);\n"
"EN_time = zeros(nSchE, nCaseE, nSetE);              % mean run time (complexity proxy)\n"
"EN_fh   = zeros(nSchE, nCaseE, nSetE);              % fused APs per user (fronthaul proxy)\n"
"EN_nf   = zeros(nSchE, nCaseE, nSetE);              % near-field link fraction\n"
"EN_seD  = nan(4, 16, nSchE, nCaseE, nSetE);         % per-user SE at ref SNR (4 detectors)")

# ---- parfor temporaries: add the new accumulators ----
rep(
"    s_nmse = zeros(nSchE, nCaseE, nSNR);\n"
"    s_be   = zeros(4, nSchE, nCaseE, nSNR);\n"
"    s_sum  = zeros(4, nSchE, nCaseE, nSNR);\n"
"    s_loc  = nan(nSchE, nCaseE, nSNR);",
"    s_nmse = zeros(nSchE, nCaseE, nSNR);\n"
"    s_be   = zeros(4, nSchE, nCaseE, nSNR);\n"
"    s_sum  = zeros(4, nSchE, nCaseE, nSNR);\n"
"    s_loc  = nan(nSchE, nCaseE, nSNR);\n"
"    s_time = zeros(nSchE, nCaseE);\n"
"    s_fh   = zeros(nSchE, nCaseE);\n"
"    s_nf   = zeros(nSchE, nCaseE);\n"
"    s_seD  = nan(4, 16, nSchE, nCaseE);")

# ---- record near-field fraction per scheme (independent of SNR/seed) ----
rep(
"            ST = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);\n",
"            ST = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);\n"
"            s_nf(scE, caE) = mean(ST.NF(:));\n")

# ---- collect time / fronthaul / per-user SE inside the SNR loop ----
rep(
"                s_nmse(scE, caE, si) = O.nmse;\n"
"                s_be(:, scE, caE, si) = O.be;\n"
"                s_sum(:, scE, caE, si) = O.sum;\n"
"                if ~isempty(O.loc), s_loc(scE, caE, si) = median(O.loc); end",
"                s_nmse(scE, caE, si) = O.nmse;\n"
"                s_be(:, scE, caE, si) = O.be;\n"
"                s_sum(:, scE, caE, si) = O.sum;\n"
"                if ~isempty(O.loc), s_loc(scE, caE, si) = median(O.loc); end\n"
"                s_time(scE, caE) = s_time(scE, caE) + O.time;\n"
"                s_fh(scE, caE)   = s_fh(scE, caE) + O.nfus;\n"
"                if si == snr_ref_idx\n"
"                    s_seD(:, 1:size(O.se, 2), scE, caE) = O.se;\n"
"                end")

# ---- reduction: store the new accumulators ----
rep(
"    EN_nmse(:, :, :, ns) = s_nmse;  EN_be(:, :, :, :, ns) = s_be;\n"
"    EN_sum(:, :, :, :, ns) = s_sum;    EN_loc(:, :, :, ns) = s_loc;",
"    EN_nmse(:, :, :, ns) = s_nmse;  EN_be(:, :, :, :, ns) = s_be;\n"
"    EN_sum(:, :, :, :, ns) = s_sum;    EN_loc(:, :, :, ns) = s_loc;\n"
"    EN_time(:, :, ns) = s_time / nSNR;  EN_fh(:, :, ns) = s_fh / nSNR;\n"
"    EN_nf(:, :, ns) = s_nf;  EN_seD(:, :, :, :, ns) = s_seD;")

# ---- replace the whole figure block (span from the figure guard to save+end) ----
START = ("if eqN_study\n"
         "    eC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.30 0.75 0.93]};\n"
         "    si0 = snr_ref_idx;")
END = ("    save('codeH_proposed_results.mat', 'SNR_dB', 'en_names', 'en_geom', 'en_N0', 'en_K', ...\n"
       "         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc');\n"
       "end")
i0 = t.index(START)
iE = t.index(END, i0) + len(END)

NEWFIG = r'''if eqN_study
    aC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56]};
    dC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56]};
    detName = {'Linear MMSE', 'Hard-SIC', 'List-SIC', 'Cross-AP List (proposed)'};
    si0 = snr_ref_idx;
    mNM = mean(EN_nmse, 4);        % nSchE x nCaseE x nSNR
    mBE = mean(EN_be, 5);          % 4 x nSchE x nCaseE x nSNR
    mSU = mean(EN_sum, 5);
    mTIME = mean(EN_time, 3);      % nSchE x nCaseE
    mFH   = mean(EN_fh, 3);
    mNF   = mean(EN_nf, 3);
    mLC = nan(nSchE, nCaseE, nSNR);
    for scE = 1:nSchE
        for caE = 1:nCaseE
            for si = 1:nSNR
                v = EN_loc(scE, caE, si, :);  v = v(isfinite(v));
                if ~isempty(v), mLC(scE, caE, si) = median(v); end
            end
        end
    end
    prop = nSchE;                  % proposed array is the last scheme
    c1   = 1;                      % antenna-rich case = main scenario

    %% ---------------- DETECTOR COMPARISON (proposed array, main scenario) ----------------
    % DET-1: BER vs SNR, four detectors
    figure('Name', 'DET-1 BER: Cross-AP List vs standard detectors', 'Position', [30 30 760 560]);
    hold on; box on; grid on;
    for d = 1:4
        semilogy(SNR_dB, max(squeeze(mBE(d, prop, c1, :)), 1e-6), '-o', 'Color', dC{d}, ...
                 'LineWidth', 1.6, 'DisplayName', detName{d});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
    title('Detector comparison (proposed array, hybrid NF/FF)'); legend('Location', 'southwest');

    % DET-2: sum rate vs SNR, four detectors
    figure('Name', 'DET-2 Sum rate: Cross-AP List vs standard detectors', 'Position', [30 30 760 560]);
    hold on; box on; grid on;
    for d = 1:4
        plot(SNR_dB, squeeze(mSU(d, prop, c1, :)), '-o', 'Color', dC{d}, ...
             'LineWidth', 1.6, 'DisplayName', detName{d});
    end
    xlabel('SNR [dB]'); ylabel('Sum rate [bit/s/Hz]');
    title('Detector comparison (proposed array)'); legend('Location', 'northwest');

    % DET-3: per-user SE CDF at reference SNR, four detectors
    figure('Name', 'DET-3 Per-user SE CDF: detectors', 'Position', [30 30 760 560]);
    hold on; box on; grid on;
    for d = 1:4
        v = EN_seD(d, :, prop, c1, :);  v = v(isfinite(v));  v = sort(v(:));
        if ~isempty(v)
            plot(v, (1:numel(v)) / numel(v), '-', 'Color', dC{d}, 'LineWidth', 1.6, 'DisplayName', detName{d});
        end
    end
    xlabel('Per-user SE [bit/s/Hz]'); ylabel('CDF');
    title(sprintf('Per-user SE at %d dB (proposed array)', SNR_dB(si0))); legend('Location', 'southeast');

    %% ---------------- ARRAY COMPARISON (Cross-AP List detector, both regimes) ----------------
    % ARR-1: BER vs SNR, four arrays, per regime
    figure('Name', 'ARR-1 BER: dense vs ECA vs std co-array vs proposed', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, max(squeeze(mBE(4, scE, caE, :)), 1e-6), '-o', 'Color', aC{scE}, ...
                     'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'southwest');
    end

    % ARR-2: sum rate vs SNR, four arrays
    figure('Name', 'ARR-2 Sum rate: four arrays', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            plot(SNR_dB, squeeze(mSU(4, scE, caE, :)), '-o', 'Color', aC{scE}, ...
                 'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        xlabel('SNR [dB]'); ylabel('Sum rate [bit/s/Hz]');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'northwest');
    end

    % ARR-3: channel NMSE vs SNR, four arrays
    figure('Name', 'ARR-3 NMSE: four arrays', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, squeeze(mNM(scE, caE, :)), '-o', 'Color', aC{scE}, ...
                     'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'southwest');
    end

    % ARR-4: computational cost (run time) and fronthaul (fused APs)
    figure('Name', 'ARR-4 Complexity and fronthaul', 'Position', [30 30 1250 480]);
    subplot(1, 2, 1); hold on; box on; grid on;
    bar(1:nSchE, mTIME);
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, 'TickLabelInterpreter', 'none');
    ylabel('Mean run time per setup [s]'); title('Computational cost (empirical)');
    legend(en_caseName, 'Location', 'northwest');
    subplot(1, 2, 2); hold on; box on; grid on;
    bar(1:nSchE, mFH);
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, 'TickLabelInterpreter', 'none');
    ylabel('Mean fused APs per user'); title('Fronthaul load (CAPE pilot fusion)');
    legend(en_caseName, 'Location', 'northwest');

    % ARR-5: per-user SE CDF and localisation CDF across arrays (ref SNR, main scenario)
    figure('Name', 'ARR-5 Per-user SE and localisation CDF', 'Position', [30 30 1250 480]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for scE = 1:nSchE
        v = EN_seD(4, :, scE, c1, :);  v = v(isfinite(v));  v = sort(v(:));
        if ~isempty(v)
            plot(v, (1:numel(v)) / numel(v), '-', 'Color', aC{scE}, 'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
    end
    xlabel('Per-user SE [bit/s/Hz]'); ylabel('CDF');
    title(sprintf('Per-user SE at %d dB (Cross-AP List)', SNR_dB(si0))); legend('Location', 'southeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    for scE = 1:nSchE
        v = EN_loc(scE, c1, :, :);  v = v(isfinite(v));  v = sort(v(:));
        if ~isempty(v)
            plot(v, (1:numel(v)) / numel(v), '-', 'Color', aC{scE}, 'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
    end
    xlabel('Localisation error [m]'); ylabel('CDF');
    title('CAPE localisation error'); legend('Location', 'southeast');

    % ARR-6: equal-N geometry, virtual sensors, near-field fraction
    figure('Name', 'ARR-6 Equal-N geometry, DOF, near-field fraction', 'Position', [30 30 1250 480]);
    subplot(1, 2, 1); hold on; box on; grid on;
    Nphys = cellfun(@(g) g.N, EN_geo(:, 1));
    Nvirt = cellfun(@(c) c.DOF, EN_CA(:, 1));
    bar(1:nSchE, [Nphys(:), Nvirt(:)]); set(gca, 'YScale', 'log');
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, 'TickLabelInterpreter', 'none');
    ylabel('count'); title(sprintf('Equal physical N=%d: antennas vs virtual sensors', en_N0(1)));
    legend({'physical antennas', 'virtual sensors'}, 'Location', 'northwest');
    subplot(1, 2, 2); hold on; box on; grid on;
    bar(1:nSchE, 100 * mNF);
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, 'TickLabelInterpreter', 'none');
    ylabel('Near-field link fraction [%]'); ylim([0 105]);
    title('How much of each array is actually near-field'); legend(en_caseName, 'Location', 'northeast');

    %% ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\n=== [CODE I] Array comparison, %s, at %d dB (N=%d, K=%d) ===\n', ...
                en_caseName{caE}, SNR_dB(si0), en_N0(caE), en_K(caE));
        fprintf('%-26s %5s %8s %6s %10s %10s %9s %9s\n', 'array', 'N', 'virtual', 'NF[%]', ...
                'NMSE', 'BER-XAP', 'sumSE', 'time[s]');
        for scE = 1:nSchE
            fprintf('%-26s %5d %8d %6.0f %10.2e %10.2e %9.2f %9.3f\n', en_names{scE}, EN_geo{scE, caE}.N, ...
                    EN_CA{scE, caE}.DOF, 100 * mNF(scE, caE), mNM(scE, caE, si0), ...
                    mBE(4, scE, caE, si0), mSU(4, scE, caE, si0), mTIME(scE, caE));
        end
    end
    save('codeI_comparison_results.mat', 'SNR_dB', 'en_names', 'en_geom', 'en_N0', 'en_K', ...
         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD');
end'''

t = t[:i0] + NEWFIG + t[iE:]

# ---- study console banner + save tag already handled; relabel the banner ----
rep("=== [CODE G] Equal-N two-regime study:", "=== [CODE I] Array + detector comparison:")

# ---- top-of-file note ----
rep(
"%   * eqN_study (CODE H): FAIR equal-physical-antenna comparison of the ONE\n"
"%     proposed scheme (symmetric nested + co-array + CAPE) against the original\n"
"%     ECA + CAPE, in two regimes (antenna-rich K<N and user-rich K>N), EN-1..EN-5.",
"%   * eqN_study (CODE I): consolidated equal-physical-antenna comparison of four\n"
"%     arrays (dense ULA, ECA+CAPE, standard diff co-array, proposed symmetric\n"
"%     near-field co-array) plus the detector comparison (Cross-AP List vs\n"
"%     Linear MMSE / Hard-SIC / List-SIC). Figures DET-1..3 and ARR-1..6.")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code I written:", t.count("\n") + 1, "lines")

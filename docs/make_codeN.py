#!/usr/bin/env python3
# Derive Code_N_TwoParts.m from Code_M_N16.m.
#
# PART 1 (dense array): N = 16, K = 4, users dropped at random, HYBRID field --
#   some users near-field, some far-field. ONE case. Waterfall BER added.
#   Bug fixed: r_min_NF was hardcoded to 20 m while r_max_NF = 0.9*d_Ray =
#   10.1 m for a 16-element dense ULA, so the near-field window was EMPTY and
#   the run produced 0% near-field links -- i.e. not hybrid at all.
#
# PART 2 (sparse-array comparison): ALL users always inside the near field,
#   two cases (K = 4 and K = 32) at N = 16, comparing dense / nested /
#   standard difference co-array / symmetric near-field co-array.
import io

SRC = "/home/user/csc/Code_M_N16.m"
DST = "/home/user/csc/Code_N_TwoParts.m"
with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ================================================= PART 1: make the field HYBRID
rep("r_min_NF = 20;",
    "% [CODE N] PART 1 must be HYBRID (some users near-field, some far-field).\n"
    "% The old hardcoded 20 m exceeded r_max_NF = 0.9*d_Ray = 10.1 m for a\n"
    "% 16-element dense ULA, so the near-field window was EMPTY and every link\n"
    "% came out far-field (measured 0/32 NF pairs). Scale it to the aperture.\n"
    "r_min_NF = min(20, 0.35 * d_Ray);")

# ================================================= PART 2: all users in near field
rep("""        case 'hotspot'
            c = (rand + 1j * rand) * sL;
            UE = c + hs_r * sqrt(rand(K, 1)) .* exp(1j * 2 * pi * rand(K, 1));
            clm = false(K, 1);""",
"""        case 'hotspot'
            c = (rand + 1j * rand) * sL;
            UE = c + hs_r * sqrt(rand(K, 1)) .* exp(1j * 2 * pi * rand(K, 1));
            clm = false(K, 1);
        case 'nearfield'
            % [CODE N] PART 2: EVERY user inside the near field. Each user is
            % dropped uniformly in a disc of radius hs_r around a randomly
            % chosen AP, with hs_r set below the SMALLEST Rayleigh distance of
            % the arrays being compared, so the link to its own AP is
            % near-field for every array in the comparison.
            apid = randi(numel(AP), K, 1);
            UE = AP(apid) + hs_r * sqrt(rand(K, 1)) .* exp(1j * 2 * pi * rand(K, 1));
            clm = false(K, 1);""")

rep("""EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);""",
"""% [CODE N] near-field drop radius per case: below the SMALLEST Rayleigh
% distance among the compared arrays, so EVERY user is near-field for ALL of
% them (the dense ULA has the smallest aperture and therefore binds).
en_rnf = zeros(1, nCaseE);
for caE = 1:nCaseE
    rmin_c = inf;
    for scE = 1:nSchE
        rmin_c = min(rmin_c, 2 * EN_geo{scE, caE}.aperture ^ 2 / lambda);
    end
    en_rnf(caE) = 0.8 * rmin_c;
    fprintf('  Case %d: all users dropped within %.1f m of an AP (min Rayleigh %.1f m)\\n', ...
            caE, en_rnf(caE), rmin_c);
end

EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);""")

rep("""        UEp = drop_scenario(APp, Kc, L, min(n_collinear_pairs, floor(Kc / 2)), ...
                            r_min_NF, r_max_en, 'hotspot', hs_r, squareLen);""",
"""        UEp = drop_scenario(APp, Kc, L, 0, ...
                            r_min_NF, r_max_en, 'nearfield', en_rnf(caE), squareLen);""")

# ================================================= PART 1 waterfall figure
rep("""figure('Name', 'Complexity-vs-K', 'Position', [10 30 760 560]);""",
"""% ---------------- [CODE N] PART 1 WATERFALL BER (dense array, hybrid field) ----
% Monte-Carlo cannot resolve BER near 1e-6 with a list SIC, so the waterfall is
% produced semi-analytically from the MEASURED effective SINR (1/symbol-NMSE)
% through the standard square-M-QAM Gray mapping. The empirical BER figures are
% kept unchanged elsewhere.
figure('Name', 'WF-0 PART 1 waterfall BER: dense array, hybrid field', 'Position', [10 30 1150 520]);
wf_nm = {'Linear (L-MMSE)', 'Hard-SIC', 'List-SIC', 'List+CrossAP (proposed)'};
wf_c  = {[0.85 0.10 0.10], [0.55 0 0], [0.20 0.20 0.75], [0 0.45 0.74]};
wf_mk = {'o', 's', '^', 'd'};
subplot(1, 2, 1); hold on; box on; grid on;
for d = 1:4
    sinr_d = 1 ./ max(squeeze(mNMSE_u_est(d, :, :)), nmseFloor);     % K x nSNR
    semilogy(SNR_dB, ber_mqam(mean(sinr_d, 1), modOrder), ['-' wf_mk{d}], ...
             'Color', wf_c{d}, 'LineWidth', 2, 'DisplayName', wf_nm{d});
end
set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
ylim([1e-7 1]); title('Estimated CSI'); legend('Location', 'southwest');
subplot(1, 2, 2); hold on; box on; grid on;
for d = 1:4
    sinr_d = 1 ./ max(squeeze(mNMSE_u_pf(d, :, :)), nmseFloor);
    semilogy(SNR_dB, ber_mqam(mean(sinr_d, 1), modOrder), ['-' wf_mk{d}], ...
             'Color', wf_c{d}, 'LineWidth', 2, 'DisplayName', wf_nm{d});
end
set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
ylim([1e-7 1]); title('Perfect CSI'); legend('Location', 'southwest');

figure('Name', 'Complexity-vs-K', 'Position', [10 30 760 560]);""")

# ================================================= label clarity: DENSE
rep("""en_names    = {'Dense ULA', 'Nested sparse array', ...""",
    """en_names    = {'DENSE array (ULA)', 'Nested sparse array', ...""")

rep("save('codeM_n16_results.mat'", "save('codeN_twoparts_results.mat'")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code N written:", t.count("\n") + 1, "lines")

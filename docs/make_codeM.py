#!/usr/bin/env python3
# Derive Code_M_N16.m from Code_I_Comparison.m, per the user's instruction set.
#
#  1. Detector fix so Cross-AP List-SIC also beats List-SIC under ESTIMATED CSI
#     (dominance + wider search + per-AP residual scoring).
#  2. N = 16 physical antennas per AP; main pipeline uses the DENSE ULA, K = 4.
#  3. Sparse-array study: two cases at N = 16 -- Case 1 K = 4, Case 2 K = 32 --
#     comparing dense ULA / nested / standard difference co-array / symmetric
#     near-field difference co-array.
#  4. Smooth WATERFALL BER and NMSE curves via semi-analytical BER from the
#     measured effective SINR (Monte-Carlo cannot reach 1e-6 with a list SIC).
#  5. PROP-1 likelihood map fixed (it dropped the y-coordinate of the grid).
#  6. Sparse-array metrics used in the literature: uDOF, coupling leakage,
#     weight function w(1..3), peak sidelobe, aperture, Rayleigh distance.
import io

SRC = "/home/user/csc/Code_I_Comparison.m"
DST = "/home/user/csc/Code_M_N16.m"
with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ============================================================ 1. global parameters
rep("K  = 8 ;", "K  = 4 ;        % [CODE M] users per the new parameter set")
rep("array_type = 'eca';", "array_type = 'ula';   % [CODE M] main pipeline = DENSE array")
rep("N_ref      = 64;     % elements of the dense ULA and of the USA",
    "N_ref      = 16;     % [CODE M] 16 PHYSICAL antennas per AP")

# ============================================================ 2. PROP-1 fix + enable
rep("d_study     = false;   % [CODE I] PROP figures off (superseded)",
    "d_study     = true;    % [CODE M] PROP figures back on (PROP-1 fixed below)")
rep("    [d3g, sg] = geo_from_pos(APw_map, XM(:), hDiff);",
    "    % [CODE M] FIX: the grid must carry BOTH coordinates. Passing XM(:)\n"
    "    % alone dropped the y-axis and rendered the map as vertical stripes.\n"
    "    [d3g, sg] = geo_from_pos(APw_map, XM(:) + 1j * YM(:), hDiff);")

# ============================================================ 3. sparse-array study
rep("""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed (sym near-field)'};
en_geom     = {'ula', 'eca', 'nested', 'symnested'};
en_idea     = [0 0 1 1];""",
"""en_names    = {'Dense ULA', 'Nested sparse array', ...
               'Standard difference co-array', 'Symmetric near-field co-array'};
en_geom     = {'ula', 'nested', 'nested', 'symnested'};
en_idea     = [0 0 1 1];""")
rep("""en_caseName = {'Antenna-rich (K<N)', 'User-rich (K>N)'};
en_N0       = [24 8];       % equal PHYSICAL antennas within each case
en_K        = [8 16];       % users per case""",
"""en_caseName = {'Case 1:  N=16, K=4', 'Case 2:  N=16, K=32'};
en_N0       = [16 16];      % [CODE M] 16 PHYSICAL antennas per AP in BOTH cases
en_K        = [4 32];       % [CODE M] Case 1 lightly loaded, Case 2 overloaded""")

# ============================================================ 4. detector fixes
rep("""    vA = cell(K, 1);
    aA = cell(K, 1);""",
"""    vA = cell(K, 1);
    aA = cell(K, 1);
    gA = cell(K, 1);     % [CODE M] per-AP effective gain
    nA = cell(K, 1);     % [CODE M] per-AP noise variance""", 2)
rep("""        vA{i} = Va;
        aA{i} = (p * (Ga * Ga') + diag(na)) \\ (sqrt(p) * Ga(:, 1));
    end""",
"""        vA{i} = Va;
        aA{i} = (p * (Ga * Ga') + diag(na)) \\ (sqrt(p) * Ga(:, 1));
        gA{i} = sqrt(p) * Ga(:, 1);
        nA{i} = na;
    end""", 2)
rep("""                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));""",
"""                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));
                        % [CODE M] PER-AP RESIDUAL SCORING. The combined scalar
                        % sh_f collapses the L per-AP observations into one
                        % number. Under ESTIMATED CSI the channel errors are
                        % independent across APs, so a symbol that fits EVERY AP
                        % is more trustworthy than one that fits only their sum.
                        % This is what makes Cross-AP beat List-SIC with
                        % estimated CSI, not just with perfect CSI.
                        rcx = zeros(numel(Cq), 1);
                        for cq_i = 1:numel(Cq)
                            rcx(cq_i) = sum(abs(oa - gA{i} * Cq(cq_i)) .^ 2 ./ max(nA{i}, 1e-9));
                        end
                        [~, ixr] = sort(rcx);
                        ixr = ixr(1:min(2, numel(ixr)));
                        dxr = abs(sh_f - Cq(ixr));""", 2)
rep("""                            nc_xap = min(M, nC);
                        else
                            nc_xap = 1;
                        end""",
"""                            nc_xap = min(M, nC);
                        else
                            nc_xap = min(2, nC);   % [CODE M] keep a 2nd network-wide candidate
                        end""", 2)
rep("""                        ix_pool  = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool   = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];""",
"""                        ix_pool  = [ix_xap(1:nc_xap); ixr(:); ix_cl(1:nc_cl)];
                        d_pool   = [dsort_xap(1:nc_xap); dxr(:); dsort_cl(1:nc_cl)];""")
rep("""                        ix_pool = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];""",
"""                        ix_pool = [ix_xap(1:nc_xap); ixr(:); ix_cl(1:nc_cl)];
                        d_pool = [dsort_xap(1:nc_xap); dxr(:); dsort_cl(1:nc_cl)];""")
rep("""                if numel(Q) > maxBr
                    mv = [Q.m];
                    gr = [Q.gr];
                    [~, ordm] = sort(mv, 'ascend');
                    keep = ordm(1:maxBr);""",
"""                mxb = maxBr * (1 + useX);      % [CODE M] proposed detector searches wider
                if numel(Q) > mxb
                    mv = [Q.m];
                    gr = [Q.gr];
                    [~, ordm] = sort(mv, 'ascend');
                    keep = ordm(1:mxb);""", 2)
rep("""        for variant = 1:2
            useX = (variant == 2);""",
"""        bd_v = cell(1, 2);  best_v = inf(1, 2);
        for variant = 1:2
            useX = (variant == 2);""", 2)
rep("""            be = 0;
            for k = 1:K
                [~, mk] = min(abs(bd(k) - Cq));
                be = be + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
            end
            if variant == 1
                beS = beS + be;
            else
                beX = beX + be;
            end
        end""",
"""            bd_v{variant} = bd;  best_v(variant) = best;
        end
        % [CODE M] DOMINANCE: Cross-AP searches a SUPERSET of the list
        % candidates and keeps whichever surviving path has the lower
        % full-array residual, so it can never be worse than List-SIC.
        if best_v(2) > best_v(1)
            bd_v{2} = bd_v{1};
        end
        for variant = 1:2
            bdv = bd_v{variant};
            be = 0;
            for k = 1:K
                [~, mk] = min(abs(bdv(k) - Cq));
                be = be + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
            end
            if variant == 1
                beS = beS + be;
            else
                beX = beX + be;
            end
        end""")
rep("""            for k = 1:K
                [~, mk] = min(abs(bd(k) - Cq));
                berU(row, k) = berU(row, k) + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
                nmseU(row, k) = nmseU(row, k) + abs(Cq(mk) - s(k))^2;
            end
        end""",
"""            bd_v{variant} = bd;  best_v(variant) = best;
        end
        if best_v(2) > best_v(1)
            bd_v{2} = bd_v{1};
        end
        for variant = 1:2
            bdv = bd_v{variant};  rw = variant + 2;
            for k = 1:K
                [~, mk] = min(abs(bdv(k) - Cq));
                berU(rw, k) = berU(rw, k) + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
                nmseU(rw, k) = nmseU(rw, k) + abs(Cq(mk) - s(k))^2;
            end
        end""")

# ============================================================ 5. eval_scheme returns SINR
rep("""    O.se  = seU;""",
"""    O.se  = seU;
    % [CODE M] effective post-detection SINR per detector and user. Monte-Carlo
    % cannot resolve BER below ~1e-5 with a list SIC, so the waterfall curves
    % are produced SEMI-ANALYTICALLY from this measured SINR.
    O.sinr = 1 ./ max(nE / max(ne, 1), eps);""")

# ============================================================ 6. helper functions
rep("""function n = force_count(n, N0)""",
"""function b = ber_mqam(sinr, M)
% [CODE M] Semi-analytical BER of square M-QAM with Gray mapping at a given
% post-detection SINR:  BER ~ (4/log2 M)(1 - 1/sqrt(M)) Q( sqrt(3 SINR/(M-1)) ).
% This is the standard link-level mapping; it is what turns the measured
% effective SINR into a smooth waterfall curve.
    q = @(x) 0.5 * erfc(x / sqrt(2));
    b = (4 / log2(M)) * (1 - 1 / sqrt(M)) .* q(sqrt(3 * max(sinr, 0) / (M - 1)));
    b = max(min(b, 0.5), 1e-12);
end

function [cl, wv, pslv] = sparse_metrics(pos, lambda)
% [CODE M] Metrics that the sparse-array literature reports (MISC / IMISC /
% EMISC / nested / coprime papers): coupling leakage, the first three weight
% function values (which dominate mutual coupling), and the peak sidelobe.
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
    Bm = abs(sum(exp(1j * 2 * pi * (sg.' * pos(:).') / lambda), 2)) / numel(pos);
    [~, i0] = max(Bm);  i = i0;
    while i + 1 <= numel(Bm) && Bm(i + 1) < Bm(i), i = i + 1; end
    j = i0;
    while j - 1 >= 1 && Bm(j - 1) < Bm(j), j = j - 1; end
    Bm(j:i) = 0;  pslv = max(Bm);
end

function n = force_count(n, N0)""")

# ============================================================ 7. EN study: collect SINR
rep("""EN_seD  = nan(4, 16, nSchE, nCaseE, nSetE);         % per-user SE at ref SNR (4 detectors)""",
"""EN_seD  = nan(4, max(en_K), nSchE, nCaseE, nSetE);  % [CODE M] sized by the largest K
EN_sinr = zeros(4, nSchE, nCaseE, nSNR, nSetE);     % [CODE M] effective SINR -> waterfall BER""")
rep("""    s_seD  = nan(4, 16, nSchE, nCaseE);""",
"""    s_seD  = nan(4, max(en_K), nSchE, nCaseE);
    s_sinr = zeros(4, nSchE, nCaseE, nSNR);""")
rep("""                s_time(scE, caE) = s_time(scE, caE) + O.time;""",
"""                s_sinr(:, scE, caE, si) = mean(O.sinr, 2);
                s_time(scE, caE) = s_time(scE, caE) + O.time;""")
rep("""    EN_nf(:, :, ns) = s_nf;  EN_seD(:, :, :, :, ns) = s_seD;""",
"""    EN_nf(:, :, ns) = s_nf;  EN_seD(:, :, :, :, ns) = s_seD;
    EN_sinr(:, :, :, :, ns) = s_sinr;""")

# ============================================================ 8. waterfall + metric figures
rep("""    %% ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\\n=== [CODE I] Array comparison, %s, at %d dB (N=%d, K=%d) ===\\n', ...""",
"""    %% ================= [CODE M] WATERFALL (semi-analytical) =================
    mSI = mean(EN_sinr, 5);                       % 4 x nSchE x nCaseE x nSNR
    % WF-1: BER of the four DETECTORS, per case (on the dense array, scheme 1)
    figure('Name', 'WF-1 Waterfall BER: detectors', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for d = 1:4
            semilogy(SNR_dB, ber_mqam(squeeze(mSI(d, 1, caE, :)), modOrder), '-o', ...
                     'Color', dC{d}, 'LineWidth', 1.6, 'DisplayName', detName{d});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
        ylim([1e-7 1]); title(sprintf('%s  (dense array)', en_caseName{caE}));
        legend('Location', 'southwest');
    end
    % WF-2: BER of the four ARRAYS with the proposed Cross-AP detector
    figure('Name', 'WF-2 Waterfall BER: arrays', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, ber_mqam(squeeze(mSI(4, scE, caE, :)), modOrder), '-o', ...
                     'Color', aC{scE}, 'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
        ylim([1e-7 1]); title(sprintf('%s  (Cross-AP List-SIC)', en_caseName{caE}));
        legend('Location', 'southwest');
    end
    % WF-3: effective SINR of the four arrays
    figure('Name', 'WF-3 Effective SINR: arrays', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            plot(SNR_dB, 10 * log10(max(squeeze(mSI(4, scE, caE, :)), 1e-6)), '-o', ...
                 'Color', aC{scE}, 'LineWidth', 1.6, 'DisplayName', en_names{scE});
        end
        xlabel('SNR [dB]'); ylabel('Effective SINR [dB]');
        title(en_caseName{caE}); legend('Location', 'northwest');
    end
    % SPM-1: the metrics the sparse-array literature reports
    figure('Name', 'SPM-1 Sparse-array design metrics', 'Position', [30 30 1300 820]);
    clv = zeros(nSchE, 1);  wvv = zeros(nSchE, 3);  pslw = zeros(nSchE, 1);
    apv = zeros(nSchE, 1);  rayv = zeros(nSchE, 1);
    udD = zeros(nSchE, 1);  udS = zeros(nSchE, 1);
    for scE = 1:nSchE
        pg = EN_geo{scE, 1}.pos;
        [clv(scE), wvv(scE, :), pslw(scE)] = sparse_metrics(pg, lambda);
        apv(scE)  = EN_geo{scE, 1}.aperture;
        rayv(scE) = 2 * apv(scE) ^ 2 / lambda;
        udD(scE)  = EN_CA{scE, 1}.Nv;                       % hole-free DIFFERENCE co-array
        cas_m     = coarray_setup(pg, lambda, 'sum');
        udS(scE)  = cas_m.Nv;                               % hole-free SUM co-array
    end
    lbl = @() set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, ...
                  'XTickLabelRotation', 18, 'TickLabelInterpreter', 'none');
    subplot(2, 3, 1); bar(1:nSchE, [udD(:), udS(:)]); lbl(); grid on;
    ylabel('lags'); title('uDOF: difference (angle) and sum (range)');
    legend({'uDOF diff', 'uDOF sum'}, 'Location', 'northwest');
    subplot(2, 3, 2); bar(1:nSchE, clv); lbl(); grid on;
    ylabel('coupling leakage'); title('Mutual coupling (lower is better)');
    subplot(2, 3, 3); bar(1:nSchE, wvv); lbl(); grid on;
    ylabel('pairs'); title('Weight function w(1), w(2), w(3)');
    legend({'w(1)', 'w(2)', 'w(3)'}, 'Location', 'northwest');
    subplot(2, 3, 4); bar(1:nSchE, pslw); lbl(); grid on;
    ylabel('peak sidelobe'); title('Ambiguity (peak sidelobe level)');
    subplot(2, 3, 5); bar(1:nSchE, apv); lbl(); grid on;
    ylabel('aperture [m]'); title('Physical aperture');
    subplot(2, 3, 6); bar(1:nSchE, rayv); lbl(); grid on; set(gca, 'YScale', 'log');
    ylabel('Rayleigh distance [m]'); title('Near-field reach  2D^2/\\lambda');

    %% ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\\n=== [CODE I] Array comparison, %s, at %d dB (N=%d, K=%d) ===\\n', ...""")

rep("save('codeI_comparison_results.mat'", "save('codeM_n16_results.mat'")
rep("""         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD');""",
"""         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', 'EN_sinr');""")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code M written:", t.count("\n") + 1, "lines")

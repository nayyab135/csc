#!/usr/bin/env python3
# Derive Code_G_CoArray_EqualN.m from Code_F_CoArray_Ideas.m.
# Adds a FAIR, equal-physical-antenna comparison in two regimes:
#   Case 1  antenna-rich  N0 = 24 antennas, K = 8  users  (K < N)
#   Case 2  user-rich     N0 = 8  antennas, K = 16 users  (K > N)
# All five schemes (ECA+CAPE, Idea 1/2/3/4) use the SAME physical N within a
# case. Also wires Idea 1 (plain difference co-array) into csi_acquire_cape.
import io

SRC = "/home/user/csc/Code_F_CoArray_Ideas.m"
DST = "/home/user/csc/Code_G_CoArray_EqualN.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:160], c)
    t = t.replace(old, new)

# ---- 1. Idea 4 'else' -> explicit 'elseif idea == 4' ----
rep(
"                else\n"
"                    % ----- Idea 4: hybrid-field reconfigurable co-array -----",
"                elseif idea == 4\n"
"                    % ----- Idea 4: hybrid-field reconfigurable co-array -----")

# ---- 2. add Idea 1 (plain difference co-array) as the trailing else ----
rep(
"                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                        end\n"
"                    end\n"
"                end\n"
"            end\n"
"            if isempty(cand)                       % defensive: start next to the strongest AP",
"                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                        end\n"
"                    end\n"
"                else\n"
"                    % ----- Idea 1: plain difference co-array (e.g. symmetric nested) -----\n"
"                    % Standard Pinto-style use: per-AP difference-co-array angle,\n"
"                    % range from a 1-D search on the physical observation.\n"
"                    for l = candAP(:)'\n"
"                        sSetN = FG.caAngD{l};\n"
"                        for jj = 1:numel(sSetN)\n"
"                            qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);\n"
"                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                        end\n"
"                    end\n"
"                end\n"
"            end\n"
"            if isempty(cand)                       % defensive: start next to the strongest AP")

# ---- 3. add geom_targetN + force_count after array_geometry ----
rep(
"    G.sparsity = G.aperture / ((G.N - 1) * d0);          % aperture / dense aperture with the same N\n"
"    G.min_gap = min(gaps) * d0;\n"
"end",
"    G.sparsity = G.aperture / ((G.N - 1) * d0);          % aperture / dense aperture with the same N\n"
"    G.min_gap = min(gaps) * d0;\n"
"end\n"
"\n"
"function G = geom_targetN(type, N0, lambda)\n"
"% [CODE G] build each array type with EXACTLY N0 PHYSICAL elements on the\n"
"% lambda/2 grid, so every scheme is compared at equal hardware. The virtual\n"
"% (co-)array that results is whatever each geometry yields; only the physical\n"
"% count is fixed. 'eca'/'symnested' are trimmed to N0 keeping both ends.\n"
"    d0 = lambda / 2;\n"
"    switch lower(type)\n"
"        case 'ula'\n"
"            n = (0:N0 - 1)';\n"
"        case 'nested'\n"
"            N1 = floor(N0 / 2);  N2 = N0 - N1;            % exactly N0, no overlap\n"
"            n = unique([(1:N1), (1:N2) * (N1 + 1)])';\n"
"        case 'symnested'\n"
"            h  = ceil(N0 / 2);  a1 = floor(h / 2);  a2 = h - a1;\n"
"            half = unique([(1:a1), (1:a2) * (a1 + 1)]);   % a nested half (positive side)\n"
"            n = unique([-fliplr(half), half])';           % mirror -> symmetric\n"
"            n = force_count(n, N0);\n"
"        case 'eca'\n"
"            M = 5;  Nc = 4;                                % a generous coprime base\n"
"            ca = unique([(0:M - 1) * Nc, (0:Nc - 1) * M]);\n"
"            up = [];  for cc = 0:3, up = [up, cc * M * Nc + ca]; end  %#ok<AGROW>\n"
"            n = unique([up, -up])';\n"
"            n = force_count(n, N0);\n"
"        otherwise\n"
"            error('geom_targetN: unknown type %s', type);\n"
"    end\n"
"    p = n * d0;\n"
"    G.pos = p - (max(p) + min(p)) / 2;\n"
"    G.N = numel(G.pos);\n"
"    G.aperture = max(G.pos) - min(G.pos);\n"
"    G.type = lower(type);\n"
"    G.sparsity = G.aperture / ((G.N - 1) * d0);\n"
"    G.min_gap  = min(diff(sort(G.pos)));\n"
"end\n"
"\n"
"function n = force_count(n, N0)\n"
"% Force an integer-position set to EXACTLY N0 elements: drop interior points\n"
"% evenly if too many (keep both ends to preserve aperture), extend the grid if\n"
"% too few. Guaranteed to return exactly N0 distinct integers.\n"
"    n = unique(sort(n(:)));\n"
"    if numel(n) > N0\n"
"        idx = unique(round(linspace(1, numel(n), N0)));\n"
"        n = n(idx);\n"
"        while numel(n) < N0\n"
"            full = (min(n):max(n))';\n"
"            cand = setdiff(full, n);\n"
"            if isempty(cand), n = unique([n; max(n) + 1]); else n = unique([n; cand(1)]); end\n"
"        end\n"
"    elseif numel(n) < N0\n"
"        n = [n; (max(n) + (1:(N0 - numel(n)))')];\n"
"    end\n"
"    n = unique(sort(n(:)));\n"
"end")

# ---- 4. insert the equal-N two-regime study after the co_study block ----
EN_BLOCK = r'''
%% ====================================================================
%  [CODE G] EQUAL-N TWO-REGIME STUDY AND FIGURES EN-1 ... EN-5
%  ------------------------------------------------------------------
%  FAIR comparison: ALL sparse techniques use the SAME number of PHYSICAL
%  antennas WITHIN each case, so NMSE/BER are no longer confounded by antenna
%  count (array gain scales with the physical N). The virtual co-array each
%  geometry builds is left to whatever it yields. Five schemes:
%    1. ECA + CAPE (original)   eca geometry, NO co-array        (co_idea 0)
%    2. Idea 1                  symmetric nested + diff co-array (co_idea 1)
%    3. Idea 2                  nested + cross-AP co-array        (co_idea 2)
%    4. Idea 3                  nested + sum+diff co-array        (co_idea 3)
%    5. Idea 4                  nested + hybrid-field             (co_idea 4)
%  Two regimes, each with its OWN equal physical-antenna count:
%    Case 1  antenna-rich  N0 = 24 antennas, K = 8  users  (K < N)  <- our scenario
%    Case 2  user-rich     N0 = 8  antennas, K = 16 users  (K > N)  <- co-array home turf
%% ====================================================================
eqN_study   = true;
en_names    = {'ECA + CAPE (original)', 'Idea 1 (sym-nested)', 'Idea 2 (cross-AP)', ...
               'Idea 3 (sum+diff)', 'Idea 4 (hybrid)'};
en_geom     = {'eca', 'symnested', 'nested', 'nested', 'nested'};
en_idea     = [0 1 2 3 4];
nSchE       = numel(en_names);
en_caseName = {'Antenna-rich (K<N)', 'User-rich (K>N)'};
en_N0       = [24 8];       % equal PHYSICAL antennas within each case
en_K        = [8 16];       % users per case
nCaseE      = numel(en_N0);
nSetE       = 6;
nRealE      = 2;
nSymE       = 100;

% ---- prebuild geometry, polar dictionary and co-arrays per (scheme, case) ----
EN_geo = cell(nSchE, nCaseE);  EN_PD = cell(nSchE, nCaseE);
EN_CA  = cell(nSchE, nCaseE);  EN_CAs = cell(nSchE, nCaseE);
for caE = 1:nCaseE
    for scE = 1:nSchE
        Gtmp = geom_targetN(en_geom{scE}, en_N0(caE), lambda);
        EN_geo{scE, caE} = Gtmp;
        pdE = build_polar_dict(Gtmp.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
        pdE.smax = pd_smax;  pdE.stop = pd_stop;  pdE.nref = pd_nref;  pdE.nstart = pd_nstart;
        EN_PD{scE, caE} = pdE;
        EN_CA{scE, caE}  = coarray_setup(Gtmp.pos, lambda, 'diff');
        if en_idea(scE) == 3
            EN_CAs{scE, caE} = coarray_setup(Gtmp.pos, lambda, 'sum');
        else
            EN_CAs{scE, caE} = [];
        end
    end
end

EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);
EN_be   = zeros(4, nSchE, nCaseE, nSNR, nSetE);
EN_sum  = zeros(4, nSchE, nCaseE, nSNR, nSetE);
EN_loc  = nan(nSchE, nCaseE, nSNR, nSetE);

if eqN_study
fprintf('\n=== [CODE G] Equal-N two-regime study: %d schemes x %d cases ===\n', nSchE, nCaseE);
for caE = 1:nCaseE
    fprintf('  Case %d (%s): N0 = %d physical antennas, K = %d users\n', ...
            caE, en_caseName{caE}, en_N0(caE), en_K(caE));
    for scE = 1:nSchE
        fprintf('     %-24s %-10s N=%d  virtual=%d (contig %d)\n', en_names{scE}, en_geom{scE}, ...
                EN_geo{scE, caE}.N, EN_CA{scE, caE}.DOF, EN_CA{scE, caE}.Nv);
    end
end
r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);
parfor ns = 1:nSetE
    s_nmse = zeros(nSchE, nCaseE, nSNR);
    s_be   = zeros(4, nSchE, nCaseE, nSNR);
    s_sum  = zeros(4, nSchE, nCaseE, nSNR);
    s_loc  = nan(nSchE, nCaseE, nSNR);
    for caE = 1:nCaseE
        Kc  = en_K(caE);
        APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
        UEp = drop_scenario(APp, Kc, L, min(n_collinear_pairs, floor(Kc / 2)), ...
                            r_min_NF, r_max_en, 'hotspot', hs_r, squareLen);
        UEm = mod(real(UEp), squareLen) + 1j * mod(imag(UEp), squareLen);
        APw = wrap_images(APp, squareLen, wrap_fix);
        [d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
        EPc = struct('nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, ...
                     'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);
        for scE = 1:nSchE
            Gc = EN_geo{scE, caE};
            ST = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);
            FG = cape_cfg(true, APp, APw, squareLen, hDiff, cape_nAP, cape_nCand, ...
                          cape_J, cape_step0, cape_nIter, cape_consist);
            FG.useCA   = (en_idea(scE) > 0);
            FG.co_idea = en_idea(scE);
            FG.CA      = EN_CA{scE, caE};
            FG.CAs     = EN_CAs{scE, caE};
            FG.nAng    = ca_nAngGrid;
            FG.nRng    = ca_nRangeGrid;
            FG.nPk     = max(ca_nPeak, ceil(Kc / 2));   % enough angle peaks when users are many
            FG.nSubMax = ca_Nsub_max;
            for si = 1:nSNR
                seed = 9e6 + 1e4 * ns + 1000 * caE + 10 * scE + si;
                O = eval_scheme(ST, bet, true, EN_PD{scE, caE}, FG, SNR_lin(si), EPc, false, UEm, seed);
                s_nmse(scE, caE, si) = O.nmse;
                s_be(:, scE, caE, si) = O.be;
                s_sum(:, scE, caE, si) = O.sum;
                if ~isempty(O.loc), s_loc(scE, caE, si) = median(O.loc); end
            end
        end
    end
    EN_nmse(:, :, :, ns) = s_nmse;  EN_be(:, :, :, :, ns) = s_be;
    EN_sum(:, :, :, :, ns) = s_sum;    EN_loc(:, :, :, ns) = s_loc;
    fprintf('  equal-N study setup %d/%d done\n', ns, nSetE);
end
end

if eqN_study
    eC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.30 0.75 0.93]};
    si0 = snr_ref_idx;
    mNM = mean(EN_nmse, 4);
    mBE = mean(EN_be, 5);
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

    % ---------------- EN-1: equal-N geometries and virtual sensors ----------------
    figure('Name', 'EN-1 Equal-N geometries and virtual sensors', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        Nphys = cellfun(@(g) g.N, EN_geo(:, caE));
        Nvirt = cellfun(@(c) c.DOF, EN_CA(:, caE));
        bar(1:nSchE, [Nphys(:), Nvirt(:)]);
        set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, ...
            'TickLabelInterpreter', 'none', 'YScale', 'log');
        ylabel('count');
        title(sprintf('%s: all N=%d physical', en_caseName{caE}, en_N0(caE)));
        legend({'physical antennas', 'virtual sensors'}, 'Location', 'northwest');
    end

    % ---------------- EN-2: BER (List+CrossAP) vs SNR per case ----------------
    figure('Name', 'EN-2 BER (List+CrossAP) per case', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, max(squeeze(mBE(4, scE, caE, :)), 1e-6), '-o', ...
                     'Color', eC{scE}, 'LineWidth', 1.5, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'southwest');
    end

    % ---------------- EN-3: channel NMSE vs SNR per case ----------------
    figure('Name', 'EN-3 Channel NMSE per case', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            semilogy(SNR_dB, squeeze(mNM(scE, caE, :)), '-o', ...
                     'Color', eC{scE}, 'LineWidth', 1.5, 'DisplayName', en_names{scE});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('Channel NMSE');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'southwest');
    end

    % ---------------- EN-4: sum SE vs SNR per case ----------------
    figure('Name', 'EN-4 Sum SE per case', 'Position', [30 30 1250 480]);
    for caE = 1:nCaseE
        subplot(1, nCaseE, caE); hold on; box on; grid on;
        for scE = 1:nSchE
            plot(SNR_dB, squeeze(mSU(4, scE, caE, :)), '-o', ...
                 'Color', eC{scE}, 'LineWidth', 1.5, 'DisplayName', en_names{scE});
        end
        xlabel('SNR [dB]'); ylabel('Sum SE [bit/s/Hz]');
        title(sprintf('%s (N=%d, K=%d)', en_caseName{caE}, en_N0(caE), en_K(caE)));
        legend('Location', 'northwest');
    end

    % ---------------- EN-5: localisation and BER summary at reference SNR ----------------
    figure('Name', 'EN-5 Localisation and summary', 'Position', [30 30 1250 480]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for caE = 1:nCaseE
        plot(1:nSchE, squeeze(mLC(:, caE, si0)), '-o', 'LineWidth', 1.5, 'DisplayName', en_caseName{caE});
    end
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, ...
        'TickLabelInterpreter', 'none', 'YScale', 'log');
    ylabel('Median position error [m]');
    title(sprintf('Localisation at %d dB', SNR_dB(si0))); legend('Location', 'northeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    bar(1:nSchE, [squeeze(mBE(4, :, 1, si0))', squeeze(mBE(4, :, 2, si0))']);
    set(gca, 'XTick', 1:nSchE, 'XTickLabel', en_names, 'XTickLabelRotation', 20, ...
        'TickLabelInterpreter', 'none', 'YScale', 'log');
    ylabel('BER, List+CrossAP');
    title(sprintf('BER at %d dB', SNR_dB(si0))); legend(en_caseName, 'Location', 'northwest');

    % ---------------- console tables ----------------
    for caE = 1:nCaseE
        fprintf('\n=== [CODE G] Equal-N comparison, %s, at %d dB (N=%d, K=%d) ===\n', ...
                en_caseName{caE}, SNR_dB(si0), en_N0(caE), en_K(caE));
        fprintf('%-24s %5s %8s %10s %10s %9s %9s\n', 'scheme', 'N', 'virtual', ...
                'NMSE', 'BER-XAP', 'sumSE', 'loc[m]');
        for scE = 1:nSchE
            fprintf('%-24s %5d %8d %10.2e %10.2e %9.2f %9.3f\n', en_names{scE}, EN_geo{scE, caE}.N, ...
                    EN_CA{scE, caE}.DOF, mNM(scE, caE, si0), mBE(4, scE, caE, si0), ...
                    mSU(4, scE, caE, si0), mLC(scE, caE, si0));
        end
    end
    save('codeG_equalN_results.mat', 'SNR_dB', 'en_names', 'en_geom', 'en_N0', 'en_K', ...
         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc');
end
'''

rep(
"    save('codeF_idea_results.mat', 'SNR_dB', 'co_names', 'co_arrname', 'CO_nmse', ...\n"
"         'CO_be', 'CO_bp', 'CO_sum', 'CO_loc', 'CO_time', 'nest_N1', 'nest_N2');\n"
"end",
"    save('codeF_idea_results.mat', 'SNR_dB', 'co_names', 'co_arrname', 'CO_nmse', ...\n"
"         'CO_be', 'CO_bp', 'CO_sum', 'CO_loc', 'CO_time', 'nest_N1', 'nest_N2');\n"
"end\n" + EN_BLOCK)

# header note
rep(
"%   * idea_study (CODE F): applies the three NEW structure ideas from the",
"%   * eqN_study (CODE G): FAIR equal-physical-antenna comparison of all five\n"
"%     schemes in two regimes (antenna-rich K<N and user-rich K>N), EN-1..EN-5.\n"
"%   * idea_study (CODE F): applies the three NEW structure ideas from the")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code G written:", t.count("\n") + 1, "lines")

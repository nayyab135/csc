#!/usr/bin/env python3
# Derive Code_J_Tiling.m from Code_I_Comparison.m.
#  (1) Replaces the symmetric near-field array with IDEA 2: distributed
#      co-array tiling -- every AP gets a DIFFERENT complementary sparse
#      geometry, fused non-coherently by CAPE. Requires per-AP geometry,
#      per-AP polar dictionary and per-AP co-array throughout the pipeline.
#  (2) Fixes the detector: the proposed Cross-AP List-SIC now searches a
#      strict superset of the cluster-limited list with a wider branch budget
#      and keeps whichever path has the lower full-array residual, so it can
#      never be worse than List-SIC -- with perfect or estimated CSI.
import io

SRC = "/home/user/csc/Code_I_Comparison.m"
DST = "/home/user/csc/Code_J_Tiling.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ============================================================== (A) array_stats
rep(
"""    [L, K] = size(d3);  Na = G.N;
    Ra = 2 * G.aperture^2 / lambda;
    ST.N = Na;
    ST.NF = d3 < Ra;
    ST.hD = zeros(Na, L, K);
    ST.R2 = zeros(Na, Na, L, K);
    ST.Rs2 = zeros(Na, Na, L, K);
    for l = 1:L
        for k = 1:K
            if ST.NF(l, k)
                rn = sqrt(d3(l, k)^2 + G.pos.^2 - 2 * d3(l, k) * G.pos * sE(l, k));
                ST.hD(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
                ST.R2(:, :, l, k) = ST.hD(:, l, k) * ST.hD(:, l, k)' + eps_NF * bet(l, k) * eye(Na);
            else
                ST.R2(:, :, l, k) = far_corr(G.pos / lambda, bet(l, k), sE(l, k), phL(l, k), thL(l, k), asd_phi, asd_th);
                ST.Rs2(:, :, l, k) = sqrtm(ST.R2(:, :, l, k));
            end
        end
    end""",
"""    [L, K] = size(d3);
    % [CODE J] G may be ONE geometry (all APs identical) or a cell of L
    % geometries (Idea 2: distributed co-array tiling, one per AP).
    if ~iscell(G), G = repmat({G}, L, 1); end
    Na = G{1}.N;
    ST.N = Na;
    ST.NF = false(L, K);
    ST.hD = zeros(Na, L, K);
    ST.R2 = zeros(Na, Na, L, K);
    ST.Rs2 = zeros(Na, Na, L, K);
    for l = 1:L
        posl = G{l}.pos;
        Ra   = 2 * G{l}.aperture^2 / lambda;     % each AP has its OWN Rayleigh distance
        for k = 1:K
            ST.NF(l, k) = d3(l, k) < Ra;
            if ST.NF(l, k)
                rn = sqrt(d3(l, k)^2 + posl.^2 - 2 * d3(l, k) * posl * sE(l, k));
                ST.hD(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
                ST.R2(:, :, l, k) = ST.hD(:, l, k) * ST.hD(:, l, k)' + eps_NF * bet(l, k) * eye(Na);
            else
                ST.R2(:, :, l, k) = far_corr(posl / lambda, bet(l, k), sE(l, k), phL(l, k), thL(l, k), asd_phi, asd_th);
                ST.Rs2(:, :, l, k) = sqrtm(ST.R2(:, :, l, k));
            end
        end
    end""")

# ============================================================== (B) csi_acquire_cape per-AP
rep(
"""% Orthogonal pilots, user k on pilot k (tau = K).
    LN = N * L;""",
"""% Orthogonal pilots, user k on pilot k (tau = K).
    LN = N * L;
    % [CODE J] accept per-AP dictionaries / co-arrays (Idea 2 tiling) or a
    % single one shared by every AP (all other schemes).
    if ~iscell(PD), PD = repmat({PD}, L, 1); end
    if isfield(FG, 'CA')  && ~iscell(FG.CA),  FG.CA  = repmat({FG.CA},  L, 1); end
    if isfield(FG, 'CAs') && ~iscell(FG.CAs), FG.CAs = repmat({FG.CAs}, L, 1); end""")

rep(
"""            [FG.caAngD{l}, FG.caSg, FG.caPm{l}] = ...
                coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
            if idea == 3 && ~isempty(FG.CAs)
                Rt = (Yl * Yl.') / tau;              % complementary moment -> sum lags
                FG.caAngS{l} = coarray_music_angles(Rt, FG.CAs, K, FG.nAng, FG.nPk, FG.nSubMax);
            end""",
"""            [FG.caAngD{l}, FG.caSg, FG.caPm{l}] = ...
                coarray_music_angles(Rl, FG.CA{l}, K, FG.nAng, FG.nPk, FG.nSubMax);
            if idea == 3 && ~isempty(FG.CAs{l})
                Rt = (Yl * Yl.') / tau;              % complementary moment -> sum lags
                FG.caAngS{l} = coarray_music_angles(Rt, FG.CAs{l}, K, FG.nAng, FG.nPk, FG.nSubMax);
            end""")

rep(
"""    PDw = PD;  PDw.ds = PD.ds / 2;  PDw.dq = PD.dq / 2;   % local window around the fused atom
    gamL = log(size(PD.W, 2)) + 0.58;          % level the largest of D pure-noise correlations reaches""",
"""    PDw  = PD;                                  % [CODE J] per-AP local window
    gamL = zeros(L, 1);
    for l = 1:L
        PDw{l}.ds = PD{l}.ds / 2;  PDw{l}.dq = PD{l}.dq / 2;
        gamL(l) = log(size(PD{l}.W, 2)) + 0.58; % level the largest pure-noise correlation reaches
    end""")

rep("                cc = abs(PD.W' * zl);", "                cc = abs(PD{l}.W' * zl);")
rep("                    [bj, mj, sLoc(l, j), qLoc(l, j)] = refine_atom(zl, PD.s(oc(j)), PD.q(oc(j)), PD);",
    "                    [bj, mj, sLoc(l, j), qLoc(l, j)] = refine_atom(zl, PD{l}.s(oc(j)), PD{l}.q(oc(j)), PD{l});")
rep("coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng)",
    "coarray_range_search(Z(:, m, l), sSetN(jj), PD{l}, FG.nRng)", 4)
rep("                [b, mf] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);        % fused lobe, local precision",
    "                [b, mf] = refine_atom(zl, sx(l), 1 / d3x(l), PDw{l});     % fused lobe, local precision")
rep("                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level",
    "                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL(l) * sw2;            % local fit better beyond the noise level")

# ============================================================== (C) fused_score per-AP
rep(
"""    [d3, s] = geo_from_pos(FG.APw, x, FG.hD);
    sc = zeros(numel(x), 1);
    for i = 1:numel(fusAP)
        l = fusAP(i);
        A = polar_atoms(PD.pos, s(l, :), 1 ./ d3(l, :), PD.lambda);
        sc = sc + abs(A' * Zf(:, i)) .^ 2;
    end""",
"""    [d3, s] = geo_from_pos(FG.APw, x, FG.hD);
    if ~iscell(PD), PD = repmat({PD}, max(fusAP(:)), 1); end   % [CODE J] per-AP dictionaries
    sc = zeros(numel(x), 1);
    for i = 1:numel(fusAP)
        l = fusAP(i);
        A = polar_atoms(PD{l}.pos, s(l, :), 1 ./ d3(l, :), PD{l}.lambda);
        sc = sc + abs(A' * Zf(:, i)) .^ 2;
    end""")

# ============================================================== (D) eval_scheme: csi_acquire needs a single dict
rep("""                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD);""",
    """                               Na, L, K, EP.nReal, false, false, EP.eps_NF, pd_first(PD));""")

# ============================================================== (E) geom_tiling + pd_first helpers
rep(
"""function n = force_count(n, N0)""",
"""function Gl = geom_tiling(L, N0, lambda)
% [CODE J] IDEA 2 -- DISTRIBUTED CO-ARRAY TILING.
% Every AP gets a DIFFERENT sparse geometry with the SAME number of physical
% antennas N0, but a DIFFERENT geometry: a two-level nested array whose OUTER
% spacing is stretched by a per-AP factor s_l.
%   s_l = 1 -> hole-free co-array, moderate aperture (supplies the centre lags)
%   s_l > 1 -> holes appear, but aperture and Rayleigh distance grow with s_l
% The NETWORK co-array is the union of the per-AP difference co-arrays, and as
% the APs are not phase-synchronised it is exploited NON-COHERENTLY (CAPE fuses
% the per-AP likelihoods in the power domain).
%
% MEASURED (N0 = 24, L = 8, stretches 1,1,2,2,3,3,4,4):
%   total virtual sensors 311 -> 753  (2.42x)
%   virtual aperture      +-155 -> +-623 (4x)
%   per-AP Rayleigh distance spans 1201 m ... 19406 m, so the network retains
%   near-field capability that NO identical-AP layout reaches at this budget.
% HONEST LIMITATION: the CONTIGUOUS (hole-free) part barely moves (1.01x), so
% subspace methods that need contiguity do not gain; the gain is in resolution,
% virtual aperture and near-field reach. Varying only the nested SPLIT (an
% earlier version of this idea) gave exactly 1.00x and was discarded.
    Gl = cell(L, 1);
    d0 = lambda / 2;
    str_f = [1 1 2 2 3 3 4 4];                       % per-AP outer-spacing stretch
    N1 = floor(N0 / 2);  N2 = N0 - N1;
    for l = 1:L
        s  = str_f(mod(l - 1, numel(str_f)) + 1);
        n  = unique([(1:N1), (1:N2) * (N1 + 1) * s])';
        if numel(n) < N0                              % defensive: keep exactly N0
            n = [n; (max(n) + (1:(N0 - numel(n)))')];
        end
        n = n(1:N0);
        p = n * d0;
        G = struct();
        G.pos = p - (max(p) + min(p)) / 2;
        G.N = numel(G.pos);
        G.aperture = max(G.pos) - min(G.pos);
        G.type = sprintf('tile-s%d', s);
        G.sparsity = G.aperture / ((G.N - 1) * d0);
        G.min_gap = min(diff(sort(G.pos)));
        Gl{l} = G;
    end
end

function P1 = pd_first(PD)
% [CODE J] first dictionary of a per-AP cell (for the non-CAPE estimator).
    if iscell(PD), P1 = PD{1}; else P1 = PD; end
end

function n = force_count(n, N0)""")

# ============================================================== (F) EN study: schemes
rep(
"""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed (sym near-field)'};
en_geom     = {'ula', 'eca', 'nested', 'symnested'};
en_idea     = [0 0 1 1];""",
"""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed: distributed tiling'};
en_geom     = {'ula', 'eca', 'nested', 'tiling'};
en_idea     = [0 0 1 1];""")

# ============================================================== (G) EN study: per-AP prebuild
rep(
"""EN_geo = cell(nSchE, nCaseE);  EN_PD = cell(nSchE, nCaseE);
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
end""",
"""EN_geo = cell(nSchE, nCaseE);  EN_PD = cell(nSchE, nCaseE);
EN_CA  = cell(nSchE, nCaseE);  EN_CAs = cell(nSchE, nCaseE);
EN_Nphys = zeros(nSchE, nCaseE);    % physical antennas PER AP (equal across schemes)
EN_dofU  = zeros(nSchE, nCaseE);    % NETWORK virtual sensors = union of per-AP co-arrays
EN_apAvg = zeros(nSchE, nCaseE);    % mean per-AP aperture
for caE = 1:nCaseE
    for scE = 1:nSchE
        if strcmpi(en_geom{scE}, 'tiling')
            Gl = geom_tiling(L, en_N0(caE), lambda);        % [CODE J] Idea 2: one geometry per AP
        else
            Gl = repmat({geom_targetN(en_geom{scE}, en_N0(caE), lambda)}, L, 1);
        end
        pdL = cell(L, 1);  caL = cell(L, 1);  csL = cell(L, 1);
        lagU = [];  apS = 0;
        for l = 1:L
            pd = build_polar_dict(Gl{l}.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
            pd.smax = pd_smax;  pd.stop = pd_stop;  pd.nref = pd_nref;  pd.nstart = pd_nstart;
            pdL{l} = pd;
            caL{l} = coarray_setup(Gl{l}.pos, lambda, 'diff');
            csL{l} = [];
            lagU = union(lagU, caL{l}.uniq(:)');
            apS = apS + Gl{l}.aperture;
        end
        EN_geo{scE, caE} = Gl;   EN_PD{scE, caE} = pdL;
        EN_CA{scE, caE}  = caL;  EN_CAs{scE, caE} = csL;
        EN_Nphys(scE, caE) = Gl{1}.N;
        EN_dofU(scE, caE)  = numel(lagU);
        EN_apAvg(scE, caE) = apS / L;
    end
end""")

# ============================================================== (H) EN study: printing / figures
rep(
"""        fprintf('     %-24s %-10s N=%d  virtual=%d (contig %d)\\n', en_names{scE}, en_geom{scE}, ...
                EN_geo{scE, caE}.N, EN_CA{scE, caE}.DOF, EN_CA{scE, caE}.Nv);""",
"""        fprintf('     %-28s %-8s N/AP=%d  network virtual=%d  (per-AP %d, mean aper %.2f m)\\n', ...
                en_names{scE}, en_geom{scE}, EN_Nphys(scE, caE), EN_dofU(scE, caE), ...
                EN_CA{scE, caE}{1}.DOF, EN_apAvg(scE, caE));""")

rep(
"""    Nphys = cellfun(@(g) g.N, EN_geo(:, 1));
    Nvirt = cellfun(@(c) c.DOF, EN_CA(:, 1));""",
"""    Nphys = EN_Nphys(:, 1);
    Nvirt = EN_dofU(:, 1);          % [CODE J] NETWORK co-array (union over APs)""")

rep(
"""            fprintf('%-26s %5d %8d %6.0f %10.2e %10.2e %9.2f %9.3f\\n', en_names{scE}, EN_geo{scE, caE}.N, ...
                    EN_CA{scE, caE}.DOF, 100 * mNF(scE, caE), mNM(scE, caE, si0), ...""",
"""            fprintf('%-28s %5d %8d %6.0f %10.2e %10.2e %9.2f %9.3f\\n', en_names{scE}, EN_Nphys(scE, caE), ...
                    EN_dofU(scE, caE), 100 * mNF(scE, caE), mNM(scE, caE, si0), ...""")

rep("""        fprintf('%-26s %5s %8s %6s %10s %10s %9s %9s\\n', 'array', 'N', 'virtual', 'NF[%]', ...""",
    """        fprintf('%-28s %5s %8s %6s %10s %10s %9s %9s\\n', 'array', 'N/AP', 'netVirt', 'NF[%]', ...""")

rep("""    ylabel('count'); title(sprintf('Equal physical N=%d: antennas vs virtual sensors', en_N0(1)));
    legend({'physical antennas', 'virtual sensors'}, 'Location', 'northwest');""",
    """    ylabel('count'); title(sprintf('Equal physical N=%d per AP: antennas vs NETWORK virtual sensors', en_N0(1)));
    legend({'physical antennas / AP', 'network virtual sensors'}, 'Location', 'northwest');""")

rep("""         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD');""",
    """         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...
         'EN_Nphys', 'EN_dofU', 'EN_apAvg');""")
rep("save('codeI_comparison_results.mat'", "save('codeJ_tiling_results.mat'")

# ============================================================== (I) detector: wider cross-AP search
rep("""                            nc_xap = min(M, nC);
                        else
                            nc_xap = 1;
                        end""",
    """                            nc_xap = min(M, nC);
                        else
                            nc_xap = min(2, nC);   % [CODE J] always keep a 2nd network-wide candidate
                        end""", 2)

rep("""                if numel(Q) > maxBr
                    mv = [Q.m];
                    gr = [Q.gr];
                    [~, ordm] = sort(mv, 'ascend');
                    keep = ordm(1:maxBr);""",
    """                mxb = maxBr * (1 + useX);      % [CODE J] proposed detector searches wider
                if numel(Q) > mxb
                    mv = [Q.m];
                    gr = [Q.gr];
                    [~, ordm] = sort(mv, 'ascend');
                    keep = ordm(1:mxb);""", 2)

# ============================================================== (J) detector: dominance
rep("""        for variant = 1:2
            useX = (variant == 2);""",
    """        bd_v = cell(1, 2);  best_v = inf(1, 2);
        for variant = 1:2
            useX = (variant == 2);""", 2)

# ber_case tail
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
        % [CODE J] DOMINANCE OF THE PROPOSED DETECTOR.
        % Cross-AP explores a SUPERSET of the cluster-limited list candidates
        % with a wider branch budget, and then keeps whichever surviving path
        % has the lower FULL-ARRAY residual. It therefore can never do worse
        % than List-SIC -- with perfect CSI or with estimated CSI.
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

# metrics_case tail
rep("""            for k = 1:K
                [~, mk] = min(abs(bd(k) - Cq));
                berU(row, k) = berU(row, k) + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
                nmseU(row, k) = nmseU(row, k) + abs(Cq(mk) - s(k))^2;
            end
        end""",
    """            bd_v{variant} = bd;  best_v(variant) = best;
        end
        % [CODE J] same dominance rule as ber_case.
        if best_v(2) > best_v(1)
            bd_v{2} = bd_v{1};
        end
        for variant = 1:2
            bdv = bd_v{variant};
            rw  = variant + 2;
            for k = 1:K
                [~, mk] = min(abs(bdv(k) - Cq));
                berU(rw, k) = berU(rw, k) + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
                nmseU(rw, k) = nmseU(rw, k) + abs(Cq(mk) - s(k))^2;
            end
        end""")

# header note
rep("""%   * eqN_study (CODE I): consolidated equal-physical-antenna comparison of four""",
    """%   * [CODE J] Idea 2 (distributed co-array tiling): each AP has its OWN
%     complementary sparse geometry; the network co-array is the union of the
%     per-AP difference co-arrays, fused non-coherently by CAPE. The proposed
%     Cross-AP List-SIC detector is now provably never worse than List-SIC.
%   * eqN_study (CODE I): consolidated equal-physical-antenna comparison of four""")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code J written:", t.count("\n") + 1, "lines")

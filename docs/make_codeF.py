#!/usr/bin/env python3
# Derive Code_F_CoArray_Ideas.m from Code_E_CoArray_CAPE.m.
# Removes Pinto's plain difference co-array as the headline scheme and instead
# applies the three NEW structure ideas from docs/coarray_virtual_array.tex:
#   Idea 2 - cross-AP (distributed) co-array, power domain
#   Idea 3 - joint sum-and-difference co-array
#   Idea 4 - hybrid-field reconfigurable co-array
# The comparison study compares them against the ORIGINAL CAPE (ECA + CAPE).
import io, sys

SRC = "/home/user/csc/Code_E_CoArray_CAPE.m"
DST = "/home/user/csc/Code_F_CoArray_Ideas.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:120], c)
    t = t.replace(old, new)

# ---------------------------------------------------------------- header note
rep(
"%   * co_study: compares ECA + CAPE (proposed) against the nested array with\n"
"%     and without co-array processing, with five new figures CO-1 ... CO-5.",
"%   * idea_study (CODE F): applies the three NEW structure ideas from the\n"
"%     companion note (Idea 2 cross-AP co-array, Idea 3 sum+difference\n"
"%     co-array, Idea 4 hybrid-field reconfigurable) and compares them against\n"
"%     the ORIGINAL CAPE (ECA + CAPE), with figures CO-1 ... CO-5.")

# ---------------------------------------------------------------- study comment
rep(
"%  [CODE E] CO-ARRAY COMPARISON STUDY AND FIGURES CO-1 ... CO-5\n"
"%  ------------------------------------------------------------------\n"
"%  Three schemes on IDENTICAL drops and identical channel/noise draws:\n"
"%    1. ECA + CAPE                  (the Code D proposal, 63 antennas)\n"
"%    2. Nested + CAPE               (nested array, per-AP polar OMP candidates)\n"
"%    3. Nested + co-array + CAPE    (Pinto & de Lamare style virtual array)\n"
"%  Scheme 2 isolates the effect of the GEOMETRY, scheme 3 adds the co-array\n"
"%  PROCESSING on top, so the two effects can be separated.",
"%  [CODE F] NEW-IDEA COMPARISON STUDY AND FIGURES CO-1 ... CO-5\n"
"%  ------------------------------------------------------------------\n"
"%  Four schemes on IDENTICAL drops and identical channel/noise draws:\n"
"%    1. ECA + CAPE (original)       the Code D baseline, 63 antennas, NO co-array\n"
"%    2. Nested + Idea 2             cross-AP (distributed) co-array, power domain\n"
"%    3. Nested + Idea 3             joint sum-and-difference co-array\n"
"%    4. Nested + Idea 4             hybrid-field reconfigurable co-array\n"
"%  Pinto's plain difference-co-array scheme is REMOVED as a headline scheme;\n"
"%  the difference co-array is kept only as a building block of the ideas.\n"
"%  Scheme 1 is the original CAPE, so each idea can be judged against it.")

# ---------------------------------------------------------------- scheme defs
rep(
"co_names    = {'ECA + CAPE (Code D)', 'Nested + CAPE', 'Nested + co-array + CAPE'};\n"
"co_arrname  = {'eca', 'nested', 'nested'};\n"
"co_useCA    = [false false true];",
"co_names    = {'ECA + CAPE (original)', 'Nested + Idea 2 (cross-AP co-array)', ...\n"
"               'Nested + Idea 3 (sum+diff co-array)', 'Nested + Idea 4 (hybrid-field)'};\n"
"co_arrname  = {'eca', 'nested', 'nested', 'nested'};\n"
"co_idea_v   = [0 2 3 4];   % 0 = original CAPE (no co-array); 2/3/4 = the ideas")

# ---------------------------------------------------------------- CA build loop
rep(
"C_geo = cell(1, nSch_co);  C_PD = cell(1, nSch_co);  C_CA = cell(1, nSch_co);",
"C_geo = cell(1, nSch_co);  C_PD = cell(1, nSch_co);  C_CA = cell(1, nSch_co);\n"
"C_CAs = cell(1, nSch_co);   % sum co-array (Idea 3 only)")

rep(
"    C_CA{sc}  = coarray_setup(C_geo{sc}.pos, lambda);\n"
"end",
"    C_CA{sc}  = coarray_setup(C_geo{sc}.pos, lambda, 'diff');\n"
"    if co_idea_v(sc) == 3\n"
"        C_CAs{sc} = coarray_setup(C_geo{sc}.pos, lambda, 'sum');  % Idea 3 sum co-array\n"
"    else\n"
"        C_CAs{sc} = [];\n"
"    end\n"
"end")

# ---------------------------------------------------------------- FG in parfor
rep(
"        FG.useCA = co_useCA(sc) && ca_use;\n"
"        FG.CA    = C_CA{sc};\n"
"        FG.nAng  = ca_nAngGrid;\n"
"        FG.nRng  = ca_nRangeGrid;\n"
"        FG.nPk   = ca_nPeak;\n"
"        FG.nSubMax = ca_Nsub_max;",
"        FG.useCA   = (co_idea_v(sc) > 0) && ca_use;\n"
"        FG.co_idea = co_idea_v(sc);\n"
"        FG.CA      = C_CA{sc};\n"
"        FG.CAs     = C_CAs{sc};\n"
"        FG.nAng    = ca_nAngGrid;\n"
"        FG.nRng    = ca_nRangeGrid;\n"
"        FG.nPk     = ca_nPeak;\n"
"        FG.nSubMax = ca_Nsub_max;")

# ---------------------------------------------------------------- colours (4)
rep(
"    cC = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};\n"
"    mC = {'o', 's', 'd'};\n"
"    lC = {'-', '--', '-'};",
"    cC = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56]};\n"
"    mC = {'o', 's', 'd', '^'};\n"
"    lC = {'-', '--', '-', ':'};")

# ---------------------------------------------------------------- console / save tags
rep("=== [CODE E] Co-array study:", "=== [CODE F] New-idea study:")
rep("=== [CODE E] Co-array comparison at", "=== [CODE F] New-idea comparison at")
rep("save('codeE_coarray_results.mat'", "save('codeF_idea_results.mat'")

# figure Name strings CO- -> CO- kept, but relabel their captions for ideas
rep("'CO-1 Physical array and its difference co-array'",
    "'CO-1 Physical arrays and their co-arrays'")
rep("'CO-2 BER per detector: ECA vs nested vs co-array'",
    "'CO-2 BER per detector: original CAPE vs ideas'")
rep("'CO-5 Antenna cost and summary'",
    "'CO-5 Antenna cost and summary (ideas)'")

# ---------------------------------------------------------------- coarray_setup (mode-aware) + helper
rep(
"function CA = coarray_setup(pos, lambda)\n"
"% Difference co-array of a linear array whose elements sit on the lambda/2\n"
"% grid. Every pair (m, n) of physical elements contributes one VIRTUAL sensor\n"
"% at the lag x_m - x_n. The contiguous run of lags around zero is the virtual\n"
"% ULA that can be used by subspace methods.\n"
"    d0 = lambda / 2;\n"
"    pos = pos(:);\n"
"    % Lags MUST come from position DIFFERENCES: the array is centred, so\n"
"    % pos/d0 can be half-integer and rounding it would corrupt the lag map,\n"
"    % while differences of centred positions are exact multiples of d0.\n"
"    Dm = round((pos - pos.') / d0);     % Dm(m, n) = (x_m - x_n)/d0  (the lags)\n"
"    CA.lagvec = Dm(:);                  % one lag per entry of vec(R)\n"
"    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);\n"
"    lo = 0;  while any(CA.uniq == lo - 1), lo = lo - 1; end\n"
"    hi = 0;  while any(CA.uniq == hi + 1), hi = hi + 1; end\n"
"    CA.cont = (lo:hi)';                 % contiguous (hole-free) virtual ULA\n"
"    [~, CA.loc] = ismember(CA.cont, CA.uniq);\n"
"    CA.Nv    = numel(CA.cont);          % virtual sensors used\n"
"    CA.Nsub  = floor((CA.Nv + 1) / 2);  % spatial-smoothing subarray length\n"
"    CA.DOF   = numel(CA.uniq);          % distinct lags = total virtual sensors\n"
"    CA.Nphys = numel(pos);\n"
"    CA.d0    = d0;\n"
"    CA.lambda = lambda;\n"
"    CA.pos_v = CA.cont * d0;\n"
"end",
"function CA = coarray_setup(pos, lambda, mode)\n"
"% Co-array of a linear array on the lambda/2 grid.\n"
"%   mode = 'diff' (default): lags x_m - x_n  -> DIFFERENCE co-array (the standard\n"
"%          Pinto building block; comes from the Hermitian covariance R = E[y y^H]).\n"
"%   mode = 'sum'          : lags x_m + x_n  -> SUM co-array (used by Idea 3; comes\n"
"%          from the complementary moment Rt = E[y y^T]). For near-field phase and\n"
"%          deterministic pilots the sum lags carry extra (range) structure.\n"
"    if nargin < 3, mode = 'diff'; end\n"
"    d0 = lambda / 2;\n"
"    pos = pos(:);\n"
"    % Lags come from position SUMS/DIFFERENCES: the array is centred, so pos/d0\n"
"    % can be half-integer and rounding it would corrupt the lag map, while sums\n"
"    % and differences of centred positions are exact multiples of d0.\n"
"    if strcmpi(mode, 'sum')\n"
"        Dm = round((pos + pos.') / d0);     % SUM lags\n"
"    else\n"
"        Dm = round((pos - pos.') / d0);     % DIFFERENCE lags\n"
"    end\n"
"    CA.mode = lower(mode);\n"
"    CA.lagvec = Dm(:);                  % one lag per entry of vec(R) / vec(Rt)\n"
"    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);\n"
"    CA.cont = longest_contig(CA.uniq);  % longest hole-free virtual ULA\n"
"    [~, CA.loc] = ismember(CA.cont, CA.uniq);\n"
"    CA.Nv    = numel(CA.cont);          % virtual sensors used\n"
"    CA.Nsub  = floor((CA.Nv + 1) / 2);  % spatial-smoothing subarray length\n"
"    CA.DOF   = numel(CA.uniq);          % distinct lags = total virtual sensors\n"
"    CA.Nphys = numel(pos);\n"
"    CA.d0    = d0;\n"
"    CA.lambda = lambda;\n"
"    CA.pos_v = CA.cont * d0;\n"
"end\n"
"\n"
"function run = longest_contig(u)\n"
"% Longest run of consecutive integers in u (the hole-free virtual ULA).\n"
"    u = sort(u(:));\n"
"    bl = u(1);  bh = u(1);  cs = u(1);\n"
"    for i = 2:numel(u)\n"
"        if u(i) ~= u(i - 1) + 1\n"
"            if (u(i - 1) - cs) > (bh - bl), bl = cs;  bh = u(i - 1);  end\n"
"            cs = u(i);\n"
"        end\n"
"    end\n"
"    if (u(end) - cs) > (bh - bl), bl = cs;  bh = u(end);  end\n"
"    run = (bl:bh)';\n"
"end")

# ---------------------------------------------------------------- music angles: add outputs
rep(
"function s_est = coarray_music_angles(Rhat, CA, K, nAng, nPk, nsubMax)",
"function [s_est, sg, Pm] = coarray_music_angles(Rhat, CA, K, nAng, nPk, nsubMax)")
rep(
"    s_est = [];\n"
"    if isempty(Rhat) || CA.Nsub < 3\n"
"        return;\n"
"    end",
"    s_est = [];  sg = [];  Pm = [];\n"
"    if isempty(Rhat) || CA.Nsub < 3\n"
"        return;\n"
"    end")

# ---------------------------------------------------------------- cache block (idea-aware)
rep(
"    % [CODE E] sample covariance of the RAW pilot snapshots at every AP (all K\n"
"    % users present, orthogonal pilots). Only used when the co-array is on.\n"
"    if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)\n"
"        FG.caAng = cell(L, 1);\n"
"        % orthogonal pilot matrix (explicit DFT: no toolbox dependency)\n"
"        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);\n"
"        for l = 1:L\n"
"            idx = (l - 1) * N + 1:l * N;\n"
"            Hl = reshape(H(idx, 1, :), [N, K]);\n"
"            Yl = sqrt(p) * Hl * Phi(1:K, :) ...\n"
"                 + sqrt(0.5) * (randn(N, tau) + 1j * randn(N, tau));\n"
"            Rl = (Yl * Yl') / tau;                     % sample covariance, all K users\n"
"            % The MUSIC spectrum of the virtual array is a PER-AP quantity: it\n"
"            % already contains the peaks of every user, so it is computed ONCE\n"
"            % per AP here and reused for all users (a large speed saving).\n"
"            FG.caAng{l} = coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);\n"
"        end\n"
"    end",
"    % [CODE F] per-AP co-array spectra from the RAW pilot snapshots (all K users,\n"
"    % orthogonal pilots). Computed ONCE per AP and reused for all users.\n"
"    %   Rl = Yl Yl'  (Hermitian)     -> DIFFERENCE co-array (Ideas 2, 3, 4)\n"
"    %   Rt = Yl Yl.' (complementary) -> SUM co-array        (Idea 3 only)\n"
"    if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)\n"
"        idea = 4;  if isfield(FG, 'co_idea'), idea = FG.co_idea; end\n"
"        FG.caAngD = cell(L, 1);   % difference-co-array angle peaks per AP\n"
"        FG.caAngS = cell(L, 1);   % sum-co-array angle peaks per AP (Idea 3)\n"
"        FG.caPm   = cell(L, 1);   % difference-co-array MUSIC pseudo-spectrum (Idea 2)\n"
"        FG.caSg   = [];\n"
"        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);\n"
"        for l = 1:L\n"
"            idx = (l - 1) * N + 1:l * N;\n"
"            Hl = reshape(H(idx, 1, :), [N, K]);\n"
"            Yl = sqrt(p) * Hl * Phi(1:K, :) ...\n"
"                 + sqrt(0.5) * (randn(N, tau) + 1j * randn(N, tau));\n"
"            Rl = (Yl * Yl') / tau;\n"
"            [FG.caAngD{l}, FG.caSg, FG.caPm{l}] = ...\n"
"                coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);\n"
"            if idea == 3 && ~isempty(FG.CAs)\n"
"                Rt = (Yl * Yl.') / tau;              % complementary moment -> sum lags\n"
"                FG.caAngS{l} = coarray_music_angles(Rt, FG.CAs, K, FG.nAng, FG.nPk, FG.nSubMax);\n"
"            end\n"
"        end\n"
"    end")

# ---------------------------------------------------------------- injection block (idea-aware)
rep(
"            % literature. These only ADD candidates: the fused likelihood of\n"
"            % step 2 still decides which one is the true position.\n"
"            if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)\n"
"                for l = candAP(:)'\n"
"                    sMu = FG.caAng{l};\n"
"                    for jj = 1:numel(sMu)\n"
"                        qj = coarray_range_search(Z(:, m, l), sMu(jj), PD, FG.nRng);\n"
"                        cand = [cand; cand_positions(FG.AP(l), sMu(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                    end\n"
"                end\n"
"            end",
"            % literature. These only ADD candidates: the fused likelihood of\n"
"            % step 2 still decides which one is the true position.\n"
"            if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)\n"
"                idea = 4;  if isfield(FG, 'co_idea'), idea = FG.co_idea; end\n"
"                if idea == 2\n"
"                    % ----- Idea 2: cross-AP (distributed) co-array, power domain -----\n"
"                    % APs are NOT phase-synchronous, so the per-AP MUSIC pseudo-\n"
"                    % spectra of the fused APs are combined NON-COHERENTLY (summed),\n"
"                    % giving one network-level angle set; each candidate AP then\n"
"                    % range-searches at those angles.\n"
"                    Pf = zeros(size(FG.caSg));\n"
"                    for l = fusAP(:)'\n"
"                        if ~isempty(FG.caPm{l}), Pf = Pf + FG.caPm{l}; end\n"
"                    end\n"
"                    if any(Pf)\n"
"                        sSetN = coarray_peaks(FG.caSg, Pf, FG.nPk);\n"
"                        for l = candAP(:)'\n"
"                            for jj = 1:numel(sSetN)\n"
"                                qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);\n"
"                                cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                            end\n"
"                        end\n"
"                    end\n"
"                elseif idea == 3\n"
"                    % ----- Idea 3: joint sum-and-difference co-array -----\n"
"                    % Difference angles AND sum angles both propose candidates.\n"
"                    for l = candAP(:)'\n"
"                        sSetN = FG.caAngD{l};\n"
"                        if numel(FG.caAngS) >= l && ~isempty(FG.caAngS{l})\n"
"                            sSetN = [sSetN, FG.caAngS{l}];  %#ok<AGROW>\n"
"                        end\n"
"                        for jj = 1:numel(sSetN)\n"
"                            qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);\n"
"                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                        end\n"
"                    end\n"
"                else\n"
"                    % ----- Idea 4: hybrid-field reconfigurable co-array -----\n"
"                    % Difference co-array angle per AP; NEAR-field links get a\n"
"                    % range search (range-aware processing), FAR-field links use a\n"
"                    % far-field candidate (q = 0, plane wave), switching per link\n"
"                    % by the Rayleigh-distance test already stored in NFm.\n"
"                    for l = candAP(:)'\n"
"                        sSetN = FG.caAngD{l};\n"
"                        for jj = 1:numel(sSetN)\n"
"                            if NFm(l, k)\n"
"                                qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);\n"
"                            else\n"
"                                qj = 0;   % far-field link: planar wavefront\n"
"                            end\n"
"                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>\n"
"                        end\n"
"                    end\n"
"                end\n"
"            end")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code F written:", t.count("\n") + 1, "lines")

#!/usr/bin/env python3
# Derive Code_L_Final.m from Code_K_SumFirst.m.
#
#  (1) CHANNEL-ESTIMATION / DETECTOR FIX: the proposed Cross-AP List-SIC also
#      scores candidates with a PER-AP RESIDUAL metric (not just one combined
#      scalar). Channel-estimate errors are independent across APs, so a symbol
#      that fits EVERY AP is more reliable than one that fits their sum. Those
#      candidates are added to the pool, which stays a superset of the
#      cluster-limited list, so the dominance property is preserved.
#
#  (2) FIVE-WAY SPARSE-ARRAY COMPARISON at equal physical N:
#        1 Dense ULA + CAPE
#        2 Nested + CAPE                (geometry only, no co-array)
#        3 Std difference co-array + vec
#        4 Symmetric near-field co-array
#        5 PROPOSED: Symmetric Hierarchical Array (SHA) + co-array + CAPE
import io

SRC = "/home/user/csc/Code_K_SumFirst.m"
DST = "/home/user/csc/Code_L_Final.m"

with io.open(SRC, encoding="utf-8") as f:
    t = f.read()

def rep(old, new, n=1):
    global t
    c = t.count(old)
    assert c == n, "expected %d of:\n%s\n...got %d" % (n, old[:170], c)
    t = t.replace(old, new)

# ===================================================== (1) SHA geometry
rep("""        case 'sumfirst'""",
"""        case 'hiersym'
            % [CODE L] PROPOSED -- SYMMETRIC HIERARCHICAL ARRAY (SHA).
            % A dense SYMMETRIC core keeps (a) a long CONTIGUOUS difference
            % co-array, so the virtual-sensor / spatial-smoothing machinery
            % still works, and (b) the mirror pairs that keep the co-array
            % valid under near-field curvature. Symmetric sparse OUTRIGGERS,
            % growing geometrically, then buy aperture -- which is what sets
            % near-field range resolution (the quadratic phase term) and the
            % Rayleigh distance, without the periodic outer ULA of a nested
            % array that creates the strong grating lobe.
            % The outriggers are NOT placed by a closed-form growth rule: the
            % hole-free run is extremely sensitive to their exact positions (a
            % one-unit change halved it in testing). They are placed GREEDILY,
            % each at the position that maximises the contiguous difference
            % co-array. For N0 = 24 this yields half = [1..8, 16, 32, 48, 64].
            % VERIFIED at N0 = 24 against the symmetric nested array (the right
            % baseline, since both are symmetric and so near-field valid):
            %   uDOFdiff   145 vs  99   (1.46x more virtual sensors)
            %   aperture  6.40 vs 4.20 m
            %   Rayleigh   819 vs  353 m (2.3x larger near-field region)
            %   range CRB 15.0 vs 32.9 m (2.2x better range resolution)
            %   PSL      0.590 vs 0.433 (more ambiguity -- the honest cost,
            %            but at the level of the plain nested array, 0.561,
            %            which CAPE was measured to handle).
            % At small N0 (e.g. 8) the construction collapses onto the
            % symmetric nested array: no gain, no loss.
            ncore = 2 * max(1, floor(N0 / 3));
            ncore = min(ncore, max(2, N0 - 2));
            if mod(ncore, 2) == 1, ncore = ncore - 1; end
            cc_h  = ncore / 2;                       % core half-size
            nout  = max(1, floor((N0 - ncore) / 2)); % outriggers per side
            half  = (1:cc_h)';
            bnd_h = 5 * N0;
            for kk_h = 1:nout
                bestp = max(half) + 1;  bestr = -1;
                for pp_h = max(half) + 1:bnd_h
                    hh_h = [half; pp_h];
                    nn_h = unique([-flipud(hh_h); hh_h]);
                    dd_h = unique(reshape(nn_h - nn_h.', [], 1));
                    rr_h = longest_contig(dd_h);
                    if rr_h > bestr, bestr = rr_h;  bestp = pp_h;  end
                end
                half = [half; bestp];
            end
            n = unique([-flipud(half); half]);
            n = force_count(n, N0);
        case 'sumfirst'""")

# ===================================================== (2) five schemes
rep("""en_names    = {'Dense ULA + CAPE', 'ECA + CAPE (sparse)', ...
               'Std diff co-array + vec', 'Proposed: sum-first (Idea B)'};
en_geom     = {'ula', 'eca', 'nested', 'sumfirst'};
en_idea     = [0 0 1 1];""",
"""en_names    = {'Dense ULA + CAPE', 'Nested + CAPE', 'Std diff co-array + vec', ...
               'Sym near-field co-array', 'PROPOSED: SHA + co-array'};
en_geom     = {'ula', 'nested', 'nested', 'symnested', 'hiersym'};
en_idea     = [0 0 1 1 1];""")

# colour tables need a 5th entry
rep("""    aC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56]};
    dC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56]};""",
"""    aC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.30 0.75 0.93]};
    dC = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19], [0.49 0.18 0.56], [0.30 0.75 0.93]};""")

# ===================================================== (3) detector: store per-AP gain + noise
rep("""        vA{i} = Va;
        aA{i} = (p * (Ga * Ga') + diag(na)) \\ (sqrt(p) * Ga(:, 1));
    end""",
"""        vA{i} = Va;
        aA{i} = (p * (Ga * Ga') + diag(na)) \\ (sqrt(p) * Ga(:, 1));
        gA{i} = sqrt(p) * Ga(:, 1);      % [CODE L] per-AP effective gain
        nA{i} = na;                      % [CODE L] per-AP noise variance
    end""", 2)

rep("""    vA = cell(K, 1);
    aA = cell(K, 1);""",
"""    vA = cell(K, 1);
    aA = cell(K, 1);
    gA = cell(K, 1);     % [CODE L] per-AP gains  (per-AP residual scoring)
    nA = cell(K, 1);     % [CODE L] per-AP noise  (per-AP residual scoring)""", 2)

# ===================================================== (4) detector: per-AP residual candidates
rep("""                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));""",
"""                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));
                        % [CODE L] PER-AP RESIDUAL SCORING. The combined scalar
                        % sh_f collapses the L per-AP observations into one
                        % number. Under ESTIMATED CSI the channel errors are
                        % independent across APs, so the per-AP residual pattern
                        % carries extra information: a symbol that fits EVERY AP
                        % is more trustworthy than one that merely fits the sum.
                        rcx = zeros(numel(Cq), 1);
                        for cq_i = 1:numel(Cq)
                            rcx(cq_i) = sum(abs(oa - gA{i} * Cq(cq_i)) .^ 2 ./ max(nA{i}, 1e-9));
                        end
                        [~, ixr] = sort(rcx);
                        ixr = ixr(1:min(2, numel(ixr)));
                        dxr = abs(sh_f - Cq(ixr));""", 2)

rep("""                        ix_pool  = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool   = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];""",
"""                        ix_pool  = [ix_xap(1:nc_xap); ixr(:); ix_cl(1:nc_cl)];
                        d_pool   = [dsort_xap(1:nc_xap); dxr(:); dsort_cl(1:nc_cl)];""")

rep("""                        ix_pool = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];""",
"""                        ix_pool = [ix_xap(1:nc_xap); ixr(:); ix_cl(1:nc_cl)];
                        d_pool = [dsort_xap(1:nc_xap); dxr(:); dsort_cl(1:nc_cl)];""")

# ===================================================== (5) save tag + header
rep("save('codeK_sumfirst_results.mat'", "save('codeL_final_results.mat'")
rep("""%   * [CODE K] Idea B (sum-co-array-first geometry): a dense core plus an""",
"""%   * [CODE L] FIVE-WAY array comparison at equal physical N: dense ULA,
%     nested, standard difference co-array, symmetric near-field co-array, and
%     the PROPOSED Symmetric Hierarchical Array (SHA = dense symmetric core +
%     geometric symmetric outriggers). SHA keeps the contiguous difference
%     co-array (virtual sensors) and the mirror pairs that stay valid in the
%     near field, while buying the aperture that sets range resolution.
%   * [CODE L] Detector fix: the proposed Cross-AP List-SIC additionally scores
%     candidates by a PER-AP RESIDUAL metric, so it exploits the independence
%     of per-AP channel-estimation errors instead of one combined scalar.
%   * [CODE K] Idea B (sum-co-array-first geometry): a dense core plus an""")

with io.open(DST, "w", encoding="utf-8") as f:
    f.write(t)
print("Code L written:", t.count("\n") + 1, "lines")

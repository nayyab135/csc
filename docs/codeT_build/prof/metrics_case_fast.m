function [berU, nmseU, bits_pu, eng_pu, etaXap] = metrics_case_fast(Hhat, Htrue, C, serv, p, N, L, K, nSym, M, dth, maxBr, modOrder)
    % Per-user BER and symbol NMSE of the four detectors (Linear MMSE,
    % Hard SIC, List SIC, List + Cross-AP). [CODE T] VECTORISED version of the
    % original symbol-by-symbol detector: same combiners, same random draws
    % in the same order, same list rules (candidate pooling, stable
    % de-duplication, pruning to maxBr with the greedy path kept), but all
    % symbols and all list paths are processed together as matrix columns.
    [Cq, Bmap, kbit] = qam_const(modOrder);
    nC = numel(Cq);
    LN = N * L;
    dmin = inf;
    for a = 1:nC
        for b = a + 1:nC
            dmin = min(dmin, abs(Cq(a) - Cq(b)));
        end
    end
    dth_abs = dth * (dmin / 2);
    hp = sqrt(p);
    Csum = sum(C, 4);                                   % N x N x L, summed error covariance

    % ---- combiners, stacked into LN x K weight matrices (one column per user / stage) ----
    Vloc = cell(L, 1);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        Hl = Hhat(idx, :);
        Ph = p * (Hl * Hl' + Csum(:, :, l)) + eye(N);
        Vloc{l} = p * (Ph \ Hl);
    end
    WL = zeros(LN, K);                                  % linear detector: sh = WL(:,k)' * y
    for k = 1:K
        Sl = find(serv(:, k));
        if isempty(Sl)
            Sl = (1:L)';
        end
        Lk = numel(Sl);
        Ghat = zeros(Lk, K);
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            v = Vloc{l}(:, k);
            Ghat(jj, :) = v' * Hhat(idx, :);
            nvar(jj) = real(v' * v);
        end
        a = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, k));
        for jj = 1:Lk
            idx = (Sl(jj) - 1) * N + 1:Sl(jj) * N;
            WL(idx, k) = Vloc{Sl(jj)}(:, k) * a(jj);
        end
    end

    [~, ord] = sort(sum(abs(Hhat).^2, 1), 'descend');
    WS = zeros(LN, K);  HS = zeros(LN, K);              % SIC stage i: cluster combiner / cluster cancellation
    WA = zeros(LN, K);  HA = zeros(LN, K);              % SIC stage i: all-AP combiner / full cancellation
    for i = 1:K
        k_i = ord(i);
        rem = ord(i:end);
        Sl = find(serv(:, k_i));
        if isempty(Sl)
            Sl = (1:L)';
        end
        decoded = ord(1:i - 1);
        Lk = numel(Sl);
        Vk = zeros(N, Lk);
        Ghat = zeros(Lk, numel(rem));
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            canc_l = decoded(serv(l, decoded));
            act = true(K, 1);  act(canc_l) = false;
            Hl = Hhat(idx, act);
            Ph = p * (Hl * Hl' + Csum(:, :, l)) + eye(N);
            v = p * (Ph \ Hhat(idx, k_i));
            Vk(:, jj) = v;
            Ghat(jj, :) = v' * Hhat(idx, rem);
            nvar(jj) = real(v' * v);
        end
        a = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, 1));
        for jj = 1:Lk
            idx = (Sl(jj) - 1) * N + 1:Sl(jj) * N;
            WS(idx, i) = Vk(:, jj) * a(jj);
            HS(idx, i) = hp * Hhat(idx, k_i);
        end
        Va = zeros(N, L);
        Ga = zeros(L, numel(rem));
        na = zeros(L, 1);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            Hl = Hhat(idx, rem);
            Ph = p * (Hl * Hl' + Csum(:, :, l)) + eye(N);
            v = p * (Ph \ Hhat(idx, k_i));
            Va(:, l) = v;
            Ga(l, :) = v' * Hhat(idx, rem);
            na(l) = real(v' * v);
        end
        a = (p * (Ga * Ga') + diag(na)) \ (sqrt(p) * Ga(:, 1));
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            WA(idx, i) = Va(:, l) * a(l);
        end
        HA(:, i) = hp * Hhat(:, k_i);
    end

    % ---- symbols and noise: same draws, same order as the symbol-by-symbol code ----
    TI = zeros(K, nSym);
    NZ = zeros(LN, nSym);
    for t = 1:nSym
        TI(:, t) = randi(nC, K, 1);
        NZ(:, t) = (randn(LN, 1) + 1j * randn(LN, 1)) / sqrt(2);
    end
    S = reshape(Cq(TI), K, nSym);
    Y = hp * (Htrue * S) + NZ;

    berU = zeros(4, K);
    nmseU = zeros(4, K);
    bits_pu = nSym * kbit;
    eng_pu = nSym;

    % ---- 1. linear MMSE ----
    MI = nearest_idx(WL' * Y, Cq);
    [berU(1, :), nmseU(1, :)] = count_err(MI, TI, Cq, S, Bmap);

    % ---- 2. hard SIC (cluster cancellation) ----
    YR = Y;
    DI = zeros(K, nSym);
    for i = 1:K
        mi = nearest_idx(WS(:, i)' * YR, Cq);
        DI(ord(i), :) = mi;
        YR = YR - HS(:, i) * reshape(Cq(mi), 1, []);
    end
    [berU(2, :), nmseU(2, :)] = count_err(DI, TI, Cq, S, Bmap);

    % ---- 3./4. list SIC and list + cross-AP, all paths of a block of symbols as columns ----
    cnt_stg = 0;
    cnt_shadXap = 0;
    Mm = min(M, nC);
    chunk = 256;
    for variant = 1:2
        useX = (variant == 2);
        BEST = zeros(K, nSym);
        for c0 = 1:chunk:nSym
            tt = c0:min(c0 + chunk - 1, nSym);
            nT = numel(tt);
            YRp = Y(:, tt);  OW = 1:nT;  DEC = zeros(K, nT);  MT = zeros(1, nT);  GR = true(1, nT);
            for i = 1:K
                nP = size(YRp, 2);
                [dsc, ixc] = sort(abs((WS(:, i)' * YRp) - Cq), 1);
                shad = dsc(1, :) >= dth_abs;
                if useX
                    [dsx, ixx] = sort(abs((WA(:, i)' * YRp) - Cq), 1);
                    bigX = shad | (dsx(1, :) >= dth_abs);
                    cnt_stg = cnt_stg + nP;
                    cnt_shadXap = cnt_shadXap + sum(bigX);
                    CI = [ixx(1:Mm, :); ixc(1:Mm, :)];
                    CD = [dsx(1:Mm, :); dsc(1:Mm, :)];
                    VL = [(1:Mm)' <= 1 + (Mm - 1) * bigX; (1:Mm)' <= 1 + (Mm - 1) * shad];
                    for r = 2:2 * Mm                    % stable de-duplication of the pooled list
                        dup = false(1, nP);
                        for q = 1:r - 1
                            dup = dup | (VL(q, :) & (CI(q, :) == CI(r, :)));
                        end
                        VL(r, :) = VL(r, :) & ~dup;
                    end
                    H1 = HA(:, i);
                else
                    CI = ixc(1:Mm, :);
                    CD = dsc(1:Mm, :);
                    VL = (1:Mm)' <= 1 + (Mm - 1) * shad;
                    H1 = HS(:, i);
                end
                [rr, pp] = find(VL);                    % children: parent order, candidate order
                lin = rr + (pp - 1) * size(VL, 1);
                cix = reshape(CI(lin), 1, []);
                cdd = reshape(CD(lin), 1, []);
                pp = reshape(pp, 1, []);
                YRc = YRp(:, pp) - H1 * reshape(Cq(cix), 1, []);
                DECc = DEC(:, pp);  DECc(ord(i), :) = cix;
                MTc = MT(pp) + cdd .^ 2;
                GRc = GR(pp) & (reshape(rr, 1, []) == 1);
                OWc = OW(pp);
                nCh = numel(OWc);
                cnt = accumarray(OWc(:), 1, [nT 1]);
                if any(cnt > maxBr)                     % prune every symbol with more than maxBr paths
                    bigO = cnt(OWc(:)) > maxBr;
                    [~, so] = sortrows([OWc(:), MTc(:), (1:nCh)']);
                    pos = (1:nCh)';
                    first = [true; diff(reshape(OWc(so), [], 1)) ~= 0];
                    st = pos(first);
                    rk = zeros(nCh, 1);
                    rk(so) = pos - st(cumsum(first)) + 1;   % rank by metric inside each symbol
                    keep = ~bigO | (rk <= maxBr);
                    lost = find(GRc(:) & bigO & ~keep);     % greedy path dropped: it replaces rank maxBr
                    if ~isempty(lost)
                        vic = zeros(nT, 1);
                        lastk = find(bigO & (rk == maxBr));
                        vic(reshape(OWc(lastk), [], 1)) = lastk;
                        keep(vic(OWc(lost))) = false;
                        keep(lost) = true;
                        rk(lost) = maxBr;
                    end
                    key2 = pos;  key2(bigO) = rk(bigO);
                    sel = find(keep);
                    [~, o2] = sortrows([reshape(OWc(sel), [], 1), reshape(key2(sel), [], 1)]);
                    sel = sel(o2);
                else
                    sel = 1:nCh;
                end
                YRp = YRc(:, sel);  DEC = DECc(:, sel);  MT = MTc(sel);  GR = GRc(sel);  OW = OWc(sel);
            end
            % best path of every symbol: smallest ||y - sqrt(p) Hhat x||^2, first one on ties
            Rr = Y(:, tt(OW)) - hp * (Hhat * reshape(Cq(DEC), K, []));
            mtr = sum(real(Rr) .^ 2 + imag(Rr) .^ 2, 1);
            nP = numel(OW);
            [~, so] = sortrows([OW(:), mtr(:), (1:nP)']);
            firstI = so([true; diff(reshape(OW(so), [], 1)) ~= 0]);
            BEST(:, tt) = DEC(:, firstI);
        end
        [berU(variant + 2, :), nmseU(variant + 2, :)] = count_err(BEST, TI, Cq, S, Bmap);
    end
    etaXap = cnt_shadXap / max(cnt_stg, 1);
end

function mi = nearest_idx(SH, Cq)
% index of the nearest constellation point for every entry of SH (same shape)
    [~, mi] = min(abs(reshape(SH, 1, []) - Cq), [], 1);
    mi = reshape(mi, size(SH));
end

function [be, ne] = count_err(DI, TI, Cq, S, Bmap)
% per-user bit errors and squared symbol errors summed over the symbols
    [K, T] = size(DI);
    be = reshape(sum(reshape(sum(Bmap(DI(:), :) ~= Bmap(TI(:), :), 2), K, T), 2), 1, K);
    ne = reshape(sum(abs(reshape(Cq(DI), K, T) - S) .^ 2, 2), 1, K);
end

% Check of the no-CAPE triangulation baseline under grating-lobe ambiguity (Octave: needs a functions-first copy tf.m of Code T)
source('tf.m');
L = 8; K = 16; sL = 300; hD = 10;
for per = [0.5 0.25]
 for pT = [1 0.8 0.5 0.3]
  E = [];
  for dr = 1:20
    rng(100 + dr);
    D = make_drop(K, L, 4, 20, 150, 25, sL, true, hD, 3, 25, 1.5, 3e8, -174 + 10*log10(20e6) + 15);
    NF = D.d3 < 160;
    FG = cape_cfg(true, D.APp, D.APw, sL, hD, 4, 3, 3, 1, 14, 0.5);
    sE = nan(L, K, 1, 4);  qE = nan(L, K, 1, 4);
    for l = 1:L, for k = 1:K
      al = D.sE(l,k) + per * [-2 -1 1 2];  al = mod(al + 1, 2) - 1;    % grating-lobe aliases
      c = al(randperm(4, 4));
      if rand < pT, c(1) = D.sE(l,k) + 0.003 * randn; end            % true lobe present with prob pT
      sE(l,k,1,:) = reshape(c, 1,1,1,4);  qE(l,k,1,:) = 0;           % angle only (no range)
    end, end
    xh = local_triangulate(sE, qE, NF, D.bet, FG);
    dx = abs(real(xh) - real(D.UEm)); dx = min(dx, sL - dx); dy = abs(imag(xh) - imag(D.UEm)); dy = min(dy, sL - dy);
    E = [E; sqrt(dx.^2 + dy.^2)];
  end
  fprintf('grating period %.2f, true lobe present w.p. %.1f: median error %.1f m, within 1 m %.0f%%\n', per, pT, median(E(isfinite(E))), 100*mean(E(isfinite(E)) < 1));
 end
end

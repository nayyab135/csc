#!/usr/bin/env python3
"""Fourth upgrade of Code T (in place): co-array comparison per parameter
 * ULA 4 lambda + co-array schemes (with / without CAPE), switch ulaw_coarray
 * two-user angular-resolution test (physical vs co-array MUSIC)
 * figures D1-D8 (one parameter per figure, one panel per array)"""
import pathlib, sys
REPO = pathlib.Path('/home/user/csc')
HERE = pathlib.Path(__file__).resolve().parent
figd = (HERE / 'figs_d.m').read_text()
ind = lambda x: ''.join(('    ' + ln if ln.strip() else ln) for ln in x.splitlines(True))

def rep(s, a, b, tag):
    if s.count(a) != 1:
        sys.exit(f'{tag}: anchor found {s.count(a)} times:\n{a[:250]}')
    return s.replace(a, b)

f = REPO / 'Code_T_ArrayComparison.m'
s = f.read_text()
s = rep(s, """%  FIGURES: A1-A4 array comparison, B1-B3 role of the difference co-array,
%  C1-C3 benefit of CAPE (T1, T2, T9 geometry; T3-T17 with fig_old = true).""",
"""%  FIGURES: A1-A4 array comparison, B1-B3 role of the difference co-array,
%  C1-C3 benefit of CAPE, D1-D8 the co-array parameter by parameter (one
%  panel per array: physical vs co-array processing of the SAME antennas)
%  (T1, T2, T9 geometry; T3-T17 with fig_old = true).""", 'header')
s = rep(s, "nSchE       = numel(en_names);\nen_caseName = {'Crowded (N = 8, K = 16)'};",
"""ulaw_coarray = true;          % also run the wide ULA (4 lambda) WITH co-array processing (with / without
                              % CAPE) so every sparse array has a co-array comparison (+ about 35% run time)
if ulaw_coarray
    en_names = [en_names, {sprintf('ULA %g lambda + co-array + CAPE', ulaw_d_lambda), ...
                           sprintf('ULA %g lambda + co-array, no CAPE', ulaw_d_lambda)}];
    en_leg   = [en_leg, {sprintf('ULA (%g\\\\lambda) + co-array', ulaw_d_lambda), ...
                         sprintf('ULA (%g\\\\lambda) + co-array, no CAPE', ulaw_d_lambda)}];
    en_geom  = [en_geom, {'ulaw', 'ulaw'}];
    en_idea  = [en_idea, 1 1];  en_cape = [en_cape, 1 0];  en_arr = [en_arr, 2 2];
end
nSchE       = numel(en_names);
en_caseName = {'Crowded (N = 8, K = 16)'};""", 'schemes')
s = rep(s, "DOA_rP = zeros(4, numel(DOA_snr));  DOA_rC = DOA_rP;  DOA_lim = zeros(4, 2);",
"""DOA_rP = zeros(4, numel(DOA_snr));  DOA_rC = DOA_rP;  DOA_lim = zeros(4, 2);
DOA_sepB = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 1 1.25 1.5];   % two-user separation [physical beamwidths lambda/aperture]
DOA_resP = zeros(4, numel(DOA_sepB));  DOA_resC = DOA_resP;   % P(both resolved): physical / co-array MUSIC""", 'doa init')
s = rep(s, """        doa_test(EN_geo{sc, 1}, EN_CA{sc, 1}, DOA_D0, 10.^(DOA_snr / 10), DOA_T, DOA_nTr, ca_Nsub_max, lambda);""",
"""        doa_test(EN_geo{sc, 1}, EN_CA{sc, 1}, DOA_D0, 10.^(DOA_snr / 10), DOA_T, DOA_nTr, ca_Nsub_max, lambda);
    [DOA_resP(a, :), DOA_resC(a, :)] = ...
        doa_resolution(EN_geo{sc, 1}, EN_CA{sc, 1}, DOA_sepB, 10^(DOA_snrD / 10), DOA_T, DOA_nTr, ca_Nsub_max, lambda);""", 'doa call')
s = rep(s, "         'DOA_pP', 'DOA_pC', 'DOA_rP', 'DOA_rC', 'DOA_lim');",
           "         'DOA_pP', 'DOA_pC', 'DOA_rP', 'DOA_rC', 'DOA_lim', 'DOA_sepB', 'DOA_resP', 'DOA_resC', 'ulaw_coarray');", 'save')
s = rep(s, "    % ---- console summary at 10 dB ----\n", ind(figd) + "\n    % ---- console summary at 10 dB ----\n", 'fig code T')
s = s.rstrip() + "\n" + r"""
function [pP, pC] = doa_resolution(G, CA, sepB, snr, T, nTr, nSubMax, lambda)
% [CODE T] Angular resolution of one array: two equal-power far-field users
% separated by sepB physical beamwidths (lambda / aperture, in sin(theta)),
% random centre inside the unambiguous sector; T snapshots. They are
% resolved if both are found within half their separation (one-to-one).
% Physical MUSIC vs difference-co-array MUSIC on the SAME snapshots.
    pos = G.pos(:);  N = numel(pos);
    per = min(2, lambda / G.unit);  w = min(per, 1.9);
    bw  = lambda / (max(pos) - min(pos));
    Ns  = min(CA.Nsub, nSubMax);
    pP = zeros(1, numel(sepB));  pC = pP;
    for ip = 1:numel(sepB)
        sep = sepB(ip) * bw;
        sg  = linspace(-w / 2, w / 2, min(20001, max(4001, ceil(40 * w / sep))));
        Ap  = exp(1j * 2 * pi * pos * sg / lambda);
        Av  = exp(1j * 2 * pi * (0:Ns - 1)' * CA.d0 * sg / lambda);
        nP = 0;  nC = 0;
        for tr = 1:nTr
            c = (rand - 0.5) * max(0.9 * w - sep, 0);
            s = c + [-sep / 2, sep / 2];
            A = exp(1j * 2 * pi * pos * s / lambda);
            Y = sqrt(snr / 2) * A * (randn(2, T) + 1j * randn(2, T)) + sqrt(0.5) * (randn(N, T) + 1j * randn(N, T));
            R = (Y * Y') / T;
            [U, Dg] = eig((R + R') / 2);  [~, o] = sort(real(diag(Dg)), 'descend');  U = U(:, o);
            Pm = 1 ./ max(sum(abs(U(:, 3:end)' * Ap) .^ 2, 1), eps);
            nP = nP + (match_peaks(s, coarray_peaks(sg, Pm, 2), sep / 2, per) == 2);
            z = R(:);  zu = accumarray(CA.ic, z) ./ max(accumarray(CA.ic, ones(size(z))), 1);
            zv = zu(CA.loc);  Rss = zeros(Ns);
            for i = 1:CA.Nv - Ns + 1
                zi = zv(i:i + Ns - 1);  Rss = Rss + zi * zi';
            end
            [U, Dg] = eig((Rss + Rss') / 2);  [~, o] = sort(real(diag(Dg)), 'descend');  U = U(:, o);
            Pc = 1 ./ max(sum(abs(U(:, 3:end)' * Av) .^ 2, 1), eps);
            nC = nC + (match_peaks(s, coarray_peaks(sg, Pc, 2), sep / 2, per) == 2);
        end
        pP(ip) = nP / nTr;  pC(ip) = nC / nTr;
    end
end
"""
f.write_text(s)

p = REPO / 'plot_codeT_figures.m'
t = p.read_text()
t = rep(t, "% ---- console summary at 10 dB ----\n", figd + "\n% ---- console summary at 10 dB ----\n", 'fig plot')
t = rep(t, """%  A1-A4: which sparse array performs better; B1-B3: role of the difference
%  co-array; C1-C3: benefit of CAPE; T1, T2, T9: geometry.""",
"""%  A1-A4: which sparse array performs better; B1-B3: role of the difference
%  co-array; C1-C3: benefit of CAPE; D1-D8: the co-array parameter by
%  parameter, one panel per array; T1, T2, T9: geometry.""", 'plot header')
p.write_text(t)
g = HERE / 'figblock_t3.m'
u = g.read_text()
g.write_text(rep(u, "% ---- console summary at 10 dB ----\n", figd + "\n% ---- console summary at 10 dB ----\n", 'figblock'))
print('ok')

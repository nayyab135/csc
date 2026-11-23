#!/usr/bin/env python3
"""Third upgrade (in place): fair no-CAPE localisation baseline (same J local
lobes per AP as CAPE) and empty-safe median/mean helpers in the figures."""
import pathlib, sys
REPO = pathlib.Path('/home/user/csc')
def rep(s, a, b, tag):
    if s.count(a) != 1:
        sys.exit(f'{tag}: anchor found {s.count(a)} times:\n{a[:200]}')
    return s.replace(a, b)
f = REPO / 'Code_T_ArrayComparison.m'
s = f.read_text()
s = rep(s, """    sEst = nan(L, K, nReal);  qEst = nan(L, K, nReal);   % local (angle, inverse range) of the first atom
    hasFG = nargin >= 17 && ~isempty(FG);""",
"""    hasFG = nargin >= 17 && ~isempty(FG);
    nJ = 1;  if hasFG && isfield(FG, 'J'), nJ = FG.J + 1; end
    sEst = nan(L, K, nReal, nJ);  qEst = nan(L, K, nReal, nJ);   % local lobes (angle, inverse range) per link""", 'init')
s = rep(s, """                    sEst(l, k, :) = reshape(sA, 1, 1, []);  qEst(l, k, :) = reshape(qA, 1, 1, []);""",
"""                    sEst(l, k, :, 1) = reshape(sA, 1, 1, []);  qEst(l, k, :, 1) = reshape(qA, 1, 1, []);
                    if nJ > 1          % [CODE T] the same J strongest local lobes CAPE starts from (baseline only)
                        cz = abs(PD.W' * z);
                        for m = 1:nReal
                            [~, oc] = sort(cz(:, m), 'descend');
                            for j = 1:nJ - 1
                                [~, ~, sEst(l, k, m, j + 1), qEst(l, k, m, j + 1)] = ...
                                    refine_atom(z(:, m), PD.s(oc(j)), PD.q(oc(j)), PD);
                            end
                        end
                    end""", 'store')
s = rep(s, """% [CODE T] Localisation WITHOUT CAPE (baseline): every one of the FG.nAP
% strongest near-field APs of a user makes its own hard decision (its local
% angle and inverse range), which gives candidate positions (two mirror
% solutions; a ray of ranges for a far-field atom). The CPU picks, among all
% candidates, the one most consistent with the other APs (smallest sum of
% distances to their nearest candidates) and averages it with those nearest
% candidates. Same APs and pilots as CAPE; only decisions are fused, not the
% pilot likelihoods.
    [~, K, M] = size(sE);""",
"""% [CODE T] Localisation WITHOUT CAPE (baseline): every one of the FG.nAP
% strongest near-field APs of a user makes its own hard decisions: its local
% estimate plus the same J strongest local lobes CAPE starts from (angle and
% inverse range), each giving candidate positions (two mirror solutions; a
% ray of ranges for a far-field atom). The CPU picks, among all candidates,
% the one most consistent with the other APs (smallest sum of distances to
% their nearest candidates) and averages it with those nearest candidates.
% Same APs, pilots and candidate lobes as CAPE: only the DECISIONS are fused
% (geometry), not the pilot likelihoods, and there is no refinement.
    [~, K, M, nJ] = size(sE);""", 'doc')
s = rep(s, """                if ~isfinite(sE(A(i), k, m)), ok(i) = false; continue; end
                P{i} = cand_positions(FG.AP(A(i)), sE(A(i), k, m), qE(A(i), k, m), sL, FG.hD);""",
"""                Pi = [];
                for j = 1:nJ
                    if ~isfinite(sE(A(i), k, m, j)), continue; end
                    Pi = [Pi; cand_positions(FG.AP(A(i)), sE(A(i), k, m, j), qE(A(i), k, m, j), sL, FG.hD)]; %#ok<AGROW>
                end
                if isempty(Pi), ok(i) = false; continue; end
                P{i} = Pi;""", 'cand')
f.write_text(s)
# empty-safe helpers in the figure code (Code T, plot script, figblock)
old_n = """tNan  = @(X) deal(arrayfun(@(i) mean(X(i, isfinite(X(i, :)))), (1:size(X, 1)).'), ...
                  arrayfun(@(i) std(X(i, isfinite(X(i, :)))) / sqrt(max(nnz(isfinite(X(i, :))), 1)), (1:size(X, 1)).'));
tMedN = @(X) arrayfun(@(i) median(X(i, isfinite(X(i, :)))), (1:size(X, 1)).');"""
new_n = """tFin  = @(v) [v(isfinite(v)), nan(1, ~any(isfinite(v)))];       % finite entries (NaN if none)
tNan  = @(X) deal(arrayfun(@(i) mean(tFin(X(i, :))), (1:size(X, 1)).'), ...
                  arrayfun(@(i) std(tFin(X(i, :))) / sqrt(max(nnz(isfinite(X(i, :))), 1)), (1:size(X, 1)).'));
tMedN = @(X) arrayfun(@(i) median(tFin(X(i, :))), (1:size(X, 1)).');"""
for p in ['Code_T_ArrayComparison.m', 'plot_codeT_figures.m', 'docs/codeT_build/figblock_t3.m', 'docs/codeT_build/figs_abc.m']:
    g = REPO / p; t = g.read_text()
    a, b = old_n, new_n
    if p == 'Code_T_ArrayComparison.m':
        a = ''.join('    ' + ln if ln.strip() else ln for ln in old_n.splitlines(True))
        b = ''.join('    ' + ln if ln.strip() else ln for ln in new_n.splitlines(True))
    g.write_text(rep(t, a, b, p))
print('ok')

#!/usr/bin/env python3
"""Figures for array_comparison_explained.pdf (Code T note)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = '/home/user/csc/docs/'
lam = 0.1
d0 = lam / 2
N0 = 32
COL = {'ula': '#595959', 'ulaw': '#007300', 'nested': '#001a8c', 'coprime': '#801a99'}
NAME = {'ula': r'ULA, $d=\lambda/2$', 'ulaw': r'ULA, $d=4\lambda$', 'nested': 'Nested', 'coprime': 'Coprime (M=8, N=17)'}


def geo(t):
    if t == 'ula':
        n = np.arange(N0)
    elif t == 'ulaw':
        n = np.arange(N0) * 8
    elif t == 'nested':
        n1 = N0 // 2
        n = np.unique(np.r_[np.arange(1, n1 + 1), np.arange(1, N0 - n1 + 1) * (n1 + 1)])
    else:
        n = np.unique(np.r_[np.arange(17) * 8, np.arange(16) * 17])
    p = n * d0
    return p - (p.max() + p.min()) / 2


def holefree(p):
    lg = np.unique(np.round((p[:, None] - p[None, :]) / d0).astype(int))
    br = np.r_[0, np.where(np.diff(lg) != 1)[0] + 1, len(lg)]
    i = np.argmax(np.diff(br))
    return lg[br[i]:br[i + 1]] * d0


types = ['ula', 'ulaw', 'nested', 'coprime']

# ---- Fig 1: layouts + hole-free co-arrays ----
fig, ax = plt.subplots(4, 1, figsize=(9, 6.2), sharex=True)
for a, t in zip(ax, types):
    p = geo(t)
    v = holefree(p)
    a.plot(v, 0.55 * np.ones_like(v), '|', color=COL[t], ms=7, alpha=0.6)
    a.plot(p, np.ones_like(p), 'o', color=COL[t], ms=4.5)
    a.set_ylim(0.3, 1.3)
    a.set_yticks([0.55, 1])
    a.set_yticklabels(['virtual\n(hole-free)', 'physical'], fontsize=7.5)
    ap = p.max() - p.min()
    a.set_title(f'{NAME[t]}:  32 antennas, aperture {ap:.2f} m, hole-free co-array {len(v)} lags',
                fontsize=9, loc='left')
    a.grid(alpha=0.3)
ax[-1].set_xlabel('position along the array [m]')
fig.tight_layout()
fig.savefig(OUT + 'codeT_fig_layouts.png', dpi=180)
plt.close(fig)

# ---- Fig 2: far-field beampatterns ----
fig, ax = plt.subplots(2, 2, figsize=(9, 5.2), sharey=True)
sg = np.linspace(-1, 1, 8001)
for a, t in zip(ax.ravel(), types):
    p = geo(t)
    B = np.abs(np.exp(1j * 2 * np.pi * np.outer(sg, p) / lam).sum(1)) / len(p)
    a.plot(sg, 20 * np.log10(np.maximum(B, 1e-4)), color=COL[t], lw=0.9)
    a.set_ylim(-40, 2)
    a.set_title(NAME[t], fontsize=9.5)
    a.grid(alpha=0.3)
    a.set_xlabel(r'$\sin\theta$', fontsize=9)
for a in ax[:, 0]:
    a.set_ylabel('beampattern [dB]', fontsize=9)
fig.tight_layout()
fig.savefig(OUT + 'codeT_fig_beams.png', dpi=180)
plt.close(fig)

# ---- Fig 3: same antennas, two processing chains ----
fig, ax = plt.subplots(figsize=(9.2, 3.6))
ax.set_xlim(0, 100)
ax.set_ylim(0, 40)
ax.axis('off')


def box(x, y, w, h, txt, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4', fc=fc, ec='#333333', lw=0.8))
    ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=7.6)


def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='->', lw=0.9, color='#333333'))


box(1, 15, 13, 10, 'SAME 32\nphysical antennas\n(nested array)', '#dfe6f7')
box(20, 27, 20, 9, 'pilot observation\nof each user\n(physical array)', '#eef2fb')
box(20, 3, 20, 10, 'co-array branch:\ncovariance -> vectorize ->\nvirtual array -> MUSIC', '#fbe9d9')
box(46, 27, 18, 9, 'candidate positions\n(local + CAPE)', '#eef2fb')
box(46, 3, 18, 10, 'EXTRA candidate\npositions\n(from co-array angles)', '#fbe9d9')
box(70, 15, 13, 10, 'SAME CAPE fused\nlikelihood picks\nthe best one', '#dff2df')
box(87, 15, 12, 10, 'channel\nestimate ->\ndetection', '#eeeeee')
arrow(14, 22, 20, 31)
arrow(14, 18, 20, 8)
arrow(40, 31, 46, 31)
arrow(40, 8, 46, 8)
arrow(64, 31, 70, 22)
arrow(64, 8, 70, 18)
arrow(83, 20, 87, 20)
ax.text(50, 39, '"Nested + CAPE" uses only the upper path.  "Nested + co-array + CAPE" adds the lower path.',
        ha='center', fontsize=8.2, style='italic')
fig.tight_layout()
fig.savefig(OUT + 'codeT_fig_processing.png', dpi=180)
plt.close(fig)
print('figures written')

#!/usr/bin/env python3
"""Build array_comparison_explained.pdf: answers to the professor's questions and Code T (crowded case)."""
import json
import pathlib
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image,
                                ListFlowable, ListItem, HRFlowable, KeepTogether)

OUT = '/home/user/csc/array_comparison_explained.pdf'
FIG = '/home/user/csc/docs/'
VERIF = pathlib.Path(FIG + 'codeT_verification.json')   # filled from the Octave check run

ss = getSampleStyleSheet()
NAVY = colors.HexColor('#12305c')
H1 = ParagraphStyle('H1', parent=ss['Heading1'], fontSize=14, spaceBefore=14, spaceAfter=6, textColor=NAVY,
                    keepWithNext=1)
H2 = ParagraphStyle('H2', parent=ss['Heading2'], fontSize=11.5, spaceBefore=10, spaceAfter=4,
                    textColor=colors.HexColor('#22548c'), keepWithNext=1)
BODY = ParagraphStyle('BODY', parent=ss['Normal'], fontSize=10, leading=14.2, alignment=TA_JUSTIFY, spaceAfter=6)
SMALL = ParagraphStyle('SMALL', parent=BODY, fontSize=8.6, leading=11.2, alignment=0, spaceAfter=0)
NOTE = ParagraphStyle('NOTE', parent=BODY, fontSize=9.4, leading=13, leftIndent=8, rightIndent=8,
                      textColor=colors.HexColor('#5a3a00'), backColor=colors.HexColor('#fff6e6'),
                      borderPadding=6, spaceBefore=6, spaceAfter=10)
KEY = ParagraphStyle('KEY', parent=BODY, fontSize=9.8, leading=13.6, leftIndent=8, rightIndent=8,
                     backColor=colors.HexColor('#eef3fb'), borderPadding=6, spaceBefore=4, spaceAfter=8)
CAP = ParagraphStyle('CAP', parent=BODY, fontSize=8.8, leading=11.8, textColor=colors.HexColor('#444444'),
                     alignment=TA_CENTER)
TITLE = ParagraphStyle('TITLE', parent=ss['Title'], fontSize=19, leading=24, textColor=NAVY)
SUB = ParagraphStyle('SUB', parent=ss['Normal'], fontSize=10.5, alignment=TA_CENTER,
                     textColor=colors.HexColor('#555555'))
REF = ParagraphStyle('REF', parent=BODY, fontSize=8.8, leading=12, leftIndent=16, firstLineIndent=-16,
                     spaceAfter=3, alignment=0)


def P(t, s=BODY):
    return Paragraph(t, s)


def bullets(items, s=BODY):
    return ListFlowable([ListItem(Paragraph(x, s), leftIndent=12) for x in items],
                        bulletType='bullet', start='•', bulletFontName='Helvetica', leftIndent=14)


def numbered(items, s=BODY):
    return ListFlowable([ListItem(Paragraph(x, s), leftIndent=14) for x in items],
                        bulletType='1', leftIndent=16)


def table(rows, widths, head=True):
    data = [[Paragraph(c, SMALL) for c in r] for r in rows]
    t = Table(data, colWidths=widths, repeatRows=1 if head else 0)
    st = [('GRID', (0, 0), (-1, -1), 0.4, colors.HexColor('#9aa7bd')),
          ('VALIGN', (0, 0), (-1, -1), 'TOP'),
          ('LEFTPADDING', (0, 0), (-1, -1), 4), ('RIGHTPADDING', (0, 0), (-1, -1), 4),
          ('TOPPADDING', (0, 0), (-1, -1), 3), ('BOTTOMPADDING', (0, 0), (-1, -1), 3)]
    if head:
        st += [('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#dbe4f3'))]
    t.setStyle(TableStyle(st))
    return t


def fig(path, w_cm, caption):
    img = Image(path)
    r = img.imageHeight / img.imageWidth
    img.drawWidth = w_cm * cm
    img.drawHeight = w_cm * cm * r
    return KeepTogether([img, Spacer(1, 3), P(caption, CAP), Spacer(1, 8)])


s = []
# ------------------------------------------------------------------ title
s += [Spacer(1, 0.6 * cm),
      P('Array Comparison for the Crowded Cell-Free Uplink', TITLE),
      Spacer(1, 4),
      P("Answers to the professor's questions, and what the new comparison code (Code T) does", SUB),
      Spacer(1, 10), HRFlowable(width='100%', thickness=1, color=NAVY), Spacer(1, 10)]

s += [P('<b>Purpose of this note.</b> The professor raised these points after looking at the Code R / Code S '
        'results: (1) the array we called a "uniform sparse array" should be called a uniform linear array; '
        '(2) the wide-spacing idea gives only a small benefit, so we should compare it properly against the '
        'normal ULA (spacing &lambda;/2) and against nested and coprime arrays; (3) the nested array and the '
        '"difference co-array" give almost the same results although they have the same antennas, so why are '
        'they treated as two things? (4) a cell-free system is a <b>loaded</b> system: few antennas per AP and '
        'many users, so we should study only the crowded case (8 antennas per AP, 16 users); and (5) BER curves '
        'must be proper waterfall curves, reaching about 10<super>-4</super> to 10<super>-5</super> around '
        '20 dB. This note answers each point in simple words, checked against our own code and the literature, '
        'and explains the new code <b>Code_T_ArrayComparison.m</b>.')]

s += [P('Short answers', H1)]
s += [table([
    ['Question', 'Short answer'],
    ['ULA with large equal spacing: "uniform linear array" or "uniform sparse array"?',
     '<b>Uniform linear array (ULA).</b> "Uniform linear" describes the geometry: antennas on a line with '
     'equal spacing. The spacing value does not change the family; large spacing only adds grating lobes. Say '
     '"a ULA with spacing d = 4&lambda;" or "a wide-spaced ULA".'],
    ['Why did the wide-spacing idea give only a small benefit?',
     'The large aperture improves resolution, but the spacing creates <b>grating lobes</b> (at d = 4&lambda;, 7 '
     'peaks as strong as the main beam). CAPE removes this ambiguity only partly.'],
    ['Nested vs "difference co-array": what is different?',
     '<b>Nothing in hardware.</b> Same antennas at the same positions. The difference co-array is a '
     '<i>virtual</i> array computed from the physical one; only the <b>processing</b> differs.'],
    ['Why are their results almost the same?',
     'The co-array only proposes a few <b>extra candidate positions</b>; the same CAPE likelihood on the '
     'physical antennas still picks the answer. More reasons in Section 3.'],
    ['Crowded system with 8 antennas: does the near field disappear?',
     'Only for the dense ULA. A sparse array has a larger aperture, so a larger Rayleigh distance. With 8 '
     'antennas, nested/coprime on the &lambda;/2 grid reach only 11-18 m, so we <b>stretch</b> them to the same '
     '~2.8 m aperture as the 4&lambda; ULA: then all three sparse arrays reach 160-180 m.'],
    ['Why were the BER curves not waterfall curves?',
     'Our "SNR" axis was really the user <b>transmit power in dBm</b>, and users had very different link '
     'qualities (22 dB spread). Code T fixes both with power control: the proposed detector now waterfalls. '
     'The simpler detectors still floor because each AP has 8 antennas for 16 users.'],
], [5.0 * cm, 11.6 * cm])]

# ------------------------------------------------------------------ 1. naming
s += [P('1. Naming: uniform linear array, not "uniform sparse array"', H1)]
s += [P('An array is named after its <b>geometry</b>, meaning where the antennas are. A <b>uniform linear '
        'array (ULA)</b> is N antennas on a straight line with the same spacing d between neighbours. The '
        'definition does not limit d. Textbooks (e.g. Van Trees, <i>Optimum Array Processing</i> [10]) explain '
        'that d &le; &lambda;/2 avoids grating lobes, and that a ULA with d &gt; &lambda;/2 is still a ULA. It '
        'is simply under-sampled in space, so it has grating lobes.')]
s += [P('The literature uses the ULA name in exactly this way:')]
s += [bullets([
    'Vaidyanathan and Pal build the <b>coprime array</b> from two <b>uniform linear arrays</b> whose spacings are '
    'M&middot;d and N&middot;d, much larger than &lambda;/2. Later papers describe a coprime array as "two ULAs '
    'that are undersampled by coprime factors" [2], [3].',
    'In near-field communications, Zhou <i>et al.</i> [6] study a uniform array with spacing larger than '
    '&lambda;/2 and call it a <b>linear sparse array (LSA)</b>, defined as a ULA with enlarged spacing. It '
    'improves near-field beam focusing but creates grating lobes as strong as the main lobe.',
    'Modular / widely spaced array papers [7], [8] describe arrays built from widely spaced <b>ULAs</b> and report '
    'the same trade-off: larger aperture and better resolution, but grating lobes that reduce the SINR of other '
    'users.',
])]
s += [P('In the sparse-array literature, "sparse array" normally means a <b>non-uniform</b> array (nested, '
        'coprime, MISC, EMISC, ...) whose positions are chosen so that the difference co-array is large and has '
        'no holes. A widely spaced ULA is "sparse" only in the loose sense of "few antennas over a long aperture".')]
s += [P('<b>Recommended wording for the paper:</b> "a uniform linear array (ULA) with inter-antenna spacing '
        'd = 4&lambda;", or "a wide-spaced ULA". Use "sparse ULA" only as a description, never as a new array '
        'type.', KEY)]

# ------------------------------------------------------------------ 2. why small gain
s += [P('2. Why the wide-spaced ULA gave only a small benefit', H1)]
s += [P('Widening the spacing from &lambda;/2 to 4&lambda; at the same number of antennas multiplies the aperture '
        'by eight. A larger aperture gives finer angle resolution and a much larger Rayleigh distance, so more '
        'users are in the near field, where range can also be resolved. The cost is the <b>grating lobes</b>: '
        'with d = 4&lambda; the array cannot tell apart directions whose sin&theta; differ by a multiple of '
        '&lambda;/d = 1/4. Figure 1 shows the far-field beampattern of the four 8-antenna arrays used in Code T: '
        'the 4&lambda; ULA has <b>7 peaks at full gain</b>.')]
s += [fig(FIG + 'codeT_fig_beams.png', 15.5,
          'Figure 1. Far-field beampattern of the four 8-antenna arrays of Code T. 4&lambda; ULA: 7 full-gain '
          'peaks. Nested and coprime, stretched to the same aperture: 3 full-gain peaks each. Dense ULA: one main '
          'beam, -13 dB sidelobes, but a 0.35 m aperture.')]
s += [P('What this meant in our earlier results (Code R, 32 antennas, spacing 4.5&lambda;, 6 random setups):')]
s += [bullets([
    '<b>Detection improved a little with many users.</b> At K = 40 and 10 dB: BER 0.100 vs 0.142 for nested, sum '
    'SE 145 vs 118 bit/s/Hz. With K = 8 the two were tied within the noise; with 6 setups the error bars were '
    'wide (about &plusmn;20 bit/s/Hz on the sum SE).',
    '<b>Localization failed in some setups</b> (median errors of 1.5 m to 40 m against centimetres for nested): a '
    'grating lobe points to a wrong position and CAPE did not always remove it.',
    '<b>Why detection survives but localization does not:</b> a grating-lobe ghost still gives a steering vector '
    'that matches the received signal well, so the channel estimate is usable for detection even when the '
    'position is wrong.',
])]
s += [P('<b>One real advantage our simulation cannot show.</b> Wide spacing removes closely spaced antennas, so '
        'mutual coupling is low (Table 2). But our channel model has <b>no mutual-coupling model</b>, so this '
        'benefit can be stated as a hardware property only; it will never appear in our BER curves.', NOTE)]

# ------------------------------------------------------------------ 3. nested vs co-array
s += [P('3. Nested array vs "difference co-array"', H1)]
s += [P('3.1 What each one is', H2)]
s += [P('The <b>physical array</b> is the set of antennas you build. The <b>difference co-array</b> is the set of '
        'all differences x<sub>m</sub> - x<sub>n</sub> between antenna positions. It appears because, for a '
        'far-field source, the covariance entry between antennas m and n depends only on that difference: '
        'R<sub>mn</sub> ~ exp(j2&pi;(x<sub>m</sub>-x<sub>n</sub>)sin&theta;/&lambda;). Rearranging the covariance '
        'entries by lag gives samples of a longer <b>virtual</b> array, on which methods such as co-array MUSIC '
        'work.')]
s += [bullets([
    '<b>Every array has a difference co-array</b>, even a plain ULA (8 antennas give lags -7 ... 7, i.e. 15 '
    'values, so nothing is gained).',
    'Nested and coprime arrays are <b>designed</b> so that the co-array is long and hole-free: with 8 antennas, '
    '39 hole-free lags for nested and 23 for coprime (Figure 2, Table 2).',
    'So "nested array" names a <b>physical geometry</b> and "difference co-array processing" names a '
    '<b>signal-processing method</b> applied to it. They are not two arrays.',
])]
s += [fig(FIG + 'codeT_fig_layouts.png', 15.5,
          'Figure 2. The four 8-antenna arrays of Code T (dots) and their hole-free difference co-arrays (ticks), '
          'on each array\'s own grid. The 4&lambda; ULA co-array is only a 15-lag ULA with 4&lambda; spacing.')]
s += [P('3.2 What our code actually does', H2)]
s += [P('In our code, "Nested + CAPE" and "Nested + co-array + CAPE" have the same geometry, the same antennas and '
        'the same CAPE estimator. The <b>only</b> difference is the flag <font face="Courier">en_idea</font>. When '
        'it is on, each AP forms the sample covariance of its pilot samples, builds the virtual array, runs MUSIC, '
        'and turns the MUSIC angles (plus a range search) into <b>extra candidate positions</b>. The final '
        'position is still chosen by the <b>same</b> CAPE fused likelihood on the <b>physical</b> antennas '
        '(Figure 3).')]
s += [fig(FIG + 'codeT_fig_processing.png', 15.5,
          'Figure 3. Same antennas, two processing chains. The co-array branch can only add candidates; the '
          'decision is always made by the same CAPE likelihood.')]
s += [P('3.3 Why the two curves were almost identical', H2)]
s += [numbered([
    '<b>The co-array can only add candidates.</b> If CAPE\'s own candidates already contain the right position, '
    'the extra candidates never win and the result is exactly the same.',
    '<b>We estimate one user at a time.</b> With orthogonal pilots each AP estimates one channel at a time. The '
    'famous co-array advantage, resolving <i>more sources than antennas</i>, matters when many sources must be '
    'separated at once from one covariance, as in DOA estimation [1], [4], [5]. Note that the crowded case '
    '(K = 16 &gt; N = 8) is exactly that regime, so if the co-array helps anywhere, it is here.',
    '<b>The co-array model assumes plane waves</b>, but near-field links have curved wavefronts, so the co-array '
    'angles are only approximate.',
    '<b>Very few samples for the covariance:</b> &tau; = K pilot samples.',
    '<b>Code R capped the virtual aperture</b> at 80 elements (3.95 m), smaller than the physical aperture. Even '
    'without the cap, spatial smoothing keeps about half of the co-array, i.e. about the physical aperture: '
    'the co-array adds degrees of freedom, not resolution beyond the physical aperture.',
    '<b>Code R comparisons were not paired:</b> every scheme had a different random seed, so the small '
    'differences between nested and co-array were Monte Carlo noise.',
])]
s += [P('<b>For the paper:</b> keep arrays and processing separate: "nested array (direct processing)" and "nested '
        'array with difference co-array processing". If the co-array gives no measurable gain, say so in one '
        'sentence with the reasons above. The professor\'s remark matches the literature: papers pick a sparse '
        'array and then apply co-array processing to it; nobody lists "nested" and "difference co-array" as two '
        'competing arrays.', KEY)]

# ------------------------------------------------------------------ 4. crowded system
s += [P('4. The crowded cell-free system (8 antennas per AP, 16 users)', H1)]
s += [P('The professor\'s point is that cell-free networks are <b>loaded</b>: many cheap APs with few antennas each, '
        'serving many users. Code T therefore keeps only the crowded case: L = 8 APs, N = 8 antennas per AP, '
        'K = 16 users (64 antennas in total for 16 users, and 16 users against 8 antennas at every AP).')]
s += [P('4.1 Does the near field disappear with 8 antennas?', H2)]
s += [P('The Rayleigh distance is 2D<super>2</super>/&lambda;, where D is the aperture. The near field matters only '
        'for users closer than this. In our layout every AP-user distance is at least about 22 m (20 m horizontal '
        'plus the 10 m AP height). With 8 antennas:')]
s += [table([
    ['Array (8 antennas)', 'Grid', 'Aperture', 'Rayleigh distance', 'Near-field links (measured)'],
    ['ULA &lambda;/2', '&lambda;/2', '0.35 m', '2.5 m', '0 %'],
    ['Coprime (2,5), classic', '&lambda;/2', '0.75 m', '11 m', '0 %'],
    ['Nested 4+4, classic', '&lambda;/2', '0.95 m', '18 m', '0-2 %'],
    ['ULA 4&lambda;', '4&lambda;', '2.80 m', '157 m', '76-80 %'],
    ['Nested 4+4, <b>stretched</b>', '1.5&lambda;', '2.85 m', '162 m', 'like the 4&lambda; ULA'],
    ['Coprime (2,5), <b>stretched</b>', '2&lambda;', '3.00 m', '180 m', 'like the 4&lambda; ULA'],
], [4.3 * cm, 1.5 * cm, 2.0 * cm, 3.0 * cm, 5.8 * cm])]
s += [Spacer(1, 6)]
s += [P('So your point is correct: a sparse array increases the Rayleigh distance. But the classic nested and '
        'coprime arrays on the &lambda;/2 grid are too short with only 8 antennas (11-18 m); their links are all '
        'far field, and CAPE (which only acts on near-field links) would do nothing for them. Code T therefore '
        '<b>stretches</b> the nested and coprime arrays to the same ~2.8 m aperture as the 4&lambda; ULA. Then '
        'all three sparse arrays have Rayleigh distances of 160-180 m, and the comparison is fair: <b>same '
        'number of antennas and same aperture</b>; only the layout differs.')]
s += [P('<b>The price of stretching.</b> If all antennas lie on a grid with spacing u &gt; &lambda;/2, the '
        'beampattern repeats with period &lambda;/u, so every stretched array has some grating lobes. But the '
        'number differs: 7 full-gain peaks for the 4&lambda; ULA, only 3 for stretched nested (1.5&lambda; grid) '
        'and stretched coprime (2&lambda; grid), see Figure 1. At equal aperture, the non-uniform layouts are '
        'clearly <b>less ambiguous</b> than the wide ULA. This is a real, explainable advantage of nested and '
        'coprime, and resolving the remaining ambiguity across APs is exactly what CAPE is for.', KEY)]

s += [P('4.2 The SNR axis: what our old plots really showed', H2)]
s += [P('In all earlier codes the path gain &beta; is already divided by the noise power (in mW), and the code '
        'sets the transmit power p equal to the "SNR" value. So the old x-axis "SNR [dB]" was really the '
        '<b>user transmit power in dBm</b>. In our diagnostic drop, the strongest-AP gain per user was -7.6 dB '
        'for the median user, -19 dB for the weakest 10 % and +2.6 dB for the best 10 %: a 22 dB spread. The '
        'average BER is dominated by the weak users, so the BER falls slowly instead of a waterfall.')]
s += [P('<b>What Code T does: uplink power control.</b> Every user scales its transmit power so that its '
        'received SNR per antenna, at its own strongest AP, equals the value on the x-axis. In the code this is '
        'one line, <font face="Courier">bet = bet ./ max(bet, [], 1)</font>, which is exactly equivalent. The '
        'x-axis is then "received SNR per antenna at the user\'s strongest AP". This is a common, honest '
        'definition in cell-free papers, and it must be stated in the paper. Weak users then transmit more '
        'power and also interfere more at other APs; that is included in the simulation.')]

s += [P('4.3 Why the old detectors had a BER floor', H2)]
s += [P('We ran a diagnostic of the crowded case in GNU Octave before changing anything (2 random drops, '
        '16-QAM, the old 3-dB clusters, classic &lambda;/2-grid nested/coprime). BER of the estimated-CSI '
        'detectors:')]
s += [table([
    ['Old transmit-power axis', '10 dBm', '20 dBm', '30 dBm', '40 dBm'],
    ['Linear MMSE, <b>perfect CSI</b> (all arrays)', '0.26-0.29', '0.23-0.24', '0.23-0.25', '0.21-0.23'],
    ['List SIC, <b>perfect CSI</b> (all arrays)', '0.20-0.24', '0.12-0.16', '0.10-0.15', '0.06-0.18'],
    ['Cross-AP List, ULA &lambda;/2', '2.0e-1', '1.1e-1', '2.3e-2', '1.6e-2'],
    ['Cross-AP List, ULA 4&lambda;', '1.8e-1', '3.5e-2', '1.4e-3', '7.8e-4'],
    ['Cross-AP List, nested (classic)', '2.0e-1', '7.3e-2', '4.3e-3', '6.3e-4'],
], [6.4 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm, 2.5 * cm])]
s += [Spacer(1, 6)]
s += [bullets([
    '<b>Linear, Hard-SIC and List-SIC stay at BER 0.1-0.25 even with perfect CSI and 40 dBm.</b> That is an '
    'interference floor, not a noise or estimation problem. Cause: these detectors first combine locally at '
    'each AP (only 8 antennas against 16 users) and the user clusters kept only APs within <b>3 dB</b> of the '
    'strongest AP, so each user was served by only 1.2-1.4 APs.',
    'Only our <b>Cross-AP List</b> detector kept falling, because it fuses all APs in its hardest stages. But '
    'it fell slowly, because of the 22 dB user spread on the transmit-power axis.',
    'With power control alone (still 16-QAM, still 3-dB clusters) the Cross-AP List curves became waterfalls '
    '(e.g. 4&lambda; ULA: 4.5e-2 at 10 dB, 7.8e-4 at 20 dB, no errors at 25 dB), while the other three '
    'detectors still floored at 0.13-0.24. This confirmed the two separate causes.',
])]
s += [P('4.4 What Code T changes for the crowded case (your choices)', H2)]
s += [numbered([
    '<b>Crowded case only:</b> N = 8, K = 16; the user sweep covers K = 8, 12, 16, 20, 24 (all K &ge; N).',
    '<b>Power control / SNR axis</b> as in Section 4.2; SNR from 0 to 28 dB in 4 dB steps.',
    '<b>16-QAM is kept.</b> QPSK was tested too, but with power control it reaches zero errors by 4-8 dB, so no '
    'waterfall is visible in the 0-20 dB range (Section 5.3).',
    '<b>Larger user clusters:</b> APs within 10 dB of the strongest AP (<font face="Courier">clu_div = 10</font>, '
    'was 2 = 3 dB), to lower the interference floor of the baseline detectors.',
    '<b>Equal-aperture sparse arrays:</b> nested and coprime stretched to ~2.8 m (Section 4.1).',
    '<b>Correct pilot overhead:</b> the sum rate now uses &tau;<sub>p</sub> = K pilots '
    '(prelog 1 - K/&tau;<sub>c</sub>). The old code always used 8 pilots, which overstated the sum rate by about '
    '4 % for K = 16 (about 20 % for K = 40 in Code R).',
    '<b>Enough bits for 10<super>-5</super>:</b> 10 drops x 4 realisations x 800 symbols x 16 users x 4 bits = '
    '2.0 x 10<super>6</super> bits per point. Fewer bits cannot show BER 10<super>-5</super> reliably.',
])]

# ------------------------------------------------------------------ 5. Code T
s += [P('5. The new code: Code_T_ArrayComparison.m', H1)]
s += [P('Code T is built from Code R (through Code S, which ran correctly on your machine). CAPE, the four '
        'detectors, the scenario generator and L = 8 APs are unchanged.')]
s += [P('<b>Table 1.</b> The six schemes (8 physical antennas per AP, all with CAPE).', BODY)]
s += [table([
    ['#', 'Scheme', 'Physical array', 'Processing'],
    ['1', 'ULA (&lambda;/2)', 'dense ULA, aperture 0.35 m', 'direct'],
    ['2', 'ULA (4&lambda;)', 'wide-spaced ULA, aperture 2.80 m', 'direct'],
    ['3', 'Nested (stretched)', 'nested 4+4 on a 1.5&lambda; grid, aperture 2.85 m', 'direct'],
    ['4', 'Nested (stretched) + co-array', '<b>same antennas as #3</b>', 'direct + co-array MUSIC candidates'],
    ['5', 'Coprime (stretched)', 'coprime (2,5) on a 2&lambda; grid, aperture 3.00 m', 'direct'],
    ['6', 'Coprime (stretched) + co-array', '<b>same antennas as #5</b>', 'direct + co-array MUSIC candidates'],
], [0.8 * cm, 4.2 * cm, 6.4 * cm, 5.2 * cm])]
s += [Spacer(1, 6)]
s += [P('<b>Table 2.</b> Geometry metrics of the four 8-antenna arrays (exact positions used in the code).', BODY)]
s += [table([
    ['Array', 'Aperture', 'Rayleigh', 'Hole-free co-array', 'w(1)', 'Coupling leakage', 'Full-gain peaks'],
    ['ULA &lambda;/2', '0.35 m', '2.5 m', '15', '7', '0.839', '1 (sidelobes -13 dB)'],
    ['ULA 4&lambda;', '2.80 m', '157 m', '15 (4&lambda; grid)', '0', '0.189', '7'],
    ['Nested, 1.5&lambda; grid', '2.85 m', '162 m', '39', '0', '0.358', '3'],
    ['Coprime, 2&lambda; grid', '3.00 m', '180 m', '23', '0', '0.243', '3'],
], [3.1 * cm, 1.8 * cm, 1.8 * cm, 2.9 * cm, 1.1 * cm, 2.6 * cm, 3.3 * cm])]
s += [P('w(1) = number of antenna pairs only &lambda;/2 apart (they couple most strongly). Coupling leakage uses '
        'the standard banded model. Co-array lengths are in units of each array\'s own grid.', SMALL),
      Spacer(1, 8)]
s += [P('5.1 Other changes compared with Code R', H2)]
s += [numbered([
    '<b>Paired comparison.</b> All schemes use the same random seed in a given drop and SNR, so they see the same '
    'users, pilot noise and data symbols; nested vs nested + co-array also see identical channels. The '
    'co-array covariance noise comes from its own random stream, so it does not shift the other random numbers. '
    'Any difference between #3 and #4 (or #5 and #6) is caused only by the co-array processing.',
    '<b>Co-array on the array\'s own grid and with no cap</b> (<font face="Courier">ca_Nsub_max = 1000</font>), so '
    'the co-array gets its best chance. On a stretched grid the virtual array spacing is also larger than '
    '&lambda;/2, so co-array MUSIC returns aliased angles; these are only candidates and CAPE decides.',
    '<b>Results are saved before plotting</b> (<font face="Courier">codeT_arraycomp_results.mat</font>, including '
    'the geometries). <font face="Courier">plot_codeT_figures.m</font> redraws all figures without rerunning.',
])]
s += [P('5.2 Figures', H2)]
s += [table([
    ['Fig.', 'Content', 'Question it answers'],
    ['T1', 'Antenna positions and hole-free co-array of each array', 'What is physically different?'],
    ['T2', 'Far-field beampattern', 'How many grating lobes does each array have?'],
    ['T3', 'Channel NMSE vs SNR', 'Which array gives the best channel estimate?'],
    ['T4', 'BER vs SNR (Cross-AP List)', 'Which array gives the best detection? (waterfall)'],
    ['T5', 'Sum SE vs SNR', 'Which array gives the best rate?'],
    ['T6', 'BER and sum SE vs number of users K at 10 dB', 'How does each array cope with more load?'],
    ['T7', 'Paired change caused by co-array processing', "The professor's co-array question"],
    ['T8', 'Median localization error vs SNR', 'Do grating lobes break localization?'],
    ['T9', 'Aperture, co-array, w(1), coupling leakage', 'Hardware cost and coupling'],
    ['T10', 'BER of all four detectors for each array', 'Do the baseline detectors still have a floor?'],
], [1.1 * cm, 7.2 * cm, 8.3 * cm])]
s += [Spacer(1, 6)]

v = json.loads(VERIF.read_text()) if VERIF.exists() else None
if v:
    s += [P('5.3 Verification run (before giving you the code)', H2)]
    s += [P(v['intro'])]
    s += [table(v['table'], v['widths'])]
    s += [Spacer(1, 4)]
    s += [P(v['comment'], NOTE)]

# ------------------------------------------------------------------ 6. how to judge
s += [P('6. How to read the results and decide', H1)]
s += [numbered([
    '<b>Use the error bars.</b> If two schemes differ by less than about two standard errors (printed in the '
    'console), treat them as equal.',
    '<b>Co-array (T7).</b> If the curves stay at 0 % and the title says the results were identical in most '
    'points, the co-array changes no decision in our system: say so in one sentence. A clear gain anywhere '
    '(most likely at low SNR, where it is the K &gt; N regime) is a real result and can be reported.',
    '<b>ULA 4&lambda; vs stretched nested vs stretched coprime (T3-T6).</b> Same antennas, same aperture: this is '
    'a clean test of <b>layout</b>. The expectation from Figure 1 is that nested and coprime (3 peaks) beat the '
    'ULA (7 peaks) in localization (T8), and possibly in BER.',
    '<b>ULA &lambda;/2 vs the sparse arrays.</b> The dense ULA has an 8-times smaller aperture and only far-field '
    'links, so CAPE cannot help it. A large gap here is mainly an aperture effect.',
    '<b>Detectors (T10).</b> In the verification run, Linear MMSE and List SIC still floor at about 0.1 with '
    '16-QAM. That is an honest message: in a loaded cell-free system local-first detectors are interference-limited, and '
    'the Cross-AP List detector is what makes the waterfall possible.',
])]

# ------------------------------------------------------------------ 7. limitations
s += [P('7. Limitations to keep in mind', H1)]
s += [bullets([
    '<b>EMISC is not included</b> (the professor allowed "coprime or EMISC"); it can be added as one more case in '
    '<font face="Courier">geom_targetN</font>.',
    '<b>No mutual coupling</b> in the channel model (Section 2).',
    '<b>Power control is channel inversion at the strongest AP.</b> It is simple and standard, but weak users '
    'transmit much more power. State it in the paper; other power-control rules would shift the curves.',
    '<b>Co-array processing uses a plane-wave virtual array</b>, as in the standard literature.',
    '<b>Run time.</b> The first version of Code T was far too slow (days). It now uses a vectorised detector '
    '(the same decisions, checked bit by bit against the old one), no unused perfect-CSI pass, a faster CAPE '
    'refinement and one parallel task per drop, scheme and SNR. One full-size task takes about 10 s in Octave '
    '(about 60 times faster); the script prints its own time estimate after the first task.',
])]

# ------------------------------------------------------------------ refs
s += [P('References', H1)]
refs = [
    '[1] P. Pal and P. P. Vaidyanathan, "Nested arrays: A novel approach to array processing with enhanced degrees '
    'of freedom," <i>IEEE Trans. Signal Process.</i>, vol. 58, no. 8, 2010.',
    '[2] P. P. Vaidyanathan and P. Pal, "Sparse sensing with co-prime samplers and arrays," <i>IEEE Trans. Signal '
    'Process.</i>, vol. 59, no. 2, 2011.',
    '[3] P. Pal and P. P. Vaidyanathan, "Coprime sampling and the MUSIC algorithm," <i>Proc. IEEE DSP/SPE '
    'Workshop</i>, 2011.',
    '[4] C.-L. Liu and P. P. Vaidyanathan, "Remarks on the spatial smoothing step in coarray MUSIC," <i>IEEE Signal '
    'Process. Lett.</i>, 2015; C.-L. Liu, PhD thesis, Caltech, 2018.',
    '[5] M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cram&eacute;r-Rao bound," <i>IEEE Trans. Signal Process.</i>, '
    '2017, arXiv:1605.03620.',
    '[6] C. Zhou, C. You, H. Zhang, L. Chen, S. Shi, "Sparse array enabled near-field communications: Beam pattern '
    'analysis and hybrid beamforming design," arXiv:2401.05690, 2024.',
    '[7] A. Kosasih, &Ouml;. T. Demir, E. Bj&ouml;rnson, "Near-field beamfocusing, localization, and channel estimation '
    'with modular linear arrays," arXiv:2505.07991, 2025.',
    '[8] Li <i>et al.</i>, modular XL-MIMO near-field analysis, arXiv:2308.11289, 2023.',
    '[9] Z. Zheng, W.-Q. Wang, Y. Kong, Y. D. Zhang, "MISC array: A new sparse array design achieving increased '
    'degrees of freedom and reduced mutual coupling effect," <i>IEEE Trans. Signal Process.</i>, 2019.',
    '[10] H. L. Van Trees, <i>Optimum Array Processing</i>, Wiley, 2002.',
    '[11] Wang <i>et al.</i>, "Exploring the advantages of sparse arrays in near-field XL-MIMO systems: Beam analysis '
    'and EDoF function," arXiv:2501.09234, 2025.',
]
s += [P(r, REF) for r in refs]


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(colors.HexColor('#777777'))
    canvas.drawString(2 * cm, 1.2 * cm, 'Array comparison note - Code T (crowded case)')
    canvas.drawRightString(19 * cm, 1.2 * cm, 'page %d' % doc.page)
    canvas.restoreState()


doc = SimpleDocTemplate(OUT, pagesize=A4, leftMargin=2.1 * cm, rightMargin=2.1 * cm,
                        topMargin=1.8 * cm, bottomMargin=2 * cm,
                        title='Array comparison: answers and Code T', author='')
doc.build(s, onFirstPage=footer, onLaterPages=footer)
print('wrote', OUT, '(verification section included)' if v else '(no verification section yet)')

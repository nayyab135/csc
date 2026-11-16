#!/usr/bin/env python3
"""Build array_comparison_explained.pdf: answers to the professor's questions and Code T description."""
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image,
                                ListFlowable, ListItem, HRFlowable, KeepTogether, PageBreak)

OUT = '/home/user/csc/array_comparison_explained.pdf'
FIG = '/home/user/csc/docs/'

ss = getSampleStyleSheet()
NAVY = colors.HexColor('#12305c')
H1 = ParagraphStyle('H1', parent=ss['Heading1'], fontSize=14, spaceBefore=14, spaceAfter=6, textColor=NAVY, keepWithNext=1)
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
                        bulletType='bullet', start='\u2022', bulletFontName='Helvetica', leftIndent=14)


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
      P('Array Comparison for the Cell-Free Near-Field Uplink', TITLE),
      Spacer(1, 4),
      P("Answers to the professor's questions, and what the new comparison code (Code T) does", SUB),
      Spacer(1, 10), HRFlowable(width='100%', thickness=1, color=NAVY), Spacer(1, 10)]

s += [P('<b>Purpose of this note.</b> The professor raised three points after looking at the Code R / Code S '
        'results: (1) the array we called a "uniform sparse array" should be called a uniform linear array; '
        '(2) the wide-spacing idea gives only a small benefit, so we should compare it properly against the '
        'normal ULA (spacing &lambda;/2) and against well-known sparse arrays (nested, coprime or EMISC); '
        '(3) the nested array and the "difference co-array" give almost the same results, and since they '
        'have the same physical antennas, why are they treated as two different things? This note answers '
        'each point in simple words, checked against our own code and the literature, and explains the new '
        'code <b>Code_T_ArrayComparison.m</b> that runs the comparison the professor asked for.')]

s += [P('Short answers', H1)]
s += [table([
    ['Question', 'Short answer'],
    ['Is a ULA with large equal spacing a "uniform linear array" or a "uniform sparse array"?',
     '<b>It is a uniform linear array (ULA).</b> "Uniform linear" describes the geometry: antennas on a '
     'line with equal spacing. The spacing value does not change the family. Large spacing only adds '
     'grating lobes. Say "a ULA with spacing d = 4&lambda;" or "a wide-spaced ULA". The professor is right.'],
    ['Why did the wide-spacing idea give only a small benefit?',
     'The large aperture improves resolution, but the spacing creates <b>grating lobes</b>: at d = 4&lambda; '
     'the beampattern has 7 peaks as strong as the main beam. CAPE removes this ambiguity only partly. '
     'Detection tolerates a wrong lobe (small BER / sum-rate gain), localization does not.'],
    ['Nested array vs "difference co-array": what is different?',
     '<b>Nothing in hardware.</b> In our code both use the same 32 antennas at the same positions. The '
     'difference co-array is not an array you build; it is a <i>virtual</i> array computed from the '
     'physical array. The only difference is the <b>processing</b>.'],
    ['Why are their results almost the same?',
     'In our code the co-array only proposes a few <b>extra candidate positions</b>; the same CAPE '
     'likelihood on the physical antennas still picks the answer. Also, the estimation is done one user at a '
     'time, the co-array model assumes plane waves (our links are near field), the covariance comes from '
     'very few pilot samples, and the co-array MUSIC was limited to a 3.95 m virtual aperture, smaller than '
     'the 13.55 m physical one.'],
    ['What does the new code do?',
     'Six schemes at equal hardware (32 antennas, all with CAPE): ULA &lambda;/2, ULA 4&lambda;, nested, '
     'nested + co-array, coprime, coprime + co-array. All schemes see the same random channels and noise '
     '(paired comparison), so any nested vs nested + co-array gap is caused only by the co-array processing.'],
], [5.0 * cm, 11.6 * cm])]

# ------------------------------------------------------------------ 1. naming
s += [P('1. Naming: uniform linear array, not "uniform sparse array"', H1)]
s += [P('An array is named after its <b>geometry</b>, meaning where the antennas are. A <b>uniform linear '
        'array (ULA)</b> is N antennas on a straight line with the same spacing d between neighbours. This '
        'definition does not put any limit on d. Textbooks (e.g. Van Trees, <i>Optimum Array Processing</i>) '
        'explain that d &le; &lambda;/2 avoids grating lobes, and that a ULA with d &gt; &lambda;/2 is still a '
        'ULA. It is simply under-sampled in space, so it has grating lobes.')]
s += [P('The literature uses the ULA name in exactly this way:')]
s += [bullets([
    'Vaidyanathan and Pal build the <b>coprime array</b> from two <b>uniform linear arrays</b> whose spacings are '
    'M&middot;d and N&middot;d, i.e. much larger than &lambda;/2. They are still called ULAs. Later papers describe a '
    'coprime array as "two ULAs that are undersampled by coprime factors" [2], [3].',
    'In near-field communications, Zhou <i>et al.</i> [6] study a uniform array with spacing larger than &lambda;/2 '
    'and call it a <b>linear sparse array (LSA)</b>. Their definition is a ULA with enlarged spacing. They show '
    'that it improves near-field beam focusing but creates grating lobes with power comparable to the main lobe.',
    'Modular / widely-spaced array papers (Kosasih, Demir and Bj&ouml;rnson [7]; Li <i>et al.</i> [8]) describe '
    'arrays built from widely spaced <b>ULAs</b>, and report the same trade-off: larger aperture, better '
    'resolution, but grating lobes that reduce the SINR of other users.',
])]
s += [P('In the sparse-array signal-processing literature, the word "sparse array" normally means a '
        '<b>non-uniform</b> array (nested, coprime, MISC, EMISC, ...) whose element positions are chosen so that '
        'the difference co-array is large and has no holes. A widely spaced ULA is "sparse" only in the loose '
        'sense of "few antennas over a long aperture". Its difference co-array has a hole after every lag, so '
        'it does not belong to that family.')]
s += [P('<b>Recommended wording for the paper:</b> "a uniform linear array (ULA) with inter-antenna spacing '
        'd = 4&lambda;", or "a wide-spaced ULA". If the word sparse is needed, use "sparse ULA" only as a '
        'description, and never present it as a new array type.', KEY)]

# ------------------------------------------------------------------ 2. why small gain
s += [P('2. Why the wide-spaced ULA gave only a small benefit', H1)]
s += [P('At the same 32 antennas, widening the spacing from &lambda;/2 to 4&lambda; multiplies the aperture by '
        'eight (1.55 m to 12.40 m). That is the attraction: a larger aperture gives finer angle resolution and a '
        'much larger Rayleigh distance (48 m to 3.1 km). This puts every user of the 300 m area inside the near '
        'field, where range can also be resolved. The cost is the <b>grating lobes</b>. With d = 4&lambda;, the '
        'array cannot tell apart directions whose sin&theta; differ by a multiple of &lambda;/d = 1/4. Figure 1 '
        'shows the far-field beampattern: the 4&lambda; ULA has <b>7 peaks at full gain</b> (peak sidelobe level '
        '0 dB), while the &lambda;/2 ULA has one main beam with -13 dB sidelobes.')]
s += [fig(FIG + 'codeT_fig_beams.png', 15.5,
          'Figure 1. Far-field beampattern of the four 32-antenna arrays. The 4&lambda; ULA has six grating lobes '
          'as strong as the main beam. Nested and coprime have one main beam and high but lower sidelobes (about '
          '-5 dB).')]
s += [P('What this meant in our earlier results (Code R, spacing 4.5&lambda;, 6 random setups):')]
s += [bullets([
    '<b>Data detection improved a little in the user-rich case.</b> At K = 40 and 10 dB the wide ULA + CAPE gave '
    'BER 0.100 vs 0.142 for nested + CAPE, and sum SE 145 vs 118 bit/s/Hz. In the antenna-rich case (K = 8) the '
    'two were tied within the noise. With only 6 setups the error bars were wide (about &plusmn;20 bit/s/Hz on '
    'the sum SE), so even the user-rich gain was not firmly established.',
    '<b>Localization failed in some setups.</b> A grating lobe points to a wrong position, and CAPE did not always '
    'remove it (median errors of 1.5 m, 4 m, 10 m and 40 m in some setups, against centimetres for nested).',
    '<b>Why detection survives but localization does not:</b> a grating-lobe ghost gives a steering vector that '
    'still matches the received signal well, so the channel estimate is usable for detection even when the '
    'position is wrong. Detection is tolerant to the ambiguity; localization is not.',
    'Code S, which chose the spacing from the full design procedure, showed the same picture: a small gain in '
    'BER and sum rate, nothing decisive.',
])]
s += [P('<b>One real advantage our simulation cannot show.</b> The 4&lambda; ULA has no closely spaced antennas '
        '(w(1) = 0) and the lowest coupling leakage of the four arrays (Table 1). That matters in hardware. But '
        'our channel model has <b>no mutual-coupling model</b>: the channel is an ideal spherical-wave or '
        'correlated-Rayleigh channel. So the low-coupling benefit can be stated as a hardware property, but it '
        'will never appear in our BER curves.', NOTE)]

# ------------------------------------------------------------------ 3. nested vs co-array
s += [P('3. Nested array vs "difference co-array"', H1)]
s += [P('3.1 What each one is', H2)]
s += [P('The <b>physical array</b> is the set of antennas you actually build. The <b>difference co-array</b> is the '
        'set of all differences x<sub>m</sub> - x<sub>n</sub> between antenna positions. It appears because the '
        'covariance entry between antennas m and n of a far-field source depends only on that difference: '
        'R<sub>mn</sub> ~ exp(j2&pi;(x<sub>m</sub>-x<sub>n</sub>)sin&theta;/&lambda;). If we rearrange the '
        'covariance entries by their lag, we obtain samples of a longer <b>virtual</b> array. Methods like '
        'co-array MUSIC (with spatial smoothing) then work on that virtual array.')]
s += [bullets([
    '<b>Every array has a difference co-array</b>, even the plain ULA (its lags are -31 ... 31, i.e. 63 values, so it '
    'gains nothing).',
    'Nested and coprime arrays are <b>designed</b> so that their co-array is long and has no holes: 543 hole-free '
    'lags for our nested array and 287 for the coprime array, from only 32 antennas (Figure 2, Table 1).',
    'So "nested array" names a <b>physical geometry</b>, and "difference co-array processing" names a '
    '<b>signal-processing method</b> applied to that geometry. They are not two arrays.',
])]
s += [fig(FIG + 'codeT_fig_layouts.png', 15.5,
          'Figure 2. The four 32-antenna arrays (dots) and their hole-free difference co-arrays (ticks). The '
          '4&lambda; ULA has no hole-free co-array: all its lags are multiples of 8 half-wavelengths.')]
s += [P('3.2 What our code actually does', H2)]
s += [P('In Code R, scheme 3 ("Nested + CAPE") and scheme 4 ("Difference co-array + CAPE") have the same geometry '
        "('nested'), the same 32 antennas and the same CAPE estimator. The <b>only</b> difference is the flag "
        '<font face="Courier">en_idea</font>. With it switched on, inside <font face="Courier">csi_acquire_cape</font> '
        'each AP forms the sample covariance of its pilot samples, builds the virtual array, runs MUSIC, and '
        'turns the MUSIC angles (plus a range search) into <b>extra candidate positions</b>. These are added to '
        "the list. The final position is chosen by the <b>same</b> CAPE fused likelihood computed on the "
        '<b>physical</b> antennas (Figure 3).')]
s += [fig(FIG + 'codeT_fig_processing.png', 15.5,
          'Figure 3. Same antennas, two processing chains. The co-array branch can only add candidates; the '
          'decision is always made by the same CAPE likelihood.')]
s += [P('3.3 Why the two curves are almost identical', H2)]
s += [numbered([
    '<b>The co-array can only add candidates.</b> If the candidates found by CAPE on the physical array already '
    'contain the correct position (and with 4 fused APs they almost always do), the extra candidates never win '
    'and the result is exactly the same.',
    '<b>Our problem is one user at a time.</b> With orthogonal pilots, each user\'s pilot is separated before '
    'estimation, so each AP estimates one channel at a time. The famous advantage of the co-array, resolving '
    '<i>more sources than antennas</i> (O(N<super>2</super>) degrees of freedom), matters when many sources must '
    'be separated at once from one covariance, as in DOA estimation [1], [4], [5]. It is not the bottleneck here.',
    '<b>The co-array model assumes plane waves.</b> The lag structure R<sub>mn</sub> = f(x<sub>m</sub>-x<sub>n</sub>) '
    'holds only for far-field sources. For our sparse arrays the Rayleigh distance is 3.1-3.7 km, so every link '
    'in the 300 m area is near field. The curved wavefront breaks the lag structure, so the co-array angles are '
    'only approximate. Our own code header already notes this.',
    '<b>Very few samples for the covariance.</b> The covariance is estimated from &tau; = K pilot samples (8 or '
    '40). Co-array MUSIC needs many snapshots to work well [5].',
    '<b>The virtual aperture was capped.</b> In Code R/S the smoothed virtual subarray was limited to 80 elements '
    '(<font face="Courier">ca_Nsub_max = 80</font>), a 3.95 m virtual aperture, which is <i>smaller</i> than the '
    '13.55 m physical aperture of the nested array. So the co-array angles were coarser than what the physical '
    'array already sees. Even without the cap, spatial smoothing keeps about half of the 543 lags, i.e. 272 '
    'elements = 13.55 m, which is the same as the physical aperture. <b>The co-array adds degrees of freedom, not '
    'resolution beyond the physical aperture.</b>',
    '<b>The earlier comparison was not paired.</b> Code R gave every scheme a different random seed (the seed '
    'included the scheme number), even though a code comment says otherwise. So nested and co-array saw '
    'different channels and noise. Their small differences (about 0.5%) were smaller than the Monte Carlo error '
    '(about &plusmn;3%), i.e. pure noise.',
])]
s += [P('3.4 How the literature uses it, and how to present it', H2)]
s += [P('The professor\'s remark matches the literature. Papers choose a sparse physical array (nested, coprime, '
        'MISC, ...) and then apply co-array processing to it. A typical comparison is either (a) different '
        'arrays processed with the same method, or (b) different methods on the same array. Nobody lists '
        '"nested" and "difference co-array" as two competing arrays, because the second one is not an array.')]
s += [P('<b>For our paper:</b> keep the arrays and the processing separate. Write "nested array (direct '
        'processing)" and "nested array with difference co-array processing". If the co-array processing gives no '
        'measurable gain (see Section 5), it is honest and stronger to say so in one sentence and explain the '
        'reasons above, rather than show two overlapping curves with different names.', KEY)]

# ------------------------------------------------------------------ 4. Code T
s += [P('4. The new comparison code: Code_T_ArrayComparison.m', H1)]
s += [P('Code T is built from Code R (through Code S, which ran correctly on your machine). The scenario, the CAPE '
        'estimator, the four detectors, the SNR grid (0 to 30 dB), L = 8 APs, the two cases (K = 8 and K = 40 with '
        '32 antennas per AP) and the user-count sweep are unchanged. Only the array set and the comparison '
        'method change.')]
s += [P('4.1 The six schemes (32 physical antennas per AP, all with CAPE)', H2)]
s += [table([
    ['#', 'Scheme', 'Physical array', 'Processing'],
    ['1', 'ULA (&lambda;/2)', 'dense ULA, d = &lambda;/2, aperture 1.55 m', 'direct (physical array)'],
    ['2', 'ULA (4&lambda;)', 'wide-spaced ULA, d = 4&lambda;, aperture 12.40 m', 'direct'],
    ['3', 'Nested', 'two-level nested (16 + 16), aperture 13.55 m', 'direct'],
    ['4', 'Nested + co-array', '<b>same antennas as #3</b>', 'direct + difference co-array MUSIC candidates'],
    ['5', 'Coprime', 'extended coprime, M = 8, N = 17, aperture 12.75 m', 'direct'],
    ['6', 'Coprime + co-array', '<b>same antennas as #5</b>', 'direct + difference co-array MUSIC candidates'],
], [0.8 * cm, 3.4 * cm, 6.6 * cm, 5.8 * cm])]
s += [Spacer(1, 6)]
s += [P('The three sparse arrays have almost the same aperture (12.4 to 13.6 m), so differences between them come '
        'from the <b>layout</b> (grating lobes vs co-array structure), not from size. Co-array processing is not '
        'applied to the 4&lambda; ULA because it has no hole-free co-array.')]
s += [P('<b>Table 1.</b> Geometry metrics (computed from the exact positions used in the code).', BODY)]
s += [table([
    ['Array', 'Aperture', 'Rayleigh dist.', 'Hole-free co-array', 'w(1)', 'Coupling leakage', 'Peak sidelobe'],
    ['ULA &lambda;/2', '1.55 m', '48 m', '63', '31', '0.865', '-13.2 dB'],
    ['ULA 4&lambda;', '12.40 m', '3075 m', '1 (holes)', '0', '0.209', '0 dB (grating lobes)'],
    ['Nested', '13.55 m', '3672 m', '543', '16', '0.772', '-5.0 dB'],
    ['Coprime (8,17)', '12.75 m', '3251 m', '287', '2', '0.432', '-5.3 dB'],
], [2.6 * cm, 1.9 * cm, 2.2 * cm, 2.6 * cm, 1.2 * cm, 2.6 * cm, 3.5 * cm])]
s += [P('w(1) = number of antenna pairs only &lambda;/2 apart (these couple most strongly). Coupling leakage uses the '
        'standard banded model (1/k for pairs k half-wavelengths apart).', SMALL), Spacer(1, 8)]

s += [P('4.2 What changed compared with Code R, and why', H2)]
s += [numbered([
    '<b>Array set</b> as requested: ULA &lambda;/2, ULA 4&lambda;, nested, coprime, plus co-array processing on '
    'nested and coprime. The spacing of the wide ULA is one parameter at the top: '
    '<font face="Courier">ulaw_d_lambda = 4</font>.',
    '<b>Paired comparison.</b> All six schemes now use the same random seed in a given setup and SNR, so they '
    'see the same user drop, the same pilot noise and the same data symbols. For nested vs nested + co-array '
    '(and coprime vs coprime + co-array) the channels are identical too, so any difference is caused only by '
    'the co-array processing. The co-array covariance noise is drawn from its own side stream, so it does not '
    'shift the random numbers used later.',
    '<b>No cap on the co-array.</b> <font face="Courier">ca_Nsub_max = 1000</font>, so co-array MUSIC uses its full '
    'smoothed virtual aperture. This gives the co-array its best chance, so a "no gain" result cannot be blamed '
    'on the cap.',
    '<b>More random setups:</b> <font face="Courier">nSetE = 10</font> (Code R had 6). Use 6 for a quick look; '
    'use at least 10 for the paper.',
    '<b>User-count sweep</b> now includes all six schemes and reuses the prebuilt geometries (same result, '
    'faster).',
    '<b>Results are saved before plotting</b> to <font face="Courier">codeT_arraycomp_results.mat</font>, including '
    'the array geometries. <font face="Courier">plot_codeT_figures.m</font> redraws all figures from the workspace or '
    'from that file, so the figure problems we had with Code S cannot happen again.',
])]

s += [P('4.3 The figures and what each one answers', H2)]
s += [table([
    ['Fig.', 'Content', 'Question it answers'],
    ['T1', 'Physical antennas and hole-free co-array of each array', 'What is physically different between the arrays?'],
    ['T2', 'Far-field beampattern', 'Where are the grating lobes (4&lambda; ULA) and sidelobes?'],
    ['T3', 'Channel NMSE vs SNR, both cases', 'Which array gives the best channel estimate?'],
    ['T4', 'BER vs SNR, both cases', 'Which array gives the best detection?'],
    ['T5', 'Sum SE vs SNR, both cases', 'Which array gives the best rate?'],
    ['T6', 'BER and sum SE vs number of users K at 10 dB', 'Does the ranking change when users outnumber antennas?'],
    ['T7', 'Paired change caused by co-array processing (with error bars and % identical)',
     "<b>The professor's question:</b> does co-array processing change anything?"],
    ['T8', 'Median localization error vs SNR', 'Do grating lobes break localization?'],
    ['T9', 'Aperture, hole-free co-array, w(1), coupling leakage', 'Hardware cost and coupling of each array'],
], [1.1 * cm, 7.4 * cm, 8.1 * cm])]
s += [Spacer(1, 6)]
s += [P('The console also prints a table at 10 dB with the mean and the standard error over the setups for every '
        'scheme, plus the percentage of (SNR, setup) points where the co-array version gave <b>exactly</b> the same '
        'BER as the direct version.')]
s += [P('<b>How it was checked before giving it to you.</b> Both files parse without errors. A reduced version (1 '
        'setup, 3 SNR points, few symbols) ran the full computation from start to end in GNU Octave (both cases, '
        'all six schemes, the user sweep) and saved the results file; all nine figures and the console table were '
        'then drawn from it with <font face="Courier">plot_codeT_figures.m</font>, which uses the same figure code. '
        'The run also confirmed that the pairing is exact: with the co-array switched on, the results were '
        'bit-identical to the plain array wherever the extra candidates did not win. Octave is not MATLAB, and '
        'such a tiny run proves only that the code works, not what the results will be. The real numbers must '
        'come from your MATLAB run.', NOTE)]

# ------------------------------------------------------------------ 5. how to judge
s += [P('5. How to read the results and decide', H1)]
s += [numbered([
    '<b>Use the error bars.</b> If two schemes differ by less than about two standard errors (printed in the '
    'console), treat them as equal. Do not claim a gain from a difference inside the noise.',
    '<b>Co-array question (Figure T7).</b> If the curves stay at 0 % and the title says the co-array version was '
    'identical in most points, co-array processing does not change any decision in our system. Then say so in '
    'one sentence and keep only the physical arrays in the main figures. If it gives a clear gain somewhere '
    '(most likely K = 40 at low SNR), that is a real, explainable result and can be reported.',
    '<b>ULA 4&lambda; vs nested vs coprime (Figures T3-T6).</b> These have nearly the same aperture, so this is a '
    'clean test of layout. If the 4&lambda; ULA is equal or better in BER/sum SE but clearly worse in localization '
    '(T8), the honest summary is: "wide spacing helps detection, grating lobes hurt localization".',
    '<b>ULA &lambda;/2 vs the sparse arrays.</b> The dense ULA has an 8-times smaller aperture and a 48 m Rayleigh '
    'distance, so many of its links are far field. A large difference here is mainly an <b>aperture</b> effect, '
    'which is expected and is not a new finding by itself.',
    '<b>Coupling.</b> Since coupling is not simulated, Table 1 / Figure T9 can support a hardware argument (for '
    'example, coprime and 4&lambda; ULA have far fewer close pairs than nested), but not a BER claim.',
])]
s += [P('<b>Expectation, to be confirmed by the run (not a result):</b> based on the mechanism above, we expect '
        'Figure T7 to stay close to zero, the three sparse arrays to be close in BER and sum SE (with the 4&lambda; '
        'ULA possibly a little ahead at K = 40), and the 4&lambda; ULA to be the worst in localization.', NOTE)]

# ------------------------------------------------------------------ 6. limitations
s += [P('6. Limitations to keep in mind', H1)]
s += [bullets([
    '<b>EMISC is not included.</b> The professor allowed "coprime or EMISC". We used the coprime array because its '
    'construction is standard and unambiguous. EMISC (or MISC [9]) can be added later as one more case in '
    '<font face="Courier">geom_targetN</font>.',
    '<b>No mutual coupling</b> in the channel model (Section 2).',
    '<b>Co-array processing uses a plane-wave virtual array</b>, as in the standard literature. Near-field-aware '
    'co-array methods exist, but they would be a new method, not part of this comparison.',
    '<b>Run time.</b> Six schemes with 10 setups and a six-scheme user sweep take roughly 1.6 to 2 times as long as '
    'Code R (Code R used 6 setups and a three-scheme sweep). Set <font face="Courier">nSetE = 6</font> for a quick '
    'first look.',
])]

# ------------------------------------------------------------------ refs
s += [P('References', H1)]
refs = [
    '[1] P. Pal and P. P. Vaidyanathan, "Nested arrays: A novel approach to array processing with enhanced degrees '
    'of freedom," <i>IEEE Trans. Signal Process.</i>, vol. 58, no. 8, 2010.',
    '[2] P. P. Vaidyanathan and P. Pal, "Sparse sensing with co-prime samplers and arrays," <i>IEEE Trans. Signal '
    'Process.</i>, vol. 59, no. 2, 2011. (Coprime array built from two ULAs.)',
    '[3] P. Pal and P. P. Vaidyanathan, "Coprime sampling and the MUSIC algorithm," <i>Proc. IEEE DSP/SPE '
    'Workshop</i>, 2011. (Extended coprime array, 2M + N - 1 sensors.)',
    '[4] C.-L. Liu and P. P. Vaidyanathan, "Remarks on the spatial smoothing step in coarray MUSIC," <i>IEEE Signal '
    'Process. Lett.</i>, 2015; and C.-L. Liu, PhD thesis, Caltech, 2018: '
    'thesis.library.caltech.edu/10970/1/Liu_Chun-Lin_2018.pdf',
    '[5] M. Wang and A. Nehorai, "Coarrays, MUSIC, and the Cram&eacute;r-Rao bound," <i>IEEE Trans. Signal Process.</i>, '
    '2017. arXiv:1605.03620 (co-array MUSIC saturates when sources exceed sensors; behaviour with snapshots).',
    '[6] C. Zhou, C. You, H. Zhang, L. Chen, S. Shi, "Sparse array enabled near-field communications: Beam pattern '
    'analysis and hybrid beamforming design," arXiv:2401.05690, 2024 (linear sparse array = ULA with spacing '
    '&gt; &lambda;/2; grating lobes; extended coprime array).',
    '[7] A. Kosasih, &Ouml;. T. Demir, E. Bj&ouml;rnson, "Near-field beamfocusing, localization, and channel estimation '
    'with modular linear arrays," arXiv:2505.07991, 2025.',
    '[8] Li <i>et al.</i>, modular XL-MIMO near-field analysis, arXiv:2308.11289, 2023 (larger aperture improves '
    'resolution, grating lobes degrade user SINR).',
    '[9] Z. Zheng, W.-Q. Wang, Y. Kong, Y. D. Zhang, "MISC array: A new sparse array design achieving increased '
    'degrees of freedom and reduced mutual coupling effect," <i>IEEE Trans. Signal Process.</i>, 2019.',
    '[10] H. L. Van Trees, <i>Optimum Array Processing</i>, Wiley, 2002 (ULA definition, grating lobes for '
    'd &gt; &lambda;/2).',
    '[11] Wang <i>et al.</i>, "Exploring the advantages of sparse arrays in near-field XL-MIMO systems: Beam analysis '
    'and EDoF function," arXiv:2501.09234, 2025.',
]
s += [P(r, REF) for r in refs]


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(colors.HexColor('#777777'))
    canvas.drawString(2 * cm, 1.2 * cm, 'Array comparison note - Code T')
    canvas.drawRightString(19 * cm, 1.2 * cm, 'page %d' % doc.page)
    canvas.restoreState()


doc = SimpleDocTemplate(OUT, pagesize=A4, leftMargin=2.1 * cm, rightMargin=2.1 * cm,
                        topMargin=1.8 * cm, bottomMargin=2 * cm,
                        title='Array comparison: answers and Code T', author='')
doc.build(s, onFirstPage=footer, onLaterPages=footer)
print('wrote', OUT)

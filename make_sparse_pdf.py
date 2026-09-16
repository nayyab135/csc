#!/usr/bin/env python3
# Build sparse_array_justification.pdf with reportlab (no LaTeX available).
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, PageBreak,
                                Table, TableStyle, ListFlowable, ListItem, HRFlowable)

OUT = "/home/user/csc/sparse_array_justification.pdf"

ss = getSampleStyleSheet()
H1 = ParagraphStyle('H1', parent=ss['Heading1'], fontSize=14, spaceBefore=14, spaceAfter=6, textColor=colors.HexColor('#1a3d6d'))
H2 = ParagraphStyle('H2', parent=ss['Heading2'], fontSize=11.5, spaceBefore=10, spaceAfter=4, textColor=colors.HexColor('#22548c'))
BODY = ParagraphStyle('BODY', parent=ss['Normal'], fontSize=9.7, leading=13.6, alignment=TA_JUSTIFY, spaceAfter=6)
EQ = ParagraphStyle('EQ', parent=ss['Normal'], fontName='Courier', fontSize=9.2, leading=13, leftIndent=18, spaceBefore=3, spaceAfter=6, textColor=colors.HexColor('#333333'))
NOTE = ParagraphStyle('NOTE', parent=BODY, fontSize=9.0, leading=12.5, leftIndent=10, textColor=colors.HexColor('#7a4a00'), backColor=colors.HexColor('#fff6e6'), borderPadding=4, spaceBefore=4, spaceAfter=8)
CAP = ParagraphStyle('CAP', parent=BODY, fontSize=8.8, textColor=colors.HexColor('#444444'), leftIndent=8)
TITLE = ParagraphStyle('TITLE', parent=ss['Title'], fontSize=19, leading=23, textColor=colors.HexColor('#12305c'))
SUB = ParagraphStyle('SUB', parent=ss['Normal'], fontSize=10.5, alignment=TA_CENTER, textColor=colors.HexColor('#555555'))
REF = ParagraphStyle('REF', parent=BODY, fontSize=8.9, leading=12, leftIndent=14, firstLineIndent=-14, spaceAfter=3)

def P(t, s=BODY): return Paragraph(t, s)
def bullets(items, s=BODY):
    return ListFlowable([ListItem(Paragraph(x, s), leftIndent=12) for x in items], bulletType='bullet', start='square', leftIndent=14)

story = []

# ---------------- Title ----------------
story += [Spacer(1, 1.2*cm),
    P("Incorporating Sparse Arrays into a Cell-Free<br/>Near-Field / Far-Field XL-MIMO Uplink", TITLE),
    Spacer(1, 6),
    P("Motivation, theory, design choices, references, and simulation evidence", SUB),
    Spacer(1, 10),
    HRFlowable(width="100%", thickness=1, color=colors.HexColor('#12305c')),
    Spacer(1, 10),
    P("<b>Scope.</b> This note explains why a sparse per-access-point (AP) array is worth adding to the "
      "existing cell-free hybrid near-field/far-field XL-MIMO uplink, what concrete benefits it brings, "
      "where the idea originates, which references support each design choice, how the advantages can be "
      "exploited, and what the five selected figures produced by the simulation code demonstrate. "
      "Closed-form (analytic) results are stated exactly; quantities that depend on the Monte-Carlo run are "
      "labelled as such so nothing is overstated.", BODY),
]

# ---------------- 1. Why ----------------
story += [P("1. Why add a sparse array?", H1),
    P("The system places <b>N = 64</b> antennas per AP at half-wavelength spacing (d = lambda/2) at "
      "f<sub>c</sub> = 3 GHz (lambda = 0.1 m). The per-AP aperture is therefore D = (N-1)&#183;lambda/2 = 3.15 m, "
      "and the Rayleigh (near-field) distance is", BODY),
    P("d_Ray = 2 D^2 / lambda = 2 &#183; (3.15)^2 / 0.1 ~ 198 m.", EQ),
    P("Two hard constraints motivate a sparse array:", BODY),
    bullets([
      "<b>The half-wavelength array cannot be made smaller.</b> Reducing N shrinks D, which shrinks d_Ray "
      "(quadratically), so the prominent near-field region on which the proposed Cross-AP list-SIC detector "
      "relies would disappear. At lambda/2 spacing, N is effectively floored.",
      "<b>Many antennas are expensive.</b> Each element needs an RF chain and contributes to fronthaul and "
      "detector complexity. For XL-MIMO this is the dominant hardware cost.",
    ]),
    P("A <b>sparse array</b> increases the element spacing to d = s&#183;lambda/2 (sparsening factor s &gt; 1). "
      "Because d_Ray depends on the <i>aperture</i>, not the element count, the same near-field distance is "
      "reached with far fewer antennas, or a much larger near-field region is obtained at the same N. The price "
      "is grating lobes, which the near-field regime makes range-selective rather than eliminating (Section 4). "
      "The remainder of this note quantifies this trade.", BODY),
]

# ---------------- 2. Origin & references ----------------
story += [P("2. Where the idea comes from", H1),
    P("Sparse arrays are a classical antenna concept (thinned/minimum-redundancy arrays, nested and coprime "
      "arrays) recently revisited for near-field XL-MIMO and &#8216;sparse MIMO&#8217; communications. The "
      "specific results used here come from:", BODY),
    bullets([
      "<b>Sparse MIMO paradigm &amp; interference suppression</b> for near- and far-field communications "
      "[R5], which shows sparse spacing raises the spatial degrees of freedom and improves interference "
      "rejection when users are in the near field.",
      "<b>Near-field beam-focusing pattern and grating-lobe characterization for modular XL-arrays</b> [R1], "
      "the origin of the 2-D (angle x range) beam-focusing analysis and the statement that the spherical "
      "wavefront makes grating lobes range-selective.",
      "<b>Advantages of sparse arrays in near-field XL-MIMO (EDoF)</b> [R2]: a closed form showing the "
      "effective degrees of freedom increase with sparsity up to a bound.",
      "<b>Near-field XL-MIMO tutorial</b> [R3] for the non-uniform spherical-wave (NUSW) model and the "
      "aperture / Rayleigh-distance relationship.",
      "<b>Reconfigurable array thinning</b> [R4] and <b>quasi-distributed (ULA to MRA) grating-lobe</b> "
      "analysis [R6], which frame the uniform-vs-non-uniform choice and the distributed setting.",
    ]),
    P("The link to the group&#8217;s own line of work is the detector: the near-field spherical wavefront is "
      "exactly what lets the Cross-AP list-SIC / IDD receiver separate users that a far-field array cannot, so "
      "anything that enlarges the near-field region feeds that receiver more of the structure it exploits.", BODY),
]

# ---------------- 3. System model / math ----------------
story += [P("3. System model and mathematical formulation", H1),
    P("<b>Array response.</b> Let p<sub>m</sub> be the position of element m along the array axis (centred). "
      "Under the plane-wave (far-field) model the steering vector depends on angle only:", BODY),
    P("a_FF,m(theta) = exp( j 2 pi p_m sin(theta) / lambda ).", EQ),
    P("Under the non-uniform spherical-wave (NUSW) near-field model it depends on angle <i>and</i> range r "
      "through the exact element-to-source distance:", BODY),
    P("d_m(r,theta) = sqrt( r^2 + p_m^2 - 2 r p_m sin(theta) ),", EQ),
    P("a_NF,m(r,theta) = ( r / d_m ) &#183; exp( - j 2 pi d_m / lambda ).", EQ),
    P("<b>Grating lobes.</b> For uniform spacing d, a beam steered to theta_0 has ambiguous replicas at", BODY),
    P("sin(theta_g) = sin(theta_0) +/- m &#183; lambda / d ,   m = 1, 2, ...", EQ),
    P("that are visible whenever |sin(theta_g)| &lt;= 1. With d = s&#183;lambda/2 this becomes "
      "sin(theta_g) = sin(theta_0) +/- 2m/s, so for s = 4 and broadside steering the grating lobes sit at "
      "sin(theta_g) = +/-0.5 and +/-1.0 (i.e. +/-30 deg and +/-90 deg). In the far field these replicas reach the "
      "full main-lobe height (0 dB): a genuine ambiguity. In the near field the non-linear phase in d_m(r,theta) "
      "prevents them from adding up coherently except near the focal range, so they become range-selective [R1].", BODY),
    P("<b>Rayleigh distance and antenna count.</b> With aperture D = (N-1)&#183;s&#183;lambda/2,", BODY),
    P("d_Ray(s) = 2 D^2 / lambda = (N-1)^2 s^2 lambda / 2   =>   d_Ray(s) = s^2 &#183; d_Ray(1).", EQ),
    P("To hold a target near-field distance d_Ray while sparsening, the required antenna count is", BODY),
    P("N(s) = 1 + sqrt( 2 d_Ray / lambda ) / s .", EQ),
    P("For d_Ray = 198 m and lambda = 0.1 m this gives N(1) = 64 (the baseline) and N(4) ~ 17, i.e. a "
      "<b>74% reduction in antennas (and RF chains) for the same near-field coverage</b> &#8212; a closed-form, "
      "assumption-free result.", BODY),
    P("<b>Collinear users.</b> Two users at the same angle theta but different ranges r1, r2 have identical "
      "far-field steering, a_FF(theta), so the 2-user channel is rank-deficient and no linear receiver can "
      "separate them: the post-MMSE SINR saturates (0 dB ceiling) at any SNR. Their near-field signatures "
      "a_NF(r1,theta), a_NF(r2,theta) differ, and a larger (sparse) aperture makes them more orthogonal, so they "
      "become separable [R5]. This is the mechanism the Cross-AP detector uses.", BODY),
]

# ---------------- 4. Benefits & exploitation ----------------
story += [PageBreak(), P("4. Benefits and how to exploit them", H1),
    P("The advantages of a sparse array can be exploited along four axes; each is tied to a reference and to a "
      "concrete way of using it in this system.", BODY)]

btab = [
 ["#", "Benefit", "Mechanism / how to exploit", "Ref"],
 ["B1", "Antenna & RF-chain reduction",
        "Same d_Ray with ~74% fewer elements (N(s)=1+sqrt(2 d_Ray/lambda)/s). Use to cut per-AP hardware, "
        "fronthaul and detector cost at equal near-field reach.", "R2,R4"],
 ["B2", "Larger near-field region",
        "d_Ray(s)=s^2 d_Ray. At fixed N, sparsening pushes more user-AP links into the near field, where the "
        "Cross-AP detector gains. Exploit by reporting the near-field link-fraction increase (Fig. S2).", "R3,R5"],
 ["B3", "Range resolution / EDoF",
        "Larger aperture raises effective DoF and resolves collinear (same-angle) users by range that the far "
        "field cannot. Exploit directly in the SIC/IDD receiver (Fig. S3).", "R2,R5"],
 ["B4", "Grating lobes made range-selective",
        "Near-field spherical wavefront defocuses replicas away from the focal range; cross-AP fusion further "
        "averages the residual across APs. This is a COST that is mitigated, not a feature (Figs. S4, S5).", "R1,R6"],
]
t = Table(btab, colWidths=[0.9*cm, 3.5*cm, 8.3*cm, 1.6*cm])
t.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#22548c')),
    ('TEXTCOLOR',(0,0),(-1,0),colors.white),
    ('FONTSIZE',(0,0),(-1,-1),8.4),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#aab7c8')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white, colors.HexColor('#eef3f9')]),
    ('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),
    ('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),
]))
story += [t, Spacer(1,8)]

story += [P("4.1 Uniform vs non-uniform sparse array", H2),
    P("<b>Recommendation: start with a uniform sparse ULA;</b> keep non-uniform (thinned / nested / coprime / "
      "minimum-redundancy) as a robustness comparator, and let the BER after cross-AP combining decide.", BODY),
    P("<b>Reasoning (reference-backed, not the manifold argument).</b> A uniform sparse array places its grating "
      "lobes at <i>predictable</i> angles (sin = 2m/s) and is exactly what the modular/sparse XL-MIMO literature "
      "models [R1,R3], so results are directly comparable to that body of work and the grating angles are "
      "analytically tractable. A non-uniform array trades the discrete grating spikes for a raised sidelobe floor "
      "everywhere &#8212; better worst-case robustness to arbitrary interferer geometry, at the cost of a higher "
      "average sidelobe and an irregular manifold [R6, nested/coprime arrays].", BODY),
    NOTE.clone('n1') if False else P("<b>Correction to a common justification.</b> The irregular-manifold penalty "
      "of non-uniform arrays matters only for codebook / beam-training / polar-domain processing. The present "
      "detector builds MMSE combiners from <i>estimated</i> channels and uses no angle-range dictionary, so the "
      "&#8216;preserve the structured manifold&#8217; argument does not apply here. The correct reasons to prefer "
      "uniform are comparability and analytical tractability &#8212; and, if a near-field channel-estimation "
      "codebook is added later, the manifold argument then becomes valid future-work justification.", NOTE),
]

# ---------------- 5. The five figures ----------------
story += [P("5. What the selected figures show", H1),
    P("The code produces five curated figures (set <font face='Courier'>sparse_study = true</font>). A sixth, the "
      "2-D near-field focusing pattern, is available behind <font face='Courier'>show_mechanism = true</font> as "
      "an optional mechanism illustration.", BODY)]

def fig(title, whatit, result):
    story.append(P("<b>"+title+"</b>", H2))
    story.append(P("<b>What it shows.</b> "+whatit, BODY))
    story.append(P("<b>Result.</b> "+result, BODY))

fig("Fig. S1 &#8212; Antenna / complexity reduction at equal near-field distance",
    "d_Ray vs N for lambda/2 vs sparse spacing, and the antenna count required to hold the baseline d_Ray as a "
    "function of the sparsening factor s.",
    "<b>Closed-form:</b> holding d_Ray = 198 m needs N ~ 17 at s = 4 versus N = 64 at lambda/2, a 74% "
    "reduction in antennas and RF chains. This is the strongest, assumption-free argument for the array.")
fig("Fig. S2 &#8212; Near-field link fraction vs sparsening (this deployment)",
    "A Monte-Carlo over the actual 300 x 300 m deployment measuring the percentage of user-AP links that are "
    "near-field (link distance &lt; d_Ray(s)) as s grows.",
    "<b>Simulation (read the printed values):</b> the near-field fraction increases monotonically with s "
    "because d_Ray(s)=s^2 d_Ray; the s=1 and s=4 percentages are printed to the console. This is the evidence "
    "that sparse arrays are <i>relevant</i> to this geometry rather than assumed to be.")
fig("Fig. S3 &#8212; Collinear-user SINR, far/near x dense/sparse",
    "Post-MMSE SINR of a desired user against one collinear interferer (same angle, delta r = 2 m) for the four "
    "combinations of {far-field, near-field} x {dense, sparse}.",
    "<b>Analytic + simulation:</b> both far-field curves saturate at the 0 dB rank-deficiency ceiling (identical "
    "steering, unrecoverable at any SNR); the near-field curves rise with SNR, and the near-field sparse array "
    "rises fastest because its larger aperture gives the finest range resolution. This is the detection-relevant "
    "benefit that connects to the Cross-AP receiver.")
fig("Fig. S4 &#8212; Grating-lobe cost and its near-field mitigation",
    "Left: the beam gain vs angle at s = 4 for the sparse array under the far-field model (full-height grating "
    "lobes) vs the near-field model (range-selective), with a dense reference. Right: the peak grating-lobe "
    "level vs element spacing, far-field vs near-field.",
    "<b>Simulation (from your run):</b> under the far-field model the grating replicas reach ~0 dB (full "
    "ambiguity); under the near-field model they are held roughly 6-9 dB lower. The honest reading is that the "
    "near field provides a ~6-9 dB range-selective mitigation, NOT immunity &#8212; a residual remains at the "
    "focal range. State it as &#8216;range-selective&#8217;, never &#8216;removed&#8217;.")
fig("Fig. S5 &#8212; Cross-AP combining vs the per-AP grating lobe (novel test)",
    "An interferer is placed exactly on AP#1&#8217;s grating angle (sin = 2/s), where AP#1&#8217;s sparse array "
    "cannot null it; the other APs see it at ordinary angles. Desired-user SINR is compared at AP#1 alone vs "
    "after cross-AP LSFD fusion, dense vs sparse.",
    "<b>Simulation (read the printed dB recovery):</b> the sparse AP#1-alone curve sits below the dense one "
    "(the grating-lobe liability is real); cross-AP fusion is expected to recover most of that loss because the "
    "remaining APs resolve the interferer. The printed &#8216;recovers X dB&#8217; line is the number that turns "
    "&#8216;sparse arrays have interesting beampatterns&#8217; into &#8216;sparse arrays work in this detector&#8217;. "
    "If fusion does NOT recover it, that is a legitimate negative result to report, not a failure.")

# ---------------- 6. Limitations ----------------
story += [P("6. Honest limitations and the open question", H1),
    bullets([
      "<b>The grating lobe is a real cost.</b> The near field mitigates it (~6-9 dB, range-selective); it does "
      "not remove it. Avoid the words &#8216;suppressed&#8217; / &#8216;immune&#8217; in the paper.",
      "<b>Benefit is conditional on the near-field fraction.</b> Sparse helps only for links that are in the "
      "near field; far-field APs gain nothing and still see fully-ambiguous grating lobes. Fig. S2 quantifies "
      "how many links are near-field, which is why it must accompany the benefit figures.",
      "<b>Figs. S1-S5 are per-array / per-scene analyses.</b> The genuinely novel, still-open claim &#8212; that "
      "cell-free cross-AP combining suppresses the residual grating lobes across APs &#8212; is only probed at the "
      "single-interferer level in Fig. S5. The decisive evidence is a full system-level BER / sum-rate comparison "
      "of {dense lambda/2, uniform sparse} under the four detectors, which requires threading the sparse geometry "
      "into the main Monte-Carlo (channel generation, near-field mask). That is the recommended next step.",
    ]),
]

# ---------------- 7. Conclusions ----------------
story += [P("7. Conclusions and recommended next steps", H1),
    P("The defensible claims, in order of strength:", BODY),
    bullets([
      "<b>(1) Hardware:</b> a uniform sparse array reaches the same near-field distance with ~74% fewer antennas "
      "(N: 64 -&gt; 17 at s = 4) &#8212; closed-form, unconditional.",
      "<b>(2) Detection:</b> the resulting large aperture resolves collinear users by range that a far-field "
      "array fundamentally cannot, directly benefiting the Cross-AP list-SIC / IDD receiver.",
      "<b>(3) Relevance:</b> on the actual 300 m deployment, sparsening moves a measurable fraction of links into "
      "the near field (Fig. S2), which is the precondition for (2).",
      "<b>(4) Cost, honestly:</b> sparsity causes grating lobes; the near field makes them range-selective "
      "(~6-9 dB), and cross-AP fusion is the mechanism proposed to handle the residual (Fig. S5).",
    ]),
    P("<b>Next step.</b> Run Figs. S1-S5. If S2 shows a meaningful near-field-fraction increase and S5 shows "
      "cross-AP fusion recovering the grating-lobe loss, proceed to the system-level BER/sum-rate integration "
      "with the sparse per-AP geometry &#8212; that figure decides whether sparse arrays belong in the paper. If "
      "S5 shows fusion does not recover the loss, report that as the finding and reconsider (e.g. non-uniform or "
      "geometry-aware detection) before investing in the full integration.", BODY),
]

# ---------------- References ----------------
story += [P("References", H1),
  P("[R1] H. Wang, Z. Xiao, Y. Zeng, <i>et al.</i>, &#8220;Near-Field Beam Focusing Pattern and Grating Lobe "
    "Characterization for Modular XL-Array,&#8221; arXiv:2305.05408, 2023.", REF),
  P("[R2] &#8220;Exploring the Advantages of Sparse Arrays in Near-Field XL-MIMO Systems: Beam Analysis and EDoF "
    "Function,&#8221; arXiv:2501.09234, 2025.", REF),
  P("[R3] Z. Wang, X. Mu, Y. Liu, <i>et al.</i>, &#8220;A Tutorial on Near-Field XL-MIMO Communications Towards "
    "6G,&#8221; arXiv:2310.11044, 2023.", REF),
  P("[R4] &#8220;Sparse Array Design for Near-Field MU-MIMO: A Reconfigurable Array Thinning Approach,&#8221; "
    "arXiv:2602.21973.", REF),
  P("[R5] &#8220;Enhancing Spatial Multiplexing and Interference Suppression for Near- and Far-Field "
    "Communications with Sparse MIMO,&#8221; arXiv:2408.01956, 2024.", REF),
  P("[R6] &#8220;Near-Field Communications with Grating Lobes for Quasi-Distributed Arrays: From ULA to MRA,&#8221; "
    "arXiv:2607.29341.", REF),
  P("[R7] P. Pal and P. P. Vaidyanathan, &#8220;Nested Arrays: A Novel Approach to Array Processing with Enhanced "
    "Degrees of Freedom,&#8221; IEEE Trans. Signal Process., 2010; P. P. Vaidyanathan and P. Pal, &#8220;Sparse "
    "Sensing with Coprime Arrays,&#8221; 2011 (non-uniform array DoF).", REF),
  P("[R8] M. Cui and L. Dai, &#8220;Channel Estimation for Extremely Large-Scale MIMO: Far-Field or Near-Field?,&#8221; "
    "IEEE Trans. Commun., 2022 (polar-domain / near-field manifold).", REF),
  Spacer(1,8),
  P("<i>Prepared as supporting material for the sparse-array extension of the cell-free NF/FF XL-MIMO uplink "
    "code. Analytic results are exact; simulation-dependent quantities should be filled from the console output "
    "of your own run.</i>", CAP),
]

doc = SimpleDocTemplate(OUT, pagesize=A4, topMargin=1.6*cm, bottomMargin=1.6*cm, leftMargin=1.7*cm, rightMargin=1.7*cm,
                        title="Sparse Arrays in Cell-Free NF/FF XL-MIMO", author="Research note")
doc.build(story)
print("WROTE", OUT)

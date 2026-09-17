#!/usr/bin/env python3
# Channel-estimation justification note (perfect vs estimated CSI, NF + FF).
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, PageBreak,
                                Table, TableStyle, ListFlowable, ListItem, HRFlowable)

OUT = "/home/user/csc/channel_estimation_note.pdf"
ss = getSampleStyleSheet()
H1  = ParagraphStyle('H1', parent=ss['Heading1'], fontSize=13.5, spaceBefore=13, spaceAfter=5, textColor=colors.HexColor('#1a3d6d'))
H2  = ParagraphStyle('H2', parent=ss['Heading2'], fontSize=11, spaceBefore=9, spaceAfter=3, textColor=colors.HexColor('#22548c'))
BODY= ParagraphStyle('BODY', parent=ss['Normal'], fontSize=9.7, leading=13.6, alignment=TA_JUSTIFY, spaceAfter=6)
EQ  = ParagraphStyle('EQ', parent=ss['Normal'], fontName='Courier', fontSize=9.1, leading=12.8, leftIndent=16, spaceBefore=2, spaceAfter=6, textColor=colors.HexColor('#333'))
NOTE= ParagraphStyle('NOTE', parent=BODY, fontSize=9.1, leading=12.8, leftIndent=8, textColor=colors.HexColor('#7a3b00'), backColor=colors.HexColor('#fff5e6'), borderPadding=5, spaceBefore=4, spaceAfter=8)
CAP = ParagraphStyle('CAP', parent=BODY, fontSize=8.7, textColor=colors.HexColor('#444'))
TITLE=ParagraphStyle('TITLE', parent=ss['Title'], fontSize=18, leading=22, textColor=colors.HexColor('#12305c'))
SUB = ParagraphStyle('SUB', parent=ss['Normal'], fontSize=10.5, alignment=TA_CENTER, textColor=colors.HexColor('#555'))
REF = ParagraphStyle('REF', parent=BODY, fontSize=8.9, leading=12, leftIndent=14, firstLineIndent=-14, spaceAfter=3)

def P(t,s=BODY): return Paragraph(t,s)
def bullets(items,s=BODY):
    return ListFlowable([ListItem(Paragraph(x,s),leftIndent=12) for x in items], bulletType='bullet', start='square', leftIndent=14)

st=[]
st += [Spacer(1,1.1*cm),
  P("Channel State Information in Near-Field Cell-Free XL-MIMO", TITLE),
  Spacer(1,5),
  P("Perfect CSI vs Estimated CSI: how each is computed, why, and with which references", SUB),
  Spacer(1,9), HRFlowable(width="100%", thickness=1, color=colors.HexColor('#12305c')), Spacer(1,9),
  P("<b>Purpose.</b> This note explains, in plain language and with the governing equations, the two channel-"
    "state-information (CSI) regimes used in the cell-free hybrid near-field / far-field (NF/FF) XL-MIMO uplink: "
    "(i) <b>perfect CSI</b> (a genie upper bound) and (ii) <b>estimated CSI</b> (pilot-based LMMSE). For each, it "
    "states how the channel is obtained, why that choice is made, who established the method, and which reference "
    "supports it &#8212; separately for the far-field and the near-field case. It also assesses what the current "
    "code actually does and flags one inconsistency to correct in the new paper.", BODY),
]

# 1. Channel model
st += [P("1. The channel that CSI refers to", H1),
  P("Each AP l has a uniform linear array of N antennas at half-wavelength spacing (Delta = lambda/2). Element n "
    "of AP l sits at position r_{n,l}. A single-antenna user k at position u_k has, to AP l, the channel vector "
    "h_{l,k} in C^N. Two regimes are distinguished by the AP&#8211;user distance d_{l,k} relative to the Rayleigh "
    "(Fraunhofer) distance d_Ray = 2 D^2 / lambda, with aperture D = (N-1) lambda/2:", BODY),
  P("&#8226; Far field (d_{l,k} &gt; d_Ray): planar wavefront. The channel is modelled as correlated Rayleigh "
    "fading,", BODY),
  P("h_{l,k} ~ CN( 0 , R_{l,k} ),", EQ),
  P("where R_{l,k} in C^{NxN} is the spatial correlation matrix set by the large-scale fading beta_{l,k} and the "
    "angular spread (the model of Demir&#8211;Bjornson&#8211;Sanguinetti [2]).", BODY),
  P("&#8226; Near field (d_{l,k} &lt;= d_Ray): spherical wavefront. Using the non-uniform spherical-wave (NUSW) "
    "model, the line-of-sight channel is deterministic given the geometry:", BODY),
  P("[h_{l,k}]_n = sqrt(beta_{l,k}) * ( d_{l,k} / d_{n,l,k} ) * exp( -j 2 pi d_{n,l,k} / lambda ),", EQ),
  P("with d_{n,l,k} = || r_{n,l} - u_k || the exact element-to-user distance. The amplitude term d/d_n and the "
    "non-linear phase are what distinguish the near field from the plane-wave model [4,5].", BODY),
  P("<b>Why this matters for CSI.</b> Perfect CSI means the receiver knows h_{l,k} exactly; estimated CSI means "
    "it forms an approximation h&#770;_{l,k} from pilots. The two cases are computed very differently in the near "
    "field than in the far field, which is the whole subject of this note.", BODY),
]

# 2. Perfect CSI
st += [P("2. Perfect CSI (the genie upper bound)", H1),
  P("<b>What it is.</b> Perfect CSI assumes the receiver knows the true channel realization h_{l,k} for every "
    "AP&#8211;user pair, with no estimation error. It is not &#8216;computed&#8217; by an algorithm: the exact "
    "channel vector generated by the model is fed directly into the combiner. In the code this means the local "
    "MMSE filter", BODY),
  P("v_{l,k} = sqrt(P) ( P sum_j h_{l,j} h_{l,j}^H + sigma^2 I )^{-1} h_{l,k}", EQ),
  P("is built from the TRUE h_{l,k} (not an estimate), and the estimation-error covariance is set to zero, "
    "C_{l,k} = 0.", BODY),
  P("<b>How it is obtained in each regime.</b>", BODY),
  bullets([
    "<b>Far field:</b> the true realization h_{l,k} ~ CN(0,R_{l,k}) drawn by the simulator is used as-is.",
    "<b>Near field (cell-free):</b> the true NUSW vector above &#8212; computed from the exact AP&#8211;element "
    "and user coordinates &#8212; is used as-is. Because the near-field LoS channel is deterministic given "
    "geometry, &#8216;perfect near-field CSI&#8217; is simply this exact spherical-wave vector.",
  ]),
  P("<b>Why it is used.</b> Perfect CSI is the standard genie benchmark in detection papers: it isolates the "
    "detector&#8217;s performance from channel-acquisition error and gives an achievable upper bound. Any "
    "estimated-CSI curve must lie below it. This is exactly the assumption your submitted paper states "
    "(&#8216;Perfect CSI is assumed at the receiver&#8217;).", BODY),
  P("<b>Who does it this way / references.</b> Perfect-CSI benchmarking is universal; for cell-free specifically "
    "it is the upper-bound convention in Demir&#8211;Bjornson&#8211;Sanguinetti [2] and Bjornson&#8211;Sanguinetti "
    "[3]. It is a modelling assumption, not a contribution, so it needs only a standard citation.", BODY),
  NOTE.clone('n2') if False else P("<b>Honest point.</b> Perfect CSI is legitimately an UPPER BOUND: it is "
    "guaranteed to be at least as good as any estimated-CSI result, because the estimated case runs the identical "
    "detector on a noisier channel with a positive error covariance. So &#8216;perfect &gt;= estimated&#8217; is "
    "physics, not a tuned outcome.", NOTE),
]

# 3. Estimated CSI - far field
st += [PageBreak(), P("3. Estimated CSI: the far-field (LMMSE) case", H1),
  P("<b>Idea in one line.</b> Users send known pilot sequences; each AP correlates the received pilot with the "
    "known sequence and applies a linear minimum-mean-square-error (LMMSE) filter to turn the noisy observation "
    "into the best linear estimate of the channel.", BODY),
  P("<b>Step 1 &#8212; pilot transmission.</b> With tau_p orthogonal pilots and a pilot-sharing set P_t (users "
    "assigned pilot t), the de-spread pilot observation at AP l is", BODY),
  P("y^p_{l,t} = sqrt(P) tau_p sum_{i in P_t} h_{l,i} + n_{l,t},   n_{l,t} ~ CN(0, sigma^2 tau_p I).", EQ),
  P("<b>Step 2 &#8212; LMMSE estimate.</b> The MMSE estimate of h_{l,k} and its error covariance are", BODY),
  P("h&#770;_{l,k} = sqrt(P) R_{l,k} Psi_{l,t}^{-1} y^p_{l,t},", EQ),
  P("Psi_{l,t} = P tau_p sum_{i in P_t} R_{l,i} + sigma^2 I,", EQ),
  P("C_{l,k} = R_{l,k} - P tau_p R_{l,k} Psi_{l,t}^{-1} R_{l,k}.", EQ),
  P("By the orthogonality principle h_{l,k} = h&#770;_{l,k} + h&#771;_{l,k} with the estimate h&#770; and the error "
    "h&#771; uncorrelated; C_{l,k} is the covariance of h&#771;. The detector then uses h&#770; in place of h in the "
    "combiner and carries C as the residual-error term.", BODY),
  P("<b>Pilot contamination.</b> When several users share pilot t (tau_p &lt; K), Psi_{l,t} contains their "
    "correlation matrices, so the estimate of user k is corrupted by its co-pilot users. This is the fundamental "
    "limit of pilot-based estimation and is what your code&#8217;s contamination mode (tau_p = K/2) reproduces.", BODY),
  P("<b>Quality metric.</b> The normalized MSE, NMSE_{l,k} = tr(C_{l,k}) / tr(R_{l,k}), quantifies estimation "
    "quality (0 = perfect, 1 = useless) and is the natural new figure for the estimation paper.", BODY),
  P("<b>Who does it this way / references.</b> This is the canonical cell-free LMMSE channel estimator of "
    "Demir&#8211;Bjornson&#8211;Sanguinetti [2, Sec. 4] and Bjornson&#8211;Sanguinetti [3]; the pilot-contamination "
    "treatment traces to Ngo et al. [1]. <b>This part of your code is correct and standard for the far field.</b>", BODY),
]

# 4. Estimated CSI - near field
st += [P("4. Estimated CSI: the near-field case (the real work)", H1),
  P("<b>Why the far-field estimator is not enough.</b> The LMMSE formulas above are correct for any correlation "
    "matrix R, but two near-field-specific issues arise:", BODY),
  bullets([
    "<b>The angular (DFT) basis fails.</b> In the far field a channel is sparse in the angle domain, which "
    "classical estimators exploit. In the near field the response depends on angle AND distance, so the DFT "
    "codebook is mismatched and angle-only estimators lose accuracy [5].",
    "<b>The correlation model changes.</b> R_{l,k} must be the near-field (spherical-wave) correlation, not the "
    "planar one; using the wrong R biases the estimate.",
  ]),
  P("<b>The established near-field solution &#8212; polar-domain estimation.</b> Cui and Dai [5] introduce a "
    "polar-domain dictionary sampled jointly in angle and distance; the near-field channel is sparse in that "
    "dictionary, and it is recovered by sparse methods (e.g., polar-domain orthogonal matching pursuit). This is "
    "the reference method for &#8216;near-field vs far-field channel estimation&#8217; and the natural estimator "
    "to cite and adapt. A lower-complexity route that fits your code with minimal change is the "
    "<b>near-field LMMSE</b>: keep equations (Step 2) but build R_{l,k} from the NUSW model so the estimator is "
    "matched to the spherical wavefront.", BODY),
  P("<b>For the cell-free / distributed setting specifically</b>, near-field combining and processing for "
    "cell-free XL-MIMO is treated by Wang et al. [6] and the hybrid near-far activity-detection line by Lei et "
    "al. [7]; the near-field XL-MIMO tutorial [4] is the umbrella reference. These justify doing estimation with "
    "a near-field-aware model in a cell-free system.", BODY),
  NOTE.clone('n4') if False else P("<b>What your current code does &#8212; and the gap.</b> For far-field links "
    "the code runs the correct LMMSE estimator (Section 3). For near-field links it does NOT estimate: it "
    "overwrites the estimate with the exact deterministic steering vector (h&#770; = h_det_all) and assigns a tiny "
    "floor covariance. That is effectively PERFECT near-field CSI, not estimation. So today the code is a hybrid: "
    "near-field perfect + far-field estimated. For a channel-estimation paper this must be replaced by a real "
    "near-field estimator (polar-domain [5] or NF-LMMSE) so that the estimated-CSI curves are genuinely "
    "estimated everywhere.", NOTE),
]

# 5. Recommendation & experiment
st += [P("5. What to compute for the new paper", H1),
  P("Present three clearly-separated CSI settings, run on the SAME system and the SAME (cross-AP) detector, so "
    "the paper is about acquisition, not detection:", BODY),
  bullets([
    "<b>Perfect CSI (upper bound):</b> true h_{l,k} everywhere, C = 0. Matches your paper&#8217;s stated "
    "assumption; use it as the benchmark curve. [2,3]",
    "<b>Proposed near-field estimation:</b> LMMSE with the NUSW correlation (or polar-domain OMP) on near-field "
    "links, standard LMMSE on far-field links, with realistic pilots and contamination. [5,6]",
    "<b>Far-field-mismatched estimation (ablation):</b> estimate every link with the planar-wave model. This "
    "shows the penalty of ignoring the near field and is the ablation that motivates the proposed estimator. [5]",
  ]),
  P("<b>New figures:</b> (a) NMSE vs SNR for the three settings; (b) BER and per-user SE vs SNR, perfect vs "
    "proposed vs mismatched, for the four detectors including cross-AP list-SIC and the IDD scheme; (c) BER vs "
    "pilot length tau_p (the estimation-overhead / contamination trade-off). Together these make &#8216;impact of "
    "channel estimation&#8217; a self-contained contribution distinct from the detector paper.", BODY),
  P("<b>Fix the inconsistency.</b> The submitted paper states perfect CSI but its results were generated with "
    "estimated (contaminated) CSI. In the new paper, label every curve by its CSI setting and generate perfect-CSI "
    "curves with C = 0 and true h, and estimated-CSI curves with the estimator of Section 4 &#8212; never mix.", BODY),
]

# 6. Reference mapping
st += [P("6. Which reference justifies which choice", H1)]
rt=[["Choice","Why / where it comes from","Ref"],
 ["Perfect CSI as an upper bound","Standard genie benchmark isolating the detector from acquisition error","[2],[3]"],
 ["Far-field LMMSE pilot estimation","Canonical cell-free estimator; error covariance C, orthogonality","[2],[3]"],
 ["Pilot contamination model","Co-pilot sharing corrupts the estimate; fundamental limit","[1],[2]"],
 ["Near-field NUSW channel model","Spherical wavefront, amplitude+non-linear phase across the array","[4],[5]"],
 ["Near-field (polar-domain) estimation","Angle-distance dictionary; DFT basis fails in the near field","[5]"],
 ["Near-field cell-free processing","Distributed combining / estimation in near-field cell-free XL-MIMO","[6],[7]"],
]
T=Table(rt, colWidths=[4.3*cm, 8.4*cm, 1.6*cm])
T.setStyle(TableStyle([
  ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#22548c')),('TEXTCOLOR',(0,0),(-1,0),colors.white),
  ('FONTSIZE',(0,0),(-1,-1),8.4),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('VALIGN',(0,0),(-1,-1),'TOP'),
  ('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#aab7c8')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#eef3f9')]),
  ('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),
]))
st += [T, Spacer(1,8),
  P("<b>Bottom line.</b> Perfect CSI is a one-line assumption backed by [2,3] and is honest as an upper bound. "
    "Estimated CSI in the far field is the standard LMMSE of [2] &#8212; correct in your code. Estimated CSI in "
    "the near field is the genuine new work: adopt a near-field-aware estimator [5] (polar-domain or NUSW-LMMSE) "
    "in place of the current deterministic shortcut, and the perfect-vs-estimated comparison then becomes a "
    "defensible conference paper on channel estimation for near-field cell-free XL-MIMO.", BODY),
]

# References
st += [P("References", H1),
 P("[1] H. Q. Ngo, A. Ashikhmin, H. Yang, E. G. Larsson, T. L. Marzetta, &#8220;Cell-free massive MIMO versus "
   "small cells,&#8221; IEEE Trans. Wireless Commun., 16(3):1834&#8211;1850, 2017.", REF),
 P("[2] O. T. Demir, E. Bjornson, L. Sanguinetti, &#8220;Foundations of user-centric cell-free massive MIMO,&#8221; "
   "Found. Trends Signal Process., 14(3-4):162&#8211;472, 2021.", REF),
 P("[3] E. Bjornson, L. Sanguinetti, &#8220;Scalable cell-free massive MIMO systems,&#8221; IEEE Trans. Commun., "
   "68(7):4247&#8211;4261, 2020.", REF),
 P("[4] H. Lu, Y. Zeng, <i>et al.</i>, &#8220;A tutorial on near-field XL-MIMO communications toward 6G,&#8221; "
   "IEEE Commun. Surveys Tuts., 26(4):2213&#8211;2257, 2024.", REF),
 P("[5] M. Cui, L. Dai, &#8220;Channel estimation for extremely large-scale MIMO: Far-field or near-field?,&#8221; "
   "IEEE Trans. Commun., 70(4):2663&#8211;2677, 2022.", REF),
 P("[6] Z. Wang, J. Zhang, <i>et al.</i>, &#8220;Low-complexity distributed combining design for near-field "
   "cell-free XL-MIMO systems,&#8221; IEEE Trans. Wireless Commun., 25:11799&#8211;11815, 2026.", REF),
 P("[7] J. Lei, <i>et al.</i>, &#8220;A unified distributed algorithm for hybrid near-far field activity detection "
   "in cell-free massive MIMO,&#8221; IEEE Trans. Wireless Commun., 2026.", REF),
 Spacer(1,6),
 P("<i>Prepared as supporting material for the follow-up paper on channel estimation for the cross-AP list-based "
   "SIC / IDD near-field cell-free XL-MIMO system. Equations follow the model and notation of the submitted "
   "paper; references [1]&#8211;[7] are from that paper&#8217;s bibliography.</i>", CAP),
]

doc=SimpleDocTemplate(OUT, pagesize=A4, topMargin=1.6*cm, bottomMargin=1.5*cm, leftMargin=1.7*cm, rightMargin=1.7*cm,
                      title="Channel Estimation Note: Perfect vs Estimated CSI")
doc.build(st)
print("WROTE", OUT)

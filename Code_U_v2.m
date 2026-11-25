clear; close all; clc;
try, delete(gcp('nocreate')); catch, end   % close any parallel pool left open by a previous run (frees memory)
%% ====================================================================
%  CODE U v2 : MOVABLE (SLIDING-MIRROR) SPARSE ARRAY vs STATIC SPARSE ARRAYS
%           crowded near-field / hybrid-field cell-free uplink (from Code T)
%  --------------------------------------------------------------------
%  Same setting as Code T: L = 8 APs, N = 8 antennas per AP, K = 16 users,
%  16-QAM, uplink power control (SNR axis = received SNR per antenna at the
%  user's strongest AP), clusters within 10 dB, hotspot users. The sparse
%  arrays all have the same ~2.8 m aperture, so most AP-user links are near
%  field; the dense lambda/2 ULA (0.35 m) is far field by physics.
%
%  Arrays (each run WITHOUT CAPE = local per-AP estimation, and WITH CAPE):
%    dense ULA | nested | nested + difference co-array | coprime |
%    coprime + difference co-array | MISC | IMISC | EMISC |
%    movable array, uniform half | movable array, sparse (increasing-gap) half
%  MOVABLE ARRAY (Idea 1): two mirror halves at +u and -u slide outward by
%  lambda/4 per frame (nPilotBlk frames). Every mirror pair gives the lag 2u
%  with the near-field range term cancelled; over all frames the lags fill a
%  hole-free virtual array (angle candidates). The pilots of all frames are
%  also stacked into a synthetic aperture for the range and gain fit. Data
%  are detected with the array parked at frame sma_sig_d.
%  FAIRNESS: every scheme gets the SAME total pilot energy, that of one Code T
%  pilot block. Static arrays use it as in Code T; the movable array splits it
%  over its nPilotBlk frames (each frame 1/nPilotBlk of the energy), so it gets
%  no extra pilot energy. Static arrays therefore run exactly as in Code T.
%  MISC / IMISC / EMISC need their inter-element spacing sets (ies_misc,
%  ies_imisc, ies_emisc below); an empty set skips that array.
%  v2: many more random drops (smooth curves; the BER noise comes from the
%  drops, not from the symbols), the run size per drop of Code T; a faster
%  movable-array estimator (its stacked lambda/4 aperture has no grating
%  lobes, so one dictionary start + the two mirror-pair angles suffice, and 6
%  refinement rounds are 200x finer than the hit tolerance); CAPE ends with a
%  full-window refinement from its chosen lobe; a no-CAPE localisation
%  baseline (angle-of-arrival triangulation from the same APs and lobes); the
%  share of users located within 1 m. Figures V1-V9 (U1-U3, U12 geometry;
%  the older U4-U11, U13 with fig_old = true).
%  fast_mode = true  -> quick check;  fast_mode = false -> full run.
%  Results saved to codeU_v2_results.mat BEFORE plotting;
%  plot_codeU_v2_figures.m redraws the figures without rerunning.
%% ====================================================================

clear;
close all;
clc;
%% ====================================================================
%  CODE E : CODE D + DIFFERENCE CO-ARRAY AND VECTORIZATION (VIRTUAL ARRAY)
%  ------------------------------------------------------------------
%  Everything in Code D is kept unchanged. What is added here:
%   * 'nested' array geometry: the two-level nested array of Pal & Vaidyanathan
%     (IEEE TSP 2010), the array used by Pinto & de Lamare (arXiv:2312.06033)
%     for super-resolution-array communications.
%   * Difference co-array + vectorization: the sample covariance of the raw
%     pilot snapshots is vectorized, duplicate lags are averaged, the
%     contiguous virtual ULA is extracted and spatially smoothed, and MUSIC on
%     that VIRTUAL array supplies high-resolution angle candidates. The range
%     is then found by a 1-D search at that angle on the physical observation,
%     which is the procedure used in the near-field co-array literature.
%   * eqN_study (CODE R): focused sparse-array comparison at 32 antennas.
%     Six schemes: Dense ULA + CAPE (ref), Nested local / Nested + CAPE,
%     Difference co-array + CAPE, Uniform sparse local / Uniform sparse + CAPE
%     (proposed). Two stories: (a) array geometry nested vs diff co-array vs
%     uniform, incl. performance vs number of users K; (b) CAPE on/off.
%     Figures R1..R9 (<=10, the important sparse-array metrics).
%   * idea_study (CODE F): applies the three NEW structure ideas from the
%     companion note (Idea 2 cross-AP co-array, Idea 3 sum+difference
%     co-array, Idea 4 hybrid-field reconfigurable) and compares them against
%     the ORIGINAL CAPE (ECA + CAPE), with figures CO-1 ... CO-5.
%
%  HONEST NOTES (please read before using the numbers)
%   1. A virtual sensor adds DEGREES OF FREEDOM and RESOLUTION, not received
%      energy. It cannot raise the detector's array gain by itself.
%   2. The standard difference co-array assumes PLANE waves. On near-field
%      links the curved wavefront breaks the clean virtual-ULA structure, so
%      the co-array angles here are approximate candidates that CAPE and the
%      local refinement then correct.
%   3. The covariance must be estimated from only tau = K pilot snapshots, so
%      it is noisy; this limits the achievable co-array gain.
%   4. With K = 8 users and N = 63 antennas per AP the O(N^2) degree-of-freedom
%      advantage is NOT the binding constraint here; the benefit to look for is
%      resolution (and antenna count), not source count.
%
%  ORIGINAL CODE D HEADER FOLLOWS.
%  ==================================================================
%  CODE D : SPARSE APs MADE USEFUL BY CROSS-AP POSITION-FUSED ESTIMATION
%           (hybrid near-/far-field cell-free XL-MIMO, cross-AP list SIC, no IDD)
%  ------------------------------------------------------------------
%  Built on the authentic-CSI code (Without_IDD_Code_Authentic_CSI.m) and
%  on the sparse-array infrastructure of Code C (array geometries,
%  geometry-aware polar dictionary, full 3 x 3 wrap-around).
%
%  WHY (what Code C showed)
%   * With the same number of antennas every array geometry has the same
%     AVERAGE inter-user interference (E[mu^2] ~ 1/N); a sparse array only
%     reduces the WORST case (co-angular / closely spaced users).
%   * With users spread over the area, the APs together already separate
%     the users (macro-diversity), so the sparse aperture adds nothing, and
%     its per-AP channel estimation is worse (narrow beams but ambiguous
%     lobes -> the local OMP picks wrong lobes on weak links). Net: the dense
%     ULA had the lower BER.
%   * With CROWDED users (hotspot) every AP sees the users inside a narrow
%     angular sector; there the large aperture is decisive (perfect-CSI BER
%     and SE clearly better), but the estimation penalty still ate the gain
%     of the cross-AP list SIC.
%
%  PROPOSED: CAPE = Cross-AP Position-fused Estimation for sparse APs
%   A user has ONE physical position and every AP's near-field channel is a
%   known function of it. The CPU fuses the pilot observations of the few
%   strongest APs in the POSITION domain (non-coherent fused likelihood),
%   which keeps only the lobe that is consistent across APs; every AP then
%   builds its atom from the fused position and refines it locally. A
%   consistency / likelihood gate keeps an AP's own local estimate when it is
%   on the same lobe and fits better, or fits better by more than the noise
%   level (so CAPE is never worse than local estimation at high SNR). The
%   sparse array's weakness (ambiguous lobes) is removed by the network, and
%   its strength (narrow beams -> sharp position fix) is kept.
%
%   Preliminary check (Octave, 6 setups, K = 8, BER of List+CrossAP,
%   estimated CSI; ULA local / ECA local / ECA + CAPE):
%     baseline drop 10 dB: 1.24e-2 / 1.72e-2 / 1.23e-2   (20 dB: all 0)
%     hotspot       10 dB: 9.2e-2  / 1.07e-1 / 8.0e-2
%     hotspot       20 dB: 1.28e-3 / 1.69e-3 / 5.5e-4   (perfect-CSI ECA: 4.2e-4)
%   Small sample: confirm with the full run of this code.
%   CAPE + sparse (ECA) APs + the cross-AP list SIC is the proposed system.
%   See csi_acquire_cape for the algorithm.
%
%  WHAT IS NEW IN THIS CODE (relative to Code C)
%   1. csi_fusion switch: NF links estimated with CAPE in every estimation
%      block of the main pipeline (old figures = proposed system).
%   2. user_scenario switch for the main pipeline: 'baseline' (original
%      drop), 'hotspot' (crowded users) or 'mixed' (half crowded).
%   3. Comparison study: 4 schemes x 3 scenarios on identical drops and
%      identical channel/noise draws:
%        ULA + local estimation (base code) | ECA + local estimation (Code C)
%        ULA + CAPE                         | ECA + CAPE (PROPOSED)
%      metrics: NMSE, localisation error, BER of the 4 detectors (estimated
%      and perfect CSI), sum and per-user SE, 5% outage SE, Jain fairness,
%      cross-AP gate rate, data and pilot fronthaul, information surplus,
%      estimation run time, and a user-crowding sweep (hotspot radius).
%   4. Figures R0..R9 (design + sparse/channel-estimation comparison).
%
%  NOVELTY BOUNDARY (state it this way in a paper): the ECA geometry is
%  from Zhou et al. (IEEE TWC 2025); location-based / near-field
%  localisation-aided estimation exists for single and modular arrays
%  (e.g. Kosasih et al. 2025). What is proposed here is cross-AP position
%  fusion used to RESOLVE SPARSE-ARRAY AMBIGUITIES in hybrid-field cell-free
%  XL-MIMO, combined with the cross-AP list SIC. Check recent literature
%  before writing "first".
%% ====================================================================
%% ====================================================================
%  INHERITED FROM CODE C (sparse-array infrastructure)
%   * array_type: AP array geometry ('ula' | 'usa' | 'eca' | 'eca_sa'); all
%     channel models, covariances, the polar dictionary and the Rayleigh
%     distance follow the true element positions. Default: 'eca' (63
%     elements, 7.6 m), the array of the proposed system.
%   * Geometry-aware polar dictionary (angular grid never coarser than the
%     aperture resolution or the dense-ULA grid; Fresnel-phase ring spacing).
%   * wrap_fix: full 3 x 3 wrap-around (the base code generated only the 3
%     diagonal images).
%   * Far-field covariance built from element positions with the same phase
%     convention as the NUSW near-field channel.
%% ====================================================================

%% ====================================================================
%  PARAMETERS
%% ====================================================================
L  = 8;
K  = 8 ;

fc     = 3e9;
c0     = 3e8;
lambda = c0 / fc;
d_ant  = lambda / 2;
d_H    = 0.5;
fc_GHz = fc / 1e9;

%% ====================================================================
%  [CODE C] AP ARRAY GEOMETRY  (non-uniform sparse array)
%  ------------------------------------------------------------------
%  array_type selects the linear array mounted at EVERY cell-free AP:
%   'ula'    : dense half-wavelength ULA with N_ref elements (base-code array).
%   'usa'    : uniform sparse array, N_ref elements, spacing usa_eta*lambda/2.
%              Grating lobes at sin(th_u) = sin(th_0) + 2u/usa_eta (Zhou,
%              Lemma 5), refocused at r_u = r_0 cos^2(th_u)/cos^2(th_0)
%              (Zhou, Lemma 6).
%   'eca'    : extended coprime array (Zhou et al.), NON-UNIFORM. eca_Lc basic
%              coprime arrays; each is an eca_M-element sparse ULA with
%              spacing eca_N*d0 plus an eca_N-element sparse ULA with spacing
%              eca_M*d0 sharing the first element; period eca_M*eca_N*d0;
%              upper and lower halves mirrored about the shared origin.
%              Q = eca_Lc*(eca_M+eca_N-1) - 1 elements, aperture
%              (eca_Lc*eca_M - 2)*eca_N*d0. Default 8/5/4: 63 elements, 7.6 m.
%   'eca_sa' : the same ECA with eca_Lc_sa basic arrays: about the SAME
%              APERTURE as the dense ULA (3.6 m) with 31 elements. Separates
%              the element-count effect from the aperture effect.
%  All geometries lie on the lambda/2 grid, so sin(theta) = -1 and +1 still
%  alias exactly in the far field and the estimator's end-fire wrap is valid.
%% ====================================================================
array_type = 'ula';   % [CODE O] PART 1 uses the DENSE array
N_ref      = 64;     % [CODE O] 64 antennas per AP in PART 1 (dense)
usa_eta    = 2;      % USA spacing in units of lambda/2
eca_M      = 5;      % coprime pair (M, N), M > N, gcd = 1
eca_N      = 4;
eca_Lc     = 8;      % number of basic coprime arrays (even)
eca_Lc_sa  = 4;      % basic arrays of the same-aperture ECA
wrap_fix   = true;   % full 3 x 3 wrap-around (see header)

%% ====================================================================
%  [CODE E] NESTED ARRAY AND DIFFERENCE CO-ARRAY (VIRTUAL ARRAY)
%  ------------------------------------------------------------------
%  Two-level nested array (Pal & Vaidyanathan, IEEE TSP 2010): an inner ULA of
%  nest_N1 elements at spacing d0 = lambda/2 and an outer ULA of nest_N2
%  elements at spacing (nest_N1 + 1) d0. Its difference co-array is a HOLE-FREE
%  virtual ULA with 2*nest_N2*(nest_N1+1) - 1 virtual sensors, i.e. about
%  N^2/2 virtual sensors from N = nest_N1 + nest_N2 physical ones.
%  The default (11, 12) gives 23 antennas over 7.2 m: about the SAME aperture
%  as the 63-element ECA with roughly one third of the antennas, which is the
%  antenna-saving comparison made by Pinto & de Lamare (arXiv:2312.06033).
%% ====================================================================
nest_N1      = 11;       % inner (dense) elements, spacing d0
nest_N2      = 12;       % outer (sparse) elements, spacing (nest_N1+1) d0
ca_use       = true;     % use the difference co-array to propose candidates
ca_nAngGrid  = 721;      % angle grid of the co-array MUSIC spectrum
ca_nRangeGrid = 60;      % inverse-range grid of the 1-D range search
ca_nPeak     = 3;        % co-array angle candidates kept per proposing AP
ca_Nsub_max  = 80;       % cap on the spatial-smoothing subarray length (speed)
APG  = array_geometry(array_type, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa, nest_N1, nest_N2);
N    = APG.N;
N_BS = L * N;

ASD_varphi = deg2rad(15);
ASD_theta  = deg2rad(15);

tau_c       = 200;
prelog_pcsi = 1;

SNR_dB  = 0:5:30;
SNR_lin = 10.^(SNR_dB / 10);
nSNR    = length(SNR_dB);
LN      = L * N;

squareLen = 300;
hDiff     = 10;
hBS       = 25;
hUT       = 1.5;

D_ap  = APG.aperture;               % [CODE C] aperture of the array in use
d_Ray = 2 * D_ap^2 / lambda;

% --- Centralized BS panel: derive a near-square Nx_BS x Ny_BS grid whose
%     element count EXACTLY equals N_BS (= L*N). This must match N_BS or the
%     R_BS covariance assignment fails (was hardcoded 32x16=512 before).
Ny_BS = floor(sqrt(N_BS));
while mod(N_BS, Ny_BS) ~= 0
    Ny_BS = Ny_BS - 1;
end
Nx_BS = N_BS / Ny_BS;
[gx_BS, gy_BS] = meshgrid(0:Nx_BS - 1, 0:Ny_BS - 1);
px_BS = gx_BS(:) - (Nx_BS - 1) / 2;
py_BS = gy_BS(:) - (Ny_BS - 1) / 2;
assert(numel(px_BS) == N_BS, 'BS panel elements (%d) must equal N_BS (%d).', numel(px_BS), N_BS);
S_ang = 300;
D_BS    = sqrt((Nx_BS - 1)^2 + (Ny_BS - 1)^2) * d_ant;
d_Ray_BS = 2 * D_BS^2 / lambda;

eps_NF = 0.01;
n_idx  = (-(N - 1) / 2:(N - 1) / 2)';

B            = 20e6;
noisePow_dBm = -174 + 10 * log10(B) + 15;

n_collinear_pairs = floor(K / 4);
r_min_NF = 20;
r_max_NF = min(0.9 * d_Ray, squareLen / 2);

squareLen_sweep = [100 200 300 500 800 1200];   % [CODE C] extended: FF links return for the large sparse aperture
nSweep          = length(squareLen_sweep);
nSetups_sweep   = 5;
% Reference SNR for the area sweep, the cluster size sweep and the CDF.
% MUST NOT use find(SNR_dB == 10, 1): if 10 is not an exact grid point that
% returns empty, p_ref becomes empty, and the area sweep dies with a
% misleading "Incorrect dimensions for matrix multiplication". Worse, the
% cluster size sweep gates on si == snr_ref_idx and would silently produce
% all zero curves. Nearest grid point is always well defined.
[~, snr_ref_idx] = min(abs(SNR_dB - 10));
snr_ref_dB = SNR_dB(snr_ref_idx);

%% ====================================================================
%  MASTER SWITCH: CONTAMINATION MODE   (see warnings W2, W3)
%  ------------------------------------------------------------------
%  false : LEGACY. tau_p = 7, tau_sym = K, tau_fig2 = K. Old figures 1 to
%          8 reproduce exactly. Pilots are effectively orthogonal at the
%          symbol level, so the proposed clustering rule reduces to
%          channel norm and the new curves coincide with the baseline.
%          This is the correct no harm behaviour, not a failure.
%  true  : STUDY. tau_p = tau_sym = tau_fig2 = K/2 with aligned pilots.
%          The only regime in which clustering can be studied. This DOES
%          change the old figures, which now share one operating point.
%% ====================================================================
contamination_mode = true;
nf_csi_consistent  = true;   % (S2) LEGACY ONLY (genie_csi = true). NF pairs: exact channel ->
                              % NUSW-floor error covariance. Has NO effect when genie_csi = false.

%% ====================================================================
%  CSI ACQUISITION MODE
%  ------------------------------------------------------------------
%  genie_csi = false : AUTHENTIC (default). No true-channel information
%      enters any estimator. The true channel is used ONLY to synthesise the
%      received pilot signal, exactly as nature would.
%        * Orthogonal pilots, tau_p = tau_sym = tau_fig2 = K (paper Sec. II-B:
%          "pilot phase of length tau_p >= K with mutually orthogonal pilots").
%        * Near-field pairs (A_lk = 1): polar-domain OMP over a polar
%          dictionary (angles uniform in sin(theta), inverse range 1/r uniform,
%          r -> inf = far-field atom) with gridless refinement of each selected
%          atom (Cui & Dai, IEEE TCOM 2022, P-SOMP / P-SIGW). The path gain is
%          LMMSE-scaled using only the large-scale coefficient beta_lk.
%        * Far-field pairs (A_lk = 0): LMMSE from the statistical covariance
%          R_lk (second-order statistics, as in Bjornson & Sanguinetti 2020).
%        * The combiners use a DATA-DRIVEN error covariance sigma_e^2 * I,
%          computed from the estimator's own residual and the known noise
%          variance (never from the true channel).
%        * A_lk (near/far-field indicator, eq. (3) of the paper) is a
%          large-scale quantity, assumed known at the CPU like beta_lk.
%  genie_csi = true  : LEGACY. Reproduces the previous figures exactly: NF
%      "estimates" come from an LMMSE whose prior R2 = h*h' + eps*beta*I is
%      built from the TRUE channel, NF links are then overwritten by the true
%      channel, and pilots follow contamination_mode. Not a real estimator.
%% ====================================================================
genie_csi = false;

% Polar-domain dictionary / OMP settings (authentic mode only)
pd_os_ang = 2;       % angular oversampling: 2N grid points in sin(theta)
pd_q_max  = 1 / 8;   % largest inverse range in the dictionary [1/m] (r >= 8 m)
pd_n_q    = 13;      % inverse-range rings, including q = 0 (far-field atom)
pd_smax   = 2;       % maximum number of atoms per link
pd_stop   = 1.3;     % stop when ||residual||^2 <= pd_stop * N * sigma_w^2
pd_nref   = 10;      % gridless refinement rounds per selected atom
pd_nstart = 3;       % refinement starts from the 3 strongest grid atoms
pd_phase_step = 1.6; % [CODE C] ring spacing: Fresnel phase change across the
                     % half-aperture between rings [rad]; rings are added only
                     % when the aperture needs more than pd_n_q of them
pd_nS_min = pd_os_ang * N_ref;  % [CODE C] never a coarser angular grid than the dense
                     % ULA of the base code (128 angles): with few elements
                     % (eca_sa) a grid at the bare aperture resolution makes
                     % the weak-link fallback unreliable at low SNR

if ~genie_csi
    % AUTHENTIC: mutually orthogonal pilots (paper Sec. II-B). Polar-domain
    % OMP cannot tell apart two users that share a pilot without genie side
    % information, so pilot sharing is not used in this mode.
    tau_p     = K;
    tau_sym   = K;
    tau_fig2  = K;
    pilot_mode = 'orthogonal';
elseif contamination_mode
    tau_p     = K / 2;
    tau_sym   = K / 2;
    tau_fig2  = K / 2;
    pilot_mode = 'aligned';     % 'aligned' or 'random'  (see W3)
else
    tau_p     = 7;
    tau_sym   = K;
    tau_fig2  = K;
    pilot_mode = 'legacy';
end
prelog = 1 - tau_p / tau_c;

% ---- startup sanity checks -----------------------------------------
% Non integer or out of range pilot lengths produce failures far from
% their cause: randn(...,tau_p) errors, mod(...,tau_p) yields fractional
% pilot indices, and find(pilotIndex == t) returns empty so entire pilot
% groups are silently skipped. Catch it here instead.
tau_p    = max(1, round(tau_p));
tau_sym  = max(1, round(tau_sym));
tau_fig2 = max(1, round(tau_fig2));
assert(tau_p <= K && tau_sym <= K && tau_fig2 <= K, ...
       'Pilot length exceeds K; there would be unused pilots.');
assert(~isempty(snr_ref_idx), 'snr_ref_idx is empty; check SNR_dB.');
assert(mod(K, 2) == 0, ...
       'K must be even for the collinear pair placement and aligned pilots.');


%% ====================================================================
%  CLUSTERING PARAMETERS
%% ====================================================================
eta_FH   = 0;
clu_div  = 2 ;   % Threshold = max/clu_div. Smaller clu_div gives SMALLER
                 % clusters, which is the fronthaul constrained regime in
                 % which both the proposed clustering and the Cross-AP
                 % detector are actually motivated. Set to 10 to restore
                 % the old near full size clusters.
rho2_OCL = 0;

%% ====================================================================
%  [CODE D] PROPOSED ESTIMATION (CAPE) AND USER SCENARIO
%% ====================================================================
csi_fusion   = true;     % true: NF links by CAPE (proposed); false: per-AP polar OMP (Code C)
cape_nAP     = 4;        % APs whose pilot observations are fused (strongest NF APs of the user)
cape_nCand   = 3;        % APs that propose candidate positions
cape_J       = 3;        % local candidates (lobes) per proposing AP
cape_step0   = 1;        % [m] initial step of the 2-D position refinement
cape_nIter   = 14;       % refinement iterations (step halves on failure)
cape_consist = 0.5;      % consistency gate: keep the local atom if |b_loc^H b_fused| >= cape_consist and it fits better
user_scenario = 'baseline';   % main pipeline drop: 'baseline' | 'hotspot' | 'mixed'
hs_r         = 25;       % [m] hotspot radius


%% ====================================================================
%  SYMBOL-LEVEL LIST DETECTION PARAMETERS
%% ====================================================================
nSym    = 400;
M_lst   = 3;
d_th    = 0.60;
maxBr   = 8;
modOrder = 16;
modOrder2 = 16;


%% ====================================================================
%  RUNTIME MODE
%% ====================================================================
fast_mode = false;
if fast_mode
    nSetups = 3;
    nReal   = 5;
    nSym    = 150;
    nSym2   = 400;
else
    nSetups = 120;
    nReal   = 10;
    nSym    = 400;
    nSym2   = 1200;
end

fprintf('=== CF-mMIMO: FF-only vs NF/FF vs Centralized + Clustering + List ===\n');
fprintf('L=%d  K=%d  N=%d  N_BS=%d  LN/K=%.0f\n', L, K, N, N_BS, LN / K);
fprintf('[CODE C] AP array: %s, %d elements, aperture %.2f m (sparsity %.2f), Rayleigh distance %.1f m\n', ...
        array_type, APG.N, APG.aperture, APG.sparsity, d_Ray);
fprintf('[CODE D] estimation: %s | user scenario (main pipeline): %s\n', ...
        ternary_str(csi_fusion, sprintf('CAPE (fuse %d APs, %d proposers x %d lobes)', cape_nAP, cape_nCand, cape_J), 'per-AP polar OMP'), user_scenario);
fprintf('d_Ray(CF)=%.2fm  d_Ray(BS)=%.0fm  squareLen=%dm\n', d_Ray, d_Ray_BS, squareLen);
fprintf('genie_csi=%d  contamination_mode=%d  pilot_mode=%s  clu_div=%d\n', ...
        genie_csi, contamination_mode, pilot_mode, clu_div);
if ~genie_csi
    fprintf(['CSI: AUTHENTIC. NF pairs -> polar-domain OMP (%d atoms max, %d refinement rounds);\n' ...
             '     FF pairs -> statistical LMMSE; orthogonal pilots tau_p=%d.\n'], pd_smax, pd_nref, tau_p);
else
    fprintf('CSI: LEGACY genie path (true-channel prior / override). Not a real estimator.\n');
end
fprintf('tau_p=%d  tau_sym=%d  tau_fig2=%d  prelog=%.3f\n', ...
        tau_p, tau_sym, tau_fig2, prelog);
fprintf('List: nSym=%d M=%d d_th=%.2f\n', nSym, M_lst, d_th);
fprintf('Reference SNR for area/cluster-size sweeps and CDF: %d dB (index %d)\n', ...
        snr_ref_dB, snr_ref_idx);
fprintf('%d setups x %d MC x %d SNR pts\n', nSetups, nReal, nSNR);
if genie_csi && ~contamination_mode
    fprintf(['\nNOTE: contamination_mode is OFF. Symbol level pilots are\n' ...
             'orthogonal, so the proposed rule reduces to channel norm by\n' ...
             'construction and its curves WILL coincide with the baseline.\n' ...
             'Set contamination_mode = true to study clustering.\n']);
end
fprintf('\n');

%% ====================================================================
%  ACCUMULATORS
%% ====================================================================
SR1_lmmse_acc = zeros(nSetups, nSNR);
BER1_lmmse_acc = zeros(nSetups, nSNR);
SR1_cnsic_acc = zeros(nSetups, nSNR);
BER1_cnsic_acc = zeros(nSetups, nSNR);

SR2_lmmse_acc = zeros(nSetups, nSNR);
BER2_lmmse_acc = zeros(nSetups, nSNR);
SR2_cnsic_acc = zeros(nSetups, nSNR);
BER2_cnsic_acc = zeros(nSetups, nSNR);

BER1_NF_acc = zeros(nSetups, nSNR);
SR1_NF_acc = zeros(nSetups, nSNR);
BER2_NF_acc = zeros(nSetups, nSNR);
SR2_NF_acc = zeros(nSetups, nSNR);
BER2_FF_acc = zeros(nSetups, nSNR);
SR2_FF_acc = zeros(nSetups, nSNR);
BER1_CL_acc = zeros(nSetups, nSNR);
SR1_CL_acc = zeros(nSetups, nSNR);
BER2_CL_acc = zeros(nSetups, nSNR);
SR2_CL_acc = zeros(nSetups, nSNR);

SR3_lmmse_acc = zeros(nSetups, nSNR);
BER3_lmmse_acc = zeros(nSetups, nSNR);
SR3_cnsic_acc = zeros(nSetups, nSNR);
BER3_cnsic_acc = zeros(nSetups, nSNR);

SR_CNcl_lmmse_acc = zeros(nSetups, nSNR);
BER_CNcl_lmmse_acc = zeros(nSetups, nSNR);
SR_CNcl_cnsic_acc = zeros(nSetups, nSNR);
BER_CNcl_cnsic_acc = zeros(nSetups, nSNR);
% [NEW, DI RENNA-STYLE ANALYTICAL RATE] accumulators for the cluster-vs
% cross-AP achievable sum-rate gap. R_clu fuses over the serving cluster
% only; R_all fuses over ALL APs (the observations cross-AP reads back in).
% Both are prelog*log2(1+Gamma_k) with Gamma_k the MMSE-combiner SINR from
% det_local_sic, averaged over realizations, following Di Renna &
% de Lamare (IEEE T-Comms 2020), Thm 1 / eqs. (31)-(33). deltaR = R_all -
% R_clu is the analytical rate value of the discarded APs. It is ZERO under
% full cooperation (Theorem 1), which is the correct reduction property.
SR_anaClu_acc = zeros(nSetups, nSNR);   % cluster-only analytical sum-rate
SR_anaAll_acc = zeros(nSetups, nSNR);   % all-AP (cross-AP) analytical sum-rate
etaList_acc = zeros(nSetups, nSNR);     % [NEW] measured List-SIC trigger rate
etaXap_acc  = zeros(nSetups, nSNR);     % [NEW] measured Cross-AP invocation rate
SR_RTcl_lmmse_acc = zeros(nSetups, nSNR);
BER_RTcl_lmmse_acc = zeros(nSetups, nSNR);
SR_RTcl_cnsic_acc = zeros(nSetups, nSNR);
BER_RTcl_cnsic_acc = zeros(nSetups, nSNR);

% [PERFECT-CSI, FIG 1/2] genie variants (true channel in the combiner, C=0)
% of the five architectures, L-MMSE and SIC, for the perfect-CSI sum-rate
% (Fig 1) and empirical BER (Fig 2) figures.
SR1_lmmse_pf_acc = zeros(nSetups, nSNR);   SR1_cnsic_pf_acc = zeros(nSetups, nSNR);
SR2_lmmse_pf_acc = zeros(nSetups, nSNR);   SR2_cnsic_pf_acc = zeros(nSetups, nSNR);
SR3_lmmse_pf_acc = zeros(nSetups, nSNR);   SR3_cnsic_pf_acc = zeros(nSetups, nSNR);
SR_CNcl_lmmse_pf_acc = zeros(nSetups, nSNR);  SR_CNcl_cnsic_pf_acc = zeros(nSetups, nSNR);
SR_RTcl_lmmse_pf_acc = zeros(nSetups, nSNR);  SR_RTcl_cnsic_pf_acc = zeros(nSetups, nSNR);


loadCN_acc = zeros(nSetups, nSNR);
loadBSR_acc = zeros(nSetups, nSNR);
clDiff_acc = zeros(nSetups, nSNR);

% rows: [Linear, Hard-SIC, List-SIC, List+CrossAP(proposed)]
eBERhf = zeros(4, nSetups, nSNR);
eBERcn = zeros(4, nSetups, nSNR);
eBERrt = zeros(4, nSetups, nSNR);

% [FIG 2 EMPIRICAL] rows [Linear, Hard-SIC] per architecture
eBER2_ff = zeros(2, nSetups, nSNR);
eBER2_hf = zeros(2, nSetups, nSNR);
eBER2_bs = zeros(2, nSetups, nSNR);
eBER2_cn = zeros(2, nSetups, nSNR);
eBER2_rt = zeros(2, nSetups, nSNR);
% [FIG 2 EMPIRICAL, PERFECT CSI] genie variants (true channel, C=0)
eBER2_ff_pf = zeros(2, nSetups, nSNR);
eBER2_hf_pf = zeros(2, nSetups, nSNR);
eBER2_bs_pf = zeros(2, nSetups, nSNR);
eBER2_cn_pf = zeros(2, nSetups, nSNR);
eBER2_rt_pf = zeros(2, nSetups, nSNR);

% [FIGS 6/7/8] per-user BER & NMSE, 4 detectors, RATE clustering
mBERrt  = zeros(4, K, nSetups, nSNR);
mNMSErt = zeros(4, K, nSetups, nSNR);

% [CSI STUDY] perfect vs estimated CSI for the 4 detectors (rate cluster).
%  - "estimated", genie_csi = false (AUTHENTIC): near-field links by
%    polar-domain OMP, far-field links by statistical LMMSE, orthogonal pilots,
%    data-driven error covariance. No true-channel information is used.
%  - "estimated", genie_csi = true (LEGACY): near-field links use an LMMSE
%    whose prior R2 = h*h' + eps*beta*I is built from the TRUE channel h, so
%    the channel direction is known in advance. Kept only to reproduce old
%    figures; it must not be reported as channel estimation.
%  - "perfect": the true channel is used in the combiner with zero error
%    covariance (a genie upper bound). Feeds the two new summary figures.
%  Set csi_study = false to skip the extra passes and recover the old runtime.
csi_study = true;
eBERrt_est = zeros(4, nSetups, nSNR);      % estimated-CSI BER, 4 detectors (near-field/hybrid)
eBERrt_pf  = zeros(4, nSetups, nSNR);      % perfect-CSI   BER, 4 detectors (near-field/hybrid)
eBERhf_est = zeros(4, nSetups, nSNR);      % [CSI STUDY] full-hybrid (all APs serve) BER, estimated
eBERhf_pf  = zeros(4, nSetups, nSNR);      % [CSI STUDY] full-hybrid (all APs serve) BER, perfect
mBERrt_est = zeros(4, K, nSetups, nSNR);   % per-user BER (for sum-rate), estimated
mBERrt_pf  = zeros(4, K, nSetups, nSNR);   % per-user BER (for sum-rate), perfect
mNMSErt_est = zeros(4, K, nSetups, nSNR);  % per-user symbol NMSE, estimated (for achievable rate)
mNMSErt_pf  = zeros(4, K, nSetups, nSNR);  % per-user symbol NMSE, perfect

% [FAR-FIELD COMPARISON] same system + same rate cluster, but ALL links use the
% far-field (planar-wave) channel model and far-field LMMSE estimation. Lets us
% compare, for each detector, near-field vs far-field channel models under both
% perfect and estimated CSI. Set ff_study = false to skip.
ff_study = true;
eBERrt_ff_est = zeros(4, nSetups, nSNR);   eBERrt_ff_pf = zeros(4, nSetups, nSNR);
mBERrt_ff_est = zeros(4, K, nSetups, nSNR);  mBERrt_ff_pf = zeros(4, K, nSetups, nSNR);
mNMSErt_ff_est = zeros(4, K, nSetups, nSNR); mNMSErt_ff_pf = zeros(4, K, nSetups, nSNR);
% [ESTIMATOR NMSE] channel-estimation quality tr(C)/tr(R), averaged over links
nmse_est_acc = zeros(nSetups, nSNR);       % near-field/hybrid estimator NMSE
nmse_ff_acc  = zeros(nSetups, nSNR);       % far-field estimator NMSE
etaXap_pf_acc = zeros(nSetups, nSNR);      % Cross-AP invocation rate under PERFECT CSI

% Polar-domain dictionary, built once (depends only on the array geometry).
PD = build_polar_dict(APG.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);   % [CODE C]
PD.smax = pd_smax;  PD.stop = pd_stop;  PD.nref = pd_nref;  PD.nstart = pd_nstart;
if ~genie_csi
    fprintf('Polar dictionary: %d angles x %d inverse-range rings = %d atoms (N=%d)\n\n', ...
            PD.nS, PD.nQ, size(PD.W, 2), N);
end

nNF_total = 0;

%% ====================================================================
%  [CODE U] SETTINGS (the Code T crowded setting, unchanged) AND SCHEMES
%% ====================================================================
fast_mode = false;    % true: quick check;  false: full run (the script prints its run-time estimate)
if fast_mode
    nSetE = 28;   nSetK = 10;
else
    nSetE = 140;  nSetK = 50;     % random drops: SNR study / user-count sweep (about the cost of Code T)
end
nRealE = 2;  nSymE = 400;         % per drop, as Code T (the BER noise comes from the drops)
Ksw = [8 12 16 20 24];            % user-count sweep (schemes with CAPE)
ulaw_d_lambda = 4;                    % spacing that sets the common aperture (as Code T)
ulaw_eta      = 2 * ulaw_d_lambda;
ca_Nsub_max   = 1000;                 % co-array MUSIC: full smoothed virtual aperture
modOrder = 16;                        % 16-QAM (as Code T)
clu_div  = 10;                        % clusters: APs within 10 dB of the strongest AP
SNR_dB   = 0:4:28;                    % received SNR per antenna at the strongest AP [dB]
SNR_lin  = 10.^(SNR_dB / 10);
nSNR     = numel(SNR_dB);
[~, snr_ref_idx] = min(abs(SNR_dB - 10));
snr_ref_dB = SNR_dB(snr_ref_idx);
N0u  = 8;                             % antennas per AP
Ku   = 16;                            % users
Dtgt = (N0u - 1) * ulaw_eta * lambda / 2;   % common aperture of the sparse arrays (2.8 m)
nPilotBlk = 24;                       % frames of the movable array (they share ONE Code T pilot block's energy)
sma_sig_d = 7;                        % data frame of the movable array (frame 7: +-2 ... +-14 lambda)
tolS = 0.5 * lambda / Dtgt;           % angle hit tolerance in sin(theta): half a beamwidth of 2.8 m
% ---- MISC family: inter-element spacing (IES) sets for N0u = 8 antennas, in units of the grid step.
%      Copy them from the papers (N0u - 1 = 7 spacings each). An empty set skips that array.
ies_misc  = [];   % MISC : Zheng, Wang, Kong, Zhang, IEEE TSP 2019
ies_imisc = [];   % IMISC: arXiv 2208.09574
ies_emisc = [];   % EMISC: Sheng, Lu, Li, de Lamare, arXiv 2309.09044

eqN_study = true;
bs_name = {'Dense ULA', 'Nested', 'Nested + co-array', 'Coprime', 'Coprime + co-array', ...
           'MISC', 'IMISC', 'EMISC', 'Movable (uniform half)', 'Movable (sparse half)'};
bs_geom = {'ula', 'nested', 'nested', 'coprime', 'coprime', 'misc', 'imisc', 'emisc', 'sma_u', 'sma_s'};
bs_kind = {'plain', 'plain', 'coarray', 'plain', 'coarray', 'plain', 'plain', 'plain', 'sma', 'sma'};
bs_ies  = {[], [], [], [], [], ies_misc, ies_imisc, ies_emisc, [], []};
keepB = true(1, numel(bs_name));
for b = 6:8
    if isempty(bs_ies{b})
        keepB(b) = false;
        fprintf('NOTE: %s skipped: its spacing set ies_%s is empty (fill it in from the paper).\n', ...
                bs_name{b}, bs_geom{b});
    end
end
bs_name = bs_name(keepB);  bs_geom = bs_geom(keepB);  bs_kind = bs_kind(keepB);  bs_ies = bs_ies(keepB);
nB = numel(bs_name);
en_base = kron(1:nB, [1 1]);          % scheme -> base scheme
en_cape = repmat([0 1], 1, nB);       % 0 = local estimation, 1 = CAPE
nSchE   = numel(en_base);
en_names = cell(1, nSchE);
for sc = 1:nSchE
    en_names{sc} = sprintf('%s, %s', bs_name{en_base(sc)}, ternary_str(en_cape(sc), 'CAPE', 'local'));
end
en_caseName = {'Crowded (N = 8, K = 16)'};
en_N0 = N0u;  en_K = Ku;  nCaseE = 1;

% ---- geometry, dictionaries and estimation set-up of every base scheme ----
EN_geo = cell(1, nB);  BPD = cell(1, nB);  BES = cell(1, nB);
fprintf('\n=== [CODE U] arrays (N = %d antennas, common aperture %.2f m, movable array %d frames) ===\n', N0u, Dtgt, nPilotBlk);
for b = 1:nB
    G = geom_u(bs_geom{b}, N0u, lambda, Dtgt, nPilotBlk, sma_sig_d, bs_ies{b});
    PDb = build_polar_dict(G.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
    PDb.smax = pd_smax;  PDb.stop = pd_stop;  PDb.nref = pd_nref;  PDb.nstart = pd_nstart;
    ES = struct('kind', bs_kind{b}, 'F', nPilotBlk, 'PDe', PDb, 'Ne', G.N, 'mir', [], 'CA', [], 'J', cape_J);
    if strcmp(bs_kind{b}, 'sma')
        PDe = build_polar_dict(G.frames(:), lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
        PDe.smax = pd_smax;  PDe.stop = pd_stop;  PDe.nref = 6;  PDe.nstart = pd_nstart;   % [v2] 6 rounds
        ES.PDe = PDe;  ES.Ne = numel(G.frames);  ES.mir = G.mir;
        ES.J = 1;          % [v2] stacked lambda/4 aperture: no grating lobes, one dictionary start
    elseif strcmp(bs_kind{b}, 'coarray')
        ES.CA = coarray_setup(G.pos, lambda, 'diff', G.unit);
    end
    EN_geo{b} = G;  BPD{b} = PDb;  BES{b} = ES;
    fprintf('  %-26s aperture %5.2f m  Rayleigh %6.0f m  smallest gap %4.2f lambda\n', bs_name{b}, ...
            G.aperture, 2 * G.aperture^2 / lambda, G.min_gap / lambda);
end

% ---- main study: one parallel task per (SNR, scheme, drop) ----
r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);
CF = struct('lambda', lambda, 'eps_NF', eps_NF, 'asd_phi', ASD_varphi, 'asd_th', ASD_theta, ...
            'clu_div', clu_div, 'sL', squareLen, 'hD', hDiff, 'nAP', cape_nAP, 'nCand', cape_nCand, ...
            'J', cape_J, 'step0', cape_step0, 'nIter', cape_nIter, 'consist', cape_consist, ...
            'nAng', ca_nAngGrid, 'nRng', ca_nRangeGrid, 'nPeak', ca_nPeak, 'nSubMax', ca_Nsub_max, ...
            'nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, 'maxBr', maxBr, ...
            'mod', modOrder, 'tau_c', tau_c, 'tolS', tolS);
rng(20261119);
DRP = cell(nSetE, 1);
for ns = 1:nSetE
    DRP{ns} = make_drop(Ku, L, n_collinear_pairs, r_min_NF, r_max_en, hs_r, squareLen, ...
                        wrap_fix, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
end
sz    = [nSNR, nSchE, nSetE];
nTask = prod(sz);
T_nmse = zeros(1, nTask);  T_be = zeros(4, nTask);  T_sum = zeros(4, nTask);  T_loc = nan(1, nTask);
T_hit  = nan(1, nTask);    T_rm = nan(1, nTask);    T_time = zeros(1, nTask);  T_fh = zeros(1, nTask);
T_nf   = zeros(1, nTask);  T_seD = nan(4, Ku, nTask);  T_loc1 = nan(1, nTask);
nW = n_workers();
% run-time estimate: one CAPE task of every array at the reference SNR
tB = zeros(1, nB);
parfor b = 1:nB
    t0 = tic;
    run_task_u(DRP{1}, EN_geo{b}, BPD{b}, BES{b}, 1, Ku, SNR_lin(snr_ref_idx), 1, CF);
    tB(b) = toc(t0);
end
fprintf('\n  task time per array [s]: %s\n', sprintf('%.1f ', tB));
fprintf(['  estimated run time with %d worker(s): main study about %.0f min, user sweep about %.0f min ' ...
         '(fast_mode = %d)\n'], nW, 2 * sum(tB) * nSNR * nSetE / nW / 60, ...
        sum(tB) * nSetK * numel(Ksw) / nW / 60, fast_mode);
fprintf('\n=== [CODE U] main study: %d schemes x %d SNR x %d drops = %d tasks ===\n', nSchE, nSNR, nSetE, nTask);
tMain = tic;
parfor tk = 1:nTask
    [si, scE, ns] = ind2sub(sz, tk);
    b = en_base(scE);
    seed = 9e6 + 1e4 * ns + 1000 + si;           % same seed for every scheme (paired; as Code T)
    [O, nfv] = run_task_u(DRP{ns}, EN_geo{b}, BPD{b}, BES{b}, en_cape(scE), Ku, SNR_lin(si), seed, CF);
    T_nmse(tk) = O.nmse;  T_be(:, tk) = O.be;  T_sum(:, tk) = O.sum;
    if ~isempty(O.loc), T_loc(tk) = median(O.loc);  T_loc1(tk) = mean(O.loc < 1); end
    T_hit(tk) = O.hit;  T_rm(tk) = O.rmse;  T_time(tk) = O.time;  T_fh(tk) = O.nfus;  T_nf(tk) = nfv;
    sd = nan(4, Ku);
    if si == snr_ref_idx, sd(:, 1:size(O.se, 2)) = O.se; end
    T_seD(:, :, tk) = sd;
end
fprintf('  main study done in %.1f min\n', toc(tMain) / 60);
EN_nmse = reshape(permute(reshape(T_nmse, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_be   = reshape(permute(reshape(T_be, [4 sz]), [1 3 2 4]), 4, nSchE, 1, nSNR, nSetE);
EN_sum  = reshape(permute(reshape(T_sum, [4 sz]), [1 3 2 4]), 4, nSchE, 1, nSNR, nSetE);
EN_loc  = reshape(permute(reshape(T_loc, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_hit  = reshape(permute(reshape(T_hit, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_rmse = reshape(permute(reshape(T_rm, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);
EN_loc1 = reshape(permute(reshape(T_loc1, sz), [2 1 3]), nSchE, 1, nSNR, nSetE);   % users within 1 m
EN_time = reshape(mean(reshape(T_time, sz), 1), nSchE, 1, nSetE);
EN_fh   = reshape(mean(reshape(T_fh, sz), 1), nSchE, 1, nSetE);
EN_nf   = reshape(mean(reshape(T_nf, sz), 1), nSchE, 1, nSetE);
seAll   = reshape(T_seD, [4, Ku, sz]);
EN_seD  = reshape(seAll(:, :, snr_ref_idx, :, :), 4, Ku, nSchE, 1, nSetE);

% ---- [CODE U] user-count sweep at the reference SNR (schemes with CAPE) ----
Ksw_sch = find(en_cape == 1);
nKsw = numel(Ksw);  nKsch = numel(Ksw_sch);
fprintf('\n=== [CODE U] user-count sweep at %d dB: K = %s ===\n', SNR_dB(snr_ref_idx), mat2str(Ksw));
rng(20261120);
DRK = cell(nSetK, nKsw);
for ns = 1:nSetK
    for ik = 1:nKsw
        DRK{ns, ik} = make_drop(Ksw(ik), L, n_collinear_pairs, r_min_NF, r_max_en, hs_r, squareLen, ...
                                wrap_fix, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    end
end
nTk  = nKsch * nKsw * nSetK;
K_be = zeros(1, nTk);  K_sum = zeros(1, nTk);
tSw  = tic;
parfor tk = 1:nTk
    [jj, ik, ns] = ind2sub([nKsch, nKsw, nSetK], tk);
    b = en_base(Ksw_sch(jj));
    O = run_task_u(DRK{ns, ik}, EN_geo{b}, BPD{b}, BES{b}, 1, Ksw(ik), SNR_lin(snr_ref_idx), ...
                   3e7 + 1e4 * ns + 100 * ik, CF);
    K_be(tk) = O.be(4);  K_sum(tk) = O.sum(4);
end
EN_Kbe  = reshape(K_be, nKsch, nKsw, nSetK);
EN_Ksum = reshape(K_sum, nKsch, nKsw, nSetK);
fprintf('  user-count sweep done in %.1f min\n', toc(tSw) / 60);

if eqN_study
    % save BEFORE plotting, so a plotting problem can never lose the run
    save('codeU_v2_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', 'Dtgt', 'tolS', ...
         'nPilotBlk', 'sma_sig_d', 'bs_name', 'bs_geom', 'bs_kind', 'en_base', 'en_cape', 'en_names', ...
         'en_caseName', 'en_N0', 'en_K', 'fast_mode', 'nSetE', 'nSetK', 'nRealE', 'nSymE', 'modOrder', 'clu_div', ...
         'EN_geo', 'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_loc1', 'EN_hit', 'EN_rmse', 'EN_time', 'EN_fh', ...
         'EN_nf', 'EN_seD', 'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\nResults saved to codeU_v2_results.mat\n');

    %% ===================== [CODE U v2] FIGURES =====================
    %  V1-V9: smooth comparison figures (see below). U1-U3, U12: geometry.
    %  The older figures U4-U11, U13: set fig_old = true.
    %  Every base scheme is run twice: WITHOUT CAPE (each AP estimates its own
    %  links) and WITH CAPE (cross-AP position fusion). Colour/marker = array.
    %  Dashed line + open marker = the same antennas with difference co-array
    %  processing. SNR axis = received SNR per antenna at each user's strongest AP.
    uNB = numel(bs_name);  uNS = numel(en_names);  uNsnr = numel(SNR_dB);
    uSet = size(EN_be, 5);
    [~, uSi] = min(abs(SNR_dB - 10));
    uDet = 4;                                      % Cross-AP List detector
    uXl  = 'SNR [dB] (received, strongest AP)';
    uCapeT = {'without CAPE (local estimation)', 'with CAPE (cross-AP fusion)'};
    uCol = zeros(uNB, 3);  uMk = cell(1, uNB);  uLs = cell(1, uNB);  uMfc = cell(1, uNB);
    for b = 1:uNB
        switch bs_geom{b}
            case 'ula',     c = [0.35 0.35 0.35];  mk = 'd';
            case 'nested',  c = [0.00 0.10 0.55];  mk = 's';
            case 'coprime', c = [0.50 0.10 0.60];  mk = '^';
            case 'misc',    c = [0.85 0.45 0.00];  mk = 'v';
            case 'imisc',   c = [0.55 0.30 0.10];  mk = '>';
            case 'emisc',   c = [0.00 0.55 0.55];  mk = '<';
            case 'sma_u',   c = [0.00 0.50 0.00];  mk = 'o';
            case 'sma_s',   c = [0.80 0.00 0.00];  mk = 'p';
            otherwise,      c = [0.00 0.00 0.00];  mk = 'x';
        end
        uCol(b, :) = c;  uMk{b} = mk;  uLs{b} = '-';  uMfc{b} = c;
        if strcmp(bs_kind{b}, 'coarray'), uLs{b} = '--';  uMfc{b} = 'none'; end
    end
    uArr = find(~strcmp(bs_kind, 'coarray'));      % one entry per physical array
    uSc  = zeros(uNB, 2);                          % scheme index of (base, without/with CAPE)
    for b = 1:uNB
        for cp = 0:1
            uSc(b, cp + 1) = find(en_base == b & en_cape == cp, 1);
        end
    end
    mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
    uTmp = EN_hit;   uCnt = sum(isfinite(uTmp), 4);  uTmp(~isfinite(uTmp)) = 0;
    mHT = 100 * sum(uTmp, 4) ./ max(uCnt, 1);  mHT(uCnt == 0) = NaN;
    uTmp = EN_rmse;  uCnt = sum(isfinite(uTmp), 4);  uTmp(~isfinite(uTmp)) = 0;
    mRM = sum(uTmp, 4) ./ max(uCnt, 1);  mRM(uCnt == 0) = NaN;
    uTmp = EN_loc;   mLC = nan(uNS, 1, uNsnr);
    for sc = 1:uNS
        for si = 1:uNsnr
            v = uTmp(sc, 1, si, :);  v = v(isfinite(v));
            if ~isempty(v), mLC(sc, 1, si) = median(v); end
        end
    end

    if ~exist('fig_old', 'var'), fig_old = false; end   % true: also draw U4-U11, U13
    % ---- U1: physical arrays (positions at the data frame) ----
    figure('Name', 'U1 Physical arrays', 'Position', [30 30 1100 120 + 95 * numel(uArr)]);
    for a = 1:numel(uArr)
        b = uArr(a);  G = EN_geo{b};  pa = G.pos(:);  ap = max(pa) - min(pa);
        subplot(numel(uArr), 1, a); hold on; box on; grid on;
        plot(pa / lambda, zeros(size(pa)), uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', uCol(b, :), 'MarkerSize', 6);
        plot([0 0], [-1 1], 'k:');
        ylim([-1 1]);  set(gca, 'YTick', []);
        title(sprintf('%s:  aperture %.2f m,  Rayleigh %.0f m,  smallest gap %.2f\\lambda', bs_name{b}, ap, ...
              2 * ap^2 / lambda, min(diff(sort(pa))) / lambda));
    end
    xlabel('position along the array [\lambda]');

    % ---- U2: the movable array (frames and mirror-pair virtual lags) ----
    uSma = find(strcmp(bs_kind, 'sma'));
    if ~isempty(uSma)
        figure('Name', 'U2 Movable array: frames and virtual lags', 'Position', [30 30 520 * numel(uSma) + 80 760]);
        for a = 1:numel(uSma)
            b = uSma(a);  G = EN_geo{b};  nh = size(G.frames, 1) / 2;  F = size(G.frames, 2);
            subplot(2, numel(uSma), a); hold on; box on; grid on;
            for sg = 1:F
                plot(G.frames(1:nh, sg) / lambda, (sg - 1) * ones(nh, 1), 'b.', 'MarkerSize', 9);
                plot(G.frames(nh + 1:end, sg) / lambda, (sg - 1) * ones(nh, 1), 'r.', 'MarkerSize', 9);
            end
            plot(G.pos / lambda, G.sig_d * ones(size(G.pos)), 'ko', 'MarkerSize', 7);
            plot([0 0], [-1 F], 'k--');
            xlabel('position [\lambda]  (blue +u, red -u, circles = data frame)');  ylabel('frame \sigma');
            title(sprintf('%s: %d frames, each half moves %.2f cm', bs_name{b}, F, 100 * (F - 1) * lambda / 4));
            subplot(2, numel(uSma), numel(uSma) + a); hold on; box on; grid on;
            bar(1:numel(G.mir.cnt), G.mir.cnt, 1, 'FaceColor', uCol(b, :));
            xlabel('mirror-pair lag 2u  [\lambda/2]');  ylabel('pairs per lag');
            title(sprintf('virtual array: lags 1..%d, %d holes', numel(G.mir.cnt), sum(G.mir.cnt == 0)));
        end
    end

    % ---- U3: far-field beampattern at the data frame ----
    figure('Name', 'U3 Beampattern at the data frame', 'Position', [30 30 1100 170 * ceil(numel(uArr) / 2) + 80]);
    sg = linspace(-1, 1, 4001);
    for a = 1:numel(uArr)
        b = uArr(a);  pa = EN_geo{b}.pos(:).';
        B = abs(exp(1j * 2 * pi * (sg.') * pa / lambda) * ones(numel(pa), 1)) / numel(pa);
        subplot(ceil(numel(uArr) / 2), 2, a);
        plot(sg, 20 * log10(max(B, 1e-4)), 'Color', uCol(b, :), 'LineWidth', 1.1);
        grid on; box on; ylim([-40 0]); xlabel('sin(\theta)'); ylabel('dB'); title(bs_name{b});
    end

    % ---- U12: geometry and hardware metrics ----
    uM = zeros(numel(uArr), 6);
    for a = 1:numel(uArr)
        b = uArr(a);  G = EN_geo{b};  pa = G.pos(:);
        Dm = round((pa - pa.') / (lambda / 2));
        U = zeros(numel(pa));  nz = (Dm ~= 0) & (abs(Dm) <= 100);
        U(nz) = 1 ./ abs(Dm(nz));  U(Dm == 0) = 1;
        if strcmp(bs_kind{b}, 'sma')
            vl = sum(G.mir.cnt > 0) * lambda / 2;               % mirror-pair virtual aperture
        else
            lg = unique(round((pa - pa.') / G.unit));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
            vl = (max(diff(bk)) - 1) * G.unit;                  % hole-free difference co-array span
        end
        ap = max(pa) - min(pa);
        uM(a, :) = [ap, 2 * ap^2 / lambda, 100 * mean(EN_nf(uSc(b, 2), 1, :)), sum(Dm(:) == 1), ...
                    norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps), vl];
    end
    figure('Name', 'U12 Geometry and hardware metrics', 'Position', [30 30 1250 720]);
    uT12 = {'aperture [m]', 'Rayleigh distance [m]', 'near-field links [%]', 'w(1): pairs \lambda/2 apart', ...
            'coupling leakage', 'virtual (hole-free) aperture [m]'};
    for m = 1:6
        subplot(2, 3, m); hold on; box on; grid on;
        for a = 1:numel(uArr)
            bar(a, uM(a, m), 0.6, 'FaceColor', uCol(uArr(a), :));
        end
        set(gca, 'XTick', 1:numel(uArr), 'XTickLabel', bs_name(uArr));
        if exist('xtickangle', 'file') || exist('xtickangle', 'builtin'), xtickangle(35); end
        title(uT12{m});
    end

    if fig_old
        % ---- U4-U6: BER, sum SE and NMSE vs SNR (without / with CAPE) ----
        uFig = {'U4 BER vs SNR', 'U5 Sum SE vs SNR', 'U6 Channel NMSE vs SNR'};
        uYl  = {'BER (Cross-AP List)', 'sum SE [bit/s/Hz]', 'channel NMSE'};
        uLoc = {'southwest', 'northwest', 'southwest'};
        for f = 1:3
            figure('Name', uFig{f}, 'Position', [30 30 1200 500]);
            for cp = 0:1
                subplot(1, 2, cp + 1); hold on; box on; grid on;
                for b = 1:uNB
                    sc = uSc(b, cp + 1);
                    if f == 1
                        y = max(squeeze(mBE(uDet, sc, 1, :)), 1e-6);
                    elseif f == 2
                        y = squeeze(mSU(uDet, sc, 1, :));
                    else
                        y = squeeze(mNM(sc, 1, :));
                    end
                    plot(SNR_dB, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                         'LineWidth', 1.4, 'DisplayName', bs_name{b});
                end
                if f ~= 2, set(gca, 'YScale', 'log'); end
                if f == 1, ylim([1e-6 1]); end
                xlabel(uXl); ylabel(uYl{f}); title(sprintf('%s, %s', en_caseName{1}, uCapeT{cp + 1}));
                legend('Location', uLoc{f});
            end
        end

        % ---- U7: angle acquisition of the near-field links (hit rate and RMSE) ----
        figure('Name', 'U7 Angle acquisition of near-field links', 'Position', [30 30 1200 760]);
        for row = 1:2
            for cp = 0:1
                subplot(2, 2, (row - 1) * 2 + cp + 1); hold on; box on; grid on;
                for b = 1:uNB
                    sc = uSc(b, cp + 1);
                    if row == 1, y = squeeze(mHT(sc, 1, :)); else, y = 1e4 * squeeze(mRM(sc, 1, :)); end
                    plot(SNR_dB, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                         'LineWidth', 1.4, 'DisplayName', bs_name{b});
                end
                xlabel(uXl);
                if row == 1
                    ylabel(sprintf('hit rate [%%]  (|\\Delta sin\\theta| < %.3f)', tolS));  ylim([0 100]);
                else
                    ylabel('RMSE of sin\theta on hits [\times 10^{-4}]');  set(gca, 'YScale', 'log');
                end
                title(uCapeT{cp + 1});  legend('Location', 'best');
            end
        end

        % ---- U8: localisation error (CAPE only: positions are estimated by the fusion) ----
        figure('Name', 'U8 Localisation error (CAPE)', 'Position', [30 30 760 500]);
        hold on; box on; grid on;
        for b = 1:uNB
            sc = uSc(b, 2);
            plot(SNR_dB, squeeze(mLC(sc, 1, :)), [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                 'LineWidth', 1.4, 'DisplayName', bs_name{b});
        end
        set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
        title(sprintf('%s, with CAPE', en_caseName{1})); legend('Location', 'northeast');

        % ---- U9: what CAPE adds, at 10 dB ----
        figure('Name', 'U9 CAPE gain at 10 dB', 'Position', [30 30 1200 480]);
        uV = {zeros(uNB, 2), zeros(uNB, 2)};
        for b = 1:uNB
            for cp = 0:1
                sc = uSc(b, cp + 1);
                uV{1}(b, cp + 1) = max(mBE(uDet, sc, 1, uSi), 1e-6);
                uV{2}(b, cp + 1) = mNM(sc, 1, uSi);
            end
        end
        uYl9 = {'BER (Cross-AP List)', 'channel NMSE'};
        for m = 1:2
            subplot(1, 2, m); hold on; box on; grid on;
            bar(1:uNB, uV{m}, 'grouped');
            set(gca, 'YScale', 'log', 'XTick', 1:uNB, 'XTickLabel', bs_name);
            if exist('xtickangle', 'file') || exist('xtickangle', 'builtin'), xtickangle(35); end
            ylabel(uYl9{m});  title(sprintf('%s at %d dB', uYl9{m}, SNR_dB(uSi)));
            legend({'without CAPE', 'with CAPE'}, 'Location', 'best');
        end

        % ---- U10: BER and sum SE vs number of users K at 10 dB (CAPE schemes) ----
        mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);
        figure('Name', 'U10 BER and sum SE vs number of users', 'Position', [30 30 1200 480]);
        for f = 1:2
            subplot(1, 2, f); hold on; box on; grid on;
            for j = 1:numel(Ksw_sch)
                b = en_base(Ksw_sch(j));
                if f == 1, y = max(mKbe(j, :), 1e-6); else, y = mKsum(j, :); end
                plot(Ksw, y, [uLs{b} uMk{b}], 'Color', uCol(b, :), 'MarkerFaceColor', uMfc{b}, ...
                     'LineWidth', 1.4, 'DisplayName', bs_name{b});
            end
            xlabel('number of users K');
            if f == 1, set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)'); else, ylabel('sum SE [bit/s/Hz]'); end
            title(sprintf('with CAPE, %d dB, N = %d per AP', SNR_dB(uSi), en_N0(1)));  legend('Location', 'northwest');
        end

        % ---- U11: does co-array processing help? (paired: same antennas, same channels) ----
        uPairs = [];  uPairNm = {};
        for g = {'nested', 'coprime'}
            bp = find(strcmp(bs_geom, g{1}) & strcmp(bs_kind, 'plain'));
            bc = find(strcmp(bs_geom, g{1}) & strcmp(bs_kind, 'coarray'));
            if ~isempty(bp) && ~isempty(bc), uPairs = [uPairs; bp bc]; uPairNm{end + 1} = bs_name{bp}; end %#ok<AGROW>
        end
        if ~isempty(uPairs)
            figure('Name', 'U11 Effect of co-array processing (paired)', 'Position', [30 30 1200 760]);
            uYl11 = {'BER change with co-array [%]  (<0 = better)', 'sum-SE change with co-array [%]  (>0 = better)'};
            for row = 1:2
                for cp = 0:1
                    subplot(2, 2, (row - 1) * 2 + cp + 1); hold on; box on; grid on;
                    plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'HandleVisibility', 'off');
                    for pr = 1:size(uPairs, 1)
                        a = uSc(uPairs(pr, 1), cp + 1);  c2 = uSc(uPairs(pr, 2), cp + 1);
                        if row == 1, X = EN_be(uDet, a, 1, :, :);  Y = EN_be(uDet, c2, 1, :, :);
                        else,        X = EN_sum(uDet, a, 1, :, :); Y = EN_sum(uDet, c2, 1, :, :); end
                        X = reshape(X, uNsnr, uSet);  Y = reshape(Y, uNsnr, uSet);  Dd = Y - X;
                        ref = max(abs(mean(X, 2)), eps);
                        hEb = errorbar(SNR_dB, 100 * mean(Dd, 2) ./ ref, 100 * std(Dd, 0, 2) / sqrt(uSet) ./ ref, ...
                                       ['-' uMk{uPairs(pr, 1)}]);
                        set(hEb, 'Color', uCol(uPairs(pr, 1), :), 'LineWidth', 1.3, 'DisplayName', uPairNm{pr});
                    end
                    xlabel(uXl); ylabel(uYl11{row}); title(uCapeT{cp + 1}); legend('Location', 'best');
                end
            end
        end

        % ---- U13: the four detectors, with CAPE, for every array ----
        uDn = {'Linear MMSE', 'Hard SIC', 'List SIC', 'Cross-AP List'};  uDs = {':', '-.', '--', '-'};
        figure('Name', 'U13 BER of the four detectors (with CAPE)', 'Position', [30 30 1100 200 * ceil(numel(uArr) / 2) + 80]);
        for a = 1:numel(uArr)
            b = uArr(a);  sc = uSc(b, 2);
            subplot(ceil(numel(uArr) / 2), 2, a); hold on; box on; grid on;
            for d = 1:4
                plot(SNR_dB, max(squeeze(mBE(d, sc, 1, :)), 1e-6), [uDs{d} uMk{b}], 'Color', uCol(b, :), ...
                     'MarkerFaceColor', uCol(b, :), 'LineWidth', 1.2, 'DisplayName', uDn{d});
            end
            set(gca, 'YScale', 'log'); ylim([1e-6 1]); xlabel(uXl); ylabel('BER'); title(bs_name{b});
            legend('Location', 'southwest');
        end

        % ---- console summary at 10 dB ----
        fprintf('\n=== [CODE U] %s at %d dB, mean +- std. error over %d drops ===\n', en_caseName{1}, SNR_dB(uSi), uSet);
        fprintf('%-38s %6s %9s %21s %15s %7s %9s %8s\n', 'scheme', 'NF[%]', 'NMSE', 'BER (Cross-AP List)', 'sum SE', 'hit[%]', 'RMSE e-4', 'loc[m]');
        for sc = 1:uNS
            be = reshape(EN_be(uDet, sc, 1, uSi, :), [], 1);  su = reshape(EN_sum(uDet, sc, 1, uSi, :), [], 1);
            fprintf('%-38s %6.0f %9.2e %10.3e +- %7.1e %7.2f +- %4.2f %7.1f %9.1f %8.3f\n', en_names{sc}, ...
                    100 * mean(EN_nf(sc, 1, :)), mNM(sc, 1, uSi), mean(be), std(be) / sqrt(uSet), ...
                    mean(su), std(su) / sqrt(uSet), mHT(sc, 1, uSi), 1e4 * mRM(sc, 1, uSi), mLC(sc, 1, uSi));
        end
    end

    % ======================================================================
    %  [CODE U v2] FIGURES V1-V9. Every curve is the mean over the drops with
    %  +- one standard error; every scheme sees the same drops, channels,
    %  pilots, noise and symbols (paired). With CAPE: solid, filled marker.
    %  Without CAPE: dotted, open marker. Co-array processing: dashed.
    % ======================================================================
    vS12  = find(SNR_dB >= 12, 1);  if isempty(vS12), vS12 = uSi; end    % summary SNR (12 dB)
    vSEm  = SNR_dB <= 16;                           % sum SE shown up to 16 dB (measured SE saturates above)
    vFl   = 1 / (max(en_K) * nRealE * nSymE * log2(modOrder) * uSet);   % one bit error over all drops
    vFin  = @(v) [v(isfinite(v)), nan(1, ~any(isfinite(v)))];
    vMS   = @(X) deal(arrayfun(@(i) mean(vFin(X(i, :))), (1:size(X, 1)).'), ...
                      arrayfun(@(i) std(vFin(X(i, :))) / sqrt(max(nnz(isfinite(X(i, :))), 1)), (1:size(X, 1)).'));
    vMed  = @(X) arrayfun(@(i) median(vFin(X(i, :))), (1:size(X, 1)).');
    vRow  = @(Q, sc) reshape(Q(sc, 1, :, :), uNsnr, uSet);
    vRowB = @(d, sc) reshape(EN_be(d, sc, 1, :, :), uNsnr, uSet);
    vRowS = @(d, sc) reshape(EN_sum(d, sc, 1, :, :), uNsnr, uSet);
    vPhys = find(~strcmp(bs_kind, 'coarray'));     % physical processing: one per array
    vSpar = vPhys(~strcmp(bs_geom(vPhys), 'ula')); % sparse arrays (incl. movable)
    vSma  = find(strcmp(bs_kind, 'sma'));
    bCop  = find(strcmp(bs_geom, 'coprime') & strcmp(bs_kind, 'plain'), 1);   % reference: best static array of Code T
    if isempty(bCop), bCop = find(strcmp(bs_geom, 'nested') & strcmp(bs_kind, 'plain'), 1); end   % (or nested)
    vRef  = lower(bs_name{bCop});
    vEB   = @(h, b, ls, mf, nm) set(h, 'LineStyle', ls, 'Marker', uMk{b}, 'Color', uCol(b, :), ...
                                    'MarkerFaceColor', mf, 'LineWidth', 1.4, 'DisplayName', nm);

    % ---------------- V1: detection of all arrays (with CAPE) ----------------
    figure('Name', 'V1 Which array: detection (Cross-AP List, with CAPE)', 'Position', [30 30 1150 820]);
    for b = vPhys
        sc = uSc(b, 2);  X = vRowB(uDet, sc);  Y = vRowS(uDet, sc);
        m = max(mean(X, 2), vFl);  s = std(X, 0, 2) / sqrt(uSet);
        subplot(2, 2, 1); hold on;  vEB(errorbar(SNR_dB, m, min(s, 0.9 * m), s), b, '-', uCol(b, :), bs_name{b});
        r = mean(X <= 1e-3, 2);  rs = sqrt(r .* (1 - r) / uSet);
        subplot(2, 2, 2); hold on;  vEB(errorbar(SNR_dB, 100 * r, 100 * rs), b, '-', uCol(b, :), bs_name{b});
        subplot(2, 2, 3); hold on;
        vEB(errorbar(SNR_dB(vSEm), mean(Y(vSEm, :), 2), std(Y(vSEm, :), 0, 2) / sqrt(uSet)), b, '-', uCol(b, :), bs_name{b});
        v = sort(max(X(vS12, :), 1e-6));
        subplot(2, 2, 4); hold on;  stairs(v, (1:uSet) / uSet, 'Color', uCol(b, :), 'LineWidth', 1.6, 'DisplayName', bs_name{b});
    end
    subplot(2, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('BER (Cross-AP List)');
    title('Mean BER \pm s.e.'); legend('Location', 'southwest');
    subplot(2, 2, 2); box on; grid on; xlabel(uXl); ylabel('drops with BER \leq 10^{-3} [%]'); ylim([0 100]);
    title('Reliability: share of drops with BER \leq 10^{-3}'); legend('Location', 'northwest');
    subplot(2, 2, 3); box on; grid on; xlabel(uXl); ylabel('sum SE [bit/s/Hz]');
    title('Sum SE \pm s.e. (up to 16 dB: above it the measured SE saturates)'); legend('Location', 'northwest');
    subplot(2, 2, 4); box on; grid on; set(gca, 'XScale', 'log'); xlim([1e-6 1]); ylim([0 1]);
    xlabel('per-drop BER (error-free drops shown at 10^{-6})'); ylabel('CDF over the drops');
    title(sprintf('Distribution of the per-drop BER at %d dB', SNR_dB(vS12))); legend('Location', 'northwest');

    % ---------------- V2: angle and localisation (with CAPE) ----------------
    figure('Name', 'V2 Which array: angle and localisation (with CAPE)', 'Position', [30 30 1500 470]);
    for b = vSpar
        sc = uSc(b, 2);
        [m, s] = vMS(vRow(EN_hit, sc));
        subplot(1, 3, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, '-', uCol(b, :), bs_name{b});
        subplot(1, 3, 2); hold on;
        plot(SNR_dB, vMed(vRow(EN_loc, sc)), '-', 'Marker', uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', uCol(b, :), ...
             'LineWidth', 1.4, 'DisplayName', bs_name{b});
        [m, s] = vMS(vRow(EN_loc1, sc));
        subplot(1, 3, 3); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, '-', uCol(b, :), bs_name{b});
    end
    subplot(1, 3, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
    title(sprintf('Angle: true lobe found (|\\Delta sin\\theta| < %.3f)', tolS)); legend('Location', 'southeast');
    subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
    title('Position error (median over the drops)'); legend('Location', 'southwest');
    subplot(1, 3, 3); box on; grid on; xlabel(uXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
    title('Users located within 1 m'); legend('Location', 'southeast');

    % ---------------- V3: the movable array resolves the angle WITHOUT cross-AP help ----------------
    figure('Name', 'V3 Movable array: angle from ONE AP (no CAPE)', 'Position', [30 30 1150 470]);
    for b = vSpar
        sc = uSc(b, 1);
        [m, s] = vMS(vRow(EN_hit, sc));
        subplot(1, 2, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ':', 'none', bs_name{b});
        R = vRow(EN_rmse, sc);  rm = arrayfun(@(i) sqrt(mean(vFin(R(i, :)) .^ 2)), (1:uNsnr).');
        subplot(1, 2, 2); hold on;
        plot(SNR_dB, rm, ':', 'Marker', uMk{b}, 'Color', uCol(b, :), 'LineWidth', 1.4, 'DisplayName', bs_name{b});
    end
    subplot(1, 2, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
    title('Each AP alone (no cross-AP fusion): true lobe found'); legend('Location', 'southeast');
    subplot(1, 2, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('angle RMSE of the hits [sin\theta]');
    title('Angle accuracy when the true lobe is found (no CAPE)'); legend('Location', 'northeast');

    % ---------------- V4: load (number of users), with CAPE ----------------
    nSK = size(EN_Kbe, 3);
    figure('Name', 'V4 Which array vs number of users (with CAPE)', 'Position', [30 30 1150 560]);
    for j = 1:numel(Ksw_sch)
        sc = Ksw_sch(j);  b = en_base(sc);
        if strcmp(bs_kind{b}, 'coarray'), continue; end
        Xb = reshape(EN_Kbe(j, :, :), numel(Ksw), nSK);  Xs = reshape(EN_Ksum(j, :, :), numel(Ksw), nSK);
        m = max(mean(Xb, 2), vFl);  s = std(Xb, 0, 2) / sqrt(nSK);
        subplot(1, 2, 1); hold on;  vEB(errorbar(Ksw, m, min(s, 0.9 * m), s), b, '-', uCol(b, :), bs_name{b});
        subplot(1, 2, 2); hold on;  vEB(errorbar(Ksw, mean(Xs, 2), std(Xs, 0, 2) / sqrt(nSK)), b, '-', uCol(b, :), bs_name{b});
    end
    subplot(1, 2, 1); box on; grid on; set(gca, 'YScale', 'log'); xlabel('number of users K'); ylabel('BER (Cross-AP List)');
    title(sprintf('BER vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southoutside', 'NumColumns', 3);
    subplot(1, 2, 2); box on; grid on; xlabel('number of users K'); ylabel('sum SE [bit/s/Hz]');
    title(sprintf('Sum SE vs load, %d dB', SNR_dB(uSi))); legend('Location', 'southoutside', 'NumColumns', 3);

    % ---------------- V5: physical and hardware metrics (benefit vs cost) ----------------
    rs0 = rng;  rng(11);
    nA = numel(vPhys);  vM = zeros(nA, 6);
    for ia = 1:nA
        b = vPhys(ia);  G = EN_geo{b};  pa = G.pos(:);
        if strcmp(bs_kind{b}, 'sma')
            pe = G.frames(:);  vM(ia, 5) = 100 * (G.F - 1) * lambda / 4;  vM(ia, 6) = G.F;
        else
            pe = pa;  vM(ia, 5) = 0;  vM(ia, 6) = 1;
        end
        vM(ia, 1) = max(pe) - min(pe);                                 % estimation aperture
        vM(ia, 2) = 100 * mean(EN_nf(uSc(b, 2), 1, :));                % near-field links
        gg = round(abs(diff(sort(unique(round(pe / (lambda / 4))))))); % grid in lambda/4 units
        gu = gg(1);  for q = gg(:).', gu = gcd(gu, q); end
        vM(ia, 3) = min(2, 2 / (gu / 2));                              % grating-lobe period of the estimation aperture
        mu = 0;
        for t = 1:300                                                  % 16 users in a 25 m hotspot, 40-150 m away
            c = (40 + 110 * rand) * exp(1j * 2 * pi * rand);
            U = c + 25 * sqrt(rand(en_K(1), 1)) .* exp(1j * 2 * pi * rand(en_K(1), 1));
            dH = abs(U);  d3 = sqrt(dH .^ 2 + 100);  s = sin(angle(U)) .* dH ./ d3;
            rn = sqrt(d3.' .^ 2 + pa .^ 2 - 2 * pa * (d3.' .* s.'));
            Am = (d3.' ./ rn) .* exp(-1j * 2 * pi * rn / lambda);  Am = Am ./ sqrt(sum(abs(Am) .^ 2, 1));
            G2 = abs(Am' * Am) .^ 2;  mu = mu + mean(G2(~eye(en_K(1))));
        end
        vM(ia, 4) = mu / 300;
    end
    rng(rs0);
    vShort = strrep(strrep(strrep(strrep(strrep(bs_name(vPhys), 'Movable (uniform half)', 'Mov-U'), 'Movable (sparse half)', 'Mov-S'), 'Dense ULA', 'ULA'), 'Coprime', 'Copr.'), 'Nested', 'Nest.');
    figure('Name', 'V5 Why the arrays differ: benefit and hardware cost', 'Position', [30 30 1500 760]);
    vY = {'estimation aperture [m]', 'near-field links [%]', 'grating-lobe period in sin\theta (2 = none)', ...
          'mean |a_i^H a_j|^2 between hotspot users', 'travel of the outer element [cm]', 'pilot frames per estimate'};
    vT = {'Aperture used for estimation', 'Links in the near field', 'Angle ambiguity of the estimation aperture', ...
          'User-channel correlation at the data positions', 'Mechanical cost: travel', 'Latency cost: frames (user must stay still)'};
    for mI = 1:6
        subplot(2, 3, mI); hold on; box on; grid on;
        for ia = 1:nA
            patch(ia + [-0.3 0.3 0.3 -0.3], [0 0 vM(ia, mI) vM(ia, mI)], uCol(vPhys(ia), :), 'EdgeColor', 'k');
        end
        set(gca, 'XTick', 1:nA, 'XTickLabel', vShort); xlim([0.4 nA + 0.6]); ylabel(vY{mI}); title(vT{mI});
    end

    % ---------------- V6: movable array vs the best static array (paired) ----------------
    figure('Name', sprintf('V6 Movable array vs %s (paired, same drops)', vRef), 'Position', [30 30 1150 950]);
    vT6 = {sprintf('BER(%s) / BER(movable)  (>1: movable better)', vRef), 'true-lobe hit rate gain [points]', ...
           'position error change [%]  (<0: movable better)', 'NMSE change [%]  (<0: movable better)'};
    for f = 1:4
        subplot(2, 2, f); hold on; box on; grid on;
        if f == 1, yr = ones(size(SNR_dB)); else, yr = zeros(size(SNR_dB)); end
        plot(SNR_dB, yr, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        for b = vSma(:).'
            for cp = [1 0]
                sA = uSc(bCop, cp + 1);  sB = uSc(b, cp + 1);
                if f == 1
                    X = vRowB(uDet, sA);  Y = vRowB(uDet, sB);  mx = mean(X, 2);  my = mean(Y, 2);
                    R = mx ./ max(my, vFl);
                    cxy = arrayfun(@(i) sum((X(i, :) - mx(i)) .* (Y(i, :) - my(i))) / (uSet - 1), (1:uNsnr).');
                    se = R .* sqrt(max(var(X, 0, 2) ./ (uSet * max(mx, vFl) .^ 2) + var(Y, 0, 2) ./ (uSet * max(my, vFl) .^ 2) ...
                         - 2 * cxy ./ (uSet * max(mx, vFl) .* max(my, vFl)), 0));
                    h = errorbar(SNR_dB, R, se);
                elseif f == 2
                    [m, s] = vMS(100 * (vRow(EN_hit, sB) - vRow(EN_hit, sA)));  h = errorbar(SNR_dB, m, s);
                else
                    if f == 3, X = vRow(EN_loc, sA);  Y = vRow(EN_loc, sB); else, X = vRow(EN_nmse, sA);  Y = vRow(EN_nmse, sB); end
                    mu = nan(uNsnr, 1);  lo = mu;  hi = mu;
                    for si = 1:uNsnr
                        ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
                        if nnz(ok) < 2, continue; end
                        rr = log(Y(si, ok) ./ X(si, ok));  mm = mean(rr);  ss = std(rr) / sqrt(nnz(ok));
                        mu(si) = 100 * (exp(mm) - 1);  lo(si) = 100 * (exp(mm - ss) - 1);  hi(si) = 100 * (exp(mm + ss) - 1);
                    end
                    h = errorbar(SNR_dB, mu, mu - lo, hi - mu);
                end
                if cp == 1, vEB(h, b, '-', uCol(b, :), [bs_name{b} ', CAPE']);
                else, vEB(h, b, ':', 'none', [bs_name{b} ', no CAPE']); end
            end
        end
        xlabel(uXl); ylabel(vT6{f}); legend('Location', 'southoutside', 'NumColumns', 2);
    end

    % ---------------- V7: co-array, paired gains (same antennas) ----------------
    vCoP = [];  vCoN = {};
    for b = find(strcmp(bs_kind, 'coarray'))
        bp = find(strcmp(bs_geom, bs_geom{b}) & strcmp(bs_kind, 'plain'), 1);
        if ~isempty(bp), vCoP = [vCoP; bp b]; vCoN{end + 1} = bs_name{bp}; end %#ok<AGROW>
    end
    if ~isempty(vCoP)                              % only when a co-array scheme is present
    figure('Name', 'V7 Co-array: paired gains vs SNR (same antennas)', 'Position', [30 30 1150 820]);
    vT7 = {'position error change [%]  (<0 = co-array better)', 'true-lobe hit rate change [points]', ...
           'NMSE change [%]  (<0 = better)', 'BER change [%] (Cross-AP List)  (no loss = around 0)'};
    for f = 1:4
        subplot(2, 2, f); hold on; box on; grid on;
        plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        for pr = 1:size(vCoP, 1)
            for cp = [1 0]
                sA = uSc(vCoP(pr, 1), cp + 1);  sB = uSc(vCoP(pr, 2), cp + 1);  b = vCoP(pr, 1);
                if f == 1 || f == 3
                    if f == 1, X = vRow(EN_loc, sA);  Y = vRow(EN_loc, sB); else, X = vRow(EN_nmse, sA);  Y = vRow(EN_nmse, sB); end
                    mu = nan(uNsnr, 1);  lo = mu;  hi = mu;
                    for si = 1:uNsnr
                        ok = isfinite(X(si, :)) & isfinite(Y(si, :)) & X(si, :) > 0 & Y(si, :) > 0;
                        if nnz(ok) < 2, continue; end
                        rr = log(Y(si, ok) ./ X(si, ok));  mm = mean(rr);  ss = std(rr) / sqrt(nnz(ok));
                        mu(si) = 100 * (exp(mm) - 1);  lo(si) = 100 * (exp(mm - ss) - 1);  hi(si) = 100 * (exp(mm + ss) - 1);
                    end
                    h = errorbar(SNR_dB, mu, mu - lo, hi - mu);
                elseif f == 2
                    [m, s] = vMS(100 * (vRow(EN_hit, sB) - vRow(EN_hit, sA)));  h = errorbar(SNR_dB, m, s);
                else
                    X = vRowB(uDet, sA);  Y = vRowB(uDet, sB);  ref = max(mean(X, 2), eps);
                    h = errorbar(SNR_dB, 100 * mean(Y - X, 2) ./ ref, 100 * std(Y - X, 0, 2) / sqrt(uSet) ./ ref);
                end
                if cp == 1, vEB(h, b, '-', uCol(b, :), [vCoN{pr} ', CAPE']);
                else, vEB(h, b, ':', 'none', [vCoN{pr} ', no CAPE']); end
            end
        end
        xlabel(uXl); ylabel(vT7{f}); legend('Location', 'best');
    end

    end

    % ---------------- V8: CAPE vs no CAPE (angle and localisation) ----------------
    figure('Name', 'V8 CAPE: angle and localisation (with vs without CAPE)', 'Position', [30 30 1500 640]);
    for b = vSpar
        for cp = [1 0]
            sc = uSc(b, cp + 1);
            if cp == 1, ls = '-';  mf = uCol(b, :);  nmc = ', CAPE'; else, ls = ':';  mf = 'none';  nmc = ', no CAPE'; end
            [m, s] = vMS(vRow(EN_hit, sc));
            subplot(1, 3, 1); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ls, mf, [bs_name{b} nmc]);
            subplot(1, 3, 2); hold on;
            plot(SNR_dB, vMed(vRow(EN_loc, sc)), ls, 'Marker', uMk{b}, 'Color', uCol(b, :), 'MarkerFaceColor', mf, ...
                 'LineWidth', 1.4, 'DisplayName', [bs_name{b} nmc]);
            [m, s] = vMS(vRow(EN_loc1, sc));
            subplot(1, 3, 3); hold on;  vEB(errorbar(SNR_dB, 100 * m, 100 * s), b, ls, mf, [bs_name{b} nmc]);
        end
    end
    subplot(1, 3, 1); box on; grid on; xlabel(uXl); ylabel('true-lobe hit rate [%]'); ylim([0 100]);
    title('Angle: true lobe found'); legend('Location', 'southoutside', 'NumColumns', 2);
    subplot(1, 3, 2); box on; grid on; set(gca, 'YScale', 'log'); xlabel(uXl); ylabel('median position error [m]');
    title('Position: CAPE vs local triangulation'); legend('Location', 'southoutside', 'NumColumns', 2);
    subplot(1, 3, 3); box on; grid on; xlabel(uXl); ylabel('users located within 1 m [%]'); ylim([0 100]);
    title('Users located within 1 m'); legend('Location', 'southoutside', 'NumColumns', 2);

    % ---------------- V9: CAPE gain per detector (List SIC vs proposed Cross-AP List) ----------------
    vG = [vSpar(:).'];
    vNc = min(3, numel(vG));                       % panels per row
    figure('Name', 'V9 CAPE gain: List SIC vs the proposed Cross-AP List', 'Position', [30 30 500 * vNc 420 * ceil(numel(vG) / vNc)]);
    vDs = {'--', '-'};  vDl = {'List SIC', 'Cross-AP List (proposed)'};
    for ig = 1:numel(vG)
        b = vG(ig);  sc = uSc(b, 2);  sn = uSc(b, 1);
        subplot(ceil(numel(vG) / vNc), vNc, ig); hold on; box on; grid on;
        plot(SNR_dB, ones(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        for id = 1:2
            X = vRowB(id + 2, sn);  Y = vRowB(id + 2, sc);  mx = mean(X, 2);  my = mean(Y, 2);  ok = mx > 0 & my > 0;
            R = nan(uNsnr, 1);  se = R;
            cxy = arrayfun(@(i) sum((X(i, :) - mx(i)) .* (Y(i, :) - my(i))) / (uSet - 1), (1:uNsnr).');
            R(ok) = mx(ok) ./ my(ok);
            se(ok) = R(ok) .* sqrt(max(var(X(ok, :), 0, 2) ./ (uSet * mx(ok) .^ 2) + var(Y(ok, :), 0, 2) ./ (uSet * my(ok) .^ 2) ...
                     - 2 * cxy(ok) ./ (uSet * mx(ok) .* my(ok)), 0));
            mf = uCol(b, :);  if id == 1, mf = 'none'; end
            vEB(errorbar(SNR_dB, R, se), b, vDs{id}, mf, vDl{id});
        end
        xlabel(uXl); ylabel('BER(no CAPE) / BER(CAPE)  (>1: CAPE helps)'); title(bs_name{b}); legend('Location', 'best');
    end

    % ---------------- console summary (V figures) ----------------
    fprintf('\n=== [CODE U v2] SUMMARY AT %d dB (mean over %d drops) ===\n', SNR_dB(vS12), uSet);
    fprintf('%-26s %10s %10s %8s %10s %9s %8s %10s\n', 'array (CAPE)', 'BER', 'BER<=1e-3', 'sum SE', 'pos. [m]', 'within1m', 'hit', 'hit noCAPE');
    for b = vPhys
        sc = uSc(b, 2);  X = vRowB(uDet, sc);  Y = vRowS(uDet, sc);
        md = vMed(vRow(EN_loc, sc));  [w1, ~] = vMS(vRow(EN_loc1, sc));  [ht, ~] = vMS(vRow(EN_hit, sc));
        [h0, ~] = vMS(vRow(EN_hit, uSc(b, 1)));
        fprintf('%-26s %10.2e %9.0f%% %8.1f %10.2f %8.0f%% %7.0f%% %9.0f%%\n', bs_name{b}, mean(X(vS12, :)), ...
                100 * mean(X(vS12, :) <= 1e-3), mean(Y(vS12, :)), md(vS12), 100 * w1(vS12), 100 * ht(vS12), 100 * h0(vS12));
    end
end

function [Hhat, C, Hpre, Cpre] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...
                                             N, L, K, nReal, genie, nfOvr, eps_NF, PD, FG)
% H     : LN x nReal x K TRUE channels. Used ONLY to synthesise the received
%         pilot signal yp, as nature would. Never passed to an estimator.
% R     : N x N x L x K statistical covariances (used for far-field pairs).
% NFm   : L x K near-field indicator A_lk (large-scale, paper eq. (3)).
% genie : true  -> LEGACY: LMMSE with prior R for every pair, then near-field
%                  pairs overwritten by hdet (the true channel). If nfOvr,
%                  their error covariance is set to eps_NF*beta*I.
%         false -> AUTHENTIC: near-field pairs by polar-domain OMP, far-field
%                  pairs by statistical LMMSE. hdet is NOT used.
% Hpre/Cpre : estimates before any legacy override (== Hhat/C if authentic).
% FG    : [CODE D] optional CAPE configuration; when FG.on (authentic mode),
%         the near-field links are estimated by csi_acquire_cape.
    if nargin >= 17 && ~isempty(FG) && FG.on && ~genie && any(NFm(:))
        assert(tau == K && isequal(pilot(:), (1:K)'), ...
               'csi_acquire: CAPE needs orthogonal pilots with user k on pilot k.');
        [Hhat, C] = csi_acquire_cape(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG);
        Hpre = Hhat;
        Cpre = C;
        return;
    end
    LN = N * L;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);
    C = zeros(N, N, L, K);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        for t = 1:tau
            ue_t = find(pilot == t)';
            if isempty(ue_t)
                continue;
            end
            yp = sqrt(p) * tau * sum(H(idx, :, ue_t), 3) + sqrt(tau) * Np(:, :, l, t);
            if genie
                Psi_t = p * tau * sum(R(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R(:, :, l, k) / Psi_t;
                    Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                    C(:, :, l, k) = R(:, :, l, k) - p * tau * RPsi * R(:, :, l, k);
                end
            else
                if numel(ue_t) > 1
                    error('csi_acquire: authentic mode needs orthogonal pilots (pilot %d shared).', t);
                end
                k = ue_t;
                if NFm(l, k)
                    % de-spread: z = h + w,  w ~ CN(0, sw2*I),  sw2 = 1/(p*tau)
                    z   = yp / (sqrt(p) * tau);
                    sw2 = 1 / (p * tau);
                    [hh, se2] = pdomp_estimate(z, PD, beta(l, k), sw2);
                    Hhat(idx, :, k) = hh;
                    C(:, :, l, k) = mean(se2) * eye(N);
                else
                    Psi_t = p * tau * R(:, :, l, k) + eye(N);
                    RPsi = R(:, :, l, k) / Psi_t;
                    Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                    C(:, :, l, k) = R(:, :, l, k) - p * tau * RPsi * R(:, :, l, k);
                end
            end
        end
    end
    Hpre = Hhat;
    Cpre = C;
    if genie
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NFm(l, k)
                    Hhat(idx, :, k) = repmat(hdet(:, l, k), [1, nReal]);
                    if nfOvr
                        C(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                    end
                end
            end
        end
    end
end

function PD = build_polar_dict(pos, lambda, os_ang, q_max, n_q_min, phase_step, nS_min)
% Polar-domain dictionary (Cui & Dai, IEEE TCOM 2022) for an ARBITRARY linear
% geometry [CODE C]. Angles uniform in s = sin(theta) over [-1, 1): os_ang*N
% points, but never fewer than the aperture's resolution, ceil(2D/lambda)
% points, nor than nS_min. Inverse range q = 1/r uniform in [0, q_max] (q = 0: far-field atom)
% with spacing chosen so the Fresnel phase across the half-aperture changes
% by phase_step rad between rings, pi (D/2)^2 dq / lambda = phase_step, and at
% least n_q_min rings (n_q_min = 13 reproduces the base code for the dense
% ULA). Atoms follow the NUSW form (amplitude r/r_n, exact spherical phase)
% and have unit norm.
    PD.lambda = lambda;
    PD.pos = pos(:);
    N = numel(PD.pos);
    D = max(PD.pos) - min(PD.pos);
    nS = max([os_ang * N, ceil(2 * D / lambda), nS_min]);
    n_q = max(n_q_min, ceil(q_max / (phase_step * lambda / (pi * (D / 2)^2))) + 1);
    s_grid = -1 + (0:nS - 1) * (2 / nS);
    q_grid = linspace(0, q_max, n_q);
    [SS, QQ] = meshgrid(s_grid, q_grid);
    PD.s = SS(:).';
    PD.q = QQ(:).';
    PD.W = polar_atoms(PD.pos, PD.s, PD.q, lambda);
    PD.ds = 2 / nS;
    PD.dq = q_grid(2) - q_grid(1);
    PD.q_max = q_max;
    PD.nS = nS;
    PD.nQ = n_q;
end

function B = polar_atoms(pos, s, q, lambda)
% Unit-norm NUSW atoms, one column per (s(j), q(j)); q = 0 gives a plane wave.
% With u = r_n / r = sqrt(1 + q^2 p^2 - 2 q p s), the path difference is
% r_n - r = (q p^2 - 2 p s) / (u + 1), which stays accurate as q -> 0.
    P = pos(:);
    S = s(:).';
    Q = q(:).';
    u  = sqrt(max(1 + (Q .^ 2) .* (P .^ 2) - 2 * Q .* (P * S), 1e-12));
    dr = ((P .^ 2) * Q - 2 * P * S) ./ (u + 1);
    B  = (1 ./ u) .* exp(-1j * 2 * pi * dr / lambda);
    B  = B ./ sqrt(sum(abs(B) .^ 2, 1));
end

function [hh, se2] = pdomp_estimate(z, PD, beta_lk, sw2)
% Polar-domain OMP with gridless refinement (P-SOMP / P-SIGW family,
% Cui & Dai 2022) applied to z = h + w, one column per realization.
% Uses only: the observation z, the known noise variance sw2, the dictionary,
% and the large-scale coefficient beta_lk. The true channel is never used.
%  * An atom is accepted only if |b^H r|^2 clears gam*sw2, the level the
%    largest of D pure-noise correlations reaches (~ln D + 0.58), plus a
%    margin. Without this test OMP latches onto noise on weak links.
%  * Path gains are Wiener-scaled with DATA-estimated path powers
%    |c_LS|^2 - sw2 (unbiased), not with any true-channel quantity.
%  * If nothing is significant (weak link), the strongest refined atom is
%    kept with a prior-based Wiener gain N*beta/(N*beta + sw2), which is
%    small but never exactly zero (an exactly-zero estimate on a serving AP
%    makes the LSFD matrix singular), and the error variance is the
%    large-scale prior beta_lk (the receiver knows nothing more).
% se2 : per-element error variance computed from the residual only.
    [N, M] = size(z);
    hh  = zeros(N, M);
    se2 = zeros(1, M);
    gam = log(size(PD.W, 2)) + 0.58 + 3;
    Cz  = PD.W' * z;
    for m = 1:M
        zm = z(:, m);
        r  = zm;
        G  = zeros(N, 0);
        for it = 1:PD.smax
            if it == 1
                cc = abs(Cz(:, m));
            else
                cc = abs(PD.W' * r);
            end
            [~, ordc] = sort(cc, 'descend');
            bBest = [];
            mBest = -inf;
            for st = 1:min(PD.nstart, numel(ordc))
                j = ordc(st);
                [b, mg] = refine_atom(r, PD.s(j), PD.q(j), PD);
                if mg > mBest
                    mBest = mg;
                    bBest = b;
                end
            end
            if mBest ^ 2 < gam * sw2
                if it == 1
                    bWeak = bBest;                         % kept for the weak-link fallback
                end
                break;                                     % no significant path left
            end
            G = [G, bBest];                                %#ok<AGROW>
            cls = G \ zm;
            r = zm - G * cls;
            if real(r' * r) <= PD.stop * N * sw2
                break;
            end
        end
        s = size(G, 2);
        if s == 0                                          % weak link: nothing significant
            hh(:, m) = bWeak * ((N * beta_lk / (N * beta_lk + sw2)) * (bWeak' * zm));
            se2(m) = beta_lk;
            continue;
        end
        cls = G \ zm;
        pw  = max(abs(cls) .^ 2 - sw2, 1e-3 * sw2);        % data-estimated path powers
        A   = G' * G + sw2 * diag(1 ./ pw);
        cl  = A \ (G' * zm);                               % Wiener (LMMSE) path gains
        hh(:, m) = G * cl;
        on  = sw2 * real(trace(G * (A \ G')));            % on-support error energy
        res = zm - G * cls;
        off = max(real(res' * res) - (N - s) * sw2, 0);   % unmodelled (off-support) energy
        se2(m) = (on + off) / N;
    end
end

function [b, mg, s, q] = refine_atom(r, s0, q0, PD)
% Gridless local search around a grid atom, run in two parameterisations,
% (s, q) and (s, kappa) with kappa = (1 - s^2) q, keeping the better fit.
% (s, kappa) decouples angle and range off end-fire; (s, q) is better
% conditioned near end-fire. s wraps around [-1, 1): with half-wavelength
% spacing sin(theta) = -1 and +1 are the same end-fire direction, and a
% clamped search started at the alias s = -1 could never reach s = +0.99.
    [b1, m1, s1, q1] = refine_param(r, s0, q0, PD, false);
    [b2, m2, s2, q2] = refine_param(r, s0, q0, PD, true);
    if m1 >= m2
        b = b1;  mg = m1;  s = s1;  q = q1;
    else
        b = b2;  mg = m2;  s = s2;  q = q2;
    end
end

function [b, mg, s0, q0] = refine_param(r, s0, q0, PD, useKappa)
% 5 x 5 window whose step halves each round, maximising |b^H r|, unit-norm b.
    ds = PD.ds;
    dq = PD.dq;
    k0 = (1 - s0 ^ 2) * q0;
    iS = ceil((1:25) / 5);  iQ = mod(0:24, 5) + 1;   % [CODE T] 5 x 5 window, meshgrid order
    for rnd = 1:PD.nref
        sc = mod(s0 + ds * (-2:2) + 1, 2) - 1;
        if useKappa
            kc = max(k0 + dq * (-2:2), 0);
            SS = sc(iS);  KK = kc(iQ);
            QQ = min(KK ./ max(1 - SS .^ 2, 1e-6), PD.q_max);
        else
            qc = min(max(q0 + dq * (-2:2), 0), PD.q_max);
            SS = sc(iS);  QQ = qc(iQ);
        end
        Bc = polar_atoms(PD.pos, SS, QQ, PD.lambda);
        [~, j] = max(abs(Bc' * r));
        s0 = SS(j);
        q0 = QQ(j);
        k0 = (1 - s0 ^ 2) * q0;
        ds = ds / 2;
        dq = dq / 2;
    end
    b = polar_atoms(PD.pos, s0, q0, PD.lambda);
    mg = abs(b' * r);
end

function [SE, BER] = det_local_sic(Hhat_mc, H_mc, C, serv, p, prelog, eta, N, L, K)
    [~, ord] = sort(sum(abs(Hhat_mc).^2, 1), 'descend');
    SE = zeros(K, 1);
    BER = zeros(K, 1);
    for i = 1:K
        k_i = ord(i);
        rem = ord(i:end);
        nr = numel(rem);
        Sl = find(serv(:, k_i));
        Lk = numel(Sl);
        if Lk == 0
            Sl = (1:L)';
            Lk = L;
        end
        Ghat = zeros(Lk, nr);
        Gtrue = zeros(Lk, nr);
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            Hl_rem = Hhat_mc(idx, rem);
            Cl = sum(C(:, :, l, rem), 4);
            Ph = p * (Hl_rem * Hl_rem' + Cl) + eye(N);
            v = p * (Ph \ Hhat_mc(idx, k_i));
            Ghat(jj, :) = v' * Hhat_mc(idx, rem);
            Gtrue(jj, :) = v' * H_mc(idx, rem);
            nvar(jj)   = real(v' * v);
        end
        Rhat = p * (Ghat * Ghat') + diag(nvar);
        a = Rhat \ (sqrt(p) * Ghat(:, 1));
        des = sqrt(p) * (a' * Gtrue(:, 1));
        sig = abs(des)^2;
        intf = 0;
        if nr > 1
            intf = p * sum(abs(a' * Gtrue(:, 2:end)).^2);
        end
        noise = real(a' * (nvar .* a));
        sinr = real(sig / ((intf + noise) * (1 + eta * Lk)));
        SE(k_i) = prelog * log2(1 + sinr);
        BER(k_i) = 0.5 * erfc(sqrt(max(sinr, 0)));
    end
end

function [SE, BER] = det_local(Hhat_mc, H_mc, C, serv, p, prelog, eta, N, L, K)
    Vloc = zeros(N, K, L);
    nv = zeros(K, L);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        Hl = Hhat_mc(idx, :);
        Ph = p * (Hl * Hl' + sum(C(:, :, l, :), 4)) + eye(N);
        Vl = p * (Ph \ Hl);
        Vloc(:, :, l) = Vl;
        nv(:, l) = real(sum(conj(Vl) .* Vl, 1)).';
    end
    SE = zeros(K, 1);
    BER = zeros(K, 1);
    for k = 1:K
        Sl = find(serv(:, k));
        Lk = numel(Sl);
        if Lk == 0
            Sl = (1:L)';
            Lk = L;
        end
        Ghat = zeros(Lk, K);
        Gtrue = zeros(Lk, K);
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            v = Vloc(:, k, l);
            Ghat(jj, :) = v' * Hhat_mc(idx, :);
            Gtrue(jj, :) = v' * H_mc(idx, :);
            nvar(jj)   = nv(k, l);
        end
        Rhat = p * (Ghat * Ghat') + diag(nvar);
        a = Rhat \ (sqrt(p) * Ghat(:, k));
        des = sqrt(p) * (a' * Gtrue(:, k));
        sig = abs(des)^2;
        oth = [1:k - 1, k + 1:K];
        intf = p * sum(abs(a' * Gtrue(:, oth)).^2);
        noise = real(a' * (nvar .* a));
        sinr = real(sig / ((intf + noise) * (1 + eta * Lk)));
        SE(k) = prelog * log2(1 + sinr);
        BER(k) = 0.5 * erfc(sqrt(max(sinr, 0)));
    end
end

%% ====================================================================
%  SYMBOL-LEVEL DETECTORS  [UNCHANGED, Cross-AP V6]
%    variant=1  ->  conventional List-SIC (cluster only)
%    variant=2  ->  List+CrossAP (proposed detector)
%  Under full cooperation Sl = 1:L, so Cross-AP reduces exactly to
%  List-SIC. That reduction is the correct safety property.
%% ====================================================================
function [beL, beH, beS, beX, bits, etaList, etaXap] = ber_case(Hhat, Htrue, C, serv, p, N, L, K, nSym, M, dth, maxBr, modOrder)
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
    % [NEW] measured gating counters (ber_case only). nStg = SIC stages;
    % cnt_shadList = stages where the cluster estimate is unreliable (List
    % trigger); cnt_shadXap = stages where cross-AP fusion is invoked.
    cnt_stg      = 0;
    cnt_shadList = 0;
    cnt_shadXap  = 0;

    Vloc = cell(L, 1);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        Hl = Hhat(idx, :);
        Ph = p * (Hl * Hl' + sum(C(:, :, l, :), 4)) + eye(N);
        Vloc{l} = p * (Ph \ Hl);
    end
    aLin = cell(K, 1);
    servL = cell(K, 1);
    for k = 1:K
        Sl = find(serv(:, k));
        if isempty(Sl)
            Sl = (1:L)';
        end
        servL{k} = Sl;
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
        aLin{k} = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, k));
    end

    [~, ord] = sort(sum(abs(Hhat).^2, 1), 'descend');
    vS = cell(K, 1);
    aS = cell(K, 1);
    sS = cell(K, 1);
    for i = 1:K
        k_i = ord(i);
        rem = ord(i:end);
        Sl = find(serv(:, k_i));
        if isempty(Sl)
            Sl = (1:L)';
        end
        sS{i} = Sl;
        decoded = ord(1:i - 1);
        Lk = numel(Sl);
        Vk = zeros(N, Lk);
        Ghat = zeros(Lk, numel(rem));
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            canc_l = decoded(serv(l, decoded));
            active_l = setdiff((1:K)', canc_l(:));
            Hl = Hhat(idx, active_l);
            Cl = sum(C(:, :, l, :), 4);
            Ph = p * (Hl * Hl' + Cl) + eye(N);
            v = p * (Ph \ Hhat(idx, k_i));
            Vk(:, jj) = v;
            Ghat(jj, :) = v' * Hhat(idx, rem);
            nvar(jj) = real(v' * v);
        end
        vS{i} = Vk;
        aS{i} = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, 1));
    end

    vA = cell(K, 1);
    aA = cell(K, 1);
    for i = 1:K
        k_i = ord(i);
        rem = ord(i:end);
        Va = zeros(N, L);
        Ga = zeros(L, numel(rem));
        na = zeros(L, 1);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            Hl = Hhat(idx, rem);
            Cl = sum(C(:, :, l, :), 4);
            Ph = p * (Hl * Hl' + Cl) + eye(N);
            v = p * (Ph \ Hhat(idx, k_i));
            Va(:, l) = v;
            Ga(l, :) = v' * Hhat(idx, rem);
            na(l) = real(v' * v);
        end
        vA{i} = Va;
        aA{i} = (p * (Ga * Ga') + diag(na)) \ (sqrt(p) * Ga(:, 1));
    end

    beL = 0;
    beH = 0;
    beS = 0;
    beX = 0;
    bits = 0;
    for t = 1:nSym
        ti = randi(nC, K, 1);
        s = Cq(ti);
        n = (randn(LN, 1) + 1j * randn(LN, 1)) / sqrt(2);
        y = sqrt(p) * (Htrue * s) + n;
        bits = bits + kbit * K;

        for k = 1:K
            Sl = servL{k};
            Lk = numel(Sl);
            o = zeros(Lk, 1);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                o(jj) = Vloc{l}(:, k)' * y(idx);
            end
            sh = aLin{k}' * o;
            [~, mi] = min(abs(sh - Cq));
            beL = beL + sum(Bmap(mi, :) ~= Bmap(ti(k), :));
        end

        yr = y;
        di = zeros(K, 1);
        for i = 1:K
            k_i = ord(i);
            Sl = sS{i};
            Lk = numel(Sl);
            o = zeros(Lk, 1);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                o(jj) = vS{i}(:, jj)' * yr(idx);
            end
            sh = aS{i}' * o;
            [~, mi] = min(abs(sh - Cq));
            di(k_i) = mi;
            xh = Cq(mi);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                yr(idx) = yr(idx) - sqrt(p) * Hhat(idx, k_i) * xh;
            end
        end
        for k = 1:K
            beH = beH + sum(Bmap(di(k), :) ~= Bmap(ti(k), :));
        end

        for variant = 1:2
            useX = (variant == 2);
            P = struct('yr', {y}, 'dec', {zeros(K, 1)}, 'm', {0}, 'gr', {true});
            for i = 1:K
                k_i = ord(i);
                Sl = sS{i};
                Lk = numel(Sl);
                maxChild = numel(P) * (2 * M + 2);
                Q = repmat(struct('yr', [], 'dec', [], 'm', 0, 'gr', false), 1, maxChild);
                qn = 0;
                for pp = 1:numel(P)
                    yrp = P(pp).yr;
                    o = zeros(Lk, 1);
                    for jj = 1:Lk
                        l = Sl(jj);
                        idx = (l - 1) * N + 1:l * N;
                        o(jj) = vS{i}(:, jj)' * yrp(idx);
                    end
                    sh = aS{i}' * o;
                    [dsort_cl, ix_cl] = sort(abs(sh - Cq));
                    shad = (dsort_cl(1) >= dth_abs);
                    if useX
                        cnt_stg = cnt_stg + 1;
                        if shad, cnt_shadList = cnt_shadList + 1; end
                    end

                    if useX
                        oa = zeros(L, 1);
                        for l = 1:L
                            idx = (l - 1) * N + 1:l * N;
                            oa(l) = vA{i}(:, l)' * yrp(idx);
                        end
                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));
                        shad_xap = (dsort_xap(1) >= dth_abs);
                        if shad || shad_xap
                            cnt_shadXap = cnt_shadXap + 1;
                            nc_xap = min(M, nC);
                        else
                            nc_xap = 1;
                        end
                        if shad
                            nc_cl = min(M, nC);
                        else
                            nc_cl = 1;
                        end
                        ix_pool  = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool   = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];
                        [ix, ui]  = unique(ix_pool, 'stable');
                        dsort    = d_pool(ui);
                        nc       = numel(ix);
                    else
                        ix = ix_cl;
                        dsort = dsort_cl;
                        if shad
                            nc = min(M, nC);
                        else
                            nc = 1;
                        end
                    end

                    for c = 1:nc
                        cq = Cq(ix(c));
                        nd = P(pp).dec;
                        nd(k_i) = cq;
                        nyr = yrp;
                        if useX
                            nyr = nyr - sqrt(p) * Hhat(:, k_i) * cq;
                        else
                            for jj = 1:Lk
                                l = Sl(jj);
                                idx = (l - 1) * N + 1:l * N;
                                nyr(idx) = nyr(idx) - sqrt(p) * Hhat(idx, k_i) * cq;
                            end
                        end
                        qn = qn + 1;
                        Q(qn).yr = nyr;
                        Q(qn).dec = nd;
                        Q(qn).m = P(pp).m + dsort(c)^2;
                        Q(qn).gr = P(pp).gr && (c == 1);
                    end
                end
                Q = Q(1:qn);
                if numel(Q) > maxBr
                    mv = [Q.m];
                    gr = [Q.gr];
                    [~, ordm] = sort(mv, 'ascend');
                    keep = ordm(1:maxBr);
                    gi = find(gr, 1);
                    if ~isempty(gi) && ~any(keep == gi)
                        keep(end) = gi;
                    end
                    Q = Q(keep);
                end
                P = Q;
            end
            best = inf;
            bd = P(1).dec;
            for pp = 1:numel(P)
                r = y - sqrt(p) * (Hhat * P(pp).dec);
                m = real(r' * r);
                if m < best
                    best = m;
                    bd = P(pp).dec;
                end
            end
            be = 0;
            for k = 1:K
                [~, mk] = min(abs(bd(k) - Cq));
                be = be + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
            end
            if variant == 1
                beS = beS + be;
            else
                beX = beX + be;
            end
        end
    end
    % [NEW] measured gating fractions for this call
    etaList = cnt_shadList / max(cnt_stg, 1);   % List-SIC list-trigger rate
    etaXap  = cnt_shadXap  / max(cnt_stg, 1);   % Cross-AP invocation rate
end

function [berU, nmseU, bits_pu, eng_pu, etaXap] = metrics_case(Hhat, Htrue, C, serv, p, N, L, K, nSym, M, dth, maxBr, modOrder)
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
    chunk = 64;   % symbol block of the list detector (memory only; results do not change)
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

function [beL, beH, bits] = ber_linsic(Hhat, Htrue, C, serv, p, N, L, K, nSym, modOrder)
    % Slim empirical BER: Linear-MMSE and Hard-SIC only.  [UNCHANGED]
    [Cq, Bmap, kbit] = qam_const(modOrder);
    nC = numel(Cq);
    LN = N * L;

    Vloc = cell(L, 1);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        Hl = Hhat(idx, :);
        Ph = p * (Hl * Hl' + sum(C(:, :, l, :), 4)) + eye(N);
        Vloc{l} = p * (Ph \ Hl);
    end
    aLin = cell(K, 1);
    servL = cell(K, 1);
    for k = 1:K
        Sl = find(serv(:, k));
        if isempty(Sl)
            Sl = (1:L)';
        end
        servL{k} = Sl;
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
        aLin{k} = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, k));
    end

    [~, ord] = sort(sum(abs(Hhat).^2, 1), 'descend');
    vS = cell(K, 1);
    aS = cell(K, 1);
    sS = cell(K, 1);
    for i = 1:K
        k_i = ord(i);
        rem = ord(i:end);
        Sl = find(serv(:, k_i));
        if isempty(Sl)
            Sl = (1:L)';
        end
        sS{i} = Sl;
        decoded = ord(1:i - 1);
        Lk = numel(Sl);
        Vk = zeros(N, Lk);
        Ghat = zeros(Lk, numel(rem));
        nvar = zeros(Lk, 1);
        for jj = 1:Lk
            l = Sl(jj);
            idx = (l - 1) * N + 1:l * N;
            canc_l = decoded(serv(l, decoded));
            active_l = setdiff((1:K)', canc_l(:));
            Hl = Hhat(idx, active_l);
            Cl = sum(C(:, :, l, :), 4);
            Ph = p * (Hl * Hl' + Cl) + eye(N);
            v = p * (Ph \ Hhat(idx, k_i));
            Vk(:, jj) = v;
            Ghat(jj, :) = v' * Hhat(idx, rem);
            nvar(jj) = real(v' * v);
        end
        vS{i} = Vk;
        aS{i} = (p * (Ghat * Ghat') + diag(nvar)) \ (sqrt(p) * Ghat(:, 1));
    end

    beL = 0;
    beH = 0;
    bits = 0;
    for t = 1:nSym
        ti = randi(nC, K, 1);
        s = Cq(ti);
        n = (randn(LN, 1) + 1j * randn(LN, 1)) / sqrt(2);
        y = sqrt(p) * (Htrue * s) + n;
        bits = bits + kbit * K;

        for k = 1:K
            Sl = servL{k};
            Lk = numel(Sl);
            o = zeros(Lk, 1);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                o(jj) = Vloc{l}(:, k)' * y(idx);
            end
            sh = aLin{k}' * o;
            [~, mi] = min(abs(sh - Cq));
            beL = beL + sum(Bmap(mi, :) ~= Bmap(ti(k), :));
        end

        yr = y;
        di = zeros(K, 1);
        for i = 1:K
            k_i = ord(i);
            Sl = sS{i};
            Lk = numel(Sl);
            o = zeros(Lk, 1);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                o(jj) = vS{i}(:, jj)' * yr(idx);
            end
            sh = aS{i}' * o;
            [~, mi] = min(abs(sh - Cq));
            di(k_i) = mi;
            xh = Cq(mi);
            for jj = 1:Lk
                l = Sl(jj);
                idx = (l - 1) * N + 1:l * N;
                yr(idx) = yr(idx) - sqrt(p) * Hhat(idx, k_i) * xh;
            end
        end
        for k = 1:K
            beH = beH + sum(Bmap(di(k), :) ~= Bmap(ti(k), :));
        end
    end
end

function [Cq, B, kbit] = qam_const(M)
    Q = round(sqrt(M));
    kbit = round(log2(M));
    kb = round(log2(Q));
    amp = zeros(Q, 1);
    lab = zeros(Q, kb);
    for i = 0:Q - 1
        amp(i + 1) = 2 * i - (Q - 1);
        g = bitxor(i, floor(i / 2));
        for b = 1:kb
            lab(i + 1, b) = bitget(g, kb - b + 1);
        end
    end
    Cq = zeros(M, 1);
    B = zeros(M, kbit);
    m = 1;
    for iI = 1:Q
        for iQ = 1:Q
            Cq(m) = amp(iI) + 1j * amp(iQ);
            B(m, :) = [lab(iI, :) lab(iQ, :)];
            m = m + 1;
        end
    end
    Cq = Cq / sqrt(mean(abs(Cq).^2));
end

%% ====================================================================
%  [CODE C] ARRAY GEOMETRY, SPARSE-ARRAY CHANNEL AND AMBIGUITY HELPERS
%% ====================================================================
function G = array_geometry(type, N_ref, lambda, eta, M, Nc, Lc, Lc_sa, nN1, nN2)
% Element positions (metres, centred on the array midpoint) of the AP array.
%  'ula'    dense ULA, N_ref elements, spacing lambda/2
%  'usa'    uniform sparse array, N_ref elements, spacing eta*lambda/2
%  'eca'    extended coprime array (Zhou et al.) with Lc basic coprime arrays
%  'eca_sa' the same ECA with Lc_sa basic arrays (about the dense aperture)
%  'nested' [CODE E] two-level nested array: inner nN1 elements at spacing d0
%           plus outer nN2 elements at spacing (nN1+1) d0 (Pal & Vaidyanathan,
%           IEEE TSP 2010). Its difference co-array is a hole-free virtual ULA.
% ECA construction: one basic coprime array is an M-element sparse ULA with
% spacing Nc*d0 and an Nc-element sparse ULA with spacing M*d0 that share
% their first element (M + Nc - 1 elements). Lc/2 copies are placed with
% period M*Nc*d0 on the upper half and mirrored onto the lower half; the
% two halves share the origin, so Q = Lc*(M + Nc - 1) - 1 elements and the
% aperture is (Lc*M - 2)*Nc*d0 (for M > Nc).
    d0 = lambda / 2;
    switch lower(type)
        case 'ula'
            n = (0:N_ref - 1)';
        case 'usa'
            n = eta * (0:N_ref - 1)';
        case 'nested'
            % [CODE E] two-level nested array (Pal & Vaidyanathan, TSP 2010):
            % inner ULA at spacing d0 plus outer ULA at spacing (N1+1) d0.
            assert(nargin >= 10 && ~isempty(nN1) && ~isempty(nN2), ...
                   'array_geometry: nested needs nest_N1 and nest_N2.');
            n = unique([(1:nN1), (1:nN2) * (nN1 + 1)])';
        case {'eca', 'eca_sa'}
            if strcmpi(type, 'eca_sa')
                Lc = Lc_sa;
            end
            assert(mod(Lc, 2) == 0, 'array_geometry: the number of basic coprime arrays must be even.');
            assert(gcd(M, Nc) == 1 && M > Nc, 'array_geometry: need coprime M > N.');
            ca = unique([(0:M - 1) * Nc, (0:Nc - 1) * M]);
            up = [];
            for c = 0:Lc / 2 - 1
                up = [up, c * M * Nc + ca];               %#ok<AGROW>
            end
            n = unique([up, -up])';
        otherwise
            error('array_geometry: unknown array type "%s".', type);
    end
    p = n * d0;
    G.pos = p - (max(p) + min(p)) / 2;
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    gaps = diff(sort(n));
    G.sparsity = G.aperture / ((G.N - 1) * d0);          % aperture / dense aperture with the same N
    G.min_gap = min(gaps) * d0;
end

function G = geom_targetN(type, N0, lambda, eta, Dtgt)
% [CODE G] build each array type with EXACTLY N0 PHYSICAL elements on the
% lambda/2 grid, so every scheme is compared at equal hardware. The virtual
% (co-)array that results is whatever each geometry yields; only the physical
% count is fixed. 'eca'/'symnested' are trimmed to N0 keeping both ends;
% [CODE T] 'ulaw' is a wide-spaced ULA (spacing eta*d0); 'coprime' is the
% extended coprime array of Pal & Vaidyanathan.
    if nargin < 4 || isempty(eta), eta = 8; end
    d0 = lambda / 2;
    u  = 1;                       % [CODE T] grid unit in half-wavelengths
    switch lower(type)
        case 'ula'
            n = (0:N0 - 1)';
        case 'nested'
            N1 = floor(N0 / 2);  N2 = N0 - N1;            % exactly N0, no overlap
            n = unique([(1:N1), (1:N2) * (N1 + 1)])';
        case 'symnested'
            h  = ceil(N0 / 2);  a1 = floor(h / 2);  a2 = h - a1;
            half = unique([(1:a1), (1:a2) * (a1 + 1)]);   % a nested half (positive side)
            n = unique([-fliplr(half), half])';           % mirror -> symmetric
            n = force_count(n, N0);
        case 'eca'
            M = 5;  Nc = 4;                                % a generous coprime base
            ca = unique([(0:M - 1) * Nc, (0:Nc - 1) * M]);
            up = [];  for cc = 0:3, up = [up, cc * M * Nc + ca]; end  %#ok<AGROW>
            n = unique([up, -up])';
            n = force_count(n, N0);
        case 'ulaw'
            % [CODE T] wide-spaced ULA: uniform spacing eta*d0 (eta = 8 ->
            % d = 4 lambda). Still a uniform linear array; the wide spacing
            % enlarges the aperture at the same N and creates grating lobes.
            n = (0:N0 - 1)';  u = eta;
        case 'coprime'
            % [CODE T] extended coprime array (Pal & Vaidyanathan 2011): Nc
            % elements at spacing Mc*d0 and 2*Mc elements at spacing Nc*d0,
            % Mc < Nc coprime, sharing the origin -> 2*Mc + Nc - 1 = N0. Among
            % the valid pairs the one with the longest hole-free co-array
            % (Mc*Nc + Mc - 1 lags each side) is used: (2, 5) for N0 = 8,
            % (8, 17) for N0 = 32.
            best = -1;  Mc = 0;  Nc = 0;
            for Mt = 2:floor(N0 / 2)
                Nt = N0 + 1 - 2 * Mt;
                if Nt > Mt && gcd(Mt, Nt) == 1 && Mt * Nt + Mt - 1 > best
                    best = Mt * Nt + Mt - 1;  Mc = Mt;  Nc = Nt;
                end
            end
            if Mc == 0, error('geom_targetN: no coprime pair for N0 = %d', N0); end
            n = unique([(0:Nc - 1) * Mc, (0:2 * Mc - 1) * Nc])';
        otherwise
            error('geom_targetN: unknown type %s', type);
    end
    % [CODE T] stretch nested/coprime onto a coarser grid so that their aperture
    % matches Dtgt (the wide-ULA aperture). All positions then lie on a grid of
    % u half-wavelengths, so for u > 1 the array also has grating lobes.
    if nargin >= 5 && ~isempty(Dtgt) && Dtgt > 0 && any(strcmpi(type, {'nested', 'coprime'}))
        u = max(1, round(Dtgt / ((max(n) - min(n)) * d0)));
    end
    p = n * u * d0;
    G.pos = p - (max(p) + min(p)) / 2;
    G.unit = u * d0;              % grid unit [m] (co-array lags are multiples of it)
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    G.sparsity = G.aperture / ((G.N - 1) * d0);
    G.min_gap  = min(diff(sort(G.pos)));
end

function [cl, wv, pslv] = sparse_metrics(pos, lambda)
% [CODE P] The design metrics reported in the sparse-array literature
% (nested / coprime / MISC / IMISC / EMISC papers):
%   cl   - coupling leakage, the standard mutual-coupling measure
%          ||U - diag(U)||_F / ||U||_F with the banded model u_k = 1/k
%   wv   - the weight function values w(1), w(2), w(3), i.e. how many antenna
%          pairs sit 1, 2 and 3 spacings apart. These dominate the coupling
%          matrix, which is why EMISC drives w(1) down to 1.
%   pslv - peak sidelobe of the beampattern (ambiguity / grating lobes)
    d0 = lambda / 2;
    Dm = round((pos(:) - pos(:).') / d0);
    wv = [sum(Dm(:) == 1), sum(Dm(:) == 2), sum(Dm(:) == 3)];
    G = 100;  Na = numel(pos);  U = zeros(Na);
    for a = 1:Na
        for b = 1:Na
            kk = abs(Dm(a, b));
            if kk == 0, U(a, b) = 1; elseif kk <= G, U(a, b) = 1 / kk; end
        end
    end
    cl = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
    sg = linspace(-1, 1, 4001);
    Bm = abs(sum(exp(1j * 2 * pi * (sg.' * pos(:).') / lambda), 2)) / Na;
    [~, i0] = max(Bm);  i = i0;
    while i + 1 <= numel(Bm) && Bm(i + 1) < Bm(i), i = i + 1; end
    j = i0;
    while j - 1 >= 1 && Bm(j - 1) < Bm(j), j = j - 1; end
    Bm(j:i) = 0;  pslv = max(Bm);
end

function n = force_count(n, N0)
% Force an integer-position set to EXACTLY N0 elements: drop interior points
% evenly if too many (keep both ends to preserve aperture), extend the grid if
% too few. Guaranteed to return exactly N0 distinct integers.
    n = unique(sort(n(:)));
    if numel(n) > N0
        idx = unique(round(linspace(1, numel(n), N0)));
        n = n(idx);
        while numel(n) < N0
            full = (min(n):max(n))';
            cand = setdiff(full, n);
            if isempty(cand), n = unique([n; max(n) + 1]); else n = unique([n; cand(1)]); end
        end
    elseif numel(n) < N0
        n = [n; (max(n) + (1:(N0 - numel(n)))')];
    end
    n = unique(sort(n(:)));
end

function R = far_corr(pos_w, beta_lk, s, phi, theta, asd_phi, asd_th)
% Far-field local-scattering covariance for an ARBITRARY linear geometry
% (the Toeplitz build of the base code only holds for uniform spacing).
% pos_w: element positions in wavelengths. Same Gaussian angular-spread
% approximation as the base code, with the phase convention of the NUSW
% near-field channel (h_n ~ exp(+j 2 pi pos_n s / lambda)):
%  [R]_mn = beta exp(j 2pi Dmn s) exp(-asd_phi^2/2 (2pi Dmn cos(phi)cos(theta))^2)
%                                 exp(-asd_th^2/2  (2pi Dmn sin(theta))^2),
%  Dmn = pos_m - pos_n.
    P = pos_w(:);
    Dm = P - P.';
    R = beta_lk * exp(1j * 2 * pi * Dm * s) ...
        .* exp(-asd_phi^2 / 2 * (2 * pi * Dm * cos(phi) * cos(theta)).^2) ...
        .* exp(-asd_th^2 / 2 * (2 * pi * Dm * sin(theta)).^2);
    R = (R + R') / 2;
end

function a = nusw_vec(pos, d3, s, lambda)
% Unit-norm spherical-wave (NUSW) array response, one column per (d3(i), s(i)).
% The common phase exp(-j 2 pi d3 / lambda) is removed (it never changes |a^H b|).
    P = pos(:);
    D = d3(:).';
    S = s(:).';
    rn = sqrt(max(D .^ 2 + P .^ 2 - 2 * P * (D .* S), 1e-12));
    a = (D ./ rn) .* exp(-1j * 2 * pi * (rn - D) / lambda);
    a = a ./ sqrt(sum(abs(a) .^ 2, 1));
end

function W = wrap_images(AP, sL, fixed)
% Wrap-around images of every AP (L x 9). fixed = true: full 3 x 3 torus.
% fixed = false reproduces the base code (3 diagonal images only).
    wr = repmat([-sL 0 sL], [3 1]);
    if fixed
        off = wr(:).' + 1j * reshape(wr.', 1, []);
    else
        off = wr(:)' + 1j * (wr(:)');
    end
    W = repmat(AP(:), [1 9]) + repmat(off, [numel(AP) 1]);
end

function [UE, clm, home] = drop_users(AP, K, L, nPairs, rmin, rmax)
% User drop of the base code: nPairs collinear (co-angular) pairs around
% APs 1..nPairs, the rest round-robin around the APs. home(k) is the AP
% each user was dropped around.
    UE = zeros(K, 1);  clm = false(K, 1);  home = zeros(K, 1);
    u = 1;
    for pr = 1:nPairs
        al = mod(pr - 1, L) + 1;
        ph = 2 * pi * rand;
        rn = rmin + (rmax / 2 - rmin) * rand;
        rf = rmax / 2 + (rmax - rmax / 2) * rand;
        UE(u) = AP(al) + rn * exp(1j * ph);
        UE(u + 1) = AP(al) + rf * exp(1j * ph);
        clm(u:u + 1) = true;
        home(u:u + 1) = al;
        u = u + 2;
    end
    ar = 1;
    while u <= K
        r = rmin + (rmax - rmin) * rand;
        UE(u) = AP(ar) + r * exp(1j * 2 * pi * rand);
        home(u) = ar;
        u = u + 1;
        ar = mod(ar, L) + 1;
    end
end

function b = large_scale_gain(d3D, d2D, fc_GHz, hBS, hUT, c0, noisePow_dBm)
% 3GPP UMi LoS/NLoS path loss with shadowing, exactly as in the base code
% (one random LoS draw and one shadowing draw), normalised by the noise.
    if d2D <= 18
        P_LoS = 1;
    else
        P_LoS = (18 / d2D) * (1 - exp(-d2D / 63)) + exp(-d2D / 63);
    end
    dBP = 4 * (hBS - 1) * (hUT - 1) * fc_GHz * 1e9 / c0;
    if d3D < dBP
        PL_LoS = 28 + 37 * log10(d3D) + 20 * log10(fc_GHz);
    else
        PL_LoS = 28 + 20 * log10(d3D) + 20 * log10(fc_GHz) - 9 * log10(dBP^2 + (hBS - hUT)^2);
    end
    if rand < P_LoS
        PL_dB = PL_LoS + 4 * randn;
    else
        PL_NLoS = 13.54 + 39.08 * log10(d3D) + 20 * log10(fc_GHz) - 0.6 * (hUT - 1.5);
        PL_dB = max(PL_LoS, PL_NLoS) + 6 * randn;
    end
    b = 10 ^ ((-PL_dB - noisePow_dBm) / 10);
end

function D = cluster_cn_rule(metric_CN, beta, NF_mask, clu_div)
% Channel-norm cluster of the base code: separate NF / FF thresholds
% max/clu_div, plus the master AP.
    [L, K] = size(beta);
    D = false(L, K);
    for k = 1:K
        nf = NF_mask(:, k);
        if any(nf)
            thr_NF = max(metric_CN(nf, k)) / clu_div;
        else
            thr_NF = inf;
        end
        thr_FF = max(beta(:, k)) / clu_div;
        for l = 1:L
            if NF_mask(l, k)
                D(l, k) = metric_CN(l, k) >= thr_NF;
            else
                D(l, k) = beta(l, k) >= thr_FF;
            end
        end
        [~, lm] = max(metric_CN(:, k));
        D(lm, k) = true;
    end
end

function D = cluster_rate_rule(R, C, p, D_CN)
% Rate cluster of the base code (symbol-level D_BSR_s): APs ranked by the
% estimate-based rate log2(1 + p (tr R - tr C) / (p tr C + 1)), cluster size
% equal to the channel-norm cluster size.
    L = size(R, 3);  K = size(R, 4);
    SR = zeros(L, K);
    for l = 1:L
        for k = 1:K
            er  = real(trace(C(:, :, l, k)));
            gmm = max(real(trace(R(:, :, l, k))) - er, 0);
            SR(l, k) = log2(1 + p * gmm / (p * er + 1));
        end
    end
    D = false(L, K);
    for k = 1:K
        [~, rk] = sort(SR(:, k), 'descend');
        D(rk(1:max(nnz(D_CN(:, k)), 1)), k) = true;
    end
end

function H = draw_channel(hD, Rs, NFm, N, L, K, nReal)
% True channel realisations: deterministic NUSW channel on NF links,
% correlated Rayleigh sqrt(0.5) Rs w on FF links.
    H = randn(L * N, nReal, K) + 1j * randn(L * N, nReal, K);
    for l = 1:L
        idx = (l - 1) * N + 1:l * N;
        for k = 1:K
            if NFm(l, k)
                H(idx, :, k) = repmat(hD(:, l, k), [1, nReal]);
            else
                H(idx, :, k) = sqrt(0.5) * Rs(:, :, l, k) * H(idx, :, k);
            end
        end
    end
end

function [aR, aG, aN] = ambiguity_areas(pos, lambda, AP, UE, sL, hD, step, th, excl, wfix)
% Ambiguity area fractions of the service area (LoS spherical model).
%  aR : per AP and user, area where the beam focused on the user has gain > th
%       inside the user's main angular lobe (|s - s_k| < lambda/D) but more
%       than excl metres away: RANGE ambiguity (beam depth).
%  aG : same, outside the main angular lobe: GRATING-LOBE / SIDELOBE ambiguity.
%  aN : network level: non-coherent all-AP gain > th more than excl away.
% Path loss weights for the network gain: d^-3.7 (no shadowing).
    xg = step / 2:step:sL;
    [X, Y] = meshgrid(xg, xg);
    Xc = X(:) + 1j * Y(:);
    nG = numel(Xc);
    L = numel(AP);  K = numel(UE);
    sml = lambda / (max(pos) - min(pos));
    APw = wrap_images(AP, sL, wfix);
    UEm = mod(real(UE), sL) + 1j * mod(imag(UE), sL);
    dxr = abs(real(Xc) - real(UEm).');  dxr = min(dxr, sL - dxr);
    dxi = abs(imag(Xc) - imag(UEm).');  dxi = min(dxi, sL - dxi);
    far = sqrt(dxr .^ 2 + dxi .^ 2) > excl;                       % nG x K
    aRl = zeros(L, K);  aGl = zeros(L, K);
    num = zeros(nG, K);  ex = zeros(nG, 1);  ek = zeros(1, K);
    for l = 1:L
        [dg, wg] = min(abs(Xc - APw(l, :)), [], 2);
        d3g = max(sqrt(hD ^ 2 + dg .^ 2), 10);
        sg  = sin(angle(Xc - APw(l, wg).')) .* cos(asin(hD ./ d3g));
        [du, wu] = min(abs(UE(:) - APw(l, :)), [], 2);
        d3u = max(sqrt(hD ^ 2 + du .^ 2), 10);
        su  = sin(angle(UE(:) - APw(l, wu).')) .* cos(asin(hD ./ d3u));
        Gm = abs(nusw_vec(pos, d3g, sg, lambda)' * nusw_vec(pos, d3u, su, lambda));   % nG x K
        inML = abs(sg - su.') < sml;
        aRl(l, :) = mean((Gm > th) & inML & far, 1);
        aGl(l, :) = mean((Gm > th) & ~inML, 1);
        bx = d3g .^ -3.7;  bk = (d3u .^ -3.7).';
        num = num + sqrt(bx .* bk) .* Gm;
        ex = ex + bx;  ek = ek + bk;
    end
    Gn = num ./ sqrt(ex .* ek);
    aR = mean(aRl(:));
    aG = mean(aGl(:));
    aN = mean(mean((Gn > th) & far, 1));
end

function [G1, Gn, nA] = beam_maps(pos, lambda, AP, UE, sL, hD, step, th)
% Maps for one user (no wrap-around): gain of AP 1 alone, non-coherent
% network gain (d^-3.7 weights), and the number of APs whose own beam has
% gain > th at each point (points with nA >= 1 away from the user are
% per-AP ambiguity points).
    xg = step / 2:step:sL;
    [X, Y] = meshgrid(xg, xg);
    Xc = X(:) + 1j * Y(:);
    L = numel(AP);
    num = zeros(numel(Xc), 1);  ex = zeros(numel(Xc), 1);  ek = 0;  nA = zeros(numel(Xc), 1);
    for l = 1:L
        d3g = max(sqrt(hD ^ 2 + abs(Xc - AP(l)) .^ 2), 10);
        sg  = sin(angle(Xc - AP(l))) .* cos(asin(hD ./ d3g));
        d3u = max(sqrt(hD ^ 2 + abs(UE - AP(l)) ^ 2), 10);
        su  = sin(angle(UE - AP(l))) * cos(asin(hD / d3u));
        Gm = abs(nusw_vec(pos, d3g, sg, lambda)' * nusw_vec(pos, d3u, su, lambda));
        if l == 1
            G1 = reshape(Gm, size(X));
        end
        nA = nA + (Gm > th);
        bx = d3g .^ -3.7;  bk = d3u ^ -3.7;
        num = num + sqrt(bx * bk) .* Gm;
        ex = ex + bx;  ek = ek + bk;
    end
    Gn = reshape(num ./ sqrt(ex * ek), size(X));
    nA = reshape(nA, size(X));
end

function g = per_ap_sinr_eval(Hh, H, p, N, L, K, nReal)
% Evaluation-only per-AP SINR: local MMSE combiners built from the
% estimates, applied to the TRUE channels (noise variance 1). Scores the
% captured information kappa_k and the information surplus Delta_k.
    g = zeros(L, K);
    for mc = 1:nReal
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            Hl = reshape(Hh(idx, mc, :), [N, K]);
            Ht = reshape(H(idx, mc, :),  [N, K]);
            V  = (p * (Hl * Hl') + eye(N)) \ Hl;
            Gc = V' * Ht;
            for k = 1:K
                v   = V(:, k);
                sig = p * abs(Gc(k, k))^2;
                itf = p * (sum(abs(Gc(k, :)).^2) - abs(Gc(k, k))^2);
                g(l, k) = g(l, k) + real(sig / (itf + real(v' * v)));
            end
        end
    end
    g = g / max(nReal, 1);
end

function [kap, del] = surplus_metrics(g, D)
% Captured-information ratio kappa_k and information surplus Delta_k
% (eq. (26) of the VTC paper) of the cluster matrix D.
    Gam  = sum(g, 1);
    GamS = sum(g .* D, 1);
    kap  = GamS ./ max(Gam, eps);
    del  = (Gam - GamS) ./ ((1 + GamS) .* (1 + Gam));
end

%% ====================================================================
%  [CODE D] CROSS-AP POSITION-FUSED CHANNEL ESTIMATION (CAPE)
%% ====================================================================
function [Hhat, C, xhat, nfus] = csi_acquire_cape(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG)
% Cross-AP position-domain estimation of the near-field (NF) links.
% A user has ONE physical position; the NF channel of every AP is a known
% function of it (NUSW model). The CPU therefore
%   1. takes the J strongest local polar-domain candidates (angle, range) of
%      the FG.nCand strongest NF APs and maps each to the two positions it is
%      consistent with (a linear array cannot tell front from back);
%   2. scores every candidate position with the non-coherent fused likelihood
%      sum_{l in F} |a_l(x)^H z_l|^2 over the FG.nAP strongest NF APs, keeps
%      the best and refines it with a 2-D pattern search;
%   3. sends the position back; every NF AP builds its atom from it and
%      refines it locally in a window of half a dictionary step. CONSISTENCY
%      GATE: the AP's own local estimate is kept if it lies on the SAME lobe
%      as the fused atom (|b_loc^H b_fus| >= FG.consist) and fits better, or
%      if it fits better by more than the max-of-noise level (ln D + 0.58)
%      sw2 (high SNR: the true lobe always wins clearly); otherwise the
%      fused atom is used (wrong-lobe picks and weak links). The path gain
%      is the prior-based Wiener gain N beta / (N beta + sw2).
% Sparse arrays have narrow beams but ambiguous lobes; the lobes of
% different APs fall at different positions, so the fused likelihood
% keeps only the true one (cf. the network ambiguity maps of Code C).
% Far-field (FF) links: statistical LMMSE, as in csi_acquire.
% H is used ONLY to synthesise the received pilots. FG.UE (true positions)
% is used ONLY to score the localisation error, never by the estimator.
% Orthogonal pilots, user k on pilot k (tau = K).
    LN = N * L;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);
    C = zeros(N, N, L, K);
    xhat = nan(K, nReal);
    nfus = zeros(K, 1);                        % APs whose pilots were fused (pilot fronthaul)
    sw2 = 1 / (p * tau);
    % [CODE F] per-AP co-array spectra from the RAW pilot snapshots (all K users,
    % orthogonal pilots). Computed ONCE per AP and reused for all users.
    %   Rl = Yl Yl'  (Hermitian)     -> DIFFERENCE co-array (Ideas 2, 3, 4)
    %   Rt = Yl Yl.' (complementary) -> SUM co-array        (Idea 3 only)
    if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)
        idea = 4;  if isfield(FG, 'co_idea'), idea = FG.co_idea; end
        FG.caAngD = cell(L, 1);   % difference-co-array angle peaks per AP
        FG.caAngS = cell(L, 1);   % sum-co-array angle peaks per AP (Idea 3)
        FG.caPm   = cell(L, 1);   % difference-co-array MUSIC pseudo-spectrum (Idea 2)
        FG.caSg   = [];
        Phi = exp(-1j * 2 * pi * (0:tau - 1)' * (0:tau - 1) / tau) / sqrt(tau);
        useSide = isfield(FG, 'caSeed');      % [CODE T] covariance noise from a separate
        if useSide, rsMain = rng; rng(FG.caSeed); end   % stream, so the rest of the run
        for l = 1:L                           % stays paired with the no-co-array scheme
            idx = (l - 1) * N + 1:l * N;
            Hl = reshape(H(idx, 1, :), [N, K]);
            Yl = sqrt(p) * Hl * Phi(1:K, :) ...
                 + sqrt(0.5) * (randn(N, tau) + 1j * randn(N, tau));
            Rl = (Yl * Yl') / tau;
            [FG.caAngD{l}, FG.caSg, FG.caPm{l}] = ...
                coarray_music_angles(Rl, FG.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
            if idea == 3 && ~isempty(FG.CAs)
                Rt = (Yl * Yl.') / tau;              % complementary moment -> sum lags
                FG.caAngS{l} = coarray_music_angles(Rt, FG.CAs, K, FG.nAng, FG.nPk, FG.nSubMax);
            end
        end
        if useSide, rng(rsMain); end          % [CODE T] back to the main stream
    end
    PDw = PD;  PDw.ds = PD.ds / 2;  PDw.dq = PD.dq / 2;   % local window around the fused atom
    gamL = log(size(PD.W, 2)) + 0.58;          % level the largest of D pure-noise correlations reaches
    for k = 1:K
        Z = zeros(N, nReal, L);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            yp = sqrt(p) * tau * H(idx, :, k) + sqrt(tau) * Np(:, :, l, k);
            Z(:, :, l) = yp / (sqrt(p) * tau);
            if ~NFm(l, k)                       % far-field link: statistical LMMSE
                Psi = p * tau * R(:, :, l, k) + eye(N);
                RPsi = R(:, :, l, k) / Psi;
                Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                C(:, :, l, k) = R(:, :, l, k) - p * tau * RPsi * R(:, :, l, k);
            end
        end
        nf = find(NFm(:, k));
        if isempty(nf)
            continue;
        end
        [~, ob] = sort(beta(nf, k), 'descend');
        candAP = nf(ob(1:min(FG.nCand, numel(nf))));
        fusAP  = nf(ob(1:min(FG.nAP, numel(nf))));
        nfus(k) = numel(fusAP);
        se2 = zeros(L, nReal);
        for m = 1:nReal
            % ---- 0. local multi-start estimate of every NF link (as the per-AP OMP) ----
            bLoc = zeros(N, L);  mLoc = zeros(L, 1);  sLoc = zeros(L, FG.J);  qLoc = zeros(L, FG.J);
            for l = nf(:)'
                zl = Z(:, m, l);
                cc = abs(PD.W' * zl);
                [~, oc] = sort(cc, 'descend');
                for j = 1:FG.J
                    [bj, mj, sLoc(l, j), qLoc(l, j)] = refine_atom(zl, PD.s(oc(j)), PD.q(oc(j)), PD);
                    if mj > mLoc(l)
                        mLoc(l) = mj;
                        bLoc(:, l) = bj;
                    end
                end
            end
            % ---- 1. candidate positions from the strongest APs ----
            cand = [];
            for l = candAP(:)'
                for j = 1:FG.J
                    cand = [cand; cand_positions(FG.AP(l), sLoc(l, j), qLoc(l, j), FG.sL, FG.hD)];  %#ok<AGROW>
                end
            end
            % ---- 1b. [CODE E] extra candidates from the DIFFERENCE CO-ARRAY ----
            % The virtual array (vectorized covariance) gives high-resolution
            % angles; the range is then found by a 1-D search at that angle on
            % the physical observation, as in the near-field co-array
            % literature. These only ADD candidates: the fused likelihood of
            % step 2 still decides which one is the true position.
            if isfield(FG, 'useCA') && FG.useCA && ~isempty(FG.CA)
                idea = 4;  if isfield(FG, 'co_idea'), idea = FG.co_idea; end
                if idea == 2
                    % ----- Idea 2: cross-AP (distributed) co-array, power domain -----
                    % APs are NOT phase-synchronous, so the per-AP MUSIC pseudo-
                    % spectra of the fused APs are combined NON-COHERENTLY (summed),
                    % giving one network-level angle set; each candidate AP then
                    % range-searches at those angles.
                    Pf = zeros(size(FG.caSg));
                    for l = fusAP(:)'
                        if ~isempty(FG.caPm{l}), Pf = Pf + FG.caPm{l}; end
                    end
                    if any(Pf)
                        sSetN = coarray_peaks(FG.caSg, Pf, FG.nPk);
                        for l = candAP(:)'
                            for jj = 1:numel(sSetN)
                                qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);
                                cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>
                            end
                        end
                    end
                elseif idea == 3
                    % ----- Idea 3: joint sum-and-difference co-array -----
                    % Difference angles AND sum angles both propose candidates.
                    for l = candAP(:)'
                        sSetN = FG.caAngD{l};
                        if numel(FG.caAngS) >= l && ~isempty(FG.caAngS{l})
                            sSetN = [sSetN, FG.caAngS{l}];  %#ok<AGROW>
                        end
                        for jj = 1:numel(sSetN)
                            qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);
                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>
                        end
                    end
                elseif idea == 4
                    % ----- Idea 4: hybrid-field reconfigurable co-array -----
                    % Difference co-array angle per AP; NEAR-field links get a
                    % range search (range-aware processing), FAR-field links use a
                    % far-field candidate (q = 0, plane wave), switching per link
                    % by the Rayleigh-distance test already stored in NFm.
                    for l = candAP(:)'
                        sSetN = FG.caAngD{l};
                        for jj = 1:numel(sSetN)
                            if NFm(l, k)
                                qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);
                            else
                                qj = 0;   % far-field link: planar wavefront
                            end
                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>
                        end
                    end
                else
                    % ----- Idea 1: plain difference co-array (e.g. symmetric nested) -----
                    % Standard Pinto-style use: per-AP difference-co-array angle,
                    % range from a 1-D search on the physical observation.
                    for l = candAP(:)'
                        sSetN = FG.caAngD{l};
                        for jj = 1:numel(sSetN)
                            qj = coarray_range_search(Z(:, m, l), sSetN(jj), PD, FG.nRng);
                            cand = [cand; cand_positions(FG.AP(l), sSetN(jj), qj, FG.sL, FG.hD)];  %#ok<AGROW>
                        end
                    end
                end
            end
            if isempty(cand)                       % defensive: start next to the strongest AP
                cand = FG.AP(candAP(1)) + 1;
            end
            % ---- 2. fused likelihood, best candidate, 2-D pattern search ----
            Zf = reshape(Z(:, m, fusAP), [N, numel(fusAP)]);
            sc = fused_score(cand, Zf, fusAP, FG, PD);
            [best, jb] = max(sc);
            x = cand(jb);
            st = FG.step0;
            nbr = exp(1j * pi / 4 * (0:7)).';
            for it = 1:FG.nIter
                xn = x + st * nbr;
                sn = fused_score(xn, Zf, fusAP, FG, PD);
                [mx, jn] = max(sn);
                if mx > best
                    best = mx;
                    x = xn(jn);
                else
                    st = st / 2;
                end
            end
            xhat(k, m) = mod(real(x), FG.sL) + 1j * mod(imag(x), FG.sL);
            % ---- 3. per-AP atom from the fused position, local refinement, Wiener gain ----
            [d3x, sx] = geo_from_pos(FG.APw, x, FG.hD);
            for l = nf(:)'
                idx = (l - 1) * N + 1:l * N;
                zl = Z(:, m, l);
                [b, mf] = refine_atom(zl, sx(l), 1 / d3x(l), PDw);        % fused lobe, local precision
                sameLobe = abs(bLoc(:, l)' * b) >= FG.consist;
                clearWin = mLoc(l) ^ 2 - mf ^ 2 > gamL * sw2;               % local fit better beyond the noise level
                if (sameLobe && mLoc(l) > mf) || clearWin                  % consistency / likelihood gate
                    b = bLoc(:, l);
                end
                c = b' * zl;
                g = N * beta(l, k) / (N * beta(l, k) + sw2);
                Hhat(idx, m, k) = b * (g * c);
                res = zl - b * c;
                off = max(real(res' * res) - (N - 1) * sw2, 0);
                se2(l, m) = (N * beta(l, k) * sw2 / (N * beta(l, k) + sw2) + off) / N;
            end
        end
        for l = nf(:)'
            C(:, :, l, k) = mean(se2(l, :)) * eye(N);
        end
    end
end

function [d3, s] = geo_from_pos(APw, x, hD)
% Per-AP 3-D distance and effective sine (array axis convention of the
% channel model) of positions x (column vector), nearest wrap image.
    L = size(APw, 1);
    x = x(:).';
    d3 = zeros(L, numel(x));
    s = zeros(L, numel(x));
    for l = 1:L
        D = abs(x - APw(l, :).');                % 9 x X
        [dH, wi] = min(D, [], 1);
        img = APw(l, wi);
        d3(l, :) = max(sqrt(hD ^ 2 + dH .^ 2), 10);
        s(l, :) = sin(angle(x - img)) .* (dH ./ d3(l, :));
    end
end

function sc = fused_score(x, Zf, fusAP, FG, PD)
% Non-coherent fused likelihood sum_l |a_l(x)^H z_l|^2 (unknown complex
% gain per AP), atoms with the NUSW amplitude profile, unit norm.
    [d3, s] = geo_from_pos(FG.APw, x, FG.hD);
    sc = zeros(numel(x), 1);
    for i = 1:numel(fusAP)
        l = fusAP(i);
        A = polar_atoms(PD.pos, s(l, :), 1 ./ d3(l, :), PD.lambda);
        sc = sc + abs(A' * Zf(:, i)) .^ 2;
    end
end

function X = cand_positions(APl, s, q, sL, hD)
% Positions consistent with an AP-local (s, q) estimate: the two mirror
% solutions phi and pi - phi; for a far-field atom (q ~ 0) the range is
% unknown and the ray is sampled. A range below the AP height (possible for
% noisy lobes, q > 1/hD) is clamped to just outside it, so the list is
% never empty.
    if q > 1 / 400
        rr = 1 / q;
    else
        rr = [20 40 60 90 130 180 250];
    end
    X = [];
    for d3 = rr
        d3 = max(d3, hD + 0.5);               % a lobe closer than the AP height: user (almost) below the AP
        dH = sqrt(d3 ^ 2 - hD ^ 2);
        sp = max(min(s * d3 / dH, 1), -1);
        p1 = asin(sp);
        for ph = [p1, pi - p1]
            x = APl + dH * exp(1j * ph);
            X = [X; mod(real(x), sL) + 1j * mod(imag(x), sL)];   %#ok<AGROW>
        end
    end
end



%% ====================================================================
%  [CODE E] DIFFERENCE CO-ARRAY AND VECTORIZATION (VIRTUAL ARRAY)
%  Pal & Vaidyanathan (IEEE TSP 2010); used for communications by
%  Pinto & de Lamare (arXiv:2312.06033).
%% ====================================================================
function CA = coarray_setup(pos, lambda, mode, unit)
% Co-array of a linear array on the lambda/2 grid.
%   mode = 'diff' (default): lags x_m - x_n  -> DIFFERENCE co-array (the standard
%          Pinto building block; comes from the Hermitian covariance R = E[y y^H]).
%   mode = 'sum'          : lags x_m + x_n  -> SUM co-array (used by Idea 3; comes
%          from the complementary moment Rt = E[y y^T]). For near-field phase and
%          deterministic pilots the sum lags carry extra (range) structure.
    if nargin < 3, mode = 'diff'; end
    if nargin < 4 || isempty(unit), unit = lambda / 2; end   % [CODE T] grid of the array
    d0 = unit;
    pos = pos(:);
    % Lags come from position SUMS/DIFFERENCES: the array is centred, so pos/d0
    % can be half-integer and rounding it would corrupt the lag map, while sums
    % and differences of centred positions are exact multiples of d0.
    if strcmpi(mode, 'sum')
        Dm = round((pos + pos.') / d0);     % SUM lags
    else
        Dm = round((pos - pos.') / d0);     % DIFFERENCE lags
    end
    CA.mode = lower(mode);
    CA.lagvec = Dm(:);                  % one lag per entry of vec(R) / vec(Rt)
    [CA.uniq, ~, CA.ic] = unique(CA.lagvec);
    CA.cont = longest_contig(CA.uniq);  % longest hole-free virtual ULA
    [~, CA.loc] = ismember(CA.cont, CA.uniq);
    CA.Nv    = numel(CA.cont);          % virtual sensors used
    CA.Nsub  = floor((CA.Nv + 1) / 2);  % spatial-smoothing subarray length
    CA.DOF   = numel(CA.uniq);          % distinct lags = total virtual sensors
    CA.Nphys = numel(pos);
    CA.d0    = d0;
    CA.lambda = lambda;
    CA.pos_v = CA.cont * d0;
end

function run = longest_contig(u)
% Longest run of consecutive integers in u (the hole-free virtual ULA).
    u = sort(u(:));
    bl = u(1);  bh = u(1);  cs = u(1);
    for i = 2:numel(u)
        if u(i) ~= u(i - 1) + 1
            if (u(i - 1) - cs) > (bh - bl), bl = cs;  bh = u(i - 1);  end
            cs = u(i);
        end
    end
    if (u(end) - cs) > (bh - bl), bl = cs;  bh = u(end);  end
    run = (bl:bh)';
end

function [s_est, sg, Pm] = coarray_music_angles(Rhat, CA, K, nAng, nPk, nsubMax)
% Vectorization -> lag averaging -> spatial smoothing -> MUSIC on the VIRTUAL
% array. Returns up to nPk high-resolution angle candidates (in sin(theta)).
%   z = vec(Rhat) behaves like one snapshot of a virtual array at the lags
%   (Pal & Vaidyanathan 2010); spatial smoothing restores the rank so that
%   MUSIC can be applied (Liu & Vaidyanathan, IEEE SPL 2015).
    s_est = [];  sg = [];  Pm = [];
    if isempty(Rhat) || CA.Nsub < 3
        return;
    end
    z = Rhat(:);
    % average the entries that share the same lag (duplicate virtual sensors)
    zs = accumarray(CA.ic, z);
    zc = accumarray(CA.ic, ones(size(z)));
    zu = zs ./ max(zc, 1);
    zv = zu(CA.loc);                       % contiguous virtual ULA, ascending lag
    Ns = min(CA.Nsub, nsubMax);             % capped subarray length (speed)
    nWin = CA.Nv - Ns + 1;                  % sliding windows over the virtual ULA
    Rss = zeros(Ns);
    for i = 1:nWin                          % spatial smoothing restores the rank
        zi = zv(i:i + Ns - 1);
        Rss = Rss + zi * zi';
    end
    Rss = Rss / nWin;
    [U, Dg] = eig((Rss + Rss') / 2);
    [~, o] = sort(real(diag(Dg)), 'descend');
    U = U(:, o);
    Kk = max(1, min(K, Ns - 2));
    En = U(:, Kk + 1:end);                  % noise subspace
    sg = linspace(-0.98, 0.98, nAng);
    Av = exp(1j * 2 * pi * (0:Ns - 1)' * CA.d0 * sg / CA.lambda);
    Pm = 1 ./ max(sum(abs(En' * Av) .^ 2, 1), eps);
    s_est = coarray_peaks(sg, Pm, nPk);
end

function sp = coarray_peaks(sg, Pm, nPk)
% The nPk strongest local maxima of a 1-D spectrum.
    Pm = Pm(:).';  sg = sg(:).';
    isPk = [false, Pm(2:end - 1) > Pm(1:end - 2) & Pm(2:end - 1) > Pm(3:end), false];
    idx = find(isPk);
    if isempty(idx)
        [~, idx] = max(Pm);
    end
    [~, oo] = sort(Pm(idx), 'descend');
    idx = idx(oo(1:min(nPk, numel(oo))));
    sp = sg(idx);
end

function q = coarray_range_search(z, s, PD, nRng)
% 1-D range search at a FIXED angle on the physical observation. This is the
% standard way the near-field co-array literature recovers the range after the
% angle has been found on the virtual array.
    qg = linspace(0, PD.q_max, nRng);
    B  = polar_atoms(PD.pos, s * ones(size(qg)), qg, PD.lambda);
    [~, j] = max(abs(B' * z));
    q = qg(j);
end

%% ====================================================================
%  [CODE D] SCENARIOS AND STUDY HELPERS
%% ====================================================================
function [UE, clm] = drop_scenario(AP, K, L, nPairs, rmin, rmax, scen, hs_r, sL)
% 'baseline': the drop of the base code (collinear pairs + round robin,
%             identical random-number order).
% 'hotspot' : all K users uniform in a disc of radius hs_r around a uniform
%             random centre (crowded users: venue, station, stadium).
% 'mixed'   : K/2 users in a hotspot, K/2 with the baseline drop.
    switch lower(scen)
        case 'baseline'
            [UE, clm] = drop_users(AP, K, L, nPairs, rmin, rmax);
        case 'hotspot'
            c = (rand + 1j * rand) * sL;
            UE = c + hs_r * sqrt(rand(K, 1)) .* exp(1j * 2 * pi * rand(K, 1));
            clm = false(K, 1);
        case 'mixed'
            Kh = floor(K / 2);
            [U1, c1] = drop_users(AP, K - Kh, L, min(nPairs, floor((K - Kh) / 2)), rmin, rmax);
            c = (rand + 1j * rand) * sL;
            U2 = c + hs_r * sqrt(rand(Kh, 1)) .* exp(1j * 2 * pi * rand(Kh, 1));
            UE = [U1; U2];
            clm = [c1; false(Kh, 1)];
        otherwise
            error('drop_scenario: unknown scenario "%s".', scen);
    end
end

function FG = cape_cfg(on, AP, APw, sL, hD, nAP, nCand, J, step0, nIter, consist)
% Configuration of the cross-AP position-fused estimator (CAPE).
    FG = struct('on', on, 'AP', AP, 'APw', APw, 'sL', sL, 'hD', hD, 'nAP', nAP, ...
                'nCand', nCand, 'J', J, 'step0', step0, 'nIter', nIter, 'consist', consist);
end

function [d3, sE, phL, thL, bet] = link_geometry(APw, UE, hD, fc_GHz, hBS, hUT, c0, noisePow_dBm)
% Distances, angles and large-scale gains of every AP-user link (nearest
% wrap image), with the path-loss model of the base code.
    L = size(APw, 1);  K = numel(UE);
    d3 = zeros(L, K);  sE = d3;  phL = d3;  thL = d3;  bet = d3;
    for l = 1:L
        for k = 1:K
            [dH, wI] = min(abs(APw(l, :) - UE(k)));
            d3(l, k)  = max(sqrt(hD^2 + dH^2), 10);
            phL(l, k) = angle(UE(k) - APw(l, wI));
            thL(l, k) = asin(hD / d3(l, k));
            sE(l, k)  = sin(phL(l, k)) * cos(thL(l, k));
            bet(l, k) = large_scale_gain(d3(l, k), max(dH, 1), fc_GHz, hBS, hUT, c0, noisePow_dBm);
        end
    end
end

function ST = array_stats(G, d3, sE, phL, thL, bet, lambda, eps_NF, asd_phi, asd_th, clu_div)
% Channel statistics of one array geometry: NUSW channel on NF links
% (d < 2D^2/lambda), local-scattering covariance on FF links, and the
% channel-norm cluster.
    [L, K] = size(d3);  Na = G.N;
    Ra = 2 * G.aperture^2 / lambda;
    ST.N = Na;
    ST.NF = d3 < Ra;
    ST.hD = zeros(Na, L, K);
    ST.R2 = zeros(Na, Na, L, K);
    ST.Rs2 = zeros(Na, Na, L, K);
    for l = 1:L
        for k = 1:K
            if ST.NF(l, k)
                rn = sqrt(d3(l, k)^2 + G.pos.^2 - 2 * d3(l, k) * G.pos * sE(l, k));
                ST.hD(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
                ST.R2(:, :, l, k) = ST.hD(:, l, k) * ST.hD(:, l, k)' + eps_NF * bet(l, k) * eye(Na);
            else
                ST.R2(:, :, l, k) = far_corr(G.pos / lambda, bet(l, k), sE(l, k), phL(l, k), thL(l, k), asd_phi, asd_th);
                ST.Rs2(:, :, l, k) = sqrtm(ST.R2(:, :, l, k));
            end
        end
    end
    mCN = bet;
    mCN(ST.NF) = Na * bet(ST.NF);
    ST.Dcn = cluster_cn_rule(mCN, bet, ST.NF, clu_div);
end

function O = eval_scheme(ST, bet, useCape, PD, FG, p, EP, wantRef, UEm, seed)
% One scheme at one SNR: draw channels (seeded, so the local and the CAPE
% variant of the same array see IDENTICAL channels, pilots and symbols),
% estimate, cluster, detect with estimated and perfect CSI.
    [L, K] = size(bet);  Na = ST.N;
    rng(seed);
    FG.caSeed = mod(seed + 5e8, 2^32 - 1);   % [CODE T] side stream for co-array covariance noise
    Ht = draw_channel(ST.hD, ST.Rs2, ST.NF, Na, L, K, EP.nReal);
    tE = tic;
    O.loc = [];
    O.nfus = 0;
    if useCape
        [Hh, Cs, xh, nfus] = csi_acquire_cape(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG);
        O.nfus = mean(nfus);
        dx = abs(real(xh) - real(UEm(:)));  dx = min(dx, FG.sL - dx);     % torus distance
        dy = abs(imag(xh) - imag(UEm(:)));  dy = min(dy, FG.sL - dy);
        e = sqrt(dx.^2 + dy.^2);
        O.loc = e(isfinite(e));
    else
        [Hh, Cs] = csi_acquire(Ht, ST.R2, ST.NF, ST.hD, bet, (1:K)', K, p, ...
                               Na, L, K, EP.nReal, false, false, EP.eps_NF, PD);
    end
    O.time = toc(tE);
    O.nmse = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);
    Drt = cluster_rate_rule(ST.R2, Cs, p, ST.Dcn);
    O.sz = mean(sum(Drt, 1));
    bE = zeros(4, K);  nE = zeros(4, K);  nb = 0;  ne = 0;  eX = 0;
    for mc = 1:EP.nReal
        Hm = reshape(Ht(:, mc, :), [L * Na, K]);
        He = reshape(Hh(:, mc, :), [L * Na, K]);
        [bu, nu, bpu, epu, ex] = metrics_case(He, Hm, Cs, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);
        bE = bE + bu;  nE = nE + nu;  nb = nb + bpu;  ne = ne + epu;  eX = eX + ex;
    end
    O.be  = sum(bE, 2) / max(K * nb, 1);
    O.bp  = nan(4, 1);   % [CODE T] perfect-CSI pass removed: never used here, it doubled the run time
    O.eta = eX / EP.nReal;
    seU   = EP.prelog * log2(1 + 1 ./ max(nE / max(ne, 1), 1 / (EP.nReal * EP.nSym)));   % 4 x K effective-SINR SE
    O.se  = seU;
    O.sum = sum(seU, 2);
    O.del = nan(1, K);
    if wantRef
        g = per_ap_sinr_eval(Hh, Ht, p, Na, L, K, EP.nReal);
        [~, O.del] = surplus_metrics(g, Drt);
    end
end

function s = ternary_str(c, a, b)
% a if c is true, b otherwise (for printing)
    if c
        s = a;
    else
        s = b;
    end
end

function D = make_drop(Kc, L, nCol, rmin, rmax, hs_r, sL, wfix, hD, fc_GHz, hBS, hUT, c0, noisePow_dBm)
% [CODE T] one random AP/user drop with uplink power control (received SNR at
% each user's strongest AP equals the SNR axis value).
    D.APp = (rand(L, 1) + 1j * rand(L, 1)) * sL;
    UEp = drop_scenario(D.APp, Kc, L, min(nCol, floor(Kc / 2)), rmin, rmax, 'hotspot', hs_r, sL);
    D.UEm = mod(real(UEp), sL) + 1j * mod(imag(UEp), sL);
    D.APw = wrap_images(D.APp, sL, wfix);
    [D.d3, D.sE, D.phL, D.thL, bet] = link_geometry(D.APw, UEp, hD, fc_GHz, hBS, hUT, c0, noisePow_dBm);
    D.bet = bet ./ max(bet, [], 1);
end

function [O, nf] = run_task(D, G, PD, CA, CAs, idea, useCape, Kc, p, seed, CF)
% [CODE T] one simulation task: one drop, one array/processing scheme, one SNR.
    ST = array_stats(G, D.d3, D.sE, D.phL, D.thL, D.bet, CF.lambda, CF.eps_NF, CF.asd_phi, CF.asd_th, CF.clu_div);
    nf = mean(ST.NF(:));
    FG = cape_cfg(true, D.APp, D.APw, CF.sL, CF.hD, CF.nAP, CF.nCand, CF.J, CF.step0, CF.nIter, CF.consist);
    FG.useCA = (idea > 0);  FG.co_idea = idea;  FG.CA = CA;  FG.CAs = CAs;
    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    EP = struct('nReal', CF.nReal, 'nSym', CF.nSym, 'M', CF.M, 'dth', CF.dth, 'maxBr', CF.maxBr, ...
                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF);
    O = eval_scheme(ST, D.bet, useCape, PD, FG, p, EP, false, D.UEm, seed);
end

function nW = n_workers()
% number of parallel workers (1 without a pool or without the Parallel Toolbox)
    nW = 1;
    try
        pl = gcp('nocreate');
        if isempty(pl), pl = gcp; end
        nW = max(1, pl.NumWorkers);
    catch
    end
end

%% ====================================================================
%  [CODE U] MOVABLE ARRAY, GEOMETRIES AND THE UNIFIED ESTIMATOR
%% ====================================================================
function G = geom_u(type, N0, lambda, Dtgt, F, sig_d, ies)
% Geometry of one array with N0 physical antennas. Sparse arrays are put on a
% grid of u half-wavelengths so that their aperture is about Dtgt.
    d0 = lambda / 2;
    switch lower(type)
        case 'ula'
            G = geom_targetN('ula', N0, lambda, 8);
        case {'nested', 'coprime'}
            G = geom_targetN(type, N0, lambda, 8, Dtgt);
        case {'misc', 'imisc', 'emisc'}
            n = [0; cumsum(ies(:))];
            assert(numel(n) == N0, 'geom_u: the %s spacing set needs N0 - 1 = %d entries.', type, N0 - 1);
            u = max(1, round(Dtgt / ((max(n) - min(n)) * d0)));
            G = finish_geom(n * u * d0, u * d0, type, d0);
        case {'sma_u', 'sma_s'}
            % MOVABLE MIRROR ARRAY: positive half at u = (b_i + sigma) lambda/4,
            % negative half at -u, frames sigma = 0 .. F-1.
            nh = N0 / 2;  q4 = lambda / 4;
            if strcmpi(type, 'sma_u')
                b = 1 + 16 * (0:nh - 1);                   % uniform half, pitch 4 lambda
            else
                b = 1 + 4 * (0:nh - 1) .* (1:nh);          % increasing gaps 2, 4, 6 ... lambda
            end
            U = (b(:) + (0:F - 1)) * q4;                   % nh x F
            G = finish_geom([U(:, sig_d + 1); -U(:, sig_d + 1)], q4, type, d0);
            G.frames = [U; -U];                            % rows 1..nh: +u_i, rows nh+1..N0: -u_i
            G.b = b;  G.F = F;  G.sig_d = sig_d;
            % mirror pairs of the stacked frames: lag 2u = (b_i + sigma) lambda/2
            [ii, ss] = ndgrid(1:nh, 1:F);
            MI.iPos = (ss(:) - 1) * N0 + ii(:);
            MI.iNeg = (ss(:) - 1) * N0 + nh + ii(:);
            MI.q    = b(ii(:)).' + ss(:) - 1;
            MI.Q    = max(MI.q);
            MI.cnt  = accumarray(MI.q, 1, [MI.Q 1]);
            MI.nfft = 2 ^ nextpow2(8 * MI.Q);
            G.mir = MI;
            if any(MI.cnt == 0)
                warning('geom_u: %s mirror lags have %d holes (increase F).', type, sum(MI.cnt == 0));
            end
        otherwise
            error('geom_u: unknown array %s', type);
    end
end

function G = finish_geom(p, unit, type, d0)
    p = p(:);
    G.pos = p - (max(p) + min(p)) / 2;
    G.unit = unit;
    G.N = numel(G.pos);
    G.aperture = max(G.pos) - min(G.pos);
    G.type = lower(type);
    G.sparsity = G.aperture / ((G.N - 1) * d0);
    G.min_gap = min(diff(sort(G.pos)));
end

function h = nf_channel(pos, d3, sE, bet, lambda)
% NUSW near-field channel (same model as array_stats) at arbitrary positions.
    [L, K] = size(d3);  pos = pos(:);  h = zeros(numel(pos), L, K);
    for l = 1:L
        for k = 1:K
            rn = sqrt(d3(l, k)^2 + pos.^2 - 2 * d3(l, k) * pos * sE(l, k));
            h(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
        end
    end
end

function sx = sma_mirror_angles(z, MI, nPk)
% Idea 1: z(+u) conj(z(-u)) ~ exp(j pi q sin(theta)) on the lag q = 2u/(lambda/2);
% the near-field range term is the same on both twins and cancels. Average the
% pairs per lag, take the spectrum of the hole-free virtual array, keep peaks.
    w = z(MI.iPos) .* conj(z(MI.iNeg));
    v = accumarray(MI.q, w, [MI.Q 1]) ./ max(MI.cnt, 1);
    sp = abs(fft([0; v], MI.nfft));                   % index k <-> sin(theta) = 2k/nfft
    pk = find(sp > circshift(sp, 1) & sp >= circshift(sp, -1));
    [~, o] = sort(sp(pk), 'descend');
    pk = pk(o(1:min(nPk, numel(o))));
    sx = mod(2 * (pk(:).' - 1) / MI.nfft + 1, 2) - 1;
end

function [Hhat, C, xhat, nfus, aErr, cSa, cQa] = csi_nf_estimate(H, R, NFm, beta, tau, p, N, L, K, nReal, PD, FG, EST, useCape)
% Unified estimator. Far-field links: statistical LMMSE (data-block pilot).
% Near-field links: single-atom polar estimation on the estimation observation
% EST.Ze (static arrays: one pilot block, as Code T; movable array: the
% frames stacked into a synthetic aperture, each frame with 1/F of the energy). Candidates per link = FG.J refined
% dictionary atoms + extra angles (co-array MUSIC or mirror-pair angles) with a
% 1-D range search. Without CAPE each AP keeps its best-fitting candidate; with
% CAPE the candidates become positions and the cross-AP fused likelihood picks
% one (as csi_acquire_cape, same consistency gate). The channel at the data
% positions is rebuilt from the chosen (angle, range).
    LN = N * L;
    PDe = EST.PDe;  Ne = EST.Ne;  sw2e = EST.sw2e;
    Np = sqrt(0.5) * (randn(N, nReal, L, tau) + 1j * randn(N, nReal, L, tau));
    Hhat = zeros(LN, nReal, K);  C = zeros(N, N, L, K);
    xhat = nan(K, nReal);  nfus = zeros(K, 1);  aErr = nan(L, K, nReal);
    PDw = PDe;  PDw.ds = PDe.ds / 2;  PDw.dq = PDe.dq / 2;
    gamL = log(size(PDe.W, 2)) + 0.58;
    nc = FG.J + EST.nExtra;
    cSa = nan(L, K, nReal, nc);  cQa = nan(L, K, nReal, nc);   % [v2] local candidate lobes (triangulation baseline)
    for k = 1:K
        for l = 1:L
            if ~NFm(l, k)                               % far-field link: statistical LMMSE
                idx = (l - 1) * N + 1:l * N;
                yp = sqrt(p) * tau * H(idx, :, k) + sqrt(tau) * Np(:, :, l, k);
                Psi = p * tau * R(:, :, l, k) + eye(N);
                RPsi = R(:, :, l, k) / Psi;
                Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                C(:, :, l, k) = R(:, :, l, k) - p * tau * RPsi * R(:, :, l, k);
            end
        end
        nf = find(NFm(:, k));
        if isempty(nf), continue; end
        [~, ob] = sort(beta(nf, k), 'descend');
        candAP = nf(ob(1:min(FG.nCand, numel(nf))));
        fusAP  = nf(ob(1:min(FG.nAP, numel(nf))));
        if useCape, nfus(k) = numel(fusAP); end
        se2 = zeros(L, nReal);
        for m = 1:nReal
            Z = reshape(EST.Ze(:, m, :, k), Ne, L);
            cS = nan(L, nc);  cQ = nan(L, nc);
            bLoc = zeros(Ne, L);  mLoc = zeros(L, 1);  sLoc = zeros(L, 1);  qLoc = zeros(L, 1);
            for l = nf(:)'
                zl = Z(:, l);
                [~, oc] = sort(abs(PDe.W' * zl), 'descend');
                for j = 1:FG.J
                    [bj, mj, cS(l, j), cQ(l, j)] = refine_atom(zl, PDe.s(oc(j)), PDe.q(oc(j)), PDe);
                    if mj > mLoc(l), mLoc(l) = mj;  bLoc(:, l) = bj;  sLoc(l) = cS(l, j);  qLoc(l) = cQ(l, j); end
                end
                if strcmp(EST.kind, 'coarray')
                    sx = EST.angCA{l};
                elseif strcmp(EST.kind, 'sma')
                    sx = sma_mirror_angles(zl, EST.mir, EST.nPkSMA);
                else
                    sx = [];
                end
                for j = 1:min(numel(sx), EST.nExtra)
                    q0 = coarray_range_search(zl, sx(j), PDe, FG.nRng);
                    [bj, mj, cS(l, FG.J + j), cQ(l, FG.J + j)] = refine_atom(zl, sx(j), q0, PDe);
                    if mj > mLoc(l)
                        mLoc(l) = mj;  bLoc(:, l) = bj;  sLoc(l) = cS(l, FG.J + j);  qLoc(l) = cQ(l, FG.J + j);
                    end
                end
            end
            cSa(:, k, m, :) = reshape(cS, L, 1, 1, nc);  cQa(:, k, m, :) = reshape(cQ, L, 1, 1, nc);
            if useCape
                cand = [];
                for l = candAP(:)'
                    for j = 1:nc
                        if isfinite(cS(l, j))
                            cand = [cand; cand_positions(FG.AP(l), cS(l, j), cQ(l, j), FG.sL, FG.hD)];  %#ok<AGROW>
                        end
                    end
                end
                if isempty(cand), cand = FG.AP(candAP(1)) + 1; end
                Zf = Z(:, fusAP);
                sc = fused_score(cand, Zf, fusAP, FG, PDe);
                [best, jb] = max(sc);
                x = cand(jb);  st = FG.step0;  nbr = exp(1j * pi / 4 * (0:7)).';
                for it = 1:FG.nIter
                    xn = x + st * nbr;
                    sn = fused_score(xn, Zf, fusAP, FG, PDe);
                    [mx, jn] = max(sn);
                    if mx > best
                        best = mx;  x = xn(jn);
                    else
                        st = st / 2;
                    end
                end
                xhat(k, m) = mod(real(x), FG.sL) + 1j * mod(imag(x), FG.sL);
                [d3x, sxp] = geo_from_pos(FG.APw, x, FG.hD);
            end
            for l = nf(:)'
                zl = Z(:, l);
                b = bLoc(:, l);  s1 = sLoc(l);  q1 = qLoc(l);
                if useCape
                    [bf, mfu, sf, qf] = refine_atom(zl, sxp(l), 1 / d3x(l), PDw);   % fused lobe
                    sameLobe = abs(bLoc(:, l)' * bf) >= FG.consist;
                    clearWin = mLoc(l) ^ 2 - mfu ^ 2 > gamL * sw2e;
                    if ~((sameLobe && mLoc(l) > mfu) || clearWin)
                        b = bf;  s1 = sf;  q1 = qf;
                    end
                    [b2, m2, s2, q2] = refine_atom(zl, s1, q1, PDe);   % [v2] full window from the chosen lobe
                    if m2 > abs(b' * zl)
                        b = b2;  s1 = s2;  q1 = q2;
                    end
                end
                c = b' * zl;
                g = Ne * beta(l, k) / (Ne * beta(l, k) + sw2e);
                idx = (l - 1) * N + 1:l * N;
                Hhat(idx, m, k) = polar_atoms(PD.pos, s1, q1, PD.lambda) * (g * c * sqrt(N / Ne));
                res = zl - b * c;
                off = max(real(res' * res) - (Ne - 1) * sw2e, 0);
                se2(l, m) = (Ne * beta(l, k) * sw2e / (Ne * beta(l, k) + sw2e) + off) / Ne;
                dS = abs(s1 - EST.sTrue(l, k));
                aErr(l, k, m) = min(dS, 2 - dS);
            end
        end
        for l = nf(:)'
            C(:, :, l, k) = mean(se2(l, :)) * eye(N);
        end
    end
end

function O = eval_scheme_u(ST, bet, useCape, PD, FG, p, EP, UEm, seed, ES)
% One scheme at one SNR: channels, pilots of all blocks / frames, estimation,
% clustering and the four detectors (estimated CSI).
    [L, K] = size(bet);  Na = ST.N;
    rng(seed);
    Ht = draw_channel(ST.hD, ST.Rs2, ST.NF, Na, L, K, EP.nReal);
    tE = tic;
    sw2 = 1 / (p * K);                                   % de-spread pilot noise of one block (Code T)
    if strcmp(ES.kind, 'sma')
        sw2e = sw2 * ES.F;                               % F frames, each with 1/F of the pilot energy
    else
        sw2e = sw2;                                      % static: one pilot block, as Code T
    end
    Ze = repmat(reshape(ES.he, ES.Ne, 1, L, K), [1, EP.nReal, 1, 1]) + ...
         sqrt(sw2e / 2) * (randn(ES.Ne, EP.nReal, L, K) + 1j * randn(ES.Ne, EP.nReal, L, K));
    EST = struct('kind', ES.kind, 'Ze', Ze, 'Ne', ES.Ne, 'sw2e', sw2e, 'PDe', ES.PDe, ...
                 'sTrue', ES.sTrue, 'nPkSMA', 2, 'nExtra', 0);
    EST.mir = ES.mir;  EST.angCA = cell(L, 1);
    if strcmp(ES.kind, 'coarray')
        rsMain = rng;  rng(mod(seed + 5e8, 2^32 - 1));     % covariance noise from a side stream (paired)
        Phi = exp(-1j * 2 * pi * (0:K - 1)' * (0:K - 1) / K) / sqrt(K);
        for l = 1:L
            idx = (l - 1) * Na + 1:l * Na;
            Hl = reshape(Ht(idx, 1, :), [Na, K]);
            Yl = sqrt(p) * Hl * Phi + sqrt(0.5) * (randn(Na, K) + 1j * randn(Na, K));   % one block, as Code T
            EST.angCA{l} = coarray_music_angles((Yl * Yl') / K, ES.CA, K, FG.nAng, FG.nPk, FG.nSubMax);
        end
        rng(rsMain);
        EST.nExtra = FG.nPk;
    elseif strcmp(ES.kind, 'sma')
        EST.nExtra = EST.nPkSMA;
    end
    [Hh, Cs, xh, nfus, aErr, cSa, cQa] = csi_nf_estimate(Ht, ST.R2, ST.NF, bet, K, p, Na, L, K, EP.nReal, PD, FG, EST, useCape);
    O.time = toc(tE);
    O.nfus = mean(nfus);
    O.loc = [];
    if ~useCape                                       % [v2] no-CAPE localisation baseline: at most 4
        nU = min(size(cSa, 4), 4);                    % lobes per AP (dictionary lobes + strongest extra),
        xh = local_triangulate(cSa(:, :, :, 1:nU), cQa(:, :, :, 1:nU), ST.NF, bet, FG);   % as in Code T
    end
    if true
        dx = abs(real(xh) - real(UEm(:)));  dx = min(dx, FG.sL - dx);
        dy = abs(imag(xh) - imag(UEm(:)));  dy = min(dy, FG.sL - dy);
        e = sqrt(dx.^2 + dy.^2);
        O.loc = e(isfinite(e));
    end
    a = aErr(isfinite(aErr));
    O.hit = NaN;  O.rmse = NaN;
    if ~isempty(a)
        hit = a < EP.tolS;
        O.hit = mean(hit);
        if any(hit), O.rmse = sqrt(mean(a(hit) .^ 2)); end
    end
    O.nmse = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);
    Drt = cluster_rate_rule(ST.R2, Cs, p, ST.Dcn);
    bE = zeros(4, K);  nE = zeros(4, K);  nb = 0;  ne = 0;
    for mc = 1:EP.nReal
        Hm = reshape(Ht(:, mc, :), [L * Na, K]);
        He = reshape(Hh(:, mc, :), [L * Na, K]);
        [bu, nu, bpu, epu] = metrics_case(He, Hm, Cs, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);
        bE = bE + bu;  nE = nE + nu;  nb = nb + bpu;  ne = ne + epu;
    end
    O.be  = sum(bE, 2) / max(K * nb, 1);
    seU   = EP.prelog * log2(1 + 1 ./ max(nE / max(ne, 1), 1 / (EP.nReal * EP.nSym)));
    O.se  = seU;
    O.sum = sum(seU, 2);
end

function [O, nf] = run_task_u(D, G, PD, ES, useCape, Kc, p, seed, CF)
% One task: one drop, one scheme, one SNR.
    ST = array_stats(G, D.d3, D.sE, D.phL, D.thL, D.bet, CF.lambda, CF.eps_NF, CF.asd_phi, CF.asd_th, CF.clu_div);
    nf = mean(ST.NF(:));
    FG = cape_cfg(true, D.APp, D.APw, CF.sL, CF.hD, CF.nAP, CF.nCand, CF.J, CF.step0, CF.nIter, CF.consist);
    FG.nAng = CF.nAng;  FG.nRng = CF.nRng;  FG.nPk = max(CF.nPeak, ceil(Kc / 2));  FG.nSubMax = CF.nSubMax;
    FG.J = ES.J;                                   % [v2] dictionary starts per link
    if strcmp(ES.kind, 'sma')
        ES.he = nf_channel(G.frames(:), D.d3, D.sE, D.bet, CF.lambda) .* reshape(ST.NF, [1, size(ST.NF)]);
    else
        ES.he = ST.hD;
    end
    ES.sTrue = D.sE;
    EP = struct('nReal', CF.nReal, 'nSym', CF.nSym, 'M', CF.M, 'dth', CF.dth, 'maxBr', CF.maxBr, ...
                'mod', CF.mod, 'prelog', 1 - Kc / CF.tau_c, 'eps_NF', CF.eps_NF, 'tolS', CF.tolS);
    O = eval_scheme_u(ST, D.bet, useCape, PD, FG, p, EP, D.UEm, seed, ES);
end

function xh = local_triangulate(sE, qE, NFm, beta, FG)
% [CODE T] Localisation WITHOUT CAPE (baseline): classical angle-of-arrival
% triangulation from hard per-AP decisions. Every one of the FG.nAP
% strongest near-field APs of a user gives its local estimate plus the same J
% strongest local lobes CAPE starts from. Each lobe gives a BEARING line (both
% mirror sides of the linear array), sampled from 20 m to 260 m, plus the
% point at its own range estimate. The CPU takes the point that is closest
% to the bearing lines of all the other APs (smallest sum of distances) and
% averages it with the nearest points of the other APs. Same APs, pilots and
% candidate lobes as CAPE: only DECISIONS are fused (geometry), not the pilot
% likelihoods, and there is no refinement on the pilots.
% Check (Octave, 20 drops): with the true angle among every AP's candidates
% the median error is 1.5-2.3 m; if the true lobe is present at only 80% /
% 50% / 30% of the APs (grating-lobe aliases otherwise) it grows to
% 64-78 / 92-122 / 101-131 m. Hard decisions break down under grating-lobe
% ambiguity; this is what the likelihood fusion of CAPE avoids.
    [~, K, M, nJ] = size(sE);
    xh = nan(K, M);
    sL = FG.sL;  rr = 20:6:260;
    td = @(a, b) sqrt(min(abs(real(a) - real(b)), sL - abs(real(a) - real(b))) .^ 2 + ...
                      min(abs(imag(a) - imag(b)), sL - abs(imag(a) - imag(b))) .^ 2);
    for k = 1:K
        nf = find(NFm(:, k));
        if isempty(nf), continue; end
        [~, ob] = sort(beta(nf, k), 'descend');
        A = nf(ob(1:min(FG.nAP, numel(nf))));
        for m = 1:M
            P = cell(numel(A), 1);  ok = true(numel(A), 1);
            for i = 1:numel(A)
                Pi = [];
                for j = 1:nJ
                    sv = sE(A(i), k, m, j);
                    if ~isfinite(sv), continue; end
                    Pi = [Pi; cand_positions(FG.AP(A(i)), sv, qE(A(i), k, m, j), sL, FG.hD); ...
                          bearing_points(FG.AP(A(i)), sv, rr, sL, FG.hD)]; %#ok<AGROW>
                end
                if isempty(Pi), ok(i) = false; continue; end
                P{i} = Pi;
            end
            P = P(ok);
            if isempty(P), continue; end
            if numel(P) == 1, xh(k, m) = P{1}(1); continue; end
            ref = vertcat(P{1:min(FG.nCand, numel(P))});
            sc = zeros(numel(ref), 1);  NB = zeros(numel(ref), numel(P));
            for i = 1:numel(P)
                [dm, jm] = min(td(ref(:), P{i}(:).'), [], 2);
                sc = sc + dm;  NB(:, i) = P{i}(jm);
            end
            [~, b] = min(sc);
            xb = [ref(b); NB(b, :).'];
            off = xb - xb(1);                              % average on the torus around ref(b)
            ox = mod(real(off) + sL / 2, sL) - sL / 2;  oy = mod(imag(off) + sL / 2, sL) - sL / 2;
            x = xb(1) + mean(ox) + 1j * mean(oy);
            xh(k, m) = mod(real(x), sL) + 1j * mod(imag(x), sL);
        end
    end
end

function X = bearing_points(APl, s, rr, sL, hD)
% Points along the bearing of an AP-local angle estimate s (effective sine of
% the channel model), at the 3-D ranges rr, on both mirror sides.
    d3 = max(rr(:), hD + 0.5);
    dH = sqrt(d3 .^ 2 - hD ^ 2);
    p1 = asin(max(min(s * d3 ./ dH, 1), -1));
    x = APl + [dH .* exp(1j * p1); dH .* exp(1j * (pi - p1))];
    X = mod(real(x), sL) + 1j * mod(imag(x), sL);
end

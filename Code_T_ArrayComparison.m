clear; close all; clc;
%% ====================================================================
%  CODE T : ARRAY COMPARISON AT EQUAL HARDWARE (built from Code R)
%  --------------------------------------------------------------------
%  Focused study requested after the Code R / Code S review. Compares, at
%  N = 32 physical antennas per AP and with CAPE channel estimation:
%    1. ULA, d = lambda/2            dense reference
%    2. ULA, d = 4 lambda            wide-spaced ULA (still a ULA)
%    3. Nested array                 physical processing
%    4. Nested array + co-array      SAME antennas as 3, co-array processing
%    5. Coprime array                physical processing
%    6. Coprime array + co-array     SAME antennas as 5, co-array processing
%  "Difference co-array" is a PROCESSING option of a sparse array, not a
%  separate array: schemes 3/4 and 5/6 use identical antennas. Every scheme
%  sees the same random drops, channels, pilots and noise in a setup
%  (paired comparison), so a 3-vs-4 or 5-vs-6 difference is caused only by
%  the co-array processing. The co-array MUSIC subarray is NOT capped here,
%  so the co-array gets its full smoothed virtual aperture.
%  Figures T1..T9. Results are saved to codeT_arraycomp_results.mat BEFORE
%  plotting; plot_codeT_figures.m redraws the figures without rerunning.
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
%  [CODE T] ARRAY SET (equal hardware: N0 = 32 physical antennas per AP)
%  --------------------------------------------------------------------
%  'ula'     dense ULA, spacing d = lambda/2 (reference)
%  'ulaw'    wide-spaced ULA, spacing d = ulaw_d_lambda * lambda. Uniform
%            spacing, so it is a ULA; wide spacing enlarges the aperture
%            with the same N but creates grating lobes (d > lambda/2).
%  'nested'  two-level nested array (Pal & Vaidyanathan, IEEE TSP 2010)
%  'coprime' extended coprime array (Pal & Vaidyanathan, DSP/SPE 2011):
%            N1 elements at spacing M*d0 plus 2M elements at spacing N1*d0,
%            M < N1 coprime, 2M + N1 - 1 = N0 (M = 8, N1 = 17 for N0 = 32)
%  Co-array processing (en_idea = 1) is only meaningful for nested/coprime:
%  a wide ULA has no hole-free co-array (all its lags are multiples of eta).
%% ====================================================================
ulaw_d_lambda = 4;                    % wide-ULA spacing in wavelengths
ulaw_eta      = 2 * ulaw_d_lambda;    % same spacing in units of d0 = lambda/2
ca_Nsub_max   = 1000;                 % no cap: full smoothed virtual aperture for the co-array

eqN_study   = true;
en_names    = {'ULA lambda/2 + CAPE', sprintf('ULA %g lambda + CAPE', ulaw_d_lambda), 'Nested + CAPE', ...
               'Nested + co-array + CAPE', 'Coprime + CAPE', 'Coprime + co-array + CAPE'};
en_leg      = {'ULA (\lambda/2)', sprintf('ULA (%g\\lambda)', ulaw_d_lambda), 'Nested', ...
               'Nested + co-array', 'Coprime', 'Coprime + co-array'};
en_geom     = {'ula', 'ulaw', 'nested', 'nested', 'coprime', 'coprime'};
en_idea     = [0 0 0 1 0 1];   % 1 = add difference co-array processing (same antennas)
en_cape     = [1 1 1 1 1 1];   % CAPE for every scheme: the array is the only variable
en_arr      = [1 2 3 3 4 4];   % physical array of each scheme (1 ULA, 2 wide ULA, 3 nested, 4 coprime)
nSchE       = numel(en_names);
en_caseName = {'Antenna-rich (K<N)', 'User-rich (K>N)'};
en_N0       = [32 32];         % 32 PHYSICAL antennas per AP, both cases
en_K        = [8 40];          % Case 1 K = 8 (K < N), Case 2 K = 40 (K > N)
nCaseE      = numel(en_N0);
nSetE       = 10;              % random AP/user drops: 6 = quick look (as Code R), >= 10 for the paper
nRealE      = 2;
nSymE       = 100;

% ---- prebuild geometry, polar dictionary and co-arrays per (scheme, case) ----
EN_geo = cell(nSchE, nCaseE);  EN_PD = cell(nSchE, nCaseE);
EN_CA  = cell(nSchE, nCaseE);  EN_CAs = cell(nSchE, nCaseE);
for caE = 1:nCaseE
    for scE = 1:nSchE
        Gtmp = geom_targetN(en_geom{scE}, en_N0(caE), lambda, ulaw_eta);
        EN_geo{scE, caE} = Gtmp;
        pdE = build_polar_dict(Gtmp.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
        pdE.smax = pd_smax;  pdE.stop = pd_stop;  pdE.nref = pd_nref;  pdE.nstart = pd_nstart;
        EN_PD{scE, caE} = pdE;
        EN_CA{scE, caE}  = coarray_setup(Gtmp.pos, lambda, 'diff');
        if en_idea(scE) == 3
            EN_CAs{scE, caE} = coarray_setup(Gtmp.pos, lambda, 'sum');
        else
            EN_CAs{scE, caE} = [];
        end
    end
end

EN_nmse = zeros(nSchE, nCaseE, nSNR, nSetE);
EN_be   = zeros(4, nSchE, nCaseE, nSNR, nSetE);
EN_sum  = zeros(4, nSchE, nCaseE, nSNR, nSetE);
EN_loc  = nan(nSchE, nCaseE, nSNR, nSetE);
EN_time = zeros(nSchE, nCaseE, nSetE);              % mean run time (complexity proxy)
EN_fh   = zeros(nSchE, nCaseE, nSetE);              % fused APs per user (fronthaul proxy)
EN_nf   = zeros(nSchE, nCaseE, nSetE);              % near-field link fraction
EN_seD  = nan(4, max(en_K), nSchE, nCaseE, nSetE);  % [CODE O] sized by the largest K

if eqN_study
fprintf('\n=== [CODE T] array comparison: %d schemes x %d cases x %d setups ===\n', nSchE, nCaseE, nSetE);
for caE = 1:nCaseE
    fprintf('  Case %d (%s): N0 = %d physical antennas, K = %d users\n', ...
            caE, en_caseName{caE}, en_N0(caE), en_K(caE));
    for scE = 1:nSchE
        fprintf('     %-24s %-10s N=%d  virtual=%d (contig %d)\n', en_names{scE}, en_geom{scE}, ...
                EN_geo{scE, caE}.N, EN_CA{scE, caE}.DOF, EN_CA{scE, caE}.Nv);
    end
end
r_max_en = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);
parfor ns = 1:nSetE
    s_nmse = zeros(nSchE, nCaseE, nSNR);
    s_be   = zeros(4, nSchE, nCaseE, nSNR);
    s_sum  = zeros(4, nSchE, nCaseE, nSNR);
    s_loc  = nan(nSchE, nCaseE, nSNR);
    s_time = zeros(nSchE, nCaseE);
    s_fh   = zeros(nSchE, nCaseE);
    s_nf   = zeros(nSchE, nCaseE);
    s_seD  = nan(4, max(en_K), nSchE, nCaseE);
    for caE = 1:nCaseE
        Kc  = en_K(caE);
        APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
        UEp = drop_scenario(APp, Kc, L, min(n_collinear_pairs, floor(Kc / 2)), ...
                            r_min_NF, r_max_en, 'hotspot', hs_r, squareLen);
        UEm = mod(real(UEp), squareLen) + 1j * mod(imag(UEp), squareLen);
        APw = wrap_images(APp, squareLen, wrap_fix);
        [d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
        EPc = struct('nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, ...
                     'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);
        for scE = 1:nSchE
            Gc = EN_geo{scE, caE};
            ST = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);
            s_nf(scE, caE) = mean(ST.NF(:));
            FG = cape_cfg(true, APp, APw, squareLen, hDiff, cape_nAP, cape_nCand, ...
                          cape_J, cape_step0, cape_nIter, cape_consist);
            FG.useCA   = (en_idea(scE) > 0);
            FG.co_idea = en_idea(scE);
            FG.CA      = EN_CA{scE, caE};
            FG.CAs     = EN_CAs{scE, caE};
            FG.nAng    = ca_nAngGrid;
            FG.nRng    = ca_nRangeGrid;
            FG.nPk     = max(ca_nPeak, ceil(Kc / 2));   % enough angle peaks when users are many
            FG.nSubMax = ca_Nsub_max;
            for si = 1:nSNR
                seed = 9e6 + 1e4 * ns + 1000 * caE + si;   % [CODE T] same seed for every scheme (paired)
                O = eval_scheme(ST, bet, en_cape(scE), EN_PD{scE, caE}, FG, SNR_lin(si), EPc, false, UEm, seed);
                s_nmse(scE, caE, si) = O.nmse;
                s_be(:, scE, caE, si) = O.be;
                s_sum(:, scE, caE, si) = O.sum;
                if ~isempty(O.loc), s_loc(scE, caE, si) = median(O.loc); end
                s_time(scE, caE) = s_time(scE, caE) + O.time;
                s_fh(scE, caE)   = s_fh(scE, caE) + O.nfus;
                if si == snr_ref_idx
                    s_seD(:, 1:size(O.se, 2), scE, caE) = O.se;
                end
            end
        end
    end
    EN_nmse(:, :, :, ns) = s_nmse;  EN_be(:, :, :, :, ns) = s_be;
    EN_sum(:, :, :, :, ns) = s_sum;    EN_loc(:, :, :, ns) = s_loc;
    EN_time(:, :, ns) = s_time / nSNR;  EN_fh(:, :, ns) = s_fh / nSNR;
    EN_nf(:, :, ns) = s_nf;  EN_seD(:, :, :, :, ns) = s_seD;
    fprintf('  equal-N study setup %d/%d done\n', ns, nSetE);
end
end

% ---- [CODE T] user-count sweep at the reference SNR, all six schemes ----
% N is fixed at en_N0(1) = 32, so the prebuilt geometry, polar dictionary and
% co-array of case 1 are reused (identical to rebuilding them).
Ksw     = [4 8 16 24 32 40 48];
Ksw_sch = 1:nSchE;
nKsw    = numel(Ksw);  nKsch = numel(Ksw_sch);
EN_Kbe  = zeros(nKsch, nKsw, nSetE);    % BER, Cross-AP detector
EN_Ksum = zeros(nKsch, nKsw, nSetE);    % sum SE, Cross-AP detector
if eqN_study
fprintf('\n=== [CODE T] user-count sweep at %d dB: K = %s ===\n', SNR_dB(snr_ref_idx), mat2str(Ksw));
parfor ns = 1:nSetE
    kbe = zeros(nKsch, nKsw);  ksum = zeros(nKsch, nKsw);
    for ik = 1:nKsw
        Kc  = Ksw(ik);
        APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
        UEp = drop_scenario(APp, Kc, L, min(n_collinear_pairs, floor(Kc / 2)), ...
                            r_min_NF, r_max_en, 'hotspot', hs_r, squareLen);
        UEm = mod(real(UEp), squareLen) + 1j * mod(imag(UEp), squareLen);
        APw = wrap_images(APp, squareLen, wrap_fix);
        [d3, sE, phL, thL, bet] = link_geometry(APw, UEp, hDiff, fc_GHz, hBS, hUT, c0, noisePow_dBm);
        EPc = struct('nReal', nRealE, 'nSym', nSymE, 'M', M_lst, 'dth', d_th, ...
                     'maxBr', maxBr, 'mod', modOrder, 'prelog', prelog, 'eps_NF', eps_NF);
        for jj = 1:nKsch
            scE = Ksw_sch(jj);
            Gc  = EN_geo{scE, 1};
            PDc = EN_PD{scE, 1};
            ST  = array_stats(Gc, d3, sE, phL, thL, bet, lambda, eps_NF, ASD_varphi, ASD_theta, clu_div);
            FG  = cape_cfg(true, APp, APw, squareLen, hDiff, cape_nAP, cape_nCand, ...
                           cape_J, cape_step0, cape_nIter, cape_consist);
            FG.useCA = (en_idea(scE) > 0);  FG.co_idea = en_idea(scE);
            FG.CA = EN_CA{scE, 1};  FG.CAs = [];
            FG.nAng = ca_nAngGrid;  FG.nRng = ca_nRangeGrid;
            FG.nPk = max(ca_nPeak, ceil(Kc / 2));  FG.nSubMax = ca_Nsub_max;
            seed = 3e7 + 1e4 * ns + 100 * ik;   % same seed for every scheme (paired)
            O = eval_scheme(ST, bet, en_cape(scE), PDc, FG, SNR_lin(snr_ref_idx), EPc, false, UEm, seed);
            kbe(jj, ik)  = O.be(4);
            ksum(jj, ik) = O.sum(4);
        end
    end
    EN_Kbe(:, :, ns) = kbe;  EN_Ksum(:, :, ns) = ksum;
    fprintf('  user-count sweep setup %d/%d done\n', ns, nSetE);
end
end

if eqN_study
    % [CODE T] save BEFORE plotting, so a plotting problem can never lose the run
    save('codeT_arraycomp_results.mat', 'SNR_dB', 'snr_ref_idx', 'lambda', 'L', ...
         'en_names', 'en_leg', 'en_geom', 'en_idea', 'en_cape', 'en_arr', 'en_caseName', ...
         'en_N0', 'en_K', 'ulaw_d_lambda', 'nSetE', 'EN_geo', ...
         'EN_nmse', 'EN_be', 'EN_sum', 'EN_loc', 'EN_time', 'EN_fh', 'EN_nf', 'EN_seD', ...
         'Ksw', 'Ksw_sch', 'EN_Kbe', 'EN_Ksum');
    fprintf('\nResults saved to codeT_arraycomp_results.mat\n');

    %% ===================== [CODE T] FIGURES T1-T9 =====================
    %  Every scheme uses CAPE. Colour = physical array. Solid line + filled
    %  marker = processing on the PHYSICAL array. Dashed line + open marker =
    %  the SAME antennas with difference co-array processing added.
    tNs   = numel(en_names);  tNc = numel(en_N0);  tNsnr = numel(SNR_dB);
    tSet  = size(EN_be, 5);
    [~, tSi] = min(abs(SNR_dB - 10));             % 10 dB reference point
    tDet  = 4;                                     % Cross-AP List detector
    tArrNm  = {'ULA (\lambda/2)', en_leg{2}, 'Nested', 'Coprime'};
    tArrSch = [1 2 3 5];                           % one scheme per physical array
    tArrCol = {[0.35 0.35 0.35], [0.00 0.45 0.00], [0.00 0.10 0.55], [0.50 0.10 0.60]};
    tCol  = tArrCol(en_arr);
    tLs   = repmat({'-'}, 1, tNs);   tLs(en_idea > 0) = {'--'};
    tMk   = {'d', 'o', 's', 's', '^', '^'};
    tMfc  = tCol;                    tMfc(en_idea > 0) = {'none'};
    tPair = [3 4; 5 6];                            % [physical, co-array] pairs
    tPairNm = {'Nested', 'Coprime'};
    mNM = mean(EN_nmse, 4);  mBE = mean(EN_be, 5);  mSU = mean(EN_sum, 5);
    mLC = nan(tNs, tNc, tNsnr);
    for sc = 1:tNs
        for ca = 1:tNc
            for si = 1:tNsnr
                v = EN_loc(sc, ca, si, :);  v = v(isfinite(v));
                if ~isempty(v), mLC(sc, ca, si) = median(v); end
            end
        end
    end

    % ---- T1: physical arrays and their hole-free difference co-arrays ----
    figure('Name', 'T1 Physical arrays and difference co-arrays', 'Position', [30 30 1100 760]);
    for a = 1:4
        pa = EN_geo{tArrSch(a), 1}.pos(:);
        lg = unique(round((pa - pa.') / (lambda / 2)));
        bk = [0; find(diff(lg) ~= 1); numel(lg)];  [nv, mi] = max(diff(bk));
        cont = lg(bk(mi) + 1:bk(mi + 1));
        ap = max(pa) - min(pa);
        subplot(4, 1, a); hold on; box on; grid on;
        plot(pa, ones(size(pa)), 'o', 'Color', tArrCol{a}, 'MarkerFaceColor', tArrCol{a}, 'MarkerSize', 5);
        plot(cont * lambda / 2, 0.6 * ones(size(cont)), '.', 'Color', tArrCol{a}, 'MarkerSize', 6);
        ylim([0.3 1.3]); set(gca, 'YTick', [0.6 1], 'YTickLabel', {'virtual', 'physical'});
        xlabel('position along the array [m]');
        title(sprintf('%s:  N = %d,  aperture = %.2f m,  hole-free co-array = %d lags,  Rayleigh = %.0f m', ...
              tArrNm{a}, numel(pa), ap, nv, 2 * ap^2 / lambda));
    end

    % ---- T2: far-field beampattern (where the grating lobes are) ----
    figure('Name', 'T2 Far-field beampattern (grating lobes)', 'Position', [30 30 1100 650]);
    sg = linspace(-1, 1, 4001);
    for a = 1:4
        pa = EN_geo{tArrSch(a), 1}.pos(:).';
        B = abs(exp(1j * 2 * pi * (sg.') * pa / lambda) * ones(numel(pa), 1)) / numel(pa);
        subplot(2, 2, a); plot(sg, 20 * log10(max(B, 1e-4)), 'Color', tArrCol{a}, 'LineWidth', 1.2);
        grid on; box on; ylim([-40 0]); xlabel('sin(\theta)'); ylabel('beampattern [dB]'); title(tArrNm{a});
    end

    % ---- T3-T5: NMSE, BER and sum SE vs SNR (all six schemes, both cases) ----
    tMet = {mNM, mBE, mSU};
    tNm  = {'T3 Channel NMSE vs SNR', 'T4 BER vs SNR', 'T5 Sum SE vs SNR'};
    tYl  = {'channel NMSE', 'BER (Cross-AP List)', 'sum SE [bit/s/Hz]'};
    tLoc = {'southwest', 'southwest', 'northwest'};
    for f = 1:3
        figure('Name', tNm{f}, 'Position', [30 30 1150 470]);
        for ca = 1:tNc
            subplot(1, 2, ca); hold on; box on; grid on;
            for sc = 1:tNs
                if f == 1
                    y = squeeze(tMet{f}(sc, ca, :));
                else
                    y = squeeze(tMet{f}(tDet, sc, ca, :));
                end
                if f == 2, y = max(y, 1e-6); end
                plot(SNR_dB, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
                     'LineWidth', 1.5, 'DisplayName', en_leg{sc});
            end
            if f < 3, set(gca, 'YScale', 'log'); end
            xlabel('SNR [dB]'); ylabel(tYl{f});
            title(sprintf('%s (K = %d), all with CAPE', en_caseName{ca}, en_K(ca)));
            legend('Location', tLoc{f});
        end
    end

    % ---- T6: BER and sum SE vs number of users K at 10 dB ----
    mKbe = mean(EN_Kbe, 3);  mKsum = mean(EN_Ksum, 3);
    figure('Name', 'T6 BER and sum SE vs number of users', 'Position', [30 30 1150 470]);
    for f = 1:2
        subplot(1, 2, f); hold on; box on; grid on;
        for j = 1:numel(Ksw_sch)
            sc = Ksw_sch(j);
            if f == 1, y = max(mKbe(j, :), 1e-6); else, y = mKsum(j, :); end
            plot(Ksw, y, [tLs{sc} tMk{sc}], 'Color', tCol{sc}, 'MarkerFaceColor', tMfc{sc}, ...
                 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
        end
        xlabel('number of users K');
        if f == 1
            set(gca, 'YScale', 'log'); ylabel('BER (Cross-AP List)');
            title(sprintf('BER vs K at %d dB, N = %d per AP', SNR_dB(tSi), en_N0(1)));
        else
            ylabel('sum SE [bit/s/Hz]');
            title(sprintf('Sum SE vs K at %d dB, N = %d per AP', SNR_dB(tSi), en_N0(1)));
        end
        legend('Location', 'northwest');
    end

    % ---- T7: does co-array processing help? paired difference, same antennas ----
    figure('Name', 'T7 Effect of co-array processing (paired)', 'Position', [30 30 1150 760]);
    tYl7 = {'BER change with co-array [%]  (<0 = better)', 'sum-SE change with co-array [%]  (>0 = better)'};
    for row = 1:2
        for ca = 1:tNc
            subplot(2, 2, (row - 1) * 2 + ca); hold on; box on; grid on;
            plot(SNR_dB, zeros(size(SNR_dB)), 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            ttl = '';
            for pr = 1:2
                a = tPair(pr, 1);  b = tPair(pr, 2);
                if row == 1
                    X = EN_be(tDet, a, ca, :, :);   Y = EN_be(tDet, b, ca, :, :);
                else
                    X = EN_sum(tDet, a, ca, :, :);  Y = EN_sum(tDet, b, ca, :, :);
                end
                X = reshape(X, tNsnr, tSet);  Y = reshape(Y, tNsnr, tSet);  Dd = Y - X;
                ref = max(abs(mean(X, 2)), eps);
                mu  = 100 * mean(Dd, 2) ./ ref;
                se  = 100 * std(Dd, 0, 2) / sqrt(tSet) ./ ref;
                hEb = errorbar(SNR_dB, mu, se, ['-' tMk{a}]);
                set(hEb, 'Color', tCol{a}, 'MarkerFaceColor', tCol{a}, 'LineWidth', 1.4, 'DisplayName', tPairNm{pr});
                ttl = [ttl sprintf('  %s %d%%', tPairNm{pr}, round(100 * mean(Dd(:) == 0)))]; %#ok<AGROW>
            end
            xlabel('SNR [dB]'); ylabel(tYl7{row});
            title({sprintf('%s (K = %d)', en_caseName{ca}, en_K(ca)), ['identical:' ttl]});
            legend('Location', 'best');
        end
    end

    % ---- T8: localisation error vs SNR ----
    figure('Name', 'T8 Localisation error vs SNR', 'Position', [30 30 1150 470]);
    for ca = 1:tNc
        subplot(1, 2, ca); hold on; box on; grid on;
        for sc = 1:tNs
            plot(SNR_dB, squeeze(mLC(sc, ca, :)), [tLs{sc} tMk{sc}], 'Color', tCol{sc}, ...
                 'MarkerFaceColor', tMfc{sc}, 'LineWidth', 1.5, 'DisplayName', en_leg{sc});
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('median position error [m]');
        title(sprintf('%s (K = %d)', en_caseName{ca}, en_K(ca))); legend('Location', 'northeast');
    end

    % ---- T9: array hardware metrics ----
    tAp = zeros(4, 1);  tUd = zeros(4, 1);  tW1 = zeros(4, 1);  tCl = zeros(4, 1);
    for a = 1:4
        pa = EN_geo{tArrSch(a), 1}.pos(:);
        Dm = round((pa - pa.') / (lambda / 2));
        lg = unique(Dm(:));  bk = [0; find(diff(lg) ~= 1); numel(lg)];
        tAp(a) = max(pa) - min(pa);  tUd(a) = max(diff(bk));  tW1(a) = sum(Dm(:) == 1);
        U = zeros(numel(pa));  nz = (Dm ~= 0) & (abs(Dm) <= 100);
        U(nz) = 1 ./ abs(Dm(nz));  U(Dm == 0) = 1;
        tCl(a) = norm(U - diag(diag(U)), 'fro') / max(norm(U, 'fro'), eps);
    end
    figure('Name', 'T9 Array hardware metrics', 'Position', [30 30 1000 720]);
    tVal = {tAp, tUd, tW1, tCl};
    tYl9 = {'aperture [m]', 'hole-free co-array length', 'w(1): pairs \lambda/2 apart', 'coupling leakage'};
    tTt9 = {'Aperture (resolution)', 'Hole-free difference co-array', 'Closest pairs', 'Mutual coupling'};
    for m = 1:4
        subplot(2, 2, m); hold on; box on; grid on;
        for a = 1:4
            bar(a, tVal{m}(a), 0.6, 'FaceColor', tArrCol{a});
        end
        set(gca, 'XTick', 1:4, 'XTickLabel', tArrNm);  xlim([0.4 4.6]);
        ylabel(tYl9{m});  title(tTt9{m});
    end

    % ---- console summary at 10 dB ----
    for ca = 1:tNc
        fprintf('\n=== [CODE T] %s, K = %d, at %d dB (mean +- std. error over %d setups) ===\n', ...
                en_caseName{ca}, en_K(ca), SNR_dB(tSi), tSet);
        fprintf('%-28s %10s %22s %18s %10s\n', 'scheme', 'NMSE', 'BER (Cross-AP List)', 'sum SE', 'loc [m]');
        for sc = 1:tNs
            be = reshape(EN_be(tDet, sc, ca, tSi, :), [], 1);
            su = reshape(EN_sum(tDet, sc, ca, tSi, :), [], 1);
            fprintf('%-28s %10.2e %11.3e +- %8.1e %8.2f +- %5.2f %10.3f\n', en_names{sc}, mNM(sc, ca, tSi), ...
                    mean(be), std(be) / sqrt(tSet), mean(su), std(su) / sqrt(tSet), mLC(sc, ca, tSi));
        end
        for pr = 1:2
            a = tPair(pr, 1);  b = tPair(pr, 2);
            Dd = EN_be(tDet, b, ca, :, :) - EN_be(tDet, a, ca, :, :);
            fprintf('  %s: co-array vs physical processing identical in %d%% of (SNR, setup) points\n', ...
                    tPairNm{pr}, round(100 * mean(Dd(:) == 0)));
        end
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
    for rnd = 1:PD.nref
        sc = mod(s0 + ds * (-2:2) + 1, 2) - 1;
        if useKappa
            kc = max(k0 + dq * (-2:2), 0);
            [SS, KK] = meshgrid(sc, kc);
            SS = SS(:).';
            QQ = min(KK(:).' ./ max(1 - SS .^ 2, 1e-6), PD.q_max);
        else
            qc = min(max(q0 + dq * (-2:2), 0), PD.q_max);
            [SS, QQ] = meshgrid(sc, qc);
            SS = SS(:).';
            QQ = QQ(:).';
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
    % Per-user BER and per-user hard-decision symbol NMSE for all four
    % detectors. Detection logic identical to ber_case; only the accounting
    % is per user. Feeds figures 6, 7 and 8.  [UNCHANGED]
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
    cnt_stg = 0;            % [CODE C] SIC stages seen by the cross-AP variant
    cnt_shadXap = 0;        % [CODE C] stages where the all-AP fusion is invoked

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

    berU = zeros(4, K);
    nmseU = zeros(4, K);
    bits_pu = nSym * kbit;
    eng_pu = nSym;
    for t = 1:nSym
        ti = randi(nC, K, 1);
        s = Cq(ti);
        n = (randn(LN, 1) + 1j * randn(LN, 1)) / sqrt(2);
        y = sqrt(p) * (Htrue * s) + n;

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
            berU(1, k) = berU(1, k) + sum(Bmap(mi, :) ~= Bmap(ti(k), :));
            nmseU(1, k) = nmseU(1, k) + abs(Cq(mi) - s(k))^2;
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
            berU(2, k) = berU(2, k) + sum(Bmap(di(k), :) ~= Bmap(ti(k), :));
            nmseU(2, k) = nmseU(2, k) + abs(Cq(di(k)) - s(k))^2;
        end

        for variant = 1:2
            useX = (variant == 2);
            row = variant + 2;
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
                        oa = zeros(L, 1);
                        for l = 1:L
                            idx = (l - 1) * N + 1:l * N;
                            oa(l) = vA{i}(:, l)' * yrp(idx);
                        end
                        sh_f = aA{i}' * oa;
                        [dsort_xap, ix_xap] = sort(abs(sh_f - Cq));
                        shad_xap = (dsort_xap(1) >= dth_abs);
                        cnt_stg = cnt_stg + 1;
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
                        ix_pool = [ix_xap(1:nc_xap); ix_cl(1:nc_cl)];
                        d_pool = [dsort_xap(1:nc_xap); dsort_cl(1:nc_cl)];
                        [ix, ui] = unique(ix_pool, 'stable');
                        dsort = d_pool(ui);
                        nc = numel(ix);
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
            for k = 1:K
                [~, mk] = min(abs(bd(k) - Cq));
                berU(row, k) = berU(row, k) + sum(Bmap(mk, :) ~= Bmap(ti(k), :));
                nmseU(row, k) = nmseU(row, k) + abs(Cq(mk) - s(k))^2;
            end
        end
    end
    etaXap = cnt_shadXap / max(cnt_stg, 1);   % [CODE C] cross-AP gate firing rate
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

function G = geom_targetN(type, N0, lambda, eta)
% [CODE G] build each array type with EXACTLY N0 PHYSICAL elements on the
% lambda/2 grid, so every scheme is compared at equal hardware. The virtual
% (co-)array that results is whatever each geometry yields; only the physical
% count is fixed. 'eca'/'symnested' are trimmed to N0 keeping both ends;
% [CODE T] 'ulaw' is a wide-spaced ULA (spacing eta*d0); 'coprime' is the
% extended coprime array of Pal & Vaidyanathan.
    if nargin < 4 || isempty(eta), eta = 8; end
    d0 = lambda / 2;
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
            n = (0:N0 - 1)' * eta;
        case 'coprime'
            % [CODE T] extended coprime array (Pal & Vaidyanathan 2011): Nc
            % elements at spacing Mc*d0 and 2*Mc elements at spacing Nc*d0,
            % Mc < Nc coprime, sharing the origin -> 2*Mc + Nc - 1 = N0. Among
            % the valid pairs the one with the longest hole-free co-array
            % (Mc*Nc + Mc - 1 lags each side) is used: (8, 17) for N0 = 32.
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
    p = n * d0;
    G.pos = p - (max(p) + min(p)) / 2;
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
function CA = coarray_setup(pos, lambda, mode)
% Co-array of a linear array on the lambda/2 grid.
%   mode = 'diff' (default): lags x_m - x_n  -> DIFFERENCE co-array (the standard
%          Pinto building block; comes from the Hermitian covariance R = E[y y^H]).
%   mode = 'sum'          : lags x_m + x_n  -> SUM co-array (used by Idea 3; comes
%          from the complementary moment Rt = E[y y^T]). For near-field phase and
%          deterministic pilots the sum lags carry extra (range) structure.
    if nargin < 3, mode = 'diff'; end
    d0 = lambda / 2;
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
    Cs0 = zeros(size(Cs));
    bE = zeros(4, K);  nE = zeros(4, K);  bP = zeros(4, K);  nb = 0;  ne = 0;  eX = 0;
    for mc = 1:EP.nReal
        Hm = reshape(Ht(:, mc, :), [L * Na, K]);
        He = reshape(Hh(:, mc, :), [L * Na, K]);
        [bu, nu, bpu, epu, ex] = metrics_case(He, Hm, Cs, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);
        bE = bE + bu;  nE = nE + nu;  nb = nb + bpu;  ne = ne + epu;  eX = eX + ex;
        bu = metrics_case(Hm, Hm, Cs0, Drt, p, Na, L, K, EP.nSym, EP.M, EP.dth, EP.maxBr, EP.mod);
        bP = bP + bu;
    end
    O.be  = sum(bE, 2) / max(K * nb, 1);
    O.bp  = sum(bP, 2) / max(K * nb, 1);
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

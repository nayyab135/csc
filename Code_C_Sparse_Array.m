
clear;
close all;
clc;
%% ====================================================================
%  CODE C : NON-UNIFORM SPARSE ARRAYS IN HYBRID NEAR-/FAR-FIELD CELL-FREE
%           XL-MIMO WITH CROSS-AP LIST-BASED SIC  (no IDD)
%  ------------------------------------------------------------------
%  Built on the authentic-CSI code (Without_IDD_Code_Authentic_CSI.m).
%  Clustering is the same as in that code (channel-norm and rate rules);
%  the P1 clustering of Code B is NOT included here.
%
%  What is new in Code C
%   1. Array geometry at every AP is selectable (array_type): dense ULA,
%      uniform sparse array (USA) or the NON-UNIFORM extended coprime array
%      (ECA) of Zhou et al. Default: ECA, 63 elements, aperture 7.6 m.
%      All channel models, covariances, the polar-domain dictionary and the
%      Rayleigh distance follow the true element positions.
%   2. Hybrid NF/FF model kept: a link is near-field when d < 2 D^2 / lambda
%      of the ARRAY IN USE, so a larger sparse aperture moves links into the
%      near field. The deployment-area sweep is extended so the far-field
%      links come back for the large sparse aperture.
%   3. Sparse-array comparison study (same drops, same fading, four arrays):
%      grating-lobe / ambiguity analysis per AP and for the whole network,
%      probability that a user pair is unresolvable at the home AP, inside
%      the serving cluster, at the OUTER (non-serving) APs and network-wide,
%      the "rescue" probability (unresolvable in the cluster but resolvable
%      at the outer APs, which is exactly what cross-AP SIC exploits),
%      estimation NMSE, BER of Linear / Hard-SIC / List-SIC / List+CrossAP
%      under estimated and perfect CSI, near-field vs far-field channel model,
%      cross-AP gate rate and fronthaul, information surplus, per-user SE.
%   4. Analytical figures: geometries, beam patterns with grating lobes and
%      their refocusing distance (Zhou Lemmas 5-6), beam depth, NF and
%      focusing fractions vs area, co-angular correlation vs range, and
%      ambiguity maps (single AP vs network).
%   New figures are named SPA-1 ... SPA-14 (SPA-14 is a neutral sparse-vs-
%   dense summary: whether sparse wins depends on the regime, see its comment); the old figures are unchanged in
%   meaning and now use the selected array.
%
%  Two corrections relative to the base code (both switchable):
%   * wrap_fix  : the base code built wrap-around images as
%                 wr(:)' + 1j*(wr(:)'), which only produces the 3 DIAGONAL
%                 images. Code C uses the full 3 x 3 torus (imaginary part
%                 from the transpose), as in Bjornson & Sanguinetti.
%   * The far-field covariance is built from element positions with the SAME
%     phase convention as the NUSW near-field channel. The base code's
%     Toeplitz build had the opposite sign (FF users appeared mirrored with
%     respect to NF users at the same AP); for a uniform dense array this
%     only mirrors the angle, for a non-uniform array Toeplitz does not hold.
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
array_type = 'eca';
N_ref      = 64;     % elements of the dense ULA and of the USA
usa_eta    = 2;      % USA spacing in units of lambda/2
eca_M      = 5;      % coprime pair (M, N), M > N, gcd = 1
eca_N      = 4;
eca_Lc     = 8;      % number of basic coprime arrays (even)
eca_Lc_sa  = 4;      % basic arrays of the same-aperture ECA
wrap_fix   = true;   % full 3 x 3 wrap-around (see header)
APG  = array_geometry(array_type, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa);
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
%  SETUP LOOP  [PARALLELIZED]
%% ====================================================================
parfor ns = 1:nSetups
    

    fprintf('Setup %d/%d  ', ns, nSetups);
    rs_bs = RandStream('mt19937ar', 'Seed', 777 + ns);

    % ------------------------------------------------------------------
    %  AP and UE placement  [UNCHANGED]
    % ------------------------------------------------------------------
    APpos = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;

    UEpos = zeros(K, 1);
    collinear_mask = false(K, 1);
    user_idx = 1;
    for pr = 1:n_collinear_pairs
        ap_l = mod(pr - 1, L) + 1;
        phi = 2 * pi * rand;
        r_near = r_min_NF + (r_max_NF / 2 - r_min_NF) * rand;
        r_far = r_max_NF / 2 + (r_max_NF - r_max_NF / 2) * rand;
        UEpos(user_idx)  = APpos(ap_l) + r_near * exp(1j * phi);
        UEpos(user_idx + 1) = APpos(ap_l) + r_far * exp(1j * phi);
        collinear_mask(user_idx) = true;
        collinear_mask(user_idx + 1) = true;
        user_idx = user_idx + 2;
    end
    ap_rr = 1;
    while user_idx <= K
        r = r_min_NF + (r_max_NF - r_min_NF) * rand;
        phi = 2 * pi * rand;
        UEpos(user_idx) = APpos(ap_rr) + r * exp(1j * phi);
        user_idx = user_idx + 1;
        ap_rr = mod(ap_rr, L) + 1;
    end

    APwrap = wrap_images(APpos, squareLen, wrap_fix);   % [CODE C] full 3 x 3 wrap-around

    % ------------------------------------------------------------------
    %  [NEW] PILOT ASSIGNMENT
    %  'legacy'  : the original mod() assignment, old figures unchanged.
    %  'aligned' : the two members of a collinear pair share a pilot, so
    %              near field range resolution has something to resolve.
    %              FAVOURABLE. Report it as such.  (see warning W3)
    %  'random'  : neutral control.
    % ------------------------------------------------------------------
    if strcmp(pilot_mode, 'orthogonal')
        pilotIndex = (1:K)';                 % one pilot per user, tau_p = K
    elseif strcmp(pilot_mode, 'aligned')
        pilotIndex = zeros(K, 1);
        nxt = 1;
        for pr = 1:n_collinear_pairs
            pilotIndex(2 * pr - 1) = nxt;
            pilotIndex(2 * pr)     = nxt;
            nxt = nxt + 1;
        end
        for k = 2 * n_collinear_pairs + 1:K
            pilotIndex(k) = nxt;
            nxt = nxt + 1;
        end
        pilotIndex = mod(pilotIndex - 1, tau_p) + 1;
    elseif strcmp(pilot_mode, 'random')
        prm = randperm(K);
        pilotIndex = zeros(K, 1);
        pilotIndex(prm) = mod((0:K - 1)', tau_p) + 1;
    else
        pilotIndex = mod((0:K - 1)', tau_p) + 1;
    end
    if contamination_mode || ~genie_csi
        pilotSym = pilotIndex;
        pilotF2  = pilotIndex;
    else
        pilotSym = mod((0:K - 1)', tau_sym) + 1;
        pilotF2  = mod((0:K - 1)', tau_fig2) + 1;
    end

    % ------------------------------------------------------------------
    %  Cell-free channel statistics (Cases 1 & 2)  [UNCHANGED]
    % ------------------------------------------------------------------
    beta     = zeros(L, K);
    R1       = zeros(N, N, L, K);
    R2       = zeros(N, N, L, K);
    h_det_all = zeros(N, L, K, 'like', 1 + 1j);
    NF_mask  = false(L, K);
    nNF_this = 0;

    for l = 1:L
        for k = 1:K
            [dH, wIdx] = min(abs(APwrap(l, :) - UEpos(k)));
            d3D = max(sqrt(hDiff^2 + dH^2), 10);
            d2D = max(dH, 1);
            phi_lk = angle(UEpos(k) - APwrap(l, wIdx));
            theta = asin(hDiff / d3D);
            sin_eff = sin(phi_lk) * cos(theta);

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
            beta(l, k) = db2pow(-PL_dB - noisePow_dBm);

            % [CODE C] covariance from the true element positions (any geometry)
            R_FF = far_corr(APG.pos / lambda, beta(l, k), sin_eff, phi_lk, theta, ASD_varphi, ASD_theta);
            R1(:, :, l, k) = R_FF;

            if d3D < d_Ray
                NF_mask(l, k) = true;
                nNF_this = nNF_this + 1;
                r_n = sqrt(d3D^2 + APG.pos.^2 - 2 * d3D * APG.pos * sin_eff);   % [CODE C] true positions
                h_det = sqrt(beta(l, k)) * (d3D ./ r_n) .* exp(-1j * 2 * pi * r_n / lambda);
                h_det_all(:, l, k) = h_det;
                R2(:, :, l, k) = h_det * h_det' + eps_NF * beta(l, k) * eye(N);
            else
                R2(:, :, l, k) = R_FF;
            end
        end
    end

    nNF_per_UE  = sum(NF_mask, 1);
    NF_dom_users = nNF_per_UE >= L / 2;
    FF_dom_users = ~NF_dom_users;
    nNF_total   = nNF_total + nNF_this;


    % ------------------------------------------------------------------
    %  CLUSTERING METRICS + CHANNEL-NORM CLUSTER  [UNCHANGED]
    % ------------------------------------------------------------------
    beta_NFeff = zeros(L, K);
    metric_CN  = zeros(L, K);
    for l = 1:L
        for k = 1:K
            if NF_mask(l, k)
                metric_CN(l, k) = real(h_det_all(:, l, k)' * h_det_all(:, l, k));
                beta_NFeff(l, k) = metric_CN(l, k) / N;
            else
                metric_CN(l, k) = beta(l, k);
                beta_NFeff(l, k) = beta(l, k);
            end
        end
    end
    beta_tilde = metric_CN;

    D_CN = false(L, K);
    for k = 1:K
        nf_l = NF_mask(:, k);
        if any(nf_l)
            thr_NF = max(metric_CN(nf_l, k)) / clu_div;
        else
            thr_NF = inf;
        end
        thr_FF = max(beta(:, k)) / clu_div;
        for l = 1:L
            if NF_mask(l, k)
                D_CN(l, k) = metric_CN(l, k) >= thr_NF;
            else
                D_CN(l, k) = beta(l, k)      >= thr_FF;
            end
        end
        [~, lm] = max(metric_CN(:, k));
        D_CN(lm, k) = true;
    end

    % ------------------------------------------------------------------
    %  CENTRALIZED BS CHANNEL STATISTICS  [UNCHANGED]
    % ------------------------------------------------------------------
    BS_pos   = (0.5 + 0.5j) * squareLen;
    n_idx_BS = (-(N_BS - 1) / 2:(N_BS - 1) / 2)';

    beta_BS = zeros(1, K);
    R_BS    = zeros(N_BS, N_BS, K);

    for k = 1:K
        d2D_BS = max(abs(BS_pos - UEpos(k)), 1);
        d3D_BS = max(sqrt(hDiff^2 + d2D_BS^2), 10);
        phi_BS = angle(UEpos(k) - BS_pos);
        theta_BS = asin(hDiff / d3D_BS);
        sin_eff_BS = sin(phi_BS) * cos(theta_BS);

        if d2D_BS <= 18
            P_LoS_BS = 1;
        else
            P_LoS_BS = (18 / d2D_BS) * (1 - exp(-d2D_BS / 63)) + exp(-d2D_BS / 63);
        end
        dBP_BS = 4 * (hBS - 1) * (hUT - 1) * fc_GHz * 1e9 / c0;
        if d3D_BS < dBP_BS
            PL_LoS_BS = 28 + 37 * log10(d3D_BS) + 20 * log10(fc_GHz);
        else
            PL_LoS_BS = 28 + 20 * log10(d3D_BS) + 20 * log10(fc_GHz) ...
          - 9 * log10(dBP_BS^2 + (hBS - hUT)^2);
        end
        if rand < P_LoS_BS
            PL_dB_BS = PL_LoS_BS + 4 * randn;
        else
            PL_NLoS_BS = 13.54 + 39.08 * log10(d3D_BS) + 20 * log10(fc_GHz) - 0.6 * (hUT - 1.5);
            PL_dB_BS = max(PL_LoS_BS, PL_NLoS_BS) + 6 * randn;
        end
        beta_BS(k) = db2pow(-PL_dB_BS - noisePow_dBm);

        ang = randn(rs_bs, 2, S_ang);
        PHI = phi_BS + ASD_varphi * ang(1, :);
        TH = theta_BS + ASD_theta * ang(2, :);
        A = exp(1j * 2 * pi * d_H * (px_BS * (sin(PHI) .* cos(TH)) + py_BS * sin(TH)));
        R_BS(:, :, k) = beta_BS(k) * (A * A') / S_ang;
    end

    fprintf('NF pairs:%d/%d(%.0f%%) | NF users:%d/%d | CL:%d | <|M_CN|>=%.1f\n', ...
            nNF_this, L * K, 100 * nNF_this / (L * K), sum(NF_dom_users), K, sum(collinear_mask), ...
            mean(sum(D_CN, 1)));

    % ------------------------------------------------------------------
    %  Precompute matrix square roots once per setup  [UNCHANGED]
    % ------------------------------------------------------------------
    Rs1 = zeros(N, N, L, K, 'like', 1 + 1j);
    Rs2 = zeros(N, N, L, K, 'like', 1 + 1j);
    for l = 1:L
        for k = 1:K
            Rs1(:, :, l, k) = sqrtm(R1(:, :, l, k));
            if ~NF_mask(l, k)
                Rs2(:, :, l, k) = sqrtm(R2(:, :, l, k));
            end
        end
    end
    Rs_BS = zeros(N_BS, N_BS, K, 'like', 1 + 1j);
    for k = 1:K
        Rs_BS(:, :, k) = sqrtm(R_BS(:, :, k));
    end


    %% ================================================================
    %  SNR LOOP
    %% ================================================================
    for si = 1:nSNR
        p = (SNR_lin(si));
        p3 = 4 * p;
        % --------------------------------------------------------------
        %  CASES 1 AND 2: cell-free processing  [UNCHANGED]
        % --------------------------------------------------------------
        for caseID = 1:2
            if caseID == 1
                R_use = R1;
                Rs_use = Rs1;
            else
                R_use = R2;
                Rs_use = Rs2;
            end

            H = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for k = 1:K
                    if caseID == 2 && NF_mask(l, k)
                        H(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                    else
                        H(idx, :, k) = sqrt(0.5) * Rs_use(:, :, l, k) * H(idx, :, k);
                    end
                end
            end

            % Case 1 is the far-field-only model (no NF pairs); Case 2 is hybrid.
            if caseID == 1
                NFm_c = false(L, K);
            else
                NFm_c = NF_mask;
            end
            [Hhat, C] = csi_acquire(H, R_use, NFm_c, h_det_all, beta, pilotIndex, tau_p, p, ...
                                    N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);

            SE_lmmse = zeros(K, nReal);
            SE_cnsic = zeros(K, nReal);
            BER_lmmse = zeros(K, nReal);
            BER_cnsic = zeros(K, nReal);
            SE_lmmse_pf = zeros(K, nReal);   % [PERFECT CSI]
            SE_cnsic_pf = zeros(K, nReal);

            servFull = true(L, K);
            C0 = zeros(size(C));             % [PERFECT CSI] zero error covariance
            for mc = 1:nReal
                Hmc    = reshape(H(:, mc, :),   [LN, K]);
                Hhat_mc = reshape(Hhat(:, mc, :), [LN, K]);

                [se_l, be_l] = det_local(Hhat_mc, Hmc, C, servFull, p, prelog, eta_FH, N, L, K);
                SE_lmmse(:, mc) = se_l;
                BER_lmmse(:, mc) = be_l;

                [se_s, be_s] = det_local_sic(Hhat_mc, Hmc, C, servFull, p, prelog, eta_FH, N, L, K);
                SE_cnsic(:, mc) = se_s;
                BER_cnsic(:, mc) = be_s;

                % [PERFECT CSI] true channel in the combiner, C = 0 (genie bound)
                SE_lmmse_pf(:, mc) = det_local(Hmc, Hmc, C0, servFull, p, prelog, eta_FH, N, L, K);
                SE_cnsic_pf(:, mc) = det_local_sic(Hmc, Hmc, C0, servFull, p, prelog, eta_FH, N, L, K);
            end

            sv = sum(mean(SE_lmmse, 2));
            sc = sum(mean(SE_cnsic, 2));
            bv = mean(BER_lmmse(:));
            bc = mean(BER_cnsic(:));
            sv_pf = sum(mean(SE_lmmse_pf, 2));   % [PERFECT CSI]
            sc_pf = sum(mean(SE_cnsic_pf, 2));

            if caseID == 1
                SR1_lmmse_acc(ns, si) = sv;
                SR1_cnsic_acc(ns, si) = sc;
                SR1_lmmse_pf_acc(ns, si) = sv_pf;
                SR1_cnsic_pf_acc(ns, si) = sc_pf;
                BER1_lmmse_acc(ns, si) = bv;
                BER1_cnsic_acc(ns, si) = bc;
                if any(NF_dom_users)
                    BER1_NF_acc(ns, si) = mean(mean(BER_lmmse(NF_dom_users, :)));
                    SR1_NF_acc(ns, si) = sum(mean(SE_lmmse(NF_dom_users, :), 2));
                end
                if any(collinear_mask)
                    BER1_CL_acc(ns, si) = mean(mean(BER_lmmse(collinear_mask, :)));
                    SR1_CL_acc(ns, si) = sum(mean(SE_lmmse(collinear_mask, :), 2));
                end
            else
                SR2_lmmse_acc(ns, si) = sv;
                SR2_cnsic_acc(ns, si) = sc;
                SR2_lmmse_pf_acc(ns, si) = sv_pf;
                SR2_cnsic_pf_acc(ns, si) = sc_pf;
                BER2_lmmse_acc(ns, si) = bv;
                BER2_cnsic_acc(ns, si) = bc;
                if any(NF_dom_users)
                    BER2_NF_acc(ns, si) = mean(mean(BER_lmmse(NF_dom_users, :)));
                    SR2_NF_acc(ns, si) = sum(mean(SE_lmmse(NF_dom_users, :), 2));
                end
                if any(FF_dom_users)
                    BER2_FF_acc(ns, si) = mean(mean(BER_lmmse(FF_dom_users, :)));
                    SR2_FF_acc(ns, si) = sum(mean(SE_lmmse(FF_dom_users, :), 2));
                end
                if any(collinear_mask)
                    BER2_CL_acc(ns, si) = mean(mean(BER_lmmse(collinear_mask, :)));
                    SR2_CL_acc(ns, si) = sum(mean(SE_lmmse(collinear_mask, :), 2));
                end
            end
            H = []; Hhat = []; C = [];
        end

        % --------------------------------------------------------------
        %  CASE 3: CENTRALIZED BS  [UNCHANGED]
        % --------------------------------------------------------------
        H_BS = randn(N_BS, nReal, K) + 1j * randn(N_BS, nReal, K);
        for k = 1:K
            H_BS(:, :, k) = sqrt(0.5) * Rs_BS(:, :, k) * H_BS(:, :, k);
        end

        Np_BS = sqrt(0.5) * (randn(N_BS, nReal, tau_p) + 1j * randn(N_BS, nReal, tau_p));
        Hhat_BS = zeros(N_BS, nReal, K);
        C_BS   = zeros(N_BS, N_BS, K);

        for t = 1:tau_p
            ue_t = find(pilotIndex == t)';
            yp_BS = sqrt(p) * tau_p * sum(H_BS(:, :, ue_t), 3) + sqrt(tau_p) * Np_BS(:, :, t);
            Psi_t_BS = p * tau_p * sum(R_BS(:, :, ue_t), 3) + eye(N_BS);
            for k = ue_t
                RPsi_BS = R_BS(:, :, k) / Psi_t_BS;
                Hhat_BS(:, :, k) = sqrt(p) * RPsi_BS * yp_BS;
                C_BS(:, :, k) = R_BS(:, :, k) - p * tau_p * RPsi_BS * R_BS(:, :, k);
            end
        end

        Psi_BS = sum(C_BS, 3);

        SE_lm3 = zeros(K, nReal);
        SE_sic3 = zeros(K, nReal);
        BER_lm3 = zeros(K, nReal);
        BER_sic3 = zeros(K, nReal);
        SE_lm3_pf = zeros(K, nReal);    % [PERFECT CSI]
        SE_sic3_pf = zeros(K, nReal);
        C_BS0 = zeros(size(C_BS));      % [PERFECT CSI] zero error covariance

        for mc = 1:nReal
            Hmc_BS  = reshape(H_BS(:, mc, :),   [N_BS, K]);
            Hhat_mc3 = reshape(Hhat_BS(:, mc, :), [N_BS, K]);

            for pass = 1:2
                if pass == 1                        % estimated CSI
                    He = Hhat_mc3;  PsiE = Psi_BS;  CE = C_BS;
                else                                % [PERFECT CSI] genie
                    He = Hmc_BS;    PsiE = zeros(N_BS);  CE = C_BS0;
                end
                Phi_BS = p3 * (He * He' + PsiE) + eye(N_BS);
                V3 = p3 * (Phi_BS \ He);
                se_lm_k = zeros(K, 1);
                for k = 1:K
                    oth = [1:k - 1, k + 1:K];
                    v = V3(:, k);
                    sinr = real(p3 * abs(v' * Hmc_BS(:, k))^2 / ...
                              (p3 * sum(abs(v' * Hmc_BS(:, oth)).^2) + norm(v)^2));
                    se_lm_k(k) = prelog * log2(1 + sinr);
                    if pass == 1, BER_lm3(k, mc) = 0.5 * erfc(sqrt(max(sinr, 0))); end
                end

                [~, ord3] = sort(sum(abs(He).^2, 1), 'descend');
                se_sic_k = zeros(K, 1);
                for i = 1:K
                    k_i = ord3(i);
                    act = ord3(i:end);
                    rem = ord3(i + 1:end);
                    Ha = He(:, act);
                    Ca = sum(CE(:, :, act), 3);
                    Ph = p3 * (Ha * Ha' + Ca) + eye(N_BS);
                    vs = p3 * (Ph \ He(:, k_i));
                    sig = p3 * abs(vs' * Hmc_BS(:, k_i))^2;
                    intf = 0;
                    if ~isempty(rem)
                        intf = p3 * sum(abs(vs' * Hmc_BS(:, rem)).^2);
                    end
                    sinr = real(sig / (intf + norm(vs)^2));
                    se_sic_k(k_i) = prelog * log2(1 + sinr);
                    if pass == 1, BER_sic3(k_i, mc) = 0.5 * erfc(sqrt(max(sinr, 0))); end
                end

                if pass == 1
                    SE_lm3(:, mc) = se_lm_k;   SE_sic3(:, mc) = se_sic_k;
                else
                    SE_lm3_pf(:, mc) = se_lm_k;  SE_sic3_pf(:, mc) = se_sic_k;
                end
            end
        end

        SR3_lmmse_acc(ns, si) = sum(mean(SE_lm3, 2));
        SR3_cnsic_acc(ns, si) = sum(mean(SE_sic3, 2));
        SR3_lmmse_pf_acc(ns, si) = sum(mean(SE_lm3_pf, 2));   % [PERFECT CSI]
        SR3_cnsic_pf_acc(ns, si) = sum(mean(SE_sic3_pf, 2));
        BER3_lmmse_acc(ns, si) = mean(BER_lm3(:));
        BER3_cnsic_acc(ns, si) = mean(BER_sic3(:));
        H_BS = []; Hhat_BS = []; C_BS = []; Np_BS = [];

        % --------------------------------------------------------------
        %  CLUSTERED DETECTION on the HYBRID NF/FF channel
        %  [EXTENDED: methods 3 and 4 are the new contamination rules]
        % --------------------------------------------------------------
        Hc = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hc(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                else
                    Hc(idx, :, k) = sqrt(0.5) * Rs2(:, :, l, k) * Hc(idx, :, k);
                end
            end
        end
        [Hhc, Cc] = csi_acquire(Hc, R2, NF_mask, h_det_all, beta, pilotIndex, tau_p, p, ...
                                N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);

        % ---- ORIGINAL rate score, LEFT EXACTLY AS WRITTEN (see W6) ----
        %  This omits the sum over j ~= k, so it is a monotone function of
        %  estimate energy and therefore a near duplicate of channel norm.
        %  It is retained unchanged so the old figures do not move. The
        %  NOTE: this score has NO sum over j ~= k, so it is a monotone
        %  function of estimate energy and therefore ranks APs almost
        %  exactly as channel norm does. That is why the CN-Cluster and
        %  Rate-Cluster curves nearly coincide. Mashdour eq. (2) measures
        %  the rate AFTER the local combiner and includes the interference
        %  the combiner leaves behind. Left as written, unchanged.
        SRlk = zeros(L, K);
        for l = 1:L
            for k = 1:K
                gmm = max(real(trace(R2(:, :, l, k))) - real(trace(Cc(:, :, l, k))), 0);
                er  = real(trace(Cc(:, :, l, k)));
                SRlk(l, k) = log2(1 + p * gmm / (p * er + 1));
            end
        end
        alpha_src = mean(SRlk(:));

        D_BSR = false(L, K);
        for k = 1:K
            [~, rk] = sort(SRlk(:, k), 'descend');
            sz_cn  = max(nnz(D_CN(:, k)), 1);
            sz_rt  = sz_cn;
            D_BSR(rk(1:sz_rt), k) = true;
        end
        loadCN_acc(ns, si) = mean(sum(D_CN, 1));
        loadBSR_acc(ns, si) = mean(sum(D_BSR, 1));
        clDiff_acc(ns, si) = mean(sum(D_CN ~= D_BSR, 1));


        for method = 1:2
            if method == 1
                Dcl = D_CN;
            else
                Dcl = D_BSR;
            end

            SEl = zeros(K, nReal);
            BEl = zeros(K, nReal);
            SEs = zeros(K, nReal);
            BEs = zeros(K, nReal);
            SEl_pf = zeros(K, nReal);   % [PERFECT CSI]
            SEs_pf = zeros(K, nReal);
            Cc0 = zeros(size(Cc));      % [PERFECT CSI] zero error covariance

            for mc = 1:nReal
                Hm = reshape(Hc(:, mc, :), [LN, K]);
                Hh = reshape(Hhc(:, mc, :), [LN, K]);
                [se_l, be_l] = det_local(Hh, Hm, Cc, Dcl, p, prelog, eta_FH, N, L, K);
                SEl(:, mc) = se_l;
                BEl(:, mc) = be_l;
                [se_s, be_s] = det_local_sic(Hh, Hm, Cc, Dcl, p, prelog, eta_FH, N, L, K);
                SEs(:, mc) = se_s;
                BEs(:, mc) = be_s;
                % [PERFECT CSI] true channel in the combiner, C = 0
                SEl_pf(:, mc) = det_local(Hm, Hm, Cc0, Dcl, p, prelog, eta_FH, N, L, K);
                SEs_pf(:, mc) = det_local_sic(Hm, Hm, Cc0, Dcl, p, prelog, eta_FH, N, L, K);
            end

            if method == 1
                SR_CNcl_lmmse_acc(ns, si) = sum(mean(SEl, 2));
                BER_CNcl_lmmse_acc(ns, si) = mean(BEl(:));
                SR_CNcl_cnsic_acc(ns, si) = sum(mean(SEs, 2));
                BER_CNcl_cnsic_acc(ns, si) = mean(BEs(:));
                SR_CNcl_lmmse_pf_acc(ns, si) = sum(mean(SEl_pf, 2));
                SR_CNcl_cnsic_pf_acc(ns, si) = sum(mean(SEs_pf, 2));
            else
                SR_RTcl_lmmse_acc(ns, si) = sum(mean(SEl, 2));
                BER_RTcl_lmmse_acc(ns, si) = mean(BEl(:));
                SR_RTcl_cnsic_acc(ns, si) = sum(mean(SEs, 2));
                BER_RTcl_cnsic_acc(ns, si) = mean(BEs(:));
                SR_RTcl_lmmse_pf_acc(ns, si) = sum(mean(SEl_pf, 2));
                SR_RTcl_cnsic_pf_acc(ns, si) = sum(mean(SEs_pf, 2));

                % [NEW, DI RENNA-STYLE ANALYTICAL RATE]
                % SEs above is the cluster-fused analytical sum-rate
                % (det_local_sic with the RATE cluster mask Dcl). Recompute
                % the SAME quantity with ALL APs serving every user to get
                % the cross-AP analytical rate, and store both. The gap is
                % the analytical value of the observations clustering
                % discards. Uses the existing det_local_sic unchanged.
                SEs_all = zeros(K, nReal);
                servAll = true(L, K);
                for mc = 1:nReal
                    Hm = reshape(Hc(:, mc, :),  [LN, K]);
                    Hh = reshape(Hhc(:, mc, :), [LN, K]);
                    [se_all, ~] = det_local_sic(Hh, Hm, Cc, servAll, p, prelog, eta_FH, N, L, K);
                    SEs_all(:, mc) = se_all;
                end
                SR_anaClu_acc(ns, si) = sum(mean(SEs, 2));      % cluster fusion
                SR_anaAll_acc(ns, si) = sum(mean(SEs_all, 2));  % all-AP fusion
            end
        end



        % --------------------------------------------------------------
        %  SYMBOL-LEVEL LIST DETECTION on HYBRID NF/FF
        %  [EXTENDED: the proposed and contamination rules are added]
        % --------------------------------------------------------------
        Hcs = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hcs(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                else
                    Hcs(idx, :, k) = sqrt(0.5) * Rs2(:, :, l, k) * Hcs(idx, :, k);
                end
            end
        end
        % [CSI STUDY] Hhs_est / Cs_est are the estimates BEFORE any legacy
        % near-field override (the "estimated CSI" of the paper figures).
        % AUTHENTIC mode: NF pairs by polar-domain OMP, FF pairs by LMMSE, and
        % there is no override, so Hhs == Hhs_est.
        % LEGACY mode: Hhs_est is the LMMSE with the true-channel prior R2,
        % and Hhs has NF links overwritten by the true channel.
        [Hhs, Cs, Hhs_est, Cs_est] = csi_acquire(Hcs, R2, NF_mask, h_det_all, beta, pilotSym, tau_sym, p, ...
                                                N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);

        SRs = zeros(L, K);
        for l = 1:L
            for k = 1:K
                gmm = max(real(trace(R2(:, :, l, k))) - real(trace(Cs(:, :, l, k))), 0);
                er = real(trace(Cs(:, :, l, k)));
                SRs(l, k) = log2(1 + p * gmm / (p * er + 1));
            end
        end
        D_BSR_s = false(L, K);
        for k = 1:K
            [~, rk] = sort(SRs(:, k), 'descend');
            szc = max(nnz(D_CN(:, k)), 1);
            D_BSR_s(rk(1:szc), k) = true;
        end


        servAll = true(L, K);
        bhf = zeros(4, 1);
        bcn = zeros(4, 1);
        brt = zeros(4, 1);
        nbit = 0;
        mB = zeros(4, K);
        mN = zeros(4, K);
        mbits = 0;
        meng = 0;
        etaList_run = 0;   % [NEW] measured List-SIC trigger rate, this point
        etaXap_run  = 0;   % [NEW] measured Cross-AP invocation rate, this point
        etaXpf_run  = 0;   % Cross-AP invocation rate under perfect CSI
        % [CSI STUDY] accumulators for perfect and proper-estimated CSI
        brt_est = zeros(4, 1);  brt_pf = zeros(4, 1);
        bhf_est = zeros(4, 1);  bhf_pf = zeros(4, 1);   % [CSI STUDY] full-hybrid (servAll)
        mB_est  = zeros(4, K);  mB_pf  = zeros(4, K);
        mN_est  = zeros(4, K);  mN_pf  = zeros(4, K);
        Cs0     = zeros(size(Cs));
        for mc = 1:nReal
            Hm = reshape(Hcs(:, mc, :), [LN, K]);
            Hh = reshape(Hhs(:, mc, :), [LN, K]);
            [e1, e2, e3, e4, bt] = ber_case(Hh, Hm, Cs, servAll, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
            bhf = bhf + [e1; e2; e3; e4];
            nbit = nbit + bt;
            [e1, e2, e3, e4] = ber_case(Hh, Hm, Cs, D_CN, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
            bcn = bcn + [e1; e2; e3; e4];
            [e1, e2, e3, e4, ~, etaL_mc, etaX_mc] = ber_case(Hh, Hm, Cs, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
            brt = brt + [e1; e2; e3; e4];
            etaList_run = etaList_run + etaL_mc;
            etaXap_run  = etaXap_run  + etaX_mc;
            [buu, nuu, bpu, epu] = metrics_case(Hh, Hm, Cs, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
            mB = mB + buu;
            mN = mN + nuu;
            mbits = mbits + bpu;
            meng = meng + epu;
            if csi_study
                % [CSI STUDY] proper ESTIMATED CSI: NF links genuinely estimated
                % (NUSW-LMMSE via Hhs_est/Cs_est), FF links LMMSE. Same rate cluster.
                Hhe = reshape(Hhs_est(:, mc, :), [LN, K]);
                [a1, a2, a3, a4] = ber_case(Hhe, Hm, Cs_est, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                brt_est = brt_est + [a1; a2; a3; a4];
                [bue, nue] = metrics_case(Hhe, Hm, Cs_est, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                mB_est = mB_est + bue;   mN_est = mN_est + nue;
                % [CSI STUDY] full-hybrid, estimated CSI (every AP serves every UE)
                [h1, h2, h3, h4] = ber_case(Hhe, Hm, Cs_est, servAll, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                bhf_est = bhf_est + [h1; h2; h3; h4];
                % [CSI STUDY] PERFECT CSI: true channel in the combiner, C = 0.
                [q1, q2, q3, q4, ~, ~, etaXq] = ber_case(Hm, Hm, Cs0, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                brt_pf = brt_pf + [q1; q2; q3; q4];
                etaXpf_run = etaXpf_run + etaXq;
                % [CSI STUDY] full-hybrid, perfect CSI (every AP serves every UE)
                [g1, g2, g3, g4] = ber_case(Hm, Hm, Cs0, servAll, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                bhf_pf = bhf_pf + [g1; g2; g3; g4];
                [bup, nup] = metrics_case(Hm, Hm, Cs0, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                mB_pf = mB_pf + bup;   mN_pf = mN_pf + nup;
            end
        end
        eBERhf(:, ns, si) = bhf / max(nbit, 1);
        eBERcn(:, ns, si) = bcn / max(nbit, 1);
        eBERrt(:, ns, si) = brt / max(nbit, 1);
        etaList_acc(ns, si) = etaList_run / max(nReal, 1);   % [NEW]
        etaXap_acc(ns, si)  = etaXap_run  / max(nReal, 1);   % [NEW]
        etaXap_pf_acc(ns, si) = etaXpf_run / max(nReal, 1);  % perfect-CSI gating rate
        mBERrt(:, :, ns, si)  = mB / max(mbits, 1);
        mNMSErt(:, :, ns, si) = mN / max(meng, 1);
        % [CSI STUDY] unconditional sliced writes (accumulators are zero when off)
        eBERrt_est(:, ns, si)    = brt_est / max(nbit, 1);
        eBERrt_pf(:, ns, si)     = brt_pf  / max(nbit, 1);
        eBERhf_est(:, ns, si)    = bhf_est / max(nbit, 1);
        eBERhf_pf(:, ns, si)     = bhf_pf  / max(nbit, 1);
        mBERrt_est(:, :, ns, si) = mB_est  / max(mbits, 1);
        mBERrt_pf(:, :, ns, si)  = mB_pf   / max(mbits, 1);
        mNMSErt_est(:, :, ns, si) = mN_est / max(meng, 1);
        mNMSErt_pf(:, :, ns, si)  = mN_pf  / max(meng, 1);
        % [ESTIMATOR NMSE] paper eq. (14). AUTHENTIC: empirical
        % sum||h - h_hat||^2 / sum||h||^2 (truth used only to SCORE the
        % estimator, never inside it). LEGACY: statistical tr(C)/tr(R).
        if ~genie_csi
            e_num = sum(abs(Hcs(:) - Hhs_est(:)).^2);
            e_den = sum(abs(Hcs(:)).^2);
        else
            e_num = 0;  e_den = 0;
            for l = 1:L
                for k = 1:K
                    e_num = e_num + real(trace(Cs_est(:, :, l, k)));
                    e_den = e_den + real(trace(R2(:, :, l, k)));
                end
            end
        end
        nmse_est_acc(ns, si) = e_num / max(e_den, eps);
        Hcs = []; Hhs = []; Cs = []; Hhs_est = []; Cs_est = [];

        % ================= [FAR-FIELD COMPARISON] =========================
        %  Same links, same rate cluster D_BSR_s, but ALL links use the
        %  far-field (planar) channel R1 and far-field LMMSE estimation, under
        %  perfect (C=0) and estimated CSI. Enables the NF-vs-FF figures.
        if ff_study
            Hcf = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for k = 1:K
                    Hcf(idx, :, k) = sqrt(0.5) * Rs1(:, :, l, k) * Hcf(idx, :, k);
                end
            end
            % All links far-field (statistical LMMSE in both modes).
            [Hhf, Cf] = csi_acquire(Hcf, R1, false(L, K), h_det_all, beta, pilotSym, tau_sym, p, ...
                                    N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);
            if ~genie_csi
                f_num = sum(abs(Hcf(:) - Hhf(:)).^2);
                f_den = sum(abs(Hcf(:)).^2);
            else
                f_num = 0;  f_den = 0;
                for l = 1:L
                    for k = 1:K
                        f_num = f_num + real(trace(Cf(:, :, l, k)));
                        f_den = f_den + real(trace(R1(:, :, l, k)));
                    end
                end
            end
            nmse_ff_acc(ns, si) = f_num / max(f_den, eps);
            brt_fe = zeros(4, 1);  brt_fp = zeros(4, 1);
            mB_fe = zeros(4, K);   mB_fp = zeros(4, K);
            mN_fe = zeros(4, K);   mN_fp = zeros(4, K);
            Cf0 = zeros(size(Cf));
            for mc = 1:nReal
                Hmf = reshape(Hcf(:, mc, :), [LN, K]);
                Hef = reshape(Hhf(:, mc, :), [LN, K]);
                [a1, a2, a3, a4] = ber_case(Hef, Hmf, Cf, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                brt_fe = brt_fe + [a1; a2; a3; a4];
                [bfe, nfe] = metrics_case(Hef, Hmf, Cf, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                mB_fe = mB_fe + bfe;  mN_fe = mN_fe + nfe;
                [q1, q2, q3, q4] = ber_case(Hmf, Hmf, Cf0, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                brt_fp = brt_fp + [q1; q2; q3; q4];
                [bfp, nfp] = metrics_case(Hmf, Hmf, Cf0, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                mB_fp = mB_fp + bfp;  mN_fp = mN_fp + nfp;
            end
            eBERrt_ff_est(:, ns, si)     = brt_fe / max(nbit, 1);
            eBERrt_ff_pf(:, ns, si)      = brt_fp / max(nbit, 1);
            mBERrt_ff_est(:, :, ns, si)  = mB_fe / max(mbits, 1);
            mBERrt_ff_pf(:, :, ns, si)   = mB_fp / max(mbits, 1);
            mNMSErt_ff_est(:, :, ns, si) = mN_fe / max(meng, 1);
            mNMSErt_ff_pf(:, :, ns, si)  = mN_fp / max(meng, 1);
            Hcf = []; Hhf = []; Cf = [];
        end

        %----------------------------------------------------------
        %  [FIGURE 2 EMPIRICAL] Monte-Carlo Linear + Hard-SIC BER
        %  [EXTENDED: the proposed clustering rule is added]
        %----------------------------------------------------------
        Hc2 = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hc2(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                else
                    Hc2(idx, :, k) = sqrt(0.5) * Rs2(:, :, l, k) * Hc2(idx, :, k);
                end
            end
        end
        [Hh2, Cc2] = csi_acquire(Hc2, R2, NF_mask, h_det_all, beta, pilotF2, tau_fig2, p, ...
                                 N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);
        SRs2 = zeros(L, K);
        for l = 1:L
            for k = 1:K
                gmm = max(real(trace(R2(:, :, l, k))) - real(trace(Cc2(:, :, l, k))), 0);
                er2 = real(trace(Cc2(:, :, l, k)));
                SRs2(l, k) = log2(1 + p * gmm / (p * er2 + 1));
            end
        end
        D_RT2 = false(L, K);
        for k = 1:K
            [~, rk] = sort(SRs2(:, k), 'descend');
            szc = max(nnz(D_CN(:, k)), 1);
            D_RT2(rk(1:szc), k) = true;
        end

        HcF2 = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                HcF2(idx, :, k) = sqrt(0.5) * Rs1(:, :, l, k) * HcF2(idx, :, k);
            end
        end
        [HhF2, CcF2] = csi_acquire(HcF2, R1, false(L, K), h_det_all, beta, pilotF2, tau_fig2, p, ...
                                   N, L, K, nReal, genie_csi, nf_csi_consistent, eps_NF, PD);

        Hc_BS2 = randn(N_BS, nReal, K) + 1j * randn(N_BS, nReal, K);
        for k = 1:K
            Hc_BS2(:, :, k) = sqrt(0.5) * Rs_BS(:, :, k) * Hc_BS2(:, :, k);
        end
        Np_BS2 = sqrt(0.5) * (randn(N_BS, nReal, tau_fig2) + 1j * randn(N_BS, nReal, tau_fig2));
        Hh_BS2 = zeros(N_BS, nReal, K);
        C_BS2 = zeros(N_BS, N_BS, 1, K);
        for t = 1:tau_fig2
            ue_t = find(pilotF2 == t)';
            yp = sqrt(p) * tau_fig2 * sum(Hc_BS2(:, :, ue_t), 3) + sqrt(tau_fig2) * Np_BS2(:, :, t);
            Psi_t = p * tau_fig2 * sum(R_BS(:, :, ue_t), 3) + eye(N_BS);
            for k = ue_t
                RPsi = R_BS(:, :, k) / Psi_t;
                Hh_BS2(:, :, k) = sqrt(p) * RPsi * yp;
                C_BS2(:, :, 1, k) = R_BS(:, :, k) - p * tau_fig2 * RPsi * R_BS(:, :, k);
            end
        end

        servAll = true(L, K);
        bff2 = zeros(2, 1);
        bhf2 = zeros(2, 1);
        bcn2 = zeros(2, 1);
        brt2 = zeros(2, 1);
        bbs2 = zeros(2, 1);
        bff2_pf = zeros(2, 1);  bhf2_pf = zeros(2, 1);   % [PERFECT CSI]
        bcn2_pf = zeros(2, 1);  brt2_pf = zeros(2, 1);  bbs2_pf = zeros(2, 1);
        CcF2_0 = zeros(size(CcF2));  Cc2_0 = zeros(size(Cc2));  C_BS2_0 = zeros(size(C_BS2));
        nb2 = 0;
        for mc = 1:nReal
            HmF = reshape(HcF2(:, mc, :), [LN, K]);
            HhFm = reshape(HhF2(:, mc, :), [LN, K]);
            [l1, h1, bt] = ber_linsic(HhFm, HmF, CcF2, servAll, p, N, L, K, nSym2, modOrder2);
            bff2 = bff2 + [l1; h1];
            nb2 = nb2 + bt;
            Hm = reshape(Hc2(:, mc, :), [LN, K]);
            Hh = reshape(Hh2(:, mc, :), [LN, K]);
            [l1, h1] = ber_linsic(Hh, Hm, Cc2, servAll, p, N, L, K, nSym2, modOrder2);
            bhf2 = bhf2 + [l1; h1];
            [l1, h1] = ber_linsic(Hh, Hm, Cc2, D_CN, p, N, L, K, nSym2, modOrder2);
            bcn2 = bcn2 + [l1; h1];
            [l1, h1] = ber_linsic(Hh, Hm, Cc2, D_RT2, p, N, L, K, nSym2, modOrder2);
            brt2 = brt2 + [l1; h1];
            HmB = reshape(Hc_BS2(:, mc, :), [N_BS, K]);
            HhB = reshape(Hh_BS2(:, mc, :), [N_BS, K]);
            [l1, h1] = ber_linsic(HhB, HmB, C_BS2, true(1, K), p, N_BS, 1, K, nSym2, modOrder2);
            bbs2 = bbs2 + [l1; h1];
            % [PERFECT CSI] true channel in the combiner, C = 0 (genie)
            [l1, h1] = ber_linsic(HmF, HmF, CcF2_0, servAll, p, N, L, K, nSym2, modOrder2);
            bff2_pf = bff2_pf + [l1; h1];
            [l1, h1] = ber_linsic(Hm, Hm, Cc2_0, servAll, p, N, L, K, nSym2, modOrder2);
            bhf2_pf = bhf2_pf + [l1; h1];
            [l1, h1] = ber_linsic(Hm, Hm, Cc2_0, D_CN, p, N, L, K, nSym2, modOrder2);
            bcn2_pf = bcn2_pf + [l1; h1];
            [l1, h1] = ber_linsic(Hm, Hm, Cc2_0, D_RT2, p, N, L, K, nSym2, modOrder2);
            brt2_pf = brt2_pf + [l1; h1];
            [l1, h1] = ber_linsic(HmB, HmB, C_BS2_0, true(1, K), p, N_BS, 1, K, nSym2, modOrder2);
            bbs2_pf = bbs2_pf + [l1; h1];
        end
        eBER2_ff(:, ns, si) = bff2 / max(nb2, 1);
        eBER2_hf(:, ns, si) = bhf2 / max(nb2, 1);
        eBER2_cn(:, ns, si) = bcn2 / max(nb2, 1);
        eBER2_rt(:, ns, si) = brt2 / max(nb2, 1);
        eBER2_bs(:, ns, si) = bbs2 / max(nb2, 1);
        eBER2_ff_pf(:, ns, si) = bff2_pf / max(nb2, 1);   % [PERFECT CSI]
        eBER2_hf_pf(:, ns, si) = bhf2_pf / max(nb2, 1);
        eBER2_cn_pf(:, ns, si) = bcn2_pf / max(nb2, 1);
        eBER2_rt_pf(:, ns, si) = brt2_pf / max(nb2, 1);
        eBER2_bs_pf(:, ns, si) = bbs2_pf / max(nb2, 1);
        Hc2 = []; Hh2 = []; Cc2 = []; HcF2 = []; HhF2 = []; CcF2 = []; Hc_BS2 = []; Hh_BS2 = []; C_BS2 = []; Np_BS2 = [];

    end % si

end % ns

fprintf('\nAvg NF pairs (CF): %.1f%%\n', 100 * nNF_total / (nSetups * L * K));

%% ====================================================================
%  AVERAGE OVER SETUPS
%% ====================================================================
SR1_lmmse = mean(SR1_lmmse_acc, 1);
SR1_cnsic = mean(SR1_cnsic_acc, 1);
BER1_lmmse = mean(BER1_lmmse_acc, 1);
BER1_cnsic = mean(BER1_cnsic_acc, 1);

SR2_lmmse = mean(SR2_lmmse_acc, 1);
SR2_cnsic = mean(SR2_cnsic_acc, 1);
BER2_lmmse = mean(BER2_lmmse_acc, 1);
BER2_cnsic = mean(BER2_cnsic_acc, 1);

SR3_lmmse = mean(SR3_lmmse_acc, 1);
SR3_cnsic = mean(SR3_cnsic_acc, 1);
BER3_lmmse = mean(BER3_lmmse_acc, 1);
BER3_cnsic = mean(BER3_cnsic_acc, 1);

SR_CNcl_lmmse = mean(SR_CNcl_lmmse_acc, 1);
BER_CNcl_lmmse = mean(BER_CNcl_lmmse_acc, 1);
SR_CNcl_cnsic = mean(SR_CNcl_cnsic_acc, 1);
BER_CNcl_cnsic = mean(BER_CNcl_cnsic_acc, 1);
SR_RTcl_lmmse = mean(SR_RTcl_lmmse_acc, 1);
BER_RTcl_lmmse = mean(BER_RTcl_lmmse_acc, 1);
SR_RTcl_cnsic = mean(SR_RTcl_cnsic_acc, 1);
% [NEW, DI RENNA-STYLE ANALYTICAL RATE] average over setups
SR_anaClu = mean(SR_anaClu_acc, 1);
SR_anaAll = mean(SR_anaAll_acc, 1);
SR_anaGap = SR_anaAll - SR_anaClu;   % analytical rate value of discarded APs
BER_RTcl_cnsic = mean(BER_RTcl_cnsic_acc, 1);

% [PERFECT CSI, FIG 1] genie sum-rate variants averaged over setups
SR1_lmmse_pf = mean(SR1_lmmse_pf_acc, 1);   SR1_cnsic_pf = mean(SR1_cnsic_pf_acc, 1);
SR2_lmmse_pf = mean(SR2_lmmse_pf_acc, 1);   SR2_cnsic_pf = mean(SR2_cnsic_pf_acc, 1);
SR3_lmmse_pf = mean(SR3_lmmse_pf_acc, 1);   SR3_cnsic_pf = mean(SR3_cnsic_pf_acc, 1);
SR_CNcl_lmmse_pf = mean(SR_CNcl_lmmse_pf_acc, 1);  SR_CNcl_cnsic_pf = mean(SR_CNcl_cnsic_pf_acc, 1);
SR_RTcl_lmmse_pf = mean(SR_RTcl_lmmse_pf_acc, 1);  SR_RTcl_cnsic_pf = mean(SR_RTcl_cnsic_pf_acc, 1);

BER1_NF = mean(BER1_NF_acc, 1);
SR1_NF = mean(SR1_NF_acc, 1);
BER2_NF = mean(BER2_NF_acc, 1);
SR2_NF = mean(SR2_NF_acc, 1);
BER2_FF = mean(BER2_FF_acc, 1);
SR2_FF = mean(SR2_FF_acc, 1);
BER1_CL = mean(BER1_CL_acc, 1);
SR1_CL = mean(SR1_CL_acc, 1);
BER2_CL = mean(BER2_CL_acc, 1);
SR2_CL = mean(SR2_CL_acc, 1);

if nSetups == 1
    aBERhf = reshape(eBERhf, 4, nSNR);
    aBERcn = reshape(eBERcn, 4, nSNR);
    aBERrt = reshape(eBERrt, 4, nSNR);
    aBER2_ff = reshape(eBER2_ff, 2, nSNR);
    aBER2_hf = reshape(eBER2_hf, 2, nSNR);
    aBER2_bs = reshape(eBER2_bs, 2, nSNR);
    aBER2_cn = reshape(eBER2_cn, 2, nSNR);
    aBER2_rt = reshape(eBER2_rt, 2, nSNR);
    aBER2_ff_pf = reshape(eBER2_ff_pf, 2, nSNR);   % [PERFECT CSI]
    aBER2_hf_pf = reshape(eBER2_hf_pf, 2, nSNR);
    aBER2_bs_pf = reshape(eBER2_bs_pf, 2, nSNR);
    aBER2_cn_pf = reshape(eBER2_cn_pf, 2, nSNR);
    aBER2_rt_pf = reshape(eBER2_rt_pf, 2, nSNR);
    mBER_u  = reshape(mBERrt, 4, K, nSNR);
    mNMSE_u = reshape(mNMSErt, 4, K, nSNR);
    aBERrt_est = reshape(eBERrt_est, 4, nSNR);      % [CSI STUDY]
    aBERrt_pf  = reshape(eBERrt_pf, 4, nSNR);
    aBERhf_est = reshape(eBERhf_est, 4, nSNR);
    aBERhf_pf  = reshape(eBERhf_pf, 4, nSNR);
    mBER_u_est = reshape(mBERrt_est, 4, K, nSNR);
    mBER_u_pf  = reshape(mBERrt_pf, 4, K, nSNR);
    mNMSE_u_est = reshape(mNMSErt_est, 4, K, nSNR);
    mNMSE_u_pf  = reshape(mNMSErt_pf, 4, K, nSNR);
    aBERrt_ff_est = reshape(eBERrt_ff_est, 4, nSNR);       % [FAR-FIELD]
    aBERrt_ff_pf  = reshape(eBERrt_ff_pf, 4, nSNR);
    mBER_u_ff_est = reshape(mBERrt_ff_est, 4, K, nSNR);
    mBER_u_ff_pf  = reshape(mBERrt_ff_pf, 4, K, nSNR);
    mNMSE_u_ff_est = reshape(mNMSErt_ff_est, 4, K, nSNR);
    mNMSE_u_ff_pf  = reshape(mNMSErt_ff_pf, 4, K, nSNR);
else
    aBERhf = squeeze(mean(eBERhf, 2));
    aBERcn = squeeze(mean(eBERcn, 2));
    aBERrt = squeeze(mean(eBERrt, 2));
    aBER2_ff = squeeze(mean(eBER2_ff, 2));
    aBER2_hf = squeeze(mean(eBER2_hf, 2));
    aBER2_bs = squeeze(mean(eBER2_bs, 2));
    aBER2_cn = squeeze(mean(eBER2_cn, 2));
    aBER2_rt = squeeze(mean(eBER2_rt, 2));
    aBER2_ff_pf = squeeze(mean(eBER2_ff_pf, 2));   % [PERFECT CSI]
    aBER2_hf_pf = squeeze(mean(eBER2_hf_pf, 2));
    aBER2_bs_pf = squeeze(mean(eBER2_bs_pf, 2));
    aBER2_cn_pf = squeeze(mean(eBER2_cn_pf, 2));
    aBER2_rt_pf = squeeze(mean(eBER2_rt_pf, 2));
    mBER_u  = squeeze(mean(mBERrt, 3));
    mNMSE_u = squeeze(mean(mNMSErt, 3));
    aBERrt_est = squeeze(mean(eBERrt_est, 2));      % [CSI STUDY]
    aBERrt_pf  = squeeze(mean(eBERrt_pf, 2));
    aBERhf_est = squeeze(mean(eBERhf_est, 2));
    aBERhf_pf  = squeeze(mean(eBERhf_pf, 2));
    mBER_u_est = squeeze(mean(mBERrt_est, 3));
    mBER_u_pf  = squeeze(mean(mBERrt_pf, 3));
    mNMSE_u_est = squeeze(mean(mNMSErt_est, 3));
    mNMSE_u_pf  = squeeze(mean(mNMSErt_pf, 3));
    aBERrt_ff_est = squeeze(mean(eBERrt_ff_est, 2));       % [FAR-FIELD]
    aBERrt_ff_pf  = squeeze(mean(eBERrt_ff_pf, 2));
    mBER_u_ff_est = squeeze(mean(mBERrt_ff_est, 3));
    mBER_u_ff_pf  = squeeze(mean(mBERrt_ff_pf, 3));
    mNMSE_u_ff_est = squeeze(mean(mNMSErt_ff_est, 3));
    mNMSE_u_ff_pf  = squeeze(mean(mNMSErt_ff_pf, 3));
end
nmse_est = mean(nmse_est_acc, 1);    % [ESTIMATOR NMSE] 1 x nSNR
nmse_ff  = mean(nmse_ff_acc, 1);



bits_sym = log2(modOrder);
SE_det   = zeros(4, nSNR);
NMSE_det = zeros(4, nSNR);
for si = 1:nSNR
    for d = 1:4
        gk = prelog * bits_sym * (1 - squeeze(mBER_u(d, :, si)));
        SE_det(d, si)   = sum(gk);
        NMSE_det(d, si) = mean(squeeze(mNMSE_u(d, :, si)));
    end
end

% [CSI STUDY] goodput sum-rate of the 4 detectors under estimated and perfect CSI
SE_det_est = zeros(4, nSNR);   SE_det_pf = zeros(4, nSNR);
SE_det_ff_est = zeros(4, nSNR);  SE_det_ff_pf = zeros(4, nSNR);   % [FAR-FIELD] goodput
NMSE_det_est = zeros(4, nSNR);  NMSE_det_pf = zeros(4, nSNR);     % [FIG 7] symbol NMSE, est/pf (hybrid)
for si = 1:nSNR
    for d = 1:4
        SE_det_est(d, si) = sum(prelog * bits_sym * (1 - squeeze(mBER_u_est(d, :, si))));
        SE_det_pf(d, si)  = sum(prelog * bits_sym * (1 - squeeze(mBER_u_pf(d, :, si))));
        SE_det_ff_est(d, si) = sum(prelog * bits_sym * (1 - squeeze(mBER_u_ff_est(d, :, si))));
        SE_det_ff_pf(d, si)  = sum(prelog * bits_sym * (1 - squeeze(mBER_u_ff_pf(d, :, si))));
        NMSE_det_est(d, si) = mean(squeeze(mNMSE_u_est(d, :, si)));
        NMSE_det_pf(d, si)  = mean(squeeze(mNMSE_u_pf(d, :, si)));
    end
end

% [NEW, FIGURE 9] Sum-rate of all four DETECTION techniques.
%
%  WHY THIS IS NOT THE SHANNON SUM-RATE OF FIGURES 1 AND 4.
%  Figures 1 and 4 evaluate prelog*log2(1+SINR_k) with SINR_k from the
%  combiner. List-SIC and Cross-AP List-SIC use the SAME combiners and
%  therefore the SAME SINR as hard SIC: they differ only in which symbol
%  is DECIDED. A closed-form Shannon curve is identical for all three by
%  construction and cannot separate the detectors. Any figure claiming
%  otherwise would be wrong.
%
%  Two honest detector-resolving rates are plotted instead:
%   (a) GOODPUT sum-rate, sum_k prelog*log2(M)*(1-BER_k). Bounded above by
%       K*prelog*log2(M) because the modulation is fixed at 16-QAM. This
%       is the same quantity as Figure 6 and is repeated here so the two
%       rate definitions sit side by side.
%   (b) EFFECTIVE-SINR sum-rate, sum_k prelog*log2(1+1/NMSE_k), where
%       NMSE_k is the measured hard-decision symbol NMSE. This is the
%       EVM-style effective SINR and it is unbounded, so it looks like a
%       conventional sum-rate curve.
%
%  FLOOR ON (b), READ BEFORE QUOTING ANY NUMBER FROM IT. NMSE is measured
%  from a finite number of symbols, so once a detector makes no symbol
%  errors the measured NMSE is 0 and log2(1+1/0) is infinite. NMSE is
%  therefore floored at 1/nSymTot, the smallest value the experiment can
%  resolve. Curves that flatten at the top of the right panel have hit
%  that floor and are NOT a real saturation: they mean the detector made
%  too few errors to measure. Raise nSym to push the floor down.
nSymTot   = nSetups * nReal * nSym;
nmseFloor = 1 / nSymTot;
SEeff_det = zeros(4, nSNR);
for si = 1:nSNR
    for d = 1:4
        nk = max(squeeze(mNMSE_u(d, :, si)), nmseFloor);
        SEeff_det(d, si) = sum(prelog * log2(1 + 1 ./ nk));
    end
end
SEeff_capped = any(squeeze(mNMSE_u(:, :, :)) < nmseFloor, 2);

% [ACHIEVABLE RATE] effective-SINR sum-rate SE_k = prelog*log2(1+1/NMSE_k). Unlike
% the goodput (which saturates at the 16-QAM alphabet bound, so all detectors meet
% at high SNR), this is unbounded and keeps growing with SNR, separating the
% detectors by their decision quality. Computed for perfect and estimated CSI,
% near-field/hybrid and far-field.
SEeff_det_est = zeros(4, nSNR);  SEeff_det_pf = zeros(4, nSNR);
SEeff_det_ff_est = zeros(4, nSNR);  SEeff_det_ff_pf = zeros(4, nSNR);
for si = 1:nSNR
    for d = 1:4
        SEeff_det_est(d, si)    = sum(prelog * log2(1 + 1 ./ max(squeeze(mNMSE_u_est(d, :, si)),    nmseFloor)));
        SEeff_det_pf(d, si)     = sum(prelog * log2(1 + 1 ./ max(squeeze(mNMSE_u_pf(d, :, si)),     nmseFloor)));
        SEeff_det_ff_est(d, si) = sum(prelog * log2(1 + 1 ./ max(squeeze(mNMSE_u_ff_est(d, :, si)), nmseFloor)));
        SEeff_det_ff_pf(d, si)  = sum(prelog * log2(1 + 1 ./ max(squeeze(mNMSE_u_ff_pf(d, :, si)),  nmseFloor)));
    end
end

si_cdf = snr_ref_idx;
seSamp = cell(4, 1);
seSamp_est = cell(4, 1);   seSamp_pf = cell(4, 1);   % [FIG 8] per-user SE samples, est/pf (hybrid)
for d = 1:4
    berk = squeeze(mBERrt(d, :, :, si_cdf));
    seSamp{d} = prelog * bits_sym * (1 - berk(:));
    berk_e = squeeze(mBERrt_est(d, :, :, si_cdf));
    seSamp_est{d} = prelog * bits_sym * (1 - berk_e(:));
    berk_p = squeeze(mBERrt_pf(d, :, :, si_cdf));
    seSamp_pf{d} = prelog * bits_sym * (1 - berk_p(:));
end

%% ====================================================================
%  [CODE C] SPARSE-ARRAY COMPARISON STUDY
%  ------------------------------------------------------------------
%  Same AP/UE drops, same large-scale fading (beta and LoS draws), same
%  orthogonal pilots and the same four detectors for every geometry in
%  arr_list; ONLY the array mounted at the APs changes. Per array and setup:
%   (a) geometry : NF fraction (d < Rayleigh distance), focusing fraction
%       (d < R cos^2(theta)/10.24, Zhou Prop. 2: beyond it no array can focus
%       in range), per-AP pair correlation mu, probability that a pair is
%       unresolvable (mu > mu_th) at the home AP, in the serving cluster, at
%       the outer (non-serving) APs and over all APs, the rescue probability
%       (unresolvable in the cluster, resolvable at the outer APs), and the
%       per-AP / network ambiguity areas;
%   (b) link     : authentic CSI (polar-domain OMP on NF links, statistical
%       LMMSE on FF links), empirical NMSE;
%   (c) detection: BER of Linear, Hard-SIC, List-SIC and List+CrossAP with
%       estimated and perfect CSI, on the hybrid model and on the all-far-
%       field model; cross-AP gate firing rate; serving-cluster size;
%   (d) at the reference SNR: captured information kappa_k, information
%       surplus Delta_k (VTC eq. (26)) and per-user effective-SINR SE.
%  Cluster-level and network correlations are NON-COHERENT across APs,
%  sum_l |h_lk^H h_lj| / sqrt(sum_l ||h_lk||^2 sum_l ||h_lj||^2), which is the
%  relevant measure for local combining + LSFD (level-3 processing): each AP
%  separates users with its own array, the CPU only weights the outputs.
%% ====================================================================
arr_study    = true;
arr_list     = {'ula', 'usa', 'eca', 'eca_sa'};
arr_names    = {'Dense ULA (64, \lambda/2)', 'Uniform sparse (64, \lambda)', ...
                'ECA non-uniform (63)', 'ECA same aperture (31)'};
nArr         = numel(arr_list);
nSetups_arr  = 20;
nReal_arr    = 4;
nSym_arr     = 200;
mu_th        = 0.7;    % pair unresolvable when its correlation exceeds mu_th: separating
                       % two users with correlation mu costs a factor 1 - mu^2 in SINR
                       % (zero-forcing view), so mu > 1/sqrt(2) ~ 0.7 is a > 3 dB loss
amb_th       = 0.5;    % ambiguity point: normalised beam gain above amb_th
amb_grid     = 2.5;    % [m] grid step of the ambiguity-area computation
amb_excl     = 5;      % [m] exclusion radius around the target (network main lobe)
lvl_names    = {'home AP', 'serving cluster', 'outer APs', 'all APs', 'rescue (clu. unres., outer res.)'};

A_geo   = cell(1, nArr);                  % geometry of every array (for figures)
PDs     = cell(1, nArr);                  % polar dictionary of every array
for a = 1:nArr
    A_geo{a} = array_geometry(arr_list{a}, N_ref, lambda, usa_eta, eca_M, eca_N, eca_Lc, eca_Lc_sa);
    PDs{a} = build_polar_dict(A_geo{a}.pos, lambda, pd_os_ang, pd_q_max, pd_n_q, pd_phase_step, pd_nS_min);
    PDs{a}.smax = pd_smax;  PDs{a}.stop = pd_stop;  PDs{a}.nref = pd_nref;  PDs{a}.nstart = pd_nstart;
end
% one UE-placement radius for every array (that of the dense reference), so
% the drops are identical and only the array differs
r_max_arr = min(0.9 * 2 * ((N_ref - 1) * d_ant)^2 / lambda, squareLen / 2);

A_nf    = zeros(nArr, nSetups_arr);          % NF link fraction
A_foc   = zeros(nArr, nSetups_arr);          % focusing-capable link fraction
A_pu    = zeros(5, 2, nArr, nSetups_arr);    % P(unresolvable): rows lvl_names, cols [all pairs, collinear pairs]
A_muA   = cell(nArr, nSetups_arr);           % home-AP correlation of every ordered pair
A_muN   = cell(nArr, nSetups_arr);           % network (all-AP) correlation of every ordered pair
A_amb   = zeros(3, nArr, nSetups_arr);       % ambiguity area fraction: [per-AP range, per-AP grating/sidelobe, network]
A_nmse  = zeros(nArr, nSNR, nSetups_arr);
A_ber_e = zeros(4, nArr, nSNR, nSetups_arr); % hybrid model, estimated CSI
A_ber_p = zeros(4, nArr, nSNR, nSetups_arr); % hybrid model, perfect CSI
A_bff_e = zeros(4, nArr, nSNR, nSetups_arr); % all-far-field model, estimated CSI
A_bff_p = zeros(4, nArr, nSNR, nSetups_arr); % all-far-field model, perfect CSI
A_eta   = zeros(nArr, nSNR, nSetups_arr);    % cross-AP gate firing rate, estimated CSI
A_size  = zeros(nArr, nSNR, nSetups_arr);    % mean serving-cluster size
A_kap   = zeros(nArr, K, nSetups_arr);       % captured information, reference SNR
A_del   = zeros(nArr, K, nSetups_arr);       % information surplus, reference SNR
A_se    = zeros(4, nArr, K, nSetups_arr);    % per-user effective-SINR SE, reference SNR
A_time  = zeros(nArr, nSetups_arr);          % run time per array and setup [s]

if arr_study
fprintf('\n=== [CODE C] Sparse-array comparison: %d setups x %d MC, arrays: %s ===\n', ...
        nSetups_arr, nReal_arr, strjoin(arr_list, ', '));
for a = 1:nArr
    fprintf('  %-8s N=%3d  aperture=%5.2f m  R=%7.1f m  focus limit (broadside)=%6.1f m  dictionary=%d atoms\n', ...
            arr_list{a}, A_geo{a}.N, A_geo{a}.aperture, 2 * A_geo{a}.aperture^2 / lambda, ...
            2 * A_geo{a}.aperture^2 / lambda / 10.24, size(PDs{a}.W, 2));
end
parfor ns = 1:nSetups_arr
    % ---------------- common drop, identical for every array ----------------
    APp = (rand(L, 1) + 1j * rand(L, 1)) * squareLen;
    [UEp, ~, homeAP] = drop_users(APp, K, L, n_collinear_pairs, r_min_NF, r_max_arr);
    APw = wrap_images(APp, squareLen, wrap_fix);
    d3 = zeros(L, K);  sE = zeros(L, K);  phL = zeros(L, K);  thL = zeros(L, K);  bet = zeros(L, K);
    for l = 1:L
        for k = 1:K
            [dH, wI] = min(abs(APw(l, :) - UEp(k)));
            d3(l, k)  = max(sqrt(hDiff^2 + dH^2), 10);
            phL(l, k) = angle(UEp(k) - APw(l, wI));
            thL(l, k) = asin(hDiff / d3(l, k));
            sE(l, k)  = sin(phL(l, k)) * cos(thL(l, k));
            bet(l, k) = large_scale_gain(d3(l, k), max(dH, 1), fc_GHz, hBS, hUT, c0, noisePow_dBm);
        end
    end
    isPair = false(K, K);                         % members of the same collinear pair
    for pr = 1:n_collinear_pairs
        isPair(2 * pr - 1, 2 * pr) = true;
        isPair(2 * pr, 2 * pr - 1) = true;
    end

    t_nf = zeros(nArr, 1);   t_foc = zeros(nArr, 1);   t_pu = zeros(5, 2, nArr);
    t_muA = cell(nArr, 1);   t_muN = cell(nArr, 1);    t_amb = zeros(3, nArr);
    t_nmse = zeros(nArr, nSNR);  t_be = zeros(4, nArr, nSNR);  t_bp = zeros(4, nArr, nSNR);
    t_fe = zeros(4, nArr, nSNR); t_fp = zeros(4, nArr, nSNR);  t_eta = zeros(nArr, nSNR);
    t_sz = zeros(nArr, nSNR);    t_kap = zeros(nArr, K);       t_del = zeros(nArr, K);
    t_se = zeros(4, nArr, K);    t_tm = zeros(nArr, 1);

    for a = 1:nArr
        tA = tic;
        G   = A_geo{a};
        Na  = G.N;
        LNa = L * Na;
        Ra  = 2 * G.aperture^2 / lambda;
        % ---------------- channel statistics for this array ----------------
        hD  = zeros(Na, L, K);
        R1a = zeros(Na, Na, L, K);
        R2a = zeros(Na, Na, L, K);
        Rs1 = zeros(Na, Na, L, K);
        Rs2 = zeros(Na, Na, L, K);
        U   = zeros(Na, L, K);                    % unit-norm link directions (hybrid model)
        NFa = d3 < Ra;
        for l = 1:L
            for k = 1:K
                R1a(:, :, l, k) = far_corr(G.pos / lambda, bet(l, k), sE(l, k), phL(l, k), thL(l, k), ASD_varphi, ASD_theta);
                Rs1(:, :, l, k) = sqrtm(R1a(:, :, l, k));
                if NFa(l, k)
                    rn = sqrt(d3(l, k)^2 + G.pos.^2 - 2 * d3(l, k) * G.pos * sE(l, k));
                    hD(:, l, k) = sqrt(bet(l, k)) * (d3(l, k) ./ rn) .* exp(-1j * 2 * pi * rn / lambda);
                    R2a(:, :, l, k) = hD(:, l, k) * hD(:, l, k)' + eps_NF * bet(l, k) * eye(Na);
                    U(:, l, k) = hD(:, l, k) / norm(hD(:, l, k));
                else
                    R2a(:, :, l, k) = R1a(:, :, l, k);
                    Rs2(:, :, l, k) = Rs1(:, :, l, k);
                    U(:, l, k) = exp(1j * 2 * pi * G.pos * sE(l, k) / lambda) / sqrt(Na);
                end
            end
        end
        t_nf(a)  = mean(NFa(:));
        t_foc(a) = mean(d3(:) < Ra * (1 - sE(:).^2) / 10.24);
        mCN = bet;
        mCN(NFa) = Na * bet(NFa);                 % ||h_det||^2 ~ N beta (as metric_CN in the main code)
        Dcn = cluster_cn_rule(mCN, bet, NFa, clu_div);

        % ---------------- (a) pair resolvability ----------------
        MU = zeros(L, K, K);
        for l = 1:L
            Ul = reshape(U(:, l, :), [Na, K]);
            MU(l, :, :) = reshape(abs(Ul' * Ul), [1, K, K]);
        end
        E = Na * bet;                             % link energies
        cntU = zeros(5, 2);  cntP = zeros(1, 2);
        muA = zeros(K * (K - 1), 1);  muN = zeros(K * (K - 1), 1);  q = 0;
        for k = 1:K
            sets = {homeAP(k), find(Dcn(:, k)), find(~Dcn(:, k)), (1:L)'};
            for j = [1:k - 1, k + 1:K]
                muS = ones(4, 1);
                for v = 1:4
                    S = sets{v};
                    if ~isempty(S)
                        muS(v) = sum(sqrt(E(S, k) .* E(S, j)) .* MU(S, k, j)) / sqrt(sum(E(S, k)) * sum(E(S, j)));
                    end
                end
                un = muS > mu_th;
                un5 = [un; un(2) && ~un(3)];
                q = q + 1;
                muA(q) = muS(1);
                muN(q) = muS(4);
                cntU(:, 1) = cntU(:, 1) + un5;   cntP(1) = cntP(1) + 1;
                if isPair(k, j)
                    cntU(:, 2) = cntU(:, 2) + un5;   cntP(2) = cntP(2) + 1;
                end
            end
        end
        t_pu(:, :, a) = cntU ./ max(cntP, 1);
        t_muA{a} = muA;
        t_muN{a} = muN;
        [aR, aG, aN] = ambiguity_areas(G.pos, lambda, APp, UEp, squareLen, hDiff, amb_grid, amb_th, amb_excl, wrap_fix);
        t_amb(:, a) = [aR; aG; aN];

        % ---------------- (b)-(d) estimation and detection ----------------
        for si = 1:nSNR
            p = SNR_lin(si);
            % hybrid model
            Ht = draw_channel(hD, Rs2, NFa, Na, L, K, nReal_arr);
            [Hh, Cs] = csi_acquire(Ht, R2a, NFa, hD, bet, (1:K)', K, p, ...
                                   Na, L, K, nReal_arr, false, false, eps_NF, PDs{a});
            t_nmse(a, si) = sum(abs(Ht(:) - Hh(:)).^2) / sum(abs(Ht(:)).^2);
            Drt = cluster_rate_rule(R2a, Cs, p, Dcn);
            t_sz(a, si) = mean(sum(Drt, 1));
            Cs0 = zeros(size(Cs));
            bE = zeros(4, K);  nE = zeros(4, K);  bP = zeros(4, K);  nbit = 0;  neng = 0;  eX = 0;
            for mc = 1:nReal_arr
                Hm = reshape(Ht(:, mc, :), [LNa, K]);
                He = reshape(Hh(:, mc, :), [LNa, K]);
                [bu, nu, bpu, epu, ex] = metrics_case(He, Hm, Cs, Drt, p, Na, L, K, nSym_arr, M_lst, d_th, maxBr, modOrder);
                bE = bE + bu;   nE = nE + nu;   nbit = nbit + bpu;   neng = neng + epu;   eX = eX + ex;
                bu = metrics_case(Hm, Hm, Cs0, Drt, p, Na, L, K, nSym_arr, M_lst, d_th, maxBr, modOrder);
                bP = bP + bu;
            end
            t_be(:, a, si) = sum(bE, 2) / max(K * nbit, 1);
            t_bp(:, a, si) = sum(bP, 2) / max(K * nbit, 1);
            t_eta(a, si)   = eX / nReal_arr;
            if si == snr_ref_idx
                g = per_ap_sinr_eval(Hh, Ht, p, Na, L, K, nReal_arr);
                [kp, dl] = surplus_metrics(g, Drt);
                t_kap(a, :) = kp;
                t_del(a, :) = dl;
                t_se(:, a, :) = reshape(prelog * log2(1 + 1 ./ max(nE / max(neng, 1), 1 / (nReal_arr * nSym_arr))), [4, 1, K]);
            end
            % all-far-field model (same cluster rule, LMMSE on every link)
            Hf = draw_channel(hD, Rs1, false(L, K), Na, L, K, nReal_arr);
            [Hfh, Cf] = csi_acquire(Hf, R1a, false(L, K), hD, bet, (1:K)', K, p, ...
                                    Na, L, K, nReal_arr, false, false, eps_NF, PDs{a});
            Df = cluster_rate_rule(R1a, Cf, p, cluster_cn_rule(bet, bet, false(L, K), clu_div));
            Cf0 = zeros(size(Cf));
            bF = zeros(4, K);  bFp = zeros(4, K);  nbf = 0;
            for mc = 1:nReal_arr
                Hm = reshape(Hf(:, mc, :), [LNa, K]);
                He = reshape(Hfh(:, mc, :), [LNa, K]);
                [bu, ~, bpu] = metrics_case(He, Hm, Cf, Df, p, Na, L, K, nSym_arr, M_lst, d_th, maxBr, modOrder);
                bF = bF + bu;  nbf = nbf + bpu;
                bu = metrics_case(Hm, Hm, Cf0, Df, p, Na, L, K, nSym_arr, M_lst, d_th, maxBr, modOrder);
                bFp = bFp + bu;
            end
            t_fe(:, a, si) = sum(bF, 2)  / max(K * nbf, 1);
            t_fp(:, a, si) = sum(bFp, 2) / max(K * nbf, 1);
        end
        t_tm(a) = toc(tA);
    end
    A_nf(:, ns) = t_nf;          A_foc(:, ns) = t_foc;       A_pu(:, :, :, ns) = t_pu;
    A_muA(:, ns) = t_muA;        A_muN(:, ns) = t_muN;       A_amb(:, :, ns) = t_amb;
    A_nmse(:, :, ns) = t_nmse;   A_ber_e(:, :, :, ns) = t_be; A_ber_p(:, :, :, ns) = t_bp;
    A_bff_e(:, :, :, ns) = t_fe; A_bff_p(:, :, :, ns) = t_fp; A_eta(:, :, ns) = t_eta;
    A_size(:, :, ns) = t_sz;     A_kap(:, :, ns) = t_kap;    A_del(:, :, ns) = t_del;
    A_se(:, :, :, ns) = t_se;    A_time(:, ns) = t_tm;
    fprintf('  array study setup %d/%d done\n', ns, nSetups_arr);
end
end

%% ====================================================================
%  DEPLOYMENT AREA SWEEP  [UNCHANGED]
%% ====================================================================
fprintf('--- Area Sweep: %s m (%d setups/pt, SNR=10dB) ---\n', ...
        mat2str(squareLen_sweep), nSetups_sweep);

NF_frac_sw = zeros(1, nSweep);
BER1_sw = zeros(1, nSweep);
BER2_sw = zeros(1, nSweep);
SR1_sw = zeros(1, nSweep);
SR2_sw = zeros(1, nSweep);

for sw = 1:nSweep
    sLen = squareLen_sweep(sw);
    rmax_sw = min(0.9 * d_Ray, sLen / 2);
    nNF_sw = 0;
    b1_sw = 0;
    b2_sw = 0;
    s1_sw = 0;
    s2_sw = 0;

    for ns = 1:nSetups_sweep
        APsw = (rand(L, 1) + 1j * rand(L, 1)) * sLen;
        UEsw = zeros(K, 1);
        ui = 1;
        for pr = 1:n_collinear_pairs
            al = mod(pr - 1, L) + 1;
            phi = 2 * pi * rand;
            rn = r_min_NF + (rmax_sw / 2 - r_min_NF) * rand;
            rf = rmax_sw / 2 + (rmax_sw - rmax_sw / 2) * rand;
            UEsw(ui) = APsw(al) + rn * exp(1j * phi);
            UEsw(ui + 1) = APsw(al) + rf * exp(1j * phi);
            ui = ui + 2;
        end
        ar = 1;
        while ui <= K
            r = r_min_NF + (rmax_sw - r_min_NF) * rand;
            phi = 2 * pi * rand;
            UEsw(ui) = APsw(ar) + r * exp(1j * phi);
            ui = ui + 1;
            ar = mod(ar, L) + 1;
        end

        APw_sw = wrap_images(APsw, sLen, wrap_fix);   % [CODE C]
        beta_sw = zeros(L, K);
        R1_sw = zeros(N, N, L, K);
        R2_sw = zeros(N, N, L, K);
        hd_sw = zeros(N, L, K, 'like', 1 + 1j);
        NF_sw = false(L, K);
        nNFs = 0;

        for l = 1:L
            for k = 1:K
                [dH, ~] = min(abs(APw_sw(l, :) - UEsw(k)));
                d3D = max(sqrt(hDiff^2 + dH^2), 10);
                d2D = max(dH, 1);
                [~, wi] = min(abs(APw_sw(l, :) - UEsw(k)));
                ph_lk = angle(UEsw(k) - APw_sw(l, wi));
                th = asin(hDiff / d3D);
                se = sin(ph_lk) * cos(th);
                if d2D <= 18
                    PL = 1;
                else
                    PL = (18 / d2D) * (1 - exp(-d2D / 63)) + exp(-d2D / 63);
                end
                dBP_sw = 4 * (hBS - 1) * (hUT - 1) * fc_GHz * 1e9 / c0;
                if d3D < dBP_sw
                    PLs = 28 + 37 * log10(d3D) + 20 * log10(fc_GHz);
                else
                    PLs = 28 + 20 * log10(d3D) + 20 * log10(fc_GHz) - 9 * log10(dBP_sw^2 + (hBS - hUT)^2);
                end
                if rand < PL
                    PL_dB = PLs + 4 * randn;
                else
                    PLn = 13.54 + 39.08 * log10(d3D) + 20 * log10(fc_GHz) - 0.6 * (hUT - 1.5);
                    PL_dB = max(PLs, PLn) + 6 * randn;
                end
                beta_sw(l, k) = db2pow(-PL_dB - noisePow_dBm);
                RFF_lk = far_corr(APG.pos / lambda, beta_sw(l, k), se, ph_lk, th, ASD_varphi, ASD_theta);   % [CODE C]
                R1_sw(:, :, l, k) = RFF_lk;
                if d3D < d_Ray
                    NF_sw(l, k) = true;
                    nNFs = nNFs + 1;
                    rn_s = sqrt(d3D^2 + APG.pos.^2 - 2 * d3D * APG.pos * se);   % [CODE C]
                    hd = sqrt(beta_sw(l, k)) * (d3D ./ rn_s) .* exp(-1j * 2 * pi * rn_s / lambda);
                    hd_sw(:, l, k) = hd;
                    R2_sw(:, :, l, k) = hd * hd' + eps_NF * beta_sw(l, k) * eye(N);
                else
                    R2_sw(:, :, l, k) = RFF_lk;
                end
            end
        end
        nNF_sw = nNF_sw + nNFs / (L * K);
        p_ref = SNR_lin(snr_ref_idx);
        piIdx = mod((0:K - 1)', tau_p) + 1;

        for cID = 1:2
            if cID == 1
                Ru = R1_sw;
            else
                Ru = R2_sw;
            end
            H_sw = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for k = 1:K
                    if cID == 2 && NF_sw(l, k)
                        H_sw(idx, :, k) = repmat(hd_sw(:, l, k), [1, nReal]);
                    else
                        H_sw(idx, :, k) = sqrt(0.5) * sqrtm(Ru(:, :, l, k)) * H_sw(idx, :, k);
                    end
                end
            end
            if cID == 1
                NFm_sw = false(L, K);
            else
                NFm_sw = NF_sw;
            end
            % Legacy area sweep overwrote NF estimates but kept the LMMSE C,
            % hence nf-consistent flag = false here (legacy reproduction).
            [Hh, Csw] = csi_acquire(H_sw, Ru, NFm_sw, hd_sw, beta_sw, piIdx, tau_p, p_ref, ...
                                    N, L, K, nReal, genie_csi, false, eps_NF, PD);
            Pl = zeros(N, N, L);
            for l = 1:L
                Pl(:, :, l) = sum(Csw(:, :, l, :), 4);
            end
            SEsw = zeros(K, nReal);
            BERsw = zeros(K, nReal);
            for mc = 1:nReal
                Hm = reshape(H_sw(:, mc, :), [LN, K]);
                Hm2 = reshape(Hh(:, mc, :), [LN, K]);
                V = zeros(LN, K);
                for l = 1:L
                    idx = (l - 1) * N + 1:l * N;
                    Hl = Hm2(idx, :);
                    Ph = p_ref * (Hl * Hl' + Pl(:, :, l)) + eye(N);
                    V(idx, :) = p_ref * (Ph \ Hl);
                end
                for k = 1:K
                    oth = [1:k - 1, k + 1:K];
                    v = V(:, k);
                    sinr = real(p_ref * abs(v' * Hm(:, k))^2 / (p_ref * sum(abs(v' * Hm(:, oth)).^2) + norm(v)^2));
                    SEsw(k, mc) = prelog * log2(1 + sinr);
                    BERsw(k, mc) = 0.5 * erfc(sqrt(max(sinr, 0)));
                end
            end
            if cID == 1
                s1_sw = s1_sw + sum(mean(SEsw, 2));
                b1_sw = b1_sw + mean(BERsw(:));
            else
                s2_sw = s2_sw + sum(mean(SEsw, 2));
                b2_sw = b2_sw + mean(BERsw(:));
            end
            clear H_sw Hh Csw;
        end
    end
    NF_frac_sw(sw) = 100 * nNF_sw / nSetups_sweep;
    BER1_sw(sw) = b1_sw / nSetups_sweep;
    BER2_sw(sw) = b2_sw / nSetups_sweep;
    SR1_sw(sw) = s1_sw / nSetups_sweep;
    SR2_sw(sw) = s2_sw / nSetups_sweep;
    fprintf('squareLen=%4dm | NF=%.0f%% | SR-gap=%.3f | BER-gap=%.2e\n', ...
            sLen, NF_frac_sw(sw), SR2_sw(sw) - SR1_sw(sw), BER1_sw(sw) - BER2_sw(sw));
end

%% ====================================================================
%  PLOT SETTINGS
%% ====================================================================
lw = 2;
ms = 7;
xt = SNR_dB(1:4:end);
c1_li = [0.00 0.45 0.74];
c1_si = [0.47 0.67 0.19];
c2_li = [0.85 0.10 0.10];
c2_si = [0.00 0.65 0.65];
c3_li = [0.49 0.18 0.56];
c3_si = [0.75 0.40 0.85];
cCN = [0.85 0.45 0.00];
cRT = [0.13 0.55 0.13];

%% ====================================================================
%  FIGURE 1 - Sum-rate   [RETAINED + proposed clustering curves added]
%% ====================================================================
figure('Name', 'SR-Comparison', 'Position', [10 520 900 560]);
hold on;
box on;
grid on;
plot(SNR_dB, SR3_lmmse_pf, '-o', 'Color', c3_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Centralized BS (far-field): L-MMSE');
plot(SNR_dB, SR3_cnsic_pf, '-s', 'Color', c3_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Centralized BS (far-field): SIC');
plot(SNR_dB, SR1_lmmse_pf, '-o', 'Color', c1_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free far-field: L-MMSE');
plot(SNR_dB, SR1_cnsic_pf, '-s', 'Color', c1_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free far-field: SIC');
plot(SNR_dB, SR2_lmmse_pf, '-o', 'Color', c2_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free hybrid NF/FF: L-MMSE');
plot(SNR_dB, SR2_cnsic_pf, '-s', 'Color', c2_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free hybrid NF/FF: SIC');
plot(SNR_dB, SR_CNcl_lmmse_pf, '--o', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid CN-clustering: L-MMSE');
plot(SNR_dB, SR_CNcl_cnsic_pf, '--s', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid CN-clustering: SIC');
plot(SNR_dB, SR_RTcl_lmmse_pf, '--^', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid IR-clustering: L-MMSE');
plot(SNR_dB, SR_RTcl_cnsic_pf, '--d', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid IR-clustering: SIC');
xlabel('SNR [dB]', 'FontSize', 13);
ylabel('Sum-rate [bps/Hz]', 'FontSize', 13);
title(sprintf(['Sum-rate, PERFECT CSI: Centralized BS (far-field) vs Cell-free (far-field & hybrid NF/FF) + Clustering\n'...
               'Channel: hybrid near-field/far-field for CF NF/FF & clustering; d_{Ray}(CF)=%.0fm  N_{BS}=%d  squareLen=%dm  L=%d  K=%d  \\eta_{FH}=%.2f'], ...
              d_Ray, N_BS, squareLen, L, K, eta_FH), 'FontSize', 10);
legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);
set(gca, 'XTick', xt);

%% ====================================================================
%  FIGURE 2 - EMPIRICAL BER  [RETAINED + proposed clustering curves added]
%% ====================================================================
figure('Name', 'BER-Comparison-Empirical', 'Position', [920 520 900 560]);
flr2 = 1e-6;
semilogy(SNR_dB, max(aBER2_bs_pf(1, :), flr2), '-o', 'Color', c3_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Centralized BS (far-field): L-MMSE');
hold on;
semilogy(SNR_dB, max(aBER2_bs_pf(2, :), flr2), '-s', 'Color', c3_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Centralized BS (far-field): SIC');
semilogy(SNR_dB, max(aBER2_ff_pf(1, :), flr2), '-o', 'Color', c1_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free far-field: L-MMSE');
semilogy(SNR_dB, max(aBER2_ff_pf(2, :), flr2), '-s', 'Color', c1_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free far-field: SIC');
semilogy(SNR_dB, max(aBER2_hf_pf(1, :), flr2), '-o', 'Color', c2_li, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free hybrid NF/FF: L-MMSE');
semilogy(SNR_dB, max(aBER2_hf_pf(2, :), flr2), '-s', 'Color', c2_si, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Cell-free hybrid NF/FF: SIC');
semilogy(SNR_dB, max(aBER2_cn_pf(1, :), flr2), '--o', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid CN-clustering: L-MMSE');
semilogy(SNR_dB, max(aBER2_cn_pf(2, :), flr2), '--s', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid CN-clustering: SIC');
semilogy(SNR_dB, max(aBER2_rt_pf(1, :), flr2), '--^', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid IR-clustering: L-MMSE');
semilogy(SNR_dB, max(aBER2_rt_pf(2, :), flr2), '--d', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'Hybrid IR-clustering: SIC');
box on;
grid on;
set(gca, 'YMinorGrid', 'on');
xlabel('SNR [dB]', 'FontSize', 13);
ylabel('Empirical BER (16-QAM)', 'FontSize', 13);
title(sprintf(['EMPIRICAL BER (Monte-Carlo 16-QAM), PERFECT CSI: Centralized BS (far-field) vs Cell-free (far-field & hybrid NF/FF) + Clustering\n'...
               'Channel: hybrid near-field/far-field for CF NF/FF & clustering; d_{Ray}(CF)=%.0fm  N_{BS}=%d  squareLen=%dm  L=%d  K=%d  nSym=%d'], ...
              d_Ray, N_BS, squareLen, L, K, nSym2), 'FontSize', 9);
legend('Location', 'southwest', 'FontSize', 8, 'NumColumns', 2);
set(gca, 'XTick', xt);
set(gca, 'YTick', 10.^(-6:0));

%% ====================================================================
%  FIGURE 3 - Deployment area sweep  [UNCHANGED]
%% ====================================================================
figure('Name', 'Area-Sweep', 'Position', [10 30 900 430]);
subplot(1, 2, 1);
yyaxis left;
plot(squareLen_sweep, NF_frac_sw, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', lw, 'MarkerSize', ms + 1);
ylabel('NF pair fraction [%]', 'FontSize', 12);
ylim([0 100]);
yyaxis right;
plot(squareLen_sweep, SR2_sw - SR1_sw, '-s', 'Color', [0.85 0.10 0.10], 'LineWidth', lw, 'MarkerSize', ms + 1);
ylabel('SR gain NF/FF over FF [bps/Hz]', 'FontSize', 11);
xlabel('squareLen [m]', 'FontSize', 12);
title(sprintf('NF fraction & SR gain vs Area\nd_{Ray}=%.0fm SNR=10dB', d_Ray), 'FontSize', 11);
grid on;
box on;
xline(300, '--k', 'LineWidth', 1.5, 'Label', '300m');
legend({'NF pair fraction', 'SR gain'}, 'Location', 'northeast', 'FontSize', 9);
set(gca, 'XTick', squareLen_sweep);
subplot(1, 2, 2);
semilogy(squareLen_sweep, BER1_sw, '-o', 'Color', [0.00 0.45 0.74], 'LineWidth', lw, 'MarkerSize', ms + 1, 'DisplayName', 'CF FF-only');
hold on;
semilogy(squareLen_sweep, BER2_sw, '-s', 'Color', [0.85 0.10 0.10], 'LineWidth', lw, 'MarkerSize', ms + 1, 'DisplayName', 'CF NF/FF');
box on;
grid on;
set(gca, 'YMinorGrid', 'on');
xlabel('squareLen [m]', 'FontSize', 12);
ylabel('BER @SNR=10dB', 'FontSize', 12);
title('BER vs Deployment Area', 'FontSize', 11);
legend('Location', 'southeast', 'FontSize', 9);
set(gca, 'YTick', 10.^(-7:0));
xline(300, '--k', 'LineWidth', 1.5, 'Label', '300m');
set(gca, 'XTick', squareLen_sweep);
sgtitle(sprintf('Area Sweep: L=%d K=%d N=%d', L, K, N), 'FontSize', 12);

% [FIGURE 4 REMOVED per new-paper figure plan -- the NF/FF-full-vs-clustering
%  analytical sum-rate is already covered by the clustering curves in Figure 1.]

%% ====================================================================
%  FIGURE 5 - EMPIRICAL BER with LIST detection, PERFECT vs ESTIMATED CSI
%  Channel model: hybrid near-field/far-field (NF users -> NUSW spherical
%  channel, FF users -> correlated Rayleigh). 16 curves:
%    full-hybrid  (all APs serve)  x {Linear,SIC,List-SIC,CrossAP} x {perfect,estimated}
%    IR-cluster   (rate clustering) x {Linear,SIC,List-SIC,CrossAP} x {perfect,estimated}
%  Solid = perfect CSI (genie, C=0), dashed = estimated CSI (LMMSE pilots).
%% ====================================================================
figure('Name', 'BER-List-CSI', 'Position', [940 620 980 600]);
flr = 1e-6;
mk5  = {'o', 's', '^', 'd'};
det5 = {'Linear', 'SIC', 'List-SIC', 'List+CrossAP (proposed)'};
cHF  = {[0.30 0.55 0.95], [0.15 0.35 0.85], [0.55 0 0], [0.30 0 0]};   % full-hybrid shades
cIR  = {[0.35 0.75 0.35], [0.20 0.60 0.20], [0 0.45 0], [0 0.25 0]};   % IR-cluster shades
hold on;
for d = 1:4
    lwd = lw + 0.4 * (d >= 3);
    semilogy(SNR_dB, max(aBERhf_pf(d, :), flr),  ['-'  mk5{d}], 'Color', cHF{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', ['Hybrid full, perfect: ' det5{d}]);
    semilogy(SNR_dB, max(aBERhf_est(d, :), flr), ['--' mk5{d}], 'Color', cHF{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', ['Hybrid full, estimated: ' det5{d}]);
    semilogy(SNR_dB, max(aBERrt_pf(d, :), flr),  [':'  mk5{d}], 'Color', cIR{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', ['IR-cluster, perfect: ' det5{d}]);
    semilogy(SNR_dB, max(aBERrt_est(d, :), flr), ['-.' mk5{d}], 'Color', cIR{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', ['IR-cluster, estimated: ' det5{d}]);
end
set(gca, 'YScale', 'log');
box on; grid on; set(gca, 'YMinorGrid', 'on');
xlabel('SNR [dB]', 'FontSize', 13);
ylabel('Empirical BER (16-QAM)', 'FontSize', 13);
title(sprintf(['Empirical BER, hybrid NF/FF channel -- perfect vs estimated CSI\n'...
               'Full-hybrid & IR-clustering, 4 detectors; \\tau_{sym}=%d, M=%d, d_{th}=%.2f, nSym=%d'], ...
              tau_sym, M_lst, d_th, nSym), 'FontSize', 10);
legend('Location', 'eastoutside', 'FontSize', 6, 'NumColumns', 1);
set(gca, 'XTick', xt);

% [FIGURE 6 REMOVED per new-paper figure plan -- goodput SE vs detector is
%  redundant with the sum-rate detector figure (Figure 9).]

cDet = {c2_li, [0.55 0 0], [0.20 0.20 0.75], [0.30 0 0]};
mk   = {'-o', '-s', '-^', '-d'};
nmDet = {'IR-Cluster: Linear (L-MMSE)', 'IR-Cluster: SIC', ...
         'IR-Cluster: List-SIC', 'IR-Cluster: List+CrossAP (proposed)'};

%% ====================================================================
%  FIGURE 7 - Symbol NMSE vs detector, PERFECT vs ESTIMATED CSI
%  Channel model: hybrid near-field/far-field, IR (rate) clustering.
%  8 curves = 4 detectors x {perfect, estimated}.
%% ====================================================================
figure('Name', 'NMSE-Detectors-CSI', 'Position', [920 520 900 560]);
hold on; box on; grid on;
set(gca, 'YScale', 'log');
for d = 1:4
    lwd = lw + 0.3 * (d >= 3);
    semilogy(SNR_dB, max(NMSE_det_pf(d, :), 1e-6),  mk{d},        'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (perfect)']);
    semilogy(SNR_dB, max(NMSE_det_est(d, :), 1e-6), ['--' mk{d}(2)], 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (estimated)']);
end
set(gca, 'YMinorGrid', 'on');
xlabel('SNR [dB]', 'FontSize', 13);
ylabel('Symbol NMSE', 'FontSize', 13);
title(sprintf(['Symbol NMSE vs detector, hybrid NF/FF channel -- perfect vs estimated CSI\n'...
               'IR clustering; NMSE_k=E|C(dec_k)-s_k|^2/E|s_k|^2, M=%d, L=%d K=%d N=%d'], modOrder, L, K, N), 'FontSize', 10);
legend('Location', 'southwest', 'FontSize', 7, 'NumColumns', 2);
set(gca, 'XTick', xt);

%% ====================================================================
%  FIGURE 8 - CDF of per-user effective SE, PERFECT vs ESTIMATED CSI
%  Channel model: hybrid near-field/far-field, IR (rate) clustering.
%  8 curves = 4 detectors x {perfect, estimated}.
%% ====================================================================
figure('Name', 'SE-CDF-CSI', 'Position', [10 30 1200 520]);
hold on; box on; grid on;
for d = 1:4
    lwd = lw + 0.3 * (d >= 3);
    xs = sort(seSamp_pf{d});   cdfv = (1:numel(xs))' / numel(xs);
    plot(xs, cdfv, mk{d}, 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms - 1, ...
         'MarkerIndices', 1:max(1, round(numel(xs) / 12)):numel(xs), 'DisplayName', [nmDet{d} ' (perfect)']);
    xe = sort(seSamp_est{d});  cdfe = (1:numel(xe))' / numel(xe);
    plot(xe, cdfe, ['--' mk{d}(2)], 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms - 1, ...
         'MarkerIndices', 1:max(1, round(numel(xe) / 12)):numel(xe), 'DisplayName', [nmDet{d} ' (estimated)']);
end
xlabel('Per-user effective SE [bps/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title(sprintf(['CDF of per-user SE @ %d dB, hybrid NF/FF channel -- perfect vs estimated CSI\n'...
               '(left tail = cell-edge)'], SNR_dB(si_cdf)), 'FontSize', 10);
legend('Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
ylim([0 1]);

%% ====================================================================
%  FIGURE 9 [NEW] - SUM-RATE of all four detection techniques
%  The sum-rate counterpart of Figure 5 (which is BER of all detectors).
%  Left  : goodput sum-rate, bounded by the 16-QAM alphabet.
%  Right : effective-SINR sum-rate from the measured symbol NMSE.
%  See the comment block above nSymTot for why the closed-form Shannon
%  rate of Figures 1 and 4 cannot separate these detectors, and for the
%  NMSE floor that limits the top of the right panel.
%% ====================================================================
figure('Name', 'SumRate-Detectors', 'Position', [10 520 1150 560]);

subplot(1, 2, 1);
hold on; box on; grid on;
for d = 1:4
    plot(SNR_dB, SE_det(d, :), mk{d}, 'Color', cDet{d}, ...
         'LineWidth', lw + 0.3 * (d >= 3), 'MarkerSize', ms, 'DisplayName', nmDet{d});
end
yline(K * prelog * log2(modOrder), '--k', 'LineWidth', 1.2, ...
      'DisplayName', 'Alphabet bound K\cdotprelog\cdotlog_2(M)');
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('Goodput sum-rate [bps/Hz]', 'FontSize', 12);
title(sprintf('Goodput sum-rate, %d-QAM\nSE_k = prelog\\cdotlog_2(M)\\cdot(1-BER_k)', modOrder), 'FontSize', 11);
legend('Location', 'southeast', 'FontSize', 8);
set(gca, 'XTick', xt);

subplot(1, 2, 2);
hold on; box on; grid on;
for d = 1:4
    plot(SNR_dB, SEeff_det(d, :), mk{d}, 'Color', cDet{d}, ...
         'LineWidth', lw + 0.3 * (d >= 3), 'MarkerSize', ms, 'DisplayName', nmDet{d});
end
yline(K * prelog * log2(1 + nSymTot), ':k', 'LineWidth', 1.2, ...
      'DisplayName', 'NMSE measurement floor');
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('Effective-SINR sum-rate [bps/Hz]', 'FontSize', 12);
title(sprintf('Effective-SINR sum-rate\nSE_k = prelog\\cdotlog_2(1+1/NMSE_k), floor 1/%d', nSymTot), 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 8);
set(gca, 'XTick', xt);
sgtitle(sprintf('Sum-rate vs detector, rate clustering, hybrid NF/FF (L=%d K=%d N=%d)', L, K, N), 'FontSize', 12);

%% ====================================================================
%  FIGURE 10 [NEW] - DETECTOR COMPARISON (BER) + DI RENNA-STYLE
%  ANALYTICAL SUM-RATE (cluster fusion vs cross-AP fusion)
%
%  LEFT  : empirical BER of the four detectors (same data as Figure 5,
%          repeated here so the detector separation sits next to the
%          rate). Cross-AP separates from List-SIC HERE, at the decision
%          level - not in the analytical rate, by construction.
%  RIGHT : analytical achievable sum-rate, Di Renna Thm 1 / eq. (31)-(33)
%          form, prelog*log2(1+Gamma_k) averaged over realizations. Two
%          curves: cluster-only fusion (Gamma_k^clu) and cross-AP fusion
%          (Gamma_k^all). The shaded gap is deltaR, the analytical rate
%          value of the observations clustering discards. It is >= 0 and
%          collapses to 0 under full cooperation (Theorem 1). This is a
%          SYSTEM-level rate figure; it does NOT separate List-SIC from
%          Cross-AP (both share the combiner SINR).
%% ====================================================================
figure('Name', 'DetectorCompare-and-AnalyticalRate', 'Position', [10 30 1200 520]);

% ---- LEFT: empirical BER of the four detectors, perfect vs estimated CSI ----
subplot(1, 2, 1);
hold on; box on; grid on;
for d = 1:4
    lwd = lw + 0.3 * (d >= 3);
    semilogy(SNR_dB, max(aBERrt_pf(d, :), 1e-6),  mk{d},          'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (perfect)']);
    semilogy(SNR_dB, max(aBERrt_est(d, :), 1e-6), ['--' mk{d}(2)], 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (estimated)']);
end
set(gca, 'YScale', 'log', 'XTick', xt);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('BER (simulated)', 'FontSize', 12);
title(sprintf('Detector comparison: empirical BER\nhybrid NF/FF channel, perfect vs estimated CSI'), 'FontSize', 10);
legend('Location', 'southwest', 'FontSize', 6, 'NumColumns', 2);
ylim([1e-5 1]);

% ---- RIGHT: Di Renna-style analytical sum-rate, cluster vs cross-AP ----
subplot(1, 2, 2);
hold on; box on; grid on;
% shade the gap deltaR first so the lines sit on top
xf = [SNR_dB, fliplr(SNR_dB)];
yf = [SR_anaClu, fliplr(SR_anaAll)];
fill(xf, yf, [0.20 0.50 0.90], 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
     'DisplayName', '\DeltaR (discarded APs)');
plot(SNR_dB, SR_anaClu, '-o', 'Color', [0.20 0.30 0.60], 'LineWidth', 1.6, ...
     'MarkerSize', ms, 'DisplayName', 'Cluster fusion  R^{clu}');
plot(SNR_dB, SR_anaAll, '-d', 'Color', [0.75 0.20 0.20], 'LineWidth', 1.9, ...
     'MarkerSize', ms, 'DisplayName', 'Cross-AP fusion  R^{all}');
set(gca, 'XTick', xt);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('Achievable sum-rate [bps/Hz]', 'FontSize', 12);
title({'Analytical sum-rate (Di Renna form), hybrid NF/FF', ...
       'R = \Sigma_k prelog\cdotlog_2(1+\Gamma_k) (combiner SINR, CSI-agnostic)'}, 'FontSize', 10);
legend('Location', 'northwest', 'FontSize', 8);
sgtitle(sprintf('Detector BER (perfect vs estimated CSI) vs analytical cross-AP rate gain, hybrid NF/FF (L=%d K=%d N=%d)', L, K, N), 'FontSize', 11);

% ---- console readout of the analytical rate gap ----
fprintf('\n--- Analytical sum-rate (Di Renna form): cluster vs cross-AP ---\n');
fprintf('%-6s %14s %14s %12s\n', 'SNR', 'R_clu [b/Hz]', 'R_all [b/Hz]', 'deltaR');
fprintf('%s\n', repmat('-', 1, 50));
for si = 1:nSNR
    fprintf('  %4d %14.3f %14.3f %12.3f\n', ...
            SNR_dB(si), SR_anaClu(si), SR_anaAll(si), SR_anaGap(si));
end
if max(SR_anaGap) < 1e-6
    fprintf(['NOTE: deltaR is ~0 at all SNR. With this cluster size the\n' ...
             'serving set already captures nearly all the SINR, so the\n' ...
             'discarded APs carry little rate. Shrink the cluster (smaller\n' ...
             'tau_p / fronthaul) to open the gap; this matches Theorem 1.\n']);
end

% [FIGURE 11 REMOVED per new-paper figure plan -- the analytical all-detector
%  sum-rate coincides for the SIC family by construction (shared combiner SINR),
%  so the detector separation is carried by BER (Figures 5 and 10) instead.]

%% ====================================================================
%  FIGURE 12 [NEW] - COMPUTATIONAL COMPLEXITY vs NUMBER OF USERS K
%
%  Complex multiplications per channel use, per detector, as K grows, with
%  L, N, M_lst, maxBr fixed at their configured values. This is the
%  "affordability" figure that pairs with BER: the classic list-detector
%  argument is near-ML accuracy at SIC-like cost (cf. Li Fig 10, Arevalo
%  Fig 1, Di Renna Fig 5 / Table II).
%
%  COUNTING MODEL (dominant complex multiplications, per user):
%   - Local MMSE combiner per AP: N x N solve ~ N^3/3, plus N^2 apply.
%     A user served by |S| APs pays |S| of these. We take |S| ~ Kbar, the
%     mean cluster size, so the count reflects the clustered architecture.
%   - Linear (L-MMSE):     cLin  = Kbar*(N^3/3 + N^2)
%   - Hard-SIC:            cSIC  = cLin + Kbar*K*N            (sequential
%                          cancellation across the K stages)
%   - List-SIC:            cList = cSIC + M_lst*maxBr*Kbar*N  (list search
%                          over M candidates and maxBr branches)
%   - Cross-AP List-SIC:   cXap  = cList + eta*(L-Kbar)*N     (the gated
%                          all-AP fusion: only the non-serving L-Kbar APs,
%                          only for the fraction eta of unreliable stages)
%
%  eta_meas is now COUNTED inside ber_case (fraction of SIC stages that
%  actually invoke cross-AP fusion) at the reference SNR - not inferred from
%  rate. All four share the same leading N^3 term, so the curves stay close;
%  cross-AP adds a sub-dominant, measured term.
%% ====================================================================
Ksweep = 2:2:12;
Kbar   = max(mean(loadBSR_acc(:)), 1);          % mean cluster size (rate rule, accumulated)
% MEASURED gating fraction at the reference SNR: the fraction of SIC stages
% where cross-AP fusion was actually invoked, counted inside ber_case.
eta_meas = mean(etaXap_acc(:, snr_ref_idx));

cLin  = zeros(size(Ksweep));
cSIC  = zeros(size(Ksweep));
cList = zeros(size(Ksweep));
cXap  = zeros(size(Ksweep));
for ii = 1:numel(Ksweep)
    Kv = Ksweep(ii);
    base   = Kbar * (N^3 / 3 + N^2);
    cLin(ii)  = base;
    cSIC(ii)  = base + Kbar * Kv * N;
    cList(ii) = cSIC(ii) + M_lst * maxBr * Kbar * N;
    cXap(ii)  = cList(ii) + eta_meas * (L - Kbar) * N;
end

figure('Name', 'Complexity-vs-K', 'Position', [10 30 760 560]);
semilogy(Ksweep, cLin,  '-o', 'Color', cDet{1}, 'LineWidth', lw, 'MarkerSize', ms, ...
         'DisplayName', 'Linear (L-MMSE)'); hold on; box on; grid on;
semilogy(Ksweep, cSIC,  '-s', 'Color', cDet{2}, 'LineWidth', lw, 'MarkerSize', ms, ...
         'DisplayName', 'Hard-SIC');
semilogy(Ksweep, cList, '-^', 'Color', cDet{3}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
         'DisplayName', 'List-SIC');
semilogy(Ksweep, cXap,  '-d', 'Color', cDet{4}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
         'DisplayName', 'List+CrossAP (proposed)');
xlabel('Number of users K', 'FontSize', 12);
ylabel('Complex multiplications per channel use', 'FontSize', 12);
title(sprintf('Detector complexity vs K  (L=%d, N=%d, \\eta_{meas}=%.2f)', L, N, eta_meas), 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9);

fprintf('\n--- Complexity vs K (Fig 12), Kbar=%.2f, eta_meas=%.3f ---\n', Kbar, eta_meas);
fprintf('%-4s %14s %14s %14s %16s\n', 'K', 'Linear', 'Hard-SIC', 'List-SIC', 'List+CrossAP');
for ii = 1:numel(Ksweep)
    fprintf('%-4d %14.3e %14.3e %14.3e %16.3e\n', ...
            Ksweep(ii), cLin(ii), cSIC(ii), cList(ii), cXap(ii));
end

%% ====================================================================
%  FIGURE 13 [NEW] - FRONTHAUL OVERHEAD vs SNR
%
%  Average scalars forwarded per user, Bbar_k = |S_k| + eta(SNR)*(L-|S_k|),
%  where eta(SNR) is the MEASURED fraction of stages whose cluster estimate
%  is unreliable (the gate that triggers cross-AP fusion). The point of the
%  figure: the cross-AP cost is not fixed - it SHRINKS as SNR rises, because
%  fewer stages are unreliable, so the gate fires less. At high SNR the
%  overhead collapses toward the cluster-only payload |S_k|. This directly
%  answers the reviewer question "what does cross-AP cost?".
%
%  eta(SNR) is MEASURED: inside ber_case we count the fraction of SIC
%  stages where cross-AP fusion is invoked (the SAC shadow event), averaged
%  over realizations and setups. This is the true gating rate, not a proxy.
%% ====================================================================
Kbar_fh = max(mean(loadBSR_acc(:)), 1);         % mean cluster size (accumulated)
% MEASURED cross-AP invocation rate per SNR, counted inside ber_case.
eta_snr = mean(etaXap_acc, 1);                   % 1 x nSNR, averaged over setups
Bbar_cluster = Kbar_fh * ones(1, nSNR);              % cluster-only payload
Bbar_xap     = Kbar_fh + eta_snr .* (L - Kbar_fh);   % gated cross-AP payload
% Same payload measured with PERFECT CSI (for perfect-CSI papers, e.g. ICASSP Fig. 2).
eta_snr_pf   = mean(etaXap_pf_acc, 1);
Bbar_xap_pf  = Kbar_fh + eta_snr_pf .* (L - Kbar_fh);

figure('Name', 'Fronthaul-vs-SNR', 'Position', [790 30 760 560]);
plot(SNR_dB, Bbar_xap, '-d', 'Color', cDet{4}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
     'DisplayName', 'Cross-AP (gated)  |S_k|+\eta(L-|S_k|)'); hold on; box on; grid on;
plot(SNR_dB, Bbar_cluster, '--s', 'Color', cDet{2}, 'LineWidth', lw, 'MarkerSize', ms, ...
     'DisplayName', 'Cluster-only  |S_k|');
yline(L, ':k', 'LineWidth', 1.2, 'DisplayName', 'All-AP (ungated) upper bound L');
set(gca, 'XTick', xt);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('Average fronthaul scalars per user', 'FontSize', 12);
title(sprintf('Fronthaul overhead vs SNR  (L=%d, mean |S_k|=%.2f)', L, Kbar_fh), 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 9);
ylim([0 L + 0.5]);

fprintf('\n--- Fronthaul overhead vs SNR (Fig 13), mean |S_k|=%.2f ---\n', Kbar_fh);
fprintf('etaList = measured List-SIC trigger rate; etaXap = measured Cross-AP invocation rate\n');
etaList_snr = mean(etaList_acc, 1);
fprintf('%-6s %10s %10s %14s %14s\n', 'SNR', 'etaList', 'etaXap', 'cluster |S_k|', 'cross-AP Bbar');
for si = 1:nSNR
    fprintf('  %4d %10.3f %10.3f %14.3f %14.3f\n', ...
            SNR_dB(si), etaList_snr(si), eta_snr(si), Bbar_cluster(si), Bbar_xap(si));
end

%% ====================================================================
%  CONSOLE SUMMARY
%% ====================================================================
fprintf('\n=== RESULTS: PERFECT-CSI sum-rate (L-MMSE), matches Figure 1 (squareLen=%dm, eta_FH=%.2f) ===\n', squareLen, eta_FH);
fprintf('%-6s|%-9s|%-9s|%-9s|%-9s|%-9s\n', 'SNR', 'SR-FF', 'SR-NF/FF', 'SR-CN-cl', 'SR-RT-cl', 'SR-Cent');
fprintf('%s\n', repmat('-', 1, 60));
for si = 1:nSNR
    fprintf('  %4d|%9.2f|%9.2f|%9.2f|%9.2f|%9.2f\n', ...
            SNR_dB(si), SR1_lmmse_pf(si), SR2_lmmse_pf(si), SR_CNcl_lmmse_pf(si), ...
            SR_RTcl_lmmse_pf(si), SR3_lmmse_pf(si));
end

fprintf('\n=== AP-SELECTION DISAGREEMENT: CHANNEL NORM vs RATE (APs per user) ===\n');
fprintf('If this is 0.000 the two clustering rules select the SAME APs and\n');
fprintf('every difference between their curves is Monte-Carlo noise.\n');
fprintf('  Rate-Cluster vs CN : %.3f\n\n', mean(clDiff_acc(:)));

fprintf('\n--- Sum-rate by detector (Fig 9): goodput | effective-SINR ---\n');
fprintf('%-6s %28s %28s\n', 'SNR', 'goodput [bps/Hz]', 'eff-SINR [bps/Hz]');
fprintf('%-6s %6s %6s %6s %6s   %6s %6s %6s %6s\n', '', ...
        'Lin', 'SIC', 'List', 'XAP', 'Lin', 'SIC', 'List', 'XAP');
for si2 = 1:nSNR
    fprintf('  %4d %6.2f %6.2f %6.2f %6.2f   %6.1f %6.1f %6.1f %6.1f\n', SNR_dB(si2), ...
            SE_det(1, si2), SE_det(2, si2), SE_det(3, si2), SE_det(4, si2), ...
            SEeff_det(1, si2), SEeff_det(2, si2), SEeff_det(3, si2), SEeff_det(4, si2));
end
if any(SEeff_capped(:))
    fprintf('NOTE: at least one detector hit the NMSE measurement floor.\n');
    fprintf('Effective-SINR values at the top of the right panel are LIMITED\n');
    fprintf('BY nSym, not by the detector. Raise nSym before quoting them.\n');
end

fprintf('\n--- Empirical BER (Hybrid full): Linear / SIC / List / CrossAP ---\n');
for si2 = 1:nSNR
    fprintf('  %4d dB | %8.4f | %8.4f | %8.4f | %8.4f  (CrossAP gain vs List %+.1f%%)\n', ...
            SNR_dB(si2), aBERhf(1, si2), aBERhf(2, si2), aBERhf(3, si2), aBERhf(4, si2), ...
            100 * (aBERhf(3, si2) - aBERhf(4, si2)) / max(aBERhf(3, si2), eps));
end

fprintf('\nCF avg NF pairs: %.1f%%\n', 100 * nNF_total / (nSetups * L * K));
fprintf('Mean cluster load: |M_CN|=%.2f  |M_BSR|=%.2f  (full L=%d)\n', ...
        mean(loadCN_acc(:)), mean(loadBSR_acc(:)), L);
fprintf('\nReport these numbers as measured. Do NOT sweep clu_div, tau_sym or\n');
fprintf('pilot_mode until a favourable sign appears; that is genie tuning.\n');
fprintf('Done.\n');

%% ====================================================================
%  NEW FIGURE 1 [FIGURE REARRANGEMENT ONLY - no technical change]
%  Left  : BER of all four detectors  (same data as the detector-BER figure)
%  Right : CDF of per-user effective SE (same data as the CDF figure)
%  Both already existed as separate figures; here they are combined into one
%  figure with two subplots.
%% ====================================================================
figure('Name', 'NewFig1-BER-and-CDF', 'Position', [60 120 1250 520]);

% ---- LEFT: empirical BER of the four detectors, perfect vs estimated CSI ----
subplot(1, 2, 1);
hold on; box on; grid on;
for d = 1:4
    lwd = lw + 0.3 * (d >= 3);
    semilogy(SNR_dB, max(aBERrt_pf(d, :), 1e-6),  mk{d},          'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (perfect)']);
    semilogy(SNR_dB, max(aBERrt_est(d, :), 1e-6), ['--' mk{d}(2)], 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms, 'DisplayName', [nmDet{d} ' (estimated)']);
end
set(gca, 'YScale', 'log', 'XTick', xt);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);
title(sprintf('BER of all four detectors\nhybrid NF/FF channel, perfect vs estimated CSI'), 'FontSize', 10);
legend('Location', 'southwest', 'FontSize', 6, 'NumColumns', 2);
ylim([1e-5 1]);

% ---- RIGHT: CDF of per-user effective SE, perfect vs estimated CSI ----
subplot(1, 2, 2);
hold on; box on; grid on;
for d = 1:4
    lwd = lw + 0.3 * (d >= 3);
    xs = sort(seSamp_pf{d});   cdfv = (1:numel(xs))' / numel(xs);
    plot(xs, cdfv, mk{d}, 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms - 1, ...
         'MarkerIndices', 1:max(1, round(numel(xs) / 12)):numel(xs), 'DisplayName', [nmDet{d} ' (perfect)']);
    xe = sort(seSamp_est{d});  cdfe = (1:numel(xe))' / numel(xe);
    plot(xe, cdfe, ['--' mk{d}(2)], 'Color', cDet{d}, 'LineWidth', lwd, 'MarkerSize', ms - 1, ...
         'MarkerIndices', 1:max(1, round(numel(xe) / 12)):numel(xe), 'DisplayName', [nmDet{d} ' (estimated)']);
end
xlabel('Per-user effective SE [bps/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title(sprintf('CDF of per-user SE @ %d dB\nhybrid NF/FF channel, perfect vs estimated CSI', SNR_dB(si_cdf)), 'FontSize', 10);
legend('Location', 'northwest', 'FontSize', 6, 'NumColumns', 2);
ylim([0 1]);

%% ====================================================================
%  NEW FIGURE 2 [FIGURE REARRANGEMENT ONLY - no technical change]
%  Left  : average fronthaul scalars per user (same data as the fronthaul figure)
%  Right : detector complexity vs number of users K (same data as the complexity figure)
%  Both already existed as separate figures; here they are combined into one
%  figure with two subplots.
%% ====================================================================
figure('Name', 'NewFig2-Fronthaul-and-Complexity', 'Position', [60 120 1250 520]);

% ---- LEFT: average fronthaul scalars per user (same as Fig 13) ----
subplot(1, 2, 1);
plot(SNR_dB, Bbar_xap, '-d', 'Color', cDet{4}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
     'DisplayName', 'Cross-AP (gated)  |S_k|+\eta(L-|S_k|)'); hold on; box on; grid on;
plot(SNR_dB, Bbar_cluster, '--s', 'Color', cDet{2}, 'LineWidth', lw, 'MarkerSize', ms, ...
     'DisplayName', 'Cluster-only  |S_k|');
yline(L, ':k', 'LineWidth', 1.2, 'DisplayName', 'All-AP (ungated) upper bound L');
set(gca, 'XTick', xt);
xlabel('SNR [dB]', 'FontSize', 12);
ylabel('Average fronthaul scalars per user', 'FontSize', 12);
title(sprintf('Fronthaul overhead vs SNR  (L=%d, mean |S_k|=%.2f)', L, Kbar_fh), 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 9);
ylim([0 L + 0.5]);

% ---- RIGHT: detector complexity vs number of users K (same as Fig 12) ----
subplot(1, 2, 2);
semilogy(Ksweep, cLin,  '-o', 'Color', cDet{1}, 'LineWidth', lw, 'MarkerSize', ms, ...
         'DisplayName', 'Linear (L-MMSE)'); hold on; box on; grid on;
semilogy(Ksweep, cSIC,  '-s', 'Color', cDet{2}, 'LineWidth', lw, 'MarkerSize', ms, ...
         'DisplayName', 'Hard-SIC');
semilogy(Ksweep, cList, '-^', 'Color', cDet{3}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
         'DisplayName', 'List-SIC');
semilogy(Ksweep, cXap,  '-d', 'Color', cDet{4}, 'LineWidth', lw + 0.3, 'MarkerSize', ms + 1, ...
         'DisplayName', 'List+CrossAP (proposed)');
xlabel('Number of users K', 'FontSize', 12);
ylabel('Complex multiplications per channel use', 'FontSize', 12);
title(sprintf('Detector complexity vs K  (L=%d, N=%d)', L, N), 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 9);

%% ====================================================================
%  NEW FIGURE A [CSI STUDY] - BER of all four detectors: ESTIMATED vs PERFECT
%  Hybrid NF/FF: near-field links use the NUSW near-field LMMSE estimate, far-
%  field links the standard LMMSE estimate; perfect CSI is the genie upper bound
%  (true channel, zero error covariance). Solid = estimated, dashed = perfect.
%% ====================================================================
if csi_study
    cDcsi = {[0.85 0.10 0.10], [0.55 0 0], [0.20 0.20 0.75], [0.30 0 0]};
    mkcsi = {'o', 's', '^', 'd'};
    nmcsi = {'Linear (L-MMSE)', 'SIC', 'List-SIC', 'List+CrossAP (proposed)'};
    flrC  = 1e-6;
    figure('Name', 'CSI-BER-Perfect-vs-Estimated', 'Position', [50 70 960 640]);
    hold on; box on; grid on;
    hEst = gobjects(1, 4);
    for d = 1:4
        hEst(d) = semilogy(SNR_dB, max(aBERrt_est(d, :), flrC), ['-' mkcsi{d}], 'Color', cDcsi{d}, ...
            'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nmcsi{d} ' - estimated CSI']);
    end
    for d = 1:4
        semilogy(SNR_dB, max(aBERrt_pf(d, :), flrC), ['--' mkcsi{d}], 'Color', cDcsi{d}, ...
            'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nmcsi{d} ' - perfect CSI']);
    end
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB);
    xlim([SNR_dB(1) SNR_dB(end)]); ylim([1e-5 1]);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Empirical BER (16-QAM)', 'FontSize', 13);
    if genie_csi
        estTxt = 'NF links: LEGACY genie-prior LMMSE';
    else
        estTxt = 'NF links: polar-domain OMP';
    end
    title(sprintf('BER of all detectors: estimated vs perfect CSI (hybrid NF/FF)\n%s;  FF links: LMMSE;  \\tau_{sym}=%d', estTxt, tau_sym), 'FontSize', 11);
    legend('Location', 'southwest', 'FontSize', 8, 'NumColumns', 2);

    %% ================================================================
    %  NEW FIGURE B [CSI STUDY] - SUM-RATE of all four detectors: EST vs PERFECT
    %% ================================================================
    figure('Name', 'CSI-SumRate-Perfect-vs-Estimated', 'Position', [80 50 960 640]);
    hold on; box on; grid on;
    for d = 1:4
        plot(SNR_dB, SE_det_est(d, :), ['-' mkcsi{d}], 'Color', cDcsi{d}, ...
            'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nmcsi{d} ' - estimated CSI']);
    end
    for d = 1:4
        plot(SNR_dB, SE_det_pf(d, :), ['--' mkcsi{d}], 'Color', cDcsi{d}, ...
            'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nmcsi{d} ' - perfect CSI']);
    end
    set(gca, 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Goodput sum-rate [bps/Hz]', 'FontSize', 13);
    title(sprintf('Sum-rate of all detectors: estimated vs perfect CSI (hybrid NF/FF)\nSE_k = prelog\\cdotlog_2(M)(1-BER_k), M=%d', modOrder), 'FontSize', 11);
    legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);

    fprintf('\n[CSI] BER at %ddB (proposed detector): estimated=%.2e  perfect=%.2e\n', ...
            SNR_dB(end), aBERrt_est(4, end), aBERrt_pf(4, end));
    fprintf('[CSI] Sum-rate at %ddB (proposed detector): estimated=%.2f  perfect=%.2f bps/Hz\n', ...
            SNR_dB(end), SE_det_est(4, end), SE_det_pf(4, end));
end

% shared style for the CSI / NF-vs-FF figures
cD4 = {[0.85 0.10 0.10], [0.55 0 0], [0.20 0.20 0.75], [0.30 0 0]};
mk4 = {'o', 's', '^', 'd'};
nm4 = {'Linear (L-MMSE)', 'SIC', 'List-SIC', 'List+CrossAP (proposed)'};
flrC = 1e-6;
cNF = [0.30 0 0];  cFF = [0 0.35 0.65];

if csi_study
    %% NEW FIGURE C - ACHIEVABLE (effective-SINR) sum-rate vs SNR ----------
    %  Unlike the goodput (which meets the 16-QAM alphabet bound at high SNR,
    %  so all detectors converge), this rate is unbounded and keeps growing,
    %  so the detectors fan out from low to high SNR.
    figure('Name', 'CSI-AchievableRate-vs-SNR', 'Position', [60 70 960 640]);
    hold on; box on; grid on;
    for d = 1:4
        plot(SNR_dB, SEeff_det_est(d, :), ['-' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2, ...
            'MarkerSize', 7, 'DisplayName', [nm4{d} ' - estimated CSI']);
    end
    for d = 1:4
        plot(SNR_dB, SEeff_det_pf(d, :), ['--' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2, ...
            'MarkerSize', 7, 'DisplayName', [nm4{d} ' - perfect CSI']);
    end
    set(gca, 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Achievable sum-rate [bps/Hz]', 'FontSize', 13);
    title(sprintf('Achievable (effective-SINR) sum-rate vs SNR, hybrid NF/FF channel\nestimated vs perfect CSI'), 'FontSize', 11);
    legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);

    %% NEW FIGURE H - ROBUSTNESS to CSI error (retained rate fraction) ------
    %  Metric: retained fraction = SE(estimated)/SE(perfect) in [0,1]. A value
    %  near 1 means the detector loses little when CSI is estimated rather than
    %  genie -- i.e. it is ROBUST to estimation error. Higher curve = more
    %  robust. This replaces the earlier absolute-loss plot, on which a detector
    %  that both starts and stays highest in absolute rate necessarily showed
    %  the largest absolute gap and looked "least robust" -- an artefact of
    %  scale, not of sensitivity. The fraction normalises that out.
    figure('Name', 'CSI-Robustness-Retained-Rate-Fraction', 'Position', [90 50 900 600]);
    hold on; box on; grid on;
    for d = 1:4
        frac = SE_det_est(d, :) ./ max(SE_det_pf(d, :), eps);
        frac = min(max(frac, 0), 1);
        plot(SNR_dB, frac, ['-' mk4{d}], 'Color', cD4{d}, ...
            'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', nm4{d});
    end
    set(gca, 'XTick', SNR_dB); ylim([0 1.02]);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Retained rate fraction  SE_{est}/SE_{perfect}', 'FontSize', 13);
    title(sprintf('Robustness to CSI error, hybrid NF/FF channel: retained goodput fraction per detector\n(closer to 1 = more robust to estimation error)'), 'FontSize', 11);
    legend('Location', 'southeast', 'FontSize', 9);

    %% NEW FIGURE G - estimator NMSE vs SNR --------------------------------
    figure('Name', 'Estimator-NMSE-vs-SNR', 'Position', [110 40 860 560]);
    if genie_csi
        nfName = 'Hybrid: legacy genie-prior LMMSE';
    else
        nfName = 'Hybrid: polar-domain OMP (NF) + LMMSE (FF)';
    end
    semilogy(SNR_dB, max(nmse_est, 1e-6), '-o', 'Color', cNF, 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', nfName); hold on; box on; grid on;
    if ff_study
        semilogy(SNR_dB, max(nmse_ff, 1e-6), '-s', 'Color', cFF, 'LineWidth', 2, 'MarkerSize', 7, ...
            'DisplayName', 'Far-field LMMSE');
    end
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB);
    if genie_csi
        xlabel('SNR [dB]', 'FontSize', 13); ylabel('Estimator NMSE = tr(C)/tr(R)', 'FontSize', 13);
        title(sprintf(['Channel-estimation quality vs SNR, LEGACY genie-prior path\n'...
                       'NF NMSE is optimistic: the prior R2 is built from the true channel']), 'FontSize', 10);
    else
        xlabel('SNR [dB]', 'FontSize', 13); ylabel('Empirical NMSE  E||h-\^h||^2 / E||h||^2', 'FontSize', 13);
        title(sprintf(['Channel-estimation quality vs SNR, \\tau_p=%d orthogonal pilots\n'...
                       'NF pairs: polar-domain OMP;  FF pairs: statistical LMMSE'], tau_sym), 'FontSize', 10);
    end
    legend('Location', 'southwest', 'FontSize', 9);
end

if ff_study
    %% NEW FIGURE D - Proposed Cross-AP: BER, near-field vs far-field ------
    figure('Name', 'CrossAP-BER-NF-vs-FF', 'Position', [60 70 900 620]);
    d = 4;
    semilogy(SNR_dB, max(aBERrt_pf(d, :), flrC),   '-d',  'Color', cNF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Near-field, perfect CSI'); hold on;
    semilogy(SNR_dB, max(aBERrt_est(d, :), flrC),  '--d', 'Color', cNF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Near-field, estimated CSI');
    semilogy(SNR_dB, max(aBERrt_ff_pf(d, :), flrC),  '-o',  'Color', cFF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Far-field, perfect CSI');
    semilogy(SNR_dB, max(aBERrt_ff_est(d, :), flrC), '--o', 'Color', cFF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Far-field, estimated CSI');
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB); ylim([1e-5 1]);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('BER (16-QAM)', 'FontSize', 13);
    title(sprintf('Proposed Cross-AP List-SIC: BER, near-field vs far-field channel model\nsolid = perfect CSI, dashed = estimated CSI'), 'FontSize', 10);
    legend('Location', 'southwest', 'FontSize', 9);

    %% NEW FIGURE E - Proposed Cross-AP: sum-rate, near-field vs far-field -
    figure('Name', 'CrossAP-SumRate-NF-vs-FF', 'Position', [90 50 900 620]);
    d = 4;
    plot(SNR_dB, SE_det_pf(d, :),     '-d',  'Color', cNF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Near-field, perfect CSI'); hold on; box on; grid on;
    plot(SNR_dB, SE_det_est(d, :),    '--d', 'Color', cNF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Near-field, estimated CSI');
    plot(SNR_dB, SE_det_ff_pf(d, :),  '-o',  'Color', cFF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Far-field, perfect CSI');
    plot(SNR_dB, SE_det_ff_est(d, :), '--o', 'Color', cFF, 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'Far-field, estimated CSI');
    set(gca, 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Goodput sum-rate [bps/Hz]', 'FontSize', 13);
    title(sprintf('Proposed Cross-AP List-SIC: sum-rate, near-field vs far-field channel model\nsolid = perfect CSI, dashed = estimated CSI'), 'FontSize', 10);
    legend('Location', 'northwest', 'FontSize', 9);

    %% NEW FIGURE F - Other detectors (Linear, SIC, List-SIC): NF vs FF ----
    figure('Name', 'OtherDetectors-NF-vs-FF', 'Position', [40 40 1220 540]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for d = 1:3
        semilogy(SNR_dB, max(aBERrt_pf(d, :), flrC),    ['-' mk4{d}],  'Color', cD4{d}, 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (near-field)']);
        semilogy(SNR_dB, max(aBERrt_ff_pf(d, :), flrC), ['--' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (far-field)']);
    end
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB); ylim([1e-5 1]);
    xlabel('SNR [dB]', 'FontSize', 12); ylabel('BER (16-QAM)', 'FontSize', 12);
    title('BER (perfect CSI)', 'FontSize', 11); legend('Location', 'southwest', 'FontSize', 8);
    subplot(1, 2, 2); hold on; box on; grid on;
    for d = 1:3
        plot(SNR_dB, SE_det_pf(d, :),    ['-' mk4{d}],  'Color', cD4{d}, 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (near-field)']);
        plot(SNR_dB, SE_det_ff_pf(d, :), ['--' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (far-field)']);
    end
    set(gca, 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 12); ylabel('Goodput sum-rate [bps/Hz]', 'FontSize', 12);
    title('Sum-rate (perfect CSI)', 'FontSize', 11); legend('Location', 'northwest', 'FontSize', 8);
    sgtitle('Other detectors (Linear, SIC, List-SIC): near-field vs far-field channel model, perfect CSI', 'FontSize', 11);
end

%% ====================================================================
%  [CODE C] SPARSE-ARRAY ANALYSIS: FIGURES SPA-1 ... SPA-14
%% ====================================================================
cA   = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.13 0.55 0.13], [0.49 0.18 0.56]};
mA   = {'o', 's', 'd', '^'};
detN = {'Linear (L-MMSE)', 'Hard-SIC', 'List-SIC', 'List+CrossAP (proposed)'};
nGeo = numel(A_geo);

% ---------------- SPA-1: array geometries ----------------
figure('Name', 'SPA-1 AP array geometries', 'Position', [40 40 900 420]);
hold on; box on; grid on;
for a = 1:nGeo
    plot(A_geo{a}.pos, a * ones(A_geo{a}.N, 1), mA{a}, 'Color', cA{a}, 'MarkerSize', 4, ...
         'MarkerFaceColor', cA{a}, 'DisplayName', sprintf('%s: N=%d, D=%.2f m, R=%.0f m', ...
         arr_names{a}, A_geo{a}.N, A_geo{a}.aperture, 2 * A_geo{a}.aperture^2 / lambda));
end
ylim([0 nGeo + 1]);
xlabel('Element position along the array [m]'); ylabel('Array');
title('Element positions (all on the \lambda/2 grid)'); legend('Location', 'southoutside');

% ---------------- SPA-2: beam pattern, grating lobes, refocusing, beam depth ----------------
r0 = 40;  s0 = 0.3;  s_gl = s0 - 2 / usa_eta;          % first USA grating lobe (Zhou Lemma 5)
r_gl = r0 * (1 - s_gl^2) / (1 - s0^2);                  % its focal distance (Zhou Lemma 6)
s_ax = linspace(-1, 1, 2001);
r_ax = linspace(5, 400, 1600);
figure('Name', 'SPA-2 Beam pattern grating lobes and beam depth', 'Position', [40 40 1300 420]);
for a = 1:3
    pos_a = A_geo{a}.pos;
    a0 = nusw_vec(pos_a, r0, s0, lambda);
    subplot(1, 3, 1); hold on; box on; grid on;
    plot(s_ax, abs(nusw_vec(pos_a, r0 * ones(size(s_ax)), s_ax, lambda)' * a0), '-', 'Color', cA{a}, ...
         'LineWidth', 1.4, 'DisplayName', arr_names{a});
    subplot(1, 3, 2); hold on; box on; grid on;
    plot(r_ax, abs(nusw_vec(pos_a, r_ax, s_gl * ones(size(r_ax)), lambda)' * a0), '-', 'Color', cA{a}, ...
         'LineWidth', 1.4, 'DisplayName', arr_names{a});
    subplot(1, 3, 3); hold on; box on; grid on;
    plot(r_ax, abs(nusw_vec(pos_a, r_ax, s0 * ones(size(r_ax)), lambda)' * a0), '-', 'Color', cA{a}, ...
         'LineWidth', 1.4, 'DisplayName', arr_names{a});
end
subplot(1, 3, 1); xlabel('sin\theta at r = r_0'); ylabel('Normalised gain |a^H a_0|');
title(sprintf('Angular cut, focus at r_0=%d m, sin\\theta_0=%.1f', r0, s0)); legend('Location', 'south');
subplot(1, 3, 2); plot([r_gl r_gl], [0 1], ':k', 'LineWidth', 1.2, 'DisplayName', 'r_u (Zhou Lemma 6)');
xlabel('Range r [m]'); ylabel('Normalised gain');
title(sprintf('Range cut at the grating-lobe angle sin\\theta_u=%.1f', s_gl)); legend('Location', 'northeast');
subplot(1, 3, 3); xlabel('Range r [m]'); ylabel('Normalised gain');
title('Range cut at the focus angle (beam depth)'); legend('Location', 'northeast');

% ---------------- SPA-3: near-field and focusing fractions vs area ----------------
area_ax = [100 200 300 500 800 1200 2000];
nDropA  = 200;
nfA  = zeros(nGeo, numel(area_ax));
focA = zeros(nGeo, numel(area_ax));
for ia = 1:numel(area_ax)
    sL = area_ax(ia);
    d3v = [];  s2v = [];
    for t = 1:nDropA
        APt = (rand(L, 1) + 1j * rand(L, 1)) * sL;
        UEt = (rand(1, K) + 1j * rand(1, K)) * sL;           % uniform users
        Aw  = wrap_images(APt, sL, wrap_fix);
        for l = 1:L
            [dH, wI] = min(abs(Aw(l, :).' - UEt), [], 1);
            dd = max(sqrt(hDiff^2 + dH.^2), 10);
            ph = angle(UEt - Aw(l, wI));
            ss = sin(ph) .* cos(asin(hDiff ./ dd));
            d3v = [d3v, dd];                                   %#ok<AGROW>
            s2v = [s2v, ss.^2];                                %#ok<AGROW>
        end
    end
    for a = 1:nGeo
        Ra = 2 * A_geo{a}.aperture^2 / lambda;
        nfA(a, ia)  = 100 * mean(d3v < Ra);
        focA(a, ia) = 100 * mean(d3v < Ra * (1 - s2v) / 10.24);
    end
end
figure('Name', 'SPA-3 NF and focusing fraction vs area', 'Position', [40 40 1100 420]);
subplot(1, 2, 1); hold on; box on; grid on;
for a = 1:nGeo
    plot(area_ax, nfA(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
end
xlabel('Square side [m]'); ylabel('Near-field links [%]  (d < 2D^2/\lambda)');
title('Hybrid regime: NF link fraction vs deployment size'); legend('Location', 'southwest');
subplot(1, 2, 2); hold on; box on; grid on;
for a = 1:nGeo
    plot(area_ax, focA(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
end
xlabel('Square side [m]'); ylabel('Focusing-capable links [%]');
title('d < R cos^2\theta / 10.24 (Zhou Prop. 2): range focusing possible'); legend('Location', 'northeast');

% ---------------- SPA-4: co-angular (same-angle) correlation vs range ----------------
r1 = 30;  r2_ax = linspace(10, 300, 600);
figure('Name', 'SPA-4 Co-angular correlation vs range', 'Position', [40 40 760 480]);
hold on; box on; grid on;
for a = 1:nGeo
    a1 = nusw_vec(A_geo{a}.pos, r1, s0, lambda);
    plot(r2_ax, abs(nusw_vec(A_geo{a}.pos, r2_ax, s0 * ones(size(r2_ax)), lambda)' * a1), '-', ...
         'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
end
plot(r2_ax, mu_th * ones(size(r2_ax)), ':k', 'LineWidth', 1.2, 'DisplayName', sprintf('\\mu_{th} = %.1f', mu_th));
xlabel('Range of the second user r_2 [m]'); ylabel('|a(r_1)^H a(r_2)|');
title(sprintf('Two users on the same ray (sin\\theta=%.1f), first user at r_1=%d m', s0, r1));
legend('Location', 'southeast');

% ---------------- SPA-5: ambiguity maps, single AP vs network ----------------
AP_map = [50+50j; 150+40j; 250+60j; 60+150j; 240+150j; 50+250j; 150+260j; 250+240j];
UE_map = 170 + 130j;
st_map = 1;          % [m] the NF beams are ~1 m wide at 100 m, a coarser grid aliases them
xg_map = st_map / 2:st_map:squareLen;
figure('Name', 'SPA-5 Ambiguity maps single AP vs network', 'Position', [40 40 1300 900]);
for a = 1:3
    [G1, Gn, nA] = beam_maps(A_geo{a}.pos, lambda, AP_map, UE_map, squareLen, hDiff, st_map, amb_th);
    subplot(3, 3, a);     imagesc(xg_map, xg_map, G1); axis xy; axis equal; colorbar;
    title(sprintf('%s: AP 1 alone', arr_names{a}));
    subplot(3, 3, 3 + a); imagesc(xg_map, xg_map, Gn); axis xy; axis equal; colorbar;
    title('Network (non-coherent, all 8 APs)');
    subplot(3, 3, 6 + a); imagesc(xg_map, xg_map, nA); axis xy; axis equal; colorbar;
    title(sprintf('# APs with gain > %.1f (ambiguous APs)', amb_th));
end

if arr_study
    % ---------------- aggregation ----------------
    m_pu   = mean(A_pu, 4);                         % 5 x 2 x nArr
    m_amb  = 100 * mean(A_amb, 3);                  % 3 x nArr  [% of area]
    m_nmse = mean(A_nmse, 3);                       % nArr x nSNR
    m_be   = mean(A_ber_e, 4);   m_bp = mean(A_ber_p, 4);
    m_fe   = mean(A_bff_e, 4);   m_fp = mean(A_bff_p, 4);
    m_eta  = mean(A_eta, 3);     m_sz = mean(A_size, 3);
    m_fh   = m_sz + m_eta .* (L - m_sz);            % fronthaul scalars per user (gated cross-AP)
    si0    = snr_ref_idx;

    % ---------------- SPA-6: probability that a pair is unresolvable ----------------
    figure('Name', 'SPA-6 Unresolvable pair probability', 'Position', [40 40 1200 440]);
    ttl6 = {'All user pairs', 'Collinear (co-angular) pairs'};
    for c = 1:2
        subplot(1, 2, c); hold on; box on; grid on;
        bar(1:5, reshape(m_pu(:, c, :), [5, nArr]));
        set(gca, 'XTick', 1:5, 'XTickLabel', {'home AP', 'cluster', 'outer APs', 'all APs', 'rescue'});
        ylabel(sprintf('P(\\mu > %.1f)', mu_th)); title(ttl6{c}); legend(arr_names, 'Location', 'northeast');
    end

    % ---------------- SPA-7: correlation CDFs and ambiguity areas ----------------
    figure('Name', 'SPA-7 Correlation CDFs and ambiguity areas', 'Position', [40 40 1300 420]);
    for v = 1:2
        subplot(1, 3, v); hold on; box on; grid on;
        for a = 1:nArr
            if v == 1, x = sort(vertcat(A_muA{a, :})); else, x = sort(vertcat(A_muN{a, :})); end
            plot(x, (1:numel(x)) / numel(x), '-', 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
        end
        plot([mu_th mu_th], [0 1], ':k', 'DisplayName', '\mu_{th}');
        xlabel('Pair correlation \mu'); ylabel('CDF');
        if v == 1, title('Home AP (single array)'); else, title('Network, all APs (non-coherent)'); end
        legend('Location', 'northwest');
    end
    subplot(1, 3, 3); hold on; box on; grid on;
    bar(1:3, m_amb);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'per-AP range', 'per-AP grating/sidelobe', 'network'});
    ylabel(sprintf('Area with gain > %.1f [%% of square]', amb_th)); title('Ambiguity area'); legend(arr_names);

    % ---------------- SPA-8: estimation NMSE ----------------
    figure('Name', 'SPA-8 Channel estimation NMSE per array', 'Position', [40 40 700 480]);
    hold on; box on; grid on;
    for a = 1:nArr
        semilogy(SNR_dB, m_nmse(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('NMSE (authentic CSI)');
    title('Polar-domain OMP (NF links) + LMMSE (FF links)'); legend('Location', 'southwest');

    % ---------------- SPA-9: BER of every detector per array ----------------
    figure('Name', 'SPA-9 BER per detector and array', 'Position', [40 40 1200 820]);
    for d = 1:4
        subplot(2, 2, d); hold on; box on; grid on;
        for a = 1:nArr
            semilogy(SNR_dB, max(squeeze(m_be(d, a, :)), 1e-6), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, ...
                     'DisplayName', [arr_names{a} ', est.']);
            semilogy(SNR_dB, max(squeeze(m_bp(d, a, :)), 1e-6), [':' mA{a}], 'Color', cA{a}, 'LineWidth', 1.2, ...
                     'DisplayName', [arr_names{a} ', perfect']);
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER'); title(detN{d});
        legend('Location', 'southwest');
    end

    % ---------------- SPA-10: role of the array in the proposed detector ----------------
    figure('Name', 'SPA-10 Cross-AP gain fronthaul and cluster size per array', 'Position', [40 40 1300 420]);
    subplot(1, 3, 1); hold on; box on; grid on;
    for a = 1:nArr
        plot(SNR_dB, squeeze(m_be(3, a, :)) ./ max(squeeze(m_be(4, a, :)), 1e-7), ['-' mA{a}], 'Color', cA{a}, ...
             'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    xlabel('SNR [dB]'); ylabel('BER_{List} / BER_{List+CrossAP}'); title('Gain of cross-AP fusion (est. CSI)');
    legend('Location', 'northwest');
    subplot(1, 3, 2); hold on; box on; grid on;
    for a = 1:nArr
        plot(SNR_dB, m_fh(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    plot(SNR_dB, L * ones(size(SNR_dB)), ':k', 'DisplayName', 'all-AP bound L');
    xlabel('SNR [dB]'); ylabel('Fronthaul scalars per user'); title('|S_k| + \rho (L - |S_k|)'); legend('Location', 'east');
    subplot(1, 3, 3); hold on; box on; grid on;
    for a = 1:nArr
        plot(SNR_dB, m_sz(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    xlabel('SNR [dB]'); ylabel('Mean |S_k|'); title('Serving-cluster size'); legend('Location', 'east');

    % ---------------- SPA-11: information surplus and gate rate ----------------
    figure('Name', 'SPA-11 Information surplus and gate rate per array', 'Position', [40 40 1100 420]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for a = 1:nArr
        x = sort(reshape(A_del(a, :, :), [], 1));
        plot(x, (1:numel(x)) / numel(x), '-', 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    xlabel('\Delta_k (information outside the cluster)'); ylabel('CDF');
    title(sprintf('Information surplus at %d dB', SNR_dB(si0))); legend('Location', 'southeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    for a = 1:nArr
        plot(SNR_dB, m_eta(a, :), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    xlabel('SNR [dB]'); ylabel('\rho'); title('Cross-AP gate firing rate'); legend('Location', 'northeast');

    % ---------------- SPA-12: per-user SE ----------------
    figure('Name', 'SPA-12 Per-user SE per array', 'Position', [40 40 1100 420]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for a = 1:nArr
        x = sort(reshape(A_se(4, a, :, :), [], 1));
        plot(x, (1:numel(x)) / numel(x), '-', 'Color', cA{a}, 'LineWidth', 1.6, 'DisplayName', arr_names{a});
    end
    xlabel('Per-user SE, List+CrossAP [bit/s/Hz]'); ylabel('CDF');
    title(sprintf('Effective-SINR SE at %d dB', SNR_dB(si0))); legend('Location', 'southeast');
    subplot(1, 2, 2); hold on; box on; grid on;
    bar(1:4, reshape(mean(mean(A_se, 4), 3), [4, nArr]));
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Linear', 'Hard-SIC', 'List-SIC', 'List+CrossAP'});
    ylabel('Mean per-user SE [bit/s/Hz]'); title('Every detector, every array'); legend(arr_names);

    % ---------------- SPA-13: near-field vs far-field channel model per array ----------------
    figure('Name', 'SPA-13 Near-field vs far-field model per array', 'Position', [40 40 1200 420]);
    dsel = [1 4];
    for c = 1:2
        subplot(1, 2, c); hold on; box on; grid on;
        for a = 1:nArr
            semilogy(SNR_dB, max(squeeze(m_be(dsel(c), a, :)), 1e-6), ['-' mA{a}], 'Color', cA{a}, 'LineWidth', 1.6, ...
                     'DisplayName', [arr_names{a} ', hybrid NF/FF']);
            semilogy(SNR_dB, max(squeeze(m_fe(dsel(c), a, :)), 1e-6), ['--' mA{a}], 'Color', cA{a}, 'LineWidth', 1.2, ...
                     'DisplayName', [arr_names{a} ', all far-field']);
        end
        set(gca, 'YScale', 'log'); xlabel('SNR [dB]'); ylabel('BER (estimated CSI)'); title(detN{dsel(c)});
        legend('Location', 'southwest');
    end

    % ---------------- SPA-14: sparse vs dense summary at the reference SNR ----------------
%  Neutral on purpose: in a preliminary run (4 setups, 300 m, K = 8, L = 8) the
%  sparse arrays won the single-AP geometry metrics (collinear pairs unresolvable
%  at the home AP: ECA 0.25 vs dense ULA 0.75; focusing-capable links 21% vs
%  0.4%) but NOT the detection BER: the serving clusters already resolved every
%  pair through macro-diversity, and the sparse arrays paid an estimation-NMSE
%  penalty. Read this figure, do not assume its sign.
    figure('Name', 'SPA-14 Sparse vs dense summary', 'Position', [40 40 1300 640]);
    sumv = {100 * mean(A_nf, 2), 100 * mean(A_foc, 2), reshape(m_pu(2, 1, :), [], 1), ...
            m_amb(1, :).', m_nmse(:, si0), squeeze(m_be(4, :, si0)).'};
    sumt = {'NF links [%]', 'Focusing-capable links [%]', 'P(unresolvable) in cluster', ...
            'Per-AP range ambiguity [% area]', sprintf('NMSE at %d dB', SNR_dB(si0)), ...
            sprintf('BER List+CrossAP at %d dB', SNR_dB(si0))};
    for v = 1:6
        subplot(2, 3, v); hold on; box on; grid on;
        bar(1:nArr, sumv{v});
        set(gca, 'XTick', 1:nArr, 'XTickLabel', arr_list); title(sumt{v});
    end

    % ---------------- console table ----------------
    fprintf('\n=== [CODE C] Sparse-array comparison at %d dB (%d setups) ===\n', SNR_dB(si0), nSetups_arr);
    fprintf('%-8s %5s %6s %6s %8s %8s %8s %8s %9s %9s %9s %9s %6s %6s %7s\n', 'array', 'N', 'NF%', 'foc%', ...
            'Pu-home', 'Pu-clu', 'Pu-out', 'rescue', 'ambR%', 'ambG%', 'NMSE', 'BER-XAP', 'rho', 'FH', 't[s]');
    for a = 1:nArr
        fprintf('%-8s %5d %6.1f %6.1f %8.3f %8.3f %8.3f %8.3f %9.2f %9.2f %9.2e %9.2e %6.3f %6.2f %7.1f\n', ...
                arr_list{a}, A_geo{a}.N, 100 * mean(A_nf(a, :)), 100 * mean(A_foc(a, :)), m_pu(1, 1, a), ...
                m_pu(2, 1, a), m_pu(3, 1, a), m_pu(5, 1, a), m_amb(1, a), m_amb(2, a), m_nmse(a, si0), ...
                m_be(4, a, si0), m_eta(a, si0), m_fh(a, si0), mean(A_time(a, :)));
    end
    save('codeC_sparse_array_results.mat', 'arr_list', 'arr_names', 'SNR_dB', 'A_nf', 'A_foc', 'A_pu', ...
         'A_amb', 'A_nmse', 'A_ber_e', 'A_ber_p', 'A_bff_e', 'A_bff_p', 'A_eta', 'A_size', 'A_kap', ...
         'A_del', 'A_se', 'A_time', 'nfA', 'focA', 'area_ax', 'mu_th', 'amb_th');
end

%% ====================================================================
%  CHANNEL ACQUISITION  (single entry point for every estimation block)
%% ====================================================================
function [Hhat, C, Hpre, Cpre] = csi_acquire(H, R, NFm, hdet, beta, pilot, tau, p, ...
                                             N, L, K, nReal, genie, nfOvr, eps_NF, PD)
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

function [b, mg] = refine_atom(r, s0, q0, PD)
% Gridless local search around a grid atom, run in two parameterisations,
% (s, q) and (s, kappa) with kappa = (1 - s^2) q, keeping the better fit.
% (s, kappa) decouples angle and range off end-fire; (s, q) is better
% conditioned near end-fire. s wraps around [-1, 1): with half-wavelength
% spacing sin(theta) = -1 and +1 are the same end-fire direction, and a
% clamped search started at the alias s = -1 could never reach s = +0.99.
    [b1, m1] = refine_param(r, s0, q0, PD, false);
    [b2, m2] = refine_param(r, s0, q0, PD, true);
    if m1 >= m2
        b = b1;
        mg = m1;
    else
        b = b2;
        mg = m2;
    end
end

function [b, mg] = refine_param(r, s0, q0, PD, useKappa)
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
function G = array_geometry(type, N_ref, lambda, eta, M, Nc, Lc, Lc_sa)
% Element positions (metres, centred on the array midpoint) of the AP array.
%  'ula'    dense ULA, N_ref elements, spacing lambda/2
%  'usa'    uniform sparse array, N_ref elements, spacing eta*lambda/2
%  'eca'    extended coprime array (Zhou et al.) with Lc basic coprime arrays
%  'eca_sa' the same ECA with Lc_sa basic arrays (about the dense aperture)
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

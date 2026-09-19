% ==========================================================================
%  Cell-Free Massive MIMO - Uplink: FF-Only vs Hybrid NF/FF vs Centralized
%  + CHANNEL-NORM and INFORMATION-RATE clustering (hybrid NF/FF case)
%  + LIST-BASED (SAC) symbol-level detection
%  + NF-AWARE CONTAMINATION-PENALISED CLUSTERING  [NEW, purely additive]
%
%  Author:     Nayyab Haider, PUC-Rio
%  Supervisor: Prof. Rodrigo C. de Lamare
%
%  ===================== NOTES ON THIS VERSION ============================
%  Figures 1 to 8 are RETAINED unchanged. Figure 9 is NEW: sum-rate of all
%  four detection techniques, the sum-rate counterpart of Figure 5.
%
%  The contamination-penalised / separability clustering study has been
%  REMOVED at the author's request: the chi separability matrix, the Psi
%  score, the BSR-fixed and contamination rules, the clustering quality
%  metrics, the cluster-size sweep, and figures C1 to C4. Clustering is
%  back to the two original rules, channel norm (D_CN) and the original
%  information-rate score (D_BSR).
%
%  ================== SETTINGS RETAINED FROM THE STUDY ====================
%  (S1) contamination_mode is still present and still TRUE, so
%       tau_p = tau_sym = tau_fig2 = K/2 with pilots aligned to the
%       collinear user pairs. This is NOT the original operating point
%       (tau_p = 7, tau_sym = K). It is kept so the results you have
%       already generated reproduce. Set contamination_mode = false to
%       return to the original pilot lengths; every figure will move.
%  (S2) nf_csi_consistent = true. NF pairs have Hhat overwritten with the
%       exact geometry channel h_det, so their error covariance is set to
%       the NUSW floor. Claiming an exact estimate while feeding the
%       combiner a large error covariance makes it over-whiten. Set false
%       to reproduce the original inconsistent path.
%  (S3) snr_ref_idx uses nearest-grid-point lookup, not exact equality,
%       so it cannot silently return empty when 10 dB is off-grid.
%
%  ===================== CROSS-AP FIX (V6) ================================
%  Full residual cancellation across all L*N segments. Under full
%  cooperation Cross-AP reduces exactly to List-SIC, which is the correct
%  reduction property, not a defect.
%  =========================================================================
%
%  REFERENCES:
%  [R1] Lu et al., IEEE COMST, vol.26, no.4, 2024.
%  [R2] Demir, Bjornson, Sanguinetti, FnT, vol.14, 2021.
%  [R3] Bjornson, Sanguinetti, IEEE TCom, vol.68, no.7, 2020.
%  [R8] 3GPP TR 38.901 V17.0.0, 2022.
%  [R9] Ngo et al., IEEE TWC, vol.16, no.3, 2017.
%  [R14] Mashdour et al. (BSR clustering), IEEE WCL, vol.13, no.7, 2024.
%  [R16] Ssettumba et al., list-based detector + AP selection, ISWCS 2022.
%  [R18] Marzetta, IEEE TWC, vol.9, no.11, 2010 (coherent contamination).
% ==========================================================================

clear;
close all;
clc;

%% ====================================================================
%  PARAMETERS
%% ====================================================================
% [OPERATING POINT] near-field detection vs interference-limited clustering.
%  false (default): the near-field study -- large per-AP aperture (N=64), so
%    the array separates users trivially, the system is NOISE-limited, and
%    CN/IR/GBSE clustering coincide (greedy is optimal). This is the regime
%    for the near-field/Cross-AP detection results.
%  true: the CLUSTERING study -- small per-AP arrays (N=4), more users and a
%    denser deployment, so local combiners are INTERFERENCE-limited and the
%    choice of serving APs (CN vs IR vs GBSE) actually matters. This sacrifices
%    the near-field aperture (d_Ray shrinks -> effectively far-field), matching
%    de Lamare's cell-free clustering setup. Run the paper as TWO operating
%    points: near-field detection (false) and clustering gains (true).
interf_preset = false;
if interf_preset
    L = 16;  K = 16;  N = 4;   squareLen_preset = 100;
else
    L = 8;   K = 8;   N = 64;  squareLen_preset = 300;
end
N_BS = L * N;

fc     = 3e9;
c0     = 3e8;
lambda = c0 / fc;
d_ant  = lambda / 2;
d_H    = 0.5;
fc_GHz = fc / 1e9;

ASD_varphi = deg2rad(15);
ASD_theta  = deg2rad(15);

tau_c       = 200;
prelog_pcsi = 1;

SNR_dB  = 0:5:30;
SNR_lin = 10.^(SNR_dB / 10);
nSNR    = length(SNR_dB);
LN      = L * N;

squareLen = squareLen_preset;
hDiff     = 10;
hBS       = 25;
hUT       = 1.5;

D_ap  = (N - 1) * d_ant;
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

squareLen_sweep = [100 200 300 400 500];
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
nf_csi_consistent  = true;   % (S2) NF pairs: exact channel -> NUSW-floor error covariance.
                              % Set false to reproduce the original inconsistent path.

% [GBSE] Graph-Based Steepest-Ascent uplink clustering (adapts Tesolin &
%  de Lamare, IEEE WCL 2026 / arXiv 2607.04074, to the uplink; see the
%  gbse_cluster function). This is the expensive study: it runs the
%  detector-rate objective over the serving-state graph. Throttle via
%  gbse_nEval / gbse_maxMoves, or set gbse_study=false to skip entirely.
gbse_study    = true;   % add the proposed GBSE clustering (method 3)
gbse_nEval    = 2;      % channel realizations per objective evaluation (raise for accuracy)
gbse_maxMoves = 6;      % max steepest-ascent swap moves
gbse_sizeSweep = 1:6;   % cluster-size budgets for the CN/IR/GBSE sweep study
                        %  (the "does clustering ever matter?" diagnostic --
                        %   greedy CN is near-optimal at size 1, so any gain
                        %   only appears at larger clusters / interference)

if contamination_mode
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
fprintf('d_Ray(CF)=%.2fm  d_Ray(BS)=%.0fm  squareLen=%dm\n', d_Ray, d_Ray_BS, squareLen);
fprintf('contamination_mode=%d  pilot_mode=%s  clu_div=%d\n', ...
        contamination_mode, pilot_mode, clu_div);
fprintf('tau_p=%d  tau_sym=%d  tau_fig2=%d  prelog=%.3f\n', ...
        tau_p, tau_sym, tau_fig2, prelog);
fprintf('List: nSym=%d M=%d d_th=%.2f\n', nSym, M_lst, d_th);
fprintf('Reference SNR for area/cluster-size sweeps and CDF: %d dB (index %d)\n', ...
        snr_ref_dB, snr_ref_idx);
fprintf('%d setups x %d MC x %d SNR pts\n', nSetups, nReal, nSNR);
if ~contamination_mode
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

% [GBSE] proposed graph-based steepest-ascent clustering (method 3)
SR_GBSE_lmmse_acc  = zeros(nSetups, nSNR);   BER_GBSE_lmmse_acc = zeros(nSetups, nSNR);
SR_GBSE_cnsic_acc  = zeros(nSetups, nSNR);   BER_GBSE_cnsic_acc = zeros(nSetups, nSNR);
SR_GBSE_lmmse_pf_acc = zeros(nSetups, nSNR); SR_GBSE_cnsic_pf_acc = zeros(nSetups, nSNR);
loadGBSE_acc = zeros(nSetups, nSNR);         % avg GBSE cluster size (fronthaul load)
gbseConv_acc = zeros(nSetups, gbse_maxMoves + 1);  % convergence trajectory @ ref SNR
% [GBSE] cluster-size sweep: sum-rate of CN/IR/GBSE vs a fixed cluster budget,
% evaluated with the coupled objective at the reference SNR (estimated CSI).
srSweep_cn_acc = zeros(nSetups, numel(gbse_sizeSweep));
srSweep_ir_acc = zeros(nSetups, numel(gbse_sizeSweep));
srSweep_gb_acc = zeros(nSetups, numel(gbse_sizeSweep));

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

% [CSI STUDY] perfect vs proper-estimated CSI for the 4 detectors (rate cluster).
%  - "estimated": far-field links use the standard LMMSE estimate; near-field
%    links use the near-field-aware (NUSW, rank-1 spherical) LMMSE estimate that
%    the code already forms via R2 -- i.e. NO deterministic override. This is the
%    genuinely estimated hybrid NF/FF case (see channel_estimation_note.pdf).
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

    wr   = repmat([-squareLen 0 squareLen], [3 1]);
    wrapL = wr(:)' + 1j * (wr(:)');
    APwrap = repmat(APpos, [1 9]) + repmat(wrapL, [L 1]);

    % ------------------------------------------------------------------
    %  [NEW] PILOT ASSIGNMENT
    %  'legacy'  : the original mod() assignment, old figures unchanged.
    %  'aligned' : the two members of a collinear pair share a pilot, so
    %              near field range resolution has something to resolve.
    %              FAVOURABLE. Report it as such.  (see warning W3)
    %  'random'  : neutral control.
    % ------------------------------------------------------------------
    if strcmp(pilot_mode, 'aligned')
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
    if contamination_mode
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

            firstRow = zeros(N, 1);
            firstRow(1) = 1;
            for m = 2:N
                dist = d_H * (m - 1);
                firstRow(m) = exp(1j * 2 * pi * dist * sin_eff) ...
                    * exp(-ASD_varphi^2 / 2 * (2 * pi * dist * cos(phi_lk) * cos(theta))^2) ...
                    * exp(-ASD_theta^2 / 2 * (2 * pi * dist * sin(theta))^2);
            end
            R_FF = beta(l, k) * toeplitz(firstRow);
            R1(:, :, l, k) = R_FF;

            if d3D < d_Ray
                NF_mask(l, k) = true;
                nNF_this = nNF_this + 1;
                r_n = sqrt(d3D^2 + (n_idx * d_ant).^2 - 2 * d3D * (n_idx * d_ant) * sin_eff);
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
    % ================================================================
    %  [GBSE] Compute the proposed graph-based steepest-ascent cluster ONCE
    %  per setup, at the REFERENCE operating point (interference-limited,
    %  where the coupled objective is most informative), then reuse it at
    %  every SNR -- clustering is an SNR-agnostic serving-state assignment,
    %  exactly like D_CN. Done BEFORE the SNR loop so the inner loop stays a
    %  plain 1:nSNR range (required for parfor sliced-variable indexing).
    %  Adapts the Hamming-graph steepest-ascent search of Tesolin & de Lamare
    %  (IEEE WCL 2026 / arXiv 2607.04074) to the uplink; see gbse_cluster.
    if gbse_study
        p_g = SNR_lin(snr_ref_idx);
        Hc_g = randn(LN, nReal, K) + 1j * randn(LN, nReal, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hc_g(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                else
                    Hc_g(idx, :, k) = sqrt(0.5) * Rs2(:, :, l, k) * Hc_g(idx, :, k);
                end
            end
        end
        Npc_g = sqrt(0.5) * (randn(N, nReal, L, tau_p) + 1j * randn(N, nReal, L, tau_p));
        Hhc_g = zeros(LN, nReal, K);
        Cc_g = zeros(N, N, L, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for t = 1:tau_p
                ue_t = find(pilotIndex == t)';
                yp = sqrt(p_g) * tau_p * sum(Hc_g(idx, :, ue_t), 3) + sqrt(tau_p) * Npc_g(:, :, l, t);
                Psi_t = p_g * tau_p * sum(R2(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R2(:, :, l, k) / Psi_t;
                    Hhc_g(idx, :, k) = sqrt(p_g) * RPsi * yp;
                    Cc_g(:, :, l, k) = R2(:, :, l, k) - p_g * tau_p * RPsi * R2(:, :, l, k);
                end
            end
        end
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hhc_g(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                    if nf_csi_consistent
                        Cc_g(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                    end
                end
            end
        end
        gbCand = false(L, K);
        for k = 1:K
            szk = max(nnz(D_CN(:, k)), 1);
            nc  = min(max(2 * szk, szk + 2), L);
            [~, ordc] = sort(metric_CN(:, k), 'descend');
            gbCand(ordc(1:nc), k) = true;
        end
        nEvalG = min(gbse_nEval, nReal);
        objG = @(Dm) gbse_obj(Dm, Hhc_g, Hc_g, Cc_g, p_g, prelog, eta_FH, N, L, K, nEvalG);
        [D_GBSE_setup, gbTraj] = gbse_cluster(D_CN, gbCand, objG, gbse_maxMoves);
        convRow = gbTraj(end) * ones(1, gbse_maxMoves + 1);
        nT = min(numel(gbTraj), gbse_maxMoves + 1);
        convRow(1:nT) = gbTraj(1:nT);
        gbseConv_acc(ns, :) = convRow;

        % ---- [GBSE] cluster-size sweep (diagnostic) ----------------------
        %  For each fixed cluster budget M, build the CN (top-M channel norm),
        %  IR (top-M rate score) and GBSE (steepest ascent from CN, candidates
        %  = top-2M) clusters and score each with the SAME coupled objective at
        %  the reference SNR. Greedy CN is near-optimal at M=1 (no choice), so
        %  any GBSE gain can only appear at larger M / interference-limited
        %  operation. Flat curves here mean the regime does not reward
        %  clustering -- that is a result, not a bug.
        SRlk_g = zeros(L, K);
        for l = 1:L
            for k = 1:K
                gmm = max(real(trace(R2(:, :, l, k))) - real(trace(Cc_g(:, :, l, k))), 0);
                er  = real(trace(Cc_g(:, :, l, k)));
                SRlk_g(l, k) = log2(1 + p_g * gmm / (p_g * er + 1));
            end
        end
        for mi = 1:numel(gbse_sizeSweep)
            Msz = min(gbse_sizeSweep(mi), L);
            Dcn_M = false(L, K);  Dir_M = false(L, K);  candM = false(L, K);
            for k = 1:K
                [~, oc] = sort(metric_CN(:, k), 'descend');
                Dcn_M(oc(1:Msz), k) = true;
                ncc = min(max(2 * Msz, Msz + 2), L);
                candM(oc(1:ncc), k) = true;
                [~, oi] = sort(SRlk_g(:, k), 'descend');
                Dir_M(oi(1:Msz), k) = true;
            end
            Dgb_M = gbse_cluster(Dcn_M, candM, objG, gbse_maxMoves);
            srSweep_cn_acc(ns, mi) = objG(Dcn_M);
            srSweep_ir_acc(ns, mi) = objG(Dir_M);
            srSweep_gb_acc(ns, mi) = objG(Dgb_M);
        end
    else
        D_GBSE_setup = false(L, K);
    end

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

            Np  = sqrt(0.5) * (randn(N, nReal, L, tau_p) + 1j * randn(N, nReal, L, tau_p));
            Hhat = zeros(LN, nReal, K);
            C = zeros(N, N, L, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for t = 1:tau_p
                    ue_t = find(pilotIndex == t)';
                    yp = sqrt(p) * tau_p * sum(H(idx, :, ue_t), 3) + sqrt(tau_p) * Np(:, :, l, t);
                    Psi_t = p * tau_p * sum(R_use(:, :, l, ue_t), 4) + eye(N);
                    for k = ue_t
                        RPsi = R_use(:, :, l, k) / Psi_t;
                        Hhat(idx, :, k) = sqrt(p) * RPsi * yp;
                        C(:, :, l, k) = R_use(:, :, l, k) - p * tau_p * RPsi * R_use(:, :, l, k);
                    end
                end
            end

            if caseID == 2
                for l = 1:L
                    idx = (l - 1) * N + 1:l * N;
                    for k = 1:K
                        if NF_mask(l, k)
                            Hhat(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                            if nf_csi_consistent
                                C(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                            end
                        end
                    end
                end
            end

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
            H = []; Hhat = []; C = []; Np = [];
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
        Npc = sqrt(0.5) * (randn(N, nReal, L, tau_p) + 1j * randn(N, nReal, L, tau_p));
        Hhc = zeros(LN, nReal, K);
        Cc = zeros(N, N, L, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for t = 1:tau_p
                ue_t = find(pilotIndex == t)';
                yp = sqrt(p) * tau_p * sum(Hc(idx, :, ue_t), 3) + sqrt(tau_p) * Npc(:, :, l, t);
                Psi_t = p * tau_p * sum(R2(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R2(:, :, l, k) / Psi_t;
                    Hhc(idx, :, k) = sqrt(p) * RPsi * yp;
                    Cc(:, :, l, k) = R2(:, :, l, k) - p * tau_p * RPsi * R2(:, :, l, k);
                end
            end
        end
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hhc(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                    if nf_csi_consistent
                        Cc(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                    end
                end
            end
        end

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

        % [GBSE] use the per-setup cluster computed before the SNR loop
        %  (cardinality-preserving swaps => fronthaul load matched to CN/IR;
        %  coupled det_local objective => sees the interference CN/IR ignore).
        if gbse_study
            D_GBSE = D_GBSE_setup;
            loadGBSE_acc(ns, si) = mean(sum(D_GBSE, 1));
        else
            D_GBSE = D_BSR;
        end

        for method = 1:3
            if method == 1
                Dcl = D_CN;
            elseif method == 3
                Dcl = D_GBSE;
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
            elseif method == 3
                % [GBSE] proposed graph-based steepest-ascent cluster
                SR_GBSE_lmmse_acc(ns, si) = sum(mean(SEl, 2));
                BER_GBSE_lmmse_acc(ns, si) = mean(BEl(:));
                SR_GBSE_cnsic_acc(ns, si) = sum(mean(SEs, 2));
                BER_GBSE_cnsic_acc(ns, si) = mean(BEs(:));
                SR_GBSE_lmmse_pf_acc(ns, si) = sum(mean(SEl_pf, 2));
                SR_GBSE_cnsic_pf_acc(ns, si) = sum(mean(SEs_pf, 2));
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
        Nps = sqrt(0.5) * (randn(N, nReal, L, tau_sym) + 1j * randn(N, nReal, L, tau_sym));
        Hhs = zeros(LN, nReal, K);
        Cs = zeros(N, N, L, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for t = 1:tau_sym
                ue_t = find(pilotSym == t)';
                yp = sqrt(p) * tau_sym * sum(Hcs(idx, :, ue_t), 3) + sqrt(tau_sym) * Nps(:, :, l, t);
                Psi_t = p * tau_sym * sum(R2(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R2(:, :, l, k) / Psi_t;
                    Hhs(idx, :, k) = sqrt(p) * RPsi * yp;
                    Cs(:, :, l, k) = R2(:, :, l, k) - p * tau_sym * RPsi * R2(:, :, l, k);
                end
            end
        end
        % [CSI STUDY] keep the proper NF-aware LMMSE estimate BEFORE the
        % deterministic near-field override below. For NF links R2 is the rank-1
        % spherical (NUSW) correlation, so Hhs_est is a genuine near-field estimate
        % (this is the "estimated CSI" used by the two new summary figures).
        Hhs_est = Hhs;  Cs_est = Cs;
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hhs(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                    if nf_csi_consistent
                        Cs(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                    end
                end
            end
        end

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
                [q1, q2, q3, q4] = ber_case(Hm, Hm, Cs0, D_BSR_s, p, N, L, K, nSym, M_lst, d_th, maxBr, modOrder);
                brt_pf = brt_pf + [q1; q2; q3; q4];
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
        % [ESTIMATOR NMSE] channel-estimation quality tr(C)/tr(R) over links (hybrid)
        e_num = 0;  e_den = 0;
        for l = 1:L
            for k = 1:K
                e_num = e_num + real(trace(Cs_est(:, :, l, k)));
                e_den = e_den + real(trace(R2(:, :, l, k)));
            end
        end
        nmse_est_acc(ns, si) = e_num / max(e_den, eps);
        Hcs = []; Hhs = []; Cs = []; Nps = []; Hhs_est = []; Cs_est = [];

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
            Npf = sqrt(0.5) * (randn(N, nReal, L, tau_sym) + 1j * randn(N, nReal, L, tau_sym));
            Hhf = zeros(LN, nReal, K);  Cf = zeros(N, N, L, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for t = 1:tau_sym
                    ue_t = find(pilotSym == t)';
                    yp = sqrt(p) * tau_sym * sum(Hcf(idx, :, ue_t), 3) + sqrt(tau_sym) * Npf(:, :, l, t);
                    Psi_t = p * tau_sym * sum(R1(:, :, l, ue_t), 4) + eye(N);
                    for k = ue_t
                        RPsi = R1(:, :, l, k) / Psi_t;
                        Hhf(idx, :, k) = sqrt(p) * RPsi * yp;
                        Cf(:, :, l, k) = R1(:, :, l, k) - p * tau_sym * RPsi * R1(:, :, l, k);
                    end
                end
            end
            f_num = 0;  f_den = 0;
            for l = 1:L
                for k = 1:K
                    f_num = f_num + real(trace(Cf(:, :, l, k)));
                    f_den = f_den + real(trace(R1(:, :, l, k)));
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
            Hcf = []; Hhf = []; Cf = []; Npf = [];
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
        Np2 = sqrt(0.5) * (randn(N, nReal, L, tau_fig2) + 1j * randn(N, nReal, L, tau_fig2));
        Hh2 = zeros(LN, nReal, K);
        Cc2 = zeros(N, N, L, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for t = 1:tau_fig2
                ue_t = find(pilotF2 == t)';
                yp = sqrt(p) * tau_fig2 * sum(Hc2(idx, :, ue_t), 3) + sqrt(tau_fig2) * Np2(:, :, l, t);
                Psi_t = p * tau_fig2 * sum(R2(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R2(:, :, l, k) / Psi_t;
                    Hh2(idx, :, k) = sqrt(p) * RPsi * yp;
                    Cc2(:, :, l, k) = R2(:, :, l, k) - p * tau_fig2 * RPsi * R2(:, :, l, k);
                end
            end
        end
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for k = 1:K
                if NF_mask(l, k)
                    Hh2(idx, :, k) = repmat(h_det_all(:, l, k), [1, nReal]);
                    if nf_csi_consistent
                        Cc2(:, :, l, k) = eps_NF * beta(l, k) * eye(N);
                    end
                end
            end
        end
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
        NpF2 = sqrt(0.5) * (randn(N, nReal, L, tau_fig2) + 1j * randn(N, nReal, L, tau_fig2));
        HhF2 = zeros(LN, nReal, K);
        CcF2 = zeros(N, N, L, K);
        for l = 1:L
            idx = (l - 1) * N + 1:l * N;
            for t = 1:tau_fig2
                ue_t = find(pilotF2 == t)';
                yp = sqrt(p) * tau_fig2 * sum(HcF2(idx, :, ue_t), 3) + sqrt(tau_fig2) * NpF2(:, :, l, t);
                Psi_t = p * tau_fig2 * sum(R1(:, :, l, ue_t), 4) + eye(N);
                for k = ue_t
                    RPsi = R1(:, :, l, k) / Psi_t;
                    HhF2(idx, :, k) = sqrt(p) * RPsi * yp;
                    CcF2(:, :, l, k) = R1(:, :, l, k) - p * tau_fig2 * RPsi * R1(:, :, l, k);
                end
            end
        end

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
        Hc2 = []; Hh2 = []; Cc2 = []; Np2 = []; HcF2 = []; HhF2 = []; CcF2 = []; NpF2 = []; Hc_BS2 = []; Hh_BS2 = []; C_BS2 = []; Np_BS2 = [];

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

% [GBSE] proposed clustering, averaged over setups
SR_GBSE_lmmse  = mean(SR_GBSE_lmmse_acc, 1);   SR_GBSE_cnsic  = mean(SR_GBSE_cnsic_acc, 1);
BER_GBSE_lmmse = mean(BER_GBSE_lmmse_acc, 1);  BER_GBSE_cnsic = mean(BER_GBSE_cnsic_acc, 1);
SR_GBSE_lmmse_pf = mean(SR_GBSE_lmmse_pf_acc, 1);  SR_GBSE_cnsic_pf = mean(SR_GBSE_cnsic_pf_acc, 1);
loadGBSE = mean(loadGBSE_acc, 1);
gbseConv = mean(gbseConv_acc, 1);   % mean convergence trajectory @ ref SNR
srSweep_cn = mean(srSweep_cn_acc, 1);   % cluster-size sweep, CN/IR/GBSE
srSweep_ir = mean(srSweep_ir_acc, 1);
srSweep_gb = mean(srSweep_gb_acc, 1);

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

        wr_sw = repmat([-sLen 0 sLen], [3 1]);
        wL_sw = wr_sw(:)' + 1j * (wr_sw(:)');
        APw_sw = repmat(APsw, [1 9]) + repmat(wL_sw, [L 1]);
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
                fR = zeros(N, 1);
                fR(1) = 1;
                for m = 2:N
                    ds = d_H * (m - 1);
                    fR(m) = exp(1j * 2 * pi * ds * se) * exp(-ASD_varphi^2 / 2 * (2 * pi * ds * cos(ph_lk) * cos(th))^2) ...
                        * exp(-ASD_theta^2 / 2 * (2 * pi * ds * sin(th))^2);
                end
                RFF_lk = beta_sw(l, k) * toeplitz(fR);
                R1_sw(:, :, l, k) = RFF_lk;
                if d3D < d_Ray
                    NF_sw(l, k) = true;
                    nNFs = nNFs + 1;
                    rn_s = sqrt(d3D^2 + (n_idx * d_ant).^2 - 2 * d3D * (n_idx * d_ant) * se);
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
            Np_sw = sqrt(0.5) * (randn(N, nReal, L, tau_p) + 1j * randn(N, nReal, L, tau_p));
            Hh = zeros(LN, nReal, K);
            Csw = zeros(N, N, L, K);
            for l = 1:L
                idx = (l - 1) * N + 1:l * N;
                for t = 1:tau_p
                    ut = find(piIdx == t)';
                    yp = sqrt(p_ref) * tau_p * sum(H_sw(idx, :, ut), 3) + sqrt(tau_p) * Np_sw(:, :, l, t);
                    Pt = p_ref * tau_p * sum(Ru(:, :, l, ut), 4) + eye(N);
                    for k = ut
                        RP = Ru(:, :, l, k) / Pt;
                        Hh(idx, :, k) = sqrt(p_ref) * RP * yp;
                        Csw(:, :, l, k) = Ru(:, :, l, k) - p_ref * tau_p * RP * Ru(:, :, l, k);
                    end
                end
            end
            if cID == 2
                for l = 1:L
                    idx = (l - 1) * N + 1:l * N;
                    for k = 1:K
                        if NF_sw(l, k)
                            Hh(idx, :, k) = repmat(hd_sw(:, l, k), [1, nReal]);
                        end
                    end
                end
            end
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
            clear H_sw Hh Csw Np_sw;
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
cGB = [0.60 0.00 0.70];   % [GBSE] proposed clustering

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
plot(SNR_dB, SR_GBSE_lmmse_pf, '-p', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 2, 'DisplayName', 'Hybrid GBSE-clustering (proposed): L-MMSE');
plot(SNR_dB, SR_GBSE_cnsic_pf, '-h', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 2, 'DisplayName', 'Hybrid GBSE-clustering (proposed): SIC');
xlabel('SNR [dB]', 'FontSize', 13);
ylabel('Sum-rate [bps/Hz]', 'FontSize', 13);
title(sprintf(['Sum-rate, PERFECT CSI: Centralized BS (far-field) vs Cell-free (far-field & hybrid NF/FF) + Clustering (CN / IR / GBSE)\n'...
               'Channel: hybrid near-field/far-field; GBSE at matched fronthaul load; d_{Ray}(CF)=%.0fm  N_{BS}=%d  squareLen=%dm  L=%d  K=%d  \\eta_{FH}=%.2f'], ...
              d_Ray, N_BS, squareLen, L, K, eta_FH), 'FontSize', 9);
legend('Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
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
    title(sprintf('BER of all detectors: estimated vs perfect CSI (hybrid NF/FF)\nNF links: near-field (NUSW) LMMSE;  FF links: LMMSE;  \\tau_{sym}=%d', tau_sym), 'FontSize', 11);
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
    semilogy(SNR_dB, max(nmse_est, 1e-6), '-o', 'Color', cNF, 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', 'Near-field (NUSW) LMMSE'); hold on; box on; grid on;
    if ff_study
        semilogy(SNR_dB, max(nmse_ff, 1e-6), '-s', 'Color', cFF, 'LineWidth', 2, 'MarkerSize', 7, ...
            'DisplayName', 'Far-field LMMSE');
    end
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Estimator NMSE = tr(C)/tr(R)', 'FontSize', 13);
    title(sprintf(['Channel-estimation quality vs SNR, hybrid model (\\tau_{sym} pilots, LMMSE)\n'...
                   'NF lower NMSE because the NUSW correlation R is rank-1 (few DoF to estimate), not a fairness gain']), 'FontSize', 10);
    legend('Location', 'southwest', 'FontSize', 9);

    %% ================================================================
    %  [FIG 20 companion] PER-DETECTOR SYMBOL NMSE under ESTIMATED CSI
    %  The estimator-NMSE figure above (tr(C)/tr(R)) is channel-estimation
    %  quality and does NOT depend on the detector, so it carries no
    %  per-detector (e.g. Cross-AP) curve. The detector-dependent quantity is
    %  the post-detection SYMBOL NMSE; here it is plotted under ESTIMATED CSI
    %  for all four detectors, including the proposed Cross-AP List-SIC.
    %  Channel model: hybrid near-field/far-field, IR (rate) clustering.
    %% ================================================================
    figure('Name', 'SymbolNMSE-Detectors-Estimated', 'Position', [130 30 880 560]);
    hold on; box on; grid on; set(gca, 'YScale', 'log');
    for d = 1:4
        semilogy(SNR_dB, max(NMSE_det_est(d, :), 1e-6), ['-' mk4{d}], 'Color', cD4{d}, ...
            'LineWidth', 2 + 0.3 * (d >= 3), 'MarkerSize', 7, 'DisplayName', nm4{d});
    end
    set(gca, 'YMinorGrid', 'on', 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 13); ylabel('Symbol NMSE', 'FontSize', 13);
    title(sprintf(['Per-detector symbol NMSE, ESTIMATED CSI (hybrid NF/FF, IR clustering)\n'...
                   'includes proposed Cross-AP List-SIC']), 'FontSize', 10);
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

    %% NEW FIGURE F - ALL FOUR detectors incl. Cross-AP: NF vs FF, ESTIMATED CSI
    figure('Name', 'AllDetectors-NF-vs-FF-Estimated', 'Position', [40 40 1220 540]);
    subplot(1, 2, 1); hold on; box on; grid on;
    for d = 1:4
        semilogy(SNR_dB, max(aBERrt_est(d, :), flrC),    ['-' mk4{d}],  'Color', cD4{d}, 'LineWidth', 2 + 0.3 * (d >= 3), 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (near-field)']);
        semilogy(SNR_dB, max(aBERrt_ff_est(d, :), flrC), ['--' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2 + 0.3 * (d >= 3), 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (far-field)']);
    end
    set(gca, 'YScale', 'log', 'YMinorGrid', 'on', 'XTick', SNR_dB); ylim([1e-5 1]);
    xlabel('SNR [dB]', 'FontSize', 12); ylabel('BER (16-QAM)', 'FontSize', 12);
    title('BER (estimated CSI)', 'FontSize', 11); legend('Location', 'southwest', 'FontSize', 7, 'NumColumns', 2);
    subplot(1, 2, 2); hold on; box on; grid on;
    for d = 1:4
        plot(SNR_dB, SE_det_est(d, :),    ['-' mk4{d}],  'Color', cD4{d}, 'LineWidth', 2 + 0.3 * (d >= 3), 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (near-field)']);
        plot(SNR_dB, SE_det_ff_est(d, :), ['--' mk4{d}], 'Color', cD4{d}, 'LineWidth', 2 + 0.3 * (d >= 3), 'MarkerSize', 7, 'DisplayName', [nm4{d} ' (far-field)']);
    end
    set(gca, 'XTick', SNR_dB);
    xlabel('SNR [dB]', 'FontSize', 12); ylabel('Goodput sum-rate [bps/Hz]', 'FontSize', 12);
    title('Sum-rate (estimated CSI)', 'FontSize', 11); legend('Location', 'northwest', 'FontSize', 7, 'NumColumns', 2);
    sgtitle('All detectors (Linear, SIC, List-SIC, Cross-AP List-SIC): near-field vs far-field channel model, ESTIMATED CSI', 'FontSize', 11);
end

%% ====================================================================
%  [GBSE] PROPOSED CLUSTERING STUDY - CN vs IR vs GBSE
%  Channel model: hybrid near-field/far-field. Demonstrates that the
%  graph-based steepest-ascent search (adapted from Tesolin & de Lamare,
%  IEEE WCL 2026 / arXiv 2607.04074) dominates the separable channel-norm
%  (CN) and information-rate (IR) rules at MATCHED fronthaul load, because
%  its inner objective is the coupled post-combiner uplink rate.
%   Panel 1: sum-rate vs SNR (perfect CSI), L-MMSE, CN/IR/GBSE.
%   Panel 2: sum-rate vs SNR (estimated CSI), L-MMSE, CN/IR/GBSE.
%   Panel 3: mean convergence trajectory of GBSE @ reference SNR (monotone).
%   Panel 4: mean cluster size CN/IR/GBSE (equal => fair comparison).
%% ====================================================================
if gbse_study
    loadCN_m  = mean(loadCN_acc(:));
    loadBSR_m = mean(loadBSR_acc(:));
    loadGBSE_m = mean(loadGBSE_acc(:));
    figure('Name', 'GBSE-Clustering-Study', 'Position', [40 40 1300 620]);

    subplot(1, 4, 1); hold on; box on; grid on;
    plot(SNR_dB, SR_CNcl_lmmse_pf, '--o', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'CN-clustering');
    plot(SNR_dB, SR_RTcl_lmmse_pf, '--^', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'IR-clustering');
    plot(SNR_dB, SR_GBSE_lmmse_pf, '-p', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 2, 'DisplayName', 'GBSE (proposed)');
    set(gca, 'XTick', xt); xlabel('SNR [dB]', 'FontSize', 11); ylabel('Sum-rate [bps/Hz]', 'FontSize', 11);
    title(sprintf('Sum-rate (perfect CSI), L-MMSE\nhybrid NF/FF, matched load'), 'FontSize', 9);
    legend('Location', 'northwest', 'FontSize', 8);

    subplot(1, 4, 2); hold on; box on; grid on;
    plot(SNR_dB, SR_CNcl_lmmse, '--o', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'CN-clustering');
    plot(SNR_dB, SR_RTcl_lmmse, '--^', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', 'IR-clustering');
    plot(SNR_dB, SR_GBSE_lmmse, '-p', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 2, 'DisplayName', 'GBSE (proposed)');
    set(gca, 'XTick', xt); xlabel('SNR [dB]', 'FontSize', 11); ylabel('Sum-rate [bps/Hz]', 'FontSize', 11);
    title(sprintf('Sum-rate (estimated CSI), L-MMSE\nhybrid NF/FF, matched load'), 'FontSize', 9);
    legend('Location', 'northwest', 'FontSize', 8);

    subplot(1, 4, 3); hold on; box on; grid on;
    plot(0:numel(gbseConv) - 1, gbseConv, '-p', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 2);
    xlabel('Steepest-ascent move', 'FontSize', 11); ylabel('Objective: uplink sum-rate [bps/Hz]', 'FontSize', 11);
    title(sprintf('GBSE convergence @ %d dB\n(monotone by construction)', snr_ref_dB), 'FontSize', 9);

    subplot(1, 4, 4); hold on; box on; grid on;
    bar([loadCN_m, loadBSR_m, loadGBSE_m], 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'CN', 'IR', 'GBSE'});
    ylabel('Mean cluster size [APs/UE]', 'FontSize', 11);
    title(sprintf('Fronthaul load (equal => fair)\nCN=%.2f IR=%.2f GBSE=%.2f', loadCN_m, loadBSR_m, loadGBSE_m), 'FontSize', 9);

    sgtitle('Proposed GBSE clustering vs CN / IR (hybrid NF/FF, matched fronthaul load)', 'FontSize', 12);

    fprintf('\n--- [GBSE] Clustering sum-rate (L-MMSE) at %d dB ---\n', snr_ref_dB);
    fprintf('   CN=%.3f  IR=%.3f  GBSE=%.3f bps/Hz (perfect CSI)\n', ...
            SR_CNcl_lmmse_pf(snr_ref_idx), SR_RTcl_lmmse_pf(snr_ref_idx), SR_GBSE_lmmse_pf(snr_ref_idx));
    fprintf('   CN=%.3f  IR=%.3f  GBSE=%.3f bps/Hz (estimated CSI)\n', ...
            SR_CNcl_lmmse(snr_ref_idx), SR_RTcl_lmmse(snr_ref_idx), SR_GBSE_lmmse(snr_ref_idx));
    fprintf('   mean cluster size: CN=%.2f IR=%.2f GBSE=%.2f APs/UE\n', loadCN_m, loadBSR_m, loadGBSE_m);

    %% ================================================================
    %  [GBSE] CLUSTER-SIZE SWEEP - the "does clustering matter?" figure.
    %  Sum-rate (coupled objective, estimated CSI, reference SNR) vs the fixed
    %  per-UE cluster budget M, for CN / IR / GBSE. If GBSE separates from CN
    %  anywhere, it is in the mid-range M (enough choice, interference binds).
    %  If the curves overlap for all M, the deployment is not interference-
    %  limited and clustering choice does not matter here -- an honest result.
    %% ================================================================
    figure('Name', 'GBSE-ClusterSize-Sweep', 'Position', [80 60 760 560]);
    hold on; box on; grid on;
    plot(gbse_sizeSweep, srSweep_cn, '--o', 'Color', cCN, 'LineWidth', lw, 'MarkerSize', ms + 1, 'DisplayName', 'CN-clustering');
    plot(gbse_sizeSweep, srSweep_ir, '--^', 'Color', cRT, 'LineWidth', lw, 'MarkerSize', ms + 1, 'DisplayName', 'IR-clustering');
    plot(gbse_sizeSweep, srSweep_gb, '-p', 'Color', cGB, 'LineWidth', lw + 0.6, 'MarkerSize', ms + 3, 'DisplayName', 'GBSE (proposed)');
    xlabel('Cluster-size budget M [APs/UE]', 'FontSize', 12);
    ylabel(sprintf('Uplink sum-rate @ %d dB [bps/Hz]', snr_ref_dB), 'FontSize', 12);
    title(sprintf(['Clustering sum-rate vs fronthaul budget (estimated CSI, hybrid NF/FF)\n'...
                   'gap opens only where choice + interference exist']), 'FontSize', 10);
    legend('Location', 'southeast', 'FontSize', 10);
    set(gca, 'XTick', gbse_sizeSweep);

    fprintf('\n--- [GBSE] cluster-size sweep, uplink sum-rate @ %d dB (estimated CSI) ---\n', snr_ref_dB);
    fprintf('%-6s %10s %10s %10s %10s\n', 'M', 'CN', 'IR', 'GBSE', 'GBSE-CN%');
    for mi = 1:numel(gbse_sizeSweep)
        fprintf('%-6d %10.3f %10.3f %10.3f %9.2f%%\n', gbse_sizeSweep(mi), ...
                srSweep_cn(mi), srSweep_ir(mi), srSweep_gb(mi), ...
                100 * (srSweep_gb(mi) - srSweep_cn(mi)) / max(srSweep_cn(mi), eps));
    end
end

%% ====================================================================
%  [GBSE] Graph-Based Steepest-Ascent uplink clustering.
%  Adapts the Hamming-graph steepest-ascent (SAHC) search of Tesolin & de
%  Lamare (IEEE WCL 2026 / arXiv 2607.04074) from downlink energy efficiency
%  to the UPLINK detection objective. The search space is the serving-state
%  graph; here the neighborhood uses CARDINALITY-PRESERVING SWAP moves (remove
%  one serving AP, add one eligible candidate), so every visited state keeps
%  the same per-user cluster size as the D_CN initializer -- the fronthaul
%  load is held equal to CN/IR by construction, which is what makes the
%  head-to-head comparison fair. objfun returns the coupled uplink sum-rate to
%  MAXIMISE. Monotone ascent over a finite state space => convergence, and the
%  returned cluster is >= the CN initialization by the ascent argument.
function [D, traj] = gbse_cluster(D0, cand, objfun, maxMoves)
    [~, K] = size(D0);
    D = logical(D0);
    best = objfun(D);
    traj = best;
    for mv = 1:maxMoves
        bestD = D; bestVal = best; improved = false;
        for k = 1:K
            served = find(D(:, k));
            cands  = find(cand(:, k) & ~D(:, k));
            for ia = 1:numel(served)
                a = served(ia);
                for ib = 1:numel(cands)
                    b = cands(ib);
                    Dt = D; Dt(a, k) = false; Dt(b, k) = true;
                    v = objfun(Dt);
                    if v > bestVal + 1e-9
                        bestVal = v; bestD = Dt; improved = true;
                    end
                end
            end
        end
        if ~improved
            break;
        end
        D = bestD; best = bestVal; traj(end + 1) = best; %#ok<AGROW>
    end
end

%  [GBSE] coupled inner objective: mean post-combiner uplink sum-rate of the
%  candidate serving matrix D, evaluated with the local L-MMSE detector under
%  ESTIMATED CSI (Hhat estimate, Htrue true channel, C error covariance).
%  Interference between co-served UEs enters det_local's SINR, so unlike the
%  separable CN/IR scores this objective sees the coupling clustering discards.
function val = gbse_obj(D, Hhat, Htrue, C, p, prelog, eta, N, L, K, nEval)
    val = 0;
    for mc = 1:nEval
        Hh = reshape(Hhat(:, mc, :),  [L * N, K]);
        Hm = reshape(Htrue(:, mc, :), [L * N, K]);
        se = det_local(Hh, Hm, C, D, p, prelog, eta, N, L, K);
        val = val + sum(se);
    end
    val = val / max(nEval, 1);
end

%% ====================================================================
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

function [berU, nmseU, bits_pu, eng_pu] = metrics_case(Hhat, Htrue, C, serv, p, N, L, K, nSym, M, dth, maxBr, modOrder)
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
                        if shad || shad_xap
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
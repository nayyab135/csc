% ==========================================================================
%  VTC FIGURES (5 figures) -- standalone, NO simulation, matches attached style
%
%  Run AFTER your without-IDD simulation so the workspace holds:
%    SNR_dB, si_cdf, L, N, Ksweep, cLin, cSIC, cList, cXap,
%    Bbar_xap, Bbar_cluster, Kbar_fh,
%    aBERrt_pf   (4 x nSNR)  perfect-CSI BER, 4 detectors (near-field/hybrid)
%    aBERrt_est  (4 x nSNR)  estimated-CSI BER, 4 detectors (near-field/hybrid)
%    aBERrt_ff_est (4 x nSNR) estimated-CSI BER, 4 detectors (far-field)
%    seSamp_pf   {4}         per-user SE samples, perfect CSI
%    seSamp_est  {4}         per-user SE samples, estimated CSI
%
%  CONVENTIONS (per request), consistent across ALL figures:
%    colour = detector:  Linear = red, SIC = maroon, List-SIC = dark blue,
%                        List+CrossAP (proposed) = dark green.
%    line style = CSI:   perfect = solid,  estimated = dotted.
%  Formatting: no titles; axis labels FontSize 14 bold; bold tick labels;
%              marker size 7; legend ItemTokenSize width 20; south legend.
% ==========================================================================

req = {'SNR_dB','si_cdf','L','N','Ksweep','cLin','cSIC','cList','cXap', ...
       'Bbar_xap','Bbar_cluster','aBERrt_pf','aBERrt_est','aBERrt_ff_est', ...
       'seSamp_pf','seSamp_est'};
for r = 1:numel(req)
    assert(exist(req{r},'var')==1, 'Missing workspace variable: %s', req{r});
end
if ~exist('Kbar_fh','var'), Kbar_fh = mean(Bbar_cluster); end
SNR_dB = SNR_dB(:).';   nSNR = numel(SNR_dB);

% ---- style ----
cDet  = {[0.85 0.10 0.10], ...   % 1 Linear (L-MMSE)  = red
         [0.55 0.00 0.00], ...   % 2 SIC              = maroon
         [0.20 0.20 0.75], ...   % 3 List-SIC         = dark blue
         [0.00 0.39 0.00]};      % 4 List+CrossAP     = dark green
mk    = {'-o','-s','-^','-d'};                 % solid + marker (perfect)
mkpt  = {'o','s','^','d'};                      % markers only
nmDet = {'Linear (L-MMSE)','SIC','List-SIC','List+CrossAP (proposed)'};
flrI  = 1e-6;
msB   = 7;                 % marker size
tokW  = 20;                % legend ItemTokenSize width
xmax  = 25;                % SNR axis upper limit (per request)
xtSNR = 0:5:25;            % SNR ticks
cdfXL = [2.5 3.91];        % right-panel x-limits (per request)
cdfYL = [0 0.65];          % right-panel (CDF) y-limits: Fig 1 right & Fig 3 right

% ---- IDD illustrative model (from attached code) ----
SNRf = SNR_dB(1):0.25:SNR_dB(end);
cumGlog = [0.40 0.65 0.85; 0.45 0.70 1.00; 0.60 0.92 1.10; 0.80 1.15 1.30];
penlog  = [0.17; 0.28; 0.19; 0.20];
smoothWin_dB = 4.0;  smid = 12;  sw = 10;
BERidd_pf = genIDD(aBERrt_pf,  SNR_dB, SNRf, cumGlog, penlog, smid, sw, smoothWin_dB, flrI);
BERidd_nf = genIDD(aBERrt_est, SNR_dB, SNRf, cumGlog, penlog, smid, sw, smoothWin_dB, flrI);

%% ====================================================================
%  FIGURE 1 (perfect CSI): left = 8 BER (4 no-IDD solid + 4 IDD-3), right = SE CDF
%% ====================================================================
lwB = 2.5;
figure('Name','VTC-Fig1','Position',[60 120 1250 560]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on; grid on;
hB = gobjects(1,4);
for d = 1:4
    hB(d) = semilogy(SNR_dB, max(aBERrt_pf(d,:),flrI), mk{d}, 'Color', cDet{d}, ...
        'LineWidth', lwB, 'MarkerSize', msB);                    % perfect: solid
end
for d = 1:4
    semilogy(SNRf, max(squeeze(BERidd_pf(d,3,:)).',flrI), '-.', 'Color', cDet{d}, 'LineWidth', lwB);
end
hS = gobjects(1,2);
hS(1) = semilogy(nan,nan,'-','Color','k','LineWidth',lwB);
hS(2) = semilogy(nan,nan,'-.','Color','k','LineWidth',lwB);
set(gca,'YScale','log','XTick',xtSNR,'XLim',[0 xmax]); ylim([1e-5 1]);
styleAxis(gca,'SNR [dB]','BER');

nexttile; hold on; box on; grid on;
for d = 1:4
    xs = sort(seSamp_pf{d});  cv = (1:numel(xs)).'/numel(xs);
    plot(xs, cv, mk{d}, 'Color', cDet{d}, 'LineWidth', lwB, 'MarkerSize', msB, ...
        'MarkerIndices', 1:max(1,round(numel(xs)/12)):numel(xs));
end
xlim(cdfXL); ylim(cdfYL);
styleAxis(gca,'Per-user effective SE [bps/Hz]','CDF');

lgd = legend([hB(:); hS(:)], [nmDet(:); {'No IDD';'IDD iter 3'}], ...
    'NumColumns', 3, 'FontSize', 9, 'Box', 'on');
lgd.Layout.Tile = 'south';  lgd.ItemTokenSize = [tokW 18];

%% ====================================================================
%  FIGURE 2 (perfect CSI): left = fronthaul (3), right = complexity (4)
%% ====================================================================
figFronthaulComplexity('VTC-Fig2', xmax, xtSNR, tokW, msB);

%% ====================================================================
%  FIGURE 3: left = 8 BER (4 perfect DOTTED + 4 estimated SOLID), right = SE CDF (est, solid)
%% ====================================================================
lwB = 2.5;
figure('Name','VTC-Fig3','Position',[60 120 1250 560]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on; box on; grid on;
hB3 = gobjects(1,4);
for d = 1:4
    hB3(d) = semilogy(SNR_dB, max(aBERrt_pf(d,:),flrI), [':' mkpt{d}], 'Color', cDet{d}, ...
        'LineWidth', lwB, 'MarkerSize', msB);                    % perfect: dotted + marker
end
for d = 1:4
    semilogy(SNR_dB, max(aBERrt_est(d,:),flrI), mk{d}, 'Color', cDet{d}, ...
        'LineWidth', lwB, 'MarkerSize', msB);                    % estimated: solid + marker
end
hS3 = gobjects(1,2);
hS3(1) = semilogy(nan,nan,':', 'Color','k','LineWidth',lwB);
hS3(2) = semilogy(nan,nan,'-','Color','k','LineWidth',lwB);
set(gca,'YScale','log','XTick',xtSNR,'XLim',[0 xmax]); ylim([1e-5 1]);
styleAxis(gca,'SNR [dB]','BER');

nexttile; hold on; box on; grid on;
for d = 1:4
    xs = sort(seSamp_est{d});  cv = (1:numel(xs)).'/numel(xs);
    plot(xs, cv, mk{d}, 'Color', cDet{d}, 'LineWidth', lwB, 'MarkerSize', msB, ...
        'MarkerIndices', 1:max(1,round(numel(xs)/12)):numel(xs));   % estimated CSI, SOLID + marker
end
xlim(cdfXL); ylim(cdfYL);
styleAxis(gca,'Per-user effective SE [bps/Hz]','CDF');

lgd = legend([hB3(:); hS3(:)], [nmDet(:); {'Perfect CSI';'Estimated CSI'}], ...
    'NumColumns', 3, 'FontSize', 9, 'Box', 'on');
lgd.Layout.Tile = 'south';  lgd.ItemTokenSize = [tokW 18];

%% ====================================================================
%  FIGURE 4: same as Figure 2 (CSI-independent cost)
%% ====================================================================
figFronthaulComplexity('VTC-Fig4', xmax, xtSNR, tokW, msB);

%% ====================================================================
%  FIGURE 5 (single figure): 12 BER curves, all estimated CSI
%    NF estimated (dotted+marker) + FF estimated (dashed) + NF IDD-3 (dash-dot)
%    colour = detector.  Line style separates the three families.
%% ====================================================================
lwB = 2.5;
figure('Name','VTC-Fig5','Position',[80 120 760 620]);
hold on; box on; grid on;
hD = gobjects(1,4);
for d = 1:4
    hD(d) = semilogy(SNR_dB, max(aBERrt_est(d,:),flrI), mk{d}, 'Color', cDet{d}, ...
        'LineWidth', lwB, 'MarkerSize', msB);                    % NF estimated: solid + marker
end
for d = 1:4
    semilogy(SNR_dB, max(aBERrt_ff_est(d,:),flrI), [':' mkpt{d}], 'Color', cDet{d}, ...
        'LineWidth', lwB, 'MarkerSize', msB);                    % FF estimated: dotted + marker
end
for d = 1:4
    semilogy(SNRf, max(squeeze(BERidd_nf(d,3,:)).',flrI), ':', 'Color', cDet{d}, ...
        'LineWidth', lwB + 1.0);                                 % NF IDD-3: large dotted, no marker, detector colour
end
hK = gobjects(1,3);
hK(1) = semilogy(nan,nan,'-o','Color','k','LineWidth',lwB,'MarkerSize',msB);   % NF est key
hK(2) = semilogy(nan,nan,':o','Color','k','LineWidth',lwB,'MarkerSize',msB);   % FF est key
hK(3) = semilogy(nan,nan,':', 'Color','k','LineWidth',lwB + 1.0);              % IDD-3 key (large dotted)
set(gca,'YScale','log','XTick',xtSNR,'XLim',[0 xmax]); ylim([1e-5 1]);
styleAxis(gca,'SNR [dB]','BER');

lgd = legend([hD(:); hK(:)], ...
    [nmDet(:); {'Near-field, estimated';'Far-field, estimated';'Near-field, IDD iter 3'}], ...
    'NumColumns', 3, 'FontSize', 9, 'Box', 'on', 'Location', 'southoutside');
lgd.ItemTokenSize = [tokW 18];

% ==========================================================================
%  LOCAL FUNCTIONS
% ==========================================================================
function styleAxis(ax, xtxt, ytxt)
    set(ax, 'FontWeight', 'bold');                               % bold tick labels
    xlabel(ax, xtxt, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel(ax, ytxt, 'FontSize', 14, 'FontWeight', 'bold');
end

function BERidd = genIDD(aBER, SNR_dB, SNRf, cumGlog, penlog, smid, sw, smoothWin_dB, flrI)
    nF = numel(SNRf);
    baseS = zeros(4,nF);
    smW = max(3, round(smoothWin_dB/(SNRf(2)-SNRf(1))));
    dS  = SNRf(2) - SNRf(1);
    for d = 1:4
        yMon = cummin(log10(max(aBER(d,:), 1e-9)));
        yf   = interp1(SNR_dB, yMon, SNRf, 'pchip');
        slL  = (yf(2)-yf(1))/dS;   slR = (yf(end)-yf(end-1))/dS;
        yl   = yf(1)  + slL*dS*(-smW:-1);
        yr   = yf(end)+ slR*dS*(1:smW);
        ypad = movmean([yl, yf, yr], smW);
        baseS(d,:) = 10.^ypad(smW+1 : smW+nF);
    end
    Sc = 0.5*(1 + tanh((SNRf - smid)/sw));
    BERidd = zeros(4,3,nF);
    for d = 1:4
        for m = 1:3
            Delta = cumGlog(d,m).*Sc - penlog(d).*(1 - Sc);
            BERidd(d,m,:) = min(max(baseS(d,:).*10.^(-Delta), flrI), 1);
        end
    end
    BERidd = cummin(BERidd, 2);
end

function figFronthaulComplexity(figName, xmax, xtSNR, tokW, msB)
    cDet = {[0.85 0.10 0.10],[0.55 0.00 0.00],[0.20 0.20 0.75],[0.00 0.39 0.00]};
    SNR_dB = evalin('base','SNR_dB'); SNR_dB = SNR_dB(:).';
    L        = evalin('base','L');       N        = evalin('base','N');
    Ksweep   = evalin('base','Ksweep');
    cLin     = evalin('base','cLin');    cSIC     = evalin('base','cSIC');
    cList    = evalin('base','cList');   cXap     = evalin('base','cXap');
    Bbar_xap = evalin('base','Bbar_xap');
    Bbar_cluster = evalin('base','Bbar_cluster');
    lwB = 3.0;

    figure('Name',figName,'Position',[60 120 1250 560]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; box on; grid on;                          % LEFT: fronthaul (0..25 dB)
    hf1 = plot(SNR_dB, L*ones(size(SNR_dB)), ':k', 'LineWidth', 1.8);
    hf2 = plot(SNR_dB, Bbar_xap, '-d', 'Color', cDet{4}, 'LineWidth', lwB, 'MarkerSize', msB);
    hf3 = plot(SNR_dB, Bbar_cluster, '-s', 'Color', cDet{3}, 'LineWidth', lwB, 'MarkerSize', msB);
    set(gca,'XTick',xtSNR,'XLim',[0 xmax]); ylim([0 L+0.5]); set(gca,'FontWeight','bold');
    xlabel('SNR [dB]','FontSize',14,'FontWeight','bold');
    ylabel('Avg fronthaul scalars per user','FontSize',14,'FontWeight','bold');

    nexttile; hold on; box on; grid on;                          % RIGHT: complexity (unchanged)
    hc1 = plot(Ksweep, cXap,  '-d', 'Color', cDet{4}, 'LineWidth', lwB, 'MarkerSize', msB);
    hc2 = plot(Ksweep, cSIC,  '-s', 'Color', cDet{2}, 'LineWidth', lwB, 'MarkerSize', msB);
    hc3 = plot(Ksweep, cLin,  '-o', 'Color', cDet{1}, 'LineWidth', lwB, 'MarkerSize', msB);
    hc4 = plot(Ksweep, cList, '-^', 'Color', cDet{3}, 'LineWidth', lwB, 'MarkerSize', msB);
    set(gca,'XTick',Ksweep); set(gca,'FontWeight','bold');
    xlabel('Number of users K','FontSize',14,'FontWeight','bold');
    ylabel('Complex mult. per channel use','FontSize',14,'FontWeight','bold');

    lgd = legend([hf1 hf2 hf3 hc2 hc3 hc4], ...
        {'All-AP bound L','List+CrossAP (proposed)','Cluster-only', ...
         'SIC','Linear (L-MMSE)','List-SIC'}, ...
        'NumColumns', 3, 'FontSize', 9, 'Box', 'on');
    lgd.Layout.Tile = 'south';  lgd.ItemTokenSize = [tokW 18];
end

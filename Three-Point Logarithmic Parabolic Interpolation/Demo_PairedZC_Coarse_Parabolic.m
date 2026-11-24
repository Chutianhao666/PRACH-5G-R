clc;
clear;
close all;

%% ============================================================
% Paired ZC coarse peak + three-point logarithmic parabola demo
%
% Uses the user's actual project functions:
%   prachSigGen.m
%   prachFrontendProcessing.m
%   pairedZCFullLagCoarseEstimator.m
%   cal_bwp_info.m
%
% Controlled single-path experiment (no TDL-D / no MP-ACR):
%   - known fractional propagation delay
%   - known CFO
%   - optional AWGN
%
% Outputs:
%   Figure 1: full correlation curves for Branch 1 and Branch 2
%   Figure 2: local three-point logarithmic parabola fits
%   Figure 3: integer-peak versus interpolated delay/CFO errors
%   CSV/MAT result files
%% ============================================================

%% 0. Add this package folder to MATLAB path
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);

requiredFunctions = { ...
    'prachSigGen', ...
    'prachFrontendProcessing', ...
    'pairedZCFullLagCoarseEstimator', ...
    'cal_bwp_info'};

fprintf('\n===== Checking required MATLAB functions =====\n');
for k = 1:numel(requiredFunctions)
    fn = requiredFunctions{k};
    fp = which(fn);
    if isempty(fp)
        error('Demo:MissingFunction', ...
            'Missing function: %s. Check package files.',fn);
    end
    fprintf('%-38s %s\n',fn,fp);
end

if exist('dsp.FIRDecimator','class') ~= 8
    error('Demo:MissingDSPToolbox', ...
        ['The supplied PRACH front end uses dsp.FIRDecimator. ', ...
         'Please install/enable DSP System Toolbox.']);
end

%% 1. Reproducible simulation parameters
rng(20261124,'twister');

simConfig = struct( ...
    'nTx',1, ...
    'nRx',2, ...
    'puschBW',100, ...
    'puschSCS',30, ...
    'prachFormat','F0', ...
    'prachZeroConfigZone',13, ...
    'prachRootSeqIdx',22, ...
    'channelType','TDL', ...
    'snrdB',Inf, ...
    'simEpoch',1, ...
    'simAlgo',1, ...
    'CarrierFrequency',2.1e9, ...
    'MobileSpeed',500/3.6, ...
    'TrainPosition',-1000, ...
    'TrackDistance',100, ...
    'BaseHeight',30, ...
    'TrainAntennaHeight',5, ...
    'TransmissionDelay',0, ...
    'c',299792458);

% Ground-truth values. Change these to explore other cases.
trueDelay_effSamples = 3.29;
trueCFO_Hz = 360;
SNR_dB = 2;

%% 2. Generate paired-ZC PRACH waveform
[txSig,prachConfig] = prachSigGen(simConfig);

fsRaw = prachConfig.bwpInfo.SamplingRate;
N = prachConfig.bwpInfo.NFFT;
NFFTraw = prachConfig.nFFT;
NCPraw = round(prachConfig.NCP_RA);
Nrep = prachConfig.Nrep;
K = prachConfig.K;
fsEffective = fsRaw/K;

if Nrep ~= 1
    error('Demo:UnexpectedRepetitions', ...
        'This demonstration is configured for PRACH F0 with Nrep=1.');
end
if abs(NFFTraw/N-K) > 1e-9
    error('Demo:SamplingRatioMismatch', ...
        'Raw/effective FFT ratio does not match prachConfig.K.');
end

fprintf('\n===== Sampling information =====\n');
fprintf('Raw Fs       : %.6f MHz\n',fsRaw/1e6);
fprintf('Effective Fs : %.6f MHz\n',fsEffective/1e6);
fprintf('Effective N  : %d\n',N);
fprintf('K            : %g\n',K);

%% 3. Construct a known single-path delayed signal
rawBurstLength = NCPraw+NFFTraw*Nrep;
if size(txSig,1) < rawBurstLength
    error('Demo:ShortWaveform','Generated TX waveform is shorter than PRACH burst.');
end

txBurst = txSig(1:rawBurstLength,1);
trueDelay_rawSamples = trueDelay_effSamples*K;
rxDelayed = fractionalDelayRawLocal( ...
    txBurst,trueDelay_rawSamples,256);

% Two identical receive branches; no channel fading or gain normalization.
rxClean = repmat(rxDelayed,1,simConfig.nRx);

%% 4. Inject the known CFO
nRaw = (0:size(rxClean,1)-1).';
cfoRotation = exp(1j*2*pi*trueCFO_Hz*nRaw/fsRaw);
rxClean = bsxfun(@times,rxClean,cfoRotation);

%% 5. Add AWGN (SNR is defined using active clean-signal support)
if isinf(SNR_dB)
    rxNoisy = rxClean;
else
    activeEnvelope = max(abs(rxClean),[],2);
    activePeak = max(activeEnvelope);
    activeIdx = find(activeEnvelope > activePeak*1e-12);
    if isempty(activeIdx)
        error('Demo:EmptySignalSupport','No active TX samples were found.');
    end
    signalPower = mean(abs(rxClean(activeIdx,:)).^2,'all');
    noiseVariance = signalPower/10^(SNR_dB/10);
    noiseSigma = sqrt(noiseVariance/2);
    noise = noiseSigma*(randn(size(rxClean))+1j*randn(size(rxClean)));
    rxNoisy = rxClean+noise;
end

%% 6. Apply the actual PRACH front end
% The same front end is used for the clean signal to derive fixed gains.
timeDataClean = prachFrontendProcessing(rxClean,prachConfig,simConfig);
timeData = prachFrontendProcessing(rxNoisy,prachConfig,simConfig);

timeFlatClean = reshape(timeDataClean,simConfig.nRx,[]);
rxPower = mean(abs(timeFlatClean).^2,2);
if any(~isfinite(rxPower)) || any(rxPower<=0)
    error('Demo:InvalidRxPower','Cannot derive positive finite RX gains.');
end

ref1 = ifft(prachConfig.fdRootSeq,N)*sqrt(N);
refPower = mean(abs(ref1).^2);
if ~isfinite(refPower) || refPower<=0
    error('Demo:InvalidReferencePower','Cannot derive positive finite reference gain.');
end
rxGain = 1./sqrt(rxPower);
refGain = 1/sqrt(refPower);

%% 7. Run the actual full-lag paired-ZC coarse estimator
opts = struct();
opts.MaxAbsCFO_Hz = 1500;
opts.MinDelaySamples = 0;
opts.MaxDelaySamples = 20;
opts.SearchGuardBins = 2;
opts.ReturnCorrelations = true;
opts.RxGain = rxGain;
opts.RefGain = refGain;

[delayParabolic,cfoParabolic,diagInfo] = ...
    pairedZCFullLagCoarseEstimator(timeData,prachConfig,opts);

%% 8. Retrieve full correlation curves and independently reproduce fit
lags = diagInfo.lags(:);
corr1 = diagInfo.branch1Correlation(:);
corr2 = diagInfo.branch2Correlation(:);

fit1 = getPeakFit(corr1,lags,diagInfo.branch1PeakIndex);
fit2 = getPeakFit(corr2,lags,diagInfo.branch2PeakIndex);

%% 9. Theoretical branch peak coordinates
OS = N/prachConfig.L_RA;
cfoBinSpacing_Hz = prachConfig.scs_RA*1e3/OS;
trueCFOBins = trueCFO_Hz/cfoBinSpacing_Hz;
trueBranch1 = trueDelay_effSamples+trueCFOBins;
trueBranch2 = trueDelay_effSamples-trueCFOBins;

%% 10. Integer-peak combined estimates
% The estimator's integer peak coordinates are combined using the same
% paired-branch delay/CFO equations as the sub-bin estimates.
delayInteger = (fit1.lagInteger+fit2.lagInteger)/2;
cfoInteger = (fit1.lagInteger-fit2.lagInteger)/2*cfoBinSpacing_Hz;

%% 11. Results and error table
fprintf('\n====================================================\n');
fprintf('TRUE VALUE VS ESTIMATED VALUE\n');
fprintf('====================================================\n');
fprintf('True delay          : %.6f effective samples\n',trueDelay_effSamples);
fprintf('Integer-peak delay  : %.6f effective samples\n',delayInteger);
fprintf('Parabolic delay     : %.6f effective samples\n',delayParabolic);
fprintf('Integer delay error : %+.6f samples\n',delayInteger-trueDelay_effSamples);
fprintf('Parabolic delay err : %+.6f samples\n',delayParabolic-trueDelay_effSamples);
fprintf('\nTrue CFO            : %+.6f Hz\n',trueCFO_Hz);
fprintf('Integer-peak CFO    : %+.6f Hz\n',cfoInteger);
fprintf('Parabolic CFO       : %+.6f Hz\n',cfoParabolic);
fprintf('Integer CFO error   : %+.6f Hz\n',cfoInteger-trueCFO_Hz);
fprintf('Parabolic CFO error : %+.6f Hz\n',cfoParabolic-trueCFO_Hz);

Branch = ["Branch 1";"Branch 2"];
TrueCoordinate = [trueBranch1;trueBranch2];
IntegerPeak = [fit1.lagInteger;fit2.lagInteger];
ParabolicPeak = [fit1.lagSub;fit2.lagSub];
IntegerError = IntegerPeak-TrueCoordinate;
ParabolicError = ParabolicPeak-TrueCoordinate;
branchTable = table(Branch,TrueCoordinate,IntegerPeak, ...
    ParabolicPeak,IntegerError,ParabolicError);

fprintf('\n===== Individual branch peak comparison =====\n');
disp(branchTable);

summaryTable = table( ...
    trueDelay_effSamples,delayInteger,delayParabolic, ...
    delayInteger-trueDelay_effSamples, ...
    delayParabolic-trueDelay_effSamples, ...
    trueCFO_Hz,cfoInteger,cfoParabolic, ...
    cfoInteger-trueCFO_Hz,cfoParabolic-trueCFO_Hz,SNR_dB);

%% 12. Figure 1: full correlation curves
fig1 = figure('Color','w','Position',[80 80 1250 760]);
tl1 = tiledlayout(fig1,2,1,'TileSpacing','compact','Padding','compact');

nexttile(tl1);
plot(lags,abs(corr1),'LineWidth',1.25); hold on;
xline(trueBranch1,'k--','LineWidth',1.5);
xline(fit1.lagInteger,'r:','LineWidth',1.5);
xline(fit1.lagSub,'g-.','LineWidth',1.5);
scatter(fit1.lagInteger,fit1.peakMagnitude,55,'r','filled');
grid on;
xlim([diagInfo.searchLagMin diagInfo.searchLagMax]);
xlabel('Lag coordinate (effective samples)');
ylabel('|Correlation|');
title('Branch 1: full correlation peak');
legend('Correlation magnitude','Theoretical branch coordinate', ...
    'Integer peak','Log-parabolic peak','Integer peak sample', ...
    'Location','best');
hold off;

nexttile(tl1);
plot(lags,abs(corr2),'LineWidth',1.25); hold on;
xline(trueBranch2,'k--','LineWidth',1.5);
xline(fit2.lagInteger,'r:','LineWidth',1.5);
xline(fit2.lagSub,'g-.','LineWidth',1.5);
scatter(fit2.lagInteger,fit2.peakMagnitude,55,'r','filled');
grid on;
xlim([diagInfo.searchLagMin diagInfo.searchLagMax]);
xlabel('Lag coordinate (effective samples)');
ylabel('|Correlation|');
title('Branch 2: full correlation peak');
legend('Correlation magnitude','Theoretical branch coordinate', ...
    'Integer peak','Log-parabolic peak','Integer peak sample', ...
    'Location','best');
hold off;
title(tl1,'Paired-ZC full-lag coarse peak search');

%% 13. Figure 2: three-point logarithmic parabola fitting
fig2 = figure('Color','w','Position',[100 100 1250 600]);
tl2 = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

nexttile(tl2);
plot(fit1.localLags,fit1.localMagnitude,'ko-', ...
    'LineWidth',1.3,'MarkerSize',7,'DisplayName','Three actual samples');
hold on;
plot(fit1.fitLags,fit1.fitMagnitude,'--','LineWidth',2, ...
    'DisplayName','Exponentiated log-parabola');
scatter(fit1.lagInteger,fit1.peakMagnitude,75,'r','filled', ...
    'DisplayName','Integer peak');
scatter(fit1.lagSub,fit1.fittedPeak,85,'g','filled', ...
    'DisplayName','Interpolated vertex');
xline(trueBranch1,'k:','LineWidth',1.5,'DisplayName','Theoretical coordinate');
grid on;
xlabel('Lag coordinate (effective samples)'); ylabel('|Correlation|');
title(sprintf('Branch 1: integer %.4f, interpolated %.4f', ...
    fit1.lagInteger,fit1.lagSub));
legend('Location','best'); hold off;

nexttile(tl2);
plot(fit2.localLags,fit2.localMagnitude,'ko-', ...
    'LineWidth',1.3,'MarkerSize',7,'DisplayName','Three actual samples');
hold on;
plot(fit2.fitLags,fit2.fitMagnitude,'--','LineWidth',2, ...
    'DisplayName','Exponentiated log-parabola');
scatter(fit2.lagInteger,fit2.peakMagnitude,75,'r','filled', ...
    'DisplayName','Integer peak');
scatter(fit2.lagSub,fit2.fittedPeak,85,'g','filled', ...
    'DisplayName','Interpolated vertex');
xline(trueBranch2,'k:','LineWidth',1.5,'DisplayName','Theoretical coordinate');
grid on;
xlabel('Lag coordinate (effective samples)'); ylabel('|Correlation|');
title(sprintf('Branch 2: integer %.4f, interpolated %.4f', ...
    fit2.lagInteger,fit2.lagSub));
legend('Location','best'); hold off;
title(tl2,'Three-point logarithmic parabolic interpolation');

%% 14. Figure 3: final delay and CFO errors
fig3 = figure('Color','w','Position',[120 120 1100 460]);
tl3 = tiledlayout(fig3,1,2,'TileSpacing','compact','Padding','compact');

nexttile(tl3);
bar([delayInteger-trueDelay_effSamples, ...
     delayParabolic-trueDelay_effSamples]);
hold on; yline(0,'k--','LineWidth',1.2); grid on;
set(gca,'XTickLabel',{'Integer peak','Log-parabolic'});
ylabel('Delay error (effective samples)');
title('Delay estimation error'); hold off;

nexttile(tl3);
bar([cfoInteger-trueCFO_Hz,cfoParabolic-trueCFO_Hz]);
hold on; yline(0,'k--','LineWidth',1.2); grid on;
set(gca,'XTickLabel',{'Integer peak','Log-parabolic'});
ylabel('CFO error (Hz)');
title('CFO estimation error'); hold off;
title(tl3,'Integer peak versus sub-bin estimation');

%% 15. Save reproducible results
writetable(branchTable,'Demo_PairedZC_BranchPeakComparison.csv');
writetable(summaryTable,'Demo_PairedZC_EstimationSummary.csv');
save('Demo_PairedZC_Coarse_Parabolic.mat', ...
    'simConfig','trueDelay_effSamples','trueCFO_Hz','SNR_dB', ...
    'delayInteger','delayParabolic','cfoInteger','cfoParabolic', ...
    'branchTable','summaryTable','diagInfo','fit1','fit2', ...
    'fsRaw','fsEffective','N','K','cfoBinSpacing_Hz');

savefig(fig1,'Demo_Fig1_FullCorrelation.fig');
exportgraphics(fig1,'Demo_Fig1_FullCorrelation.png','Resolution',220);
savefig(fig2,'Demo_Fig2_ParabolicFit.fig');
exportgraphics(fig2,'Demo_Fig2_ParabolicFit.png','Resolution',220);
savefig(fig3,'Demo_Fig3_EstimationError.fig');
exportgraphics(fig3,'Demo_Fig3_EstimationError.png','Resolution',220);

fprintf('\nDemo completed. Figures, CSV and MAT files were saved here:\n%s\n',projectDir);

%% ============================================================
% Local function: extract/reproduce the estimator's three-point fit
%% ============================================================
function out = getPeakFit(c,lags,peakIndex)
mag = abs(c(:));
lags = lags(:);
if peakIndex<=1 || peakIndex>=numel(c)
    error('Demo:PeakAtCorrelationBoundary', ...
        'Peak lies at full-correlation boundary; three-point fit is unavailable.');
end

pMinus = mag(peakIndex-1);
pZero  = mag(peakIndex);
pPlus  = mag(peakIndex+1);

yMinus = log(max(pMinus,realmin));
yZero  = log(max(pZero,realmin));
yPlus  = log(max(pPlus,realmin));

den = yMinus-2*yZero+yPlus;
if isfinite(den) && abs(den)>1e-12
    delta = 0.5*(yMinus-yPlus)/den;
    delta = max(-0.5,min(0.5,delta));
else
    delta = 0;
end

lagInteger = lags(peakIndex);
lagSub = lagInteger+delta;

% Quadratic in local coordinate x = -1,0,+1:
% log|C| = a*x^2+b*x+c
% Its vertex is x = -b/(2a).
a = 0.5*den;
b = 0.5*(yPlus-yMinus);
c = yZero;

xFit = linspace(-1,1,301);
fitLog = a*xFit.^2+b*xFit+c;
fitMagnitude = exp(fitLog);
fittedPeak = exp(a*delta^2+b*delta+c);

out.lagInteger = lagInteger;
out.lagSub = lagSub;
out.delta = delta;
out.peakMagnitude = pZero;
out.fittedPeak = fittedPeak;
out.localLags = lags(peakIndex-1:peakIndex+1);
out.localMagnitude = mag(peakIndex-1:peakIndex+1);
out.fitLags = lagInteger+xFit;
out.fitMagnitude = fitMagnitude;
out.logPoints = [yMinus,yZero,yPlus];
end

%% ============================================================
% Local function: apply a fractional raw-sample delay
%% ============================================================
function y = fractionalDelayRawLocal(x,delaySamples,guardSamples)
x = x(:);
if ~isscalar(delaySamples) || ~isfinite(delaySamples) || delaySamples<0
    error('Demo:InvalidFractionalDelay','Delay must be finite and nonnegative.');
end

integerDelay = floor(delaySamples);
fractionalPart = delaySamples-integerDelay;
xInteger = [zeros(integerDelay,1,'like',x);x];
if fractionalPart<1e-12
    y = xInteger;
    return;
end

guardSamples = max(guardSamples,32);
xPad = [zeros(guardSamples,1,'like',xInteger); ...
        xInteger;zeros(guardSamples,1,'like',xInteger)];
nFFT = 2^nextpow2(2*numel(xPad));
k = [0:nFFT/2,-nFFT/2+1:-1].';
X = fft(xPad,nFFT);
Y = ifft(X.*exp(-1j*2*pi*k*fractionalPart/nFFT));
first = guardSamples+1;
last = first+numel(x)+integerDelay;
y = Y(first:last);
end

clc;
clear;
close all;
clear functions;
rehash;

%% ============================================================
% V13.6 - TDL-D MULTIPATH VALIDATION OF THE V13.5 FULL-LAG COARSE
%         ESTIMATOR
%
% Purpose:
%   Validate the unchanged V13.5 full-lag/sub-bin coarse estimator under
%   a TDL-D multipath channel before connecting the MP-ACR fine estimator.
%
% Important modeling separation:
%   1) nrTDLChannel uses TDL-D with MaximumDopplerShift = 0 in this stage.
%      This deliberately freezes the multipath taps during one PRACH burst,
%      so that multipath peak-selection behavior can be isolated.
%   2) The railway LOS Doppler and residual oscillator offset are applied
%      as one known common CFO after the TDL-D channel:
%           CFO = geometry Doppler + residual CFO.
%   3) This is a multipath/coarse-estimator validation stage, NOT yet the
%      final time-varying high-speed fading validation. A later stage can
%      introduce time-varying path Doppler without double-counting CFO.
%
% Channel:
%   DelayProfile           = TDL-D
%   DelaySpread            = 100 ns (RMS delay spread)
%   KFactorScaling         = false (retain standard TDL-D first-tap K)
%   MaximumDopplerShift    = 0 Hz (static taps for this isolation stage)
%   NumTransmitAntennas    = 1
%   NumReceiveAntennas     = 2
%
% Signal:
%   - Same paired-ZC waveform generation as V13.5
%   - Same prachFrontendProcessing function
%   - Same pairedZCFullLagCoarseEstimator function, unchanged
%   - No time_freq_recovery / V6 / prachPowerDelayProfileGen
%   - No transmit waveform energy normalization
%   - No per-algorithm received-power normalization
%
% TDL-D channel-filter delay:
%   The integer implementation delay reported by info(nrTDLChannel) is
%   removed before applying the known railway propagation delay. The
%   physical TDL path delays remain in the received waveform.
%
% Delay calibration:
%   An independent set of noise-free TDL-D zero-propagation-delay trials
%   estimates one fixed delay-coordinate bias. Both raw and calibrated
%   delay errors are exported. Calibration data are NOT included in the
%   reported Monte-Carlo performance trials.
%
% 80-point one-dimensional fit experiment:
%   CFO mode:
%       80 distinct CFO values, fixed propagation delay = 4 effective samples
%       3 independent TDL-D channel/noise realizations per CFO point
%       240 evaluation trials total
%   Delay mode:
%       300 distinct propagation-delay values, fixed total CFO = railway Doppler
%       (zero additional CFO), 3 independent TDL-D channel/noise realizations
%       per delay point, 240 evaluation trials total
%
% Each scatter figure displays exactly 80 markers. Each marker is the mean
% estimate from 3 Monte-Carlo repetitions at that true parameter value.
% CFO and delay are swept independently; no 80 x 300 Cartesian grid is used.
%
% Required:
%   5G Toolbox (nrTDLChannel)
%   Communications/DSP functionality used by the supplied front end
%
% Outputs:
%   V13_6_TDLD_CFO80_Scatter_Raw.csv
%   V13_6_TDLD_CFO80_Scatter_Metrics.csv
%   V13_6_TDLD_CFO80_Scatter.mat
%   V13_6_TDLD_CFO80_Scatter_Performance.png/.fig
%   V13_6_TDLD_CFO80_Scatter_Accuracy.png/.fig
%% ============================================================

fprintf('\n============================================================\n');
fprintf('V13.6 TDL-D 80-POINT CFO ESTIMATION SCATTER\n');
fprintf('============================================================\n');

%% 0. Configuration

cfg.CarrierFrequency_Hz = 2.1e9;
cfg.MobileSpeed_kmh = 500;
cfg.nTx = 1;
cfg.nRx = 2;

cfg.puschBW = 100;
cfg.puschSCS = 30;
cfg.prachFormat = 'F0';
cfg.prachZeroConfigZone = 13;
cfg.prachRootSeqIdx = 22;

% Railway LOS Doppler at 500 km/h, 2.1 GHz, initial geometry.
cfg.TrainPosition_m = -1000;
cfg.TrackDistance_m = 100;
cfg.BaseHeight_m = 30;
cfg.TrainAntennaHeight_m = 5;
cfg.c_mps = 299792458;

[fd,rangeAtT0_m,propagationDelayAtT0_s] = railwayDopplerModel( ...
    cfg.TrainPosition_m,cfg.MobileSpeed_kmh/3.6, ...
    cfg.TrackDistance_m, ...
    cfg.BaseHeight_m-cfg.TrainAntennaHeight_m, ...
    cfg.CarrierFrequency_Hz,0);

cfg.geometryDoppler_Hz = fd;
cfg.rangeAtT0_m = rangeAtT0_m;
cfg.propagationDelayAtT0_s = propagationDelayAtT0_s;

% CFO-only sweep: 300 values across railway Doppler +/- 300 Hz.
% Fix propagation delay at 4 effective samples for this independent sweep.
cfg.CFOList_Hz = linspace(-fd-300, fd+300, 80);
cfg.delayList_effSamples = 4.0;
cfg.fitMode = 'CFO';

% Representative SNR retained from the previous V13.6 accuracy plot.
cfg.snrList_dB = 0;

% TDL-D configuration. DelaySpread is RMS delay spread, not maximum delay.
cfg.tdlDelayProfile = 'TDL-D';
cfg.tdlDelaySpread_s = 100e-9;

% Stage A isolates multipath with a static channel realization per burst.
% Common railway Doppler is explicitly included in CFOList_Hz above.
cfg.tdlMaximumDopplerShift_Hz = 0;

% Independent channel realizations.
cfg.nChannelRealizations = 3;
cfg.nCalibrationChannelRealizations = 3;

cfg.randomSeed = 20261121;
cfg.fractionalDelayGuardSamples = 256;

% The search guard is set after reading the actual scaled TDL-D path delays.
cfg.searchGuardBinsExtra = 2;

% Output controls: retain only the two principal result figures by default.
cfg.makeFigures = true;
cfg.makeDiagnosticFigures = false;

% Success thresholds are diagnostic metrics, not detection decisions.
cfg.CFOAcceptanceThreshold_Hz = 50;
cfg.DelayAcceptanceThreshold_effSamples = 0.1;

rng(cfg.randomSeed,'twister');

%% 1. Verify active functions and required toolbox

requiredFunctions = { ...
    'prachSigGen', ...
    'prachFrontendProcessing', ...
    'pairedZCFullLagCoarseEstimator', ...
    'railwayDopplerModel'};

fprintf('\n===== Active MATLAB function paths =====\n');
for k = 1:numel(requiredFunctions)
    p = which(requiredFunctions{k});
    if isempty(p)
        error('V13_6:MissingFunction', ...
            'Required function not found: %s',requiredFunctions{k});
    end
    fprintf('%-40s %s\n',requiredFunctions{k},p);
end

if exist('nrTDLChannel','file') ~= 2 && exist('nrTDLChannel','class') ~= 8
    error('V13_6:Missing5GToolbox', ...
        ['nrTDLChannel was not found. This V13.6 script requires ', ...
         '5G Toolbox.']);
end

fprintf('\nEstimators intentionally NOT called:\n');
fprintf('  time_freq_recovery\n');
fprintf('  time_freq_recovery_V6\n');
fprintf('  prachPowerDelayProfileGen\n');

%% 2. Generate the same paired-ZC waveform as V13.5

simConfig = struct( ...
    'nTx',cfg.nTx, ...
    'nRx',cfg.nRx, ...
    'puschBW',cfg.puschBW, ...
    'puschSCS',cfg.puschSCS, ...
    'prachFormat',cfg.prachFormat, ...
    'prachZeroConfigZone',cfg.prachZeroConfigZone, ...
    'prachRootSeqIdx',cfg.prachRootSeqIdx, ...
    'channelType','TDL', ...
    'snrdB',Inf, ...
    'simEpoch',1, ...
    'simAlgo',1, ...
    'CarrierFrequency',cfg.CarrierFrequency_Hz, ...
    'MobileSpeed',cfg.MobileSpeed_kmh/3.6, ...
    'TrainPosition',cfg.TrainPosition_m, ...
    'TrackDistance',cfg.TrackDistance_m, ...
    'BaseHeight',cfg.BaseHeight_m, ...
    'TrainAntennaHeight',cfg.TrainAntennaHeight_m, ...
    'TransmissionDelay',0, ...
    'c',cfg.c_mps);

[txSig,prachConfig] = prachSigGen(simConfig);

fsRaw = prachConfig.bwpInfo.SamplingRate;
NIFFT = prachConfig.bwpInfo.NFFT;
NFFTraw = prachConfig.nFFT;
NCPraw = round(prachConfig.NCP_RA);
Nrep = prachConfig.Nrep;
K = prachConfig.K;
fsEffective = fsRaw/K;

if Nrep ~= 1
    error('V13_6:UnexpectedNrep', ...
        'This validation expects PRACH F0 with Nrep=1.');
end
if abs(NFFTraw/NIFFT-K)>1e-9
    error('V13_6:SamplingRatioMismatch', ...
        'Raw/effective FFT ratio does not match prachConfig.K.');
end

rawBurstLength = NCPraw+NFFTraw*Nrep;
if size(txSig,1)<rawBurstLength
    error('V13_6:ShortTX', ...
        'Generated TX waveform is shorter than the active PRACH burst.');
end

txBurst = txSig(1:rawBurstLength,1);
txBurst = txBurst(:);

OS = NIFFT/prachConfig.L_RA;
cfoBinSpacing_Hz = prachConfig.scs_RA*1e3/OS;

maxAbsCFO_Hz = max(abs(cfg.CFOList_Hz));
maxCFOBins = maxAbsCFO_Hz/cfoBinSpacing_Hz;

%% 3. Configure TDL-D and inspect its actual scaled delay profile

tdlChan = makeTDLChannelLocal(cfg,fsRaw,cfg.randomSeed);
tdlInfo = info(tdlChan);

if ~isfield(tdlInfo,'ChannelFilterDelay') || ...
        ~isfield(tdlInfo,'MaximumChannelDelay') || ...
        ~isfield(tdlInfo,'PathDelays')
    error('V13_6:UnexpectedTDLInfo', ...
        'nrTDLChannel info structure lacks required delay fields.');
end

channelFilterDelay_raw = round(tdlInfo.ChannelFilterDelay);
tdlMaxExcessDelay_raw = max(tdlInfo.PathDelays)*fsRaw;
tdlMaxExcessDelay_eff = max(tdlInfo.PathDelays)*fsEffective;

% Extend the signed-lag search window to include the TDL-D excess-delay
% support in addition to CFO and absolute propagation delay.
cfg.searchGuardBins = ceil(tdlMaxExcessDelay_eff) + ...
    cfg.searchGuardBinsExtra;

lagSearchMin = floor(min(cfg.delayList_effSamples)- ...
    maxCFOBins-cfg.searchGuardBins);
lagSearchMax = ceil(max(cfg.delayList_effSamples)+ ...
    maxCFOBins+cfg.searchGuardBins);

tailLength_raw = max(64,ceil(tdlInfo.MaximumChannelDelay)+32);

fprintf('\n===== Waveform configuration =====\n');
fprintf('Format / Root / ZCZ        : %s / %d / %d\n', ...
    cfg.prachFormat,cfg.prachRootSeqIdx,prachConfig.N_CS);
fprintf('NFFT / Nrep / K            : %d / %d / %d\n',NIFFT,Nrep,K);
fprintf('Raw / effective Fs         : %.6f / %.6f MHz\n', ...
    fsRaw/1e6,fsEffective/1e6);
fprintf('CFO bin spacing            : %.9f Hz/bin\n',cfoBinSpacing_Hz);
fprintf('CFO range                  : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));
fprintf('Absolute delay range       : [%.4f, %.4f] effective samples\n', ...
    min(cfg.delayList_effSamples),max(cfg.delayList_effSamples));
fprintf('TX burst energy            : %.12e (not normalized)\n', ...
    sum(abs(txBurst).^2));

fprintf('\n===== TDL-D channel configuration =====\n');
fprintf('Profile                    : %s\n',cfg.tdlDelayProfile);
fprintf('RMS delay spread           : %.3f ns\n',cfg.tdlDelaySpread_s*1e9);
fprintf('Maximum Doppler in TDL     : %.3f Hz\n', ...
    cfg.tdlMaximumDopplerShift_Hz);
fprintf('Railway common CFO         : %+.6f Hz\n',fd);
fprintf('Train-BS range at t=0      : %.3f m\n',rangeAtT0_m);
fprintf('Geometric propagation delay: %.6f us\n', ...
    propagationDelayAtT0_s*1e6);
fprintf('Channel filter delay       : %d raw samples\n', ...
    channelFilterDelay_raw);
fprintf('Maximum TDL excess delay   : %.3f raw / %.4f effective samples\n', ...
    tdlMaxExcessDelay_raw,tdlMaxExcessDelay_eff);
fprintf('First-tap K factor         : %.3f dB\n',tdlInfo.KFactorFirstTap);
fprintf('Search guard               : %d effective samples\n', ...
    cfg.searchGuardBins);
fprintf('Signed-lag search          : [%d, %d]\n', ...
    lagSearchMin,lagSearchMax);
fprintf('Tail padding               : %d raw samples\n',tailLength_raw);

maxTotalDelayRaw = ...
    max(cfg.delayList_effSamples)*K + tdlMaxExcessDelay_raw;
if maxTotalDelayRaw >= NCPraw
    warning('V13_6:DelayExceedsCP', ...
        ['Absolute plus TDL excess delay (%.2f raw samples) is not ', ...
         'smaller than the PRACH CP (%d raw samples).'], ...
         maxTotalDelayRaw,NCPraw);
end

%% 4. Estimator options (same estimator and same coordinate model as V13.5)

estimatorOpts = struct();
estimatorOpts.MaxAbsCFO_Hz = maxAbsCFO_Hz;
estimatorOpts.MinDelaySamples = min(cfg.delayList_effSamples);
estimatorOpts.MaxDelaySamples = max(cfg.delayList_effSamples);
estimatorOpts.SearchGuardBins = cfg.searchGuardBins;
estimatorOpts.ReturnCorrelations = false;

%% 5. Independent TDL-D zero-delay calibration

fprintf('\n===== Independent TDL-D delay calibration =====\n');

calibrationDelayEstimates = nan( ...
    cfg.nCalibrationChannelRealizations*numel(cfg.CFOList_Hz),1);

calIdx = 0;

for iCal = 1:cfg.nCalibrationChannelRealizations

    calChan = makeTDLChannelLocal( ...
        cfg,fsRaw,cfg.randomSeed+1000+iCal);
    rxFaded = generateTDLRealizationLocal( ...
        calChan,txBurst,tailLength_raw,channelFilterDelay_raw);

    for iCFO = 1:numel(cfg.CFOList_Hz)

        currentCFO = cfg.CFOList_Hz(iCFO);

        rxDelayed = fractionalDelayRawLocal( ...
            rxFaded,0,cfg.fractionalDelayGuardSamples);

        nRaw = (0:size(rxDelayed,1)-1).';
        cfoRot = exp(1j*2*pi*currentCFO*nRaw/fsRaw);
        rxClean = bsxfun(@times,rxDelayed,cfoRot);

        tdClean = prachFrontendProcessing( ...
            rxClean,prachConfig,simConfig);

        [rxGain,refGain] = deriveSharedGainsLocal( ...
            tdClean,prachConfig,cfg.nRx);

        optsNow = estimatorOpts;
        optsNow.RxGain = rxGain;
        optsNow.RefGain = refGain;

        [delayCal,~,~] = pairedZCFullLagCoarseEstimator( ...
            tdClean,prachConfig,optsNow);

        calIdx = calIdx+1;
        calibrationDelayEstimates(calIdx) = delayCal;
    end

    fprintf('Calibration channel realization %d/%d complete.\n', ...
        iCal,cfg.nCalibrationChannelRealizations);
end

delayBiasCalibrationSamples = median(calibrationDelayEstimates,'omitnan');
delayBiasMAD_samples = median( ...
    abs(calibrationDelayEstimates-delayBiasCalibrationSamples),'omitnan');

if ~isfinite(delayBiasCalibrationSamples)
    error('V13_6:InvalidDelayCalibration', ...
        'TDL-D zero-delay calibration produced an invalid bias.');
end

fprintf('TDL-D calibration bias       : %+.9f effective samples\n', ...
    delayBiasCalibrationSamples);
fprintf('TDL-D calibration MAD        : %.9f effective samples\n', ...
    delayBiasMAD_samples);

if delayBiasMAD_samples > 0.1
    warning('V13_6:UnstableDelayCalibration', ...
        ['TDL-D zero-delay calibration MAD is %.4f samples. ', ...
         'A single fixed delay correction may not be reliable. ', ...
         'Use raw delay metrics as the primary result.'], ...
         delayBiasMAD_samples);
end

%% 6. Preallocate Monte-Carlo results

nCFO = numel(cfg.CFOList_Hz);
nDelay = numel(cfg.delayList_effSamples);
nSNR = numel(cfg.snrList_dB);
nRep = cfg.nChannelRealizations;

nTrials = nRep*nCFO*nDelay*nSNR;

TrialID = (1:nTrials).';
ChannelRealization = nan(nTrials,1);
TrueCFO_Hz = nan(nTrials,1);
TrueDelay_effSamples = nan(nTrials,1);
TrueDelay_rawSamples = nan(nTrials,1);
SNR_dB = nan(nTrials,1);

Branch1LagInteger = nan(nTrials,1);
Branch1LagSubBin = nan(nTrials,1);
Branch2LagInteger = nan(nTrials,1);
Branch2LagSubBin = nan(nTrials,1);
Branch1PeakValue = nan(nTrials,1);
Branch2PeakValue = nan(nTrials,1);
Branch2ToBranch1Peak_dB = nan(nTrials,1);
Branch1PeakToMedian_dB = nan(nTrials,1);
Branch2PeakToMedian_dB = nan(nTrials,1);

CFOEstimate_Hz = nan(nTrials,1);
DelayEstimateRaw_effSamples = nan(nTrials,1);
DelayEstimateCal_effSamples = nan(nTrials,1);
CFOError_Hz = nan(nTrials,1);
DelayErrorRaw_effSamples = nan(nTrials,1);
DelayErrorCal_effSamples = nan(nTrials,1);

Branch1LagError_effSamples = nan(nTrials,1);
Branch2LagError_effSamples = nan(nTrials,1);
Peak1AtSearchBoundary = false(nTrials,1);
Peak2AtSearchBoundary = false(nTrials,1);
NoiseVariance = nan(nTrials,1);
RxSignalPower = nan(nTrials,1);

%% 7. Main TDL-D Monte-Carlo loop

fprintf('\n===== Main TDL-D Monte-Carlo run =====\n');

trial = 0;

for iRep = 1:nRep

    % One independent static TDL-D realization per Monte-Carlo epoch.
    % All CFO/delay/SNR points within this epoch reuse the same channel.
    evalChan = makeTDLChannelLocal( ...
        cfg,fsRaw,cfg.randomSeed+10000+iRep);
    rxFaded = generateTDLRealizationLocal( ...
        evalChan,txBurst,tailLength_raw,channelFilterDelay_raw);

    for iCFO = 1:nCFO

        currentCFO = cfg.CFOList_Hz(iCFO);

        for iDelay = 1:nDelay

            currentDelayEff = cfg.delayList_effSamples(iDelay);
            currentDelayRaw = currentDelayEff*K;

            rxDelayed = fractionalDelayRawLocal( ...
                rxFaded,currentDelayRaw,cfg.fractionalDelayGuardSamples);

            nRaw = (0:size(rxDelayed,1)-1).';
            cfoRot = exp(1j*2*pi*currentCFO*nRaw/fsRaw);
            rxClean = bsxfun(@times,rxDelayed,cfoRot);

            % Fixed gains from the clean waveform. No waveform-energy
            % normalization is performed; the two receive antennas use
            % one gain each, shared by both paired-ZC branches.
            tdClean = prachFrontendProcessing( ...
                rxClean,prachConfig,simConfig);

            [rxGain,refGain] = deriveSharedGainsLocal( ...
                tdClean,prachConfig,cfg.nRx);

            optsNow = estimatorOpts;
            optsNow.RxGain = rxGain;
            optsNow.RefGain = refGain;

            activeEnvelope = max(abs(rxClean),[],2);
            activePeak = max(activeEnvelope);
            activeIdx = find(activeEnvelope>activePeak*1e-12);

            if isempty(activeIdx)
                error('V13_6:EmptyActiveSupport', ...
                    'The TDL-D received waveform has no active samples.');
            end

            pRxSignal = mean(abs(rxClean(activeIdx,:)).^2,'all');

            if ~isfinite(pRxSignal) || pRxSignal<=0
                error('V13_6:InvalidRXPower', ...
                    'Invalid clean TDL-D receive power.');
            end

            for iSNR = 1:nSNR

                currentSNR = cfg.snrList_dB(iSNR);

                if isinf(currentSNR)
                    noiseVariance = 0;
                else
                    noiseVariance = pRxSignal/10^(currentSNR/10);
                end

                noiseSigma = sqrt(noiseVariance/2);

                noise = noiseSigma*( ...
                    randn(size(rxClean))+1j*randn(size(rxClean)));

                rxNoisy = rxClean+noise;

                timeData = prachFrontendProcessing( ...
                    rxNoisy,prachConfig,simConfig);

                [delayRaw,cfoEstimate,diagInfo] = ...
                    pairedZCFullLagCoarseEstimator( ...
                    timeData,prachConfig,optsNow);

                trial = trial+1;

                TrueCFO_Hz(trial) = currentCFO;
                TrueDelay_effSamples(trial) = currentDelayEff;
                TrueDelay_rawSamples(trial) = currentDelayRaw;
                SNR_dB(trial) = currentSNR;
                ChannelRealization(trial) = iRep;

                Branch1LagInteger(trial) = diagInfo.branch1LagInteger;
                Branch1LagSubBin(trial) = diagInfo.branch1LagSubBin;
                Branch2LagInteger(trial) = diagInfo.branch2LagInteger;
                Branch2LagSubBin(trial) = diagInfo.branch2LagSubBin;
                Branch1PeakValue(trial) = diagInfo.branch1PeakValue;
                Branch2PeakValue(trial) = diagInfo.branch2PeakValue;
                Branch2ToBranch1Peak_dB(trial) = ...
                    diagInfo.branch2ToBranch1Peak_dB;
                Branch1PeakToMedian_dB(trial) = ...
                    diagInfo.branch1PeakToMedian_dB;
                Branch2PeakToMedian_dB(trial) = ...
                    diagInfo.branch2PeakToMedian_dB;

                CFOEstimate_Hz(trial) = cfoEstimate;
                DelayEstimateRaw_effSamples(trial) = delayRaw;
                DelayEstimateCal_effSamples(trial) = ...
                    delayRaw-delayBiasCalibrationSamples;

                CFOError_Hz(trial) = cfoEstimate-currentCFO;
                DelayErrorRaw_effSamples(trial) = ...
                    delayRaw-currentDelayEff;
                DelayErrorCal_effSamples(trial) = ...
                    DelayEstimateCal_effSamples(trial)-currentDelayEff;

                % Expected branch coordinates under the V13.5 model:
                % d1 = delay + CFO/binSpacing
                % d2 = delay - CFO/binSpacing
                expectedD1 = currentDelayEff+currentCFO/cfoBinSpacing_Hz;
                expectedD2 = currentDelayEff-currentCFO/cfoBinSpacing_Hz;

                Branch1LagError_effSamples(trial) = ...
                    diagInfo.branch1LagSubBin-expectedD1;
                Branch2LagError_effSamples(trial) = ...
                    diagInfo.branch2LagSubBin-expectedD2;

                Peak1AtSearchBoundary(trial) = ...
                    diagInfo.branch1PeakAtSearchBoundary;
                Peak2AtSearchBoundary(trial) = ...
                    diagInfo.branch2PeakAtSearchBoundary;

                NoiseVariance(trial) = noiseVariance;
                RxSignalPower(trial) = pRxSignal;
            end
        end
    end

    fprintf('Completed TDL-D channel realization %d/%d (%d/%d trials).\n', ...
        iRep,nRep,trial,nTrials);
end

if trial ~= nTrials
    error('V13_6:TrialCountMismatch', ...
        'Generated %d trials but preallocated %d.',trial,nTrials);
end

%% 8. Raw results table

resultsTable = table( ...
    TrialID,ChannelRealization,TrueCFO_Hz, ...
    TrueDelay_effSamples,TrueDelay_rawSamples,SNR_dB, ...
    Branch1LagInteger,Branch1LagSubBin, ...
    Branch2LagInteger,Branch2LagSubBin, ...
    Branch1PeakValue,Branch2PeakValue, ...
    Branch2ToBranch1Peak_dB, ...
    Branch1PeakToMedian_dB,Branch2PeakToMedian_dB, ...
    CFOEstimate_Hz,DelayEstimateRaw_effSamples, ...
    DelayEstimateCal_effSamples,CFOError_Hz, ...
    DelayErrorRaw_effSamples,DelayErrorCal_effSamples, ...
    Branch1LagError_effSamples,Branch2LagError_effSamples, ...
    Peak1AtSearchBoundary,Peak2AtSearchBoundary, ...
    NoiseVariance,RxSignalPower);

rawFile = 'V13_6_TDLD_CFO80_Scatter_Raw.csv';
writetable(resultsTable,rawFile);

%% 9. Per-SNR performance metrics

metricSNR = strings(0,1);
metricCount = zeros(0,1);
CFO_RMSE_Hz = zeros(0,1);
CFO_MAE_Hz = zeros(0,1);
CFO_P95AbsError_Hz = zeros(0,1);
CFO_SuccessRate = zeros(0,1);
Delay_RMSE_Raw = zeros(0,1);
Delay_RMSE_Calibrated = zeros(0,1);
Delay_MAE_Calibrated = zeros(0,1);
Delay_SuccessRate = zeros(0,1);
Branch1Lag_RMSE = zeros(0,1);
Branch2Lag_RMSE = zeros(0,1);
BoundaryRate = zeros(0,1);
MeanBranchPeakRatio_dB = zeros(0,1);

for iSNR = 1:nSNR

    if isinf(cfg.snrList_dB(iSNR))
        mask = isinf(resultsTable.SNR_dB);
        label = "NoNoise";
    else
        mask = resultsTable.SNR_dB==cfg.snrList_dB(iSNR);
        label = string(cfg.snrList_dB(iSNR))+"dB";
    end

    eCFO = resultsTable.CFOError_Hz(mask);
    eDelayRaw = resultsTable.DelayErrorRaw_effSamples(mask);
    eDelayCal = resultsTable.DelayErrorCal_effSamples(mask);

    metricSNR(end+1,1) = label; %#ok<AGROW>
    metricCount(end+1,1) = sum(mask); %#ok<AGROW>

    CFO_RMSE_Hz(end+1,1) = sqrt(mean(eCFO.^2,'omitnan')); %#ok<AGROW>
    CFO_MAE_Hz(end+1,1) = mean(abs(eCFO),'omitnan'); %#ok<AGROW>
    CFO_P95AbsError_Hz(end+1,1) = ...
        localPercentile(abs(eCFO),95); %#ok<AGROW>
    CFO_SuccessRate(end+1,1) = mean( ...
        abs(eCFO)<=cfg.CFOAcceptanceThreshold_Hz,'omitnan'); %#ok<AGROW>

    Delay_RMSE_Raw(end+1,1) = ...
        sqrt(mean(eDelayRaw.^2,'omitnan')); %#ok<AGROW>
    Delay_RMSE_Calibrated(end+1,1) = ...
        sqrt(mean(eDelayCal.^2,'omitnan')); %#ok<AGROW>
    Delay_MAE_Calibrated(end+1,1) = ...
        mean(abs(eDelayCal),'omitnan'); %#ok<AGROW>
    Delay_SuccessRate(end+1,1) = mean( ...
        abs(eDelayCal)<=cfg.DelayAcceptanceThreshold_effSamples, ...
        'omitnan'); %#ok<AGROW>

    Branch1Lag_RMSE(end+1,1) = sqrt(mean( ...
        resultsTable.Branch1LagError_effSamples(mask).^2,'omitnan')); %#ok<AGROW>
    Branch2Lag_RMSE(end+1,1) = sqrt(mean( ...
        resultsTable.Branch2LagError_effSamples(mask).^2,'omitnan')); %#ok<AGROW>

    BoundaryRate(end+1,1) = mean( ...
        resultsTable.Peak1AtSearchBoundary(mask) | ...
        resultsTable.Peak2AtSearchBoundary(mask)); %#ok<AGROW>

    MeanBranchPeakRatio_dB(end+1,1) = ...
        mean(resultsTable.Branch2ToBranch1Peak_dB(mask),'omitnan'); %#ok<AGROW>
end

metricsTable = table(metricSNR,metricCount, ...
    CFO_RMSE_Hz,CFO_MAE_Hz,CFO_P95AbsError_Hz,CFO_SuccessRate, ...
    Delay_RMSE_Raw,Delay_RMSE_Calibrated, ...
    Delay_MAE_Calibrated,Delay_SuccessRate, ...
    Branch1Lag_RMSE,Branch2Lag_RMSE,BoundaryRate, ...
    MeanBranchPeakRatio_dB);

metricsFile = 'V13_6_TDLD_CFO80_Scatter_Metrics.csv';
writetable(metricsTable,metricsFile);

%% 10. Save complete experiment data

matFile = 'V13_6_TDLD_CFO80_Scatter.mat';

save(matFile,'cfg','simConfig','prachConfig','tdlInfo', ...
    'resultsTable','metricsTable', ...
    'calibrationDelayEstimates','delayBiasCalibrationSamples', ...
    'delayBiasMAD_samples','fsRaw','fsEffective','NIFFT','K','OS', ...
    'cfoBinSpacing_Hz','maxAbsCFO_Hz','maxCFOBins', ...
    'channelFilterDelay_raw','tdlMaxExcessDelay_raw', ...
    'tdlMaxExcessDelay_eff','lagSearchMin','lagSearchMax', ...
    'tailLength_raw','-v7.3');

%% 11. CFO true-versus-estimated scatter (80 displayed markers)

if cfg.makeFigures

    maskFit = resultsTable.SNR_dB == cfg.snrList_dB(1);
    xRaw = resultsTable.TrueCFO_Hz(maskFit);
    yRaw = resultsTable.CFOEstimate_Hz(maskFit);

    valid = isfinite(xRaw) & isfinite(yRaw);
    xRaw = xRaw(valid);
    yRaw = yRaw(valid);

    [xUnique,~,groupIdx] = unique(xRaw);
    meanEstimate = accumarray(groupIdx,yRaw,[],@(v) mean(v,'omitnan'));
    stdEstimate = accumarray(groupIdx,yRaw,[],@(v) std(v,0,'omitnan'));
    sampleCount = accumarray(groupIdx,1);

    if numel(xUnique) ~= 80
        error('V13_6:CFOFitPointCount', ...
            'Expected 80 unique CFO points, found %d.',numel(xUnique));
    end

    RMSE_AllRuns = sqrt(mean((yRaw-xRaw).^2,'omitnan'));
    RMSE_MeanPoints = sqrt(mean((meanEstimate-xUnique).^2,'omitnan'));

    fitPointsTable = table(xUnique,meanEstimate,stdEstimate,sampleCount, ...
        meanEstimate-xUnique, ...
        'VariableNames',{'TrueCFO_Hz','MeanEstimatedCFO_Hz', ...
        'StdEstimatedCFO_Hz','MonteCarloCount','MeanError_Hz'});
    fitPointsFile = 'V13_6_TDLD_CFO80_Scatter_FitPoints.csv';
    writetable(fitPointsTable,fitPointsFile);

    % Plot the 80 empirical mean points only. No regression/smoothing curve.
    % Draw the ideal y=x reference first so it stays behind the data markers.
    figFit = figure('Color','w','Position',[100 100 1050 760]);
    refMin = min([xUnique; meanEstimate]);
    refMax = max([xUnique; meanEstimate]);
    plot([refMin refMax],[refMin refMax],'k--','LineWidth',1.3);
    hold on;
    plot(xUnique,meanEstimate,'b+','MarkerSize',7.0,'LineWidth',1.2);
    grid on; box on;
    xlabel('True CFO (Hz)');
    ylabel('Mean estimated CFO (Hz)');
    title(sprintf('CFO estimation: 80 empirical points, SNR = %g dB', ...
        cfg.snrList_dB(1)));
    legend('Ideal y = x','Monte-Carlo mean (3 runs/point)', ...
        'Location','northwest');

    statText = sprintf(['RMSE (all 240 runs) = %.3f Hz\\n', ...
        'RMSE (80 mean points) = %.3f Hz'], ...
        RMSE_AllRuns,RMSE_MeanPoints);
    text(0.03,0.97,statText,'Units','normalized', ...
        'VerticalAlignment','top','BackgroundColor','w','Margin',6);
    axis equal;
    xlim([refMin refMax]);
    ylim([refMin refMax]);
    hold off;

    savefig(figFit,'V13_6_TDLD_CFO80_Scatter_Accuracy.fig');
    exportgraphics(figFit,'V13_6_TDLD_CFO80_Scatter_Accuracy.png','Resolution',250);
end

%% 12. Optional diagnostics

if cfg.makeDiagnosticFigures

    figDiag = figure('Color','w','Position',[100 100 1300 800]);
    td = tiledlayout(figDiag,2,2, ...
        'TileSpacing','compact','Padding','compact');

    nexttile(td);
    scatter(resultsTable.Branch1LagError_effSamples, ...
        resultsTable.Branch2LagError_effSamples,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('Branch 1 lag-coordinate error (samples)');
    ylabel('Branch 2 lag-coordinate error (samples)');
    title('Branch peak coordinate errors');

    nexttile(td);
    histogram(resultsTable.Branch2ToBranch1Peak_dB,40);
    grid on;
    xlabel('Branch 2 / Branch 1 peak ratio (dB)');
    ylabel('Count');
    title('Paired branch peak-strength distribution');

    nexttile(td);
    scatter(resultsTable.TrueDelay_effSamples, ...
        resultsTable.DelayErrorRaw_effSamples,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True propagation delay (samples)');
    ylabel('Raw delay error (samples)');
    title('Raw delay error versus propagation delay');

    nexttile(td);
    scatter(resultsTable.TrueCFO_Hz, ...
        resultsTable.CFOError_Hz,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True CFO (Hz)');
    ylabel('CFO error (Hz)');
    title('CFO error versus true CFO');

    title(td,'V13.6 TDL-D diagnostic plots');
    savefig(figDiag,'V13_6_TDLD_CFO80_Scatter_Diagnostics.fig');
    exportgraphics(figDiag, ...
        'V13_6_TDLD_CFO80_Scatter_Diagnostics.png','Resolution',250);
end

%% 13. Summary

fprintf('\n============================================================\n');
fprintf('V13.6 TDL-D 80-POINT SCATTER VALIDATION FINISHED\n');
fprintf('============================================================\n');
fprintf('Evaluation channel realizations : %d\n',nRep);
fprintf('Distinct CFO / delay points      : %d / %d\n',nCFO,nDelay);
fprintf('SNR used for 80-point scatter       : %g dB\n',cfg.snrList_dB(1));
fprintf('Independent calibration channels: %d\n', ...
    cfg.nCalibrationChannelRealizations);
fprintf('Total evaluation trials         : %d\n',height(resultsTable));
fprintf('TDL-D delay spread               : %.3f ns RMS\n', ...
    cfg.tdlDelaySpread_s*1e9);
fprintf('TDL-D maximum Doppler            : %.3f Hz (static-tap stage)\n', ...
    cfg.tdlMaximumDopplerShift_Hz);
fprintf('Railway common CFO range         : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));
fprintf('TDL-D delay calibration bias     : %+.9f samples\n', ...
    delayBiasCalibrationSamples);
fprintf('TDL-D delay calibration MAD      : %.9f samples\n', ...
    delayBiasMAD_samples);

fprintf('\nPer-SNR metrics:\n');
disp(metricsTable);

fprintf('\nOutput files:\n');
fprintf('  %s\n',rawFile);
fprintf('  %s\n',metricsFile);
fprintf('  %s\n',matFile);
if cfg.makeFigures
    fprintf('  V13_6_TDLD_CFO80_Scatter_Performance.png / .fig\n');
    fprintf('  V13_6_TDLD_CFO80_Scatter_Accuracy.png / .fig\n');
end
if cfg.makeDiagnosticFigures
    fprintf('  V13_6_TDLD_CFO80_Scatter_Diagnostics.png / .fig\n');
end

fprintf('\nThis is the static-tap TDL-D multipath isolation stage.\n');
fprintf('The time-varying Doppler-spread stage is intentionally separate.\n');
fprintf('No legacy recovery estimator or PDP generator was called.\n');
fprintf('No transmit waveform energy normalization was applied.\n');

%% ============================================================
% LOCAL FUNCTION: construct a seeded independent TDL-D channel
%% ============================================================

function tdlChan = makeTDLChannelLocal(cfg,fsRaw,seed)

tdlChan = nrTDLChannel;
tdlChan.DelayProfile = cfg.tdlDelayProfile;
tdlChan.DelaySpread = cfg.tdlDelaySpread_s;
tdlChan.MaximumDopplerShift = cfg.tdlMaximumDopplerShift_Hz;
tdlChan.SampleRate = fsRaw;
tdlChan.NumTransmitAntennas = cfg.nTx;
tdlChan.NumReceiveAntennas = cfg.nRx;
tdlChan.RandomStream = 'mt19937ar with seed';
tdlChan.Seed = seed;
tdlChan.ChannelFiltering = true;

% Keep the standard TDL-D first-tap LOS/Rician specification.
if isprop(tdlChan,'KFactorScaling')
    tdlChan.KFactorScaling = false;
end

end

%% ============================================================
% LOCAL FUNCTION: filter one burst through one seeded TDL-D realization
%% ============================================================

function rxFaded = generateTDLRealizationLocal( ...
    tdlChan,txBurst,tailLength,channelFilterDelay)

txWithTail = [txBurst; ...
    zeros(tailLength,size(txBurst,2),'like',txBurst)];

rxFull = tdlChan(txWithTail);

first = channelFilterDelay+1;
last = first+size(txBurst,1)-1;

if size(rxFull,1)<last
    error('V13_6:TDLOutputTooShort', ...
        ['TDL-D output has %d samples; %d are required after ', ...
         'channel-filter-delay compensation.'],size(rxFull,1),last);
end

rxFaded = rxFull(first:last,:);

end

%% ============================================================
% LOCAL FUNCTION: fractional raw-sample delay for multiple antennas
%% ============================================================

function y = fractionalDelayRawLocal(x,delaySamples,guardSamples)

if ~isnumeric(x) || isempty(x) || ndims(x)~=2
    error('V13_6:InvalidDelayInput', ...
        'Fractional-delay input must be a nonempty 2-D matrix.');
end

if ~isfinite(delaySamples) || delaySamples<0
    error('V13_6:InvalidDelay', ...
        'Delay must be finite and nonnegative.');
end

nAnt = size(x,2);
yCell = cell(1,nAnt);

for ia = 1:nAnt
    yCell{ia} = fractionalDelayRawVectorLocal( ...
        x(:,ia),delaySamples,guardSamples);
end

outLen = cellfun(@numel,yCell);
if any(outLen~=outLen(1))
    error('V13_6:DelayLengthMismatch', ...
        'Fractional-delay output lengths differ across antennas.');
end

y = complex(zeros(outLen(1),nAnt,'like',x));
for ia = 1:nAnt
    y(:,ia) = yCell{ia};
end

end

function y = fractionalDelayRawVectorLocal(x,delaySamples,guardSamples)

x = x(:);

integerDelay = floor(delaySamples);
fractionalPart = delaySamples-integerDelay;
xInteger = [zeros(integerDelay,1,'like',x);x];

if fractionalPart<1e-12
    y = xInteger;
    return;
end

guardSamples = max(guardSamples,32);

xPad = [zeros(guardSamples,1,'like',xInteger); ...
        xInteger; ...
        zeros(guardSamples,1,'like',xInteger)];

nFFT = 2^nextpow2(2*numel(xPad));

if mod(nFFT,2)==0
    k = [0:nFFT/2,-nFFT/2+1:-1].';
else
    k = [0:(nFFT-1)/2,-(nFFT-1)/2:-1].';
end

X = fft(xPad,nFFT);
Y = ifft(X.*exp(-1j*2*pi*k*fractionalPart/nFFT));

first = guardSamples+1;
last = first+numel(x)+integerDelay;
y = Y(first:last);

end

%% ============================================================
% LOCAL FUNCTION: shared gains, consistent with V13.5
%% ============================================================

function [rxGain,refGain] = deriveSharedGainsLocal(timeData,prachConfig,nRx)

N = prachConfig.bwpInfo.NFFT;
timeFlat = reshape(timeData,nRx,[]);
rxPower = mean(abs(timeFlat).^2,2);

ref = ifft(prachConfig.fdRootSeq,N)*sqrt(N);
refPower = mean(abs(ref).^2);

if any(~isfinite(rxPower)) || any(rxPower<=0) || ...
        ~isfinite(refPower) || refPower<=0
    error('V13_6:InvalidGainSource', ...
        'Unable to derive finite positive shared gains.');
end

rxGain = 1./sqrt(rxPower);
refGain = 1/sqrt(refPower);

end

%% ============================================================
% LOCAL FUNCTION: percentile without requiring prctile
%% ============================================================

function q = localPercentile(x,p)

x = sort(x(isfinite(x)));
if isempty(x)
    q = NaN;
    return;
end

if p<=0
    q = x(1);
elseif p>=100
    q = x(end);
else
    pos = 1+(numel(x)-1)*p/100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo==hi
        q = x(lo);
    else
        q = x(lo)+(pos-lo)*(x(hi)-x(lo));
    end
end

end

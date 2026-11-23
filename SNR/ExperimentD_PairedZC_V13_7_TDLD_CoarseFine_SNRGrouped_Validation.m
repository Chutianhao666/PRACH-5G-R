clc;
clear;
close all;
clear functions;
rehash;

%% ============================================================
% V13.7 - TDL-D FULL-LAG COARSE + MP-ACR FINE ESTIMATION
%         INTEGRATED VALIDATION
%
% Purpose:
%   Connect the unchanged V13.6 full-lag/sub-bin coarse estimator to the
%   current HSR_MP_ACR fine estimator and compare coarse/final accuracy.
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
%   - HSR_MP_ACR performs local multipath-aware CFO/timing refinement
%   - No time_freq_recovery / V6 / prachPowerDelayProfileGen
%   - No transmit waveform energy normalization
%   - No per-algorithm received-power normalization
%
% TDL-D channel-filter delay:
%   The integer implementation delay reported by info(nrTDLChannel) is
%   removed before applying the known railway propagation delay. The
%   physical TDL path delays remain in the received waveform.
%
% Coarse/fine interface:
%   The raw coarse delay is used to circularly advance the received window.
%   The raw coarse CFO is passed to HSR_MP_ACR for residual-frequency
%   refinement. No external delay-bias calibration is applied before MP-ACR;
%   the fine estimator must correct the coarse residual itself.
%   Both first-path and strongest-path fine delays are retained. The
%   first-path estimate is compared with the known geometric propagation
%   delay; strongest-path delay is retained as a multipath diagnostic.
%
% Run modes:
%   smoke   : 1 channel realization  (20 trials/SNR)
%   pilot   : 5 channel realizations (100 trials/SNR)
%   full    : 20 channel realizations (400 trials/SNR)
%   snr200  : 10 channel realizations (200 trials/SNR)
%
% SNR study used in this customized version:
%   SNR = +9 : -1 : -15 dB
%   5 CFO points x 4 delay points x 10 channel realizations
%   => 200 Monte-Carlo trials per SNR
%   => 25 SNR points, 5,000 total trials
%
% Accuracy figures:
%   The script groups the SNRs three at a time:
%   [+9,+8,+7], [+6,+5,+4], ..., [-12,-13,-14], [-15]
%   Each grouped figure uses one row per SNR and two columns:
%       left  = CFO true-vs-estimated scatter (coarse/fine)
%       right = delay true-vs-estimated scatter (coarse/fine)
%
% Required:
%   5G Toolbox (nrTDLChannel)
%   Communications/DSP functionality used by the supplied front end
%
% Outputs:
%   V13_7_TDLD_CoarseFine_SNRGrouped_Raw.csv
%   V13_7_TDLD_CoarseFine_SNRGrouped_Metrics.csv
%   V13_7_TDLD_CoarseFine.mat
%   V13_7_TDLD_CoarseFine_SNRGrouped_Performance.png/.fig
%   V13_7_TDLD_CoarseFine_SNRGrouped_Accuracy.png/.fig
%% ============================================================

fprintf('\n============================================================\n');
fprintf('V13.7 TDL-D COARSE + FINE SNR-GROUPED VALIDATION\n');
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

% Geometry Doppler +/- residual oscillator CFO of up to 300 Hz.
cfg.CFOList_Hz = unique(sort([ ...
    -fd-300, -fd, 0, fd, fd+300]));

% Absolute propagation-delay test points on the effective sample grid.
% 17.460638 samples corresponds to approximately 1 km propagation range
% at the current effective sampling rate. TDL-D excess delays are additional.
cfg.delayList_effSamples = [0, 0.4, 4, 17.460638];

% Include noise-free and low-SNR stress points.
cfg.snrList_dB = 9:-1:-15;

% TDL-D configuration. DelaySpread is RMS delay spread, not maximum delay.
cfg.tdlDelayProfile = 'TDL-D';
cfg.tdlDelaySpread_s = 100e-9;

% Stage A isolates multipath with a static channel realization per burst.
% Common railway Doppler is explicitly included in CFOList_Hz above.
cfg.tdlMaximumDopplerShift_Hz = 0;

% Select a short smoke test first, then use 'snr200' for the requested
% 200-trials-per-SNR configuration.
% smoke = 1 channel (20 trials/SNR)
% pilot = 5 channel realizations (100 trials/SNR)
% full  = 20 channel realizations (400 trials/SNR)
% snr200 = 10 channel realizations (200 trials/SNR)
cfg.runMode = 'snr200';  % 'smoke' | 'pilot' | 'full' | 'snr200'
switch lower(cfg.runMode)
    case 'smoke'
        cfg.nChannelRealizations = 1;
    case 'pilot'
        cfg.nChannelRealizations = 5;
    case 'full'
        cfg.nChannelRealizations = 20;
    case 'snr200'
        cfg.nChannelRealizations = 10;
    otherwise
        error('V13_7:InvalidRunMode','Unknown run mode: %s',cfg.runMode);
end

cfg.randomSeed = 20261121;
cfg.fractionalDelayGuardSamples = 256;

% The search guard is set after reading the actual scaled TDL-D path delays.
cfg.searchGuardBinsExtra = 2;

% Output controls: retain only the two principal result figures by default.
cfg.makeFigures = true;
cfg.makeDiagnosticFigures = false;

rng(cfg.randomSeed,'twister');

%% 1. Verify active functions and required toolbox

requiredFunctions = { ...
    'prachSigGen', ...
    'prachFrontendProcessing', ...
    'pairedZCFullLagCoarseEstimator', ...
    'HSR_MP_ACR', ...
    'railwayDopplerModel'};

fprintf('\n===== Active MATLAB function paths =====\n');
for k = 1:numel(requiredFunctions)
    p = which(requiredFunctions{k});
    if isempty(p)
        error('V13_7:MissingFunction', ...
            'Required function not found: %s',requiredFunctions{k});
    end
    fprintf('%-40s %s\n',requiredFunctions{k},p);
end

if exist('nrTDLChannel','file') ~= 2 && exist('nrTDLChannel','class') ~= 8
    error('V13_7:Missing5GToolbox', ...
        ['nrTDLChannel was not found. This V13.7 script requires ', ...
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
    error('V13_7:UnexpectedNrep', ...
        'This validation expects PRACH F0 with Nrep=1.');
end
if abs(NFFTraw/NIFFT-K)>1e-9
    error('V13_7:SamplingRatioMismatch', ...
        'Raw/effective FFT ratio does not match prachConfig.K.');
end

rawBurstLength = NCPraw+NFFTraw*Nrep;
if size(txSig,1)<rawBurstLength
    error('V13_7:ShortTX', ...
        'Generated TX waveform is shorter than the active PRACH burst.');
end

txBurst = txSig(1:rawBurstLength,1);
txBurst = txBurst(:);

% Build the paired-ZC local reference through the SAME PRACH front end.
% This preserves the downconversion/decimation/FIR response seen by RX.
txRefData = prachFrontendProcessing( ...
    repmat(txSig,1,cfg.nRx),prachConfig,simConfig);
txRef = reshape(txRefData(1,:,1),[],1);
if numel(txRef)~=NIFFT || any(~isfinite(txRef))
    error('V13_7:InvalidFineReference', ...
        'The front-end-processed paired-ZC reference is invalid.');
end

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
    error('V13_7:UnexpectedTDLInfo', ...
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
    warning('V13_7:DelayExceedsCP', ...
        ['Absolute plus TDL excess delay (%.2f raw samples) is not ', ...
         'smaller than the PRACH CP (%d raw samples).'], ...
         maxTotalDelayRaw,NCPraw);
end

%% 4. Coarse and fine estimator options

% Keep the V13.6 full-lag coarse estimator unchanged.
estimatorOpts = struct();
estimatorOpts.MaxAbsCFO_Hz = maxAbsCFO_Hz;
estimatorOpts.MinDelaySamples = min(cfg.delayList_effSamples);
estimatorOpts.MaxDelaySamples = max(cfg.delayList_effSamples);
estimatorOpts.SearchGuardBins = cfg.searchGuardBins;
estimatorOpts.ReturnCorrelations = false;

% MP-ACR fine-estimator configuration.
% The measured V13.6 coarse CFO error is below 35 Hz in the full 2,000-trial set;
% +/-150 Hz therefore leaves a comfortable residual-search margin.
% The coarse delay residual is well inside +/-1.5 effective samples.
% A +/-4-sample path-search window covers the dominant early TDL-D taps.
fineOpts = struct();
fineOpts.inputCFOCompensated = false;
fineOpts.freqRadiusHz = 150;
fineOpts.delayRadiusSamples = 1.5;
fineOpts.maxIter = 4;
fineOpts.freqGridN = 81;
fineOpts.delayGridN = 61;
fineOpts.maxPaths = 4;
fineOpts.pathSearchSpanSamples = 4.0;
fineOpts.fractionalDelayGridStep = 0.05;
fineOpts.pathRelThreshold = 0.12;
fineOpts.minPathSeparationSamples = 0.20;
fineOpts.ridge = 1e-8;
fineOpts.pathOffsetsSamples = [];
fineOpts.delayInitSamples = [];

fprintf('\n===== Fine estimator configuration =====\n');
fprintf('Fine estimator              : HSR_MP_ACR\n');
fprintf('Residual CFO search radius  : +/- %.1f Hz\n',fineOpts.freqRadiusHz);
fprintf('Residual delay search radius: +/- %.2f samples\n', ...
    fineOpts.delayRadiusSamples);
fprintf('Path search span            : +/- %.2f samples\n', ...
    fineOpts.pathSearchSpanSamples);
fprintf('Maximum fitted paths        : %d\n',fineOpts.maxPaths);
fprintf('Path grid step              : %.3f samples\n', ...
    fineOpts.fractionalDelayGridStep);
fprintf('Run mode / channel epochs   : %s / %d\n', ...
    cfg.runMode,cfg.nChannelRealizations);

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

CFOEstimateCoarse_Hz = nan(nTrials,1);
DelayEstimateCoarse_effSamples = nan(nTrials,1);
CFOEstimateFine_Hz = nan(nTrials,1);
DelayEstimateFineFirst_effSamples = nan(nTrials,1);
DelayEstimateFineStrongest_effSamples = nan(nTrials,1);
CFOErrorCoarse_Hz = nan(nTrials,1);
CFOErrorFine_Hz = nan(nTrials,1);
DelayErrorCoarse_effSamples = nan(nTrials,1);
DelayErrorFineFirst_effSamples = nan(nTrials,1);
DelayErrorFineStrongest_effSamples = nan(nTrials,1);
FineDeltaFreq_Hz = nan(nTrials,1);
FineDeltaTauSamples = nan(nTrials,1);
FineNumPaths = nan(nTrials,1);
FineNormalizedResidual = nan(nTrials,1);

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
                error('V13_7:EmptyActiveSupport', ...
                    'The TDL-D received waveform has no active samples.');
            end

            pRxSignal = mean(abs(rxClean(activeIdx,:)).^2,'all');

            if ~isfinite(pRxSignal) || pRxSignal<=0
                error('V13_7:InvalidRXPower', ...
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

                % Convert front-end output to samples-by-Rx for MP-ACR.
                rxWin = reshape(timeData,cfg.nRx,[]).';
                if size(rxWin,1)~=NIFFT || size(rxWin,2)~=cfg.nRx
                    error('V13_7:InvalidFineInputWindow', ...
                        'Unexpected effective-rate receive-window dimensions.');
                end

                % Coarse timing alignment: advance by the RAW coarse delay.
                % Do not apply the V13.6 post-hoc delay-bias calibration here;
                % MP-ACR must refine the actual coarse residual itself.
                rxAligned = advanceCircular(rxWin,delayRaw);

                fineOut = HSR_MP_ACR( ...
                    rxAligned,txRef,fsEffective, ...
                    delayRaw,cfoEstimate,fineOpts);

                requiredFineFields = { ...
                    'freqFineHz','tauFirstPathSamples', ...
                    'tauStrongestPathSamples','deltaFreqHz', ...
                    'deltaTauSamples','numPaths','normalizedResidual'};
                for iField = 1:numel(requiredFineFields)
                    if ~isfield(fineOut,requiredFineFields{iField})
                        error('V13_7:MissingFineOutput', ...
                            'HSR_MP_ACR output is missing field %s.', ...
                            requiredFineFields{iField});
                    end
                end
                if any(~isfinite([fineOut.freqFineHz, ...
                        fineOut.tauFirstPathSamples, ...
                        fineOut.tauStrongestPathSamples, ...
                        fineOut.normalizedResidual]))
                    error('V13_7:InvalidFineOutput', ...
                        'HSR_MP_ACR returned a non-finite estimate.');
                end

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

                CFOEstimateCoarse_Hz(trial) = cfoEstimate;
                DelayEstimateCoarse_effSamples(trial) = delayRaw;
                CFOEstimateFine_Hz(trial) = fineOut.freqFineHz;
                DelayEstimateFineFirst_effSamples(trial) = ...
                    fineOut.tauFirstPathSamples;
                DelayEstimateFineStrongest_effSamples(trial) = ...
                    fineOut.tauStrongestPathSamples;

                CFOErrorCoarse_Hz(trial) = cfoEstimate-currentCFO;
                CFOErrorFine_Hz(trial) = fineOut.freqFineHz-currentCFO;
                DelayErrorCoarse_effSamples(trial) = ...
                    delayRaw-currentDelayEff;
                DelayErrorFineFirst_effSamples(trial) = ...
                    fineOut.tauFirstPathSamples-currentDelayEff;
                DelayErrorFineStrongest_effSamples(trial) = ...
                    fineOut.tauStrongestPathSamples-currentDelayEff;

                FineDeltaFreq_Hz(trial) = fineOut.deltaFreqHz;
                FineDeltaTauSamples(trial) = fineOut.deltaTauSamples;
                FineNumPaths(trial) = fineOut.numPaths;
                FineNormalizedResidual(trial) = fineOut.normalizedResidual;

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
    error('V13_7:TrialCountMismatch', ...
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
    CFOEstimateCoarse_Hz,DelayEstimateCoarse_effSamples, ...
    CFOEstimateFine_Hz,DelayEstimateFineFirst_effSamples, ...
    DelayEstimateFineStrongest_effSamples, ...
    CFOErrorCoarse_Hz,CFOErrorFine_Hz, ...
    DelayErrorCoarse_effSamples,DelayErrorFineFirst_effSamples, ...
    DelayErrorFineStrongest_effSamples, ...
    FineDeltaFreq_Hz,FineDeltaTauSamples,FineNumPaths, ...
    FineNormalizedResidual, ...
    Branch1LagError_effSamples,Branch2LagError_effSamples, ...
    Peak1AtSearchBoundary,Peak2AtSearchBoundary, ...
    NoiseVariance,RxSignalPower);

rawFile = 'V13_7_TDLD_CoarseFine_SNRGrouped_Raw.csv';
writetable(resultsTable,rawFile);

%% 9. Per-SNR coarse/fine performance metrics

metricSNR = strings(0,1);
metricCount = zeros(0,1);
CFO_RMSE_Coarse_Hz = zeros(0,1);
CFO_RMSE_Fine_Hz = zeros(0,1);
CFO_MAE_Coarse_Hz = zeros(0,1);
CFO_MAE_Fine_Hz = zeros(0,1);
Delay_RMSE_CoarseRaw_effSamples = zeros(0,1);
Delay_RMSE_FineFirst_effSamples = zeros(0,1);
Delay_RMSE_FineStrongest_effSamples = zeros(0,1);
Delay_MAE_CoarseRaw_effSamples = zeros(0,1);
Delay_MAE_FineFirst_effSamples = zeros(0,1);
FineResidual_Mean = zeros(0,1);
FineNumPaths_Mean = zeros(0,1);
FineEstimatorSuccessRate = zeros(0,1);
BoundaryRate = zeros(0,1);

for iSNR = 1:nSNR

    if isinf(cfg.snrList_dB(iSNR))
        mask = isinf(resultsTable.SNR_dB);
        label = "NoNoise";
    else
        mask = resultsTable.SNR_dB==cfg.snrList_dB(iSNR);
        label = string(cfg.snrList_dB(iSNR))+"dB";
    end

    eCFOCoarse = resultsTable.CFOErrorCoarse_Hz(mask);
    eCFOFine = resultsTable.CFOErrorFine_Hz(mask);
    eDelayCoarse = resultsTable.DelayErrorCoarse_effSamples(mask);
    eDelayFineFirst = resultsTable.DelayErrorFineFirst_effSamples(mask);
    eDelayFineStrongest = resultsTable.DelayErrorFineStrongest_effSamples(mask);

    metricSNR(end+1,1) = label; %#ok<AGROW>
    metricCount(end+1,1) = sum(mask); %#ok<AGROW>

    CFO_RMSE_Coarse_Hz(end+1,1) = ...
        sqrt(mean(eCFOCoarse.^2,'omitnan')); %#ok<AGROW>
    CFO_RMSE_Fine_Hz(end+1,1) = ...
        sqrt(mean(eCFOFine.^2,'omitnan')); %#ok<AGROW>
    CFO_MAE_Coarse_Hz(end+1,1) = ...
        mean(abs(eCFOCoarse),'omitnan'); %#ok<AGROW>
    CFO_MAE_Fine_Hz(end+1,1) = ...
        mean(abs(eCFOFine),'omitnan'); %#ok<AGROW>

    Delay_RMSE_CoarseRaw_effSamples(end+1,1) = ...
        sqrt(mean(eDelayCoarse.^2,'omitnan')); %#ok<AGROW>
    Delay_RMSE_FineFirst_effSamples(end+1,1) = ...
        sqrt(mean(eDelayFineFirst.^2,'omitnan')); %#ok<AGROW>
    Delay_RMSE_FineStrongest_effSamples(end+1,1) = ...
        sqrt(mean(eDelayFineStrongest.^2,'omitnan')); %#ok<AGROW>
    Delay_MAE_CoarseRaw_effSamples(end+1,1) = ...
        mean(abs(eDelayCoarse),'omitnan'); %#ok<AGROW>
    Delay_MAE_FineFirst_effSamples(end+1,1) = ...
        mean(abs(eDelayFineFirst),'omitnan'); %#ok<AGROW>

    FineResidual_Mean(end+1,1) = ...
        mean(resultsTable.FineNormalizedResidual(mask),'omitnan'); %#ok<AGROW>
    FineNumPaths_Mean(end+1,1) = ...
        mean(resultsTable.FineNumPaths(mask),'omitnan'); %#ok<AGROW>
    FineEstimatorSuccessRate(end+1,1) = mean( ...
        isfinite(resultsTable.CFOEstimateFine_Hz(mask)) & ...
        isfinite(resultsTable.DelayEstimateFineFirst_effSamples(mask))); %#ok<AGROW>

    BoundaryRate(end+1,1) = mean( ...
        resultsTable.Peak1AtSearchBoundary(mask) | ...
        resultsTable.Peak2AtSearchBoundary(mask)); %#ok<AGROW>
end

metricsTable = table(metricSNR,metricCount, ...
    CFO_RMSE_Coarse_Hz,CFO_RMSE_Fine_Hz, ...
    CFO_MAE_Coarse_Hz,CFO_MAE_Fine_Hz, ...
    Delay_RMSE_CoarseRaw_effSamples, ...
    Delay_RMSE_FineFirst_effSamples, ...
    Delay_RMSE_FineStrongest_effSamples, ...
    Delay_MAE_CoarseRaw_effSamples,Delay_MAE_FineFirst_effSamples, ...
    FineResidual_Mean,FineNumPaths_Mean,FineEstimatorSuccessRate, ...
    BoundaryRate);

metricsFile = 'V13_7_TDLD_CoarseFine_SNRGrouped_Metrics.csv';
writetable(metricsTable,metricsFile);

%% 10. Save complete experiment data

matFile = 'V13_7_TDLD_CoarseFine.mat';

save(matFile,'cfg','simConfig','prachConfig','tdlInfo', ...
    'estimatorOpts','fineOpts','txRef', ...
    'resultsTable','metricsTable', ...
    'fsRaw','fsEffective','NIFFT','K','OS', ...
    'cfoBinSpacing_Hz','maxAbsCFO_Hz','maxCFOBins', ...
    'channelFilterDelay_raw','tdlMaxExcessDelay_raw', ...
    'tdlMaxExcessDelay_eff','lagSearchMin','lagSearchMax', ...
    'tailLength_raw','-v7.3');

%% 11. Main coarse-versus-fine performance figures

if cfg.makeFigures

    finiteRows = ~(metricsTable.metricSNR=="NoNoise");
    finiteLabels = erase(metricsTable.metricSNR(finiteRows),"dB");
    finiteSNR = str2double(finiteLabels);
    [finiteSNR,order] = sort(finiteSNR);
    metricRows = find(finiteRows);
    metricRows = metricRows(order);

    figPerf = figure('Color','w','Position',[80 80 1300 560]);
    tl = tiledlayout(figPerf,1,2, ...
        'TileSpacing','compact','Padding','compact');

    nexttile(tl);
    plot(finiteSNR,metricsTable.CFO_RMSE_Coarse_Hz(metricRows),'-o', ...
        'LineWidth',1.4,'MarkerSize',6);
    hold on;
    plot(finiteSNR,metricsTable.CFO_RMSE_Fine_Hz(metricRows),'-s', ...
        'LineWidth',1.4,'MarkerSize',6);
    grid on;
    xlabel('SNR (dB)');
    ylabel('CFO RMSE (Hz)');
    title('CFO: coarse versus MP-ACR fine');
    legend('Full-lag coarse','MP-ACR fine','Location','best');
    hold off;

    nexttile(tl);
    plot(finiteSNR, ...
        metricsTable.Delay_RMSE_CoarseRaw_effSamples(metricRows),'-o', ...
        'LineWidth',1.4,'MarkerSize',6);
    hold on;
    plot(finiteSNR, ...
        metricsTable.Delay_RMSE_FineFirst_effSamples(metricRows),'-s', ...
        'LineWidth',1.4,'MarkerSize',6);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Delay RMSE (effective samples)');
    title('Propagation delay: coarse versus fine first path');
    legend('Full-lag coarse','MP-ACR fine first path','Location','best');
    hold off;

    title(tl,sprintf( ...
        ['V13.7 TDL-D coarse/fine SNR sweep, %.0f ns, static taps, ', ...
         'SNR = %+d:%+d dB, 200 trials/SNR'], ...
        cfg.tdlDelaySpread_s*1e9,max(cfg.snrList_dB),min(cfg.snrList_dB)));
    savefig(figPerf,'V13_7_TDLD_CoarseFine_SNRGrouped_Performance.fig');
    exportgraphics(figPerf, ...
        'V13_7_TDLD_CoarseFine_SNRGrouped_Performance.png','Resolution',250);

    % Group the SNRs three at a time for accuracy scatter figures.
    groupSize = 3;
    nGroups = ceil(numel(cfg.snrList_dB)/groupSize);

    for iGroup = 1:nGroups

        firstIdx = (iGroup-1)*groupSize + 1;
        lastIdx = min(iGroup*groupSize,numel(cfg.snrList_dB));
        snrGroup = cfg.snrList_dB(firstIdx:lastIdx);
        nRows = numel(snrGroup);

        figAcc = figure('Color','w','Position',[100 80 1250 330*nRows]);
        tla = tiledlayout(figAcc,nRows,2, ...
            'TileSpacing','compact','Padding','compact');

        for ir = 1:nRows
            thisSNR = snrGroup(ir);
            maskNow = resultsTable.SNR_dB == thisSNR;

            % -------- CFO panel --------
            nexttile(tla);
            scatter(resultsTable.TrueCFO_Hz(maskNow), ...
                resultsTable.CFOEstimateCoarse_Hz(maskNow),20,'filled');
            hold on;
            scatter(resultsTable.TrueCFO_Hz(maskNow), ...
                resultsTable.CFOEstimateFine_Hz(maskNow),22,'filled');
            cfoAll = [resultsTable.TrueCFO_Hz(maskNow); ...
                      resultsTable.CFOEstimateCoarse_Hz(maskNow); ...
                      resultsTable.CFOEstimateFine_Hz(maskNow)];
            lo = min(cfoAll); hi = max(cfoAll);
            plot([lo hi],[lo hi],'k--','LineWidth',1.2);
            grid on;
            xlabel('True CFO (Hz)');
            ylabel('Estimated CFO (Hz)');
            title(sprintf('CFO estimate at %+d dB',thisSNR));
            legend('Coarse','Fine','y=x','Location','best');
            hold off;

            % -------- Delay panel --------
            nexttile(tla);
            scatter(resultsTable.TrueDelay_effSamples(maskNow), ...
                resultsTable.DelayEstimateCoarse_effSamples(maskNow),20,'filled');
            hold on;
            scatter(resultsTable.TrueDelay_effSamples(maskNow), ...
                resultsTable.DelayEstimateFineFirst_effSamples(maskNow),22,'filled');
            dAll = [resultsTable.TrueDelay_effSamples(maskNow); ...
                    resultsTable.DelayEstimateCoarse_effSamples(maskNow); ...
                    resultsTable.DelayEstimateFineFirst_effSamples(maskNow)];
            lo = min(dAll); hi = max(dAll);
            plot([lo hi],[lo hi],'k--','LineWidth',1.2);
            grid on;
            xlabel('True propagation delay (effective samples)');
            ylabel('Estimated delay (effective samples)');
            title(sprintf('First-path delay estimate at %+d dB',thisSNR));
            legend('Coarse','Fine first path','y=x','Location','best');
            hold off;
        end

        title(tla,sprintf( ...
            ['V13.7 TDL-D full-lag coarse + MP-ACR fine accuracy, ', ...
             'SNR group %s to %s'], ...
            snrTagLocal(snrGroup(1)),snrTagLocal(snrGroup(end))));

        figBase = sprintf('V13_7_TDLD_CoarseFine_SNRGrouped_Accuracy_%s_to_%s', ...
            snrTagLocal(snrGroup(1)),snrTagLocal(snrGroup(end)));
        savefig(figAcc,[figBase '.fig']);
        exportgraphics(figAcc,[figBase '.png'],'Resolution',250);
    end
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
        resultsTable.DelayErrorCoarse_effSamples,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True propagation delay (samples)');
    ylabel('Raw delay error (samples)');
    title('Raw delay error versus propagation delay');

    nexttile(td);
    scatter(resultsTable.TrueCFO_Hz, ...
        resultsTable.CFOErrorCoarse_Hz,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True CFO (Hz)');
    ylabel('CFO error (Hz)');
    title('CFO error versus true CFO');

    title(td,'V13.7 TDL-D coarse/fine diagnostic plots');
    savefig(figDiag,'V13_7_TDLD_CoarseFine_SNRGrouped_Diagnostics.fig');
    exportgraphics(figDiag, ...
        'V13_7_TDLD_CoarseFine_SNRGrouped_Diagnostics.png','Resolution',250);
end

%% 13. Summary

fprintf('\n============================================================\n');
fprintf('V13.7 TDL-D COARSE + FINE SNR-GROUPED VALIDATION FINISHED\n');
fprintf('============================================================\n');
fprintf('Run mode                        : %s\n',cfg.runMode);
fprintf('Evaluation channel realizations : %d\n',nRep);
fprintf('Total evaluation trials         : %d\n',height(resultsTable));
fprintf('Trials per SNR                  : %d\n',nRep*nCFO*nDelay);
fprintf('TDL-D delay spread               : %.3f ns RMS\n', ...
    cfg.tdlDelaySpread_s*1e9);
fprintf('TDL-D maximum Doppler            : %.3f Hz (static-tap stage)\n', ...
    cfg.tdlMaximumDopplerShift_Hz);
fprintf('Railway common CFO range         : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));

fprintf('\nPer-SNR metrics:\n');
disp(metricsTable);

fprintf('\nOutput files:\n');
fprintf('  %s\n',rawFile);
fprintf('  %s\n',metricsFile);
fprintf('  %s\n',matFile);
if cfg.makeFigures
    fprintf('  V13_7_TDLD_CoarseFine_SNRGrouped_Performance.png / .fig\n');
    fprintf('  V13_7_TDLD_CoarseFine_SNRGrouped_Accuracy_<SNR1>_to_<SNR2>.png / .fig\n');
end
if cfg.makeDiagnosticFigures
    fprintf('  V13_7_TDLD_CoarseFine_SNRGrouped_Diagnostics.png / .fig\n');
end

fprintf('\nThis is the static-tap TDL-D multipath isolation stage.\n');
fprintf('The time-varying Doppler-spread stage is intentionally separate.\n');
fprintf('No legacy recovery estimator or PDP generator was called.\n');
fprintf('No transmit waveform energy normalization was applied.\n');

%% ============================================================
% LOCAL FUNCTION: fractional circular timing advance
%% ============================================================

function y = advanceCircular(x,tau)

if isvector(x)
    x = x(:);
end

N = size(x,1);
k = (0:N-1).';
k(k>floor(N/2)) = k(k>floor(N/2))-N;

% Positive tau advances the received sequence by tau effective samples.
y = ifft(bsxfun(@times,fft(x,[],1), ...
    exp(1j*2*pi*k*tau/N)),[],1);

end

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
    error('V13_7:TDLOutputTooShort', ...
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
    error('V13_7:InvalidDelayInput', ...
        'Fractional-delay input must be a nonempty 2-D matrix.');
end

if ~isfinite(delaySamples) || delaySamples<0
    error('V13_7:InvalidDelay', ...
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
    error('V13_7:DelayLengthMismatch', ...
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
    error('V13_7:InvalidGainSource', ...
        'Unable to derive finite positive shared gains.');
end

rxGain = 1./sqrt(rxPower);
refGain = 1/sqrt(refPower);

end


%% ============================================================
% LOCAL FUNCTION: compact SNR tag for titles/file names
%% ============================================================

function tag = snrTagLocal(v)

if ~isfinite(v)
    tag = 'NoNoise';
elseif v >= 0
    tag = sprintf('p%d',round(v));
else
    tag = sprintf('m%d',abs(round(v)));
end

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

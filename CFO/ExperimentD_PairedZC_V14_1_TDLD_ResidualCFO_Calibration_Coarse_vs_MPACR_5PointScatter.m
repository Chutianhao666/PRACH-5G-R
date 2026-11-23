clc;
clear;
close all;
clear functions;
rehash;

%% ============================================================
% V14.1 - GEOMETRY-REFERENCED RANDOM RESIDUAL CFO CALIBRATION
%         FULL-LAG COARSE + MP-ACR FINE
%
% Purpose:
%   Validate coarse/fine CFO calibration over a residual/random CFO range
%   of -300 to +300 Hz while retaining the existing 500 km/h railway model.
%
% Fixed conditions:
%   SNR                         = +5 dB
%   Additional propagation delay = 0 effective samples
%   Geometry Doppler            = retained from railwayDopplerModel
%   TDL-D RMS DelaySpread       = 100 ns
%   TDL-D MaximumDopplerShift   = 0 Hz (static taps)
%
% Additional residual-CFO sweep:
%   -300:30:+300 Hz, i.e. 21 equally spaced offset levels.
%   The extra term represents residual/random frequency offset on top of
%   the existing railway-geometry Doppler. It is NOT another full railway
%   Doppler term and must not be confused with the 500 km/h geometric shift.
%   With the current geometry Doppler of about +967.768 Hz, the total
%   applied CFO therefore spans about +667.768 to +1267.768 Hz.
%
% Monte-Carlo:
%   5 independent TDL-D channel realizations per residual-CFO level.
%   Each channel realization is reused across the full residual-CFO sweep;
%   each point receives a fresh noise realization. Coarse/fine estimators
%   see the same channel and noise realization (paired comparison).
%   21 frequency levels x 5 paired trials = 105 trials.
%
% Compared CFO estimators:
%   1) Full-lag paired-ZC coarse estimator with 3-point sub-bin interpolation
%   2) HSR_MP_ACR fine estimator
%
% Important scope:
%   TDL-D taps remain static. The swept Doppler is applied as a common
%   complex exponential after the channel, i.e. a common CFO-like Doppler
%   component. This experiment calibrates common frequency offset; it does
%   NOT yet model independent time-varying Doppler on every TDL path.
%
% Required:
%   5G Toolbox (nrTDLChannel)
%   Communications/DSP functionality used by the supplied PRACH front end
%
% Outputs:
%   V14_1_TDLD_ResidualCFOCalibration_Raw.csv
%   V14_1_TDLD_ResidualCFOCalibration_Metrics.csv
%   V14_1_TDLD_ResidualCFOCalibration.mat
%   V14_1_TDLD_ResidualCFOCalibration_Calibration.png/.fig
%   V14_1_TDLD_ResidualCFOCalibration_Error.png/.fig
%% ============================================================

fprintf('\n============================================================\n');
fprintf('V14.1 TDL-D RESIDUAL-CFO CALIBRATION VALIDATION\n');
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

% Residual/random CFO range selected for this sensitivity experiment.
% Scan -300 to +300 Hz in 30 Hz steps (21 levels). The railway-geometry
% Doppler is retained separately and added exactly once to obtain total CFO.
cfg.extraDopplerMax_Hz = 300;
cfg.extraDopplerStep_Hz = 30;
cfg.extraDopplerList_Hz = ...
    -cfg.extraDopplerMax_Hz:cfg.extraDopplerStep_Hz:cfg.extraDopplerMax_Hz;
cfg.CFOList_Hz = cfg.geometryDoppler_Hz + cfg.extraDopplerList_Hz;

% No additional propagation delay in this experiment.
cfg.extraDelay_effSamples = 0;

% Fixed SNR, matching the V13.9 delay-calibration experiment.
cfg.snrList_dB = 5;

% TDL-D configuration. DelaySpread is RMS delay spread, not maximum delay.
cfg.tdlDelayProfile = 'TDL-D';
cfg.tdlDelaySpread_s = 100e-9;

% Stage A isolates multipath with a static channel realization per burst.
% Common railway Doppler is explicitly included in CFOList_Hz above.
cfg.tdlMaximumDopplerShift_Hz = 0;

% smoke = 2; five = 5 (five visible raw points per residual-CFO level);
% pilot = 20; full = 200.
% Five independent TDL-D channel realizations are used at every
% residual-CFO level. Coarse and fine estimators share the same
% channel/noise trial (paired comparison).
cfg.runMode = 'pilot';  % 'smoke' | 'five' | 'pilot' | 'full'
switch lower(cfg.runMode)
    case 'smoke'
        cfg.nChannelRealizations = 2;
    case 'five'
        cfg.nChannelRealizations = 5;
    case 'pilot'
        cfg.nChannelRealizations = 20;
    case 'full'
        cfg.nChannelRealizations = 200;
    otherwise
        error('V14_1:InvalidRunMode','Unknown run mode: %s',cfg.runMode);
end

cfg.randomSeed = 20261122;
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
        error('V14_1:MissingFunction', ...
            'Required function not found: %s',requiredFunctions{k});
    end
    fprintf('%-40s %s\n',requiredFunctions{k},p);
end

if exist('nrTDLChannel','file') ~= 2 && exist('nrTDLChannel','class') ~= 8
    error('V14_1:Missing5GToolbox', ...
        ['nrTDLChannel was not found. This V14.1 script requires ', ...
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

% Keep the geometric propagation delay fixed; no additional delay is swept.
cfg.geometryDelay_effSamples = ...
    cfg.propagationDelayAtT0_s*fsEffective;
cfg.delayList_effSamples = cfg.geometryDelay_effSamples;

if numel(cfg.extraDopplerList_Hz) ~= 21
    error('V14_1:UnexpectedDopplerGrid', ...
        'The configured residual-CFO sweep must contain exactly 21 points.');
end

if Nrep ~= 1
    error('V14_1:UnexpectedNrep', ...
        'This validation expects PRACH F0 with Nrep=1.');
end
if abs(NFFTraw/NIFFT-K)>1e-9
    error('V14_1:SamplingRatioMismatch', ...
        'Raw/effective FFT ratio does not match prachConfig.K.');
end

rawBurstLength = NCPraw+NFFTraw*Nrep;
if size(txSig,1)<rawBurstLength
    error('V14_1:ShortTX', ...
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
    error('V14_1:InvalidFineReference', ...
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
    error('V14_1:UnexpectedTDLInfo', ...
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
fprintf('Total applied CFO range    : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));
fprintf('Extra residual-CFO sweep   : [%+.1f, %+.1f] Hz (%d points)\n', ...
    min(cfg.extraDopplerList_Hz),max(cfg.extraDopplerList_Hz), ...
    numel(cfg.extraDopplerList_Hz));
fprintf('Fixed geometric delay      : %.6f effective samples\n', ...
    cfg.geometryDelay_effSamples);
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
fprintf('Geometric delay baseline    : %.6f effective samples\n', ...
    cfg.geometryDelay_effSamples);
fprintf('Additional delay            : %.3f effective samples (fixed)\n', ...
    cfg.extraDelay_effSamples);
fprintf('Extra residual-CFO range     : [%+.1f, %+.1f] Hz\n', ...
    min(cfg.extraDopplerList_Hz),max(cfg.extraDopplerList_Hz));
fprintf('Applied total CFO range      : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));
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
    warning('V14_1:DelayExceedsCP', ...
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
fprintf('Trials per residual-CFO level: %d\n',cfg.nChannelRealizations);
fprintf('Total expected trials       : %d\n', ...
    cfg.nChannelRealizations*numel(cfg.extraDopplerList_Hz));

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
TrueExtraDoppler_Hz = nan(nTrials,1);
GeometryDoppler_Hz = nan(nTrials,1);
GeometryDelay_effSamples = nan(nTrials,1);
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
DelayEstimateFineStrongest_effSamples = nan(nTrials,1);
CFOErrorCoarse_Hz = nan(nTrials,1);
CFOErrorFine_Hz = nan(nTrials,1);
DelayErrorCoarse_effSamples = nan(nTrials,1);
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
    % All additional-Doppler points reuse this same channel realization.
    evalChan = makeTDLChannelLocal( ...
        cfg,fsRaw,cfg.randomSeed+10000+iRep);
    rxFaded = generateTDLRealizationLocal( ...
        evalChan,txBurst,tailLength_raw,channelFilterDelay_raw);

    for iCFO = 1:nCFO

        currentExtraDoppler = cfg.extraDopplerList_Hz(iCFO);
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
                error('V14_1:EmptyActiveSupport', ...
                    'The TDL-D received waveform has no active samples.');
            end

            pRxSignal = mean(abs(rxClean(activeIdx,:)).^2,'all');

            if ~isfinite(pRxSignal) || pRxSignal<=0
                error('V14_1:InvalidRXPower', ...
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

                % Convert front-end output to samples-by-Rx for the fine estimator.
                rxWin = reshape(timeData,cfg.nRx,[]).';
                if size(rxWin,1)~=NIFFT || size(rxWin,2)~=cfg.nRx
                    error('V14_1:InvalidFineInputWindow', ...
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
                    'freqFineHz','tauStrongestPathSamples','deltaFreqHz', ...
                    'deltaTauSamples','numPaths','normalizedResidual'};
                for iField = 1:numel(requiredFineFields)
                    if ~isfield(fineOut,requiredFineFields{iField})
                        error('V14_1:MissingFineOutput', ...
                            'HSR_MP_ACR output is missing field %s.', ...
                            requiredFineFields{iField});
                    end
                end
                if any(~isfinite([fineOut.freqFineHz, ...
                        fineOut.tauStrongestPathSamples, ...
                        fineOut.normalizedResidual]))
                    error('V14_1:InvalidFineOutput', ...
                        'HSR_MP_ACR returned a non-finite estimate.');
                end

                trial = trial+1;

                TrueCFO_Hz(trial) = currentCFO;
                TrueExtraDoppler_Hz(trial) = currentExtraDoppler;
                GeometryDoppler_Hz(trial) = cfg.geometryDoppler_Hz;
                TrueDelay_effSamples(trial) = currentDelayEff;
                TrueDelay_rawSamples(trial) = currentDelayRaw;
                GeometryDelay_effSamples(trial) = cfg.geometryDelay_effSamples;
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
                DelayEstimateFineStrongest_effSamples(trial) = ...
                    fineOut.tauStrongestPathSamples;

                CFOErrorCoarse_Hz(trial) = cfoEstimate-currentCFO;
                CFOErrorFine_Hz(trial) = fineOut.freqFineHz-currentCFO;
                DelayErrorCoarse_effSamples(trial) = ...
                    delayRaw-currentDelayEff;
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
    error('V14_1:TrialCountMismatch', ...
        'Generated %d trials but preallocated %d.',trial,nTrials);
end

%% 8. Raw results table

CFOEstimateCoarseExtra_Hz = ...
    CFOEstimateCoarse_Hz-cfg.geometryDoppler_Hz;
CFOEstimateFineExtra_Hz = ...
    CFOEstimateFine_Hz-cfg.geometryDoppler_Hz;

resultsTable = table( ...
    TrialID,ChannelRealization,TrueExtraDoppler_Hz, ...
    GeometryDoppler_Hz,TrueCFO_Hz, ...
    GeometryDelay_effSamples,TrueDelay_effSamples,TrueDelay_rawSamples,SNR_dB, ...
    Branch1LagInteger,Branch1LagSubBin, ...
    Branch2LagInteger,Branch2LagSubBin, ...
    Branch1PeakValue,Branch2PeakValue, ...
    Branch2ToBranch1Peak_dB, ...
    Branch1PeakToMedian_dB,Branch2PeakToMedian_dB, ...
    CFOEstimateCoarse_Hz,CFOEstimateCoarseExtra_Hz, ...
    DelayEstimateCoarse_effSamples, ...
    CFOEstimateFine_Hz,CFOEstimateFineExtra_Hz, ...
    DelayEstimateFineStrongest_effSamples, ...
    CFOErrorCoarse_Hz,CFOErrorFine_Hz, ...
    DelayErrorCoarse_effSamples,DelayErrorFineStrongest_effSamples, ...
    FineDeltaFreq_Hz,FineDeltaTauSamples,FineNumPaths, ...
    FineNormalizedResidual, ...
    Branch1LagError_effSamples,Branch2LagError_effSamples, ...
    Peak1AtSearchBoundary,Peak2AtSearchBoundary, ...
    NoiseVariance,RxSignalPower);

rawFile = 'V14_1_TDLD_ResidualCFOCalibration_Raw.csv';
writetable(resultsTable,rawFile);

%% 9. Per-residual-CFO calibration metrics

ExtraDoppler_Hz = cfg.extraDopplerList_Hz(:);
TrueTotalCFO_Hz = cfg.geometryDoppler_Hz+ExtraDoppler_Hz;
metricCount = zeros(numel(ExtraDoppler_Hz),1);

Coarse_Bias_Hz = nan(size(ExtraDoppler_Hz));
Fine_Bias_Hz = nan(size(ExtraDoppler_Hz));
Coarse_RMSE_Hz = nan(size(ExtraDoppler_Hz));
Fine_RMSE_Hz = nan(size(ExtraDoppler_Hz));
Coarse_MAE_Hz = nan(size(ExtraDoppler_Hz));
Fine_MAE_Hz = nan(size(ExtraDoppler_Hz));
Coarse_Std_Hz = nan(size(ExtraDoppler_Hz));
Fine_Std_Hz = nan(size(ExtraDoppler_Hz));
Coarse_P95AbsError_Hz = nan(size(ExtraDoppler_Hz));
Fine_P95AbsError_Hz = nan(size(ExtraDoppler_Hz));
Fine_RMSE_Reduction_pct = nan(size(ExtraDoppler_Hz));
FineEstimatorSuccessRate = nan(size(ExtraDoppler_Hz));
BoundaryRate = nan(size(ExtraDoppler_Hz));

for iDoppler = 1:numel(ExtraDoppler_Hz)

    mask = resultsTable.TrueExtraDoppler_Hz == ExtraDoppler_Hz(iDoppler);

    eCoarse = resultsTable.CFOErrorCoarse_Hz(mask);
    eFine = resultsTable.CFOErrorFine_Hz(mask);

    metricCount(iDoppler) = sum(mask);

    Coarse_Bias_Hz(iDoppler) = mean(eCoarse,'omitnan');
    Fine_Bias_Hz(iDoppler) = mean(eFine,'omitnan');

    Coarse_RMSE_Hz(iDoppler) = sqrt(mean(eCoarse.^2,'omitnan'));
    Fine_RMSE_Hz(iDoppler) = sqrt(mean(eFine.^2,'omitnan'));

    Coarse_MAE_Hz(iDoppler) = mean(abs(eCoarse),'omitnan');
    Fine_MAE_Hz(iDoppler) = mean(abs(eFine),'omitnan');

    Coarse_Std_Hz(iDoppler) = std(eCoarse,'omitnan');
    Fine_Std_Hz(iDoppler) = std(eFine,'omitnan');

    Coarse_P95AbsError_Hz(iDoppler) = ...
        localPercentile(abs(eCoarse),95);
    Fine_P95AbsError_Hz(iDoppler) = ...
        localPercentile(abs(eFine),95);

    if isfinite(Coarse_RMSE_Hz(iDoppler)) && ...
            Coarse_RMSE_Hz(iDoppler)>0
        Fine_RMSE_Reduction_pct(iDoppler) = 100*(1- ...
            Fine_RMSE_Hz(iDoppler)/Coarse_RMSE_Hz(iDoppler));
    end

    FineEstimatorSuccessRate(iDoppler) = mean( ...
        isfinite(resultsTable.CFOEstimateFine_Hz(mask)));

    BoundaryRate(iDoppler) = mean( ...
        resultsTable.Peak1AtSearchBoundary(mask) | ...
        resultsTable.Peak2AtSearchBoundary(mask));
end

metricsTable = table( ...
    ExtraDoppler_Hz,TrueTotalCFO_Hz,metricCount, ...
    Coarse_Bias_Hz,Fine_Bias_Hz, ...
    Coarse_RMSE_Hz,Fine_RMSE_Hz, ...
    Coarse_MAE_Hz,Fine_MAE_Hz, ...
    Coarse_Std_Hz,Fine_Std_Hz, ...
    Coarse_P95AbsError_Hz,Fine_P95AbsError_Hz, ...
    Fine_RMSE_Reduction_pct, ...
    FineEstimatorSuccessRate,BoundaryRate);

metricsFile = 'V14_1_TDLD_ResidualCFOCalibration_Metrics.csv';
writetable(metricsTable,metricsFile);

%% 10. Save complete experiment data

matFile = 'V14_1_TDLD_ResidualCFOCalibration.mat';

save(matFile,'cfg','simConfig','prachConfig','tdlInfo', ...
    'estimatorOpts','fineOpts','txRef', ...
    'resultsTable','metricsTable', ...
    'fsRaw','fsEffective','NIFFT','K','OS', ...
    'cfoBinSpacing_Hz','maxAbsCFO_Hz','maxCFOBins', ...
    'channelFilterDelay_raw','tdlMaxExcessDelay_raw', ...
    'tdlMaxExcessDelay_eff','lagSearchMin','lagSearchMax', ...
    'tailLength_raw','-v7.3');

%% 11. Main residual-CFO calibration figures

if cfg.makeFigures

    xExtra = metricsTable.ExtraDoppler_Hz;

    % Figure 1: every Monte-Carlo CFO estimate is shown as an individual dot.
    % No line connects individual trials; only the ideal y=x line is drawn.
    figCal = figure('Color','w','Position',[80 80 1320 820]);
    ax = axes(figCal);
    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    plot(ax,[min(xExtra) max(xExtra)],[min(xExtra) max(xExtra)], ...
        'k--','LineWidth',1.8,'DisplayName','Ideal calibration: y=x');

    coarseColor = [0.0000 0.4470 0.7410];
    fineColor   = [0.8500 0.3250 0.0980];

    % Filled-circle markers scale consistently with scatter SizeData.
    coarseScatterSize = 58;
    fineScatterSize   = 62;

    allEstimatedExtra = [];
    uniqueX = unique(xExtra);
    if numel(uniqueX)>1
        dopplerStep = min(diff(uniqueX));
    else
        dopplerStep = 1;
    end

    for iPlotDoppler = 1:numel(xExtra)
        pointMask = resultsTable.TrueExtraDoppler_Hz == xExtra(iPlotDoppler);
        pointRows = resultsTable(pointMask,:);
        pointRows = sortrows(pointRows,'ChannelRealization');
        nPoint = height(pointRows);

        if nPoint == 0
            continue;
        end

        % Display-only horizontal offsets separate the five repeated trials.
        % The y-values are the actual, unmodified estimated frequency offsets.
        coarseOffset = linspace(-0.14,-0.035,nPoint).'*dopplerStep;
        fineOffset   = linspace( 0.035, 0.14,nPoint).'*dopplerStep;
        xCoarsePlot = xExtra(iPlotDoppler)+coarseOffset;
        xFinePlot   = xExtra(iPlotDoppler)+fineOffset;

        yCoarsePlot = pointRows.CFOEstimateCoarseExtra_Hz;
        yFinePlot = pointRows.CFOEstimateFineExtra_Hz;
        allEstimatedExtra = [allEstimatedExtra; yCoarsePlot; yFinePlot]; %#ok<AGROW>

        if iPlotDoppler == 1
            scatter(ax,xCoarsePlot,yCoarsePlot,coarseScatterSize, ...
                coarseColor,'.','filled', ...
                'DisplayName','Full-lag coarse + 3-point parabola');
            scatter(ax,xFinePlot,yFinePlot,fineScatterSize, ...
                fineColor,'.','filled','DisplayName','MP-ACR fine');
        else
            scatter(ax,xCoarsePlot,yCoarsePlot,coarseScatterSize, ...
                coarseColor,'.','filled','HandleVisibility','off');
            scatter(ax,xFinePlot,yFinePlot,fineScatterSize, ...
                fineColor,'.','filled','HandleVisibility','off');
        end
    end

    xlabel(ax,'True additional residual CFO (Hz)');
    ylabel(ax,'Estimated additional residual CFO (Hz)');
    title(ax,sprintf([ ...
        'Individual residual-CFO estimates: SNR = %+g dB, additional delay = %.1f samples\\n', ...
        'TDL-D RMS delay spread = %.0f ns; %d paired trials per frequency point'], ...
        cfg.snrList_dB(1),cfg.extraDelay_effSamples,cfg.tdlDelaySpread_s*1e9, ...
        cfg.nChannelRealizations));
    text(ax,0.015,0.985, ...
        'Horizontal offsets are for visibility only; estimated values are unchanged.', ...
        'Units','normalized','VerticalAlignment','top', ...
        'FontSize',9,'Color',[0.25 0.25 0.25], ...
        'BackgroundColor','w','Margin',3);
    legend(ax,'Location','northwest');
    xPad = max(30,0.025*(max(xExtra)-min(xExtra)));
    xlim(ax,[min(xExtra)-xPad max(xExtra)+xPad]);
    if ~isempty(allEstimatedExtra)
        yLow = min([xExtra(:); allEstimatedExtra(:)]);
        yHigh = max([xExtra(:); allEstimatedExtra(:)]);
        yPad = max(30,0.04*(yHigh-yLow));
        ylim(ax,[yLow-yPad yHigh+yPad]);
    end
    hold(ax,'off');

    savefig(figCal,'V14_1_TDLD_ResidualCFOCalibration_Calibration.fig');
    exportgraphics(figCal, ...
        'V14_1_TDLD_ResidualCFOCalibration_Calibration.png','Resolution',300);

    % Figure 2: RMSE and systematic bias versus additional Doppler offset.
    figErr = figure('Color','w','Position',[80 80 1450 620]);
    te = tiledlayout(figErr,1,2,'TileSpacing','compact','Padding','compact');

    nexttile(te);
    hold on;
    plot(xExtra,metricsTable.Coarse_RMSE_Hz,'-o', ...
        'LineWidth',1.6,'MarkerSize',6);
    plot(xExtra,metricsTable.Fine_RMSE_Hz,'-s', ...
        'LineWidth',1.7,'MarkerSize',6);
    grid on;
    xlabel('True additional residual CFO (Hz)');
    ylabel('CFO RMSE (Hz)');
    title('CFO estimation error versus residual CFO');
    legend('Full-lag coarse + 3-point parabola', ...
        'MP-ACR fine','Location','best');
    hold off;

    nexttile(te);
    hold on;
    yline(0,'k--','LineWidth',1.2,'DisplayName','Zero bias');
    plot(xExtra,metricsTable.Coarse_Bias_Hz,'-o', ...
        'LineWidth',1.6,'MarkerSize',6, ...
        'DisplayName','Full-lag coarse + 3-point parabola');
    plot(xExtra,metricsTable.Fine_Bias_Hz,'-s', ...
        'LineWidth',1.7,'MarkerSize',6,'DisplayName','MP-ACR fine');
    grid on;
    xlabel('True additional residual CFO (Hz)');
    ylabel('CFO bias (estimate minus truth, Hz)');
    title('Systematic CFO calibration bias');
    legend('Location','best');
    hold off;

    title(te,'V14.1 TDL-D residual-CFO calibration');
    savefig(figErr,'V14_1_TDLD_ResidualCFOCalibration_Error.fig');
    exportgraphics(figErr, ...
        'V14_1_TDLD_ResidualCFOCalibration_Error.png','Resolution',300);
end

%% 12. Optional diagnostics

if cfg.makeDiagnosticFigures

    figDiag = figure('Color','w','Position',[100 100 1300 800]);
    td = tiledlayout(figDiag,2,2, ...
        'TileSpacing','compact','Padding','compact');

    nexttile(td);
    scatter(resultsTable.Branch1LagError_effSamples, ...
        resultsTable.Branch2LagError_effSamples,8, ...
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
        resultsTable.DelayErrorCoarse_effSamples,8, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True propagation delay (samples)');
    ylabel('Raw delay error (samples)');
    title('Fixed-delay residual under the residual-CFO sweep');

    nexttile(td);
    scatter(resultsTable.TrueCFO_Hz, ...
        resultsTable.CFOErrorCoarse_Hz,8, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True CFO (Hz)');
    ylabel('CFO error (Hz)');
    title('CFO error versus true total CFO');

    title(td,'V14.1 TDL-D coarse/fine Doppler diagnostic plots');
    savefig(figDiag,'V14_1_TDLD_ResidualCFOCalibration_Diagnostics.fig');
    exportgraphics(figDiag, ...
        'V14_1_TDLD_ResidualCFOCalibration_Diagnostics.png','Resolution',250);
end

%% 13. Summary

fprintf('\\n============================================================\\n');
fprintf('V14.1 TDL-D RESIDUAL-CFO CALIBRATION FINISHED\\n');
fprintf('============================================================\\n');
fprintf('Run mode                         : %s\\n',cfg.runMode);
fprintf('Independent channel realizations: %d\\n',nRep);
fprintf('Residual-CFO test levels          : %d\\n',nCFO);
fprintf('Total evaluation trials           : %d\\n',height(resultsTable));
fprintf('Fixed SNR                         : %+g dB\\n',cfg.snrList_dB(1));
fprintf('Fixed additional delay            : %.3f samples\\n', ...
    cfg.extraDelay_effSamples);
fprintf('Applied geometry Doppler          : %+.6f Hz\\n', ...
    cfg.geometryDoppler_Hz);
fprintf('Extra residual-CFO range          : [%+.1f, %+.1f] Hz\\n', ...
    min(cfg.extraDopplerList_Hz),max(cfg.extraDopplerList_Hz));
fprintf('Total applied CFO range            : [%+.3f, %+.3f] Hz\\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));
fprintf('TDL-D RMS delay spread             : %.3f ns\\n', ...
    cfg.tdlDelaySpread_s*1e9);
fprintf('TDL-D maximum Doppler              : %.3f Hz (static taps)\\n', ...
    cfg.tdlMaximumDopplerShift_Hz);

fprintf('\\nPer-residual-CFO calibration metrics:\\n');
disp(metricsTable);

fprintf('\\nOutput files:\\n');
fprintf('  %s\\n',rawFile);
fprintf('  %s\\n',metricsFile);
fprintf('  %s\\n',matFile);
if cfg.makeFigures
    fprintf('  V14_1_TDLD_ResidualCFOCalibration_Calibration.png / .fig\\n');
    fprintf('  V14_1_TDLD_ResidualCFOCalibration_Error.png / .fig\\n');
end
if cfg.makeDiagnosticFigures
    fprintf('  V14_1_TDLD_ResidualCFOCalibration_Diagnostics.png / .fig\\n');
end

fprintf('\\nThe experiment sweeps a common additional Doppler/CFO term on top of railway geometry.\\n');
fprintf('Compared methods: full-lag coarse + 3-point parabola, and MP-ACR fine CFO estimation.\\n');
fprintf('TDL-D taps remain static; this is not a per-path time-varying Doppler-spread test.\\n');

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
    error('V14_1:TDLOutputTooShort', ...
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
    error('V14_1:InvalidDelayInput', ...
        'Fractional-delay input must be a nonempty 2-D matrix.');
end

if ~isfinite(delaySamples) || delaySamples<0
    error('V14_1:InvalidDelay', ...
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
    error('V14_1:DelayLengthMismatch', ...
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
    error('V14_1:InvalidGainSource', ...
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

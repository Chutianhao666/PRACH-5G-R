clc;
clear;
close all;
clear functions;
rehash;

%% ============================================================
% V13.8 - TDL-D + OPPOSING-TRAIN INTERFERENCE
%         FULL-LAG COARSE + MP-ACR FINE ESTIMATION
%
% Purpose:
%   Extend the validated V13.7 paired-ZC coarse/fine receiver by adding
%   one opposing-train interferer. The receiver chain is unchanged:
%
%       paired-ZC RX -> V13.6 Full-Lag coarse -> HSR_MP_ACR fine
%
% Signal / geometry:
%   Desired train:
%       Paired-conjugate ZC, RootIndex = 22, unnormalized.
%       Position = -1000 m, speed = +500 km/h, track distance = 100 m.
%
%   Opposing interferer:
%       Traditional ZC, RootIndex = 0, unnormalized.
%       Same longitudinal snapshot, speed = -500 km/h,
%       adjacent track distance = 105 m (5 m separation).
%       Imposed relative PRACH arrival offset = +8 effective samples.
%
% Channel:
%   Independent TDL-D channels for target and interferer.
%   DelaySpread = 100 ns; MaximumDopplerShift = 0 Hz (static taps).
%   Railway LOS Doppler is applied externally:
%       target CFO = current V13.7 CFO test value
%       interferer CFO = target CFO + (fd_interferer - fd_target)
%
% Power convention (retained from the user's A2 setup):
%   SNR and SIR are referenced to the received traditional-ZC Root-22
%   target power over the target active PRACH support.
%   The paired target waveform and Root-0 interferer are NOT normalized.
%   Thus the paired waveform retains its original energy advantage.
%   SIR = Inf means no interferer; SNR = Inf means no AWGN.
%
% Important:
%   - V13.7 coarse estimator and HSR_MP_ACR are not modified.
%   - The coarse search window remains the target acquisition window;
%     it is not widened to deliberately search for the interferer.
%   - The +8-sample relative interferer delay is the imposed A2 offset.
%     The small geometry-derived propagation-delay difference is logged
%     but not separately added, matching the earlier A2 experiment.
%   - FineTargetLock / FineInterfererLock mean that the estimate lies
%     within the configured MP-ACR local capture window (±150 Hz,
%     ±1.5 effective samples); they are association diagnostics, not P_D.
%   - Per-trial MP-ACR exceptions/non-finite outputs are recorded rather
%     than aborting the complete Monte-Carlo run.
%   - This remains a static-TDL-tap stage, not time-varying TDL fading.
%   - No legacy recovery estimator / PDP detector / threshold is used.
%
% Run modes:
%   smoke : 1 channel realization (8,820 trials)
%   pilot : 5 channel realizations (44,100 trials)
%   full  : 20 channel realizations (176,400 trials)
%
% Outputs:
%   V13_8_TDLD_OpposingTrain_CoarseFine_Raw.csv
%   V13_8_TDLD_OpposingTrain_CoarseFine_Metrics.csv
%   V13_8_TDLD_OpposingTrain_CoarseFine.mat
%   V13_8_TDLD_OpposingTrain_CFO_CoarseFine_3D.png/.fig
%   V13_8_TDLD_OpposingTrain_Delay_Strongest_CoarseFine_3D.png/.fig
%% ============================================================

fprintf('\n============================================================\n');
fprintf('V13.8 TDL-D OPPOSING-TRAIN COARSE + FINE VALIDATION\n');
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
cfg.interfererRootIndex = 0;

% Desired train: same V13.7 geometry.
cfg.TrainPosition_m = -1000;
cfg.TrackDistance_m = 100;
cfg.BaseHeight_m = 30;
cfg.TrainAntennaHeight_m = 5;
cfg.c_mps = 299792458;

% Opposing train: same longitudinal snapshot, opposite direction,
% adjacent parallel track separated by 5 m.
cfg.InterfererPosition_m = cfg.TrainPosition_m;
cfg.InterfererSpeed_kmh = -cfg.MobileSpeed_kmh;
cfg.TrackSeparation_m = 5;
cfg.InterfererTrackDistance_m = cfg.TrackDistance_m + cfg.TrackSeparation_m;

% Target railway LOS Doppler at t=0.
[fdTarget0,rangeTarget0_m,delayTarget0_s] = railwayDopplerModel( ...
    cfg.TrainPosition_m,cfg.MobileSpeed_kmh/3.6, ...
    cfg.TrackDistance_m, ...
    cfg.BaseHeight_m-cfg.TrainAntennaHeight_m, ...
    cfg.CarrierFrequency_Hz,0);

% Opposing-train LOS Doppler at t=0.
[fdInterferer0,rangeInterferer0_m,delayInterferer0_s] = railwayDopplerModel( ...
    cfg.InterfererPosition_m,cfg.InterfererSpeed_kmh/3.6, ...
    cfg.InterfererTrackDistance_m, ...
    cfg.BaseHeight_m-cfg.TrainAntennaHeight_m, ...
    cfg.CarrierFrequency_Hz,0);

cfg.geometryDopplerTarget_Hz = fdTarget0;
cfg.geometryDopplerInterferer_Hz = fdInterferer0;
cfg.relativeGeometryDoppler_Hz = fdInterferer0-fdTarget0;
cfg.rangeTargetAtT0_m = rangeTarget0_m;
cfg.rangeInterfererAtT0_m = rangeInterferer0_m;
cfg.propagationDelayTargetAtT0_s = delayTarget0_s;
cfg.propagationDelayInterfererAtT0_s = delayInterferer0_s;
cfg.geometryPropagationDelayDifference_s = delayInterferer0_s-delayTarget0_s;

% Keep the exact V13.7 target CFO test points and absolute target delays.
cfg.CFOList_Hz = unique(sort([ ...
    -fdTarget0-300, -fdTarget0, 0, fdTarget0, fdTarget0+300]));
cfg.delayList_effSamples = [0, 0.4, 4, 17.460638];

% Requested dense 1 dB-step SNR/SIR robustness grid; no Inf cases.
% SIR is defined against conventional Root-22 target reference power.
cfg.snrList_dB = -15:3:9;
cfg.sirList_dB = 10:-5:-15;
cfg.interferenceDelayEffSamples = 8;

cfg.tdlDelayProfile = 'TDL-D';
cfg.tdlDelaySpread_s = 100e-9;
cfg.tdlMaximumDopplerShift_Hz = 0;

% Select a short smoke test first, then pilot/full.
cfg.runMode = 'full';  % 'smoke' | 'pilot' | 'full'
switch lower(cfg.runMode)
    case 'smoke'
        cfg.nChannelRealizations = 1;
    case 'pilot'
        cfg.nChannelRealizations = 5;
    case 'full'
        cfg.nChannelRealizations = 20;
    otherwise
        error('V13_8:InvalidRunMode','Unknown run mode: %s',cfg.runMode);
end

cfg.randomSeed = 20261121;
cfg.fractionalDelayGuardSamples = 256;
cfg.searchGuardBinsExtra = 2;

% Final figures only; diagnostics remain off by default.
cfg.makeFigures = true;
cfg.makeDiagnosticFigures = false;

% Accuracy scatter condition: a meaningful interference stress point.
cfg.accuracySNR_dB = -5;
cfg.accuracySIR_dB = -5;

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
        error('V13_8:MissingFunction', ...
            'Required function not found: %s',requiredFunctions{k});
    end
    fprintf('%-40s %s\n',requiredFunctions{k},p);
end

if exist('nrTDLChannel','file') ~= 2 && exist('nrTDLChannel','class') ~= 8
    error('V13_8:Missing5GToolbox', ...
        ['nrTDLChannel was not found. This V13.8 script requires ', ...
         '5G Toolbox.']);
end

fprintf('\nEstimators intentionally NOT called:\n');
fprintf('  time_freq_recovery\n');
fprintf('  time_freq_recovery_V6\n');
fprintf('  prachPowerDelayProfileGen\n');

%% 2. Generate desired paired waveform, power reference and interferer

simConfig = struct( ...
    'nTx',cfg.nTx, ...
    'nRx',cfg.nRx, ...
    'puschBW',cfg.puschBW, ...
    'puschSCS',cfg.puschSCS, ...
    'prachFormat',cfg.prachFormat, ...
    'prachZeroConfigZone',cfg.prachZeroConfigZone, ...
    'prachRootSeqIdx',cfg.prachRootSeqIdx, ...
    'channelType','TDL', ...
    'snrdB',cfg.snrList_dB, ...
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

% Desired signal: paired-conjugate ZC, RootIndex 22, unchanged/unscaled.
cfgPaired = simConfig;
cfgPaired.simAlgo = 1;
cfgPaired.prachRootSeqIdx = cfg.prachRootSeqIdx;
[txPairedSig,prachConfig] = prachSigGen(cfgPaired);

% Power reference only: conventional ZC, same RootIndex 22.
% This waveform is not used as the desired signal in this experiment.
cfgTraditionalReference = simConfig;
cfgTraditionalReference.simAlgo = 0;
cfgTraditionalReference.prachRootSeqIdx = cfg.prachRootSeqIdx;
[txTraditionalReferenceSig,prachConfigTraditional] = ...
    prachSigGen(cfgTraditionalReference);

% Opposing-train interferer: conventional ZC, RootIndex 0.
cfgInterferer = simConfig;
cfgInterferer.simAlgo = 0;
cfgInterferer.prachRootSeqIdx = cfg.interfererRootIndex;
[txInterfererSig,prachConfigInterferer] = prachSigGen(cfgInterferer);

fsRaw = prachConfig.bwpInfo.SamplingRate;
NIFFT = prachConfig.bwpInfo.NFFT;
NFFTraw = prachConfig.nFFT;
NCPraw = round(prachConfig.NCP_RA);
Nrep = prachConfig.Nrep;
K = prachConfig.K;
fsEffective = fsRaw/K;

% All three waveforms must share the same F0 time/frequency configuration.
if prachConfigTraditional.NCP_RA~=prachConfig.NCP_RA || ...
        prachConfigTraditional.nFFT~=prachConfig.nFFT || ...
        prachConfigTraditional.Nrep~=prachConfig.Nrep || ...
        prachConfigTraditional.K~=prachConfig.K || ...
        prachConfigInterferer.NCP_RA~=prachConfig.NCP_RA || ...
        prachConfigInterferer.nFFT~=prachConfig.nFFT || ...
        prachConfigInterferer.Nrep~=prachConfig.Nrep || ...
        prachConfigInterferer.K~=prachConfig.K
    error('V13_8:PRACHConfigurationMismatch', ...
        'Target, traditional reference and interferer PRACH configurations differ.');
end

if Nrep ~= 1
    error('V13_8:UnexpectedNrep', ...
        'This validation expects PRACH F0 with Nrep=1.');
end
if abs(NFFTraw/NIFFT-K)>1e-9
    error('V13_8:SamplingRatioMismatch', ...
        'Raw/effective FFT ratio does not match prachConfig.K.');
end

rawBurstLength = NCPraw+NFFTraw*Nrep;
if size(txPairedSig,1)<rawBurstLength || ...
        size(txTraditionalReferenceSig,1)<rawBurstLength || ...
        size(txInterfererSig,1)<rawBurstLength
    error('V13_8:ShortTX', ...
        'One or more generated waveforms are shorter than the active PRACH burst.');
end

txBurstPaired = txPairedSig(1:rawBurstLength,1);
txBurstTraditionalReference = txTraditionalReferenceSig(1:rawBurstLength,1);
txBurstInterferer = txInterfererSig(1:rawBurstLength,1);

pairedEnergy = sum(abs(txBurstPaired).^2);
traditionalReferenceEnergy = sum(abs(txBurstTraditionalReference).^2);
interfererEnergy = sum(abs(txBurstInterferer).^2);
pairedEnergyGain_dB = 10*log10(pairedEnergy/(traditionalReferenceEnergy+eps));
interfererEnergyOffset_dB = 10*log10(interfererEnergy/(traditionalReferenceEnergy+eps));

rawInterferenceDelay = round( ...
    cfg.interferenceDelayEffSamples*fsRaw/fsEffective);

% Build the paired-ZC local reference through the SAME PRACH front end.
txRefData = prachFrontendProcessing( ...
    repmat(txPairedSig,1,cfg.nRx),prachConfig,simConfig);
txRef = reshape(txRefData(1,:,1),[],1);
if numel(txRef)~=NIFFT || any(~isfinite(txRef))
    error('V13_8:InvalidFineReference', ...
        'The front-end-processed paired-ZC reference is invalid.');
end

OS = NIFFT/prachConfig.L_RA;
cfoBinSpacing_Hz = prachConfig.scs_RA*1e3/OS;
maxAbsCFO_Hz = max(abs(cfg.CFOList_Hz));
maxCFOBins = maxAbsCFO_Hz/cfoBinSpacing_Hz;

%% 3. Configure independent static TDL-D target/interferer channels

tdlChan = makeTDLChannelLocal(cfg,fsRaw,cfg.randomSeed);
tdlInfo = info(tdlChan);

if ~isfield(tdlInfo,'ChannelFilterDelay') || ...
        ~isfield(tdlInfo,'MaximumChannelDelay') || ...
        ~isfield(tdlInfo,'PathDelays')
    error('V13_8:UnexpectedTDLInfo', ...
        'nrTDLChannel info structure lacks required delay fields.');
end

channelFilterDelay_raw = round(tdlInfo.ChannelFilterDelay);
tdlMaxExcessDelay_raw = max(tdlInfo.PathDelays)*fsRaw;
tdlMaxExcessDelay_eff = max(tdlInfo.PathDelays)*fsEffective;

% Preserve V13.7's target acquisition search bounds.
% The coarse estimator is not widened to search for the later interferer.
cfg.searchGuardBins = ceil(tdlMaxExcessDelay_eff)+cfg.searchGuardBinsExtra;
lagSearchMin = floor(min(cfg.delayList_effSamples)- ...
    maxCFOBins-cfg.searchGuardBins);
lagSearchMax = ceil(max(cfg.delayList_effSamples)+ ...
    maxCFOBins+cfg.searchGuardBins);

tailLength_raw = max(64,ceil(tdlInfo.MaximumChannelDelay)+32);

fprintf('\n===== V13.8 waveform / power setup =====\n');
fprintf('Target sequence               : paired-conjugate ZC, root %d\n', ...
    cfg.prachRootSeqIdx);
fprintf('Power reference               : traditional ZC, root %d\n', ...
    cfg.prachRootSeqIdx);
fprintf('Interferer sequence           : traditional ZC, root %d\n', ...
    cfg.interfererRootIndex);
fprintf('Paired TX energy gain         : %.4f dB (NOT normalized)\n', ...
    pairedEnergyGain_dB);
fprintf('Interferer TX energy offset   : %.4f dB vs traditional target\n', ...
    interfererEnergyOffset_dB);
fprintf('SNR list (traditional ref.)   : %s dB\n',mat2str(cfg.snrList_dB));
fprintf('SIR list (traditional ref.)   : %s dB\n',mat2str(cfg.sirList_dB));
fprintf('Relative interferer delay    : %.3f effective samples (%d raw)\n', ...
    cfg.interferenceDelayEffSamples,rawInterferenceDelay);
fprintf('Raw / effective Fs            : %.6f / %.6f MHz\n', ...
    fsRaw/1e6,fsEffective/1e6);
fprintf('CFO bin spacing               : %.9f Hz/bin\n',cfoBinSpacing_Hz);
fprintf('Target CFO range              : [%+.3f, %+.3f] Hz\n', ...
    min(cfg.CFOList_Hz),max(cfg.CFOList_Hz));

fprintf('\n===== Opposing-train geometry =====\n');
fprintf('Target speed                  : %+8.3f km/h\n',cfg.MobileSpeed_kmh);
fprintf('Interferer speed              : %+8.3f km/h\n',cfg.InterfererSpeed_kmh);
fprintf('Target / interferer tracks    : %.2f / %.2f m\n', ...
    cfg.TrackDistance_m,cfg.InterfererTrackDistance_m);
fprintf('Track separation              : %.2f m\n',cfg.TrackSeparation_m);
fprintf('Target fd at t=0              : %+.6f Hz\n',fdTarget0);
fprintf('Interferer fd at t=0          : %+.6f Hz\n',fdInterferer0);
fprintf('Relative fd at t=0            : %+.6f Hz\n', ...
    cfg.relativeGeometryDoppler_Hz);
fprintf('Target range / delay          : %.3f m / %.6f us\n', ...
    rangeTarget0_m,delayTarget0_s*1e6);
fprintf('Interferer range / delay      : %.3f m / %.6f us\n', ...
    rangeInterferer0_m,delayInterferer0_s*1e6);
fprintf('Geometry delay difference     : %.6f us (logged, not separately imposed)\n', ...
    cfg.geometryPropagationDelayDifference_s*1e6);

fprintf('\n===== TDL-D channel configuration =====\n');
fprintf('Profile                       : %s\n',cfg.tdlDelayProfile);
fprintf('RMS delay spread              : %.3f ns\n',cfg.tdlDelaySpread_s*1e9);
fprintf('Maximum Doppler in TDL        : %.3f Hz (static taps)\n', ...
    cfg.tdlMaximumDopplerShift_Hz);
fprintf('Channel filter delay          : %d raw samples\n', ...
    channelFilterDelay_raw);
fprintf('Maximum TDL excess delay      : %.3f raw / %.4f effective samples\n', ...
    tdlMaxExcessDelay_raw,tdlMaxExcessDelay_eff);
fprintf('Signed-lag target search      : [%d, %d]\n', ...
    lagSearchMin,lagSearchMax);
fprintf('Run mode / channel epochs     : %s / %d\n', ...
    cfg.runMode,cfg.nChannelRealizations);

maxTotalDelayRaw = ...
    (max(cfg.delayList_effSamples)+cfg.interferenceDelayEffSamples)*K + ...
    tdlMaxExcessDelay_raw;
if maxTotalDelayRaw >= NCPraw
    warning('V13_8:InterfererDelayExceedsCP', ...
        ['Target plus relative interferer delay and TDL excess delay ', ...
         '(%.2f raw samples) is not smaller than the PRACH CP (%d raw samples).'], ...
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

% Keep the V13.7 MP-ACR search settings unchanged.
% Explicitly select the strongest-path timing output for this comparison.
% The opposing interferer is intentionally NOT included in the fine
% estimator's search bounds; this measures the receiver's robustness
% without giving it an oracle interferer-specific initialization.
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
% Explicitly make the estimator's primary timing output the strongest path.
% The result table and 3D delay plot also read tauStrongestPathSamples.
fineOpts.timingPath = 'strongest';

fprintf('\n===== Fine estimator configuration =====\n');
fprintf('Fine estimator                 : HSR_MP_ACR\n');
fprintf('Residual CFO search radius     : +/- %.1f Hz\n',fineOpts.freqRadiusHz);
fprintf('Residual delay search radius   : +/- %.2f samples\n', ...
    fineOpts.delayRadiusSamples);
fprintf('Path search span               : +/- %.2f samples\n', ...
    fineOpts.pathSearchSpanSamples);
fprintf('Maximum fitted paths           : %d\n',fineOpts.maxPaths);

%% 6. Preallocate Monte-Carlo results

nCFO = numel(cfg.CFOList_Hz);
nDelay = numel(cfg.delayList_effSamples);
nSNR = numel(cfg.snrList_dB);
nSIR = numel(cfg.sirList_dB);
nRep = cfg.nChannelRealizations;

nTrials = nRep*nCFO*nDelay*nSNR*nSIR;

TrialID = (1:nTrials).';
ChannelRealization = nan(nTrials,1);
TrueCFO_Hz = nan(nTrials,1);
TrueDelay_effSamples = nan(nTrials,1); % nominal propagation delay
TrueDelay_rawSamples = nan(nTrials,1);
TrueFirstPathDelay_effSamples = nan(nTrials,1);
TrueStrongestPathDelay_effSamples = nan(nTrials,1);
TargetStrongestPathIndex = nan(nTrials,1);
TargetStrongestPathPower = nan(nTrials,1);
InterfererCFOReference_Hz = nan(nTrials,1);
InterfererDelay_effSamples = nan(nTrials,1);
RelativeDopplerMean_Hz = nan(nTrials,1);
SNR_dB = nan(nTrials,1);
SIR_dB = nan(nTrials,1);

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
DelayErrorCoarseStrongest_effSamples = nan(nTrials,1);
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

CoarseTargetLock = false(nTrials,1);
CoarseInterfererLock = false(nTrials,1);
FineEstimatorValid = false(nTrials,1);
FineTargetLock = false(nTrials,1);
FineInterfererLock = false(nTrials,1);
FineErrorID = strings(nTrials,1);
FineErrorMessage = strings(nTrials,1);

NoiseVariance = nan(nTrials,1);
TargetTraditionalReferencePower = nan(nTrials,1);
TargetPairedPower = nan(nTrials,1);
InterfererPowerBeforeScale = nan(nTrials,1);
InterferenceScale = nan(nTrials,1);
ActualSNR_dB = nan(nTrials,1);
ActualSIR_dB = nan(nTrials,1);
ActualPairedSNR_dB = nan(nTrials,1);
ActualPairedSIR_dB = nan(nTrials,1);
CompositePower = nan(nTrials,1);

%% 7. Main opposing-train TDL-D Monte-Carlo loop

fprintf('\n===== V13.8 opposing-train Monte-Carlo run =====\n');
fprintf('Expected total trials: %d\n',nTrials);

trial = 0;

for iRep = 1:nRep

    % Independent static TDL-D realizations for desired and interfering trains.
    % The paired target and its traditional power-reference waveform use
    % separate channel objects with the SAME seed, ensuring matched taps.
    targetSeed = cfg.randomSeed+10000+iRep;
    interfererSeed = cfg.randomSeed+20000+iRep;

    targetChanPaired = makeTDLChannelLocal(cfg,fsRaw,targetSeed);
    targetChanTraditional = makeTDLChannelLocal(cfg,fsRaw,targetSeed);
    interfererChan = makeTDLChannelLocal(cfg,fsRaw,interfererSeed);

    [rxTargetPairedStatic,targetPathGains] = generateTDLRealizationLocal( ...
        targetChanPaired,txBurstPaired,tailLength_raw,channelFilterDelay_raw);
    rxTargetTraditionalStatic = generateTDLRealizationLocal( ...
        targetChanTraditional,txBurstTraditionalReference, ...
        tailLength_raw,channelFilterDelay_raw);
    rxInterfererStatic = generateTDLRealizationLocal( ...
        interfererChan,txBurstInterferer,tailLength_raw,channelFilterDelay_raw);

    % Determine the strongest target TDL tap from this exact realization.
    % The channel-filter delay is removed inside generateTDLRealizationLocal,
    % so the path-delay reference below contains only the TDL tap delay.
    [targetStrongestPathIndex,targetPathPower] = ...
        getStrongestTDLPathLocal(targetPathGains,numel(tdlInfo.PathDelays));
    targetFirstPathDelay_effSamples = ...
        min(tdlInfo.PathDelays)*fsEffective;
    targetStrongestPathDelay_effSamples = ...
        tdlInfo.PathDelays(targetStrongestPathIndex)*fsEffective;

    for iCFO = 1:nCFO

        currentCFO = cfg.CFOList_Hz(iCFO);

        for iDelay = 1:nDelay

            currentDelayEff = cfg.delayList_effSamples(iDelay);
            currentDelayRaw = currentDelayEff*K;

            % Ground-truth path delays include the imposed propagation
            % delay plus the actual first/strongest tap of this TDL-D draw.
            currentTrueFirstPathDelayEff = ...
                currentDelayEff+targetFirstPathDelay_effSamples;
            currentTrueStrongestPathDelayEff = ...
                currentDelayEff+targetStrongestPathDelay_effSamples;

            currentInterfererDelayEff = ...
                currentDelayEff+cfg.interferenceDelayEffSamples;
            currentInterfererDelayRaw = currentDelayRaw+rawInterferenceDelay;

            % Apply the same absolute target delay to paired and traditional
            % target reference. Interferer arrives 8 effective samples later.
            rxTargetPairedDelayed = fractionalDelayRawLocal( ...
                rxTargetPairedStatic,currentDelayRaw, ...
                cfg.fractionalDelayGuardSamples);
            rxTargetTraditionalDelayed = fractionalDelayRawLocal( ...
                rxTargetTraditionalStatic,currentDelayRaw, ...
                cfg.fractionalDelayGuardSamples);
            rxInterfererDelayed = fractionalDelayRawLocal( ...
                rxInterfererStatic,currentInterfererDelayRaw, ...
                cfg.fractionalDelayGuardSamples);

            nComposite = max([size(rxTargetPairedDelayed,1), ...
                size(rxTargetTraditionalDelayed,1), ...
                size(rxInterfererDelayed,1)]);

            rxTargetPairedBase = padRowsLocal(rxTargetPairedDelayed,nComposite);
            rxTargetTraditionalBase = padRowsLocal( ...
                rxTargetTraditionalDelayed,nComposite);
            rxInterfererBase = padRowsLocal(rxInterfererDelayed,nComposite);

            tRaw = (0:nComposite-1).'/fsRaw;

            % Apply target CFO exactly as in V13.7. The interferer is
            % frequency-shifted relative to the target receiver reference
            % by its own opposing-train geometry Doppler.
            [fdTargetVector,~,~] = railwayDopplerModel( ...
                cfg.TrainPosition_m,cfg.MobileSpeed_kmh/3.6, ...
                cfg.TrackDistance_m, ...
                cfg.BaseHeight_m-cfg.TrainAntennaHeight_m, ...
                cfg.CarrierFrequency_Hz,tRaw);
            [fdInterfererVector,~,~] = railwayDopplerModel( ...
                cfg.InterfererPosition_m,cfg.InterfererSpeed_kmh/3.6, ...
                cfg.InterfererTrackDistance_m, ...
                cfg.BaseHeight_m-cfg.TrainAntennaHeight_m, ...
                cfg.CarrierFrequency_Hz,tRaw);

            relativeDoppler = fdInterfererVector(:)-fdTargetVector(:);
            interfererInitialPhase = randomPhaseLocal( ...
                cfg.randomSeed+50000+100*iRep+iDelay);
            relativePhase = 2*pi*cumtrapz(tRaw,relativeDoppler) + ...
                interfererInitialPhase;

            targetCFOPhase = exp(1j*2*pi*currentCFO*tRaw);
            interfererCFOPhase = exp(1j*( ...
                2*pi*currentCFO*tRaw+relativePhase));

            rxTargetPairedClean = bsxfun(@times, ...
                rxTargetPairedBase,targetCFOPhase);
            rxInterfererClean = bsxfun(@times, ...
                rxInterfererBase,interfererCFOPhase);

            % Active support and power reference follow the previous A2
            % convention: traditional Root-22 target received power.
            activeIdx = signalSupportWindowLocal(rxTargetPairedBase);
            if isempty(activeIdx)
                error('V13_8:EmptyTargetSupport', ...
                    'No active target support at realization %d, delay %.6f.', ...
                    iRep,currentDelayEff);
            end

            pTargetTraditionalReference = mean( ...
                abs(rxTargetTraditionalBase(activeIdx,:)).^2,'all');
            pTargetPaired = mean(abs(rxTargetPairedBase(activeIdx,:)).^2,'all');
            pInterfererReference = mean(abs(rxInterfererClean(activeIdx,:)).^2,'all');

            if ~isfinite(pTargetTraditionalReference) || ...
                    pTargetTraditionalReference<=0 || ...
                    ~isfinite(pTargetPaired) || pTargetPaired<=0 || ...
                    ~isfinite(pInterfererReference) || pInterfererReference<=0
                error('V13_8:InvalidLinkPower', ...
                    'Target/interferer active-window power is invalid.');
            end

            relativeDopplerMean = mean(relativeDoppler(activeIdx));
            interfererCFOReference = currentCFO+relativeDopplerMean;

            % Keep V13.7's common per-Rx gain derivation from the CLEAN
            % desired paired signal only. Interference does not redefine AGC.
            tdClean = prachFrontendProcessing( ...
                rxTargetPairedClean,prachConfig,simConfig);
            [rxGain,refGain] = deriveSharedGainsLocal( ...
                tdClean,prachConfig,cfg.nRx);

            optsNow = estimatorOpts;
            optsNow.RxGain = rxGain;
            optsNow.RefGain = refGain;

            for iSNR = 1:nSNR

                currentSNR = cfg.snrList_dB(iSNR);

                if isinf(currentSNR)
                    noiseVariance = 0;
                else
                    noiseVariance = pTargetTraditionalReference / ...
                        10^(currentSNR/10);
                end

                noiseSeed = cfg.randomSeed + 1000000 + ...
                    iRep*10000+iCFO*1000+iDelay*100+iSNR;
                noise = makeComplexNoiseLocal( ...
                    [nComposite,cfg.nRx],noiseVariance,noiseSeed);

                pNoiseMeasured = mean(abs(noise(activeIdx,:)).^2,'all');
                if pNoiseMeasured<=0
                    actualSNR = Inf;
                else
                    actualSNR = 10*log10( ...
                        pTargetTraditionalReference/(pNoiseMeasured+eps));
                end
                if pNoiseMeasured<=0
                    actualPairedSNR = Inf;
                else
                    actualPairedSNR = 10*log10( ...
                        pTargetPaired/(pNoiseMeasured+eps));
                end

                for iSIR = 1:nSIR

                    currentSIR = cfg.sirList_dB(iSIR);

                    if isinf(currentSIR)
                        interferenceScale = 0;
                        rxInterfererScaled = zeros(size(rxInterfererClean), ...
                            'like',rxInterfererClean);
                        actualSIR = Inf;
                        actualPairedSIR = Inf;
                    else
                        interferenceScale = sqrt( ...
                            pTargetTraditionalReference / ...
                            (pInterfererReference*10^(currentSIR/10)));
                        rxInterfererScaled = ...
                            rxInterfererClean.*interferenceScale;

                        pInterfererScaled = mean( ...
                            abs(rxInterfererScaled(activeIdx,:)).^2,'all');
                        actualSIR = 10*log10( ...
                            pTargetTraditionalReference/(pInterfererScaled+eps));
                        actualPairedSIR = 10*log10( ...
                            pTargetPaired/(pInterfererScaled+eps));
                    end

                    rxComposite = rxTargetPairedClean + ...
                        rxInterfererScaled + noise;

                    compositePower = mean(abs(rxComposite(activeIdx,:)).^2,'all');

                    timeData = prachFrontendProcessing( ...
                        rxComposite,prachConfig,simConfig);

                    [delayRaw,cfoEstimate,diagInfo] = ...
                        pairedZCFullLagCoarseEstimator( ...
                        timeData,prachConfig,optsNow);

                    rxWin = reshape(timeData,cfg.nRx,[]).';
                    if size(rxWin,1)~=NIFFT || size(rxWin,2)~=cfg.nRx
                        error('V13_8:InvalidFineInputWindow', ...
                            'Unexpected effective-rate receive-window dimensions.');
                    end

                    % Exactly the V13.7 coarse-to-fine handoff.
                    rxAligned = advanceCircular(rxWin,delayRaw);

                    fineIsValid = false;
                    fineErrorIDNow = "";
                    fineErrorMessageNow = "";
                    fineOut = struct( ...
                        'freqFineHz',NaN, ...
                        'tauFirstPathSamples',NaN, ...
                        'tauStrongestPathSamples',NaN, ...
                        'deltaFreqHz',NaN, ...
                        'deltaTauSamples',NaN, ...
                        'numPaths',NaN, ...
                        'normalizedResidual',NaN);

                    try
                        fineCandidate = HSR_MP_ACR( ...
                            rxAligned,txRef,fsEffective, ...
                            delayRaw,cfoEstimate,fineOpts);

                        requiredFineFields = { ...
                            'freqFineHz','tauFirstPathSamples', ...
                            'tauStrongestPathSamples','deltaFreqHz', ...
                            'deltaTauSamples','numPaths','normalizedResidual'};
                        for iField = 1:numel(requiredFineFields)
                            if ~isfield(fineCandidate,requiredFineFields{iField})
                                error('V13_8:MissingFineOutput', ...
                                    'HSR_MP_ACR output is missing field %s.', ...
                                    requiredFineFields{iField});
                            end
                        end

                        if all(isfinite([fineCandidate.freqFineHz, ...
                                fineCandidate.tauFirstPathSamples, ...
                                fineCandidate.tauStrongestPathSamples, ...
                                fineCandidate.deltaFreqHz, ...
                                fineCandidate.deltaTauSamples, ...
                                fineCandidate.numPaths, ...
                                fineCandidate.normalizedResidual]))
                            fineOut = fineCandidate;
                            fineIsValid = true;
                        else
                            fineErrorIDNow = "V13_8:InvalidFineOutput";
                            fineErrorMessageNow = ...
                                "HSR_MP_ACR returned one or more non-finite fields.";
                        end
                    catch MEfine
                        fineErrorIDNow = string(MEfine.identifier);
                        fineErrorMessageNow = string(MEfine.message);
                    end

                    trial = trial+1;

                    TrialID(trial) = trial;
                    ChannelRealization(trial) = iRep;
                    TrueCFO_Hz(trial) = currentCFO;
                    TrueDelay_effSamples(trial) = currentDelayEff;
                    TrueDelay_rawSamples(trial) = currentDelayRaw;
                    TrueFirstPathDelay_effSamples(trial) = ...
                        currentTrueFirstPathDelayEff;
                    TrueStrongestPathDelay_effSamples(trial) = ...
                        currentTrueStrongestPathDelayEff;
                    TargetStrongestPathIndex(trial) = targetStrongestPathIndex;
                    TargetStrongestPathPower(trial) = ...
                        targetPathPower(targetStrongestPathIndex);
                    InterfererCFOReference_Hz(trial) = interfererCFOReference;
                    InterfererDelay_effSamples(trial) = currentInterfererDelayEff;
                    RelativeDopplerMean_Hz(trial) = relativeDopplerMean;
                    SNR_dB(trial) = currentSNR;
                    SIR_dB(trial) = currentSIR;

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
                    CFOErrorCoarse_Hz(trial) = cfoEstimate-currentCFO;
                    DelayErrorCoarseStrongest_effSamples(trial) = ...
                        delayRaw-currentTrueStrongestPathDelayEff;

                    if fineIsValid
                        CFOEstimateFine_Hz(trial) = fineOut.freqFineHz;
                        DelayEstimateFineFirst_effSamples(trial) = ...
                            fineOut.tauFirstPathSamples;
                        DelayEstimateFineStrongest_effSamples(trial) = ...
                            fineOut.tauStrongestPathSamples;

                        CFOErrorFine_Hz(trial) = fineOut.freqFineHz-currentCFO;
                        DelayErrorFineFirst_effSamples(trial) = ...
                            fineOut.tauFirstPathSamples-currentTrueFirstPathDelayEff;
                        DelayErrorFineStrongest_effSamples(trial) = ...
                            fineOut.tauStrongestPathSamples-currentTrueStrongestPathDelayEff;

                        FineDeltaFreq_Hz(trial) = fineOut.deltaFreqHz;
                        FineDeltaTauSamples(trial) = fineOut.deltaTauSamples;
                        FineNumPaths(trial) = fineOut.numPaths;
                        FineNormalizedResidual(trial) = fineOut.normalizedResidual;
                    end

                    expectedD1 = currentTrueStrongestPathDelayEff+ ...
                        currentCFO/cfoBinSpacing_Hz;
                    expectedD2 = currentTrueStrongestPathDelayEff- ...
                        currentCFO/cfoBinSpacing_Hz;
                    Branch1LagError_effSamples(trial) = ...
                        diagInfo.branch1LagSubBin-expectedD1;
                    Branch2LagError_effSamples(trial) = ...
                        diagInfo.branch2LagSubBin-expectedD2;

                    Peak1AtSearchBoundary(trial) = ...
                        diagInfo.branch1PeakAtSearchBoundary;
                    Peak2AtSearchBoundary(trial) = ...
                        diagInfo.branch2PeakAtSearchBoundary;

                    FineEstimatorValid(trial) = fineIsValid;
                    FineErrorID(trial) = fineErrorIDNow;
                    FineErrorMessage(trial) = fineErrorMessageNow;

                    CoarseTargetLock(trial) = ...
                        abs(cfoEstimate-currentCFO)<=fineOpts.freqRadiusHz && ...
                        abs(delayRaw-currentTrueStrongestPathDelayEff)<= ...
                            fineOpts.delayRadiusSamples;
                    CoarseInterfererLock(trial) = ...
                        abs(cfoEstimate-interfererCFOReference)<=fineOpts.freqRadiusHz && ...
                        abs(delayRaw-currentInterfererDelayEff)<= ...
                            fineOpts.delayRadiusSamples;

                    if fineIsValid
                        FineTargetLock(trial) = ...
                            abs(fineOut.freqFineHz-currentCFO)<=fineOpts.freqRadiusHz && ...
                            abs(fineOut.tauStrongestPathSamples- ...
                                currentTrueStrongestPathDelayEff)<= ...
                                fineOpts.delayRadiusSamples;
                        FineInterfererLock(trial) = ...
                            abs(fineOut.freqFineHz-interfererCFOReference)<= ...
                                fineOpts.freqRadiusHz && ...
                            abs(fineOut.tauFirstPathSamples-currentInterfererDelayEff)<= ...
                                fineOpts.delayRadiusSamples;
                    end

                    NoiseVariance(trial) = noiseVariance;
                    TargetTraditionalReferencePower(trial) = ...
                        pTargetTraditionalReference;
                    TargetPairedPower(trial) = pTargetPaired;
                    InterfererPowerBeforeScale(trial) = pInterfererReference;
                    InterferenceScale(trial) = interferenceScale;
                    ActualSNR_dB(trial) = actualSNR;
                    ActualSIR_dB(trial) = actualSIR;
                    ActualPairedSNR_dB(trial) = actualPairedSNR;
                    ActualPairedSIR_dB(trial) = actualPairedSIR;
                    CompositePower(trial) = compositePower;
                end
            end
        end
    end

    fprintf('Completed target/interferer TDL-D realization %d/%d (%d/%d trials).\n', ...
        iRep,nRep,trial,nTrials);
end

if trial~=nTrials
    error('V13_8:TrialCountMismatch', ...
        'Generated %d trials but preallocated %d.',trial,nTrials);
end

%% 8. Raw results table

resultsTable = table( ...
    TrialID,ChannelRealization,TrueCFO_Hz,TrueDelay_effSamples, ...
    TrueDelay_rawSamples,TrueFirstPathDelay_effSamples, ...
    TrueStrongestPathDelay_effSamples,TargetStrongestPathIndex, ...
    TargetStrongestPathPower,InterfererCFOReference_Hz, ...
    InterfererDelay_effSamples,RelativeDopplerMean_Hz,SNR_dB,SIR_dB, ...
    Branch1LagInteger,Branch1LagSubBin,Branch2LagInteger,Branch2LagSubBin, ...
    Branch1PeakValue,Branch2PeakValue,Branch2ToBranch1Peak_dB, ...
    Branch1PeakToMedian_dB,Branch2PeakToMedian_dB, ...
    CFOEstimateCoarse_Hz,DelayEstimateCoarse_effSamples, ...
    CFOEstimateFine_Hz,DelayEstimateFineFirst_effSamples, ...
    DelayEstimateFineStrongest_effSamples, ...
    CFOErrorCoarse_Hz,CFOErrorFine_Hz, ...
    DelayErrorCoarseStrongest_effSamples, ...
    DelayErrorFineFirst_effSamples,DelayErrorFineStrongest_effSamples, ...
    FineDeltaFreq_Hz,FineDeltaTauSamples,FineNumPaths, ...
    FineNormalizedResidual,Branch1LagError_effSamples, ...
    Branch2LagError_effSamples,Peak1AtSearchBoundary, ...
    Peak2AtSearchBoundary,CoarseTargetLock,CoarseInterfererLock, ...
    FineEstimatorValid,FineTargetLock,FineInterfererLock, ...
    FineErrorID,FineErrorMessage, ...
    NoiseVariance,TargetTraditionalReferencePower,TargetPairedPower, ...
    InterfererPowerBeforeScale,InterferenceScale,ActualSNR_dB, ...
    ActualSIR_dB,ActualPairedSNR_dB,ActualPairedSIR_dB,CompositePower);

rawFile = 'V13_8_TDLD_OpposingTrain_CoarseFine_Raw.csv';
writetable(resultsTable,rawFile);

%% 9. Per-SNR / per-SIR coarse-fine metrics

metricSNR = strings(0,1);
metricSIR = strings(0,1);
metricCount = zeros(0,1);

CFO_RMSE_Coarse_Hz = [];
CFO_RMSE_Fine_Hz = [];
CFO_MAE_Coarse_Hz = [];
CFO_MAE_Fine_Hz = [];
CFO_P95AbsError_Fine_Hz = [];
Delay_RMSE_CoarseStrongest_effSamples = [];
Delay_RMSE_FineFirst_effSamples = [];
Delay_RMSE_FineStrongest_effSamples = [];
Delay_MAE_FineFirst_effSamples = [];
Delay_P95AbsError_FineFirst_effSamples = [];
FineResidual_Mean = [];
FineNumPaths_Mean = [];
FineEstimatorReturnRate = [];
CoarseTargetLockRate = [];
CoarseInterfererLockRate = [];
FineTargetLockRate = [];
FineInterfererLockRate = [];
BoundaryRate = [];
ActualSNR_Mean_dB = [];
ActualSIR_Mean_dB = [];
ActualPairedSNR_Mean_dB = [];
ActualPairedSIR_Mean_dB = [];

for iSNR = 1:nSNR

    if isinf(cfg.snrList_dB(iSNR))
        snrMask = isinf(resultsTable.SNR_dB);
        snrLabel = "NoNoise";
    else
        snrMask = resultsTable.SNR_dB==cfg.snrList_dB(iSNR);
        snrLabel = string(cfg.snrList_dB(iSNR))+"dB";
    end

    for iSIR = 1:nSIR

        if isinf(cfg.sirList_dB(iSIR))
            sirMask = isinf(resultsTable.SIR_dB);
            sirLabel = "NoInterference";
        else
            sirMask = resultsTable.SIR_dB==cfg.sirList_dB(iSIR);
            sirLabel = string(cfg.sirList_dB(iSIR))+"dB";
        end

        mask = snrMask & sirMask;

        eCFOCoarse = resultsTable.CFOErrorCoarse_Hz(mask);
        eCFOFine = resultsTable.CFOErrorFine_Hz(mask);
        eDelayCoarseStrongest = ...
            resultsTable.DelayErrorCoarseStrongest_effSamples(mask);
        eDelayFineFirst = resultsTable.DelayErrorFineFirst_effSamples(mask);
        eDelayFineStrongest = resultsTable.DelayErrorFineStrongest_effSamples(mask);

        metricSNR(end+1,1) = snrLabel; %#ok<AGROW>
        metricSIR(end+1,1) = sirLabel; %#ok<AGROW>
        metricCount(end+1,1) = sum(mask); %#ok<AGROW>

        CFO_RMSE_Coarse_Hz(end+1,1) = ...
            sqrt(mean(eCFOCoarse.^2,'omitnan')); %#ok<AGROW>
        CFO_RMSE_Fine_Hz(end+1,1) = ...
            sqrt(mean(eCFOFine.^2,'omitnan')); %#ok<AGROW>
        CFO_MAE_Coarse_Hz(end+1,1) = ...
            mean(abs(eCFOCoarse),'omitnan'); %#ok<AGROW>
        CFO_MAE_Fine_Hz(end+1,1) = ...
            mean(abs(eCFOFine),'omitnan'); %#ok<AGROW>
        CFO_P95AbsError_Fine_Hz(end+1,1) = ...
            localPercentile(abs(eCFOFine),95); %#ok<AGROW>

        Delay_RMSE_CoarseStrongest_effSamples(end+1,1) = ...
            sqrt(mean(eDelayCoarseStrongest.^2,'omitnan')); %#ok<AGROW>
        Delay_RMSE_FineFirst_effSamples(end+1,1) = ...
            sqrt(mean(eDelayFineFirst.^2,'omitnan')); %#ok<AGROW>
        Delay_RMSE_FineStrongest_effSamples(end+1,1) = ...
            sqrt(mean(eDelayFineStrongest.^2,'omitnan')); %#ok<AGROW>
        Delay_MAE_FineFirst_effSamples(end+1,1) = ...
            mean(abs(eDelayFineFirst),'omitnan'); %#ok<AGROW>
        Delay_P95AbsError_FineFirst_effSamples(end+1,1) = ...
            localPercentile(abs(eDelayFineFirst),95); %#ok<AGROW>

        FineResidual_Mean(end+1,1) = ...
            mean(resultsTable.FineNormalizedResidual(mask),'omitnan'); %#ok<AGROW>
        FineNumPaths_Mean(end+1,1) = ...
            mean(resultsTable.FineNumPaths(mask),'omitnan'); %#ok<AGROW>
        FineEstimatorReturnRate(end+1,1) = ...
            mean(resultsTable.FineEstimatorValid(mask)); %#ok<AGROW>
        CoarseTargetLockRate(end+1,1) = ...
            mean(resultsTable.CoarseTargetLock(mask)); %#ok<AGROW>
        CoarseInterfererLockRate(end+1,1) = ...
            mean(resultsTable.CoarseInterfererLock(mask)); %#ok<AGROW>
        FineTargetLockRate(end+1,1) = ...
            mean(resultsTable.FineTargetLock(mask)); %#ok<AGROW>
        FineInterfererLockRate(end+1,1) = ...
            mean(resultsTable.FineInterfererLock(mask)); %#ok<AGROW>
        BoundaryRate(end+1,1) = mean( ...
            resultsTable.Peak1AtSearchBoundary(mask) | ...
            resultsTable.Peak2AtSearchBoundary(mask)); %#ok<AGROW>

        ActualSNR_Mean_dB(end+1,1) = ...
            mean(resultsTable.ActualSNR_dB(mask),'omitnan'); %#ok<AGROW>
        ActualSIR_Mean_dB(end+1,1) = ...
            mean(resultsTable.ActualSIR_dB(mask),'omitnan'); %#ok<AGROW>
        ActualPairedSNR_Mean_dB(end+1,1) = ...
            mean(resultsTable.ActualPairedSNR_dB(mask),'omitnan'); %#ok<AGROW>
        ActualPairedSIR_Mean_dB(end+1,1) = ...
            mean(resultsTable.ActualPairedSIR_dB(mask),'omitnan'); %#ok<AGROW>
    end
end

metricsTable = table(metricSNR,metricSIR,metricCount, ...
    CFO_RMSE_Coarse_Hz,CFO_RMSE_Fine_Hz, ...
    CFO_MAE_Coarse_Hz,CFO_MAE_Fine_Hz,CFO_P95AbsError_Fine_Hz, ...
    Delay_RMSE_CoarseStrongest_effSamples, ...
    Delay_RMSE_FineFirst_effSamples, ...
    Delay_RMSE_FineStrongest_effSamples, ...
    Delay_MAE_FineFirst_effSamples,Delay_P95AbsError_FineFirst_effSamples, ...
    FineResidual_Mean,FineNumPaths_Mean,FineEstimatorReturnRate, ...
    CoarseTargetLockRate,CoarseInterfererLockRate, ...
    FineTargetLockRate,FineInterfererLockRate,BoundaryRate, ...
    ActualSNR_Mean_dB,ActualSIR_Mean_dB, ...
    ActualPairedSNR_Mean_dB,ActualPairedSIR_Mean_dB);

metricsFile = 'V13_8_TDLD_OpposingTrain_CoarseFine_Metrics.csv';
writetable(metricsTable,metricsFile);

%% 10. Save complete experiment data

matFile = 'V13_8_TDLD_OpposingTrain_CoarseFine.mat';
save(matFile,'cfg','simConfig','prachConfig','prachConfigTraditional', ...
    'prachConfigInterferer','tdlInfo','estimatorOpts','fineOpts','txRef', ...
    'resultsTable','metricsTable','pairedEnergyGain_dB', ...
    'interfererEnergyOffset_dB','fdTarget0','fdInterferer0', ...
    'rangeTarget0_m','rangeInterferer0_m','delayTarget0_s', ...
    'delayInterferer0_s','fsRaw','fsEffective','NIFFT','K','OS', ...
    'cfoBinSpacing_Hz','maxAbsCFO_Hz','maxCFOBins', ...
    'channelFilterDelay_raw','tdlMaxExcessDelay_raw', ...
    'tdlMaxExcessDelay_eff','lagSearchMin','lagSearchMax', ...
    'rawInterferenceDelay','tailLength_raw','-v7.3');

%% 11. Separate 3D coarse-versus-fine RMSE figures

if cfg.makeFigures

    % Metrics are pooled across the configured CFO and nominal-delay test
    % points for each SNR/SIR pair. The two estimators use identical trials.
    % Rows = SNR, columns = SIR.
    mapCFOCoarse = reshape( ...
        metricsTable.CFO_RMSE_Coarse_Hz,nSIR,nSNR).';
    mapCFOFine = reshape( ...
        metricsTable.CFO_RMSE_Fine_Hz,nSIR,nSNR).';

    mapDelayCoarse = reshape( ...
        metricsTable.Delay_RMSE_CoarseStrongest_effSamples,nSIR,nSNR).';
    mapDelayFineStrongest = reshape( ...
        metricsTable.Delay_RMSE_FineStrongest_effSamples,nSIR,nSNR).';

    % Use physical dB values on the axes, not categorical index numbers.
    [sirGrid,snrGrid] = meshgrid(cfg.sirList_dB,cfg.snrList_dB);

    coarseColor = [0.0000 0.3500 0.7200];
    fineColor   = [0.9000 0.3000 0.0800];

    % ---------------- Figure 1: CFO RMSE ----------------
    figCFO = figure('Color','w','Position',[80 80 1120 820]);
    ax1 = axes(figCFO);
    hold(ax1,'on');

    hCFOCoarse = surf(ax1,sirGrid,snrGrid,mapCFOCoarse, ...
        'FaceColor','none', ...
        'EdgeColor',coarseColor, ...
        'LineStyle','-', ...
        'LineWidth',1.8, ...
        'DisplayName','Coarse: full-lag + 3-point parabolic');

    hCFOFine = surf(ax1,sirGrid,snrGrid,mapCFOFine, ...
        'FaceColor','none', ...
        'EdgeColor',fineColor, ...
        'LineStyle','--', ...
        'LineWidth',1.8, ...
        'DisplayName','Fine: HSR\_MP\_ACR');

    scatter3(ax1,sirGrid(:),snrGrid(:),mapCFOCoarse(:),62, ...
        coarseColor,'filled', ...
        'Marker','o', ...
        'MarkerFaceColor',coarseColor, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.55, ...
        'HandleVisibility','off');

    scatter3(ax1,sirGrid(:),snrGrid(:),mapCFOFine(:),62, ...
        fineColor,'filled', ...
        'Marker','d', ...
        'MarkerFaceColor',fineColor, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.55, ...
        'HandleVisibility','off');

    hold(ax1,'off');
    grid(ax1,'on');
    box(ax1,'on');
    xlabel(ax1,'SIR (dB)  [ +10 \rightarrow -10: interference increases ]');
    ylabel(ax1,'SNR (dB)');
    zlabel(ax1,'CFO RMSE (Hz)');
    title(ax1,{ ...
        'CFO RMSE: Coarse Three-Point Parabolic vs MP-ACR Fine', ...
        'TDL-D opposing-train interference; both methods use the same trials'});
    legend(ax1,[hCFOCoarse,hCFOFine], ...
        {'Coarse: full-lag + 3-point parabolic','Fine: HSR\_MP\_ACR'}, ...
        'Location','northeast');
    xlim(ax1,[-10 10]);
    ylim(ax1,[-15 5]);
    xticks(ax1,-10:5:10);
    yticks(ax1,-15:5:5);
    set(ax1,'XDir','reverse','FontSize',11,'LineWidth',1.0);
    view(ax1,[-38 27]);
    set(ax1,'Projection','perspective');
    savefig(figCFO,'V13_8_TDLD_OpposingTrain_CFO_CoarseFine_3D.fig');
    exportgraphics(figCFO, ...
        'V13_8_TDLD_OpposingTrain_CFO_CoarseFine_3D.png', ...
        'Resolution',300);

    % ---------------- Figure 2: strongest-path delay RMSE ----------------
    figDelay = figure('Color','w','Position',[80 80 1120 820]);
    ax2 = axes(figDelay);
    hold(ax2,'on');

    hDelayCoarse = surf(ax2,sirGrid,snrGrid,mapDelayCoarse, ...
        'FaceColor','none', ...
        'EdgeColor',coarseColor, ...
        'LineStyle','-', ...
        'LineWidth',1.8, ...
        'DisplayName','Coarse: full-lag + 3-point parabolic');

    hDelayFine = surf(ax2,sirGrid,snrGrid,mapDelayFineStrongest, ...
        'FaceColor','none', ...
        'EdgeColor',fineColor, ...
        'LineStyle','--', ...
        'LineWidth',1.8, ...
        'DisplayName','Fine: HSR\_MP\_ACR strongest path');

    scatter3(ax2,sirGrid(:),snrGrid(:),mapDelayCoarse(:),62, ...
        coarseColor,'filled', ...
        'Marker','o', ...
        'MarkerFaceColor',coarseColor, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.55, ...
        'HandleVisibility','off');

    scatter3(ax2,sirGrid(:),snrGrid(:),mapDelayFineStrongest(:),62, ...
        fineColor,'filled', ...
        'Marker','d', ...
        'MarkerFaceColor',fineColor, ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.55, ...
        'HandleVisibility','off');

    hold(ax2,'off');
    grid(ax2,'on');
    box(ax2,'on');
    xlabel(ax2,'SIR (dB)  [ +10 \rightarrow -10: interference increases ]');
    ylabel(ax2,'SNR (dB)');
    zlabel(ax2,'Strongest-path delay RMSE (effective samples)');
    title(ax2,{ ...
        'Strongest-Path Delay RMSE: Coarse Three-Point Parabolic vs MP-ACR Fine', ...
        'Ground truth uses the strongest tap in each actual TDL-D channel realization'});
    legend(ax2,[hDelayCoarse,hDelayFine], ...
        {'Coarse: full-lag + 3-point parabolic', ...
         'Fine: HSR\_MP\_ACR strongest path'}, ...
        'Location','northeast');
    xlim(ax2,[-10 10]);
    ylim(ax2,[-15 5]);
    xticks(ax2,-10:5:10);
    yticks(ax2,-15:5:5);
    set(ax2,'XDir','reverse','FontSize',11,'LineWidth',1.0);
    view(ax2,[-38 27]);
    set(ax2,'Projection','perspective');
    savefig(figDelay, ...
        'V13_8_TDLD_OpposingTrain_Delay_Strongest_CoarseFine_3D.fig');
    exportgraphics(figDelay, ...
        'V13_8_TDLD_OpposingTrain_Delay_Strongest_CoarseFine_3D.png', ...
        'Resolution',300);
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
    scatter(resultsTable.TrueStrongestPathDelay_effSamples, ...
        resultsTable.DelayErrorCoarseStrongest_effSamples,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True strongest-path delay (effective samples)');
    ylabel('Coarse strongest-path delay error (samples)');
    title('Coarse strongest-path delay error versus true delay');

    nexttile(td);
    scatter(resultsTable.TrueCFO_Hz, ...
        resultsTable.CFOErrorCoarse_Hz,14, ...
        resultsTable.SNR_dB,'filled');
    colorbar;
    grid on;
    xlabel('True CFO (Hz)');
    ylabel('CFO error (Hz)');
    title('CFO error versus true CFO');

    title(td,'V13.8 TDL-D coarse/fine diagnostic plots');
    savefig(figDiag,'V13_8_TDLD_OpposingTrain_CoarseFine_Diagnostics.fig');
    exportgraphics(figDiag, ...
        'V13_8_TDLD_OpposingTrain_CoarseFine_Diagnostics.png','Resolution',250);
end

%% 13. Summary

fprintf('\n============================================================\n');
fprintf('V13.8 TDL-D OPPOSING-TRAIN COARSE + FINE VALIDATION FINISHED\n');
fprintf('============================================================\n');
fprintf('Run mode                        : %s\n',cfg.runMode);
fprintf('Evaluation channel realizations : %d\n',nRep);
fprintf('Total evaluation trials         : %d\n',height(resultsTable));
fprintf('TDL-D delay spread               : %.3f ns RMS\n', ...
    cfg.tdlDelaySpread_s*1e9);
fprintf('TDL internal maximum Doppler     : %.3f Hz (static taps)\n', ...
    cfg.tdlMaximumDopplerShift_Hz);
fprintf('Target / interferer speeds      : %+g / %+g km/h\n', ...
    cfg.MobileSpeed_kmh,cfg.InterfererSpeed_kmh);
fprintf('Target / interferer RootIndex   : %d / %d\n', ...
    cfg.prachRootSeqIdx,cfg.interfererRootIndex);
fprintf('Relative geometry Doppler       : %+.6f Hz at t=0\n', ...
    cfg.relativeGeometryDoppler_Hz);
fprintf('Imposed relative delay          : %.3f effective samples\n', ...
    cfg.interferenceDelayEffSamples);
fprintf('SNR list                        : %s dB\n',mat2str(cfg.snrList_dB));
fprintf('SIR list                        : %s dB\n',mat2str(cfg.sirList_dB));

fprintf('\nPer-SNR / per-SIR metrics:\n');
disp(metricsTable);

fprintf('\nOutput files:\n');
fprintf('  %s\n',rawFile);
fprintf('  %s\n',metricsFile);
fprintf('  %s\n',matFile);
if cfg.makeFigures
    fprintf('  V13_8_TDLD_OpposingTrain_CFO_CoarseFine_3D.png / .fig\n');
    fprintf('  V13_8_TDLD_OpposingTrain_Delay_Strongest_CoarseFine_3D.png / .fig\n');
end
if cfg.makeDiagnosticFigures
    fprintf('  V13_8_TDLD_OpposingTrain_CoarseFine_Diagnostics.png / .fig\n');
end

fprintf('\nModeling notes:\n');
fprintf('  Independent TDL-D static-tap channels were used for both trains.\n');
fprintf('  Opposing-train geometry Doppler was applied externally.\n');
fprintf('  SIR = 10log10(Ptarget/Pinterference): lower SIR means stronger interference.\n');
fprintf('  SNR grid is +5: -1: -15 dB; SIR grid is +10: -1: -10 dB; no Inf cases.\n');
fprintf('  SNR/SIR references use conventional Root-22 target power.\n');
fprintf('  Paired target and Root-0 interferer were not energy-normalized.\n');
fprintf('  Coarse method: full-lag cross-correlation + log-magnitude 3-point parabola.\n');
fprintf('  Fine method: HSR_MP_ACR; delay RMSE is evaluated on the strongest path.\n');
fprintf('  Strongest-path truth is selected from each TDL-D realization pathGains.\n');
fprintf('  No legacy recovery estimator, PDP detector or threshold was used.\n');

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

function [rxFaded,pathGains] = generateTDLRealizationLocal( ...
    tdlChan,txBurst,tailLength,channelFilterDelay)

txWithTail = [txBurst; ...
    zeros(tailLength,size(txBurst,2),'like',txBurst)];

[rxFull,pathGains] = tdlChan(txWithTail);

first = channelFilterDelay+1;
last = first+size(txBurst,1)-1;

if size(rxFull,1)<last
    error('V13_8:TDLOutputTooShort', ...
        ['TDL-D output has %d samples; %d are required after ', ...
         'channel-filter-delay compensation.'],size(rxFull,1),last);
end

rxFaded = rxFull(first:last,:);

end

%% ============================================================
% LOCAL FUNCTION: select strongest path from the actual TDL realization
%% ============================================================

function [strongestIndex,pathPower] = ...
    getStrongestTDLPathLocal(pathGains,nPaths)

if isempty(pathGains) || ~isnumeric(pathGains)
    error('V13_8:InvalidTDLPathGains', ...
        'TDL path-gain output is empty or nonnumeric.');
end

% nrTDLChannel path-gain dimensions:
% time samples x paths x transmit antennas x receive antennas.
if size(pathGains,2)~=nPaths
    error('V13_8:UnexpectedTDLPathGainDimensions', ...
        ['Expected pathGains dimension 2 to contain %d paths; got %d.'], ...
        nPaths,size(pathGains,2));
end

pgPermuted = permute(pathGains,[1 3 4 2]);
pgMatrix = reshape(pgPermuted,[],nPaths);
pathPower = mean(abs(pgMatrix).^2,1);
pathPower = real(pathPower);

if any(~isfinite(pathPower)) || all(pathPower<=0)
    error('V13_8:InvalidTDLPathPower', ...
        'TDL path powers are invalid.');
end

[~,strongestIndex] = max(pathPower);

end

%% ============================================================
% LOCAL FUNCTION: fractional raw-sample delay for multiple antennas
%% ============================================================

function y = fractionalDelayRawLocal(x,delaySamples,guardSamples)

if ~isnumeric(x) || isempty(x) || ndims(x)~=2
    error('V13_8:InvalidDelayInput', ...
        'Fractional-delay input must be a nonempty 2-D matrix.');
end

if ~isfinite(delaySamples) || delaySamples<0
    error('V13_8:InvalidDelay', ...
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
    error('V13_8:DelayLengthMismatch', ...
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
    error('V13_8:InvalidGainSource', ...
        'Unable to derive finite positive shared gains.');
end

rxGain = 1./sqrt(rxPower);
refGain = 1/sqrt(refPower);

end

%% ============================================================
% LOCAL FUNCTION: pad a samples-by-Rx matrix to a required row count
%% ============================================================

function y = padRowsLocal(x,nRows)

if size(x,1)>nRows
    error('V13_8:PadRowsTooShort','Requested row count is shorter than input.');
end

if size(x,1)==nRows
    y = x;
else
    y = [x; zeros(nRows-size(x,1),size(x,2),'like',x)];
end

end

%% ============================================================
% LOCAL FUNCTION: locate target active support
%% ============================================================

function idx = signalSupportWindowLocal(x)

if ~isnumeric(x) || isempty(x) || ndims(x)~=2
    error('V13_8:InvalidSupportInput', ...
        'Support input must be a nonempty samples-by-Rx matrix.');
end

env = max(abs(x),[],2);
peak = max(env,[],'omitnan');

if isempty(peak) || ~isfinite(peak) || peak<=0
    idx = [];
    return;
end

mask = isfinite(env) & env>peak*1e-12;
first = find(mask,1,'first');
last = find(mask,1,'last');

if isempty(first) || isempty(last) || last<first
    idx = [];
else
    idx = (first:last).';
end

end

%% ============================================================
% LOCAL FUNCTION: deterministic unit-variance complex AWGN
%% ============================================================

function noise = makeComplexNoiseLocal(signalSize,noiseVariance,seed)

if ~isscalar(noiseVariance) || ~isfinite(noiseVariance) || noiseVariance<0
    error('V13_8:InvalidNoiseVariance','Noise variance must be finite and nonnegative.');
end

if noiseVariance==0
    noise = complex(zeros(signalSize));
    return;
end

oldState = rng;
restoreRng = onCleanup(@()rng(oldState)); %#ok<NASGU>
rng(seed,'twister');

z = (randn(signalSize)+1i*randn(signalSize))/sqrt(2);
noise = sqrt(noiseVariance).*z;

end

%% ============================================================
% LOCAL FUNCTION: reproducible independent interferer phase
%% ============================================================

function phi = randomPhaseLocal(seed)

oldState = rng;
restoreRng = onCleanup(@()rng(oldState)); %#ok<NASGU>
rng(seed,'twister');
phi = 2*pi*rand;

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

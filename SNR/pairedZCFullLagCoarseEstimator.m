function [delayEffSamples, cfoHz, diagInfo] = ...
    pairedZCFullLagCoarseEstimator(timeData, prachConfig, opts)
% ============================================================
% PAIRED ZC FULL-LAG SUB-BIN COARSE ESTIMATOR
%
% V13.5 independent coarse estimator
%
% 功能：
%   1. 两个分支分别进行完整有符号时延互相关
%   2. 两个分支均允许正、负时延搜索
%   3. 使用对数幅度三点抛物线插值
%   4. 联合计算 CFO 和传播时延
%   5. 输出完整诊断信息
%
% 输入：
%   timeData   : 接收天线数 × 有效采样点数
%   prachConfig: prachSigGen生成的配置
%   opts       : 可选参数结构体
%
% 输出：
%   delayEffSamples : 时延估计（有效采样点）
%   cfoHz           : CFO估计（Hz）
%   diagInfo        : 诊断信息
%
% 坐标关系：
%   d1 = delay + CFO_bins
%   d2 = delay - CFO_bins
%
%   delay = (d1+d2)/2
%   CFO   = (d1-d2)/2
%
% 注意：
%   不调用 time_freq_recovery
%   不调用 time_freq_recovery_V6
%   不调用 prachPowerDelayProfileGen
% ============================================================

if nargin < 3 || isempty(opts)
    opts = struct();
end

%% 1. 输入检查

if ~isnumeric(timeData) || isempty(timeData) || ndims(timeData) > 3
    error('pairedZCFullLagCoarseEstimator:InvalidTimeData', ...
        'timeData must be a nonempty numeric receive-by-samples array.');
end

requiredFields = {'fdRootSeq','fdRootSeq_', ...
    'L_RA','scs_RA','N_CS','bwpInfo'};

for k = 1:numel(requiredFields)

    if ~isfield(prachConfig,requiredFields{k})

        error('pairedZCFullLagCoarseEstimator:MissingConfigField', ...
            'prachConfig is missing field %s.', ...
            requiredFields{k});

    end

end

nRx = size(timeData,1);

timeFlat = reshape(timeData,nRx,[]);

N = prachConfig.bwpInfo.NFFT;

if size(timeFlat,2) ~= N

    error('pairedZCFullLagCoarseEstimator:SampleCountMismatch', ...
        'Expected %d effective samples per Rx, received %d.', ...
        N,size(timeFlat,2));

end

%% 2. 参数设置

opts = setDefault(opts,'MaxAbsCFO_Hz',1500);

opts = setDefault(opts,'MinDelaySamples',0);

opts = setDefault(opts,'MaxDelaySamples',20);

opts = setDefault(opts,'SearchGuardBins',2);

opts = setDefault(opts,'ReturnCorrelations',false);

if ~isscalar(opts.MaxAbsCFO_Hz) || ...
        opts.MaxAbsCFO_Hz < 0 || ...
        ~isfinite(opts.MaxAbsCFO_Hz)

    error('pairedZCFullLagCoarseEstimator:InvalidCFOBound', ...
        'MaxAbsCFO_Hz must be a finite nonnegative scalar.');

end

if ~isscalar(opts.MinDelaySamples) || ...
        ~isscalar(opts.MaxDelaySamples) || ...
        ~isfinite(opts.MinDelaySamples) || ...
        ~isfinite(opts.MaxDelaySamples) || ...
        opts.MaxDelaySamples < opts.MinDelaySamples

    error('pairedZCFullLagCoarseEstimator:InvalidDelayBounds', ...
        'Invalid MinDelaySamples / MaxDelaySamples.');

end

%% 3. 计算频偏坐标间隔

OS = N/prachConfig.L_RA;

cfoBinSpacing_Hz = ...
    prachConfig.scs_RA*1e3/OS;

if ~isfinite(cfoBinSpacing_Hz) || cfoBinSpacing_Hz <= 0

    error('pairedZCFullLagCoarseEstimator:InvalidBinSpacing', ...
        'Unable to derive a valid PRACH CFO-bin spacing.');

end

maxCFOBins = opts.MaxAbsCFO_Hz/cfoBinSpacing_Hz;

%% 4. 完整有符号时延搜索窗口

lagSearchMin = floor( ...
    opts.MinDelaySamples ...
    - maxCFOBins ...
    - opts.SearchGuardBins);

lagSearchMax = ceil( ...
    opts.MaxDelaySamples ...
    + maxCFOBins ...
    + opts.SearchGuardBins);

lagSearchMin = max(-(N-1),lagSearchMin);

lagSearchMax = min(N-1,lagSearchMax);

if lagSearchMin >= lagSearchMax

    error('pairedZCFullLagCoarseEstimator:EmptySearchRegion', ...
        'The configured signed-lag search region is empty.');

end

%% 5. 生成两个配对ZC参考序列

ref1 = ifft(prachConfig.fdRootSeq,N)*sqrt(N);

ref2 = ifft(prachConfig.fdRootSeq_,N)*sqrt(N);

ref1 = ref1(:);

ref2 = ref2(:);

%% 6. 接收信号增益处理

if isfield(opts,'RxGain') && ~isempty(opts.RxGain)

    rxGain = opts.RxGain(:);

    if isscalar(rxGain)

        rxGain = repmat(rxGain,nRx,1);

    end

    if numel(rxGain) ~= nRx || ...
            any(~isfinite(rxGain)) || ...
            any(rxGain <= 0)

        error('pairedZCFullLagCoarseEstimator:InvalidRxGain', ...
            'RxGain must be a positive scalar or one positive value per Rx.');

    end

else

    rxPower = mean(abs(timeFlat).^2,2);

    if any(~isfinite(rxPower)) || any(rxPower <= 0)

        error('pairedZCFullLagCoarseEstimator:InvalidRxPower', ...
            'Cannot derive finite positive receive gains.');

    end

    rxGain = 1./sqrt(rxPower);

end

%% 7. 参考序列增益处理

if isfield(opts,'RefGain') && ~isempty(opts.RefGain)

    refGain = opts.RefGain;

else

    refPower = mean(abs(ref1).^2);

    if ~isfinite(refPower) || refPower <= 0

        error('pairedZCFullLagCoarseEstimator:InvalidRefPower', ...
            'Cannot derive a finite positive reference gain.');

    end

    refGain = 1/sqrt(refPower);

end

if ~isscalar(refGain) || ...
        ~isfinite(refGain) || refGain <= 0

    error('pairedZCFullLagCoarseEstimator:InvalidRefGain', ...
        'RefGain must be a finite positive scalar.');

end

ref1g = ref1*refGain;

ref2g = ref2*refGain;

%% 8. 完整互相关计算

c1Sum = complex(zeros(2*N-1,1));

c2Sum = complex(zeros(2*N-1,1));

for ir = 1:nRx

    rx = timeFlat(ir,:).'*rxGain(ir);

    % Branch 1
    [c1,lags1] = xcorr(rx,ref1g);

    % Branch 2
    [c2,lags2] = xcorr(rx,ref2g);

    c1 = c1(:);

    c2 = c2(:);

    lags1 = lags1(:);

    lags2 = lags2(:);

    if numel(c1) ~= 2*N-1 || ...
            numel(c2) ~= 2*N-1 || ...
            any(lags1 ~= lags2)

        error('pairedZCFullLagCoarseEstimator:UnexpectedCorrelationLength', ...
            'Unexpected full cross-correlation length or lag alignment.');

    end

    c1Sum = c1Sum+c1;

    c2Sum = c2Sum+c2;

    lags = lags1;

end

%% 9. Branch 1完整有符号峰值搜索

[d1Int,d1Sub,p1,idx1,boundary1] = ...
    peakInSignedRange( ...
    c1Sum,lags,lagSearchMin,lagSearchMax);

%% 10. Branch 2完整有符号峰值搜索

[d2Int,d2Sub,p2,idx2,boundary2] = ...
    peakInSignedRange( ...
    c2Sum,lags,lagSearchMin,lagSearchMax);

%% 11. 双分支联合估计

cfoBins = (d1Sub-d2Sub)/2;

delayBins = (d1Sub+d2Sub)/2;

cfoHz = cfoBins*cfoBinSpacing_Hz;

delayEffSamples = delayBins;

%% 12. 计算峰值质量

searchMask = ...
    lags >= lagSearchMin & ...
    lags <= lagSearchMax;

searchMag1 = abs(c1Sum(searchMask));

searchMag2 = abs(c2Sum(searchMask));

%% 13. 输出诊断信息

diagInfo = struct();

% Branch 1
diagInfo.branch1LagInteger = d1Int;

diagInfo.branch1LagSubBin = d1Sub;

diagInfo.branch1PeakValue = p1;

diagInfo.branch1PeakIndex = idx1;

% Branch 2
diagInfo.branch2LagInteger = d2Int;

diagInfo.branch2LagSubBin = d2Sub;

diagInfo.branch2PeakValue = p2;

diagInfo.branch2PeakIndex = idx2;

% Branch power comparison
diagInfo.branch2ToBranch1Peak_dB = ...
    20*log10((p2+eps)/(p1+eps));

% Boundary check
diagInfo.branch1PeakAtSearchBoundary = boundary1;

diagInfo.branch2PeakAtSearchBoundary = boundary2;

diagInfo.anyPeakAtSearchBoundary = ...
    boundary1 || boundary2;

% Final coordinates
diagInfo.cfoBins = cfoBins;

diagInfo.delayBins = delayBins;

diagInfo.cfoBinSpacing_Hz = cfoBinSpacing_Hz;

diagInfo.oversampling = OS;

% Search region
diagInfo.searchLagMin = lagSearchMin;

diagInfo.searchLagMax = lagSearchMax;

diagInfo.maxCFOBins = maxCFOBins;

% Gain information
diagInfo.rxGain = rxGain;

diagInfo.refGain = refGain;

% Peak-to-median ratio
diagInfo.branch1PeakToMedian_dB = ...
    20*log10((p1+eps)/(median(searchMag1)+eps));

diagInfo.branch2PeakToMedian_dB = ...
    20*log10((p2+eps)/(median(searchMag2)+eps));

%% 14. 可选输出完整相关曲线

if opts.ReturnCorrelations

    diagInfo.lags = lags;

    diagInfo.branch1Correlation = c1Sum;

    diagInfo.branch2Correlation = c2Sum;

end

end

%% ============================================================
% LOCAL FUNCTION 1: 默认参数
%% ============================================================

function s = setDefault(s,name,value)

if ~isfield(s,name) || isempty(s.(name))

    s.(name) = value;

end

end

%% ============================================================
% LOCAL FUNCTION 2: 完整有符号峰值搜索 + 亚采样插值
%% ============================================================

function [lagInt,lagSub,peakValue,peakIndex,isBoundary] = ...
    peakInSignedRange(c,lags,lagMin,lagMax)

%% 1. 建立完整有符号搜索区域

mask = lags >= lagMin & lags <= lagMax;

idx = find(mask);

if isempty(idx)

    error('pairedZCFullLagCoarseEstimator:EmptyLagMask', ...
        'No correlation samples fall inside the configured lag interval.');

end

%% 2. 搜索最大相关峰

mag = abs(c(idx));

[peakValue,k] = max(mag);

peakIndex = idx(k);

lagInt = lags(peakIndex);

isBoundary = ...
    (lagInt == lagMin) || ...
    (lagInt == lagMax);

lagSub = lagInt;

%% 3. 对数幅度三点抛物线亚采样插值

if peakIndex > 1 && peakIndex < numel(c)

    ym = log(max(abs(c(peakIndex-1)),realmin));

    y0 = log(max(abs(c(peakIndex)),realmin));

    yp = log(max(abs(c(peakIndex+1)),realmin));

    den = ym-2*y0+yp;

    if isfinite(den) && abs(den) > 1e-12

        delta = 0.5*(ym-yp)/den;

        % 限制插值范围
        delta = max(-0.5,min(0.5,delta));

        lagSub = lagInt+delta;

    end

end

end
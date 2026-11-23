
function out = HSR_MP_ACR(rxWin,txRef,fs,tauCoarseSamples,fCoarseHz,opts)

if nargin < 6 || isempty(opts)
    opts = struct();
end

%% 1. 默认参数

opts = setDefault(opts,'inputCFOCompensated',false);
opts = setDefault(opts,'freqRadiusHz',1500);
opts = setDefault(opts,'delayRadiusSamples',1.5);
opts = setDefault(opts,'maxIter',4);
opts = setDefault(opts,'freqGridN',81);
opts = setDefault(opts,'delayGridN',61);

opts = setDefault(opts,'maxPaths',4);
opts = setDefault(opts,'pathSearchSpanSamples',1.5);
opts = setDefault(opts,'fractionalDelayGridStep',0.05);
opts = setDefault(opts,'pathRelThreshold',0.12);
opts = setDefault(opts,'minPathSeparationSamples',0.20);
opts = setDefault(opts,'ridge',1e-8);

opts = setDefault(opts,'pathOffsetsSamples',[]);
opts = setDefault(opts,'delayInitSamples',[]);

% 兼容旧版参数：
% 旧代码可能将pathSearchSpanSamples设为64。
% 当前假设输入窗口已粗对齐，因此限制在局部搜索范围。
if opts.pathSearchSpanSamples > 4
    warning(['pathSearchSpanSamples原值大于4，', ...
        '当前局部精估计将其调整为1.5个采样点。']);
    opts.pathSearchSpanSamples = 1.5;
end

if opts.pathSearchSpanSamples <= 0
    error('pathSearchSpanSamples必须大于0。');
end

if opts.fractionalDelayGridStep <= 0
    error('fractionalDelayGridStep必须大于0。');
end

if opts.minPathSeparationSamples < 0
    error('minPathSeparationSamples不能小于0。');
end

if opts.maxPaths < 1
    error('maxPaths必须至少为1。');
end

%% 2. 输入检查

y = rxWin;

if isvector(y)
    y = y(:);
end

x = txRef(:);

if isempty(y) || isempty(x)
    error('接收信号和参考信号不能为空。');
end

if size(y,1) ~= length(x)
    error('rxWin与txRef必须具有相同的采样点数。');
end

if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
    error('采样率fs必须为正数。');
end

if any(~isfinite(y(:))) || any(~isfinite(x(:)))
    error('输入信号包含NaN或Inf。');
end

if ~isscalar(tauCoarseSamples) || ~isscalar(fCoarseHz)
    error('粗时延和粗频偏必须为标量。');
end

N = length(x);
Nr = size(y,2);
n = (0:N-1).';

%% 3. 检测分数采样多径

if ~isempty(opts.pathOffsetsSamples)

    % 外部提供相对路径偏移时，跳过自动路径检测
    pathOffsets = sort(opts.pathOffsetsSamples(:).');

    if isempty(opts.delayInitSamples)
        tauCenter = 0;
    else
        tauCenter = opts.delayInitSamples;
    end

    pathLagsDetected = tauCenter + pathOffsets;

else

    [pathOffsets,tauCenter,pathLagsDetected] = ...
        detectFractionalPathOffsets( ...
        y,x,n,fs,fCoarseHz,opts);

end

if isempty(pathOffsets)
    pathOffsets = 0;
    tauCenter = 0;
    pathLagsDetected = 0;
end

pathOffsets = sort(pathOffsets(:).');
L = length(pathOffsets);

%% 4. 交替精化频偏和时延

% deltaF：相对于粗频偏的修正量
% deltaTau：参考路径相对于粗对齐窗口的剩余时延
deltaF = 0;
deltaTau = tauCenter;

history = zeros(opts.maxIter,3);

tauLower = tauCenter - opts.delayRadiusSamples;
tauUpper = tauCenter + opts.delayRadiusSamples;

for it = 1:opts.maxIter

    %% 4.1 固定时延，搜索频偏

    fRad = opts.freqRadiusHz / (2^(it-1));

    fLo = max(-opts.freqRadiusHz,deltaF-fRad);
    fHi = min( opts.freqRadiusHz,deltaF+fRad);

    if fHi > fLo

        deltaF = gridThenContinuousMin( ...
            @(df)objective(deltaTau,df), ...
            fLo,fHi,opts.freqGridN);

    end

    %% 4.2 固定频偏，搜索时延

    tRad = opts.delayRadiusSamples / (2^(it-1));

    tLo = max(tauLower,deltaTau-tRad);
    tHi = min(tauUpper,deltaTau+tRad);

    if tHi > tLo

        deltaTau = gridThenContinuousMin( ...
            @(dt)objective(dt,deltaF), ...
            tLo,tHi,opts.delayGridN);

    end

    %% 4.3 边界检查

    boundaryTol = 1e-3;

    if abs(deltaTau-tauLower) < boundaryTol || ...
       abs(deltaTau-tauUpper) < boundaryTol

        fprintf(['[时延边界警告] iter=%d, ', ...
            'deltaTau=%.6f, lower=%.6f, upper=%.6f\n'], ...
            it,deltaTau,tauLower,tauUpper);

    end

    [costNow,~] = objective(deltaTau,deltaF);

    history(it,:) = [deltaF,deltaTau,costNow];

end

[finalCost,hFinal] = objective(deltaTau,deltaF);

%% 5. 计算首径和最强径

% pathOffsets以检测阶段的最强参考路径为零点。
% 路径系数hFinal的每一行对应一条路径。

finalPathPower = sum(abs(hFinal).^2,2);

[~,strongIdxFinal] = max(finalPathPower);

tauCenterFine = tauCoarseSamples + deltaTau;

tauFirst = tauCenterFine + min(pathOffsets);

tauStrongest = tauCenterFine + pathOffsets(strongIdxFinal);

%% 6. 输出结果

out = struct();

% 为兼容旧版主程序，tauFineSamples定义为最终最强径时延
out.tauFineSamples = tauStrongest;

out.tauCenterSamples = tauCenterFine;
out.tauFirstPathSamples = tauFirst;
out.tauStrongestPathSamples = tauStrongest;

out.freqFineHz = fCoarseHz + deltaF;

out.deltaTauSamples = deltaTau;
out.deltaFreqHz = deltaF;

out.pathOffsetsSamples = pathOffsets;

% 初始路径检测结果，相对于粗对齐窗口
out.pathLagsDetectedSamples = pathLagsDetected;

% 精化后的各路径总时延估计
out.pathLagsTotalSamples = ...
    tauCenterFine + pathOffsets;

out.pathCoefficients = hFinal;

% 归一化残差为每根接收天线残差的平均值
out.normalizedResidual = finalCost;

out.history = history;
out.numReceiveAntennas = Nr;
out.numPaths = L;

out.delayGridStepSamples = ...
    opts.fractionalDelayGridStep;

%% 嵌套目标函数
% 最小化所有接收天线上的多径拟合残差。
% 嵌套函数可直接访问本主函数中的变量。

    function [cost,hMat] = objective(tauTry,dfTry)

        if opts.inputCFOCompensated
            % 输入信号已补偿粗频偏，只需补偿剩余频偏
            fRemove = dfTry;
        else
            % 输入信号未补偿粗频偏
            fRemove = fCoarseHz + dfTry;
        end

        rot = exp(-1j*2*pi*fRemove*n/fs);

        z = bsxfun(@times,y,rot);

        S = complex(zeros(N,L));

        for ell = 1:L

            shiftNow = tauTry + pathOffsets(ell);

            S(:,ell) = fractionalShift(x,n,shiftNow);

        end

        % 联合估计各条路径的复增益
        G = S'*S;

        regVal = opts.ridge * ...
            max(real(trace(G))/max(L,1),eps);

        Greg = G + regVal*eye(L);

        hMat = Greg \ (S'*z);

        err = z-S*hMat;

        cost = 0;

        for rr = 1:Nr

            cost = cost + ...
                sum(abs(err(:,rr)).^2) / ...
                (sum(abs(z(:,rr)).^2)+eps);

        end

        cost = real(cost/Nr);

    end

end


%% ============================================================
% 分数采样时延多径检测：匹配追踪 + 联合最小二乘
% =============================================================
function [offsets,centerLag,selectedLags] = ...
    detectFractionalPathOffsets(y,x,n,fs,fCoarseHz,opts)

N = length(x);

%% 1. 粗频偏补偿

if opts.inputCFOCompensated
    fRemove = 0;
else
    fRemove = fCoarseHz;
end

rot = exp(-1j*2*pi*fRemove*n/fs);

z = bsxfun(@times,y,rot);

%% 2. 构造分数时延模板

span = opts.pathSearchSpanSamples;
step = opts.fractionalDelayGridStep;

delayGrid = (-span:step:span).';

if isempty(delayGrid) || delayGrid(end) < span-1e-10
    delayGrid = [delayGrid;span];
end

K = length(delayGrid);

D = complex(zeros(N,K));

for k = 1:K
    D(:,k) = fractionalShift(x,n,delayGrid(k));
end

templateEnergy = sum(abs(D).^2,1).';

%% 3. 匹配追踪

R = z;

selectedIdx = [];
firstScore = [];

for p = 1:opts.maxPaths

    corrMat = D'*R;

    score = sum(abs(corrMat).^2,2) ./ ...
        (templateEnergy+eps);

    % 避免重复选择相邻网格点
    if ~isempty(selectedIdx)

        alreadySelected = delayGrid(selectedIdx);

        for q = 1:length(alreadySelected)

            nearMask = ...
                abs(delayGrid-alreadySelected(q)) < ...
                opts.minPathSeparationSamples;

            score(nearMask) = -Inf;

        end

    end

    [bestScore,bestIdx] = max(score);

    if ~isfinite(bestScore) || bestScore <= 0
        break;
    end

    if isempty(firstScore)

        firstScore = bestScore;

    else

        if bestScore < opts.pathRelThreshold*firstScore
            break;
        end

    end

    selectedIdx(end+1) = bestIdx; %#ok<AGROW>

    % 每增加一条路径，联合重新拟合路径系数
    S = D(:,selectedIdx);

    G = S'*S;

    regVal = opts.ridge * ...
        max(real(trace(G))/size(G,1),eps);

    h = (G+regVal*eye(size(G,1))) \ (S'*z);

    R = z-S*h;

end

%% 4. 兜底：至少保留一条路径

if isempty(selectedIdx)

    score = sum(abs(D'*z).^2,2) ./ ...
        (templateEnergy+eps);

    [~,bestIdx] = max(score);

    selectedIdx = bestIdx;

end

%% 5. 按时延排序

[selectedLags,order] = sort(delayGrid(selectedIdx));

selectedIdx = selectedIdx(order);

S = D(:,selectedIdx);

G = S'*S;

regVal = opts.ridge * ...
    max(real(trace(G))/size(G,1),eps);

h = (G+regVal*eye(size(G,1))) \ (S'*z);

%% 6. 选择最强路径作为参考路径

pathPower = sum(abs(h).^2,2) .* ...
    sum(abs(S).^2,1).';

[~,strongestIdx] = max(pathPower);

centerLag = selectedLags(strongestIdx);

% 以最强路径为零点，构造相对路径偏移
offsets = (selectedLags-centerLag).';

end


%% ============================================================
% 分数采样时延模板
% =============================================================
function s = fractionalShift(x,n,tau)

% 模拟：
% s[n] = x[n-tau]
%
% 使用三次样条插值，支持非整数采样点时延。
% 超出信号边界的部分补零。

s = interp1(n,x,n-tau,'spline',0);

s = s(:);

end


%% ============================================================
% 网格搜索 + 连续优化
% =============================================================
function xBest = gridThenContinuousMin(fun,lo,hi,nGrid)

if hi <= lo
    xBest = lo;
    return;
end

nGrid = max(3,round(nGrid));

xx = linspace(lo,hi,nGrid);

yy = zeros(size(xx));

for k = 1:length(xx)
    yy(k) = fun(xx(k));
end

[yGridBest,idx] = min(yy);

xBest = xx(idx);

if idx == 1

    a = xx(1);
    b = xx(2);

elseif idx == length(xx)

    a = xx(end-1);
    b = xx(end);

else

    a = xx(idx-1);
    b = xx(idx+1);

end

if b > a

    opt = optimset('Display','off','TolX',1e-5);

    [xCont,yCont] = fminbnd(fun,a,b,opt);

    if yCont < yGridBest
        xBest = xCont;
    end

end

end


%% ============================================================
% 设置默认参数
% =============================================================
function s = setDefault(s,fieldName,value)

if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = value;
end

end

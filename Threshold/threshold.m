clc;
clear;
close all;


%% ================= 参数配置 =================

simConfig = struct('nTx',1,...
    'nRx',2,...
    'puschBW',100,...
    'puschSCS',30,...
    'prachFormat','F0',...
    'prachZeroConfigZone',13,...
    'prachRootSeqIdx',22,...
    'channelType','TDL',...
    'snrdB',0:-1:-15,...
    'simEpoch',1,...
    'simAlgo',0,...
    'CarrierFrequency',2.1e9,...
    'MobileSpeed',500/3.6,...
    'TrainPosition',-1000,...
    'TrackDistance',100,...
    'BaseHeight',30,...
    'TrainAntennaHeight',5,...
    'TransmissionDelay',0,...
    'c',physconst('lightspeed'));



%% 阈值因子搜索范围

thrFactor = 0.1:0.05:0.8;


falseAlarm=zeros(length(thrFactor),1);

missDetection=zeros(length(thrFactor),1);



%% ================= 开始仿真 =================


for thrIdx=1:length(thrFactor)


    currentFactor=thrFactor(thrIdx);


    false_alarm_snr=zeros(length(simConfig.snrdB),1);

    miss_detection_snr=zeros(length(simConfig.snrdB),1);



    for snrIdx=1:length(simConfig.snrdB)


        currentSNR=simConfig.snrdB(snrIdx);


        fa=[];

        md=[];



        for epoch=1:simConfig.simEpoch



            %% PRACH生成

            [txSig,prachConfig]=prachSigGen(simConfig);



            %% =====================
            % 铁路Doppler模型
            %% =====================

            fs=prachConfig.bwpInfo.SamplingRate;

            t=(0:length(txSig)-1)'/fs;


            [fd,~,~]=railwayDopplerModel(...
                simConfig.TrainPosition,...
                simConfig.MobileSpeed,...
                simConfig.TrackDistance,...
                simConfig.BaseHeight-simConfig.TrainAntennaHeight,...
                simConfig.CarrierFrequency,...
                t);


            mobileMaxDoppler=max(abs(fd));



            %% =====================
            % TDL-D信道
            %% =====================


            tdlChannel=nrTDLChannel;


            tdlChannel.DelayProfile='TDL-D';

            tdlChannel.DelaySpread=100e-9;


            tdlChannel.SampleRate=...
                prachConfig.bwpInfo.SamplingRate;


            tdlChannel.MaximumDopplerShift=...
                mobileMaxDoppler;


            tdlChannel.MIMOCorrelation='Low';


            tdlChannel.NumTransmitAntennas=1;

            tdlChannel.NumReceiveAntennas=2;



            [tdlOut,~,~]=tdlChannel(txSig);



            %% 接收信号

            rxSig=awgn(tdlOut,currentSNR);



            %% 噪声信号

            noiseSig=zeros(size(tdlOut));

            rxNoiseSig=awgn(noiseSig,currentSNR);



            %% =====================
            % PRACH处理
            %% =====================


            timeData=prachFrontendProcessing(...
                rxSig,prachConfig,simConfig);


            noiseData=prachFrontendProcessing(...
                rxNoiseSig,prachConfig,simConfig);



            delayProfile=prachPowerDelayProfileGen(...
                timeData,...
                prachConfig,...
                simConfig);


            noiseProfile=prachPowerDelayProfileGen(...
                noiseData,...
                prachConfig,...
                simConfig);



            %% 检测

            detResults=prachThr4F0(...
                currentFactor,...
                currentSNR,...
                delayProfile,...
                prachConfig);


            falResults=fa_cal(...
                currentFactor,...
                currentSNR,...
                noiseProfile,...
                prachConfig);



            md=[md;detResults.miss_detection];

            fa=[fa;falResults.false_alarm];



        end



        false_alarm_snr(snrIdx)=...
            sum(fa)/simConfig.simEpoch;


        miss_detection_snr(snrIdx)=...
            sum(md)/simConfig.simEpoch;



    end



    falseAlarm(thrIdx)=...
        mean(false_alarm_snr);


    missDetection(thrIdx)=...
        mean(miss_detection_snr);



    fprintf('Finish HSR factor = %.2f\n',...
        currentFactor);



end



%% ================= 结果 =================


figure;

plot(thrFactor,falseAlarm,'r-o',...
    'LineWidth',1.5,...
    'MarkerSize',6);

hold on;

plot(thrFactor,missDetection,'b-s',...
    'LineWidth',1.5,...
    'MarkerSize',6);


grid on;

xlabel('Threshold factor');
ylabel('Probability');

legend('False alarm probability',...
       'Miss detection probability',...
       'Location','northwest');


set(gca,'FontSize',12);

title('Optimization of Threshold Factor for HSR PRACH Detection');




save HSR_threshold_factor.mat...
    thrFactor...
    falseAlarm...
    missDetection

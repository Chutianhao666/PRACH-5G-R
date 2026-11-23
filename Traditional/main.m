clc;
clear;
close all;

% 使用主序列+互补序列

%% Parameter Configuration
simConfig = struct('nTx', 1, ...
    'nRx', 2, ...
    'puschBW', 100, ...
    'puschSCS', 30, ...
    'prachFormat', 'F0', ...
    'prachZeroConfigZone', 13, ...
    'prachRootSeqIdx', 22, ...
    'channelType', 'TDL', ...
    'snrdB', 0:-1:-15, ...
    'simEpoch', 1, ...
    'simAlgo', 1, ...
    'CarrierFrequency', 2.1e9, ...
    'MobileSpeed', 500/3.6, ...
    'TrainPosition',-1000,...
    'TrackDistance',100,...
    'BaseHeight',30,...
    'TrainAntennaHeight',5,...
    'TransmissionDelay', 0, ...
    'c', physconst('lightspeed'));




%% Simulation Starts

detection_snr = zeros(length(simConfig.snrdB),1);

for snrIdx = 1:length(simConfig.snrdB)
    currentSNR = simConfig.snrdB(snrIdx);


    de = [];

    for epoch = 1:simConfig.simEpoch
        [txSig, prachConfig] = prachSigGen(simConfig);
        fs=prachConfig.bwpInfo.SamplingRate;
        t=(0:length(txSig)-1)'/fs;

        [fd,rangeM,propagationDelay]=railwayDopplerModel(...
            simConfig.TrainPosition,...
            simConfig.MobileSpeed,...
            simConfig.TrackDistance,...
            simConfig.BaseHeight-simConfig.TrainAntennaHeight,...
            simConfig.CarrierFrequency,...
            t);


        maxDoppler=max(abs(fd));

        txSig = [zeros(simConfig.TransmissionDelay,simConfig.nTx); txSig(1:end-simConfig.TransmissionDelay,:)];


        % pass through wireless channel
        if strcmp(simConfig.channelType, 'TDL')
            mobileMaxDoppler = simConfig.MobileSpeed*simConfig.CarrierFrequency/simConfig.c;
            % railwayDopplerShift = railwayDopplerModel( ...
            %     simConfig.ElevationAngle, simConfig.SatelliteAltitude, ...
            %     simConfig.MobileAltitude, simConfig.CarrierFrequency);
     
            tdlChannel = nrTDLChannel;
            tdlChannel.DelayProfile = 'TDL-D';
            tdlChannel.DelaySpread = 100e-9;
            tdlChannel.SampleRate = prachConfig.bwpInfo.SamplingRate;
            tdlChannel.MaximumDopplerShift = maxDoppler;
            %tdlChannel.SatelliteDopplerShift = railwayDopplerShift;

            tdlChannel.MIMOCorrelation = 'Low';
            tdlChannel.NumTransmitAntennas = 1;
            tdlChannel.NumReceiveAntennas = 2;
            % ntnChannel.RandomStream = 'mt19937ar with seed';
            % ntnChannel.Seed = simConfig.seed;


            % rng(simConfig.seed);
            tdlChannelInfo = info(tdlChannel);
            [tdlOut, tdlPathGains, tdlSampleTimes] = tdlChannel(txSig);

            rxSig = awgn(tdlOut,currentSNR);
            
            % rxSig = tdlOut;
            % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));



        elseif strcmp(simConfig.channelType, 'AWGN')
            rxSig = txSig;
            for nrx = 1:simConfig.nRx
                rxSig(:,nrx) = sum(txSig,2);
            end

            rxSig = awgn(rxSig,currentSNR);
            % figure; subplot(211); plot(abs(rxSig(:,1))); subplot(212); plot(abs(rxSig(:,2)));

        end


        timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig);

        delayProfile = prachPowerDelayProfileGen(timeData, prachConfig, simConfig);

        detResults = prachDetection(0.33, currentSNR, delayProfile, prachConfig);


        de = [de; detResults.detection];

        if epoch == simConfig.simEpoch
            detection_snr(snrIdx) = sum(de)/simConfig.simEpoch;
            fprintf('SNR = %2.1f dB, PRACH_Detection_Rate: %1.3f\n',...
                         currentSNR, detection_snr(snrIdx));
        end

        if epoch == simConfig.simEpoch && snrIdx == length(simConfig.snrdB)
            figure; plot(simConfig.snrdB, detection_snr, 'r-*');
        end

    end % end of epoch
end % end of snr





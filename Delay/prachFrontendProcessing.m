function timeData = prachFrontendProcessing(rxSig, prachConfig, simConfig)
%PRACHFRONTPROCESSING Remove CP, frequency-shift and decimate PRACH samples.
% Input:  samples-by-receive-antennas raw waveform.
% Output: receive-antennas-by-effective-samples-by-repetitions.
%
% The FIR group delay is expressed on the DECIMATED sample grid. If a
% short burst does not contain enough trailing samples to flush the FIR,
% zeros are appended BEFORE filtering (not after filtering), so the
% retained samples are actual filter outputs rather than artificial zeros.

nRx = simConfig.nRx;
NCP_RA = prachConfig.NCP_RA;
L_RA = prachConfig.L_RA;
Nrep = prachConfig.Nrep;
nFFT = prachConfig.nFFT;
K = prachConfig.K;
k1 = prachConfig.k1;
k_bar = prachConfig.k_bar;

if ~isnumeric(rxSig) || isempty(rxSig) || ~ismatrix(rxSig)
    error('prachFrontendProcessing:InvalidInput', ...
        'rxSig must be a nonempty numeric samples-by-antennas matrix.');
end
if isvector(rxSig) && nRx == 1
    rxSig = rxSig(:);
end
if size(rxSig,2) ~= nRx
    error('prachFrontendProcessing:ReceiveAntennaMismatch', ...
        'simConfig.nRx=%d but rxSig has %d columns.',nRx,size(rxSig,2));
end
if size(rxSig,1) <= NCP_RA
    error('prachFrontendProcessing:InsufficientSamples', ...
        'rxSig does not contain samples beyond the configured PRACH CP.');
end
if ~isfinite(nFFT) || nFFT <= 0 || mod(nFFT,Nrep) ~= 0
    error('prachFrontendProcessing:InvalidFFTSize', ...
        'Configured PRACH FFT size/repetition count is invalid.');
end

% Work in antennas-by-time layout.
rx_prach_sig = rxSig.';
rx_prach_sig = rx_prach_sig(:,NCP_RA+1:end);

% Translate the PRACH allocation to baseband.
nRaw = size(rx_prach_sig,2);
mix = exp(-1i*2*pi*(K*k1+k_bar)*(0:nRaw-1)/nFFT);
rx_prach_sig_fs = bsxfun(@times,rx_prach_sig,mix);

if L_RA == 839
    effectiveFFT = prachConfig.bwpInfo.NFFT;
else
    effectiveFFT = nFFT;
end
timeData = complex(zeros(nRx,effectiveFFT,Nrep));

for nrx = 1:nRx
    if L_RA == 839
        % A fresh decimator object is created for each antenna so that
        % System-object filter state cannot leak from one antenna to another.
        decFilter = cascadeDecimatorGen(simConfig,prachConfig);

        switch prachConfig.mu
            case 0  % 15 kHz PUSCH numerology
                if prachConfig.scs_RA == 1.25
                    D = 12;
                    decimator = decFilter.cascadeDecimator_12;
                elseif prachConfig.scs_RA == 5
                    D = 3;
                    decimator = decFilter.cascadeDecimator_3;
                else
                    error('prachFrontendProcessing:UnsupportedRACHSpacing', ...
                        'Unsupported PRACH SCS %.6g for mu=0.',prachConfig.scs_RA);
                end

            case 1  % 30 kHz PUSCH numerology
                if prachConfig.scs_RA == 1.25
                    D = 24;
                    decimator = decFilter.cascadeDecimator_24;
                elseif prachConfig.scs_RA == 5
                    D = 6;
                    decimator = decFilter.cascadeDecimator_6;
                else
                    error('prachFrontendProcessing:UnsupportedRACHSpacing', ...
                        'Unsupported PRACH SCS %.6g for mu=1.',prachConfig.scs_RA);
                end

            otherwise
                error('prachFrontendProcessing:UnsupportedNumerology', ...
                    'Unsupported PUSCH numerology mu=%d.',prachConfig.mu);
        end

        if mod(nFFT,D) ~= 0
            error('prachFrontendProcessing:NonintegerDecimationLength', ...
                'PRACH FFT size %d is not divisible by decimation factor %d.',nFFT,D);
        end

        nExpected = (nFFT/D)*Nrep;
        if (nFFT/D) ~= prachConfig.bwpInfo.NFFT
            error('prachFrontendProcessing:EffectiveFFTSizeMismatch', ...
                'Decimated FFT size %d does not match BWP NFFT %d.', ...
                nFFT/D,prachConfig.bwpInfo.NFFT);
        end
        idxStart = decFilter.delay + 1;
        idxEnd = decFilter.delay + nExpected;

        % Ensure the FIR has enough raw input samples to emit every sample
        % through idxEnd after group-delay removal. Padding is done BEFORE
        % the decimator so its transient response is handled correctly.
        rawRow = rx_prach_sig_fs(nrx,:);
        minRawLength = D*idxEnd + 1;
        if numel(rawRow) < minRawLength
            rawRow = [rawRow, zeros(1,minRawLength-numel(rawRow),'like',rawRow)];
        end

        decimated = decimator(rawRow(:));
        decimated = decimated(:).';

        if numel(decimated) < idxEnd
            error('prachFrontendProcessing:DecimatorOutputTooShort', ...
                ['Decimator produced %d samples, but %d are required ', ...
                 'after FIR group-delay compensation.'],numel(decimated),idxEnd);
        end

        prach_sig_dec = decimated(idxStart:idxEnd);

    elseif L_RA == 139
        % The short-sequence formats do not use this decimation chain.
        nExpected = nFFT*Nrep;
        rawRow = rx_prach_sig_fs(nrx,:);
        if numel(rawRow) < nExpected
            rawRow = [rawRow, zeros(1,nExpected-numel(rawRow),'like',rawRow)];
        end
        prach_sig_dec = rawRow(1:nExpected);

    else
        error('prachFrontendProcessing:UnsupportedSequenceLength', ...
            'Unsupported PRACH sequence length L_RA=%d.',L_RA);
    end

    if numel(prach_sig_dec) ~= (size(timeData,2)*Nrep)
        error('prachFrontendProcessing:OutputLengthMismatch', ...
            'Front-end returned %d samples; expected %d.', ...
            numel(prach_sig_dec),size(timeData,2)*Nrep);
    end

    timeData(nrx,:,:) = reshape(prach_sig_dec,[],Nrep);
end

end

%% Auxiliary Function
function [decFilter] = cascadeDecimatorGen(simConfig,prachConfig)
%
% Design of Cascade Filter, "FlexRAN L1 Algorithm Description" 5.5.3.
%
% cascadeDecimator_24: decimation 24, 4*3*2
% cascadeDecimator_12: decimation 12, 4*3
% cascadeDecimator_6: decimation 6, 3*2
% cascadeDecimator_3: decimation 3
%
% 1st decimator
% Fs = 122.88;                % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 30.18;              % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff1 = firpm(N, Fo, Ao, W, {dens});
% firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
% 
% % 2nd decimator
% Fs = 30.72;                 % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 9.7;                % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff2 = firpm(N, Fo, Ao, W, {dens});
% firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
% 
% % 3rd decimator
% Fs = 10.24;                 % Sampling Frequency
% Fpass = 0.54;               % Passband Frequency
% Fstop = 2.56;               % Stopband Frequency
% Dpass = 0.0057563991496;    % Passband Ripple
% Dstop = 0.001;              % Stopband Attenuation
% dens  = 20;                 % Density Factor
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
% coeff3 = firpm(N, Fo, Ao, W, {dens});
% firdecim3 = dsp.FIRDecimator(2,'Numerator',coeff3);
% 
% cascadeDecimator.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
% cascadeDecimator.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
% 
% decFilter.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
% decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
% decFilter.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
% decFilter.cascadeDecimator_3 = cascade(firdecim2);

samplingrate = prachConfig.bwpInfo.SamplingRate/10^6;
Fpass = 0.54;               % Passband Frequency 839-->864-->864*1.25KHz/2=1.08MHz/2=0.54MHz
Dpass = 0.0057563991496;    % Passband Ripple
Dstop = 0.001;              % Stopband Attenuation
dens  = 20;                 % Density Factor

L_RA = prachConfig.L_RA; 
scs_RA = prachConfig.scs_RA;
switch prachConfig.mu
    case 0  % 15KHz, sampling rate can be [15.360       30.720       61.440]MHz, 12x, 3x downsampling
%         % 1st decimator
%         Fs = 122.88;                % Sampling Frequency
%         Fpass = 0.54;               % Passband Frequency
%         Fstop = Fs/4-0.54;              % Stopband Frequency
%         Dpass = 0.0057563991496;    % Passband Ripple
%         Dstop = 0.001;              % Stopband Attenuation
%         dens  = 20;                 % Density Factor
%         [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
%         coeff1 = firpm(N, Fo, Ao, W, {dens});
%         firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
%         
%         % 2nd decimator
%         Fs = Fs/4;                 % Sampling Frequency
%         Fpass = 0.54;               % Passband Frequency
%         Fstop = Fs/3-0.54;                % Stopband Frequency
%         Dpass = 0.0057563991496;    % Passband Ripple
%         Dstop = 0.001;              % Stopband Attenuation
%         dens  = 20;                 % Density Factor
%         [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
%         coeff2 = firpm(N, Fo, Ao, W, {dens});
%         firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
%         
%         decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
%         decFilter.cascadeDecimator_3 = cascade(firdecim2);
        if samplingrate == 61.44   % 4x, 3x,  FFT size =4096
            % 1st decimator
            bold =4;
            Fs = samplingrate;          % Sampling Frequency, 61.44
            targetFs = Fs/4;            % target sampling rate, 15.36MHz
            Fstop = targetFs-Fpass;     % Stopband Frequency, 14.82MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
            
             % 2nd decimator
            bold =3;
            Fs = targetFs;                % Sampling Frequency, 15.36MHz
            targetFs = Fs/3;              % target sampling rate,5.12MHz
            Fstop = targetFs/2;           % Stopband Frequency, 2.56MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
            decFilter.cascadeDecimator_3 = cascade(firdecim2); 
        elseif  samplingrate ==30.72   % 4x, 3x, FFT size = 2048
            % 1st decimator
            bold =4;
            Fs = samplingrate;          % Sampling Frequency
            targetFs = Fs/4;            % target sampling rate, 7.68MHz
            Fstop = targetFs-Fpass;     % Stopband Frequency, 7.14MHz
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
           
             % 2nd decimator
            bold =3;
            Fs = targetFs;          % Sampling Frequency, 7.68MHz
            targetFs = Fs/3;        % target sampling rate, 2.56MHz
            Fstop = targetFs/2;     % Stopband Frequency, 1.28MHz
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            decFilter.cascadeDecimator_12 = cascade(firdecim1,firdecim2);
            decFilter.cascadeDecimator_3 = cascade(firdecim2); 
        else
            error('cascadeDecimatorGen:UnsupportedSamplingRate', ...
                'Unsupported mu=0 sampling rate %.3f MHz.',samplingrate);
        end
        
        Fs1 = samplingrate;       % MHz
        Fs2 = Fs1/4;               % MHz

        % Linear-phase FIR group delay in microseconds.
        delay12_us = (length(coeff1)-1)/2/Fs1 + ...
                     (length(coeff2)-1)/2/Fs2;
        delay3_us = (length(coeff2)-1)/2/Fs2;

        if scs_RA == 1.25
            decimationFactor = 12;
            delayUs = delay12_us;
        else
            decimationFactor = 3;
            delayUs = delay3_us;
        end

        outputRateMHz = samplingrate/decimationFactor;
        decFilter.delayExactSamples = delayUs*outputRateMHz;
        decFilter.delay = max(0,round(decFilter.delayExactSamples));
    case 1  % 30KHz
        if samplingrate == 122.88
            % 1st decimator  4x
            Fs = 122.88;                % Sampling Frequency
            targetFs = Fs/4;            % target sampling rate, 30.72MHz
            Fpass = 0.54;               % Passband Frequency
            Fstop = targetFs - Fpass;    % Stopband Frequency,  30.18MHz
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff1 = firpm(N, Fo, Ao, W, {dens});
            firdecim1 = dsp.FIRDecimator(4,'Numerator',coeff1);
            
            % 2nd decimator    3x
            Fs = 30.72;                 % Sampling Frequency
            Fpass = 0.54;               % Passband Frequency
            Fstop = 9.7;                % Stopband Frequency
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff2 = firpm(N, Fo, Ao, W, {dens});
            firdecim2 = dsp.FIRDecimator(3,'Numerator',coeff2);
            
            % 3rd decimator   2x
            Fs = 10.24;                 % Sampling Frequency
            Fpass = 0.54;               % Passband Frequency
            Fstop = 2.56;               % Stopband Frequency
%             Dpass = 0.0057563991496;    % Passband Ripple
%             Dstop = 0.001;              % Stopband Attenuation
%             dens  = 20;                 % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
            coeff3 = firpm(N, Fo, Ao, W, {dens});
            firdecim3 = dsp.FIRDecimator(2,'Numerator',coeff3);
            
            decFilter.cascadeDecimator_24 = cascade(firdecim1,firdecim2,firdecim3);
            decFilter.cascadeDecimator_6 = cascade(firdecim2,firdecim3);
        else
            error('cascadeDecimatorGen:UnsupportedSamplingRate', ...
                'Unsupported mu=1 sampling rate %.3f MHz.',samplingrate);
        end
        Fs1 = samplingrate;       % MHz
        Fs2 = Fs1/4;               % MHz
        Fs3 = Fs2/3;               % MHz

        % Linear-phase FIR group delay in microseconds.
        delay24_us = (length(coeff1)-1)/2/Fs1 + ...
                     (length(coeff2)-1)/2/Fs2 + ...
                     (length(coeff3)-1)/2/Fs3;
        delay6_us = (length(coeff2)-1)/2/Fs2 + ...
                    (length(coeff3)-1)/2/Fs3;

        if scs_RA == 1.25
            decimationFactor = 24;
            delayUs = delay24_us;
        else
            decimationFactor = 6;
            delayUs = delay6_us;
        end

        outputRateMHz = samplingrate/decimationFactor;
        decFilter.delayExactSamples = delayUs*outputRateMHz;
        decFilter.delay = max(0,round(decFilter.delayExactSamples));

end
end


clear
close all

% We're here to filter the data

load EMG_Validation.mat

% Field names index
Names = fieldnames(EMG);

% Creating a structure to store the filtered data
EMGfilt = EMG;

% Pulling the sampling rate and setting up a time vector
Fs = EMG.HO.posture.Fs(1,1);

% Filters design
% Highpass, cutoff 30Hz

Fp = 30;    % Passband freq in Hz
Fst = 25;   % Stopband freq in Hz
Ap = 1;     % Passband ripple in dB
Ast = 95;   % Stopband attenuation dB

df = designfilt('highpassfir','PassbandFrequency',Fp, ...
    'StopbandFrequency',Fst,'PassbandRipple',Ap, ...
    'StopbandAttenuation',Ast,'SampleRate',Fs);

% Check filter response I'm not doing that for all these
%hfvt = fvtool(df,'Fs',Fs,'FrequencyScale','log','FrequencyRange', ...
%    'Specify freq. vector','FrequencyVector',F);

% Lowpass, cutoff 480Hz

Fp2 = 480;    % Passband freq in Hz
Fst2 = 500;   % Stopband freq in Hz
Ap2 = 1;     % Passband ripple in dB
Ast2 = 95;   % Stopband attenuation dB

df2 = designfilt('lowpassfir','PassbandFrequency',Fp2, ...
    'StopbandFrequency',Fst2,'PassbandRipple',Ap2, ...
    'StopbandAttenuation',Ast2,'SampleRate',Fs);

% Check filter response
%hfvt2 = fvtool(df2,'Fs',Fs,'FrequencyScale','log','FrequencyRange', ...
%    'Specify freq. vector','FrequencyVector',F);

% I don't need to plot this for every posture so I'm commenting it out
[P,F] = pwelch(EMG.HO.posture.Data(3,:),ones(8192,1),8192/2,8192,Fs,'power');

figure(1);
y1 = loglog(F,P);
ylim([0,1.5*10^-8])
title('hand open: EMG1')
hold on



for j = 1:11
    figure('Name',Names{j})
    for i=1:8 
        t = EMG.(Names{j}).posture.Time;

        % Applying the HP filter
        HPdata = filter(df,EMG.(Names{j}).posture.Data(i,:));
        
        % Applying the Lowpass filter
        LPdata = filter(df2,HPdata);

        % Plotting unfiltered and filtered signals
        subplot(4,2,i)
        plot(t,EMG.(Names{j}).posture.Data(i,:),t,LPdata)
        EMGfilt.(Names{j}).posture.Data(i,:)=LPdata;
        ylim([-0.003, 0.003])
    end
end

figure(1)
[P2,F2] = pwelch(EMGfilt.HO.posture.Data(3,:),ones(8192,1),8192/2,8192,Fs,'power');
y1 = loglog(F2,P2);
ylim([0,1.5*10^-8])
title('hand open: EMG1')

% Plotting the signal back after filtering it though doesn't show the kind
% of power drop I was expecting, maybe that's just because most of the
% power is already in the region I'm passing or maybe I made a mistake
% somewher, but I don't know what it would be

% I still need to put loops around the bits that apply the filters, and
% store the filtered data in a structure so I can call it later, but I've
% got working HP and LP filters I can apply to a signal


% After looking at a few of the data streams I don't see spikes at 60 Hz or
% any of the harmonics, so I'm leaving the notch filters out, if it becomes
% a problem I'll bring it back and put it in a loop to make a few of them 
% and apply them to the data

% Notch @ 60 and harmonics (120, 180, 240, 300, 360, 420 up to pass band
% 
% Fp3 = 50;    % Passband freq in Hz
% Fst3 = 55;   % Stopband freq in Hz
% Fst4 = 65;   % Stopband freq 2
% Fp4 = 70;    % Passband freq 2
% Ap3 = 1;     % Passband ripple in dB
% Ast3 = 50;   % Stopband attenuation dB
% 
% df3 = designfilt('bandstopfir','PassbandFrequency1',Fp3, ...
%     'StopbandFrequency1',Fst3, 'StopbandFrequency2', ...
%     Fst4, 'Passbandfrequency2', Fp4, 'PassbandRipple1',Ap3, ...
%     'StopbandAttenuation',Ast3,'SampleRate',Fs);
% 
% % Check filter response
% hfvt3 = fvtool(df3,'Fs',Fs,'FrequencyScale','log','FrequencyRange', ...
%     'Specify freq. vector','FrequencyVector',F);

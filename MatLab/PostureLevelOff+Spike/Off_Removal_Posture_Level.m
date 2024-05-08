clear
close all

% We're here to remove "off" sections from data

load Filtered_EMG.mat
Sum1 = zeros(8,167761);
% Field names index

% Pulling the sampling rate and setting up a time vector
Fs = EMGfilt.HO.posture.Fs(1,1);
t = EMGfilt.HO.posture.Time(1,:);


Names = fieldnames(EMGfilt);

for j=1:11
    if j == 4
        Sum1(j,:) = sum(abs(EMGfilt.(Names{j}).posture.Data))-abs(EMGfilt.(Names{j}).posture.Data(1,:))-abs(EMGfilt.(Names{j}).posture.Data(3,:))-abs(EMGfilt.(Names{j}).posture.Data(5,:));
        % Posture 4 (RD) is a weak signal with an antagonist and some noisy 
        % channels that drown out the "On" signal, So I'm subtracting those
        % two for posture level off/spike removal
    else
        Sum1(j,:) = sum(abs(EMGfilt.(Names{j}).posture.Data),1);
        % Summing up all the signals for each posture
    end
end

SmoothBINsize = 2000; % bin size in ms for smoothing the data
Smoothed = zeros(1, 167761);
OnOffIndex = ones(11,167761);

for k=2:11
    % To run off removal and spike removal at the posture level
    % This version skips posture 1 (off) entirely
    Smoothed = smoothdata(Sum1(k,:),'gaussian',(SmoothBINsize*1000/Fs));

    % Off removal and spike detection at the poosture level requires using
    % mean and standard deviation for detection instead of just setting
    % values

    Mean = mean(Smoothed(1,:));
    Std = std(Smoothed(1,:));

    for i = 1: 167761

% setting a threshold for when the combined signal is "on" vs. "off", and a
% threshold for spike removal

       if Smoothed(1,i) > Mean + 2.2*Std
            OnOffIndex(k,i)=0;
       elseif Smoothed(1,i) > Mean + .25*Std
           OnOffIndex(k,i)=1;
       else
            OnOffIndex(k,i)=0;
       end
    end
    % plotting the smoothed data and OnOffIndex for each posture
    figure('name',Names{k})
    plot(t,Smoothed,t,OnOffIndex(k,:)/500)
end

% plotting all 8 channels of emg data for each posture with its
% corresponding OnOffIndex
% for i = 1:11
%     figure('name',Names{i})
%     for j =1:8
%         subplot(4,2,j)
%         plot(t,EMGfilt.(Names{i}).posture.Data(j,:),t,OnOffIndex(i,:)/10000)
%     end
% end



EMG_on=EMGfilt;
%Making a copy of the filtered data to remove the "off" sections for
%feature extraction


% Removing the sections of data we don't want
for j=2:11
    for i=1:8
    EMG_on.(Names{j}).posture.Data(i,:)=EMG_on.(Names{j}).posture.Data(i,:).*OnOffIndex(j,:);
    % Multiplying everything in the "off" sections by zero
    end
    EMG_on.(Names{j}).posture.Data(:,all(EMG_on.(Names{j}).posture.Data == 0)) = [];
    % Removing columns of data that are all zeros
end


% Plotting all 8 emg signals for each posture after removing off sections
% and spikes
% for j=1:11
%     figure('name',Names{j})
%     for i=1:8
%     subplot(4,2,i)
%     plot(EMG_on.(Names{j}).posture.Data(i,:))
%     title(['EMG', num2str(i)])
%     ylim ([-.002,.002])
%     hold on
%     end
% end

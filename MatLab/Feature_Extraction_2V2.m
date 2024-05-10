clear
close all
clc

% Binning data and extracting MAV with various bin sizes

load("EMG_Off_Removed.mat")
Names = fieldnames(EMG_on);

% Sampling frequency and time vector
Fs = EMG_on.HC.posture.Fs(1);

ylabels = {'50ms', '100 ms', '150ms', '200ms'};

l = length(EMG_on.(Names{10}).posture.Data(1,:)); %length of data for HC

for k = 1:4
    binms = k*50; % binsize in ms
    window = round(binms*Fs/1000); % active window in #of data points
    NumBins = ceil(l/window); % 69852/window rounded up to the next whole number
    
    bin = zeros([8,window]);
    Features = zeros([32,NumBins]);
    t = linspace(0,(l/Fs),NumBins);

    for j = 1:NumBins
            first = 1+j*window-window;
            last=j*window;
        
            if last > l
                last=l;
            end
        
            for i = 1:(last-first)
                bin(:,i) = EMG_on.HC.posture.Data(:,first+i);
            end
            
            for i = 1:8
                % Mean Absolute Value
                MABS = mean(abs(bin),2);
                mean(abs(bin),2);
                Features((4*i)-3,j)=MABS(i,1);
            end
    end
    figure(1)
    subplot(4,1,k)
    plot(t,Features(1,:))
    lgd = legend(ylabels{k})
    lgd.Location = 'best'
    xlabel('time [s]')
    ylabel('[V]')
    sgtitle('Mean Absolute Value, HC Posture')
end


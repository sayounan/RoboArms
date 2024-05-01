clear
close all

% Binning data and extracting MAV with various bin sizes

load("EMG_Off_Removed.mat")

% Sampling frequency and time vector
Fs = EMG_on.HC.posture.Fs(1);

ylabels = {'50ms', '100 ms', '150ms', '200ms'};

for k = 1:4
    binms = k*50; % binsize in ms
    window = round(binms*Fs/1000); % active window in #of data points
    NumBins = ceil(69852/window); % 69852/window rounded up to the next whole number
    
    bin = zeros([8,window]);
    Features = zeros([32,NumBins]);
    t = linspace(0,(69852/Fs),NumBins);

    for j = 1:NumBins
            first = 1+j*window-window;
            last=j*window;
        
            if last > 69852
                last=69852;
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
    ylabel(ylabels{k})
end
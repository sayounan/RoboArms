clear
close all

% Binning data and extracting all four features at with 50ms bins

load("EMG_Off_Removed.mat")

Names = fieldnames(EMG_on);
EMG_Feat = EMG_on;

% Sampling frequency and time vector
Fs = EMG_on.HC.posture.Fs(1);
t = linspace(0,(69852/Fs),69852);

binms = 50; % binsize in ms
window = round(binms*Fs/1000); % active window in #of data points
NumBins = ceil(69852/window); % 69852/window rounded up to the next whole number

bin = zeros([8,window]);
Features = zeros([32,315]);

for p = 1:11
    EMG_Feat.(Names{p}).posture.Time = linspace(0,(69852/Fs),NumBins);
    for j = 1:NumBins
        first = 1+j*window-window;
        last=j*window;
    
        if last > 69852
            last=69852;
        end
    
        for i = 1:(last-first)
            bin(:,i) = EMG_on.(Names{p}).posture.Data(:,first+i);
        end
        
        for i = 1:8
            % Mean Absolute Value
            MABS = mean(abs(bin),2);
            mean(abs(bin),2);
            Features((4*i)-3,j)=MABS(i,1);
    
            % Zero Crossings
            ZC = 0;
            for pts=1:window-1
                if bin(i,pts)*bin(i,pts+1) < 0
                    ZC = ZC+1;
                end
            end
            Features((4*i)-2,j) = ZC;
        
            % Slope Sign Changes
            SSC = 0;
            for pts=1:window-2
                S1 = bin(i,pts+1)-bin(i,pts);
                    S2 = bin(i,pts+2)-bin(i,pts+1);
                if S1*S2 < 0
                    SSC = SSC +1;
                end
            end
            Features((4*i)-1,j) = SSC;
    
            % Waveform Length
            WL = 0;
            for pts=1:window-1
                dist = bin(i,pts)-bin(i,pts+1);
                WL = WL + dist;
            end
            Features((4*i),j) = WL;

            % The features are organized thusly with N bins forming a 32XN
            % matrix
            % MAV1.1 MAV1.2 ... MAV1.N
            % ZC1.1  ZC1.2      ZC1.N
            % SSC1.1 SSC1.2     SSC1.N
            % WL1.1  WL1.2      WL1.N
            % MAV2.1 MAV2.2     MAV2.N
            % ZC2.1  ZC2.2      ZC2.N
            % ...
            % SSC8.1 SSC8.2     SSC8.N
            % WL8.1  WL8.2  ... WL8.N
        end
    end
    EMG_Feat.(Names{p}).posture.Data = Features;
end

figure(1)
subplot(4,1,1)
plot(EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(1,:),EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(25,:))
title('Features Hand Closed: EMG1 & EMG 7');
ylabel('MAV')
subplot(4,1,2)
plot(EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(2,:),EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(26,:))
ylabel('ZC')
subplot(4,1,3)
plot(EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(3,:),EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(27,:))
ylabel('SSC')
subplot(4,1,4)
plot(EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(4,:),EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(28,:))
ylabel('WL')


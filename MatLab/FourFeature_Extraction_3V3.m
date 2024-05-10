clear
close all
clc

% Binning data and extracting all four features at with 50ms bins

load("EMG_Off_Removed_Val.mat")

Names = fieldnames(EMG_on);
EMG_Feat = EMG_on;
Fs = EMG_on.HC.posture.Fs(1);
binms = 50; % binsize in ms
window = round(binms*Fs/1000); % active window in #of data points
bin = zeros([8,window]);


for p = 1:11
    % Sampling frequency and time vector
    Size = size(EMG_on.(Names{p}).posture.Data,2);

    t = linspace(0,(Size/Fs),Size);
    NumBins = ceil(Size/window); % data size/active window rounded up to the next whole number
    EMG_Feat.(Names{p}).posture.Time = linspace(0,(Size/Fs),NumBins);
    %Creating a new time vector and saving it to the binned data structure
    Features = zeros([32,NumBins]);

    for j = 1:NumBins
        first = 1+j*window-window;
        last=j*window;
    
        if last > Size
            last=Size;
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
                    if abs(bin(i,pts)-bin(i,pts+1)) > 1e-5 %ZC threshold
                    ZC = ZC+1;
                    end
                end
            end
            Features((4*i)-2,j) = ZC;
        
            % Slope Sign Changes
            SSC = 0;
            for pts=1:window-2
                S1 = bin(i,pts+1)-bin(i,pts);
                S2 = bin(i,pts+2)-bin(i,pts+1);
                if S1*S2 < 0
                    if abs(S2-S1) > 1e-6 %SSC threshold
                    SSC = SSC +1;
                    end
                end
            end
            Features((4*i)-1,j) = SSC;
    
            % Waveform Length
            WL = 0;
            for pts=1:window-1
                %dist = bin(i,pts)-bin(i,pts+1); %dV only
                dist = ((bin(i,pts)-bin(i,pts+1))^2 + (4.5001e-04)^2)^0.5; %Euclidean distance
                WL = WL + abs(dist);
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


%% RD re-thresh
for p = 4
    % Sampling frequency and time vector
    Size = size(EMG_on.(Names{p}).posture.Data,2);

    t = linspace(0,(Size/Fs),Size);
    NumBins = ceil(Size/window); % data size/active window rounded up to the next whole number
    EMG_Feat.(Names{p}).posture.Time = linspace(0,(Size/Fs),NumBins);
    %Creating a new time vector and saving it to the binned data structure
    Features4 = zeros([32,NumBins]);

    for j = 1:NumBins
        first = 1+j*window-window;
        last=j*window;
    
        if last > Size
            last=Size;
        end
    
        for i = 1:(last-first)
            bin(:,i) = EMG_on.(Names{p}).posture.Data(:,first+i);
        end
       
        for i = 1:8

            % Mean Absolute Value
            MABS = mean(abs(bin),2);
            mean(abs(bin),2);
            Features4((4*i)-3,j)=MABS(i,1);

            % Zero Crossings
            ZC = 0;
            for pts=1:window-1
                if bin(i,pts)*bin(i,pts+1) < 0
                    if abs(bin(i,pts)-bin(i,pts+1)) > 1e-9 %ZC no threshold
                    ZC = ZC+1;
                    end
                end
            end
            Features4((4*i)-2,j) = ZC;
        
            % Slope Sign Changes
            SSC = 0;
            for pts=1:window-2
                S1 = bin(i,pts+1)-bin(i,pts);
                S2 = bin(i,pts+2)-bin(i,pts+1);
                if S1*S2 < 0
                    if abs(S2-S1) > 1e-9 %SSC threshold
                    SSC = SSC +1;
                    end
                end
            end
            Features4((4*i)-1,j) = SSC;

            % Waveform Length
            WL = 0;
            for pts=1:window-1
                %dist = bin(i,pts)-bin(i,pts+1); %dV only
                dist = ((bin(i,pts)-bin(i,pts+1))^2 + (4.5001e-04)^2)^0.5; %Euclidean distance
                WL = WL + abs(dist);
            end
            Features4((4*i),j) = WL;

        end
    end
    EMG_Feat.(Names{4}).posture.Data = Features4;
end

%% TAB re-thresh
for p = 9
    % Sampling frequency and time vector
    Size = size(EMG_on.(Names{p}).posture.Data,2);

    t = linspace(0,(Size/Fs),Size);
    NumBins = ceil(Size/window); % data size/active window rounded up to the next whole number
    EMG_Feat.(Names{p}).posture.Time = linspace(0,(Size/Fs),NumBins);
    %Creating a new time vector and saving it to the binned data structure
    Features9 = zeros([32,NumBins]);

    for j = 1:NumBins
        first = 1+j*window-window;
        last=j*window;
    
        if last > Size
            last=Size;
        end
    
        for i = 1:(last-first)
            bin(:,i) = EMG_on.(Names{p}).posture.Data(:,first+i);
        end
       
        for i = 1:8

            % Mean Absolute Value
            MABS = mean(abs(bin),2);
            mean(abs(bin),2);
            Features9((4*i)-3,j)=MABS(i,1);

            % Zero Crossings
            ZC = 0;
            for pts=1:window-1
                if bin(i,pts)*bin(i,pts+1) < 0
                    %if abs(bin(i,pts)-bin(i,pts+1)) > 1e-9 %ZC no threshold
                    ZC = ZC+1;
                    %end
                end
            end
            Features9((4*i)-2,j) = ZC;
        
            % Slope Sign Changes
            SSC = 0;
            for pts=1:window-2
                S1 = bin(i,pts+1)-bin(i,pts);
                S2 = bin(i,pts+2)-bin(i,pts+1);
                if S1*S2 < 0
                    if abs(S2-S1) > 1e-9 %SSC threshold
                    SSC = SSC +1;
                    end
                end
            end
            Features9((4*i)-1,j) = SSC;

            % Waveform Length
            WL = 0;
            for pts=1:window-1
                %dist = bin(i,pts)-bin(i,pts+1); %dV only
                dist = ((bin(i,pts)-bin(i,pts+1))^2 + (4.5001e-04)^2)^0.5; %Euclidean distance
                WL = WL + abs(dist);
            end
            Features9((4*i),j) = WL;

        end
    end
    EMG_Feat.(Names{9}).posture.Data = Features9;
end

%%

figure(1)
subplot(4,1,1)
plot(EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(1,:),EMG_Feat.HC.posture.Time,EMG_Feat.HC.posture.Data(25,:))
title('Features Hand Closed,50 ms bins: EMG1 & EMG 7');
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
xlabel('time [s]')
%% Another Posture Features plot

figure
subplot(4,1,1)
plot(EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(1,:),EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(25,:))
title((Names{9}), ' Features Hand Closed,50 ms bins: EMG1 & EMG 7');
ylabel('MAV')
subplot(4,1,2)
plot(EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(2,:),EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(26,:))
ylabel('ZC')
subplot(4,1,3)
plot(EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(3,:),EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(27,:))
ylabel('SSC')
subplot(4,1,4)
plot(EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(4,:),EMG_Feat.(Names{9}).posture.Time,EMG_Feat.(Names{9}).posture.Data(28,:))
ylabel('WL')
xlabel('time [s]')

figure
subplot(4,1,1)
plot(EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(1,:),EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(25,:))
title((Names{4}), ' Features Hand Closed,50 ms bins: EMG1 & EMG 7');
ylabel('MAV')
subplot(4,1,2)
plot(EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(2,:),EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(26,:))
ylabel('ZC')
subplot(4,1,3)
plot(EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(3,:),EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(27,:))
ylabel('SSC')
subplot(4,1,4)
plot(EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(4,:),EMG_Feat.(Names{4}).posture.Time,EMG_Feat.(Names{4}).posture.Data(28,:))
ylabel('WL')
xlabel('time [s]')

%%
EMG_Feat_Val = EMG_Feat;  %Run these two lines if extracting from Val data, else comment out
clear EMG_Feat

clear ans bin binms dist EMG_on Features first Fs i j last MABS Names NumBins p pts S1 S2 Size SSC t window WL ZC Features4 Features9

save('EMG_Features_Val.mat')
clear
close all

% Loading the EMG data streams
D1 = load ("off_train.mat");
D2 = load("WF_train.mat");
D3 = load("WE_train.mat");
D4 = load("RD_train.mat");
D5 = load("UD_train.mat");
D6 = load("WP_train.mat");
D7 = load("WS_train.mat");
D8 = load("Tad_train.mat");
D9 = load("Tab_train.mat");
D10 = load("HC_train.mat");
D11 = load("HO_train.mat");

% Moving the data streams into a single structure 
EMG=struct('off',D1,'WF',D2,'WE',D3,'RD',D4,'UD',D5,'WP',D6,'WS',D7, ...
    'TAD',D8,'TAB',D9,'HC',D10,'HO',D11);

% Fieled names index
Names = fieldnames(EMG);

for j=1:11
    figure('name',Names{j})
    for i=1:8
        subplot(4,2,i)
        plot(EMG.(Names{j}).posture.Time,EMG.(Names{j}).posture.Data(i,:))
        ylim([-0.003, 0.003])
        title(['Raw EMG', num2str(i)])
        xlabel('time [s]')
    end
end

% The structure from this was saved as "EMG_Training.mat" So that's what I
% call in the filter finction. I did also run this loading the validation
% data, so I have that as well saved as "EMG_Validation.mat" for when I
% need that
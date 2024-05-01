clear
close all

% Moving towards the actual LDA stuff

load EMG_Features.mat

% Pulling Sampling frequency and time vector, idk if I need them but might
% as well

Fs = EMG_Feat.HC.posture.Fs(1,1);
t = EMG_Feat.HC.posture.Time;

Names = (fieldnames(EMG_Feat));

%EMG_Feat_Mean_Vecs = EMG_Feat;

% Computing class mean vectors

M = zeros(32,1);
for j=1:11
    for i=1:32
    M(i,:) = mean(EMG_Feat.(Names{j}).posture.Data(i,:));
    end
    EMG_Feat.(Names{j}).posture.Class_Mean = M;
end

% Computing the All class mean matrix

AllClassMean = zeros(32,11);
for j=1:11
    AllClassMean(:,j) = EMG_Feat.(Names{j}).posture.Class_Mean;
end

% Because of how the cov function is set up, i need to put in the transpose
% of my matrices

BetweenClassScatter = cov(AllClassMean.');

% All class mean vector

AllClassMeanVector = zeros(11, 1);
WithinClassScatter = zeros(32,32);

for i=1:11
AllClassMeanVector(i,:) = mean(AllClassMean(i,:));

% Class scatter matrix for each posture
ScatterIntermediate = zeros(32,630);

for j=1:630
    ScatterIntermediate(:,j) = EMG_Feat.(Names{i}).posture.Data(:,j)-EMG_Feat.(Names{i}).posture.Class_Mean;
end
EMG_Feat.(Names{i}).posture.ClassScatter = cov(ScatterIntermediate.');

WithinClassScatter = WithinClassScatter + EMG_Feat.(Names{i}).posture.ClassScatter;
end

% Optimization Matrix
Optimization = BetweenClassScatter/WithinClassScatter;

% Getting Eigenvectors/Values from the optimization matrix
[Vectors,Values]= eig(Optimization);
% Vectors' columns are the eigenvectors


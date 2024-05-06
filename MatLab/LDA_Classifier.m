clear
close all

% Moving towards the actual LDA stuff
% Fai, Erin Radcliffe, Sari Younan, & Tim Young

load EMG_Features2.mat

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
    % THIS ISN'T MECESSARY, BUT DOING IT DOESN'T SEEM TO EFFECT PREDICTIONS
    %ScatterIntermediate = zeros(32,590);

    % for j=1:590
    %     ScatterIntermediate(:,j) = EMG_Feat.(Names{i}).posture.Data(:,j)-EMG_Feat.(Names{i}).posture.Class_Mean;
    % end
    %EMG_Feat.(Names{i}).posture.ClassScatter = cov(ScatterIntermediate.');

    EMG_Feat.(Names{i}).posture.ClassScatter = cov(EMG_Feat.(Names{i}).posture.Data.');
    WithinClassScatter = WithinClassScatter + EMG_Feat.(Names{i}).posture.ClassScatter;
end

% Optimization Matrix
Optimization = BetweenClassScatter/WithinClassScatter;

% Getting Eigenvectors/Values from the optimization matrix
[Vectors,Values]= eig(Optimization);
% Vectors' columns are the eigenvectors

% We need the first 10 eigenvectors, because the Eigenvalues drop off
% prodigiously
%  0.6811 + 0.0000i 1
%  0.3931 + 0.0000i 2
%  0.2581 + 0.0000i 3
%  0.2464 + 0.0000i 4
%  0.1957 + 0.0000i 5 
%  0.1316 + 0.0000i 6
%  0.0610 + 0.0000i 7
%  0.0492 + 0.0000i 8
%  0.0314 + 0.0000i 9
%  0.0297 + 0.0000i 10
%  0.0000 + 0.0000i 11
% ...
% -0.0000 + 0.0000i 32

  W = zeros(32,10);

  for i=1: 10
      W(:,i) = Vectors(:, i);
  end

Class_Mean_Transform = AllClassMean.' * W;

load EMG_Features_Validation2.mat

Size = size(EMG_Feat_Val.HO.posture.Data,2);

for i = 1:11
    Validation_Transform = EMG_Feat_Val.(Names{i}).posture.Data.' * W;

    EMG_Feat_Val.(Names{i}).posture.Transformed = Validation_Transform;
end

P = zeros(11, 11);
Pred = zeros(1,11);
I = zeros(Size,11);

for h = 1:11
    A = EMG_Feat_Val.(Names{h}).posture.Transformed;
    distance = zeros(Size,11);

    for k =1:Size
    d2 = zeros(1,11);
        for j = 1:11
            for i = 1:10
                d2(1,j) = d2(1,j) + (Class_Mean_Transform(j,i)-A(k,i))^2;
            end
            distance(k,j) = sqrt(d2(1,j));
        end
    end
    [~,Index] = min(distance,[],2);
    I(:,h) = Index;
    for i = 1:11

        Pred(1,i) = sum(Index == i);
    end
    EMG_Feat_val.(Names{h}).posture.Class_Predictions = Pred;
    P(h,:) = Pred;
end

classes = {'off','WF','WE','RD','UD','WP','WS','TAD','TAB','HC','HO'};
H = heatmap(classes, classes, (P/Size*100));
H.Title = 'LDA Predictions';
H.YLabel = 'Expected Class';
H.XLabel = 'Predicted Class';

Correct_Pred = 0;

for i= 1:11
    Correct_Pred = Correct_Pred + P(i,i);
    % This is for identifying patterns in misclassifications i.e. large
    % sections of misclassification or incorrect matrix sizes 
    %figure(i+1)
    %plot(I(:,i), LineStyle='none', Marker="o")
end

% Reporting important data
AllClassMean
BetweenClassScatter
WithinClassScatter
Optimization
Vectors
Values
W
Classifier_Accuracy = Correct_Pred/(Size*11)*100

clear
close all

% Moving towards the actual LDA stuff
% Fai, Erin Radcliffe, Sari Younan, & Tim Young

load EMG_Features7.mat

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
invWCS = inv(WithinClassScatter);
Optimization = BetweenClassScatter*invWCS;

% Getting Eigenvectors/Values from the optimization matrix
[Vectors,Values]= eig(Optimization);
% Vectors' columns are the eigenvectors

% We need the first 10 eigenvectors, because the Eigenvalues drop off
% prodigously

% Sorting the eigenvectors and values, this really doesn't seem to be
% necessary, if you want to use it comment out lines 88-90, and uncomment
% lines 69-71 and 92-94
% [v,ind] = sort(diag(Values));
% Values_Sort = Values(ind,ind);
% Vectors_sort = Vectors(:,ind);

  W = zeros(32,10);

  for i=1: 10
      W(:,i) = Vectors(:, i);
  end

  % for i=1: 10
  %   W(:,i) = Vectors_sort(:, 33-i);
  % end

Class_Mean_Transform = AllClassMean.' * W;

load EMG_Features_Validation7.mat

for i = 1:11
    Validation_Transform = EMG_Feat_Val.(Names{i}).posture.Data.' * W;

    EMG_Feat_Val.(Names{i}).posture.Transformed = Validation_Transform;
end

P = zeros(11, 11);
Pred = zeros(1,11);


for h = 1:11
    Size = size(EMG_Feat_Val.(Names{h}).posture.Data,2);
    A = EMG_Feat_Val.(Names{h}).posture.Transformed;
    distance = zeros(Size,11);
    I = zeros(Size,11);

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
    EMG_Feat_Val.(Names{h}).posture.Class_Predictions = Pred;
    P(h,:) = Pred/Size*100; 

    % This is for identifying patterns in misclassifications i.e. large
    % sections of misclassification or incorrect matrix sizes 
    figure(h)
    plot(I(:,h), LineStyle='none', Marker="o")
end

figure(12)
classes = {'off','WF','WE','RD','UD','WP','WS','TAD','TAB','HC','HO'};
H = heatmap(classes, classes, P);
H.Title = 'LDA Predictions';
H.YLabel = 'Expected Class';
H.XLabel = 'Predicted Class';

Correct_Pred = 0;

for i= 1:11
    Correct_Pred = Correct_Pred + P(i,i);
end

% Reporting important data
AllClassMean;
BetweenClassScatter;
WithinClassScatter;
Optimization;
Vectors;
Values;
W;
Classifier_Accuracy = Correct_Pred/(sum(sum(P)))*100

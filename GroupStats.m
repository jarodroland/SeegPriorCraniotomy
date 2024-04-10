flags.saveFigures = false; %true; %
flags.closeFigures = true; %false; %
addpath(genpath('./Violinplot-Matlab'))

projDir = GetProjectPath();
outDir = fullfile(projDir, 'figures');

% define which group of data to use
dataGroup = 'SLCH';

maxNumContacts = 16;        %MAGICNUMBER: longest electrode has 16 contacts

% load subject list
data = load(['subjectListPriorCran' dataGroup '.mat']);
subjListNaiveBone = data.subjListNaiveBone;
subjListPriorCran = data.subjListPriorCran;

%% Native Bone Subjects
numSubj = length(subjListNaiveBone);
naiveFirstContactErrors = [];
naiveLastContactErrors = [];
naiveAllContactErrors = [];
allContactErrors = [];
for subjIdx = 1:numSubj
    try
        [electrodeData, axes3DFig, axes2DFig] = CalcSEEGTrajectory(subjListNaiveBone{subjIdx}, 'dataGroup', dataGroup, 'isFigVisible', ~flags.closeFigures);  % if not saving figures, then don't show them

        % save figures
        if(flags.saveFigures)
            outFileRoot = fullfile(outDir, [subjListNaiveBone{subjIdx} '_ErrorContactsToBoltTrajectory']);
            saveas(axes3DFig.Parent, [outFileRoot '-3dCoords.png']);
            savefig(axes3DFig.Parent, [outFileRoot '-3dCoords.fig']);
            saveas(axes2DFig.Parent, [outFileRoot '-2dErrors.png']);
            savefig(axes2DFig.Parent, [outFileRoot '-2dErrors.fig']);
        end %if(flags.saveFigures)
    
        % close figures
        if(flags.closeFigures)
            close(axes3DFig.Parent);
            close(axes2DFig.Parent);
        end

        % collect first and last contact error measurements
        for elecIdx = 1:length(electrodeData)
            naiveFirstContactErrors = [naiveFirstContactErrors, electrodeData(elecIdx).contactsError(1)];       %MAGICNUMBER: 1st contact is the furthest from the bolt
            naiveLastContactErrors = [naiveLastContactErrors, electrodeData(elecIdx).contactsError(end)];       %MAGICNUMBER: last contact is the furthest from the bolt
            naiveAllContactErrors(end+1, length(electrodeData(elecIdx).contactsError):-1:1) = electrodeData(elecIdx).contactsError;
            allContactErrors(end+1, length(electrodeData(elecIdx).contactsError):-1:1) = electrodeData(elecIdx).contactsError;
        end
    catch exception
        disp(['Error on subject ' subjListNaiveBone{subjIdx}]);
        disp(getReport(exception));
    end
end

%% Prior Craniotomy Subjects
numSubj = length(subjListPriorCran);
cranFirstContactErrors = [];
cranLastContactErrors = [];
cranAllContactErrors = [];
for subjIdx = 1:numSubj
    try
        [electrodeData, axes3DFig, axes2DFig] = CalcSEEGTrajectory(subjListPriorCran{subjIdx}, 'dataGroup', dataGroup, 'isFigVisible', ~flags.closeFigures);  % if not saving figures, then don't show them

        % save figures
        if(flags.saveFigures)
            outFileRoot = fullfile(outDir, [subjListPriorCran{subjIdx} '_ErrorContactsToBoltTrajectory']);
            saveas(axes3DFig.Parent, [outFileRoot '-3dCoords.png']);
            savefig(axes3DFig.Parent, [outFileRoot '-3dCoords.fig']);
            saveas(axes2DFig.Parent, [outFileRoot '-2dErrors.png']);
            savefig(axes2DFig.Parent, [outFileRoot '-2dErrors.fig']);
        end %if(flags.saveFigures)
    
        % close figures
        if(flags.closeFigures)
            close(axes3DFig.Parent);
            close(axes2DFig.Parent);
        end

        % collect first and last contact error measurements
        for elecIdx = 1:length(electrodeData)
%             cranFirstContactErrors = arrayfun(@(x) [cranFirstContactErrors; x.contactsError(1)], electrodeData, 'UniformOutput', true);
            cranFirstContactErrors = [cranFirstContactErrors, electrodeData(elecIdx).contactsError(1)];     %MAGICNUMBER: 1st contact is the furthest from the bolt
            cranLastContactErrors = [cranLastContactErrors, electrodeData(elecIdx).contactsError(end)];     %MAGICNUMBER: last contact is the closet to the bolt
            cranAllContactErrors(end+1, length(electrodeData(elecIdx).contactsError):-1:1) = electrodeData(elecIdx).contactsError;
            allContactErrors(end+1, length(electrodeData(elecIdx).contactsError):-1:1) = electrodeData(elecIdx).contactsError;
        end
    catch exception
        disp(['Error on subject ' subjListPriorCran{subjIdx}]);
        disp(getReport(exception));
    end
end

%% Violin plots
%% plot distribution of errors for contact position for both groups lumped together
allContactErrors(allContactErrors==0) = NaN;
cMapContactPosition = hsv(maxNumContacts);
figViolinContactPosition = figure();
positionLabels = cellstr(num2str((16:-1:1)'));
% violinContactErr = violinplot(fliplr(allContactErrors), flip(positionLabels));
violinContactErr = violinplot(allContactErrors, flip(positionLabels));
xlim([0.5, maxNumContacts + 0.5]);
for idx = 1:maxNumContacts
    subViolinPlot = violinContactErr(idx);
    BeautifyViolinPlot(subViolinPlot)
    subViolinPlot.BoxWidth = 0.05;
    subViolinPlot.ScatterPlot.MarkerFaceColor = cMapContactPosition(idx, :);
    %subViolinPlot.BoxPlot.EdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.ScatterPlot.MarkerEdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.MeanPlot.Color = cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.EdgeColor = [0, 0, 0]; %cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.LineWidth = 1;
    subViolinPlot.ScatterPlot.LineWidth = 1.0;
    subViolinPlot.ScatterPlot.MarkerFaceAlpha = 0.9;
    %subViolinPlot.ScatterPlot.SizeData = 10;
    subViolinPlot.ViolinColor = {cMapContactPosition(idx,:)};
end
title('Error Distribution by Contact Position');
ylabel('Error (mm)');
xlabel('Contact Position (reverse order)');

%% plot both groups as Left and Right halves of violin plot
violinGap = 0.07;

% Naive distribution of errors for contact position
naiveAllContactErrors(naiveAllContactErrors==0) = NaN;
cMapContactPosition = hsv(maxNumContacts);
figNaiveAndCranViolinContactPosition = figure();
positionLabels = cellstr(num2str((16:-1:1)'));
% violinContactErr = violinplot(fliplr(allContactErrors), flip(positionLabels));
violinNaiveContactErr = violinplot(naiveAllContactErrors, flip(positionLabels), 'QuartileStyle', 'shadow', 'HalfViolin', 'left', 'ShowWhiskers', false);
xlim([0.5, maxNumContacts + 0.5]);

figNaiveAndCranViolinContactPosition.Position = [0, 0, 1400, 400];
figNaiveAndCranViolinContactPosition.CurrentAxes.FontSize = 20;

colorNative = [0.3, 0.3, 1.0];
colorPrior = [1.0, 0.0, 0.0];
colorNativeMarker = [0.7, 0.7, 1.0];
colorPriorMarker = [1.0, 0.7, 0.7];

for idx = 1:maxNumContacts
    subViolinPlot = violinNaiveContactErr(idx);
    BeautifyViolinPlot(subViolinPlot)
    % subViolinPlot.ViolinColor = {colorNative}; %{cMapContactPosition(idx,:)};
    subViolinPlot.BoxWidth = 0.05;
    subViolinPlot.ScatterPlot.MarkerFaceColor = colorNativeMarker; %cMapContactPosition(idx, :);
    subViolinPlot.ViolinColor = {colorNative}; %{cMapContactPosition(idx,:)};   % NOTE: Must set ViolinColor AFGER ScatterPlot.MarkerFaceFolor
    %subViolinPlot.BoxPlot.EdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.ScatterPlot.MarkerEdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.MeanPlot.Color = cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.EdgeColor = [0, 0, 0]; %cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.LineWidth = 0.5;
    subViolinPlot.ScatterPlot.LineWidth = 0.5;
    subViolinPlot.ScatterPlot.MarkerEdgeAlpha = 0.5;
    subViolinPlot.ScatterPlot.XData = subViolinPlot.ScatterPlot.XData - violinGap;
    %subViolinPlot.ScatterPlot.SizeData = 10;
    subViolinPlot.ViolinPlot.XData = subViolinPlot.ViolinPlot.XData - violinGap;
    subViolinPlot.MedianPlot.MarkerEdgeColor = [0 0 0];
    subViolinPlot.MedianPlot.XData = subViolinPlot.MedianPlot.XData - violinGap;
    subViolinPlot.ViolinPlotQ.XData = subViolinPlot.ViolinPlotQ.XData - violinGap;
    subViolinPlot.ViolinPlotQ.LineWidth = 2.0;
    subViolinPlot.ViolinPlotQ.EdgeColor = [0 0 0];
    subViolinPlot.ViolinPlotQ.FaceAlpha = 0.4;
end

% Cran distribution of errors for contact position
cranAllContactErrors(cranAllContactErrors==0) = NaN;
cMapContactPosition = hsv(maxNumContacts);
% figViolinContactPosition = figure();
positionLabels = cellstr(num2str((16:-1:1)'));
% violinContactErr = violinplot(fliplr(allContactErrors), flip(positionLabels));
violinCranContactErr = violinplot(cranAllContactErrors, flip(positionLabels), 'QuartileStyle', 'shadow', 'HalfViolin', 'right', 'ShowWhiskers', false);
xlim([0.5, maxNumContacts + 0.5]);
for idx = 1:maxNumContacts
    subViolinPlot = violinCranContactErr(idx);
    BeautifyViolinPlot(subViolinPlot)
    subViolinPlot.BoxWidth = 0.05;
    subViolinPlot.ScatterPlot.MarkerFaceColor = colorPriorMarker; %cMapContactPosition(idx, :);
    subViolinPlot.ViolinColor = {colorPrior}; %{cMapContactPosition(idx,:)};
    %subViolinPlot.BoxPlot.EdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.ScatterPlot.MarkerEdgeColor = cMapContactPosition(idx, :);
    %subViolinPlot.MeanPlot.Color = cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.EdgeColor = [0, 0, 0]; %cMapContactPosition(idx, :);
    subViolinPlot.ViolinPlot.LineWidth = 0.5;
    subViolinPlot.ScatterPlot.LineWidth = 1.0;
    subViolinPlot.ScatterPlot.MarkerEdgeAlpha = 0.5;
    subViolinPlot.ScatterPlot.XData = subViolinPlot.ScatterPlot.XData + violinGap;
    %subViolinPlot.ScatterPlot.SizeData = 10;
    subViolinPlot.ViolinPlot.XData = subViolinPlot.ViolinPlot.XData + violinGap;
    subViolinPlot.MedianPlot.XData = subViolinPlot.MedianPlot.XData + violinGap;
    subViolinPlot.MedianPlot.MarkerEdgeColor = [0 0 0];
    subViolinPlot.ViolinPlotQ.XData = subViolinPlot.ViolinPlotQ.XData + violinGap;
    subViolinPlot.ViolinPlotQ.LineWidth = 2.0;
    subViolinPlot.ViolinPlotQ.EdgeColor = [0 0 0];
    subViolinPlot.ViolinPlotQ.FaceAlpha = 0.4;
end
title('Native (left) and Prior Craniotomy (right) Error Distribution by Contact Position');
ylabel('Error (mm)');
xlabel('Contact Position (reverse order)');

%%
% plot distribution of first and last contacts as grouped by Native and Prior
grpLabels = vertcat( repmat({'Native Bone'}, length(naiveFirstContactErrors), 1), repmat({'Prior Cran'}, length(cranFirstContactErrors), 1) );

figFirstErr = figure();
pltViolinFirstErr = violinplot([naiveFirstContactErrors, cranFirstContactErrors]', grpLabels);
BeautifyViolinPlot(pltViolinFirstErr(1));
BeautifyViolinPlot(pltViolinFirstErr(2));
pltViolinFirstErr(1).ScatterPlot.MarkerFaceColor = colorNativeMarker;
pltViolinFirstErr(1).ViolinPlot.FaceColor = colorNative;
pltViolinFirstErr(2).ScatterPlot.MarkerFaceColor = colorPriorMarker;
pltViolinFirstErr(2).ViolinPlot.FaceColor = colorPrior;
pltViolinFirstErr(1).ScatterPlot.LineWidth = 1.0;
pltViolinFirstErr(2).ScatterPlot.LineWidth = 1.0;
pltViolinFirstErr(1).ScatterPlot.MarkerFaceAlpha = 0.5;
pltViolinFirstErr(2).ScatterPlot.MarkerFaceAlpha = 0.5;
pltViolinFirstErr(1).ViolinColor = {colorNative};
pltViolinFirstErr(2).ViolinColor = {colorPrior};
ylim([0, 20]);
ylabel('Contact-to-Trajectory Error (mm)');
title('Error at First Contact');
figFirstErr.CurrentAxes.FontSize = 20;
% figFirstErr.CurrentAxes.FontName = 'Hack Nerd Font';

figLastErr = figure();
pltViolinLastErr = violinplot([naiveLastContactErrors, cranLastContactErrors]', grpLabels);
BeautifyViolinPlot(pltViolinLastErr(1));
BeautifyViolinPlot(pltViolinLastErr(2));
pltViolinLastErr(1).ScatterPlot.MarkerFaceColor = colorNativeMarker;
pltViolinLastErr(2).ScatterPlot.MarkerFaceColor = colorPriorMarker;
pltViolinLastErr(1).ViolinPlot.FaceColor = colorNative;
pltViolinLastErr(2).ViolinPlot.FaceColor = colorPrior;
pltViolinLastErr(1).ScatterPlot.LineWidth = 1.0;
pltViolinLastErr(2).ScatterPlot.LineWidth = 1.0;
pltViolinLastErr(1).ScatterPlot.MarkerFaceAlpha = 0.5;
pltViolinLastErr(2).ScatterPlot.MarkerFaceAlpha = 0.5;
pltViolinLastErr(1).ViolinColor = {colorNative};
pltViolinLastErr(2).ViolinColor = {colorPrior};

ylim([0, 20]);
ylabel('Contact-to-Trajectory Error (mm)');
title('Error at Last Contact');
figLastErr.CurrentAxes.FontSize = 20;
% figLastErr.CurrentAxes.FontName = 'Hack Nerd Font';

%% stats
errorThreshold = 6.0; %5.5; %3.5; %6.5; %

numNaive = length(naiveFirstContactErrors);
numCran = length(cranFirstContactErrors);

fprintf('Using an error threshold = %1.1f\n', errorThreshold);

% first contact error stats
% chiSquareTable = table([1;2],[3;4],'VariableNames',{'Experimental','Control'},'RowNames',{'Disease','NoDisease'})
outlierNaive = length(naiveFirstContactErrors(naiveFirstContactErrors >= errorThreshold));
outlierCran = length(cranFirstContactErrors(cranFirstContactErrors >= errorThreshold));
contTable = [outlierNaive, numNaive-outlierNaive; outlierCran, numCran-outlierCran];
[pValue, chiSquareValue] = chisquarecont(contTable);
[ksHypothesisFirst, ksPvalueFirst, ks2statFirst] = kstest2(naiveFirstContactErrors, cranFirstContactErrors);
fprintf('First electrode errors: %d native bone, %d prior cran, %d native outliers, %d cran outliers\n', numNaive, numCran, outlierNaive, outlierCran);
fprintf('Chi-Square test p-value = %1.10f, chi-square value = %1.3f\n', pValue, chiSquareValue);
fprintf('K-S test p-value=%1.5f, K-S test statistic = %1.3f, are distributions non-equal = %d\n', ksPvalueFirst, ks2statFirst, ksHypothesisFirst);

fprintf('-----\n');

% first contact error stats
outlierNaive = length(naiveLastContactErrors(naiveLastContactErrors >= errorThreshold));
outlierCran = length(cranLastContactErrors(cranLastContactErrors >= errorThreshold));
contTable = [outlierNaive, numNaive-outlierNaive; outlierCran, numCran-outlierCran];
[pValue, chiSquareValue] = chisquarecont(contTable);
[ksHypothesisLast, ksPvalueLast, ks2statLast] = kstest2(naiveLastContactErrors, cranLastContactErrors);
fprintf('Last electrode errors: %d native bone, %d prior cran, %d native outliers, %d cran outliers\n', numNaive, numCran, outlierNaive, outlierCran);
fprintf('Chi-Square test p-value = %1.3f, chi-square value = %1.3f\n', pValue, chiSquareValue);
fprintf('K-S test p-value=%1.3f, K-S test statistic = %1.3f, are distributions non-equal = %d\n', ksPvalueLast, ks2statLast, ksHypothesisLast);

fprintf('-----\n\n\n');

%% Plot smoothing spline that is used to find peak difference in outliers at varying thresholds

listThresholds = 2.0:0.1:10.0;
numThresholds = length(listThresholds);
numOverThresholdNaive = zeros(0, numThresholds);
numOverThresholdCran = zeros(0, numThresholds);
for idx = 1:numThresholds
    numOverThresholdNaive(idx) = length(naiveFirstContactErrors(naiveFirstContactErrors >= listThresholds(idx)));
    numOverThresholdCran(idx) = length(cranFirstContactErrors(cranFirstContactErrors >= listThresholds(idx)));
end
deltaOutliers = (numOverThresholdCran - numOverThresholdNaive);

% code generated by Matlab Curve Fitter App
[xData, yData] = prepareCurveData( listThresholds, deltaOutliers );
ft = fittype( 'smoothingspline' );  % Set up fittype and options.
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.9;

[fitresult, gof] = fit( xData, yData, ft, opts ); % Fit model to data.

figure( 'Name', 'SpineFitOutliers' ); % Plot fit with data.
h = plot( fitresult, xData, yData);
legend( h, 'Number Over Threshold', 'Spline Model', 'Location', 'NorthEast', 'Interpreter', 'none' );
xlabel( 'Thresholds (mm)', 'Interpreter', 'none' );
ylabel( 'Differnce in Number of Outliers Between Cran and Naive', 'Interpreter', 'none' );
grid on


function BeautifyViolinPlot(pltViolin)

    % beautify violin plots
%     color = [0.1, 0.1, 0.1];
%     pltViolin.BoxPlot.EdgeColor = color;
%     width = 0.8;
%     pltViolin.BoxPlot.LineWidth = width;
%     width = 0.02;%0.0175;
%     pltViolin.BoxWidth = width;
%     alpha = 0.5;
%     pltViolin.BoxPlot.FaceAlpha = alpha;
    marker = 'o';%'.';
    pltViolin.ScatterPlot.Marker = marker;
    color = [0, 0, 0];
    pltViolin.ScatterPlot.MarkerEdgeColor = color;
    alpha = 0.5;
    pltViolin.ScatterPlot.MarkerEdgeAlpha = alpha;
    width = 2.0;% 1.5;% 
    pltViolin.ScatterPlot.LineWidth = width;
    alpha = 0.5;
    pltViolin.ScatterPlot.MarkerFaceAlpha = alpha;
    color = [0, 0, 0];
    pltViolin.MeanPlot.Color = color;
    width = 1.5;
    pltViolin.MeanPlot.LineWidth = width;
    color = [0, 0, 0];
    pltViolin.ViolinPlot.EdgeColor = color;
    width = 1.5;
    pltViolin.ViolinPlot.LineWidth = width;

end %function BeutifyViolinPlot

function [electrodeData, axes3DSimple, axes3DModel, axes2DErrors] = CalcSEEGTrajectory(subjID, varargin)
% [electrodeData, axes3DSimple, axes2DFig] = CalcSEEGTrajectory(subjID, varargin) 
% 
% Description:
%   Estimates a trajectory from points on the SEEG bolt, then calculates the error from the trajectory for each electrode contact
%   The data for each electrode and a figure handle is returned.
%
% Usage:
%   >> [electrodeData, axesPlot] = CalcSEEGTrajectory('SUBJECT001');
%
% Output:
%   electrodeData - structure with member variables label, contactsCoords, contactsError, and boltPts
%   axes3DSimple - axes object for the plot of electrode with a simple circle marker indicating the 3D contact points
%   axes3DFig - axes handle of the 2D contact error plots
%   
% Required Parameters:
%   subjID - The subject ID (with or without a 'sub-' prefix)
%
% Optional Parameters:
%   dataGroup - The group by which subject data is organized into seperate bids directories under the data/bids/ directory
%   useAxes - (axes hanlde) pass an axes handle on which to plot
%   isFigVisible - (bool) Indicates if figues should be visible (default) or hidden
%   isSaveFigures - (bool) Indicates if figues should be saved to the default output directory
%   
% Author:
%   Jarod L Roland
%   Department of Neurosurgery
%   Washington University in St. Louis
%
params = inputParser;
addRequired(params, 'subjID', @(x) true);
addParameter(params, 'dataGroup', 'slch', @ischar);   % slch, ucsf, 'vera-oldAMCandBJH' 
addParameter(params, 'useAxes', 0, @(x) isa(x, 'matlab.graphics.axis.Axes'));
addParameter(params, 'isFigVisible', true, @islogical);
addParameter(params, 'isSaveFigures', false, @islogical);
parse(params, subjID, varargin{:});

%% init vars
dataGroup = lower(params.Results.dataGroup);
axes3DSimple = params.Results.useAxes;
isFigVisible = params.Results.isFigVisible;
flags.saveFigures = params.Results.isSaveFigures;

projDir = GetProjectPath();
bidsDir = fullfile(projDir, 'data', 'bids', dataGroup);
outDir = fullfile(projDir, 'figures');

% append 'sub-' prefix if not already present
if(length(subjID) < 4 || ~strcmp(subjID(1:4), 'sub-'))
    subjID = ['sub-' subjID];
end

% create a figure unless a set of axes is already passed in
if(axes3DSimple == 0)
    figAllContactsWithTrajectory = figure('visible', isFigVisible);
    axes();
    axes3DSimple = figAllContactsWithTrajectory.CurrentAxes;
else
    figAllContactsWithTrajectory = axes3DSimple.Parent;
end
axis equal;
axis vis3d;
title(subjID);
hold on

% create another figure for the 3D model
fig3DModel = figure('DoubleBuffer', 'on');
axis();
axis equal;
axis vis3d
title(subjID);
hold on;
axes3DModel = fig3DModel.CurrentAxes;


%% load electrode coordinates
electrodeDir = fullfile(bidsDir, 'derivatives', 'trajectories', subjID, [subjID '_electrodes-all']);
assert(~isempty(electrodeDir), ['Error: No electrodes directory found for ' subjID]);

electrodeFileList = dir( fullfile(electrodeDir, [subjID '_electrodes-sortedcontacts-*.dat']) );
assert(~isempty(electrodeFileList), ['Error: No electrodes files found for ' subjID]);

% get list of all electrodes files
numElectrodes = length(electrodeFileList);
electrodeData(numElectrodes).contactsCoords = [];
electrodeData(numElectrodes).boltPts = [];
electrodeData(numElectrodes).label = '';
electrodeData(numElectrodes).contactsError = [];

%% loop through all electrode files
for electrodeIdx = 1:numElectrodes
    electrodeFilepath = fullfile(electrodeFileList(electrodeIdx).folder, electrodeFileList(electrodeIdx).name);
    
    % read electrode contacts
    fid = fopen(electrodeFilepath);
%     assert(fid, ['Error: Error opening electrode contacts file (' electrodeFilepath ')']);
    numLines = 1;
    strLine = fgetl(fid);
    while(~strcmp(strtrim(strLine), 'info') && ~feof(fid))
        strLine = fgetl(fid);
        numLines = numLines + 1;
    end
    numContacts = fscanf(fid, 'numpoints %d', 1);
    % assert(isnumeric(numContacts) && (numContacts > 0), ['Error: Invalid numContacts read from electrodes file (' electrodeFilepath ')']);
    if(isnumeric(numContacts) && any(numContacts > 0))
        % SLCH data files have info and numpoints
        fseek(fid, 0, 'bof');   % now go back to begining of file and read all contact coordinates
        contactsCoords = zeros(numContacts, 3);
        for contactIdx = 1:numContacts
            contactsCoords(contactIdx, :) = fscanf(fid, '%f %f %f\n', [1, 3]);
        end
        fclose(fid);
    else
        error('Error reading electrode contacts dat file');
    end

    %NOTE: x coordinates are reversed positive/negative relative to Matlab's coordinate system
    contactsCoords(:, 1) = -contactsCoords(:, 1);    % invert X coordinates for plotting in Matlab

    % read bolt coordinates
    dashIndices = strfind(electrodeFileList(electrodeIdx).name, '-');
    dotIndices = strfind(electrodeFileList(electrodeIdx).name, '.');
    electrodeData(electrodeIdx).label = electrodeFileList(electrodeIdx).name(dashIndices(end)+1:dotIndices(end)-1);
    boltFilename = [subjID '_electrodes-bolt' electrodeFileList(electrodeIdx).name(dashIndices(end):end)];
    boltFilepath = fullfile(electrodeFileList(electrodeIdx).folder, boltFilename);
    assert(isfile(boltFilepath), ['Error: Corresponding bolt points file not found (' boltFilepath ')']);
    fid = fopen(boltFilepath);
%     assert(fid, ['Error: Error opening bolt points file (' boltFilepath ')']);
    strLine = fgetl(fid);
    numLines = 0;
    while(~strcmp(strLine, 'info') && ~feof(fid))
        strLine = fgetl(fid);
        numLines = numLines + 1;
    end
    numBoltPts = fscanf(fid, 'numpoints %d', 1);
    % assert(isnumeric(numBoltPts) && (numBoltPts > 0), ['Error: Invalid numBoltPts read from bolts file (' boltFilepath ')']);
    if(isnumeric(numBoltPts) && any(numBoltPts > 0))
        % SLCH data files have info and numpoints
        fseek(fid, 0, 'bof');   % now go back to begining of file and read all contact coordinates
        boltPts = zeros(numBoltPts, 3);
        for i = 1:numBoltPts
            boltPts(i, :) = fscanf(fid, '%f %f %f\n', [1, 3]);
        end
        fclose(fid);
    else
        error('Error reading bolt coordinates dat file');
    end

    %NOTE: x coordinates are reversed positive/negative relative to Matlab's coordinate system
    boltPts(:, 1) = -boltPts(:, 1);    % invert X coordinates for plotting in Matlab


    %% calculate a straight line through bolt and find distance of each electrode away from line
    % ref: https://www.mathworks.com/matlabcentral/answers/406619-3d-coordinates-line-of-fit
    meanBoltPts = mean(boltPts, 1);
    residBoltPts = bsxfun(@minus, boltPts, meanBoltPts);
    residContactsToBoltPts = bsxfun(@minus, contactsCoords, meanBoltPts);
    covarianceMat = (residBoltPts' * residBoltPts) / (numContacts - 1);
    [svdU, ~] = svd(covarianceMat, 0);
    
    projResidBoltPts = residBoltPts * svdU(:, 1);     % project residual on to svdU(:, 1)
    projResidBoltMin = min(projResidBoltPts);
    projResidBoltMax = max(projResidBoltPts);
    projBoltStart = projResidBoltMin * svdU(:, 1)' + meanBoltPts;
%     projBoltEnd = projResidBoltMax * svdU(:, 1)' + meanBoltPts;
    
    projResidContactsToBolt = residContactsToBoltPts(1, :) * svdU(:, 1);       %MAGICNUMBER residContactsToBoltPts(1, :) is the first contact, which is furthest from the bolt
%     projResidContactsToBolt = residContactsToBoltPts(idxFurthestContact, :) * svdU(:, 1);
    projBoltToLastCoord = projResidContactsToBolt * svdU(:, 1)' + meanBoltPts;
    
    % rename for easy reference
    trajectoryStart = projBoltStart;
    trajectoryEnd = projBoltToLastCoord;

    % for each contact project onto line of trajectory and draw a tangential line from trajectory to each contact representing measured error
    contactsError = zeros(numContacts, 1);
    for i = 1:numContacts
        iContact = contactsCoords(i, :);
        projection = trajectoryStart + ( (trajectoryStart - trajectoryEnd) / norm(trajectoryStart - trajectoryEnd) ) * dot( (iContact - trajectoryStart), (trajectoryStart - trajectoryEnd) ) / norm(trajectoryStart - trajectoryEnd);  % projection of iContact onto trajectory
        plot3(axes3DSimple, [projection(1), iContact(1)], [projection(2), iContact(2)], [projection(3), iContact(3)], 'r-', 'LineWidth', 3);% ); %

        contactsError(i) = norm(iContact - projection);
    end

    % project bolt coordinates onto trajectory for 3D model visualization
    boltPtsProjected = zeros(size(boltPts));
    for boltIdx = 1:numBoltPts
        iBoltPt = boltPts(boltIdx, :);
        boltPtsProjected(boltIdx, :) = trajectoryStart + ( (trajectoryStart - trajectoryEnd) / norm(trajectoryStart - trajectoryEnd) ) * dot( (iBoltPt - trajectoryStart), (trajectoryStart - trajectoryEnd) ) / norm(trajectoryStart - trajectoryEnd);  % projection of iContact onto trajectory
        
        % % TESTING: plot green X on bolt points projected onto trajectory line
        % plot3(axes3DSimple, boltPtsProjected(boltIdx, 1), boltPtsProjected(boltIdx, 2), boltPtsProjected(boltIdx, 3), 'gx', 'LineWidth', 3);
    end

    %% store data in a structure
%     electrodeData(fileIdx).label = electrodeFilepath(strfind(electrodeFilepath, [electrodeContactsSet '-'])+length(electrodeContactsSet)+1:strfind(electrodeFilepath, '.dat')-1);
    electrodeData(electrodeIdx).boltPts = boltPts;
    electrodeData(electrodeIdx).contactsCoords = contactsCoords;
    electrodeData(electrodeIdx).contactsError = contactsError;

    %% Plot contacts and trajectory figure
    title(strrep(subjID, '_', '\_'));
    
    % plot bolt line
    % plot3(axes3DSimple, [projBoltStart(1), projBoltEnd(1)], [projBoltStart(2), projBoltEnd(2)], [projBoltStart(3), projBoltEnd(3)], '-k', 'LineWidth', 3);
    
    % trajectory line
    plot3(axes3DSimple, [trajectoryStart(1), trajectoryEnd(1)], [trajectoryStart(2), trajectoryEnd(2)], [trajectoryStart(3), trajectoryEnd(3)], '-.k', 'LineWidth', 1);
    
    % plot bolt points as simple circles
    for i=1:numBoltPts
        plot3(axes3DSimple, boltPts(i, 1), boltPts(i, 2), boltPts(i, 3), 'ko', 'LineWidth', 2);
    end

    % plot electrode contacts as simple circles
    % plot3(axes3DSimple, contactsCoords(1, 1), contactsCoords(1, 2), contactsCoords(1, 3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);        % red circle for 1st contact
    % plot3(axes3DSimple, contactsCoords(1, 1), contactsCoords(1, 2), contactsCoords(1, 3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);        % blue circle for 1st contact
    for i=1:numContacts 
        plot3(axes3DSimple, contactsCoords(i, 1), contactsCoords(i, 2), contactsCoords(i, 3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);    % blue circle for all subsequent contacts
    end

    %% Plot 3D model of electrode
    [contactsModels, electrodeModel, boltModel, ~] = PlotElectrodeModel(contactsCoords, boltPtsProjected, 'useAxes', axes3DModel);
%     [contactsModels, electrodeModel, boltModel, ~] = PlotElectrodeModel(contactsCoords, boltPtsProjected, 'useAxes', axes3DSimple);
    electrodeData(electrodeIdx).contactsModels = contactsModels;
    electrodeData(electrodeIdx).electrodeModel = electrodeModel;
    electrodeData(electrodeIdx).boltModel = boltModel;

end %for fileIdx = 1:length(fileList)

backgroundColor = [1, 1, 1]; %[0, 0, 0]; %
axes3DSimple.Parent.Color = backgroundColor;
axes3DSimple.Color = backgroundColor;
axis(axes3DSimple, 'tight');
grid(axes3DSimple, 'on');
axes3DSimple.LineWidth = 2;

% give same view angle to both 3D plots
axes3DSimple.View = [45, 45];   % view from 45 degree oblique angle
axes3DModel.View = [45, 45];


%% Plot electrode errors on 2D graph
cMap = hsv(numElectrodes);
figErrors = figure('visible', isFigVisible);
for elecIdx = 1:numElectrodes
    hold on;
    %NOTE: use flipud to reorder from nearest to furthest to plot from Left to Right on the 2D plot
    plot(figErrors.CurrentAxes, flipud(electrodeData(elecIdx).contactsError), 'o:', 'MarkerSize', 10, 'LineWidth', 1, 'Color', cMap(elecIdx, :), 'MarkerEdgeColor', 'black', 'MarkerFaceColor', cMap(elecIdx, :));
end
legend(figErrors.CurrentAxes, electrodeData.label, 'Location', 'Northwest', 'NumColumns', 1)
ylim(figErrors.CurrentAxes, [0, 20]);
title(strrep(subjID, '_', '\_'));
ylabel('Contact to Trajectory Error (mm)');
xlabel('Contact Position');
axes2DErrors = figErrors.CurrentAxes;

%% save figure
if(flags.saveFigures)
    outFileRoot = fullfile(outDir, [subjID '_ErrorContactsToBoltTrajectory']);
    saveas(figAllContactsWithTrajectory, [outFileRoot '.png']);
    savefig(figAllContactsWithTrajectory, [outFileRoot '.fig']);
    close(figAllContactsWithTrajectory);
end %if(flags.saveFigures)


end % function
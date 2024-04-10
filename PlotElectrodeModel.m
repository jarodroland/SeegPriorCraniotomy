function [contactsModels, electrodeModel, boltModel, axes3DModel] = PlotElectrodeModel(contactCoordinates, boltPts, varargin)
% axes3DModel = PlotElectrodeModel(contactCoordinates, varargin)
% 
% Description:
%   Plots a 3D model of an SEEG electrode given a Nx3 matrix of coordinates for the N contact positions
%
% Usage:
%   >> PlotElectrodeMode(contactCoordinates, boltPots);
%
% Output:
%   axes3DModel - axes handle for the rendered model
%   
% Required Parameters:
%   contactCoordinates  - N x 3 matrix of 3-dimensional spacial coordinates for the N contacts that make up the electrode
%   boltPts             - N x 3 matrix of 3-dimensional spacial coordinates for the N control points that approximate the bolt 
%
% Optional Parameters:
%   useAxes - (axes hanlde) pass an axes handle on which to plot the 3D electrode model
%   
% Author:
%   Jarod L Roland
%   Department of Neurosurgery
%   Washington University in St. Louis
%
params = inputParser;
addRequired(params, 'contactCoordinates', @(x) size(x, 2) == 3); % must be an N x 3 matrix
addRequired(params, 'boltPts', @(x) size(x, 2) == 3); % must be an N x 3 matrix
addParameter(params, 'useAxes', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));

parse(params, contactCoordinates, boltPts, varargin{:});

%% initialize variables
axes3DModel = params.Results.useAxes;

electrodeDiameter = 0.7; %0.3;
backgroundColor = [1, 1, 1]; %[0, 0, 0]; %
electrodeColor = [0.5, 0.5, 0.5];
contactColor = [0.9, 0.9, 0.9];
boltColor = [0, 0, 0];

cylResolutionLength = 10; %100;              % resolution of cylinder along its length (higher number generates smoother mesh)
cylResolutionCircumference = 10; %30;        % resolution of cylinder along its diameter (higher number generates smoother mesh)

% calculate a vector point from the Nth to the N-1th contact
vectorAdjacentElectrodes = diff(contactCoordinates);
% meanDistanceBetweenElectrodes = mean(vecnorm(vectorAdjacentElectrodes, 2, 2)); % 2-norm across dimension 2. N.B. some SEEG electrodes have different gaps between contacts

numElectrodes = length(contactCoordinates);
numBoltPts = size(boltPts, 1);
numControlPts = 1 + numElectrodes + 1 + numBoltPts;     % one in front (distal) + all electrodes (middle) + all bolts (proximal)

vecTowardsBolt = contactCoordinates(end, :) - boltPts(end, :);
vecTowardsBolt = vecTowardsBolt ./ norm(vecTowardsBolt);
vecTowardsBolt = vecTowardsBolt .* norm(vectorAdjacentElectrodes(end, :));
proximalControlPoint = contactCoordinates(end, :) - vecTowardsBolt;


% setup control points for spline interpolation
controlCoordinates = zeros(numControlPts, 3);
controlCoordinates(2:length(contactCoordinates)+1, :) = contactCoordinates;
controlCoordinates(1, :) = contactCoordinates(1,:) - vectorAdjacentElectrodes(1, :);    % add a control point just beyond the 1st (most distal) contact
% cappedCoordinates(end, :) = contactCoordinates(end,:) + vectorAdjacentElectrodes(end, :);
controlCoordinates(length(contactCoordinates)+2, :) = proximalControlPoint;             % add a control point that is one interElectrodeDistance toward the bolt
controlCoordinates(length(contactCoordinates)+3:end, :) = boltPts;                      % add controls points for electrode going through bolt

%% interpolate a spline through control points
splineSteps = linspace(1, numControlPts, (numControlPts - 1) * cylResolutionLength + 1);
splineX = interp1(1:numControlPts, controlCoordinates(:, 1), splineSteps, 'spline');
splineY = interp1(1:numControlPts, controlCoordinates(:, 2), splineSteps, 'spline');
splineZ = interp1(1:numControlPts, controlCoordinates(:, 3), splineSteps, 'spline');

% create mesh of electrode as a cylinder around the spline path
[modelElectrodeX, modelElectrodeY, modelElectrodeZ] = gencyl( [splineX; splineY; splineZ], ones(length(splineX), 1) * electrodeDiameter, cylResolutionLength, cylResolutionCircumference);
modelElectrodeX = modelElectrodeX((3 * cylResolutionLength):end, :);     %MAGICNUMBER: 3 x cylResolutionLength is an empiric length of the electrode past the last contact
modelElectrodeY = modelElectrodeY((3 * cylResolutionLength):end, :);
modelElectrodeZ = modelElectrodeZ((3 * cylResolutionLength):end, :);
[electrodeModel.faces, electrodeModel.vertices, ~] = surf2patch(modelElectrodeX, modelElectrodeY, modelElectrodeZ, 'triangles');

%% plot the 3D model of the electrode
if(isempty(axes3DModel))
    figure();
    axes3DModel = axis();
end
axes3DModel.Parent.Color = backgroundColor;
axes3DModel.Color = backgroundColor;
surf(axes3DModel, modelElectrodeX, modelElectrodeY, modelElectrodeZ, 'FaceColor', electrodeColor,'FaceLighting','gouraud','EdgeColor','none');
axis(axes3DModel, 'tight');
grid(axes3DModel, 'on');
axes3DModel.LineWidth = 2;
hold on;
material metal;

% % TESTING: plot red circle at contact points
% plot3(contactCoordinates(:, 1), contactCoordinates(:, 2), contactCoordinates(:, 3), 'ro', 'MarkerSize', 20)

% % TESTING: plot green circle at bolt points
% plot3(boltPts(:, 1), boltPts(:, 2), boltPts(:, 3), 'go', 'MarkerSize', 20)

% % TESTING: plot blue circle at control points
% plot3(controlCoordinates(:, 1), controlCoordinates(:, 2), controlCoordinates(:, 3), 'bo', 'MarkerSize', 30)

% % TESTING: plot black circle at proximal control point (virtual step in the spline between proximal contact and bolt)
% plot3(proximalControlPoint(:, 1), proximalControlPoint(:, 2), proximalControlPoint(:, 3), 'ko', 'MarkerSize', 35);

% plot the contacts as a similar model with different color and slightly larger diameter than the rest of the electrode
contactsModels(numElectrodes).faces = [];
contactsModels(numElectrodes).vertices = [];
% contactsModels(numElectrodes).colors = [];
for contactNum = 1:numElectrodes
    idx = contactNum * cylResolutionLength + 1; % skip the first and last index of spline control points as these are simply references points to define the spline past the first and last electrode contact
    contactHalfLength = floor(cylResolutionLength * 0.35);  %MAGNICNUMBER: contact segment length equal 30% x 2 = 60% (i.e., 30% forward + 30% backward) of the distance bewteen adjacent contacts (cylResolution defines the number of cylinder vertices are drawn bewteen adjacent spline nodes)
    [contactX, contactY, contactZ] = gencyl( [splineX(idx-contactHalfLength:idx+contactHalfLength); splineY(idx-contactHalfLength:idx+contactHalfLength); splineZ(idx-contactHalfLength:idx+contactHalfLength)], ones(contactHalfLength*2+1, 1) * electrodeDiameter * 1.1, cylResolutionLength, cylResolutionCircumference); %MAGICNUMBER: 1.1 - expand the diameter by 10% to visualize this segment over the rest of the electrode
%     [contactsModels(contactNum).faces, contactsModels(contactNum).vertices, contactsModels(contactNum).contactColor] = surf2patch(contactX, contactY, contactZ);
    [contactsModels(contactNum).faces, contactsModels(contactNum).vertices, ~] = surf2patch(contactX, contactY, contactZ, 'triangles');
    surf(axes3DModel, contactX, contactY, contactZ, 'FaceColor', contactColor, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
    material('shiny');
end

% plot the bolts
boltModel.faces = [];
boltModel.vertices = [];
% contactsModels(numElectrodes).colors = [];
startIdx = (numElectrodes + 2) * cylResolutionLength;
endIdx = (numControlPts - 1) * cylResolutionLength;

[boltX, boltY, boltZ] = gencyl( [splineX(startIdx:endIdx); splineY(startIdx:endIdx); splineZ(startIdx:endIdx)], ones((endIdx-startIdx)+1, 1) * electrodeDiameter * 1.1, cylResolutionLength, cylResolutionCircumference); %MAGICNUMBER: 1.1 - expand the diameter by 10% to visualize this segment over the rest of the electrode
[boltModel.faces, boltModel.vertices, ~] = surf2patch(boltX, boltY, boltZ, 'triangles');
surf(axes3DModel, boltX, boltY, boltZ, 'FaceColor', boltColor, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
material('metal');


end %function PlotElectrodeModel(electrodeCoordinates)
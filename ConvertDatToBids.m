% Convert .dat files as saved from FreeView when manually localizing
% electrode contacts to a BIDS format electrode definition file

% init variables
dataGroup = 'slch';
projDir = GetProjectPath();
bidsDir = fullfile(projDir, 'data', 'bids', dataGroup);

% load list of subjects
data = open('subjectListPriorCranSLCH.mat');
subjListAll = vertcat(data.subjListNaiveBone, data.subjListPriorCran);
numSubj = length(subjListAll);

% loop through all subjects
for idxSubj = 1:numSubj
    subjID = subjListAll{idxSubj};
    
    % get electrode data for subject
    [electrodeData, axes3DSimple, axes3DModel, axes2DErrors] = CalcSEEGTrajectory(subjID, 'isFigVisible', false);
    close(axes3DSimple.Parent, axes3DModel.Parent, axes2DErrors.Parent);
    
    % build paths and open bids files for writing
    pathSubjDir = fullfile(bidsDir, ['sub-' subjID]);
    pathIeegDir = fullfile(pathSubjDir, 'ieeg');
    if(not(isfolder(pathIeegDir)))  % make sure directory exists
        mkdir(pathIeegDir);
    end
    pathBidsElectrodesFile = fullfile(pathIeegDir, ['sub-' subjID '_electrodes.tsv']);
    fidElectrodes = fopen(pathBidsElectrodesFile, 'w+');
    assert(fidElectrodes, ['Error: Error opening BIDS electrodes file (' pathBidsElectrodesFile ')']);
    fprintf(fidElectrodes, 'name\tx\ty\tz\tsize\ttype\n');

    % loop through all contacts on all electrodes
    numElectrodes = length(electrodeData);
    for idxElectrode = 1:numElectrodes
        numContacts = length(electrodeData(idxElectrode).contactsCoords);        
        for idxContacts = 1:numContacts
            fprintf(fidElectrodes, '%s%d\t%.3f\t%.3f\t%.3f\tunknown\tSEEG\n', ...
                electrodeData(idxElectrode).label, idxContacts, ...
                electrodeData(idxElectrode).contactsCoords(idxContacts, 1), ...
                electrodeData(idxElectrode).contactsCoords(idxContacts, 2), ...
                electrodeData(idxElectrode).contactsCoords(idxContacts, 3));
        end %for idxContacts = 1:numContacts
        
    end %for electrodeIdx = 1:numElectrodes
    fclose(fidElectrodes);
    fidElectrodes = 0;
    
    pathBidsCoordSystemFile = fullfile(pathIeegDir, ['sub-' subjID '_coordsystem.json']);
    fidCoordSystem = fopen(pathBidsCoordSystemFile, 'w+');
    assert(fidCoordSystem, ['Error: Error opening BIDS coordsystem file (' pathBidsCoordSystemFile ')']);
    fprintf(fidCoordSystem, ...
        ['{\n' ...
        '\t"iEEGCoordinateSystem": "Other",\n' ...
        '\t"iEEGCoordinateUnits": "mm",\n' ...
        '\t"iEEGCoordinateSystemDescription": "Coordinates from native post-op CT image space. Not aligned to any other imaging or atlas."\n', ...
        '}\n'] ...
        );
    fclose(fidCoordSystem);
    fidCoordSystem = 0;

end %for subjIdx = 1:numSubj

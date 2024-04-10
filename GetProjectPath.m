function projPath = GetProjectPath()
% Returns the root of the project directory depending on which machine it's being run on. 
% ASSUMPTION: assumes this script is being *run* from a 'code' folder within the project folder. e.g.: working directory = /home/user/projects/projectName/code/

% get project name from assumed current working directory path structure
pathParts = split(pwd, filesep);
projectDirectory = pathParts(end-1);
projectName = projectDirectory{1};

% generate path based on hostname
[~, hostname] = system('hostname');
hostname = strtrim(hostname);               % remove the Line Feed (ascii 10) at the end of the string
projPath = fullfile(getenv('HOME'), 'projects', projectName);


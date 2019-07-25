% Load previous session

function handles=MODAload

[filename,pathname] = uigetfile('*.*');  % Request file location from user
name = fullfile(pathname,filename);

handles = load(name); % Load data
function [t,Ffoam_top,Ffoam_bottom,Fx_tail, Fy_tail, T_tail,Fx_rebound, Fy_rebound, Fx_contact, Fy_contact, Fy_fricTop, Fy_fricBottom,F_hardstopTop, F_hardstopBottom] = BerkeleyImpact_import3_forces(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [T,FFOAM_TOP,FFOAM_BOTTOM,FX_TAIL,T_TAIL,FMAG_REBOUND] =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   [T,FFOAM_TOP,FFOAM_BOTTOM,FX_TAIL,T_TAIL,FMAG_REBOUND] =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [t,Ffoam_top,Ffoam_bottom,Fx_tail,T_tail,Fmag_rebound] =
%   importfile('BerkeleyImpact.3',5, 98);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/07/12 13:59:33

%% Initialize variables.
if nargin<=2
    startRow = 5;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%	column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
t = dataArray{:, 1};
Ffoam_top = dataArray{:, 2};
Ffoam_bottom = dataArray{:, 3};
Fx_tail = dataArray{:, 4};
Fy_tail = dataArray{:, 5};
T_tail = dataArray{:, 6};
Fx_rebound = dataArray{:, 7};
Fy_rebound = dataArray{:, 8};
Fx_contact = dataArray{:, 9};
Fy_contact = dataArray{:, 10};
Fy_fricTop = dataArray{:, 11};
Fy_fricBottom = dataArray{:, 12};
F_hardstopTop = dataArray{:, 13};
F_hardstopBottom = dataArray{:, 14};


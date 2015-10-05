function [vx,vy,vgamma,vphi,ax,ay,agamma,aphi] = BerkeleyImpact_import6_vel_acc(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [VX,VY,VGAMMA,VPHI,AX,AY,AGAMMA,APHI,VARNAME9] = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   [VX,VY,VGAMMA,VPHI,AX,AY,AGAMMA,APHI,VARNAME9] = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   [vx,vy,vgamma,vphi,ax,ay,agamma,aphi,VarName9] =
%   importfile('BerkeleyImpact_TallBody.6',5, 205);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/09/25 12:51:31

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
%   column7: double (%f)
%	column8: double (%f)
%   column9: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14f%15f%15f%15f%15f%15f%15f%15f%s%[^\n\r]';

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
vx = dataArray{:, 1};
vy = dataArray{:, 2};
vgamma = dataArray{:, 3};
vphi = dataArray{:, 4};
ax = dataArray{:, 5};
ay = dataArray{:, 6};
agamma = dataArray{:, 7};
aphi = dataArray{:, 8};
VarName9 = dataArray{:, 9};


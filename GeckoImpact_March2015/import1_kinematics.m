function [t,x,y,x_arm,y_arm,x_tail,y_tail,gamma1,theta,phi,VarName11] = import1_kinematics(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [T,X,Y,X_ARM,Y_ARM,X_TAIL,Y_TAIL,GAMMA1,THETA,PHI,VARNAME11] =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   [T,X,Y,X_ARM,Y_ARM,X_TAIL,Y_TAIL,GAMMA1,THETA,PHI,VARNAME11] =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [t,x,y,x_arm,y_arm,x_tail,y_tail,gamma1,theta,phi,VarName11] =
%   importfile('GeckoImpact.1',5, 205);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/03/26 14:16:42

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
%   column9: double (%f)
%	column10: double (%f)
%   column11: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14f%15f%15f%15f%15f%15f%15f%15f%15f%15f%s%[^\n\r]';

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
x = dataArray{:, 2};
y = dataArray{:, 3};
x_arm = dataArray{:, 4};
y_arm = dataArray{:, 5};
x_tail = dataArray{:, 6};
y_tail = dataArray{:, 7};
gamma1 = dataArray{:, 8};
theta = dataArray{:, 9};
phi = dataArray{:, 10};
VarName11 = dataArray{:, 11};

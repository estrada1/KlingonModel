function [x_FoamTop,y_FoamTop,x_FoamBottom,y_FoamBottom,AttachPt_x_WorldFrame,AttachPt_y_WorldFrame, x_hardstop, x_Ccm, y_Ccm, x_top, y_top, x_bottom, y_bottom ] = BerkeleyImpact_import2_kinematics(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [X_FOAMTOP,Y_FOAMTOP,X_FOAMBOTTOM,Y_FOAMBOTTOM,ATTACHPT_X_WORLDFRAME,ATTACHPT_Y_WORLDFRAME]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [X_FOAMTOP,Y_FOAMTOP,X_FOAMBOTTOM,Y_FOAMBOTTOM,ATTACHPT_X_WORLDFRAME,ATTACHPT_Y_WORLDFRAME]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [x_FoamTop,y_FoamTop,x_FoamBottom,y_FoamBottom,AttachPt_x_WorldFrame,AttachPt_y_WorldFrame]
%   = importfile('BerkeleyImpact.2',5, 52);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2015/07/12 12:57:34

%% Initialize variables.
if nargin<=2
    startRow = 5;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%14s%15s%15s%15s%15s%17s%15s%15s%15s%15s%15s%15s%[^\n\r]';

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

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9, 10, 11, 12, 13]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Allocate imported array to column variable names
x_FoamTop = cell2mat(raw(:, 1));
y_FoamTop = cell2mat(raw(:, 2));
x_FoamBottom = cell2mat(raw(:, 3));
y_FoamBottom = cell2mat(raw(:, 4));
AttachPt_x_WorldFrame = cell2mat(raw(:, 5));
AttachPt_y_WorldFrame = cell2mat(raw(:, 6));
x_hardstop = cell2mat(raw(:, 7));
x_Ccm = cell2mat(raw(:, 8));
y_Ccm = cell2mat(raw(:, 9));
x_top = cell2mat(raw(:, 10));
y_top =  cell2mat(raw(:, 11));
x_bottom =  cell2mat(raw(:, 12));
y_bottom = cell2mat(raw(:, 13));


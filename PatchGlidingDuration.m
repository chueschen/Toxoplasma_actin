% Calculating patch gliding duration and number of reversals

clearvars
%select directory, starting from chosen base path, for tracking results. 
% Tracking results are in a .csv file whose name starts with "Results..." These are "spot tracks"
% exports from the ImageJ Manual Tracking plugin https://imagej.net/ij/plugins/track/track.html
folder_name = uigetdir('/Users/...patchDuration/reversals'); 
cd(folder_name);

%% 1. Import tracking data and metadata
% need to get: timestamps (save)

% this struct saves info about reversal tracking results and metadata
dirFileListTracks = dir([folder_name, '/R*.csv']);
dirFileListMeta = dir([folder_name, '/CH0*.txt']);

% import velocity tracking data into struct
for i = 1:size(dirFileListTracks,1)
    data(i).tracks = readmatrix(dirFileListTracks(i).name);
end
% columns are track no, frame no, x, y, distance, velocity (pixels/frame), pixel value.

% Pull out time stamps from text file metadata for all movies
for k = 1:length(dirFileListMeta)
    filetext = fileread(dirFileListMeta(k).name);
    expr = 'ElapsedTime-ms":\s*\d\d\d....'; % note that some time stamps are 3 digits and some up to 7! So final '.' takes any digit, which might be a comma
    matches = regexp(filetext,expr,'match'); % pull out matches to the expression expr
    for i = 1:length(matches)
        times(i,k) = str2double(matches{i}(17:23)); % pull out numbers from matches (in ms) into columns of time for each movie
    end
end
times = times ./ 1000; % convert times to seconds

%% Calculate duration of patch gliding for each event (cell-movie)
%singPix = 0.1; % pixel size in um for TIRF 100x Andor

for i = 1:size(data,2) % for struct, apparently dim 2 is num cells
  for j = 1:size(data(i).tracks,1) % loop over tracking points
 % Real time in seconds - add as column 8
    data(i).tracks(j,8) = times(data(i).tracks(j,2),i);
  end
 % Time (s) between reversals - add as column 9
   for j = 2:size(data(i).tracks,1) % loop over intervals, start at 2
    data(i).tracks(j,9) = data(i).tracks(j,8) - data(i).tracks(j-1,8);
   end
   
   % Pull out total duration into independent vector
    totalDur(i) = data(i).tracks(end,8) - data(i).tracks(1,8);
   % Calculate number of reversals for each patch gliding event - new vector
    numRev(i) = size(data(i).tracks,1) - 2; % NOTE THAT THIS IS 1 LESS THAN BODY LENGTH GLIDES
end

%% Annotate gliding events that extend beyond movie, so total duration is longer than measured
% note: needed to open imported metadata struct and make sure to match movies
for i = [2,4,8:10,12:17] % cell-movies manually annotated as extending beyond recording
    data(i).longer = 1; %mark with 1 in data struct
end


%% Plot scatter of total duration vs number of reversals
close all

figure; hold on
for i = 1:size(data,2)
    if data(i).longer == 1 % plot extended movies with +
        plot(numRev(i), totalDur(i), '+', 'MarkerSize', 10) %, 'LineWidth', 2)
    else % plot others with 0
        plot(numRev(i), totalDur(i), 'o', 'MarkerSize', 10) %, 'LineWidth', 2)
    end
end
xlabel('number of direction reversals')
ylabel('duration of patch gliding event (s)')




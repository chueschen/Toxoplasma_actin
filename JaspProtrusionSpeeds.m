% Calculating and Plotting Speeds of Jasplakinolide-Induced Actin
% Protrusions ("horns")

% I. Calculate and plot instantaneous velocity in cell long-axis vs. time
% (over 2 cycles) for 3 examples (Fig. S4A)

% II. Calculate and plot mean speed for 10 cells (Fig. S4B)

%% %%%%% PART I: Calculate and plot instantaneous velocity in cell long-axis vs. time (over 2 cycles) for 3 examples  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
%select directory where tracking results are saved, starting from chosen base
%path. Tracking results are in a .csv file whose name starts with "Results..." These are "spot tracks"
% exports from the ImageJ Manual Tracking plugin https://imagej.net/ij/plugins/track/track.html

folder_name = uigetdir('/Users/.../narwhal_horn_speeds/3examples'); 
cd(folder_name);

%% Import tracking data and metadata

% this struct saves info about tracking results and metadata
dirFileListTracks = dir([folder_name, '/R*.csv']);
dirFileListMeta = dir([folder_name, '/CH0*.txt']);

% import velocity tracking data into one array
tracks = readmatrix(dirFileListTracks(1).name);
temp = readmatrix(dirFileListTracks(2).name);
tracks = cat(1, tracks, temp); clear temp
% columns are track no, frame no, x, y, distance, velocity (pixels/frame), pixel value.
% each cell has 4 track IDs: horn/protrusion base, center of cell, posterior pole, apical pole

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


%% Organize tracking data

% Pull out and organize tracking data into structs. For each example, the
% first object tracked was the horn/protrusion, then the center of the cell
% nucleus, then one cell pole, then the other.
%set indices
    z=1;
    k=1; % start counter for # of cells
    %find TrackID
    currTrackID=tracks(z,1);
    %pulls out vector of row #s for that trackID
    currIndex=find(tracks(:,1)==currTrackID);
%    
    while z <= size(tracks,1)
        currTrackID=tracks(z,1);
        currIndex=find(tracks(:,1)==currTrackID);
        data(k).horn = tracks(currIndex(1):currIndex(end),:);
       %find next trackID
        z=currIndex(end)+1; 
        currTrackID=tracks(z,1);
        currIndex=find(tracks(:,1)==currTrackID);
        data(k).cell = tracks(currIndex(1):currIndex(end),:); % cell track
       %find next trackID - posterior pole
        z=currIndex(end)+1; 
        currTrackID=tracks(z,1);
        currIndex=find(tracks(:,1)==currTrackID);
        data(k).p = tracks(currIndex(1):currIndex(end),:); % posterior pole
       %find next trackID - apical pole
        z=currIndex(end)+1; 
        currTrackID=tracks(z,1);
        currIndex=find(tracks(:,1)==currTrackID);
        data(k).a = tracks(currIndex(1):currIndex(end),:);
       % increase cell count
        k = k+1;
       %find next trackID - prep for next horn
        z=currIndex(end)+1; 
    end


%% Calculate horn velocities, horn-cell velocities, cell long axis, velocity in long axis
singPix = 0.1; % pixel size in um for TIRF 100x Andor

% First, calculate horn positions - cell positions as columns 8 and 9
for i = 1:size(data,2) % for struct, dim 2 is num cells
  for j = 1:size(data(i).horn,1)
    data(i).horn(j,8) = data(i).horn(j,3) - data(i).cell(j,3); % horn x - cell x position
    data(i).horn(j,9) = data(i).horn(j,4) - data(i).cell(j,4); % horn y - cell y position
    
    % I'm also going to pull out time as column 15, for convenience, and
    % have them start at 0 on row 2, when velocity starts. Not used until plotting.
    clear tempTime0
    tempTime0 = times(data(i).horn(2,2),i);
    data(i).horn(j,15) = times(data(i).horn(j,2),i) - tempTime0;
  end
end

% Now, digging into the main stuff
for i = 1:size(data,2) % for struct, apparently dim 2 is num cells
  for j = 2:size(data(i).horn,1)
    
    % 1. Horn velocity in lab frame in um/s - add as column 10
    data(i).horn(j,10) = data(i).horn(j,6) * singPix / ( times(data(i).horn(j,2),i)-times(data(i).horn(j,2)-1,i) );

    % 2. Horn velocity in cell frame in um/s - add as column 11
    data(i).horn(j,11) = sqrt((data(i).horn(j,8)-data(i).horn(j-1,8))^2+(data(i).horn(j,9)-data(i).horn(j-1,9))^2)...
        * singPix / ( times(data(i).horn(j,2),i)-times(data(i).horn(j,2)-1,i) ); % horn velocity in um/s
    
    % 3. Calculate unit vector in direction of cell long axis, posterior-to-apical. 
    % x component as column 12, y component as column 13
    data(i).horn(j,12) = data(i).a(j,3) - data(i).p(j,3); % x of a minus x of p
    data(i).horn(j,13) = data(i).a(j,4) - data(i).p(j,4); % y of a minus y of p
    clear tempVecMag
    tempVecMag = sqrt(data(i).horn(j,12)^2 + data(i).horn(j,13)^2);
    data(i).horn(j,12) = data(i).horn(j,12) / tempVecMag; % norm to get unit vector
    data(i).horn(j,13) = data(i).horn(j,13) / tempVecMag; % norm to get unit vector
    
    % 4. Calculate one-dimension horn velocity in a-p long axis -- column 14
    % USING HORN - CELL VELOCITIES
    tempXcomp = (data(i).horn(j,8)-data(i).horn(j-1,8)) * data(i).horn(j,12); % dot product to get component in a-p
    tempYcomp = (data(i).horn(j,9)-data(i).horn(j-1,9)) * data(i).horn(j,13);
    data(i).horn(j,14) = tempXcomp + tempYcomp; % value of vector in a-p axis.
  end
end


%% Plots!

% First one-dim velocity vs. time
figure
hold on
for i = 1:size(data,2)
    plot(data(i).horn(2:end,15), data(i).horn(2:end,14),'.', 'MarkerSize', 15)
end
xlabel('time (s)')
ylabel('velocity in cell long-axis (\mum/s)')
legend('cell # 1','cell # 2', 'cell # 3')
set(gca,'FontSize',18)
ylim([-15,15])


%% %%%%% NOW, PART II: Calculate and plot mean speed for all cells (Fig. S4B) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
%select directory, starting from chosen base path
folder_name = uigetdir('/Users/...narwhal_horn_speeds/meanSpeeds'); 

%% Import tracking data and metadata
% need to get: timestamps (save)

% this struct saves info about horn tracking results and metadata
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

%% Calculate speeds
singPix = 0.1; % pixel size in um for TIRF 100x Andor

for i = 1:size(data,2) % for struct, apparently dim 2 is num cells
  for j = 2:size(data(i).tracks,1)
 % Horn velocity in lab frame in um/s - add as column 8
 % note - make sure time interval is calling from correct frames, NOT
 % subsequent (frames were skipped in tracking)
    data(i).tracks(j,8) = data(i).tracks(j,5) * singPix / ( times(data(i).tracks(j,2),i)-times(data(i).tracks(j-1,2),i) );
  end
end

% Calculate mean speeds for each cell
for i = 1:size(data,2) % for struct, dim 2 is num cells
    meanSpeeds(i) = mean(data(i).tracks(2:end,8));
end


%% Plotting all speeds
speeds = [];

for i = 1:size(data,2)
   tempSpeeds = data(i).tracks(2:end,8);
   speeds = [speeds; tempSpeeds];
   clear tempSpeeds
end
 dataPlot = {speeds};
 figure
 plotSpread(dataPlot, 'xNames', {''}, ...
'distributionMarkers', {'o'},'distributionColors',{'m'}, 'showMM',5); %"showMM',5 adds stdv (4 = sem)

ylabel('mean speed of recirculating actin bundle (\mum/s)')
set(gca,'FontSize',18)


%% Mean and stdev of pooled speeds
meanSpeedPooled = mean(speeds);
stdevSpeedPooled = std(speeds);



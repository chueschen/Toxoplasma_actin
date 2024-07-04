%%%% Analyzing fast, directional actin speckle tracks %%%%

folder_name = uigetdir('/Volumes/'); % navigate to folder ("ManualActinTrackingActive") containing a numbered subfolder 
% for every movie. In each subfolder, there's a movie of the actin speckle
% channel and one of the bulk labeled cell, plus manual tracking results in
% a .csv file whose name starts with "Results..." These are "spot tracks"
% exports from the ImageJ Manual Tracking plugin https://imagej.net/ij/plugins/track/track.html
cd(folder_name)

%% 1. Import manual tracking output and save info, dp vectors, theta at each time step
% directory location, frame #s, time step, track coordinates, etc.

%clearvars
singPix = 0.1; % pixel size in um
myMovies = 1:21; % input('How many movies do you have?'); %currently 10 actin, 11 myosin. Need to recheck that full actin movies were tracked

for i = myMovies
    cd(folder_name) % reset to tracking directory
    cd(num2str(i)); % move to directory for first movie
    dirFileList = dir('R*.csv'); % using dir even though one file b/c it allows wildcards
    currFileName = dirFileList.name; %find filename
    spotTables(i).Tables = readtable(currFileName);  
end

%% 2. Organize spotTables into one matrix

for i = myMovies

    %set indices
    z=1;
    j=1; % start counter for # of cells/tracks in this one movie
    %find TrackID
    currTrackID=spotTables(i).Tables.TrackN_(z);
    %pulls out vector of row #s for that trackID
    currIndex=find(spotTables(i).Tables.TrackN_==currTrackID);
    
    while currIndex(end) < size(spotTables(i).Tables,1)
        allTracks(i).tracks(j).x = spotTables(i).Tables.X(currIndex(1):currIndex(end)); % pull out x coords
        allTracks(i).tracks(j).y = spotTables(i).Tables.Y(currIndex(1):currIndex(end)); % pull out y coords
        allTracks(i).tracks(j).f = spotTables(i).Tables.SliceN_(currIndex(1):currIndex(end)); % pull out frame numbers
        allTracks(i).tracks(j).u = spotTables(i).Tables.Velocity(currIndex(2):currIndex(end)); % pull out plugin-calculated velocities
        tempFrameInt = (spotTables(i).Tables.SliceN_(currIndex(2):currIndex(end))... % find interval between positions, in # of frames
            - spotTables(i).Tables.SliceN_(currIndex(1):currIndex(end-1)));
        allTracks(i).tracks(j).u = allTracks(i).tracks(j).u ./ tempFrameInt; % **fixing velocity calculation problem from plugin** accounting for longer time during gaps
        
        %find next trackID
        j=j+1;
        z=currIndex(end)+1; 
        currTrackID=spotTables(i).Tables.TrackN_(z);
        %pulls out vector of row #s for that trackID
        clear currIndex;
        currIndex=find(spotTables(i).Tables.TrackN_==currTrackID);
    end
    
    if currIndex(end) == size(spotTables(i).Tables,1)
        allTracks(i).tracks(j).x = spotTables(i).Tables.X(currIndex(1):currIndex(end)); % pull out x coords
        allTracks(i).tracks(j).y = spotTables(i).Tables.Y(currIndex(1):currIndex(end)); % pull out y coords
        allTracks(i).tracks(j).f = spotTables(i).Tables.SliceN_(currIndex(1):currIndex(end)); % pull out frame numbers
        allTracks(i).tracks(j).u = spotTables(i).Tables.Velocity(currIndex(2):currIndex(end)); % pull out plugin-calculated velocities
        allTracks(i).frameInterval = spotTables(i).Tables.Distance(currIndex(end)) ./ spotTables(i).Tables.Velocity(currIndex(end)); % pull out dt
        tempFrameInt = (spotTables(i).Tables.SliceN_(currIndex(2):currIndex(end))... % find interval between positions, in # of frames
            - spotTables(i).Tables.SliceN_(currIndex(1):currIndex(end-1)));
        allTracks(i).tracks(j).u = allTracks(i).tracks(j).u ./ tempFrameInt; % **fixing velocity calculation problem from plugin** accounting for longer time during gaps
 
    end
    
    % Calculate dp vectors (instantaneous velocity) and theta (angle of dp from -pi to pi)
    for j = 1:size(allTracks(i).tracks, 2)
    trackLength = length(allTracks(i).tracks(j).x);
        for k = 1:trackLength - 1 %loop through all intervals between frames of track  
            % Calculate velocity vector for each dt
            p1(1) = allTracks(i).tracks(j).x(k); % x-coord of first point, time t
            p1(2) = allTracks(i).tracks(j).y(k); % y-coord of first point
            p2(1) = allTracks(i).tracks(j).x(k+1); % x-coord of 2nd point, time t+1
            p2(2) = allTracks(i).tracks(j).y(k+1); % y-coord of 2nd point
            allTracks(i).tracks(j).dp(k,:) = p2-p1; % save velocity for later as vectors [dx, dy] (first column dx = v_x/dt)
            % Calculate theta, orientation of velocity vector for each dt, in "standard" coordinates where positive y is up
            allTracks(i).tracks(j).theta(k,:) = atan2( -1 * allTracks(i).tracks(j).dp(k,2) , (allTracks(i).tracks(j).dp(k,1)) ); % atan of -dy, dx. NOTE! flipping dy to match real space coordinate system
        end
    end
end


%% 3. Read in frames of bulk-labeled movie showing cell shape, and define polarity.
% Use function segmentAlign to choose tail-to-head vector and use that to rotate image. Segment cell and use center of mass
% to translate image. Also, pull out cropped and aligned images of bulk cell and store in imAlignedBulk.

cd(folder_name) % reset to actin tracking directory

if exist('imAlignedBulk', 'var') == 0
    imAlignedBulk = [];
end

for i = myMovies
    [allTracks, imAlignedBulk] = segmentAlignManualTracks(allTracks, i, imAlignedBulk);
end
% Note that rotation vector and translation vector are stored in allTracks
% Currently, cell (ellipse) centroids are at position (50,50)
% Note that each bulk image has been cropped to 11:89, so centroid is now at 50,50 and image is 79 x 79

%% Normalize bulk images in imAlignedBulk
for i = myMovies
    imAlignedBulk(:,:,i) = imAlignedBulk(:,:,i) ./ max(max(imAlignedBulk(:,:,i)));
end


%% 4. Rotate and translate tracks.

% Rotate and translate x and y coords from u-track output
for i = myMovies
    numTracks = size(allTracks(i).tracks, 2);
    for j = 1:numTracks %loop through all tracks
        trackLength = length(allTracks(i).tracks(j).x);
        for k = 1:trackLength %loop through all frames of track  
        point = [allTracks(i).tracks(j).y(k) ; allTracks(i).tracks(j).x(k)]; % pull each pt (x,y) into [i;j] = [y;x]
        rotatedPoint = allTracks(i).tracks(j).rotMatrix*(point-allTracks(i).tracks(j).centerImOG)+allTracks(i).tracks(j).centerRotIm; % new coordinates of point after rotation
        translatedPoint = rotatedPoint + allTracks(i).tracks(j).translVector;
        allTracks(i).tracks(j).alignedTracksYXF(1,k) = translatedPoint(1) - 11; %keeping convention where first row is i = y-coord. 11 is currently hard-coded as crop box start
        allTracks(i).tracks(j).alignedTracksYXF(2,k) = translatedPoint(2) - 11; %keeping convention where second row is j = x-coord
        allTracks(i).tracks(j).alignedTracksYXF(3,k) = allTracks(i).tracks(j).f(k); %store matching movie frame # for convenience, in case I make movies later
        end
    end
end


%% Plots tracks on top of reference images

imMeanL = zeros(size(imAlignedBulk,1),size(imAlignedBulk,2),1); 
imMeanR = zeros(size(imAlignedBulk,1),size(imAlignedBulk,2),1); 
imMeanB = zeros(size(imAlignedBulk,1),size(imAlignedBulk,2),1);

%calculate mean reference images
for i = myMovies
    for j = 1:size(allTracks(i).tracks, 2)

        if allTracks(i).tracks(j).LRB == 'l'
            imMeanL = imMeanL + imAlignedBulk(:,:,allTracks(i).tracks(j).count);
        end
        if allTracks(i).tracks(j).LRB == 'b'
            imMeanB = imMeanB + imAlignedBulk(:,:,allTracks(i).tracks(j).count);
        end
        if allTracks(i).tracks(j).LRB == 'r'
            imMeanR = imMeanR + imAlignedBulk(:,:,allTracks(i).tracks(j).count);
        end
    end
end

% plot left tracks
figure()
imshow(imMeanL, []); title('reference image: left side'); % show image to check
hold on
for i = myMovies
    for j = 1:size(allTracks(i).tracks, 2)
        if allTracks(i).tracks(j).LRB == 'l'
            col = rand(1,3); % generate random RGB color to use for this track
            for k = 1:length(allTracks(i).tracks(j).x) - 1 %loop through all intervals between frames of track  
                p1(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k); % x-coord of first point, time t
                p1(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k); % y-coord of first point
                p2(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k+1); % x-coord of 2nd point, time t+1
                p2(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k+1); % y-coord of 2nd point
                dp = p2-p1; % vector of delta_x, delta_y
                allTracks(i).tracks(j).dpAligned(k,:) = dp; % save velocity for later as vectors [dx, dy] (first column dx = v_x/dt)
                quiver(p1(1),p1(2),dp(1),dp(2), 0, 'MaxHeadSize',3/norm(dp),'color',col) % draw arrows between points 
                    % ^ 0 means no scaling, so will actually connect points. norm() finds length of vector dp.
                    % Defining MaxHeadSize is not quite working to get same-size arrowheads
            end
        end
    end
end
set(gcf,'Position',[10 10 500 500])
hold off
pause(1)

% plot right tracks
figure()
imshow(imMeanR, []); title('reference image: right side'); % show image to check
hold on
for i = myMovies
    for j = 1:size(allTracks(i).tracks, 2)
        if allTracks(i).tracks(j).LRB == 'r'
            col = rand(1,3); % generate random RGB color to use for this track
            for k = 1:length(allTracks(i).tracks(j).x) - 1 %loop through all intervals between frames of track  
                p1(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k); % x-coord of first point, time t
                p1(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k); % y-coord of first point
                p2(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k+1); % x-coord of 2nd point, time t+1
                p2(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k+1); % y-coord of 2nd point
                dp = p2-p1; % vector of delta_x, delta_y
                allTracks(i).tracks(j).dpAligned(k,:) = dp; % save velocity for later as vectors [dx, dy] (first column dx = v_x/dt)
                quiver(p1(1),p1(2),dp(1),dp(2), 0, 'MaxHeadSize',3/norm(dp),'color',col) % draw arrows between points 
                    % ^ 0 means no scaling, so will actually connect points. norm() finds length of vector dp.
                    % Defining MaxHeadSize is not quite working to get same-size arrowheads
            end
        end
    end
end
set(gcf,'Position',[10 10 500 500])
hold off
pause(1)

% plot back tracks
figure()
imshow(imMeanB, []); title('reference image: back side'); % show image to check
hold on
for i = myMovies
    for j = 1:size(allTracks(i).tracks, 2)
        if allTracks(i).tracks(j).LRB == 'b'
            col = rand(1,3); % generate random RGB color to use for this track
            for k = 1:length(allTracks(i).tracks(j).x) - 1 %loop through all intervals between frames of track  
                p1(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k); % x-coord of first point, time t
                p1(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k); % y-coord of first point
                p2(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k+1); % x-coord of 2nd point, time t+1
                p2(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k+1); % y-coord of 2nd point
                dp = p2-p1; % vector of delta_x, delta_y
                allTracks(i).tracks(j).dpAligned(k,:) = dp; % save velocity for later as vectors [dx, dy] (first column dx = v_x/dt)
                quiver(p1(1),p1(2),dp(1),dp(2), 0, 'MaxHeadSize',3/norm(dp),'color',col) % draw arrows between points 
                    % ^ 0 means no scaling, so will actually connect points. norm() finds length of vector dp.
                    % Defining MaxHeadSize is not quite working to get same-size arrowheads
            end
        end
    end
end
set(gcf,'Position',[10 10 500 500])
hold off
pause(1)


%% Optional: Plot bulk image from each movie separately, with all aligned tracks on top

cd(folder_name)
for i = myMovies
    % Load all frames from bulk-labeled movie.
    cd(num2str(i)); % move to directory for first movie

  imBulkMean = mean(imAlignedBulk(:,:,allTracks(i).tracks(1).count:allTracks(i).tracks(end).count),3); % average mean aligned cell images from each track of a given movie
  
  figure 
  imshow(imBulkMean, []);
  hold on
  for j = 1:size(allTracks(i).tracks, 2)
    trackLength = length(allTracks(i).tracks(j).x);
    col = 'm'; %rand(1,3); % generate random RGB color to use for this track
        for k = 1:trackLength - 1 %loop through all intervals between frames of track  
            p1(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k); % x-coord of first point, time t
            p1(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k); % y-coord of first point
            p2(1) = allTracks(i).tracks(j).alignedTracksYXF(2,k+1); % x-coord of 2nd point, time t+1
            p2(2) = allTracks(i).tracks(j).alignedTracksYXF(1,k+1); % y-coord of 2nd point
            dp = p2-p1; % vector of delta_x, delta_y
            quiver(p1(1),p1(2),dp(1),dp(2), 0, 'MaxHeadSize',3/norm(dp),'color',col) % draw arrows between points 
        end
  end
    set(gcf,'Position',[10 10 500 500])
    hold off
    pause(1)
    clear bulkMovie; clear imBulkMean; clear framesForBulk; clear count

 cd ../ %move back up to tracking directory in prep for next movie
end


%% 6. Calculate speeds

allInstSpeeds = []; 

for i = myMovies 
    for j = 1:size(allTracks(i).tracks, 2)
        for k = 1:size(allTracks(i).tracks(j).u, 1)
            allInstSpeeds = [allInstSpeeds allTracks(i).tracks(j).u(k,1)];
        end
    end
end


%% Plot histogram of speeds

figure()
h1=histogram(allInstSpeeds);
hold on
h1.FaceColor = 'b';
h1.BinWidth = 2;
%h1.Normalization = 'probability';

xlabel('speed (\mum/sec)')
ylabel ('number of events')
title('Histogram of Speeds for All Tracks')
set(gca,'FontSize',18)
xlim([0,25])


%% 7. Rose plot
% Plot polar histogram of angle (with respect to cell "head tail axis") - using dpAligned !!!!

% Note: set up for possible filtering by persistence or speed of tracks,
% but what's currently here and in manuscript is NO filtering, all tracks
minSpeed = 0; %minimum speed to include in histogram, in um/sec
minPersist = -1; %0.35;

allHighInstThetas = [];   
for i = myMovies 
    for j = 1:size(allTracks(i).tracks, 2)
        for k = 1:size(allTracks(i).tracks(j).dp, 1)
            if allTracks(i).tracks(j).u(k,1) > minSpeed && mean(allTracks(i).tracks(j).velCorr(2:end,2),1) > minPersist % remember this cut off is in real units of um/sec, if you've done right
                tempTheta = atan2( -1 * allTracks(i).tracks(j).dpAligned(k,2) , (allTracks(i).tracks(j).dpAligned(k,1)) ); % atan of -dy, dx. NOTE! flipping dy to match real space coordinate system
                allTracks(i).tracks(j).thetaAligned(k) = tempTheta; % save aligned theta for later
                allHighInstThetas = [allHighInstThetas tempTheta];
            end
        end
    end
end

figure
polarhistogram(allHighInstThetas, 14)
title(['Rose plot of angle of instantaneous speeds GREATER than ', num2str(minSpeed), ' um/s and minPersist ', num2str(minPersist)])



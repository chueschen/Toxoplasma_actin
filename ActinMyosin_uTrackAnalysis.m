%%%% Analysis of uTrack (Jaqaman et al., Nature Methods 2008) automatic tracking results for actin and myosin speckles %%%%
 

% 0. Add uTrack code to path and navigate to location of movies for analysis

cd('/Users/Christina/Documents/GitHub/Miscellany/u-track'); % navigate to directory containing uTrack code
addpath(genpath(pwd)); % add folder and subfolders to path

folder_name = uigetdir('/Volumes/');
cd(folder_name)

%% 1. Run uTrack (Danuser lab; Jaqaman et al., Nature Methods 2008) on each movie. 
% Each folder e.g. "Movie1" contains .tif movies of both the speckle/single molecule
% channel (SM) and the cell-body-labeling (bulk), as well as the uTrack outputs.
% We used release 2.2.0.

movieSelectorGUI % opens the uTrack graphical interface for uTrack analysis.
% Documentation for uTrack available here: https://github.com/DanuserLab/u-track?tab=readme-ov-file


%% 2. Import ALL TRACKS of tracking output from all movies and save info.

% Import uTrack tracking output and save info, displacement (dp) vectors, orientation of velocity (theta) 
% at each time step, directory location, frame #s, time step, track coordinates, etc.

%clearvars
singPix = 0.1; % pixel size in um
myMovies = 11; % number of movies tracked for this condition

for i = myMovies
    cd(folder_name)
    cd(num2str(i)); % move to directory for first movie
    loadw('*.mat'); % load MovieData
    allMovieData(i) = MD; %save MovieData for later
    cd('TrackingPackage/tracks'); % trying this cuz weird mistake loading wrong .mat
    load('Channel_1_tracking_result.mat'); %load tracking results
    
    % Pull out all tracks from this movie
    for j = 1:size(tracksFinal, 1)
        allTracks(i).tracks(j).tracksCoordAmpCG = tracksFinal(j).tracksCoordAmpCG; %save all tracks for later
            % 8 columns = x-coordinate, y-coordinate, z-coordinate (0 if 2D),amplitude, x-coordinate standard deviation, 
            % y- coordinate standard deviation, z-coordinate standard deviation (0 if 2D) and amplitude standard deviation.
        allTracks(i).tracks(j).frameStart = tracksFinal(j).seqOfEvents(1); %pull out frame number at which the track starts
    end
    
    % Pull out time interval between frames
    allTracks(i).frameInterval = allMovieData(1,i).timeInterval_; % in this case, time interval info was correctly saved during uTrack tracking.
    % If not, import metadata here or correct definition of frameInternval
    
    % Calculate dp vectors (instantaneous velocity) and theta (angle of dp from -pi to pi)
    for j = 1:size(allTracks(i).tracks, 2)
    trackLength = round(length(allTracks(i).tracks(j).tracksCoordAmpCG) ./ 8);
        for k = 1:trackLength - 1 %loop through all intervals between frames of track  
            % Calculate velocity vector for each dt
            p1(1) = allTracks(i).tracks(j).tracksCoordAmpCG(8*(k-1) + 1); % x-coord of first point, time t
            p1(2) = allTracks(i).tracks(j).tracksCoordAmpCG(8*(k-1) + 2); % y-coord of first point
            p2(1) = allTracks(i).tracks(j).tracksCoordAmpCG(8*(k) + 1); % x-coord of 2nd point, time t+1
            p2(2) = allTracks(i).tracks(j).tracksCoordAmpCG(8*(k) + 2); % y-coord of 2nd point
            allTracks(i).tracks(j).dp(k,:) = p2-p1; % save velocity for later as vectors [dx, dy] (first column dx = v_x/dt)
            % calculate speeds in microns/sec
            allTracks(i).tracks(j).u(k,:) = norm(allTracks(i).tracks(j).dp(k,:)) .* singPix ./ allTracks(i).frameInterval; %um/s
            % Calculate theta, orientation of velocity vector for each dt, in "standard" coordinates where positive y is up
            allTracks(i).tracks(j).theta(k,:) = atan2( -1 * allTracks(i).tracks(j).dp(k,2) , (allTracks(i).tracks(j).dp(k,1)) ); % atan of -dy, dx. NOTE! flipping dy to match real space coordinate system
        end
    end
    
    clear tracksFinal; clear MD; clear p1; clear p2; clear frameInterval; clear trackLength
end


%% 3. Calculate speeds

allInstSpeeds = []; 

for i = myMovies 
    for j = 1:size(allTracks(i).tracks, 2)
        for k = 1:size(allTracks(i).tracks(j).dp, 1)
            allInstSpeeds = [allInstSpeeds allTracks(i).tracks(j).u(k,1)];
        end
    end
end


%% 4. Plot histogram of speeds
figure()
h1=histogram(allInstSpeeds);
hold on
h1.FaceColor = 'm';
h1.BinWidth = 0.5;
% h1.Normalization = 'probability';

xlabel('Speed (\mum/sec)')
ylabel ('Probability')
title('Histogram of Speeds for All Tracks')
set(gca,'FontSize',18)


%% Comparing histograms of speeds
% "allInstSpeeds" renamed to myoSpeeds, actinSpeeds, and fixedSpeeds for each
% relevant condition.
figure
hold on

h3=histogram(myoSpeeds);
h3.FaceColor = 'c';
h3.BinWidth = 0.125;
h3.Normalization = 'probability';

h2=histogram(actinSpeeds);
hold on
h2.FaceColor = 'm';
h2.BinWidth = 0.125;
h2.Normalization = 'probability';

h1=histogram(fixedSpeeds);
h1.FaceColor = 'k';
h1.BinWidth = 0.125;
h1.Normalization = 'probability';

xlabel('Speed (\mum/sec)')
ylabel ('Probability')
title('Histogram of Speeds')
legend( 'Live Myosin', 'Live Actin', 'Fixed Actin')
set(gca,'FontSize',18)
xlim([0,10])


%% Plot CDF 

figure
hold on
cdfplot(fixedSpeeds)
cdfplot(myoSpeeds)
cdfplot(actinSpeeds)

xlabel('Speed (\mum/sec)')
ylabel ('CDF')
title('CDF of Speeds')
legend(  'Fixed Control', 'Myosin (MLC1)', 'Actin')
set(gca,'FontSize',18)
xlim([0,10])


%% Compute CDFs and compare with KS test

% first do this myself, with my own CDF with equal bin numbers

binWidth = 0.1;
edges = (0:binWidth:10); % set bins to go from 0 um/s to 10 um/s
binCenters = edges(2:end)-binWidth/2;

cdfFixed = histcounts(fixedSpeeds, edges, 'Normalization', 'cdf'); % compute empirical CDF with set number of bins
cdfMyo = histcounts(myoSpeeds, edges, 'Normalization', 'cdf'); 
cdfActin = histcounts(actinSpeeds, edges, 'Normalization', 'cdf'); 

figure
hold on
plot(binCenters, cdfFixed)
plot(binCenters, cdfMyo)
plot(binCenters, cdfActin)
legend('fixed','myosin','actin')

myKSfixedVsMyo = max(abs(cdfFixed-cdfMyo));
myKSfixedVsActin = max(abs(cdfFixed-cdfActin));
myKSmyoVsActin = max(abs(cdfMyo-cdfActin));

% second, matlab's build-in two-sample K-S test on raw speeds:
[ksFixedVsMyo, ksFixedVsMyoP, ksFixedVsMyoD] = kstest2(fixedSpeeds, myoSpeeds, 'Alpha', 0.05);
[ksFixedVsActin, ksFixedVsActinP, ksFixedVsActinD] = kstest2(fixedSpeeds, actinSpeeds, 'Alpha', 0.05);
[ksMyoVsActin, ksMyoVsActinP, ksMyoVsActinD] = kstest2(myoSpeeds, actinSpeeds, 'Alpha', 0.05);





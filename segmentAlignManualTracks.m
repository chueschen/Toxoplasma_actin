function [allTracks, imAlignedBulk] = segmentAlignManualTracks(allTracks, myMovie, imAlignedBulk)

% Read in all frames of bulk movie and define polarity.
% Choose tail-to-head vector and use that to rotate image. Segment cell and use center of mass
% to translate image.

for i = myMovie 
    
% Load all frames from bulk-labeled movie.
    cd(num2str(i)); % move to directory for first movie
    for j = 1:size(allTracks(i).tracks, 2) 
        
        clear imBulkMean; clear imFilt; clear imNorm; clear imThresh; clear imClear; clear imLab; clear props; clear rotatedIm; clear translatedIm
    
        disp(['movie'  num2str(i)  ', track'  num2str(j)])
        framesForBulk = (allTracks(i).tracks(j).f(1):1:allTracks(i).tracks(j).f(end));
        for k = framesForBulk
            bulkMovie(:,:,k) = imread(getfield(dir('C*bulk.tif'),'name'), k ); %import bulk movie frames
        end
        imBulkMean = mean(bulkMovie,3); % average all frames of bulk movie to get "average" image for polarity choice
        imBulkMean = mat2gray(imBulkMean); % NOTE! Normalizing here so max intensity = 1, so imAlignedBulk images will all be scaled.
        %%%%%%%% next gen: can segment cell at every frame and use totranslate/rotate each track segment into cell frame %%%%%%%% 
        clear bulkMovie
        figure(10); imshow(imBulkMean, []); title('imBulkMean OG'); % display "average" cell image

    %%%%%%%%% Segment cell and fit ellipse.%%%%%%%%%%%%

        % Apply a median filter to reduce the noise in the image. 
        imFilt = medfilt2(imBulkMean); %sets each pixel to median pixel of the 3x3 box centered there
        imNorm = mat2gray(imFilt); %scales values (intensity) so that min = 0 and max = 1 %NOTE - duplicate for now since I'm trying out normalization above.

        figure; histogram(imNorm);
    %     check = 'n';
    %     while (check == 'n')

        thresh = graythresh(imNorm); % use graythresh to find a good threshold 
            % "The graythresh function uses Otsu's method, which chooses the threshold
            % to minimize the intraclass variance of the black and white pixels"

        imThresh = imNorm > thresh; %keep only pixels above threshold

        % Clear the border and label. 
            imClear = imclearborder(imThresh); %supresses structures lighter than surroundings and touching a border
            imLab = bwlabel(imClear); %label connected components in 2D binary image

        % Apply an area threshold. 
            props = regionprops(imLab, 'Area'); %uses regionprops to find area of a binary array, imLab
            %             tempAreas(:) = props.Area;
            approvedObj = zeros(size(imBulkMean));
            for k=1:size(props,1)
                tempAreas(k,1) = props(k).Area;
            end

        %instead of throwing out things below a certain area, I'll just
        %keep the biggest single object
            tempIndex = find(tempAreas == max(tempAreas)); %find index of biggest object
            approvedObj = approvedObj + (imLab==tempIndex);
            imLab = bwlabel(approvedObj); %relabel only these approved objects

            figure(11); imshow(imLab)
    %         check = input('Are you happy with this segmentation of the CELL? y/n','s');


        props = regionprops(imLab, 'Area', 'Centroid','Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity'); %use regionprops 
        cellCentroidX(i) = props(1).Centroid(1); %pulls centroids into separate matrix
        cellCentroidY(i) = props(1).Centroid(2); %pulls centroids into separate matrix
        cellCentroidInt(2,i) = round(cellCentroidX(i),0);
        cellCentroidInt(1,i) = round(cellCentroidY(i),0); % Y-coords will be first row, X-coords 2nd row %%%%!!!! so coords are i,j
        orient(i) = props.Orientation; %pulls into separate matrix . in DEGREES.

        figure(10) % show bulk sum image to check where nucleus is. Note: before starting this analysis, 
        % carefully annotate cell polarity using microtubule imaging or
        % nearby-in-time twirling events. Confirm that nuclear position and
        % this second polarity indicator agree.
        topOrBottom = input('Is the cell pointing apical end up? y/n','s');
        if topOrBottom == 'y' && orient(i) < 0
            polarity(i) = orient(i) + 180; % Get polarity vector into correct quadrant, tail to apical end
        elseif topOrBottom == 'n' && orient(i) >= 0
            polarity(i) = orient(i) - 180;
        else
            polarity(i) = orient(i); % in degrees, since that's what imrotate takes
        end

    %%%%%%%%% Rotate and translate cell image to align %%%%%%%%%%%%
    alpha = 90 - polarity(i);   % angle for rotation to get apical end up
    rotatedIm = imrotate(imBulkMean, alpha);   % rotation of the image
    figure(14); imshow(rotatedIm, []); title('rotatedIm');
    leftRightBack = input('Is the cell on its left side, right side, or back? l/r/b','s');
    rotMatrixTemp = [cosd(alpha) -sind(alpha); sind(alpha) cosd(alpha)]; 
    imCenterA = (size(imBulkMean)/2)';         % Center of the main image. First row is i=y, second row is size of j=x
    imCenterB = (size(rotatedIm )/2)';  % Center of the transformed image
    rotatedCentroid = rotMatrixTemp*(cellCentroidInt(1:2,i)-imCenterA)+imCenterB; % new coordinates of centroid after rotation
    allTracks(i).tracks(j).LRB = leftRightBack;
    allTracks(i).tracks(j).rotMatrix = rotMatrixTemp;
    allTracks(i).tracks(j).centerImOG = imCenterA;
    allTracks(i).tracks(j).centerRotIm = imCenterB;

    % Now use centroid to translate and align to center of mass.
    % Assuming movies of around 100 x 100 pixels, align centroids to position 50, 50
    translationVector = [round(50-rotatedCentroid(2,1))  round(50-rotatedCentroid(1,1))]; % translation vector (x,y) = (j,i)
    translatedIm = imtranslate(rotatedIm, translationVector, 'OutputView','full'); 
    figure(12); imshow(translatedIm, [])
    allTracks(i).tracks(j).translVector(1,1) = translationVector(2); %saving translation vector in form [ i ; j ] for compatibility
    allTracks(i).tracks(j).translVector(2,1) = translationVector(1);

    % save, and crop so that 50 is center and all images will be same size
    translateCropBox = -1 .* translationVector;
    translateCropBox = max(translateCropBox, 0); %set to zero unless original value was negative
    startCropI = 11+translateCropBox(2); allTracks(i).tracks(j).cropIJ(1,1) = startCropI; % save range of crop box for later
    endCropI = 89+translateCropBox(2); allTracks(i).tracks(j).cropIJ(2,1) = endCropI; % [ startCropI/Y   startCropJ/X ;
    startCropJ = 11+translateCropBox(1); allTracks(i).tracks(j).cropIJ(1,2) = startCropJ;%  endCropI/Y   endCropJ/X ]
    endCropJ = 89+translateCropBox(1); allTracks(i).tracks(j).cropIJ(2,2) = endCropJ;
    
    % Figure out what number this track is, of all total tracks (not i - that's movies!)
    count = 0;
    for l = 1:i-1   
       count = count + length(allTracks(l).tracks); %count all tracks in preceding movies
    end
    count = count + j; % now add tracks to date from this movie
    allTracks(i).tracks(j).count = count; % save this index for later
    
    imAlignedBulk(:,:,count) = translatedIm(startCropI:endCropI, startCropJ:endCropJ); % save final image in imAlignedBulk
    
   end

cd ../ %move back up to tracking directory in prep for next movie

end
% Outputs for when I make this a function: imAlignedBulk, allTracks

end
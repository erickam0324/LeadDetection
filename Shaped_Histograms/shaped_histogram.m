% Function to transform an object in the "foreground" of an image into another image
% that has a histogram that looks like the foreground itself.
% I provide three examples with 3 demo images:
%   I change the image of a city skyline so that the histogram of the new image looks like the skyline.
%   I change the image of a car so that the histogram of the new image looks like the car shape.
%   I change the image of a woman so that the histogram of the new image looks like the shape of the woman.
% By Image Analyst
% September 2010
function shaped_histogram()
try
	% Change the current folder to the folder of this m-file.
	% (The line of code below is from Brett Shoelson of The Mathworks.)
	if(~isdeployed)
		cd(fileparts(which(mfilename)));
	end
	
	clc;	% Clear command window.
	clear;	% Delete all variables.
	close all;	% Close all figure windows except those created by imtool.
	imtool close all;	% Close all figure windows created by imtool.
	workspace;	% Make sure the workspace panel is showing.
	fontSize = 14;
	
	% Introduce the demo, and ask user if they want to continue or exit.
	message = sprintf('This demo will create an image with that has a histogram shaped like the object in the image, or with a flat histogram (TRULY histogram equalized).\nIt requires the Image Processing Toolbox.\nMake your selection.');
	reply = questdlg(message, 'Run Demo?', 'Shaped Histogram', 'Flat Histogram', 'Cancel', 'Shaped Histogram');
	if strcmpi(reply, 'Cancel')
		% User canceled so exit.
		return;
	elseif strcmpi(reply, 'Shaped Histogram')
		useFlatHistogram = false;
	else
		% Flat histogram.
		useFlatHistogram = true;
	end

	% Check that user has the Image Processing Toolbox installed.
	versionInfo = ver; % Capture their toolboxes in the variable.
	hasIPT = false;
	for k = 1:length(versionInfo)
		if strcmpi(versionInfo(k).Name, 'Image Processing Toolbox') > 0
			hasIPT = true;
		end
	end
	if ~hasIPT
		% User does not have the toolbox installed.
		message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
		reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
		if strcmpi(reply, 'No')
			% User said No, so exit.
			return;
		end
	end

	message = sprintf('Which demo image do you want to use?');
	selectedImage = questdlg(message, 'Which Demo Image?', 'Chicago Skyline', 'Ford Shelby Car', 'Woman on Beach', 'Chicago Skyline');
	if strcmp(selectedImage, 'Chicago Skyline')
	 	baseFileName = 'chicago_skyline_small.png';
		demoImage = 1;
	elseif strcmp(selectedImage, 'Ford Shelby Car')
		baseFileName = '2011_Ford_Shelby_GT500.png';
		demoImage = 2;
	else
		baseFileName = 'Beach_Woman.png';
		demoImage = 3;
	end
	
	% Read in demo image.
	rgbImage = imread(baseFileName);
	[rows columns numberOfColorBands] = size(rgbImage);
	f1 = figure;
	subplot(2, 3, 1);
	imshow(rgbImage, [0 255]);
	title('Original Color Image', 'FontSize', fontSize, 'Color', 'b');
	set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
	set(gcf, 'name', 'Demo (page 1 of 2) by ImageAnalyst', 'numbertitle', 'off') 
	drawnow; % Force it to display immediately.

	% Plot its histograms
	PlotRGBHistogram(rgbImage);

	% Get a gray version of the image that we will histogram and transform to get a new, shaped histogram.
	% Get a binary version of the image where the object along the bottom part of the photo (the "foreground") is black.
	% IMPORTANT NOTE:  Each image has a different background and foreground, so they must each have a
	% customized way of extracting the foreground.  Even though they all use thresholding,
	% they do it on different color space channels and with different threshold values.
	switch demoImage
		case 1
			% This image uses segmentation of the Value channel in HSV color space.
			[grayImage, binaryImage] = GetSkylineBinaryImage(rgbImage);
		case 2
			% This image uses segmentation of the Red channel in RGB color space.
			[grayImage binaryImage] = GetCarBinaryImage(rgbImage);
		otherwise
			% This image uses segmentation of the Hue channel in HSV color space.
			[grayImage binaryImage] = GetWomanBinaryImage(rgbImage);
	end
	% Display the initial gray version of the image.
	figure(f1);	% Switch to the first figure (in case any of the above routines created an extra figure).
	subplot(2, 3, 3);
	imshow(grayImage, []);
	caption = sprintf('Initial Gray Image (used to derive\nbinary image below).');
	title(caption, 'FontSize', fontSize, 'Color', 'b');
	
	% Display the binary image.
	subplot(2, 3, 6);
	imshow(binaryImage, []);
	caption = sprintf('Binary Image of the "Foreground" Shape\n(This is an image, not a histogram.)');
	title(caption, 'FontSize', fontSize, 'Color', 'b');
	axis square;

	% Now, let's get its histogram.
	[pixelCount grayLevels] = imhist(grayImage);
	subplot(2, 3, 2); 
	bar(pixelCount, 'BarWidth', 1); 
	caption = sprintf('Histogram\nof Initial Gray Image');
	title(caption, 'FontSize', fontSize, 'Color', 'b');
	xlim([0 grayLevels(end)]); % Scale x axis manually.
	drawnow; % Force it to display immediately.

	% Now let's transform the gray image so that it has a histogram like the binary image.
	[outputImage desiredPixelsInEachBin] = TransformImage(grayImage, binaryImage, useFlatHistogram);
	outputImage = uint16(outputImage);

	% Now, let's get the histogram of the transformed image.
	% It should look like the skyline.
	numberOfOutputBins = columns;
	% Get the histogram, one bin per gray level.
	[pixelCountT grayLevelsT] = imhist(outputImage, 65536);
	% We don't need all those gray levels, just the same number of gray levels
	% as there are columns in the binary image.  So take just those.
% 	pixelCountT(1) = 0;
	pixelCountT = pixelCountT(1:numberOfOutputBins);
	
	% If you want to do an accuracy check, you can subtract the actual
	% histogram from the desired one.
	misMatches = pixelCountT-desiredPixelsInEachBin';
	% This array should have 0 in every location.
	
	% Display the images and plots larger.
	% Display the original gray image in the upper left.
	figure;
	subplot(2, 2, 1);
	imshow(grayImage, []);
	title('Initial Gray Image', 'FontSize', fontSize, 'Color', 'b');
	% Display the transformed image in the upper right.
	subplot(2, 2, 2);
	imshow(outputImage, []);
	title('Transformed Gray Image', 'FontSize', fontSize, 'Color', 'b');
	axis on;
	% Display the original histogram in the lower left.
	subplot(2, 2, 3);
	bar(pixelCount, 'BarWidth', 1); 
	grid on;
% 	axis square;
	caption = sprintf('Histogram of\nInitial Gray Image');
	title(caption, 'FontSize', fontSize, 'Color', 'b');
	% Display the transformed histogram in the lower right.
	subplot(2, 2, 4);
	bar(pixelCountT, 'BarWidth', 1); 
% 	axis square;
	grid on;
	caption = sprintf('Histogram of\nTransformed Gray Image');
	title(caption, 'FontSize', fontSize, 'Color', 'b');
	xlim([0 numberOfOutputBins]); % Scale x axis manually.
	% If it's a flat histogram, let's make the max Y height a little bigger than the height of the bars.
	% Otherwise, for some histograms (like the woman demo image), it has the bars take up 100% of the y axis and it just looks
	% like a big black block.  Extend the axis a little bit so we can actually see the axis and the top of the histogram.
	if useFlatHistogram
		maxBarheight = max(pixelCountT);
		yAxisLength = maxBarheight * 1.2;
		ylim([0 yAxisLength]);
	else
		% Adjust the height of the bar chart to match the aspect ratio of the binary image.
		newYMax = AdjustBarChartHeight(binaryImage, pixelCountT);
		ylim([0 newYMax]);
	end
	
	% Enlarge figure to full screen.
	set(gcf, 'Position', get(0,'Screensize')); 
	set(gcf, 'name', 'Demo (page 2 of 2) by ImageAnalyst', 'numbertitle', 'off');

	% Ask if user wants to save the transformed image to a new file.
	message = sprintf('Done with the demo!\nCheck out both of the figure windows that were created.\n\nDo you want to save the transformed image to a new file?');
	reply = questdlg(message, 'Done.  Save Image?', 'Save Image', 'Exit', 'Exit');
	if strcmpi(reply, 'Save Image')
		filterSpec = '*.*';
		defaultName = sprintf('UniqueHistogramImage.png');
		dialogTitle = 'Save Tranformed Image?';
		[baseFileName, folder] = uiputfile(filterSpec, dialogTitle, defaultName);
		if baseFileName == 0
			% User clicked cancel.
			return; % Bail out.
		end
		fullFileName = fullfile(folder, baseFileName);
		% If the number of bins is less than 256, convert it to a uint8 image
		% for ease of getting histograms in other software packages.
		% Otherwise it will be a uint16 image and the gray levels will be in the low end of the 32768 gray levels,
		% so other software may compress the histogram and gray levels (because they want to do the whole range)
		% and this will make it difficult to see unless you do special things to see them.
		if numberOfOutputBins <= 256
			% Our demo images are NOT this because they all have more than 256 columns
			% so we won't get here for our demos.
			outputImage = uint8(outputImage);
		end
		% Save the array out to a disk file.
		imwrite(outputImage, fullFileName);
		message = sprintf('Image file has been saved.\n%s\n\nNOTE: To make this 16 bit image look right in Photoshop\nyou have to do Image->Adjustments->Auto Levels', fullFileName);
		msgbox(message);
	end

catch ME
	errorMessage = sprintf('Error in shaped_histogram():\n\nError Message:\n%s', ME.message);
	uiwait(warndlg(errorMessage));
end
return; % from shaped_histogram, the main routine.

%=====================================================================
% Get a binary version of the skyline, where the skyline at the bottom of the image is black.
% Of course, this algorithm may need to be modified for individual images.
% I know this algorithm works fairly well for the Chicago skyline demo image I supplied.
function [grayImage binaryImage] = GetSkylineBinaryImage(rgbImage)
	try
		% Convert to floating point so it does the calculations correctly.
		% Also needs to be normalized to 0-1.
		rgbFloating = double(rgbImage) / 255.0;

		% Compute hsv image
		hsvImage = rgb2hsv(rgbFloating);
		% V image:
		vImage = hsvImage(:,:,3);
		% Binarize it.
		binaryImage = vImage > 0.7;
		[rows columns] = size(binaryImage);
		% Extract the skyline.
		% Convert to columns, filling in any holes between the top of the silhouette and the bottom of the image.
		for col = 1 : columns
			oneColumn = binaryImage(:, col);
			startingRow = find(oneColumn==0, 1, 'first');
			binaryImage(startingRow:end, col) = 0;
		end
		% Prepare which image we're going to supply as the monochrome image to use in transforming and histogramming.
		grayImage = rgb2gray(rgbImage);
	catch ME
		errorMessage = sprintf('Error in GetSkylineBinaryImage():\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	return; % from GetSkylineBinaryImage

%=====================================================================
% Get a binary version of the skyline, where the skyline at the bottom of the image is black.
% Of course, this algorithm may need to be modified for individual images.
% I know this algorithm works fairly well for the Ford car demo image I supplied.
function [grayImage binaryImage] = GetCarBinaryImage(rgbImage)
	try
		% Prepare which image we're going to supply as the monochrome image to use in transforming and histogramming.
		grayImage = rgbImage(:, :, 1); % Extract the red channel.
		binaryImage = grayImage < 50;  % Red less than 50 is the background.
		[rows columns] = size(binaryImage);
		% Convert to columns, filling in any holes between the top of the silhouette and the bottom of the image.
		for col = 1 : columns
			oneColumn = binaryImage(:, col);
			startingRow = find(oneColumn==0, 1, 'first');
			binaryImage(startingRow:end, col) = 0;
		end
	catch ME
		errorMessage = sprintf('Error in GetSkylineBinaryImage():\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	return; % from GetSkylineBinaryImage

%=====================================================================
% Get a binary version of the skyline, where the skyline at the bottom of the image is black.
% Of course, this algorithm may need to be modified for individual images.
% I know this algorithm works fairly well for the woman demo image I supplied.
function [grayImage binaryImage] = GetWomanBinaryImage(rgbImage)
	try
		hsvImage = rgb2hsv(rgbImage);
		% Compute hsv image
	% 	hsvImage = rgb2hsv(rgbFloating);
		hsvImage = rgb2hsv(rgbImage);
		% H image:
		hImage = hsvImage(:,:, 1);
		% S image:
		sImage = hsvImage(:,:, 2);
		% V image:
		vImage = hsvImage(:,:, 3);
		clear('hsvImage');
		
		% Get the thresholds for Hue, Saturation, and value channels.
		lowThresholdH = 0.0;
		highThresholdH = 0.11;
		lowThresholdS = 0.16;
		highThresholdS = 0.75;
		lowThresholdV = 0.0;
		highThresholdV = 0.83;
		
		% Ask user to confirm.
% 		[lowThresholdH highThresholdH] = threshold(lowThresholdH, highThresholdH, hImage);
% 		[lowThresholdS highThresholdS] = threshold(lowThresholdS, highThresholdS, sImage);
% 		[lowThresholdV highThresholdV] = threshold(lowThresholdV, highThresholdV, vImage);

		% Combine the binary images from each channel to get an overall binary image.
		% Find the water and sky
		binaryImageH = (hImage > lowThresholdH) & (hImage < highThresholdH);
		binaryImageS = (sImage > lowThresholdS) & (sImage < highThresholdS);
		binaryImageV = (vImage > lowThresholdV) & (vImage < highThresholdV);
		% Actually it looks pretty good with using the hue channel ONLY, so let's just do that.
		binaryImage = binaryImageH; % & binaryImageS & binaryImageV;
		% Ask user to confirm.
% 		sumImage = binaryImageH + binaryImageS + binaryImageV;
% 		[lowThresholdV highThresholdV] = threshold(0, 2, sumImage);
		[rows columns] = size(binaryImage);
		
		% Get rid of small specks, anything smaller than 10% of the image.
		smallestAllowableObject = 0.1 * rows * columns;
		binaryImage = bwareaopen(binaryImage, smallestAllowableObject, 4);

		% Optional: Display the HSV images and binary images.
		showResults = true;  % Set = true if you want to see these images, false if you don't.
		if showResults
			% This is helpful in determining the threshold values to use.
			figure;
			set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
			set(gcf, 'name', 'HSV images of Woman', 'numbertitle', 'off') 
			subplot(3, 3, 1);
			imshow(rgbImage, []);
			title('RGB Image', 'FontSize', 20);
			subplot(3, 3, 2);
			imshow(binaryImage, []);
			title('Binary Image', 'FontSize', 20);
			axis on;
			% Show HSV binarized images.
			subplot(3, 3, 4);
			imshow(hImage, []);
			colorbar;
			title('H Image', 'FontSize', 20);
			axis on;
			subplot(3, 3, 5);
			imshow(sImage, []);
			colorbar;
			title('S Image', 'FontSize', 20);
			axis on;
			subplot(3, 3, 6);
			imshow(vImage, []);
			colorbar;
			title('V Image', 'FontSize', 20);
			axis on;
			% Show HSV binarized images.
			subplot(3, 3, 7);
			imshow(binaryImageH, []);
			title('H Image Binarized', 'FontSize', 20);
			axis on;
			subplot(3, 3, 8);
			imshow(binaryImageS, []);
			title('S Image Binarized', 'FontSize', 20);
			axis on;
			subplot(3, 3, 9);
			imshow(binaryImageV, []);
			title('V Image Binarized', 'FontSize', 20);
			axis on;
		end
		% binaryImage is now the woman.  But we want the woman to be black so that it looks like her silhouette.
		% So invert it.
		binaryImage = ~binaryImage;
		% Convert to columns, filling in any holes between the top of the silhouette and the bottom of the image.
		for col = 1 : columns
			oneColumn = binaryImage(:, col);
			startingRow = find(oneColumn==0, 1, 'first');
			binaryImage(startingRow:end, col) = 0;
		end
		% Show the final result.
		subplot(3, 3, 3);
		imshow(binaryImage, []);
		title('Final Binary Image', 'FontSize', 20);
		% Prepare which image we're going to supply as the monochrome image to use in transforming and histogramming.
		grayImage = rgb2gray(rgbImage);
	catch ME
		errorMessage = sprintf('Error in GetSkylineBinaryImage():\n\nError Message:\n%s', ME.message);
		fprintf(1, '\n%s\n', errorMessage);
		uiwait(warndlg(errorMessage));
	end
	return; % from GetSkylineBinaryImage

%=====================================================================
function PlotRGBHistogram(rgbImage)
	try
		fontSize = 20;
		redPlane = rgbImage(:, :, 1);
		greenPlane = rgbImage(:, :, 2);
		bluePlane = rgbImage(:, :, 3);

		% Let's get its histograms.
		[pixelCountR grayLevelsR] = imhist(redPlane);
		subplot(2, 3, 4:5); 
		plot(pixelCountR, 'r', 'LineWidth', 2); 
		caption = sprintf('Histograms of\nRed, Green, & Blue Color Channels');
		title(caption, 'FontSize', fontSize);
		xlim([0 grayLevelsR(end)]); % Scale x axis manually.

		[pixelCountG grayLevelsG] = imhist(greenPlane);
		hold on;
		plot(pixelCountG, 'g', 'LineWidth', 2); 

		[pixelCountB grayLevelsB] = imhist(bluePlane);
		xlim([0 grayLevelsB(end)]); % Scale x axis manually.
		plot(pixelCountB, 'b', 'LineWidth', 2); 
	catch ME
		errorMessage = sprintf('Error in PlotRGBHistogram():\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	return; % from PlotRGBHistogram()
	
	
%=====================================================================
function [outputImage desiredPixelsInEachBin] = TransformImage(grayImage, binaryImage, useFlatHistogram)
	try
	% Find the dimensions the image.
	[rows columns] = size(binaryImage);
	numberOfPixelsInImage = double(rows) * double(columns);
	
	% A gray level image has only 256 gray levels, but the binary image
	% can have more of less than that number of columns in it.
	% To do this accurately, we'll need to convert the uint8 image so that
	% it will have the same number of intensity levels as the binary image has columns.
	dblImage = double(grayImage);
	% I know this sounds strange but you need some randomness to do this,
	% otherwise, how can you split up the number of pixels at one gray level
	% (and thus all lie in one bin) into a different number of bins?
	% You may have to think about that a bit, but if you do, you'll see that you
	% can do it two ways.  Let's take an example.  Let's say that you have 1000 pixels
	% in bin 50, and let's say the binary image is 512 columns across.
	% So we need to take the pixels that are in bin 50 (which may be scattered
	% across the image) and put them into two bins, which may be at locations
	% 71 and 72, for example.  So which of those 1000 pixels gets mapped to gray level
	% 71 and which of the 1000 pixels gets mapped to gray level 72?  See - that's why
	% you need randomness.  One way to do it is to just find the 1000 pixels 
	% in the image and randomly take right proportion of them to put in each bin.
	% The way I'm going to do it is to add up to +/- 0.5 gray levels of noise to the pixel values of our
	% floating point image.  This is valid because we have a quantized image.  For example
	% take a pixel at gray level 50.  Now you know that it could have had a
	% number of photons coming in that corresponded to a gray level anywhere in between
	% 49.5 gray levels and 50.5 gray levels.  Anything in that range would have been
	% quantized to 50.  You wouldn't be able to tell.  So we're basically just undoing the process.
	% We're taking everything at 50 and redistributing it back over 49.5 to 50.5.
	% Then, when we sort the image, and take some fractions of the pixels with gray level 50,
	% and remap them to another gray level (e.g. 71 or 72), they will come from random locations all over the image.
	% Pretty clever, huh?
	noisyImage = dblImage + rand(rows, columns) - 0.5;
	minGL = min(min(noisyImage));
	maxGL = max(max(noisyImage));
	slope = 255.0 / (maxGL - minGL);
	noisyImage = slope * (noisyImage - minGL);
	% Now there may be some pixels in the -0.5 - 0 range.
	% We don't want that to mess up our bin positions so let's put all negative
	% values in the 0-0.5 bin by simply negating them.
	noisyImage(noisyImage<0) = -noisyImage(noisyImage<0);
	% Same thing for pixels in the 255 - 255.5 range.
	% Move them back into the 254.5 - 255 range by simply subtracting 0.5.
	noisyImage(noisyImage>255) = noisyImage(noisyImage>255) - 0.5;
	% Check to make sure it did it right.
	minGL = min(min(noisyImage));
	maxGL = max(max(noisyImage));
	
	% Let's get its histogram.
	% Have the number of bins be equal to the number of columns.
	% Find out the counts that each bin ACTUALLY holds.
	[pixelCountsActual grayLevelsActual] = imhist(noisyImage, columns);
	numberOfBins = length(pixelCountsActual);
	
	% Find out the percentage that each bin SHOULD hold.
	binPercents = GetBinPercentages(binaryImage);
	% Find out the number of pixels that each bin SHOULD hold.
	desiredPixelsInEachBin = round(binPercents * numberOfPixelsInImage);
	% If they wanted a flat histogram, change it to a flat histogram.
	if useFlatHistogram
		binPercents = ones(1, numberOfBins) / numberOfBins;
		pixelsPerBin = round(numel(grayImage) / numberOfBins);
		desiredPixelsInEachBin(:) = pixelsPerBin;
	end
	
	lineImage = reshape(noisyImage, [1, rows*columns]);
	% Sort this
	[sortedValues, sortedIndexes] = sort(lineImage);
	% Now here's the tough, tricky part.
	% We need to go through every intensity of the output image (0 - #columns in image), and:
	% 1. Find the number of the pixels of that intensity that should go in that bin.
	% 2. Start going along the sorted input image and find out where those gray levels occur
	% 3. Assign those locations in the output image the intensity that we're processing.
	% 4. Then move on to the next intensity.
	% By the way, if they asked for a flat histogram, this will give a much flatter histogram than 
	% the function adapthisteq(), in fact, it will give a perfectly flat histogram.
	outputImage = -1 * ones(1, length(lineImage));
	% Make a pointer to our sorted list of input pixels.
	% We're going to march this along until we've transformed
	% each and every one of them and stored the new value in our output image.
	pointer = 0;
	maxOutputGrayLevel = length(grayLevelsActual);
	for outputBin = 1 : maxOutputGrayLevel
% 		fprintf(1, 'outputBin = %d,   pointer = %d\n', outputBin, pointer);
		desiredPixelsInThisBin = desiredPixelsInEachBin(outputBin);
		% Start marching along all input pixels...
		for p = 1:desiredPixelsInThisBin
			if pointer+p > numberOfPixelsInImage
				% Bail out if rounding errors cause this to be greater
				% than the array length.
				break;
			end
% 			originalGrayLevel = sortedValues(pointer+p);
			% Find out the linear index where this sorted pixel originally lived.
			originalIndex = sortedIndexes(pointer+p);
			% Replace the pixel at that original location
			% with the gray level that we need it to have.
			outputImage(originalIndex) = outputBin - 1; % New, desired gray level
		end
		% Move the pointer along to the next input pixel that we will transform.
		pointer = pointer + desiredPixelsInThisBin;
		% Now go on to the next gray level, if we need to, by continuing the loop.
	end
	
	% Sometimes, due to rounding, there are still a very few pixels (~10 or so)
	% that did not get assigned.  These will have the value -1 and should be
	% assigned the maximum gray level.  Unless fixed, they show up a dark
	% anomolous dark specks in the bright parts of the image.  Let's fix them.
	outputImage(outputImage == -1) = maxOutputGrayLevel;
	
	% Now outputImage is good - it has been transformed.
	% But it's the wrong shape - it's still in a line.
	% Reshape it to the size of the original image.
	outputImage = reshape(outputImage, [rows columns]);
	
	catch ME
		errorMessage = sprintf('Error in TransformImage():\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	
	return; % from TransformImage
	
	
%=====================================================================
function binPercents = GetBinPercentages(binaryImage)
	try
		% Find the number of pixels in the image (and counts in the histogram).
		% Find the dimensions the image.
		[rows columns] = size(binaryImage);
		numberOfPixels = double(rows) * double(columns);
		% Find out how many black pixels there are in each column.
		blackPixelsPerColumn =  rows - sum(binaryImage, 1);
		numberOfBlackPixelsTotal = sum(blackPixelsPerColumn);
		binPercents = blackPixelsPerColumn ./ numberOfBlackPixelsTotal;
	
	catch ME
		errorMessage = sprintf('Error in GetBinPercentages():\n\nError Message:\n%s', ME.message);
		uiwait(warndlg(errorMessage));
	end
	return; % from GetBinPercentages
	
%=====================================================================
% Determines what the max Y height of the histogram bar chart should be
% such that the bar chart has the same aspect ratio as the binary image shape.
function newYMax = AdjustBarChartHeight(binaryImage, pixelCountT)
	% Get the size of the binary image.  We're interested in the height (rows).
	rows = size(binaryImage, 1);
	verticalRange = min(binaryImage, [], 2);  % Look for black pixels.
	% Get the tallest column.
	firstYValue = find(verticalRange, 1, 'last');
	imageAspectRatio = (rows - firstYValue + 1) / rows;
	% Get the max bar height
	maxBarHeight = max(pixelCountT);
	% Calculate the new ylim value.
	newYMax = maxBarHeight / imageAspectRatio;
	return; % from AdjustBarChartHeight()

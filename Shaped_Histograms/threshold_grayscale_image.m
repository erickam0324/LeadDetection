%  File and function name : threshold_grayscale_image
%  thresholded_image = threshold_grayscale_image(original_image,min_threshold,max_threshold)
%            (0 <= min_threshold, max_threshold <= 255)
%
%  Inputs        :   original_image  - original monochrome image (0 - black, 255 - white)
%                    min_threshold   - minimum value that should be recognized as a feature (inclusive).
%                    max_threshold   - maximum value that should be recoginized as a feature (inclusive).
%                    Inputs can be any data type.
%
%  Outputs       :   thresholded_image - binary thresholded image of type uint8
%                         Outside threshold = black(0),   Within threshold = white(1)
%
%  Description   :   If the pixel in the original image has a value between the min_threshold 
%                    and max_threshold (inclusive), then it will be assigned as a white pixel in the binary 
%                    image. Otherwise it will be a black pixel.
%
%  Sample Code to run:
%  monoImageArray = int16(imread('coins.png'));
%  thresholded_image = threshold_grayscale_image(original_image, 83, 255)
%            (0 <= min_threshold, max_threshold <= 255)

function thresholded_image = threshold_grayscale_image(original_image, min_threshold, max_threshold)
	% Initialize the value of white pixels.
	white_pixel = uint8(1);
	
	% Ensure that the threshold inputs are acceptable
	if min_threshold > max_threshold
		error('threshold_grayscale_image : min_threshold is greater then max_threshold');
	end
	
	%Initialize the image to be all black (non-feature)
	thresholded_image = uint8(zeros(size(original_image)));

	% You can only multiply integers if they are of the same type.
	% Convert the type of thresholded_image to the same data type as original_image.
% 	strDataType = class(original_image); % Get data type of original_image.
% 	thresholded_image = eval([strDataType '(thresholded_image)']); % Convert type of thresholded_image to same type as original_image.

	% Find all pixels which are in the threshold selection range.
	selected = (min_threshold <= original_image) & (original_image <= max_threshold);
	% Assign the pixel values, of the pixels in the threshold selection range, to white.
	thresholded_image(selected) = white_pixel;
	return;
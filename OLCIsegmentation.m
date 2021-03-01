function imout = OLCIsegmentation(inputimagename, outputname)
%% image segmentation for OLCI tryout from code by Computerphile on youtube 

    % This script is the one mentioned during the Computerphile Image
    % Segmentation video. I chose matlab because it's a popular tool for
    % quickly prototyping things. Matlab licenses are pricey, if you don't have
    % one (or, like me, work for an organisation that does) try Octave as a
    % good free alternative. This code should work in Octave too.

    % Load in an input image
    %im = imread('/Users/ericka/Desktop/Thesis/OLCIData/olci.png');
    im = imread(inputimagename);
    
    % In matlab, K-means operates on a 2D array, where each sample is one row,
    % and the features are the columns. We can use the reshape function to turn
    % the image into this format, where each pixel is one row, and R,G and B
    % are the columns. We are turning a W,H,3 image into W*H,3

    % We also cast to a double array, because K-means requires it in matlab
    imflat = double(reshape(im, size(im,1) * size(im,2), 3));

    % I specify that initialisation shuold sample points at
    % random, rather than anything complex like kmeans++ initialisation.
    % Kmeans++ takes a long time if you are using 256 classes.

    % Perform k-means. This function returns the class IDs assigned to each
    % pixel, and in this case we also want the mean values for each class -
    % what colour is each class. This can take a long time if the value for K
    % is large, I've used the sampling start strategy to speed things up.

    % While KMeans is running, it will show you the iteration count, and the
    % number of pixels that have changed class since last iteration. This
    % number should get lower and lower, as the means settle on appropriate
    % values. For large K, it's unlikely that we will ever reach zero movement
    % (convergence) within 150 iterations.
    K = 3; 
    [kIDs, kC] = kmeans(imflat, K, 'Display', 'iter', 'MaxIter', 150, 'Start', 'sample') ;

    % Matlab can output paletted images, that is, grayscale images where the
    % colours are stored in a separate array. This array is kC, and kIDs are
    % the grayscale indices.
    colormap = kC / 256 ; % Scale 0-1, since this is what matlab wants

    % Reshape kIDs back into the original image shape
    imout = reshape(uint8(kIDs), size(im,1), size(im,2)) ;

    % Save file out, you need to subtract 1 from the image classes, since once
    % stored in the file the values should go from 0 - 255, not 1 - 256 like matlab
    % would do.
    output = fullfile( '/Users/ericka/Desktop/Thesis/output', outputname + '.jpg') ; 
    imwrite(imout - 1, colormap, output) ;


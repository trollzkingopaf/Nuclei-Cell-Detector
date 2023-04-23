% Reading and displaying the nuclei image
nuclei = imread("StackNinja3.bmp");
figure;
imshow(nuclei);
title("Nuclei image");

% Brightening low-light areas
contrastAdjustment = imlocalbrighten(nuclei, "AlphaBlend", true);

% Convert RGB to CIE L*a*b* colour space
lab = rgb2lab(contrastAdjustment);

% Define thresholds for Lightness channel based on histogram settings 0 100
LMin = 0.000;
LMax = 100.000;

% Define thresholds for a* channel based on histogram settings
aMin = -128.000;
aMax = -1.000;

% Define thresholds for b* channel based on histogram settings
bMin = 0.000;
bMax = 127.000;

% Create binary mask based on chosen histogram thresholds
sliderBW = (lab(:, :, 1) >= LMin ) & (lab(:, :, 1) <= LMax) & ...
           (lab(:, :, 2) >= aMin ) & (lab(:, :, 2) <= aMax) & ...
           (lab(:, :, 3) >= bMin ) & (lab(:, :, 3) <= bMax);
BW = sliderBW;

% Initialize the image of green nuclei cells based on contrast adjusted nuclei image
greenCells= contrastAdjustment;

% Set green nuclei cells background pixels to zero and overlay onto the contrast adjusted nuclei image
greenCells(repmat(~BW, [1 1 3])) = 0;

% Extracting red, green and blue channels from the image of green nuclei cells
redCh = greenCells(:, :, 1);
greenCh = greenCells(:, :, 2);
blueCh = greenCells(:, :, 3);

% Computing greenness based on the formula of greenness
greenness = abs(greenCh - (redCh + blueCh) / 2);

% Bilateral filtering
greenness2 = imclearborder(greenness);
bilateral = imbilatfilt(greenness2);

% Defining the value for adaptive thresholding
bwThreshold = adaptthresh(bilateral) + 0.03;

% Binarization based on the value of adaptive thresholding
greennessBW = imbinarize(bilateral, bwThreshold);

% Thresholding to minimise neighbourhoods size < 4
thresholdCC = bwareafilt(greennessBW, [4 524288], 8);

% Filtering larger connected components
extractLarge = bwareafilt(thresholdCC, 10, "largest", 8);

% Disk structuring element with radius of 2
se2 = strel("disk", 2);

% Opening performed on larger connected components
openLarge = imopen(extractLarge, se2);

% Watershed Segmentation performed on larger connected components
% Complementing the distance transformation
largeDistance = -bwdist(~openLarge, "euclidean");

% Computing the watershed distance transformation
Ld = watershed(largeDistance);

% Segmenting using watershed ridge lines
watershedBW2 = openLarge;
watershedBW2(Ld == 0) = 0;

% Filter out tiny local minima to modify the distance transformation to minimise oversegmentation
watershedMask = imextendedmin(largeDistance, 2, 8);

% Modify the distance transformation so it only has minima at desired regions
largeDistance2 = imimposemin(largeDistance, watershedMask);
Ld2 = watershed(largeDistance2);
watershedBW3 = openLarge;
watershedBW3(Ld2 == 0) = 0;
fillwatershedBW3 = imfill(watershedBW3, "holes");

% Filtering smaller connected components
extractSmall = thresholdCC - extractLarge;

% Disk structuring element with radius of 1
se1 = strel("disk", 1);

% First iteration of opening performed on smaller connected components
openSmall = imopen(extractSmall, se1);
fillSmall = imfill(openSmall, "holes");

% Thresholding to remove disjoint connected components 6
thresholdSmall = bwareaopen(fillSmall, 15, 4);

% Second iteration of opening performed on smaller connected components
openSmall2 = imopen(thresholdSmall, se2);
fillSmall2 = imfill(openSmall2, "holes");

% Concatenate larger and smaller connected components
bwMorp = fillSmall2 + fillwatershedBW3;
imshow(bwMorp);
title("Finalised Morphology");

% Initialize image for the accuracy of green nucleis based on finalised morphological image.
greenMorp = contrastAdjustment;

% Set background pixels to 0.
greenMorp(repmat(~bwMorp, [1 1 3])) = 0;

% Labeling connected components in the binary image and returns a binary matrix containing labels from objects of the connected neighbourhood
labelNuclei = bwlabel(bwMorp, 4);

% Generates a random colour map with as many label rows in the matrix
numLabels = max(labelNuclei(:));
rng('shuffle');
colourmap = rand(numLabels, 3);

% Overlaying and displaying the colourisation result on top of the final binarization image
colourisation = labeloverlay(double(bwMorp), labelNuclei, "Colormap", colourmap);
figure;
imshow(colourisation);
title("Colourisation result");
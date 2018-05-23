function MyNewColorMap = jet_white
%Creates a colormap that follows the "jet" colormap, but sets zero
%values to white and linearly blends the white to the blue. This
%is equivalent to opening the colormapeditor, double-clicking the
%minimum color (darkest blue) and setting it to white in the popup
%palette.
%
%Since the "jet" colormap matrix can vary in size depending on the
%plot that it will be applied to, this function should be executed
%after any plot commands.
%
%Example:
%
%  grayImage = imread('tire.tif');
%  imshow(grayImage, []);
%  cMap = jet_white;
%  colormap(cMap)
%
% Author: Jim Kossin (NOAA's Center for Weather and Climate)
% Date: 28 October 2015
    
    MyNewColorMap = jet;
    iblue = find(ismember(MyNewColorMap, [0 0 1], 'rows'));
    yy = [1 0]; xx = [1 iblue];
    blend = interp1(xx, yy, 1:iblue);
    for i = 1:iblue
	MyNewColorMap(i,:) = [blend(i) blend(i) 1.0000];
    end
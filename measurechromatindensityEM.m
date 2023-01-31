function [data, FigH] = measurechromatindensityEM(tiffile)
% THis script measures the 'width' of the (dense) heterochromatin abutting
% the nuclear envelope.
% A first figure appears where the user is tasked to delimit the nuclear
% membrane.
%
% INPUT             -tiffile, string name of the image file (tif format)
%
%OUTPUT             -data, n-by-2 numerical array with the first column
%                   corresponding to the normalized location of th datapoint along the
%                   perimeter of the nucleus and the second column
%                   corresponds to the width of the dense heterochromatin
%                   (in pixels length)

%% Read and prepare file
I = imread(tiffile);
I = cast(I,'uint8');
Ic = imresize(I,0.2);

%% User defining nuclear envelope and creating mask
figure, imshow(Ic)
roi = drawfreehand(gca);
NUCMask = createMask(roi);
SMASK = regionprops(imdilate(NUCMask,ones(20)),'BoundingBox');
rect = SMASK.BoundingBox;
close;

%% Dense chromatin segmentation
Ica = imadjust(Ic);
Icac = imcomplement(Ica);
Icact = imtophat(Icac,strel('disk',25));
BW = imbinarize(Icact);
BWO5 = ordfilt2(BW,1,true(5));
BWA = bwareaopen(BWO5,1000);
f = 8; %define domain for second order minimum value
BWO = ordfilt2(BW,1,true(f));
RMask = imreconstruct(NUCMask,BWO) | NUCMask;
BWO(~RMask) = false;
BWOn = BWO & BWA;

%% Measure length of chromatin in a line orthogonal to the nuclear membrane
BWOHper = bwboundaries(RMask);
Sch = regionprops(RMask,'Area','Centroid','EquivDiameter');
diam = Sch.EquivDiameter / 5;
PixH = BWOHper{1}; PixH(end,:) = [];
np = size(PixH,1);
k = 1:10:np;
ThickInd = zeros(numel(k),2);
count = 1;
for i = 1:10:np %sampling every 10 pixels
    ThickInd(count,1) = i / np;
    [C, D] = perpcoord(PixH(1,:), PixH(10,:), diam);
    CD = [C; D];
    cline = improfile(BWOn,CD(:,2),CD(:,1));
    cline(isnan(cline)) = 0;
    [w, initcross] = pulsewidth(cline);
    
    if ~isempty(initcross)
        [m, idx] = min(abs(initcross-diam));
        if m < diam
            ThickInd(count,2) = w(idx);
        end
    end
    PixH = circshift(PixH,10,1);
    count = count+1;
end
data = ThickInd;

O = imoverlay(Ic,bwperim(BWOn),'r');
O = imcrop(O,rect);
FigH = figure; imshow(O);
%figure, area(data(:,1),data(:,2),'LineStyle','none');

end

function [C, D] = perpcoord(A, B, Len)

% A = [3 3]; %[x,y]
% B = [12,8]; %[x,y]
% Len = 3; % distance off the line AB that C will lie
% Call AB the vector that points in the direction
% from A to B.
AB = B - A;
% Normalize AB to have unit length
AB = AB / norm(AB);
% compute the perpendicular vector to the line
% because AB had unit norm, so will ABperp
ABperp = AB * [0 -1;1 0];
% midpoint between A and B
ABmid = (A + B) / 2;
% Compute new points C and D, each at a ditance
% Clen off the line. Note that since ABperp is
% a vector with unit eEuclidean norm, if I
% multiply it by Clen, then it has length Clen.
C = ABmid + Len * ABperp;
D = ABmid - Len * ABperp;
% plot them all
end


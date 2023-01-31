
function CW=drawCircle2(x, y,Centroids,r)
% drawCircle2 creates a binary image of size y-by-x containing full circles
% defined by Centroids and radii r
% INPUT:    -x -y size of MASK image
%           -Centroids numcentroids*2 vector or 2-by-numCentroids matrix
%           containing xy coordinates of circle centroids
%           -r scalar or vector same size of numcentroids defining radius
%           of circles
%
% OUTPUT:   -CW binary image containing circles described by centroids and
%           radii

[X, Y]=meshgrid(1:x,1:y);

BWgrid=false(y,x);
nc=numel(Centroids)/2;
count=1;
if numel(r)==1
    r=repmat(r,nc,1);
elseif numel(r)~=nc
    error('r must be a scalar or the same size of centroid points.')
end
if any(size(Centroids)==1)
    for i=1:nc


        Xc = Centroids(count);
        Yc = Centroids(count+1);

        radius = r(i);

        bwtemp = sqrt((X-Xc).^2 + (Y-Yc).^2) <= radius;
        %bw2 = sqrt((x-center2).^2 + (y-center2).^2) <= radius;
        %bw = bw1 | bw2;

        BWgrid= BWgrid | bwtemp;
        count = count+2;
    end

else
    if size(Centroids,1)==2
        Centroids = Centroids';
    end

    for i=1:nc


        Xc = Centroids(i,1);
        Yc = Centroids(i,2);
        %center2 = -center1;
        %dist = sqrt(2*(2*center1)^2);
        radius = r(i);
        % lims = [floor(center1-1.2*radius) ceil(center1+1.2*radius)];
        % [x,y] = meshgrid(lims(1):lims(2));
        bwtemp = sqrt((X-Xc).^2 + (Y-Yc).^2) <= radius;
        %bw2 = sqrt((x-center2).^2 + (y-center2).^2) <= radius;
        %bw = bw1 | bw2;

        BWgrid = BWgrid | bwtemp;
        count = count+2;
    end
end

CW = BWgrid;
%  figure
%  imshow(BWgrid)
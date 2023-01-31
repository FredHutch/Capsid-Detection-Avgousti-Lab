function [tdata,scDat, FigH]=detectcapsidsEM(imfile)
% This function identifies capsids within a user defined area and measures
% their distance from the nuclear membrane. It produces a figure showing
% the region of interest with the perimeter of the user defined ROI
% highlighted in yellow, capsids within 200nm of the nuclear envelope
% highlighted in red, and other capsids highlighted in blue.
%
% INPUT             -imfile: string name of the image file
%
% OUTPUT            -tdata: 1-by-5 table that summarized data with
%                           variables:
%                       --Nuc.Area, Area of the nuclear mask in um
%                       --Nuc.Perim., perimeter of the nuclear mask, in um
%                       --Total#Capsids, total number of counted capsids
%                       within the nuclear mask
%                       --#Capsids200nm, number of capsids within 200nm
%                       range of the nuclear membrane
%
%                   -scdata: nCapsids-by-4 table containing single capsids
%                   data with 4 variables:
%                      --Area, the area of the capsid
%                      --DistancefromEnveloppe, distance of the capsid from
%                      the nuclear membrane, in um
%                      --EmptynessFactor, built-in factor aimed at
%                      characterizing the empty nature of the capsid, one
%                      tending to completely empty
%                      --MeanIntensity, mean pixel intensity of the capsid
%                      (from complement image)
%
%                   -FigH: a handle of the generated figure
%
% Julien Dubrulle, 2023

%% Load image and get physical pixel size

IS = bfopen(imfile);
I = cast(IS{1,1}{1,1},'uint8');
if ndims(I) == 3
    I = rgb2gray(I);
end
omeMeta = IS{1,4};
if ~isempty(omeMeta.getPixelsPhysicalSizeX(0))
    Xvoxel = omeMeta.getPixelsPhysicalSizeX(0).value();
    VoxX = Xvoxel.doubleValue(); %express length in microns
else
    VoxX = (100 / 60) / 1000; %assuming detected diameter is 50-60 pixels and capsid diameters is 100nm (expressed  in microns)
end


%% Create mask of nucleus of interest

figure, imshow(I)
roi = drawfreehand(gca);
NUCMask = createMask(roi);
close;
S = regionprops(NUCMask,'Area','Perimeter','BoundingBox');
if numel(S) > 1
    disp('NUCLEAR MASK CONTAINS MORE THAN ONE OBJECT!')
    return
else
    I = imcrop(I,S.BoundingBox);
    NUCMask = imcrop(NUCMask,S.BoundingBox);
    [yy,xx] = size(I);
end

%% Find capsids

Ic = imcomplement(I);
Ica = imadjust(Ic);
Icaf = medfilt2(Ica,[10 10]);
%Icafr=imresize(Icaf,0.5);
Icafr = imresize(Icaf,1);
[xc, yc] = size(Icafr);
[centers,radii] = imfindcircles(Icafr,[25 30],'Sensitivity',0.89, 'Method','TwoStage');
CW = drawCircle2(yc,xc,centers,radii);
CW = imresize(CW,[yy, xx]);
MW = CW & NUCMask;

BWD = bwdist(~NUCMask);


%% Categorize capsids (A-empty; B-Intermediate; C-Full)


StatsJ = regionprops(MW,Icaf,'Area','Centroid','MeanIntensity','PixelValues','PixelIdxList');

aj = [StatsJ.Area]';
tfe = aj < 1650 | aj > 3000;
StatsJ(tfe) = [];


PID = {StatsJ.PixelIdxList}';
PV = {StatsJ.PixelValues}';
PV = cellfun(@double, PV, 'UniformOutput', false);
iiqqrr = cellfun(@iqr,PV);

MI = [StatsJ.MeanIntensity]';
FullIndex = iiqqrr./ MI;

CapsD = cellfun(@(x) mean(BWD(x)),PID);

CapsArea = [StatsJ.Area]' * VoxX^2;
CapsD = CapsD * VoxX;
tfDist = CapsD < 0.2; %less than 200nm
scDat = table(CapsArea,CapsD,FullIndex,MI,'VariableNames',{'Area','DistancefromEnveloppe','EmptynessFactor','MeanIntensity'});




%% plotting



CWAll = false(size(CW));
PixAll = vertcat(StatsJ(:).PixelIdxList);
CWAll(PixAll) = true;
CWAll = imresize(CWAll,size(I));
BAll = bwboundaries(CWAll);

CWClose = false(size(CW));
PixClose = vertcat(StatsJ(tfDist).PixelIdxList);
CWClose(PixClose) = true;
CWClose = imresize(CWClose,size(I));
BClose = bwboundaries(CWClose);




NBound = bwboundaries(NUCMask);
AreaNuc = numel(find(NUCMask)) * VoxX^2;
PerNuc = S.Perimeter * VoxX;

data = [AreaNuc PerNuc numel(StatsJ) sum(tfDist) sum(tfDist) / PerNuc];
tdata = array2table(data,'VariableNames',{'Nuc.Area','Nuc.Perim.','Total # Capsids','#Capsids200nm','ratio'});

FigH = figure('Name',imfile); imshow(I); hold on;

for i = 1:length(BAll)
    bb = BAll{i};
    plot(gca,bb(:,2),bb(:,1),'-b','LineWidth',2);
end

for i = 1:length(BClose)
    bb = BClose{i};
    plot(gca,bb(:,2),bb(:,1),'-r','LineWidth',2);
end


for i = 1:length(NBound)
    nn = NBound{i};
    plot(gca,nn(:,2),nn(:,1),'-y','LineWidth',1);
end






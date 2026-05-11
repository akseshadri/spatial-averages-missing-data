function [latu, lonu, datmat] = vect2matgeophys(lat,lon,dat)
% VECT2MATGEOPHYS  Reshape vector-indexed geophysical data to a lat/lon matrix.
%
%   [latu, lonu, datmat] = vect2matgeophys(lat, lon, dat) takes data stored
%   as a vector indexed by station/grid-point number and reshapes it into a
%   2-D matrix indexed by unique latitude and longitude values. This is
%   needed for contour plotting with contourm.
%
%   Inputs:
%     lat — Vector of latitudes for each observation point (n x 1).
%     lon — Vector of longitudes for each observation point (n x 1).
%     dat — Data vector to be reshaped (n x 1), e.g. optimal weights.
%
%   Outputs:
%     latu   — Sorted unique latitudes (nlat x 1).
%     lonu   — Sorted unique longitudes (nlon x 1).
%     datmat — Data matrix (nlon x nlat), with NaN where no observation
%              exists. Oriented for use with contourm(latu, lonu, datmat').
%
%   Used in PFig2_3.m, PFig4.m, and PFig6abc.m for spatial plotting.

latu = unique(lat);
lonu = unique(lon);

nilat = NaN(size(lat));
nilon = NaN(size(lon));

for i = 1:numel(latu)
    latui = latu(i);
    indexi = ismember(lat,latui);
    nilat(indexi) = i;
end

for i = 1:numel(lonu)
    lonui = lonu(i);
    indexi = ismember(lon,lonui);
    nilon(indexi) = i;
end

datmat = NaN(numel(lonu), numel(latu));
for i = 1:numel(lat)
    datmat(nilon(i),nilat(i)) = dat(i);
end
    
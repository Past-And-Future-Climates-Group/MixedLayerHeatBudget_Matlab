% Function to perform averaging of a three-dimensional field over a mixed
% layer depth taken from an input field
% February 7, 2013
% Sam Stevenson
% Updated September 19, 2013 to use a time-varying mixed layer depth
% (appropriate for ROMS model output)

% Inputs:
% mld   -   matrix containing mixed-layer depth field, dimensions
%           nt x nlat x nlon, units of m
% field -   matrix containing desired field to integrate, dimensions nt x
%           nz x nlat x nlon
% time  -   array containing measurement times for field, dimensions nt
%           units should be in days since 0000-01-01 for compatibility with
%           the Matlab "datenum" function
% z     -   array of depths, units of m

% Outputs:
% fldint -  matrix containing the values of the input field averaged over
%           the mixed layer, dimensions nt x nlat x nlon

function [fldint]=mldavg_varytime(mld,field,time,z)

    for tt=1:size(field,1)
        for la=1:size(field,3)
            for lo=1:size(field,4)
                myz=find(z(:,la,lo) > mld(tt,la,lo));
                field(tt,myz,la,lo)=0/0;
            end
        end
    end
    
    % Average
    fldint=squeeze(nanmean(field,2));
end
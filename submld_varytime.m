% Function to find the value of a three-dimensional field just below a mixed
% layer depth taken from an input field
% March 2015
% Sam Stevenson

% Inputs:
% mld   -   matrix containing mixed-layer depth field, dimensions
%           nt x nlat x nlon, units of m
% field -   matrix containing desired field to integrate, dimensions nt x
%           nz x nlat x nlon
% time  -   array containing measurement times for field, dimensions nt
%           units should be in days since 0000-01-01 for compatibility with
%           the Matlab "datenum" function
% z     -   array of depths, units of m
% type  -   type of searching to do on depth array. For CESM/POP: should
%           use 'first', for ROMS 'last' since depths are ordered backwards

% Outputs:
% fldint -  matrix containing the values of the input field averaged over
%           the mixed layer, dimensions nt x nlat x nlon

function [fldint]=submld_varytime(mld,field,time,z,type)
    fldint=zeros(size(field,1),size(field,3),size(field,4));

    for tt=1:size(field,1)
        for la=1:size(field,3)
            for lo=1:size(field,4)
                myz=find(z(:,la,lo) > mld(tt,la,lo),1,type);
                if ~isempty(myz)
                    fldint(tt,la,lo)=field(tt,myz,la,lo);
                end
            end
        end
    end
end
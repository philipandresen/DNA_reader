function [ POI chunksize ] = lookforgoodchunks( xmatch,threshold )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
currentsize=0;
chunksize=0;
POI=0;
ind=1;
chunk_number=1;
%FIND THE POINTS OF INTEREST
while ind+currentsize<length(xmatch);
    if xmatch(ind)==1;
        while xmatch(ind+currentsize)==1;
            if ind+currentsize>=length(xmatch);
                break;
            else
                currentsize=currentsize+1;
            end;
        end;
        %currentsize=currentsize-1;
    end;
    if currentsize>threshold;
        POI(chunk_number)=ind;
        chunksize(chunk_number)=currentsize;
        chunk_number=chunk_number+1;
        %disp(currentsize)
    end;
    %if currentsize<0; currentsize=0; end;
    ind=ind+currentsize+1;
    currentsize=0;
end;
% if POI~=1;
%     POIS=[POIS POI];
% end;
%END POINTS OF INTEREST
end


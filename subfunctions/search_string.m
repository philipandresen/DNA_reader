function [ location instances inframe] = search_string(instring, searchstr,...
    number)
%search_string by Philip Andresen (Version 23:AUGUST:2011)
%INTENDED CALLER: ANY
%PURPOSE: This function searches a string of characters for the nth instance of a
%   number of sequential characters. This is basically strfind but it is a
%   strfind that will return the nth instance of the string location rather
%   than all instances or just the first.
%INPUTS:
%   instring: The string to be searched within
%   searchstr: The string for which to search.
%   number: The location of searchstr in instr to return. (1st, 2nd, 3rd...)
%OUTPUTS:
%   Location: The position in which the string was found.
%   instances: The total number of instances found.
%   inframe: Whether or not the string was found in-frame (DNA nucleotides)
%CHANGELOG:
%   Changes have not been logged as of (23:AUGUST:2011)
%External function dependencies:
%   None
%SPECIAL NOTES: 
%   This is intended primarily for searching for nucleotide sequences
%   within other nucleotide sequences.
found=strfind(upper(instring),upper(searchstr));
instances=length(found);
location=1;
inframe=1;
if instances>0
    if number<=instances;
        location=floor((2+found(1,number))/3);
    end;
    if mod(2+found(1,number),3)>0;
        inframe=0;
    end;
end;
end


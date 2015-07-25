function [CODE] = formatcode(CODE)
%formatcode by Philip Andresen (Version 23:AUGUST:2011)
%INTENDED CALLER: ANY
%PURPOSE: This program simply takes input strings and removes all
%   characters that are not 'A','G','C','T', or '_'. It also capitalizes
%   the string.
%INPUTS:
%   CODE: The code to be formatted
%OUTPUTS:
%   CODE: The fromatted code.
%CHANGELOG:
%   Changes have not been logged as of (23:AUGUST:2011)
%External function dependencies:
%   None
%SPECIAL NOTES: 
%   Meant to completely replace addcode.m, this code was written to reduce
%   the size of DNA_reader.m by compacting code that is called many times
%   in many places.

CODE=upper(CODE);
for i=1:20; CODE=strrep(CODE,' ',''); end;
IND=sort([strfind(CODE,'A') strfind(CODE,'G')...
    strfind(CODE,'C') strfind(CODE,'T') strfind(CODE,'_')]);
CODE=CODE(IND);



end
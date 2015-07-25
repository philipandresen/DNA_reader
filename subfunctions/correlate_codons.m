function [ cellarray ] = correlate_codons( handles )
%correlate_codons by Philip Andresen (Version 23:AUGUST:2011)
%INTENDED CALLER: DNA_reader.m
%PURPOSE: This program takes the handles structure from DNA_reader and
%   returns a cell matrix containing 'good' 'No correlation' and 'silent
%   mutation' that will ultimately wind up in the third column of the
%   program's text boxes.
%INPUTS:
%   handles: The handles structure from DNA_reader
%OUTPUTS:
%   cellarray: A cell matrix with correlation quality words.
%CHANGELOG:
%   Changes have not been logged as of (23:AUGUST:2011)
%External function dependencies:
%   codonify.m
%SPECIAL NOTES: 
%   none.
index=1;
A=codonify(handles.data1,1,'none');
B=codonify(handles.data2,1,'none');
C=codonify(handles.data1,1,'long');
D=codonify(handles.data2,1,'long');
%cellarray=zeros(size(A));
cellarray(index)={'No Correlation'};
while index<=length(A) && index<=length(B);
    cellarray(index)={'No Correlation'};
    if strcmp(A{index},B{index}); 
        cellarray(index)={'Good'}; 
    end;
    if strcmp(C{index},D{index}) && ~strcmp(A{index},B{index});
        cellarray(index)={'Silent Mut.'};
    end;
    index=index+1;
end;

end


function [ reference_sequence corrected_sequence ] = align_DNA(reference_sequence,target_sequence,requirement)
%align_DNA by Philip Andresen (Version 23:AUGUST:2011)
%INTENDED CALLER: DNA_reader.m
%PURPOSE: Ths function takes two DNA sequences and attempts to correct for
%   shifting caused by insertions and deletions. The method by which this
%   is done depends on the fitting of significantly large chunks of DNA as
%   the two strands are shifted past one another. Significance is defined
%   by the "requirement" input variable, which is a percentage (0.00-1.00)
%   of the total sequence length required to define a chunk. Chunks below
%   this requirement are not counted in the alignment process.
%INPUTS:
%   reference_sequence: The DNA to which the user wishes to fit.
%   target_sequence: The DNA that will be shifted to fit the ref. seq.
%   requirement: A minimum percentage (0-1) of the DNA length that is
%       required to define fitting chunks
%OUTPUTS:
%   reference_sequence: Will be the same as the input but may have a 1 or
%       0 appended to the end during the alignment process.
%   corrected_sequence: The target_sequence that has been shifted,
%       inserted, and deleted to fit sequentially with the reference 
%       sequence.
%CHANGELOG:
%   Changes have not been logged as of (23:AUGUST:2011)
%External function dependencies:
%   lookforgoodchunks.m
%clc;
original_reference=reference_sequence;
hwait=waitbar(0,'Aligning code sequences');
%NOTE! X and A are reference DNA
%Y and B are to be modified!
X=reference_sequence;
Y=target_sequence;
while length(X)>length(Y); Y=[Y '0']; end; %#ok<*AGROW>
while length(X)<length(Y); X=[X '1']; end;
resultant_values=zeros(1,length(X));
resultant_shifts=resultant_values;
saved_chunks='';
A=X;
B=Y;
Z(1:length(Y))=1;%Z is a logic mask for Y to prevent codon reuse
furthestYmatch=0;
FRAMESHIFTLOG{1,1}='No frame shifts found';
numshifts=0;
matchinggraph{1}=X==X;
titleofgraph{1}='Reference DNA, 100% self correlation';
numgraphs=1;
%Shifting B left
for index=0:length(X)-1;
    if mod(index,floor(2*length(X)/50))==1;
        waitbar(index/(2*length(X)),hwait,'Aligning code sequences');
    end;
    %disp([num2str(sum(A==B)) ' at ' num2str(index)])
    resultant_values(1,length(X)-index)=sum((A==B))+((index))/4-length(X)/4;
    resultant_shifts(1,length(X)-index)=length(X)-index;
    %disp([num2str(resultant_values(1,length(X)-index)) ' seen at ' num2str(index)])
    if resultant_values(1,length(X)-index)>length(X)/16; 
        xmatch=(A==B);
        %disp(['Xmatch is (top)' num2str(sum(xmatch))]);
        [POIS chunksizes]=lookforgoodchunks(xmatch,requirement*length(X)/100);
        %figure('Name',['Offset of ' num2str(index)]);
        if sum(chunksizes)>0; %dont make a graph if there aren't POIS
            %disp(['POIs found: ' num2str(POIS) '--' num2str(POIS+chunksizes) ...
            %' at shift ' num2str(index)])
            numgraphs=numgraphs+1;
            titleofgraph{numgraphs}=[num2str(sum(chunksizes)) ' correlations with offset of ' num2str(index)];
            matchinggraph{numgraphs}=xmatch(:);
        end;
        %bar(xmatch);
        for COUNTER=1:length(POIS) %THis will allow for multiple chunks
            POI=POIS(COUNTER);
            chunksize=chunksizes(COUNTER);
            if chunksize==0; continue; end; %if no points of interest were found
            try
                totalshift=index;
                if POI+chunksize>length(Y); chunksize=length(Y)-POI; end;%changed from >= to just =
                %masked=Y(Y & Z); %Eliminate already-correlated nucleotides
                masked=Y; %this used to be masked=Y(Y&Z) but that removed blank space.
                masked(~(Y&Z))='_';
                if masked(POI:POI+chunksize)==Y(POI:POI+chunksize); %If no pre-used nucleotides in this chunk
                    saved_chunks(POI-totalshift:POI+chunksize-totalshift)=masked(POI:POI+chunksize);
                    Z(POI:POI+chunksize-1)=0;
                    if POI+chunksize>furthestYmatch; furthestYmatch=POI+chunksize; end;
                    numshifts=numshifts+1;
                    if POI~=1; FRAMESHIFTLOG{numshifts,1}=['Found insertion at ' num2str(POI)]; end;
                else
                    disp('PRE USED DETECTED 1')
                end;
            catch exception
                disp(exception.message)
                disp(exception.stack)
            end;
        end;
    end;
    A=['1' A];
    B=[B '0'];
end;
%Shifting B right
A=X;
B=Y;

for index=length(X):2*length(X);
    if mod(index,floor(2*length(X)/50))==1;
        waitbar(index/(2*length(X)),hwait,'Aligning code sequences');
    end;
    resultant_values(1,index)=sum((A==B))+(index-length(X))/4-length(X)/4;
    resultant_shifts(1,index)=index;
    if resultant_values(1,index)>length(X)/16; 
        xmatch=(A==B);
        [POIS chunksizes]=lookforgoodchunks(xmatch,requirement*length(X)/100);
        %figure('Name',['Offset of ' num2str(length(X)-index)]);
        if sum(chunksizes)>0; %dont make a graph if there aren't POIS
            %disp(['POIS found: ' num2str(POIS) '--' num2str(POIS+chunksizes) ...
            %' at shift ' num2str(length(X)-index)])
            numgraphs=numgraphs+1;
            titleofgraph{numgraphs}=[num2str(sum(chunksizes)) ' correlations with offset of ' num2str(length(X)-index)];
            matchinggraph{numgraphs}=xmatch(:);
        end;
        %bar(xmatch);
        for COUNTER=1:length(POIS) %THis will allow for multiple chunks
            POI=POIS(COUNTER);
            chunksize=chunksizes(COUNTER);
            if chunksize==0; continue; end; %if no points of interest were found
            try
                totalshift=length(X)-index;
                if POI+chunksize+totalshift>length(Y); chunksize=length(Y)-POI-totalshift; end;%changed from >= to just =
                %if POI+totalshift<1; POI=1-totalshift; end;
                masked=Y; %this used to be masked=Y(Y&Z) but that removed blank space.
                masked(~(Y&Z))='_';
                if masked(POI+totalshift:POI+chunksize+totalshift)==Y(POI+totalshift:POI+chunksize+totalshift);
                    saved_chunks(POI:POI+chunksize)=masked(POI+totalshift:POI+chunksize+totalshift);
                    Z(POI+totalshift:POI+chunksize+totalshift-1)=0;
                    if POI+totalshift+chunksize>furthestYmatch; furthestYmatch=POI+totalshift+chunksize; end;
                    numshifts=numshifts+1;
                    FRAMESHIFTLOG{numshifts,1}=['Found deletion at ' num2str(POI+totalshift)];
                else
                    disp('PRE USED DETECTED 2')
                    %disp(masked)
                    %disp(Y)
                end;
                %disp(length(Z))
                %disp(length(Y))
            catch exception
                disp(exception.message)
                disp(exception.stack)
            end;
        end;
    end;
    A=[A '1'];
    B=['0' B];
end
while length(X)>length(saved_chunks); saved_chunks=[saved_chunks '0']; end;
while length(X)<length(saved_chunks); X=[X '1']; end;
sucalignments=sum(X==saved_chunks);
figure;
bar(resultant_values)
disp(['Correlation threshold for this sequence was ' num2str(requirement*length(X)/100) ' nucleotides'])
%disp(['New correlation: ' num2str(sum(saved_chunks==X))])
%disp(X)
%disp(saved_chunks)
disp(['Sequences were synchronized with ' num2str(sum(X==saved_chunks))...
    ' nucleotides out of ' num2str(length(X))])

titleofgraph{1}=['Correlation after correction- Threshold:' num2str(requirement*length(X)/100)...
    ' with ' num2str(sum(X==saved_chunks)) '/' num2str(length(X)) ' correlated'];
matchinggraph{1}=X==saved_chunks;
saved_chunks(saved_chunks==0)='_';
saved_chunks(saved_chunks=='0')='';
X(X=='1')='';
% try
%      saved_chunks(end:length(Z))=Y(length(saved_chunks)+1:length(Z));
%  catch %#ok<CTCH> %Just try this, not super important
%  end;
% disp('saved chunks:')
% disp(saved_chunks);
% disp(deblank(num2str(Z)))
% disp('Y:')
% disp(Y);
% disp('X')
% disp(X)
% disp(furthestYmatch)
saved_chunks(end+1:end+(length(Y)-furthestYmatch))=Y(furthestYmatch+1:end); %This was added to stick uncorrelated regions onto the end of the return string.
% disp(saved_chunks)
%saved_chunks(Z==1)=Y(Z==1);
n=length(matchinggraph);
if n>15; n=15; end; %Make sure it doesn't try to draw too much.
if sucalignments==0; 
    warndlg('Correlation was not possible. Returning uncorrected code.','Failed');
    corrected_sequence=target_sequence;
    reference_sequence=original_reference;
else
    reference_sequence=X;
    corrected_sequence=saved_chunks;
    figure('Name','Correlation graphs (for first 15 or fewer shifts)');
    for randomindex=1:n
        subplot(n,1,randomindex)
        bar(matchinggraph{randomindex})
        title(titleofgraph{randomindex})
        set(gca,'XTickLabel','');
        set(gca,'YTickLabel','');
    end
end;
close(hwait);
disp(FRAMESHIFTLOG)

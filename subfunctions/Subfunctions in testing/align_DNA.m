function [ reference_sequence corrected_sequence ] = align_DNA(reference_sequence,target_sequence)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clc;
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
matchinggraph{1}=X==X;
titleofgraph{1}='Reference DNA, 100% self correlation';
numgraphs=1;
%Shifting B left
for index=0:length(X)-1;
    resultant_values(1,length(X)-index)=sum(A==B)+((index))/4-length(X)/4;
    resultant_shifts(1,length(X)-index)=length(X)-index;
    if resultant_values(1,length(X)-index)>length(X)/16; 
        xmatch=A==B;
        [POIS chunksizes]=lookforgoodchunks(xmatch,length(X)/32);
        disp(['POIs found: ' num2str(POIS) '--' num2str(POIS+chunksizes) ...
            ' at shift ' num2str(index)])
        %figure('Name',['Offset of ' num2str(index)]);
        numgraphs=numgraphs+1;
        titleofgraph{numgraphs}=['Correlation with offset of ' num2str(index)];
        matchinggraph{numgraphs}=xmatch(:);
        %bar(xmatch);
        for COUNTER=1:length(POIS) %THis will allow for multiple chunks
            POI=POIS(COUNTER);
            chunksize=chunksizes(COUNTER);
            try
                totalshift=index;
                if POI+chunksize+totalshift>=length(Y); chunksize=length(Y)-POI-totalshift; end;
                saved_chunks(POI:POI+chunksize)=Y(POI+totalshift:POI+chunksize+totalshift);
            catch exception
                disp(exception)
                disp(POI+totalshift)
                disp(POI+chunksize)
                disp(length(Y))
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
    resultant_values(1,index)=sum(A==B)+(index-length(X))/4-length(X)/4;
    resultant_shifts(1,index)=index;
    if resultant_values(1,index)>length(X)/16; 
        xmatch=A==B;
        [POIS chunksizes]=lookforgoodchunks(xmatch,length(X)/32);
        disp(['POIS found: ' num2str(POI) '--' num2str(POIS+chunksizes) ...
            ' at shift ' num2str(length(X)-index)])
        %figure('Name',['Offset of ' num2str(length(X)-index)]);
        numgraphs=numgraphs+1;
        titleofgraph{numgraphs}=['Correlation with offset of ' num2str(length(X)-index)];
        matchinggraph{numgraphs}=xmatch(:);
        %bar(xmatch);
        for COUNTER=1:length(POIS) %THis will allow for multiple chunks
            POI=POIS(COUNTER);
            chunksize=chunksizes(COUNTER);
            try
                totalshift=length(X)-index;
                %if POI+chunksize>=length(Y); chunksize=length(Y)-POI; end;
                %if POI+totalshift<1; POI=1-totalshift; end;
                saved_chunks(POI:POI+chunksize)=Y(POI+totalshift:POI+chunksize+totalshift);
            catch exception
                disp(exception)
                %disp(POI-totalshift)
                %disp(POI+chunksize-totalshift)
            end;
        end;
    end;
    A=[A '1'];
    B=['0' B];
end
while length(X)>length(saved_chunks); saved_chunks=[saved_chunks '0']; end;
while length(X)<length(saved_chunks); X=[X '1']; end;
figure;
bar(resultant_values)
disp(['Correlation threshold for this sequence was ' num2str(ceil(length(X))) ' nucleotides'])
%disp(['New correlation: ' num2str(sum(saved_chunks==X))])
disp(X)
disp(saved_chunks)
disp(['Sequences were synchronized with ' num2str(sum(X==saved_chunks))...
    ' nucleotides out of ' num2str(length(X))])

titleofgraph{1}=['Correlation after correction- Threshold:' num2str(length(X)/16)...
    ' with ' num2str(sum(X==saved_chunks)) '/' num2str(length(X)) ' correlated'];
matchinggraph{1}=X==saved_chunks;

%saved_chunks
reference_sequence=X;
corrected_sequence=saved_chunks;
n=length(matchinggraph);
figure;
for randomindex=1:n
    subplot(n,1,randomindex) 
    bar(matchinggraph{randomindex})
    title(titleofgraph{randomindex})
end


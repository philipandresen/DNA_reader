function [ DNAmatrix ] = codonify( raw_sequence,cellwidth,translate)
%Codonify by Philip Andresen (Version 23:AUGUST:2011)
%INTENDED CALLER: ANY
%PURPOSE: This program takes a raw sequence of characters and separates the
%   data into codons in a 'cellwidth' by n array where n is however many  
%   rows it takes to represent all the data. It also will translate codons
%   into amino acid long/short names. The input string can be anything and
%   the codonify program will eliminate all non-(A,G,C,T,_) characters.
%INPUTS:
%   raw_sequence: The nucleotide sequence to be formatted
%   cellwidth: The number of codon triplets wide the output should be
%   translate: ('none'/'short'/'long') decides to translate the codons into
%       amino acid short/long names or does not translate at all.
%OUTPUTS:
%   DNAmatrix: This is the formatted/segmented code in cell array form.
%CHANGELOG:
%   Changes have not been logged as of (23:AUGUST:2011)
%External function dependencies:
%   None
%SPECIAL NOTES: The first few lines may be reduntant from the formatcode.m
%   file. This should be carefully investigated and replaced.
%   Translation Source:
%   http://www.colorado.edu/chemistry/bioinfo/AminoAciddata.htm
for i=1:20; raw_sequence=strrep(raw_sequence,' ',''); end;
IND=sort([strfind(raw_sequence,'A') strfind(raw_sequence,'G')...
    strfind(raw_sequence,'C') strfind(raw_sequence,'T') strfind(raw_sequence,'_')]);
raw_sequence=upper(raw_sequence(IND));
codon=''; %#ok<NASGU>
index1=1;
index2=1;
currentchar=1;
datasize=length(raw_sequence)+1;
DNAmatrix=num2cell(zeros(floor((datasize/3)/cellwidth),cellwidth));
while currentchar+3<=datasize
    codon=raw_sequence(currentchar:currentchar+2);
    if strcmp(translate,'short') || strcmp(translate,'long');
        switch codon
            case 'TTT'; codon='Phe';
            case 'TTC'; codon='Phe';
            case 'TTA'; codon='Leu';
            case 'TTG'; codon='Leu';
            case 'CTT'; codon='Leu';
            case 'CTC'; codon='Leu';
            case 'CTA'; codon='Leu';
            case 'CTG'; codon='Leu';
            case 'ATT'; codon='Ile';
            case 'ATC'; codon='Ile';
            case 'ATA'; codon='Ile';
            case 'ATG'; codon='Met';
            case 'GTT'; codon='Val';
            case 'GTC'; codon='Val';
            case 'GTA'; codon='Val';
            case 'GTG'; codon='Val';
                
            case 'TCT'; codon='Ser';
            case 'TCC'; codon='Ser';
            case 'TCA'; codon='Ser';
            case 'TCG'; codon='Ser';
            case 'CCT'; codon='Pro';
            case 'CCC'; codon='Pro';
            case 'CCA'; codon='Pro';
            case 'CCG'; codon='Pro';
            case 'ACT'; codon='Thr';
            case 'ACC'; codon='Thr';
            case 'ACA'; codon='Thr';
            case 'ACG'; codon='Thr';
            case 'GCT'; codon='Ala';
            case 'GCC'; codon='Ala';
            case 'GCA'; codon='Ala';
            case 'GCG'; codon='Ala';
                
            case 'TAT'; codon='Tyr';
            case 'TAC'; codon='Tyr';
            case 'TAA'; codon='STOP';
            case 'TAG'; codon='STOP';
            case 'CAT'; codon='His';
            case 'CAC'; codon='His';
            case 'CAA'; codon='Gln';
            case 'CAG'; codon='Gln';
            case 'AAT'; codon='Asn';
            case 'AAC'; codon='Asn';
            case 'AAA'; codon='Lys';
            case 'AAG'; codon='Lys';
            case 'GAT'; codon='Asp';
            case 'GAC'; codon='Asp';
            case 'GAA'; codon='Glu';
            case 'GAG'; codon='Glu';
                
            case 'TGT'; codon='Cys';
            case 'TGC'; codon='Cys';
            case 'TGA'; codon='STOP';
            case 'TGG'; codon='Trp';
            case 'CGT'; codon='Arg';
            case 'CGC'; codon='Arg';
            case 'CGA'; codon='Arg';
            case 'CGG'; codon='Arg';
            case 'AGT'; codon='Ser';
            case 'AGC'; codon='Ser';
            case 'AGA'; codon='Arg';
            case 'AGG'; codon='Arg';
            case 'GGT'; codon='Gly';
            case 'GGC'; codon='Gly';
            case 'GGA'; codon='Gly';
            case 'GGG'; codon='Gly';
        
        end;
    end;
    if strcmp(translate,'long');
        codon=strrep(codon,'Phe','Phenylalanine');
        codon=strrep(codon,'Leu','Leucine');
        codon=strrep(codon,'Ile','Isoleucine');
        codon=strrep(codon,'Met','Methionine');
        codon=strrep(codon,'Val','Valine');
        codon=strrep(codon,'Ser','Serine');
        codon=strrep(codon,'Pro','Proline');
        codon=strrep(codon,'Thr','Threonine');
        codon=strrep(codon,'Ala','Alanine');
        codon=strrep(codon,'Tyr','Tyrosine');
        codon=strrep(codon,'His','Histidine');
        codon=strrep(codon,'Glu','Glutamic acid');
        codon=strrep(codon,'Gln','Glutamine');
        codon=strrep(codon,'Lys','Lysine');
        codon=strrep(codon,'Asp','Aspartic acid');
        codon=strrep(codon,'Asn','Asparagine');
        codon=strrep(codon,'Cys','Cysteine');
        codon=strrep(codon,'Trp','Tryptophan');
        codon=strrep(codon,'Arg','Arginine');
        codon=strrep(codon,'Gly','Glycine');
    end;
    DNAmatrix(index2,index1)={codon}; 
    currentchar=currentchar+3;
    index1=index1+1;
    if index1>cellwidth; 
        index1=1;
        index2=index2+1;
    end;
end


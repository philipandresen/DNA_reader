function [handles ] = addcode(handles,CODE,whichcolumn )
%Addcode.m by Philip Andresen. Version (23:AUGUST:2011)
%INTENDED CALLER: DNA_reader
%PURPOSE: To simplify recurring assignments and formatting made in the
%   parent program. Essentially to make better organized, cleaner code.
%NOTICE:
%   This code has been made obsolete by formatcode.m, which performs the
%   same task but more officintly and without errors.
switch whichcolumn
    case 1
        set(handles.sourcetext1,'string',CODE)
        CODE=upper(get(handles.sourcetext1,'string'));
        for i=1:20; CODE=strrep(CODE,' ',''); end;
        IND=sort([strfind(CODE,'A') strfind(CODE,'G')...
            strfind(CODE,'C') strfind(CODE,'T') strfind(CODE,'_')]);
        CODE=CODE(IND);
        set(handles.codetext1,'string',codonify(CODE,1,handles.mode))
        set(handles.codetext1,'value',1)
    case 2
        set(handles.sourcetext2,'string',CODE)
        CODE=upper(get(handles.sourcetext2,'string'));
        for i=1:20; CODE=strrep(CODE,' ',''); end;
        IND=sort([strfind(CODE,'A') strfind(CODE,'G')...
            strfind(CODE,'C') strfind(CODE,'T') strfind(CODE,'_')]);
        CODE=CODE(IND);
        set(handles.codetext2,'string',codonify(CODE,1,handles.mode))
        set(handles.codetext2,'value',1)
end;

handles.result=correlate_codons(handles);
set(handles.resultslistbox,'string',handles.result)


end


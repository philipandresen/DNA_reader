function varargout = DNA_reader(varargin)
% DNA_READER by Philip Andresen (Version 24:August:2011)
%INTENDED CALLER: USER
%PURPOSE: DNA_reader takes and correlates two strings of DNA bases, translating
%     them into corresponding amino acids and identifying mutations. The
%     reader can also search for base sequences within larger sequences,
%     and can attempt to correlate frame shifted data with a quantified
%     degree of accuracy given afterwards (Histogram and %alignment)
%INPUTS:
%   N/A
%OUTPUTS:
%   N/A
%CHANGELOG: Last modified 8-24-2011 by Philip Andresen
%   In: subfunctions\align_DNA.m
%      Change: Completed the Z masking on the input data to prevent codon
%      re-use in the alignment algorithm. Also added a feature where
%      unaligned code at the end of the input sequence is no longer deleted
%      so that the data is easier to analyse.
%External function dependencies:
%   subfunctions\align_DNA.m
%   subfunctions\code_database.m
%   subfunctions\codonify.m
%   subfunctions\correlate_codons.m
%   subfunctions\find_badpoint.m
%   subfunctions\formatcode.m
%   subfunctions\lookforgoodchunks.m
%   subfunctions\search_string.m
%SPECIAL NOTES: 
%   This is a GUI and is intedned to be run explicitly by the user. More
%   information on fitting routines and functionality may be found by
%   referencing the help text of the dependent functions (Above).

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DNA_reader_OpeningFcn, ...
    'gui_OutputFcn',  @DNA_reader_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before DNA_reader is made visible.
function DNA_reader_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.output = hObject;
handles.mode='long';
handles.data2='';
handles.data1='';
handles.result={''};
handles.badpointnumber=0;
handles.searchstr='';
handles.searchstringindex=1;
handles.instances=0;
addpath('subfunctions')
presets={'Presets:','Hemagglutinin (HA) (H3N2) EF614248.1','Pa-mCherry',...
    'Actin (Beta) (Human) NM_001101.3','Dendra2-N',...
    'Neuraminidase (NA) (X-31, H3N2)','Cofilin (Human)','Pa-mKate',...
    'Profilin (Human)','TM1 (Human) NM_001018007.1',...
    'TM5 NM1 (Human) NM_153649 ','EGFP - U57609.1'};
set(handles.presetsmenu1,'string',presets)
set(handles.presetsmenu2,'string',presets)
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = DNA_reader_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on selection change in codetext1.
function codetext1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
set(handles.resultslistbox,'value',get(handles.codetext1,'value'))
set(handles.codetext2,'value',get(handles.codetext1,'value'))
%pos=get(hObject,'Value')*3;
%handles.data1(pos-2:pos)

function codetext1_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in codetext2.
function codetext2_Callback(hObject, eventdata, handles)
set(handles.resultslistbox,'value',get(handles.codetext2,'value'))
set(handles.codetext1,'value',get(handles.codetext2,'value'))

function codetext2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sourcetext1_Callback(hObject, eventdata, handles)
DATA=upper(get(hObject,'string'));
handles.data1=formatcode(DATA);
set(handles.sourcetext1,'string',handles.data1);
set(handles.codetext1,'string',codonify(handles.data1,1,handles.mode));
if get(handles.codetext1,'value')>=floor(length(get(handles.sourcetext1,'string'))/3);
    set(handles.codetext1,'value',1);
end;
handles.result=correlate_codons(handles);
set(handles.resultslistbox,'string',handles.result);

len=length(get(handles.codetext1,'string'));
if len<length(get(handles.codetext2,'string')); len=length(get(handles.codetext2,'string')); end;
set(handles.matchqualitytext,'string',...
    [num2str(sum(strcmp(upper(handles.result),'GOOD'))/len*100)...
    '% match between sequences']) %#ok<*STCI>
set(handles.presetsmenu1,'value',1)
guidata(hObject,handles);

function sourcetext1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sourcetext2_Callback(hObject, eventdata, handles)
DATA=upper(get(hObject,'string'));
handles.data2=formatcode(DATA);
set(handles.sourcetext2,'string',handles.data2);
set(handles.codetext2,'string',codonify(handles.data2,1,handles.mode));
if get(handles.codetext2,'value')>=floor(length(get(handles.sourcetext2,'string'))/3);
    set(handles.codetext2,'value',1);
end;
handles.result=correlate_codons(handles);
set(handles.resultslistbox,'string',handles.result);

len=length(get(handles.codetext1,'string'));
if len<length(get(handles.codetext2,'string')); len=length(get(handles.codetext2,'string')); end;
set(handles.matchqualitytext,'string',...
    [num2str(sum(strcmp(upper(handles.result),'GOOD'))/len*100)...
    '% match between sequences']) %#ok<*STCI>
set(handles.presetsmenu2,'value',1)
guidata(hObject,handles);

function sourcetext2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag')
    case 'codonbutton';
        handles.mode='no';
    case 'shortbutton';
        handles.mode='short';
    case 'longbutton';
        handles.mode='long';
end;
set(handles.codetext1,'string',codonify(handles.data1,1,handles.mode));
set(handles.codetext2,'string',codonify(handles.data2,1,handles.mode));
guidata(hObject,handles);

% --- Executes on selection change in resultslistbox.
function resultslistbox_Callback(hObject, eventdata, handles)
set(handles.codetext1,'value',get(hObject,'value'))
set(handles.codetext2,'value',get(hObject,'value'))

function resultslistbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
[location]=...
    find_badpoint(handles.result,get(handles.resultslistbox,'value'),'down');
set(handles.resultslistbox,'value',location)
set(handles.codetext1,'value',location)
set(handles.codetext2,'value',location)
guidata(hObject,handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
[location]=...
    find_badpoint(handles.result,get(handles.resultslistbox,'value'),'up');
set(handles.resultslistbox,'value',location)
set(handles.codetext1,'value',location)
set(handles.codetext2,'value',location)
guidata(hObject,handles)

function searchtext_Callback(hObject, eventdata, handles)
handles.searchstr=get(hObject,'String');
handles.searchstringindex=1;
handles.instances=0;
guidata(hObject,handles)

function searchtext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in searchbutton.
function searchbutton_Callback(hObject, eventdata, handles)
data=handles.data2;
location=get(handles.codetext2,'value');
if get(handles.radiobutton4,'value')==1; 
    data=handles.data1; 
    location=get(handles.codetext1,'value');
end;
if handles.instances>0;
    handles.searchstringindex=...
        mod(handles.searchstringindex,handles.instances+1);
    if handles.searchstringindex==0; handles.searchstringindex=1; end;
end;
handles.instances=1;
tries=0;
inframe=0;
currentloc=location;
while inframe==0 && handles.instances>0 && tries<=length(data);
    [location handles.instances inframe]=...
        search_string(data,handles.searchstr,...
        handles.searchstringindex);
    handles.searchstringindex=...
        mod(handles.searchstringindex+1,handles.instances+1);
    tries=tries+1;
    if handles.searchstringindex==0; handles.searchstringindex=1; end;
    if get(handles.inframebutton,'value')==0; break; end;
end;
if tries>length(data); 
    handles.instances=0; 
    location=currentloc;
end;

if inframe==1; set(handles.textout,'string','Found, in frame.'); end;
if inframe==0; set(handles.textout,'string','Found, but not in frame.'); end;
if handles.instances==0; set(handles.textout,'string','Not Found.'); end;
set(handles.resultslistbox,'value',location)
set(handles.codetext1,'value',location)
set(handles.codetext2,'value',location)

guidata(hObject,handles)

% --- Executes on selection change in presetsmenu1.
function presetsmenu1_Callback(hObject, eventdata, handles)
if get(hObject,'Value')>1;
    DATA=code_database(get(hObject,'Value'));
    handles.data1=formatcode(DATA);
    set(handles.sourcetext1,'string',handles.data1);
    set(handles.codetext1,'string',codonify(handles.data1,1,handles.mode));
    if get(handles.codetext1,'value')>=floor(length(get(handles.sourcetext1,'string'))/3);
        set(handles.codetext1,'value',1);
    end;
    handles.result=correlate_codons(handles);
    set(handles.resultslistbox,'string',handles.result);
    
    len=length(get(handles.codetext1,'string'));
    if len<length(get(handles.codetext2,'string')); len=length(get(handles.codetext2,'string')); end;
    set(handles.matchqualitytext,'string',...
        [num2str(sum(strcmp(upper(handles.result),'GOOD'))/len*100)...
        '% match between sequences']) %#ok<*STCI>
    guidata(hObject,handles);
end;

function presetsmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in presetsmenu2.
function presetsmenu2_Callback(hObject, eventdata, handles)
if get(hObject,'Value')>1;
    DATA=code_database(get(hObject,'Value'));
    handles.data2=formatcode(DATA);
    set(handles.sourcetext2,'string',handles.data2);
    set(handles.codetext2,'string',codonify(handles.data2,1,handles.mode));
    if get(handles.codetext2,'value')>=floor(length(get(handles.sourcetext2,'string'))/3);
        set(handles.codetext2,'value',1);
    end;
    handles.result=correlate_codons(handles);
    set(handles.resultslistbox,'string',handles.result);
    
    len=length(get(handles.codetext1,'string'));
    if len<length(get(handles.codetext2,'string')); len=length(get(handles.codetext2,'string')); end;
    set(handles.matchqualitytext,'string',...
        [num2str(sum(strcmp(upper(handles.result),'GOOD'))/len*100)...
        '% match between sequences']) %#ok<*STCI>
    guidata(hObject,handles);
end;

function presetsmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in inframebutton.
function inframebutton_Callback(hObject, eventdata, handles)

% --- Executes on button press in sequence1align.
function sequence1align_Callback(hObject, eventdata, handles)

% --- Executes on button press in sequence2align.
function sequence2align_Callback(hObject, eventdata, handles)

% --- Executes on button press in alignDNA.
function alignDNA_Callback(hObject, eventdata, handles)
if get(handles.sequence1align,'value')==1;
    required=str2double(get(handles.requiredfittext,'string'));
    [handles.data1 handles.data2]=align_DNA(handles.data1,handles.data2,required);
    set(handles.presetsmenu2,'value',1)
else
    required=str2double(get(handles.requiredfittext,'string'));
    [handles.data2 handles.data1]=align_DNA(handles.data2,handles.data1,required);
    set(handles.presetsmenu1,'value',1)
end;
%Addcode takes care of all the technical stuff.
    handles.data2=formatcode(handles.data2);
    set(handles.sourcetext2,'string',handles.data2);
    set(handles.codetext2,'string',codonify(handles.data2,1,handles.mode));
    set(handles.codetext2,'value',1);
    handles.result=correlate_codons(handles);
    set(handles.resultslistbox,'string',handles.result);

    handles.data1=formatcode(handles.data1);
    set(handles.sourcetext1,'string',handles.data1);
    set(handles.codetext1,'string',codonify(handles.data1,1,handles.mode));
    set(handles.codetext1,'value',1);
    handles.result=correlate_codons(handles);
    set(handles.resultslistbox,'string',handles.result);

    len=length(get(handles.codetext1,'string'));
    if len<length(get(handles.codetext2,'string')); len=length(get(handles.codetext2,'string')); end;
    set(handles.matchqualitytext,'string',...
        [num2str(sum(strcmp(upper(handles.result),'GOOD'))/len*100)...
        '% match between sequences']) %#ok<*STCI>

guidata(hObject,handles);

function requiredfittext_Callback(hObject, eventdata, handles)
try %Reset if the user enters invalid rubbish.
    str2double(get(handles.requiredfittext,'string'))
catch %#ok<*CTCH>
    set(handles.requiredfittext,'string','3.125')
end;

% --- Executes during object creation, after setting all properties.
function requiredfittext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function varargout = digital_filters(varargin)
% DIGITAL_FILTERS MATLAB code for digital_filters.fig
%      DIGITAL_FILTERS, by itself, creates a new DIGITAL_FILTERS or raises the existing
%      singleton*.
%
%      H = DIGITAL_FILTERS returns the handle to a new DIGITAL_FILTERS or the handle to
%      the existing singleton*.
%
%      DIGITAL_FILTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIGITAL_FILTERS.M with the given input arguments.
%
%      DIGITAL_FILTERS('Property','Value',...) creates a new DIGITAL_FILTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before digital_filters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to digital_filters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help digital_filters


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @digital_filters_OpeningFcn, ...
                   'gui_OutputFcn',  @digital_filters_OutputFcn, ...
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


% --- Executes just before digital_filters is made visible.
function digital_filters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to digital_filters (see VARARGIN)

%%%%%%%%%%%%%%%%%%%%%%%%%% //Global Variables// %%%%%%%%%%%%%%%%%%%%%%%%%%
global Data_Set;
global y_fft;
%%%%%%%%%%%%%%%%%%%%%%%%%% //Draw unit Circle// %%%%%%%%%%%%%%%%%%%%%%%%%%
c_draw(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%% //Browsing For Data// %%%%%%%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile();
filename = fullfile(PathName,FileName);
Data_Set = importdata(filename);
%plot the original signal in time domain
axes(handles.axes4)
plot(Data_Set(1,:))
xlabel('Signal in T before filteration')

%plot the original signal in freq domain
y_fft=fft(Data_Set(1,:));
axes(handles.axes6)
plot(abs(y_fft))
xlabel('Signal in F before filteration')
%%%%%%%%%%%%%%%%%%%%%%%%%% //End Browsing// %%%%%%%%%%%%%%%%%%%%%%%%%%

handles.p=[];   % Array holds poles in complex formula
handles.z=[];   % Array holds zeros in complex formula

handles.A=[];   % Array holds poles points in the listBox in string formula
handles.B=[];   % Array holds zeros points in the listBox in string formula

% Choose default command line output for digital_filters
handles.output = hObject;
%%%%%%%%%%%%%%%%%%%%%%%%%% //Live Sound// %%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.axes12);
handles.state = 0;
handles.Fs = 8000;

global nBits;
nBits = 24;
global recObj;

recObj = audiorecorder(handles.Fs,nBits,1);
set(recObj,'TimerPeriod',0.05,'TimerFcn',{@audioTimerCallback,handles});

xlabel(handles.axes12,'Time');
ylabel(handles.axes12, 'Amplitude');

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = digital_filters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button Add Poles
function PP_Callback(hObject, eventdata, handles)
% hObject    handle to PP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to add a pole point in the unit circle
[x,y]=ginput(1);

%push a point and its conjugate to the poles array 
handles.p(length(handles.p)+1)=x+1j*y;
% handles.p(length(handles.p)+1)=x+1j*(-y);

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);
%push a point and its conjugate to the poles listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
% handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);


% --- Executes on button Add Zeros
function ZZ_Callback(hObject, eventdata, handles)
% hObject    handle to ZZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to add a zero point in the unit circle

[x,y]=ginput(1);
%push a point and its conjugate to the zeros array 
handles.z(length(handles.z)+1)=x+1j*y;
% handles.z(length(handles.z)+1)=x+1j*(-y);

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%push a point and its conjugate to the zeros listbox array 
handles.B{length(handles.B)+1}=['(',num2str(x),' , ',num2str(y),')'];
% handles.B{length(handles.B)+1}=['(',num2str(x),' , ',num2str(-y),')'];

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button Clear All.
function nA_Callback(hObject, eventdata, handles)
% hObject    handle to nA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear all arrays
handles.p=[];
handles.z=[];
handles.A=[];
handles.B=[];

%clear all plots
cla(handles.axes1,'reset')
cla(handles.axes2,'reset')
cla(handles.axes3,'reset')
cla(handles.axes4,'reset')
cla(handles.axes5,'reset')
cla(handles.axes6,'reset')
cla(handles.axes7,'reset')
cla(handles.axes12,'reset')


c_draw(hObject, eventdata, handles)

%update the listBoxes after clearing all
set(handles.listbox1,'String',handles.p);
set(handles.listbox2,'String',handles.z);
% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%% //Draw unit Circle// %%%%%%%%%%%%%%%%%%%%%%%%%%
function c_draw(hObject, eventdata, handles)
%draw unit circle
circle_1 = exp(1i*(0:63)*2*pi/64); 
axes(handles.axes1)
% datacursormode on
plot(real(circle_1),imag(circle_1),'.');
axis([-2 2 -2 2]); 
axis('equal'); 
hold on
grid on
plot( [0 0], [1.5 -1.5], '-')
plot( [1.5 -1.5], [0 0], '-')
xlim([-1.5 1.5])
ylim([-1.5 1.5]) 
hold off;
% Update handles structure
guidata(hObject, handles);

function freq_plot(hObject, eventdata, handles)

%clear the unit circuit axes
cla(handles.axes1,'reset');
axes(handles.axes1)
c_draw(hObject, eventdata, handles);
hold on

%plot poles and zeros markers
plot_p=plot(real(handles.p),imag(handles.p),'X');
plot_z=plot(real(handles.z),imag(handles.z),'O');
set(plot_p,'markersize',8,'linewidth',2);
set(plot_z,'markersize',8,'linewidth',2);
hold off;

%Get the transfer function coeffecients
global a;
global b;
[b,a]=zp2tf(handles.z',handles.p,1);

%Get the frequency response 
global y_fft;
[h,w] = freqz(b,a,length(y_fft));

%plot the frequency response mag 
cla(handles.axes2,'reset');

axes(handles.axes2)
plot(w/pi,20*log10(abs(h)))
xlabel('Normalized Frequency')
ylabel('Magnitude (dB)')
grid on;

%plot the frequency response phase
axes(handles.axes3)
plot(w/pi,20*log10(angle(h)))
xlabel('Normalized Frequency')
ylabel('Phase (dB)')
grid on;

%apply the filter in the orignial signal
filter=h'.*(y_fft);
axes(handles.axes5)
plot(real(ifft(filter)))
xlabel('Signal in T after filteration')

%plot the filtered signal in time and freq domains
axes(handles.axes7)
plot(abs(filter))
xlabel('fft after filteration')
xlabel('Signal in F after filteration')
  
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%% //Delete Poles&Zeros// %%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button in Delete Pole
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to btnRemovePNU_FromUse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PNUList1 = handles.listbox1;
indexedPNU = get(PNUList1,'value');
PNUnames = get(PNUList1,'String');
numResults = size(PNUnames,1);
PNUnames(indexedPNU) = [];
handles.p(indexedPNU)=[];
handles.A(indexedPNU)=[];
% and change the list sting to <empty>
if isequal(numResults,length(indexedPNU)),
    PNUnames = {'<empty>'};
    indexedPNU = 1;
end
% Ensure that list box Value is valid, then reset Value and String
indexedPNU = min(indexedPNU,size(PNUnames,1));
set(PNUList1,'String',PNUnames, 'Value', 1);
freq_plot(hObject, eventdata, handles)


% --- Executes on button press in Delete Zero
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PNUList2 = handles.listbox2;
indexedPNU2 = get(PNUList2,'value');
PNUnames2 = get(PNUList2,'String');
numResults2 = size(PNUnames2,1);
PNUnames2(indexedPNU2) = [];
handles.z(indexedPNU2)=[];
handles.B(indexedPNU2)=[];
% and change the list sting to <empty>
if isequal(numResults2,length(indexedPNU2)),
    PNUnames2 = {'<empty>'};
    indexedPNU2 = 1;
end
% Ensure that list box Value is valid, then reset Value and String
indexedPNU2 = min(indexedPNU2,size(PNUnames2,1));
set(PNUList2,'String',PNUnames2, 'Value', 1);
freq_plot(hObject, eventdata, handles)


% --- Executes on button press in Remove Poles
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to remove a pole point in the unit circle
[x,y]=ginput(1);

%search for the selected pole from the poles array and remove it
% temp=find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) );
temp=find(real(handles.p <(x+0.1)))
% handles.p(find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) ))=[];
handles.p(find(real(handles.p <(x+0.1))))=[];
%remove the selected pole from  the listBox
handles.A(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox1,'String',handles.A);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Remove Zeros
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to remove a zero point in the unit circle
[x,y]=ginput(1);

%search for the selected zero from the zeros array and remove it
% temp=find(real(handles.z <(x+0.1)) & real(handles.z>(x-0.1)) );
% handles.z(find(real(handles.z <(x+0.1)) & real(handles.z>(x-0.1)) ))=[];
temp=find(real(handles.z <(x+0.1)));
handles.z(find(real(handles.z <(x+0.1))))=[];
%remove the selected zero from  the listBox
handles.B(temp)=[];

% Update handles structure
guidata(hObject, handles);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%show the listbox after removing  points
set(handles.listbox2,'String',handles.B);

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%% //End of Delete Poles&Zeros// %%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button Browse
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data_Set;
global y_fft;
[FileName,PathName] = uigetfile();
filename = fullfile(PathName,FileName);
Data_Set = importdata(filename);
y_fft=[];
%plot the original signal in time domain
axes(handles.axes4)
plot(Data_Set(1,:))
xlabel('Signal in T before filteration')

%plot the original signal in freq domain
y_fft=fft(Data_Set(1,:));
axes(handles.axes6)
plot(abs(y_fft))
xlabel('Signal in F before filteration')

%%%%%%%%%%%%%%%%%%%%%%%%%% //Live Sound// %%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button Record/Plot
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%% Record your voice for 5 seconds.
zoom on;
recObj = audiorecorder;
%set(handles.text3,'TALK for 5 Sec')
g=msgbox('Start Talking');
n=0;
disp('Start speaking.')
recordblocking(recObj, 5);
while rand> 1e-3; 
    n=n+1;
end;
delete(g);
disp('End of Recording.');

% Store data in double-precision array.
global a;
global b;
myRecording = getaudiodata(recObj);
myRecording_fft = fft(myRecording(1,:));
[j,k] = freqz(b,a,length(myRecording_fft));
filtered=j'.*(myRecording_fft);
% Plot the waveform.
cla(handles.axes12,'reset');
axes(handles.axes12)
% plot(myRecording);
plot(real(ifft(filtered)))
title('DataPlot') %Graph Title
xlabel('Time') % x-axis label
ylabel('Amplitude') %y-axis label

%%%%%%%%%%%%%%%%%%%%%%%%%% //Move Points// %%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button move Zeros
function pushbutton19_Callback(hObject, eventdata, handles)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to move a zero point in the unit circle
[x,y]=ginput(1);

%search for the selected zero from the zero array and move it
temp=find(real(handles.z <(x+0.1)) & real(handles.z>(x-0.1)) );
handles.z(find(real(handles.z <(x+0.1)) & real(handles.z>(x-0.1)) ))=[];
%remove the selected zero from  the listBox
handles.B(temp)=[];
% [r,l]=ginput()
% x=r;
% y=l;
h = impoint.empty;
h = impoint
position = getPosition(h)
x=position(1);
y=position(2);

handles.z(length(handles.z)+1)=x+1j*y;
handles.z(length(handles.z)+1)=x+1j*(-y);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%push a point and its conjugate to the zeros listbox array 
handles.B{length(handles.B)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.B{length(handles.B)+1}=['(',num2str(x),' , ',num2str(-y),')'];

%show the listbox with all added points
set(handles.listbox2,'String',(handles.B));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button move Poles.
function pushbutton21_Callback(hObject, eventdata, handles)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%click to move a zero point in the unit circle
[x,y]=ginput(1);

%search for the selected zero from the zero array and move it
temp=find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) );
handles.p(find(real(handles.p <(x+0.1)) & real(handles.p>(x-0.1)) ))=[];
%remove the selected zero from  the listBox
handles.A(temp)=[];

h = impoint.empty;
h = impoint;
position = getPosition(h);
x=position(1);
y=position(2);

handles.p(length(handles.p)+1)=x+1j*y;
handles.p(length(handles.p)+1)=x+1j*(-y);

%plot the freq response and its effect in the original signal
freq_plot(hObject, eventdata, handles);

%push a point and its conjugate to the zeros listbox array 
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(y),')'];
handles.A{length(handles.A)+1}=['(',num2str(x),' , ',num2str(-y),')'];

%show the listbox with all added points
set(handles.listbox1,'String',(handles.A));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button Exit
function pushbutton16_Callback(hObject, eventdata, handles)
cl = questdlg('Are You Handsome?','EXIT',...
            'Yes','No','No');
switch cl
    case 'No'
        close();
        clear all;
        return;
    case 'Yes'
        quit cancel;
end 

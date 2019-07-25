function varargout = MODA(varargin)
% MODA MATLAB code for MODA.fig
%      MODA, by itself, creates a new MODA or raises the existing
%      singleton*.
%
%      H = MODA returns the handle to a new MODA or the handle to
%      the existing singleton*.
%
%      MODA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODA.M with the given input arguments.
%
%      MODA('Property','Value',...) creates a new MODA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MODA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MODA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MODA

% Last Modified by GUIDE v2.5 03-Oct-2017 10:58:42

% Begin initialization code - DO NOT EDIT
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MODA_OpeningFcn, ...
                   'gui_OutputFcn',  @MODA_OutputFcn, ...
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


% --- Executes just before MODA is made visible.
function MODA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MODA (see VARARGIN)

% Choose default command line output for MODA
handles.output = hObject;
movegui(gcf,'center')
axes(handles.logo)
matlabImage = imread('frontbanner.png');
image(matlabImage)
axis off
%axis image

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MODA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MODA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in TFA.
function TFA_Callback(~,~,~)

TimeFrequencyAnalysis


% --- Executes on button press in WPC.
function WPC_Callback(~,~,~)

CoherenceMulti


% --- Executes on button press in filt.
function filt_Callback(~,~,~)

Filtering


% --- Executes on button press in bisp.
function bisp_Callback(~,~,~)

Bispectrum


% --- Executes on button press in bayesian.
function bayesian_Callback(~,~,~)

Bayesian

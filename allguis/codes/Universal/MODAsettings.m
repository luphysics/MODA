% MODA GUI settings function

function handles=MODAsettings(hObject,handles)

% Positioning
screensize = get( groot, 'Screensize' );
x=screensize(3);
y=screensize(4);
if x<1600 || y<860
    GUIsize=[x y];
    set(gcf,'units','pixels','position',[0 0 GUIsize(1) GUIsize(2)])
else
    GUIsize=[1600 860];
    x2=x-GUIsize(1);
    y2=y-GUIsize(2);
    set(gcf,'units','pixels','position',[x2/2 y2/2 GUIsize(1) GUIsize(2)])
end

% Colours
load('cmap.mat')
handles.cmap=cmap;
handles.linecol= cmap([1,18,40,50,60,64,15],:);
handles.line2width=2;

% Logos
axes(handles.logo);
matlabImage = imread('physicslogo.png');
image(matlabImage)
axis off
axis image

axes(handles.nbmplogo);
matlabImage = imread('MODAbanner5.png');
image(matlabImage)
axis off
axis image

% Fonts
h = findall(0,'Type','uicontrol');
set(h,'FontUnits','points');
set(h,'FontSize',8);
set(h,'FontUnits','normalized');

% Default calculation types
handles.calc_type = 1;
handles.plot_type = 2;

handles.it=0; % Number of times data has been loaded




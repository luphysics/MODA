% MODA data loading check function

function [handles,A]=MODAreadcheck(handles)

handles.it=handles.it+1;
A=0;
if handles.it>1
choice = questdlg('Loading new data will erase unsaved data. Continue?', ...
        'Data Import','Yes','No','default');
   switch choice
       case 'Yes'
           A=1;
             
       case 'No'
           
           A=0;
        
   end
   
else
    A=1;
end
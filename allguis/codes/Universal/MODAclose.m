% MODA GUI closing function
function MODAclose(hObject,handles)
if isfield(handles,'sig')
choice = questdlg('Are you sure you want to close? Unsaved data will be lost.', ...
        'Exit','Yes','No','default');
    switch choice
        case 'Yes'
            if isfield(handles,'h')
                delete(handles.h)
            else
            end
            delete(hObject);            
            return;            
            
        case 'No'        
    end
    
    if isempty(choice)
         return;
    end
else
    delete(hObject);  
end

% Save current session

function MODAsave(handles)
try
[FileName,PathName] = uiputfile('.mat','Save data as');
save_location = strcat(PathName,FileName);

save(save_location,'handles');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


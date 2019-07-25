function clear_axes_lines( h )
%CLEAR_AXES_LINES Deletes lines on an axes
child_handles = allchild(h);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line')) 
            delete(child_handles(i))
    end
end


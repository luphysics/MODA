function clear_axes_points( h )
%CLEAR_AXES_POINTS Deletes points on an axes
child_handles = allchild(h);
for i = 1:size(child_handles,1)    
    if(strcmp(get(child_handles(i),'Type'),'line'))
        xdat = get(child_handles(i),'XData');        
        if(length(xdat) == 1)
            delete(child_handles(i))
        end
    elseif(strcmp(get(child_handles(i),'Type'),'text'))
        delete(child_handles(i))
    end
end


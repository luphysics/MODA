function handles=MODAbayes_intdelete(hObject,eventdata,handles)
% Function to deleta data from handles

%switch eventdata.Key
    %case 'delete' % Excutes when user presses delete
        interval_selected = get(handles.interval_list_1,'Value');
        n=1:handles.c;
        
        ne=n(1:end ~=interval_selected);
        if isfield(handles,'pinput')
        else
        handles.int1=handles.int1(ne,:); % Delete data
        handles.int2=handles.int2(ne,:); 
        end
        handles.winds=handles.winds(ne);
        handles.pr=handles.pr(ne);
        handles.ovr=handles.ovr(ne);
        handles.forder=handles.forder(ne);
        handles.ns=handles.ns(ne);
        handles.confidence_level=handles.confidence_level(ne);
        
        if isfield(handles,'tm') && size(handles.tm,2)==handles.c           
        handles.tm=handles.tm(:,ne);
        handles.cc=handles.cc(:,ne);
        handles.e=handles.e(:,ne);
        handles.cpl1=handles.cpl1(:,ne);
        handles.cpl2=handles.cpl2(:,ne);
        handles.cf1=handles.cf1(:,ne);
        handles.cf2=handles.cf2(:,ne);
        handles.mcf1=handles.mcf1(:,ne);
        handles.mcf2=handles.mcf2(:,ne);
        handles.surr_cpl1=handles.surr_cpl1(:,ne);
        handles.surr_cpl2=handles.surr_cpl2(:,ne);
        handles.p1=handles.p1(:,ne);
        handles.p2=handles.p2(:,ne);
        else
        end
        
            
        if min(interval_selected)>1
            set(handles.interval_list_1,'Value',min(interval_selected)-1);
        else
            set(handles.interval_list_1,'Value',1);
        end
        
        list = get(handles.interval_list_1,'String');
        list(interval_selected,:) = [];
        set(handles.interval_list_1,'String',list);
        
        interval_selected = get(handles.interval_list_2,'Value');
        if min(interval_selected)>1
            set(handles.interval_list_2,'Value',min(interval_selected)-1);
        else
            set(handles.interval_list_2,'Value',1);
        end
        list = get(handles.interval_list_2,'String');
        list(interval_selected,:) = [];
        set(handles.interval_list_2,'String',list);
        if isfield(handles,'bands')
            handles.bands(:,interval_selected) = [];
            guidata(hObject,handles);
        else            
        end        
        drawnow;
        handles.c=handles.c-1;
        guidata(hObject,handles);
        
%         handles.int1(interval_selected,:)=[]; % Delete data
%         handles.int2(interval_selected,:)=[];       
%         handles.winds(interval_selected,:)=[];
%         handles.pr(interval_selected,:)=[];
%         handles.ovr(interval_selected,:)=[];
%         handles.forder(interval_selected,:)=[];
%         handles.ns(interval_selected,:)=[];           
%             
%         if min(interval_selected)>1
%             set(handles.interval_list_1,'Value',min(interval_selected)-1);
%         else
%             set(handles.interval_list_1,'Value',1);
%         end
%         
%         list = get(handles.interval_list_1,'String');
%         list(interval_selected,:) = [];
%         set(handles.interval_list_1,'String',list);
%         
%         interval_selected = get(handles.interval_list_2,'Value');
%         if min(interval_selected)>1
%             set(handles.interval_list_2,'Value',min(interval_selected)-1);
%         else
%             set(handles.interval_list_2,'Value',1);
%         end
%         list = get(handles.interval_list_2,'String');
%         list(interval_selected,:) = [];
%         set(handles.interval_list_2,'String',list);
%         if isfield(handles,'bands')
%             handles.bands(:,interval_selected) = [];
%             guidata(hObject,handles);
%         else            
%         end        
%         drawnow;
%         
%         guidata(hObject,handles);
%end
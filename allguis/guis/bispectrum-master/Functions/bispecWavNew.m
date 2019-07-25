% update for MODA 20.07.2018
% code compatible with Dmytro Iatsenko's wt.m
% author: Aleksandra Pidde a.pidde@gmail.com, a.pidde@lancaster.ac.uk

function [Bisp, freq, opt,wt1,wt2] = bispecWavNew(sig1, sig2, fs, varargin)
% function calculating wavelet bispectrum
%
% INPUT:
% sig1, sig2: 	signals
% fs: 			sampling frequency
% optional:
% opt:			structure of optimal parameters returned by wt.m
%
% OUTPUT:
% Bisp: 		bispectrum, 2 dim matrix of complex values
% freq: 		frequencies
% opt:			structure of optimal parameters returned by wt.m
try
    
bstype={'111';'222';'122';'211'};
if nargin == 4
	opt = varargin{1};
	[wt1, freq, opt] = wt(sig1, fs, opt);
	[wt2, freq, opt] = wt(sig2, fs, opt);
else
	[wt1, freq, opt] = wt(sig1, fs);
	[wt2, freq, opt] = wt(sig2, fs);
end

dt = 1 / fs;
nfreq = length(freq);
Bisp = NaN * zeros(nfreq, nfreq);
auto = false;
if compareMatrix(wt1, wt2)
    auto = true;
end
wbar=0;

if nargin >= 2 + 3
    for i = 1 : 2 : nargin - 3
        switch varargin{i}
            
            case 'handles'
                handles=varargin{i+1};
            case 'hObject'
                hObject=varargin{i+1};
            case 'num'
                numb=varargin{i+1};
            case 'wbar'
                wbar=varargin{i+1};
            
        end
    end
end

if wbar==1
handles.h = waitbar(0,'Calculating bispectrum...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles); 
else
end

for j = 1 : nfreq
    if wbar==1
    if getappdata(handles.h,'canceling')
        guidata(hObject,handles);
                break;
    end
    else
    end
    kstart = 1;
    if auto
        kstart = j;
    end

        if wbar==1
        if getappdata(handles.h,'canceling')
            guidata(hObject,handles);
                break;
        end
        else
        end
    for k = kstart : nfreq
        f3 = freq(j) + freq(k);
        bigger = max([j k]);
        idx3 = find(freq >= f3, 1);
        if (f3 <= freq(end)) && (freq(idx3 - 1) > freq(bigger)) 
            WTdat = wtAtf2(sig2, fs, f3, opt);
            WTdat = WTdat(:).'; % make sure it is horizontal vector
            xx = wt1(j, :) .* wt2(k, :) .* conj(WTdat);
            %xx = TFR1(j, :) .* TFR2(k, :) * transpose(conj(TFR2(idx3, :)));
            ss = nanmean(xx);
            Bisp(j, k) = ss;
        end
    end
    if wbar==1
    waitbar(j / length(freq),handles.h,sprintf(['Calculating Bispectrum ',bstype{numb},' (%d/%d)'],j,length(freq)));
    else
    end
end
if wbar==1
delete(handles.h); 
else
end

catch e
    errordlg(e.message,'Error'); 
    if wbar==1
    delete(handles.h); 
    else
    end
    rethrow(e)

end



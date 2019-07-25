function w = modwt(x,varargin)
%MODWT Maximal overlap discrete wavelet transform.
%   W = MODWT(X) computes the maximal overlap discrete wavelet transform of
%   a 1-D real-valued, double-precision input signal, X. The signal can be
%   a row or column vector and must contain at least two samples. By
%   default, the maximal overlap discrete wavelet transform is computed
%   down to level floor(log2(length(X))) using the Daubechies
%   least-asymmetric wavelet with 4 vanishing moments ('sym4') and periodic
%   boundary handling. W is a LEV+1-by-N matrix containing the wavelet
%   coefficients and final-level scaling coefficients. LEV is the level of
%   the MODWT. The m-th row of the matrix contains the wavelet (detail)
%   coefficients for scale 2^m. The LEV+1-th row of the matrix contains the
%   scaling coefficients for scale 2^LEV.
%
%   W = MODWT(X,WNAME) computes the MODWT using the wavelet, WNAME. WNAME
%   is a character vector denoting the name of an orthogonal wavelet.
%   Orthogonal wavelets are designated as type 1 wavelets in the wavelet
%   manager. Valid built-in orthogonal wavelet families begin with
%   'haar','dbN', 'fkN', 'coifN', or 'symN' where N is the number of
%   vanishing moments for all families except 'fk'. For 'fk', N is the
%   number of filter coefficients. You can determine valid values for N by
%   using waveinfo. For example, waveinfo('db'). You can check if your
%   wavelet is orthogonal by using wavemngr('type',wname) to see if a 1 is
%   returned. For example, wavemngr('type','db2').
%
%   W = MODWT(X,Lo,Hi) computes the maximal overlap discrete wavelet
%   transform using the scaling filter, Lo, and the wavelet filter, Hi. Lo
%   and Hi are even-length row or column vectors. These filters must
%   satisfy the conditions for an orthogonal wavelet. You cannot specify
%   both WNAME and a filter pair, Lo and Hi. 
%
%   W = MODWT(...,LEV) computes the maximal overlap discrete wavelet
%   transform down to the level, LEV. LEV is a positive integer that cannot
%   exceed floor(log2(length(X))). If unspecified, LEV defaults to
%   floor(log2(length(X))).
%
%   W = MODWT(...,'reflection') uses reflection boundary handling by
%   extending the signal symmetrically at the right boundary to twice the
%   signal length,[x flip(x)], before computing the wavelet transform. The
%   number of wavelet and scaling coefficients returned are twice the
%   length of the input signal. By default, the signal is extended
%   periodically. You must enter the entire character vector 'reflection'.
%   If you added a wavelet named 'reflection' using the wavelet manager,
%   you must rename that wavelet prior to using this option. 'reflection'
%   may be placed in any position in the input argument list after X.
%
%
%   % Example 1:
%   %   Obtain the maximal overlap discrete wavelet transform of the Nile 
%   %   river minimum water level data. The data is 663 samples in length 
%   %   sampled yearly. Use the Haar wavelet and transform the data down 
%   %   to level 8. Plot the level-3 wavelet coefficients.
%    
%   load nileriverminima;
%   w = modwt(nileriverminima,'haar',8);
%   plot(w(3,:)); title('Level-3 Wavelet Coefficients');
%
%   % Example 2:
%   %   Check that the maximal overlap discrete wavelet transform 
%   %   partitions the variance of the signal by scale.
%
%   load noisdopp;
%   [~,~,Lo,Hi] = wfilters('sym8');
%   w = modwt(noisdopp,Lo,Hi,10);
%   varbylev = var(w,1,2);
%   sum(varbylev)
%   var(noisdopp,1)
%
%   See also IMODWT MODWTMRA MODWTCORR MODWTXCORR


% Check number of input arguments
narginchk(1,5);

% Validate that data is real, 1-D double-precision
% with no NaNs or Infs

validateattributes(x,{'double'},{'real','nonnan','finite'});

    

%Input must be 1-D
if (~isrow(x) && ~iscolumn(x))
    error(message('Wavelet:modwt:OneD_Input'));
end

%Input must contain at least two samples
if (numel(x)<2)
    error(message('Wavelet:modwt:LenTwo'));
end


% Convert data to row vector
x = x(:)';
% Record original data length
datalength = length(x);

%Parse input arguments
params = parseinputs(datalength,varargin{:});

%Check that the level of the transform does not exceed floor(log2(numel(x))
J = params.J;
Jmax = floor(log2(datalength));
if (J <= 0) || (J > Jmax) || (J ~= fix(J))
    error(message('Wavelet:modwt:MRALevel'));
end

boundary = params.boundary;
if (~isempty(boundary) && ~strcmpi(boundary,'reflection'))
    error(message('Wavelet:modwt:Invalid_Boundary'));
end

% increase signal length if 'reflection' is specified
if strcmpi(boundary,'reflection')
    x = [x flip(x)];
end

% obtain new signal length if needed
siglen = length(x);
Nrep = siglen;


% If wavelet specified as a string, ensure that wavelet is orthogonal
if (isfield(params,'wname') && ~isfield(params,'Lo'))
    [~,~,Lo,Hi] = wfilters(params.wname);
    wtype = wavemngr('type',params.wname);
    if (wtype ~= 1)
        error(message('Wavelet:modwt:Orth_Filt'));
    end
end

%If scaling and wavelet filters are specified as vectors, ensure they
%satisfy the orthogonality conditions

if (isfield(params,'Lo') && ~isfield(params,'wname'))
    filtercheck = CheckFilter(params.Lo,params.Hi);
    if ~filtercheck
        error(message('Wavelet:modwt:Orth_Filt'));
    end
    Lo = params.Lo;
    Hi = params.Hi;
end

% Scale the scaling and wavelet filters for the MODWT
Lo = Lo./sqrt(2);
Hi = Hi./sqrt(2);

% Ensure Lo and Hi are row vectors
Lo = Lo(:)';
Hi = Hi(:)';

% If the signal length is less than the filter length, need to 
% periodize the signal in order to use the DFT algorithm

if (siglen<numel(Lo))
    x = [x repmat(x,1,ceil(numel(Lo)-siglen))];
    Nrep = numel(x);
end

% Allocate coefficient array
w = zeros(J+1,Nrep);

% Obtain the DFT of the filters
G = fft(Lo,Nrep);
H = fft(Hi,Nrep);

%Obtain the DFT of the data
Vhat = fft(x);

% Main MODWT algorithm
for jj = 1:J
    [Vhat,What] = modwtdec(Vhat,G,H,jj);
    w(jj,:) = ifft(What);
end

w(J+1,:) = ifft(Vhat);

% Truncate data to length of boundary condition
w = w(:,1:siglen);


%----------------------------------------------------------------------
function [Vhat,What] = modwtdec(X,G,H,J)
% [Vhat,What] = modwtfft(X,G,H,J)

N = length(X);
upfactor = 2^(J-1);
Gup = G(1+mod(upfactor*(0:N-1),N));
Hup = H(1+mod(upfactor*(0:N-1),N));
Vhat = Gup.*X;
What = Hup.*X;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = CheckFilter(Lo,Hi)
% For a user-supplied scaling and wavelet filter, check that
% both correspond to an orthogonal wavelet
Lo = Lo(:);
Hi = Hi(:);
Lscaling = length(Lo);
Lwavelet = length(Hi);
evenlengthLo = 1-rem(Lscaling,2);
evenlengthHi = 1-rem(Lwavelet,2);
if all([evenlengthLo evenlengthHi])
    evenlength = 1;
else
    evenlength = 0;
end
 
if (Lscaling ~= Lwavelet)
    equallen = 0;
else
    equallen = 1;
end

normLo = norm(Lo,2);
sumLo = sum(Lo);
normHi = norm(Hi,2);
sumHi = sum(Hi);
tol = 1e-7;

if (abs(normLo-1) > tol && abs(normHi -1) > tol)
    unitnorm = 0;
else
    unitnorm = 1;
     
end

if (abs(sumLo-sqrt(2)) > tol && abs(sumHi) > tol)
    sumfilters = 0;
else
    sumfilters = 1;
end


zeroevenlags = 1;

if (Lscaling > 2)
    L = Lscaling;
    xcorrHi = conv(Hi,flipud(Hi));
    xcorrLo = conv(Lo,flipud(Lo));
    xcorrLo = xcorrLo(L+2:2:end);
    xcorrHi = xcorrHi(L+2:2:end);
    zeroevenlagsLo = 1-any(abs(xcorrLo>tol));
    zeroevenlagsHi = 1-any(abs(xcorrHi>tol));
    zeroevenlags = max(zeroevenlagsLo,zeroevenlagsHi);
end
out = all([evenlength equallen unitnorm sumfilters zeroevenlags]);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = parseinputs(siglen,varargin)
 % Parse varargin and check for valid inputs
 
% Assign defaults 
params.boundary = [];
params.J = floor(log2(siglen));
params.wname = 'sym4';
 
% Check for 'reflection' boundary      
tfbound = strcmpi(varargin,'reflection');
 
% Determine if 'reflection' boundary is specified
if any(tfbound)
    params.boundary = varargin{tfbound>0};
    varargin(tfbound>0) = [];
end
 
 % If boundary is the only input in addition to the data, return with
 % defaults 
if isempty(varargin)
    return;
end
   
 % Only remaining char variable must be wavelet name
tfchar = cellfun(@ischar,varargin);
if (nnz(tfchar) == 1)
    params.wname = varargin{tfchar>0};
end
        
% Only scalar input must be the level
tfscalar = cellfun(@isscalar,varargin); 
 
% Check for numeric inputs 
tffilters = cellfun(@isnumeric,varargin);

% At most 3 numeric inputs are supported
if nnz(tffilters)>3
    error(message('Wavelet:modwt:Invalid_Numeric'));
end

% There's one numeric argument and it's not a wavelet level
if (nnz(tffilters)==1) && (nnz(tfscalar) == 0)
    error(message('Wavelet:FunctionInput:InvalidLoHiFilters'));
end

% If there are at least two numeric inputs, the first two must be the
% scaling and wavelet filters
if (nnz(tffilters)>1)
    idxFilt = find(tffilters,2,'first');
    params.Lo = varargin{idxFilt(1)};
    params.Hi = varargin{idxFilt(2)};
    params = rmfield(params,'wname');
    
    if (length(params.Lo) < 2 || length(params.Hi) < 2)
        error(message('Wavelet:modwt:Invalid_Filt_Length'));
    end
    
end
 
% Any scalar input must be the level
if any(tfscalar)
    params.J = varargin{tfscalar>0};
end
 
% If the user specifies a filter, use that instead of default wavelet
if (isfield(params,'Lo') && any(tfchar))
     error(message('Wavelet:FunctionInput:InvalidWavFilter'));
end
 
 
 
 
 
 

 
 
 
 
     
 
  
    
    
    
    
     
       

    
         
        



% [EOF] modwt.m

    

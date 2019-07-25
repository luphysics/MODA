function xrec = imodwt(w,varargin)
%IMODWT Inverse maximal overlap discrete wavelet transform.
%   XREC = IMODWT(W) returns the reconstructed signal based on the maximal
%   overlap discrete wavelet transform (MODWT) coefficients in W. W is a
%   LEV+1-by-N matrix which is the MODWT of an N-point input signal down to
%   level LEV. By default, IMODWT assumes that you used the 'sym4' wavelet
%   with periodic boundary handling to obtain the MODWT. If you do not
%   modify the coefficients, XREC is a perfect reconstruction of the
%   signal.
%
%   XREC = IMODWT(W,WNAME) reconstructs the signal using the wavelet WNAME.
%   WNAME must be the same wavelet used in the analysis of the signal with
%   MODWT.
%
%   XREC = IMODWT(W,Lo,Hi) reconstructs the signal using the scaling
%   filter Lo and the wavelet filter Hi. You cannot specify both WNAME and
%   Lo and Hi. Lo and Hi must be the same filters used in the analysis with
%   MODWT.
%
%   XREC = IMODWT(...,LEV) reconstructs the signal up to level LEV. XREC is
%   a projection onto the scaling space at level LEV. LEV is a nonnegative
%   integer between 0 and strictly less than size(W,1)-1. The default is
%   LEV = 0, which results in perfect reconstruction if you do not modify
%   the coefficients. 
%
%   XREC = IMODWT(...,'reflection') uses the 'reflection' boundary
%   condition in the reconstruction. If you specify 'reflection', IMODWT
%   assumes that the length of the original signal is 1/2 the number of
%   columns in the input coefficient matrix W. You must enter the entire
%   string 'reflection'. If you added a wavelet named 'reflection' using
%   the wavelet manager, you must rename that wavelet prior to using this
%   option. 'reflection' may be placed in any position in the input
%   argument list after W. By default both MODWT and IMODWT assume periodic
%   signal extension at the boundary.
%
%   %Example 1:
%   %   Demonstrate perfect reconstruction of ECG data with the MODWT.
%
%   load wecg;
%   wecg = wecg(1:end-1);
%   w = modwt(wecg,'sym4',10);
%   xrec = imodwt(w);
%   max(abs(xrec-wecg'))
%   subplot(2,1,1)
%   plot(wecg); title('Original Signal');
%   subplot(2,1,2)
%   plot(xrec); title('Reconstructed Signal');
%
%   %EXAMPLE 2:
%   %   Reconstruct a signal approximation based on the level-3 and
%   %   level-4 wavelet coefficients.
%
%   load wecg;
%   wecg = wecg(1:end-1);
%   w = modwt(wecg,'db2',10);
%   idx = 3:4;
%   wnew = zeros(size(w));
%   wnew = w(idx,:);
%   xrec = imodwt(wnew);
%   subplot(211)
%   plot(wecg); title('Original Signal');
%   subplot(212)
%   plot(xrec);
%   title('Reconstruction from Level-3 and Level-4 Wavelet Coefficients');
%
%   See also modwt, modwtmra, modwtvar, modwtcorr, modwtxcorr

% Check number of input arguments
narginchk(1,5);

% Input cannot be a row or column vector. IMODWT expects at least a two row
% matrix
if (isrow(w) || iscolumn(w))
    error(message('Wavelet:modwt:InvalidCFSSize'));
end

% Input must be real-value, finite, and double precision
validateattributes(w,{'double'},{'real','nonnan','finite'});

% Parse input arguments
params = parseinputs(varargin{:});

% Get the original input size
% Get the level of the MODWT
N = size(w,2);
Nrep = N;
J = size(w,1)-1;


boundary = params.boundary;
if (~isempty(boundary) && ~strcmpi(boundary,'reflection'))
    error(message('Wavelet:modwt:Invalid_Boundary'));
end

% Adjust final output length if MODWT obtained with 'reflection'
if strcmpi(boundary,'reflection')
    N = N/2;
end

% If the wavelet is specified as a string, obtain filters from wavemngr
if (isfield(params,'wname') && ~isfield(params,'Lo'))
    [~,~,Lo,Hi] = wfilters(params.wname);
    wtype = wavemngr('type',params.wname);
    if (wtype ~= 1)
        error(message('Wavelet:modwt:Orth_Filt'));
    end
end

%If scaling and wavelet filters specified as vectors, ensure they
%satisfy the orthogonality conditions

if (isfield(params,'Lo') && ~isfield(params,'wname'))
    filtercheck = CheckFilter(params.Lo,params.Hi);
    if ~filtercheck
        error(message('Wavelet:modwt:Orth_Filt'));
    end
    Lo = params.Lo;
    Hi = params.Hi;
end


% Scale scaling and wavelet filters for MODWT
Lo = Lo./sqrt(2);
Hi = Hi./sqrt(2);


% If the number of samples is less than the length of the scaling filter
% we have to periodize the data and then truncate.
if (Nrep<numel(Lo))
    w = [w repmat(w,1,ceil(numel(Lo)-Nrep))];
    Nrep = size(w,2);
end

% Get the DFTs of the scaling and wavelet filters
G = fft(Lo,Nrep);
H = fft(Hi,Nrep);

% Get the level of the reconstruction
lev = params.lev+1;

% Error
if (lev>J)
    error(message('Wavelet:modwt:Incorrect_ReconLevel'));
end

vin = w(J+1,:);

% IMODWT algorithm

for jj = J:-1:lev
    vout = imodwtrec(vin,w(jj,:),G,H,jj);
    vin = vout;
end

% Return proper output length
xrec = vout(1:N);




%-----------------------------------------------------------------
function Vout = imodwtrec(Vin,Win,G,H,J)
N = length(Vin);
Vhat = fft(Vin);
What = fft(Win);
upfactor = 2^(J-1);
Gup = conj(G(1+mod(upfactor*(0:N-1),N)));
Hup = conj(H(1+mod(upfactor*(0:N-1),N)));
Vout = ifft(Gup.*Vhat+Hup.*What);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = parseinputs(varargin)
% Parse varargin and check for valid inputs

% Assign default inputs
params.boundary = [];
params.lev = 0;
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

% If there are at least two numeric inputs, the first two must be the
% scaling and wavelet filters

if (nnz(tffilters)>1)
    idxFilt = find(tffilters,2,'first');
    params.Lo = varargin{idxFilt(1)};
    params.Hi = varargin{idxFilt(2)};
    if (length(params.Lo) < 2 || length(params.Hi) < 2)
        error(message('Wavelet:modwt:Invalid_Filt_Length'));
    end
    params = rmfield(params,'wname');
end


% Any scalar input must be the level
if any(tfscalar)
    params.lev = varargin{tfscalar>0};
end

% If the user specifies a filter, use that instead of default wavelet
if (isfield(params,'Lo') && any(tfchar))
    error(message('Wavelet:FunctionInput:InvalidWavFilter'));
end

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

if (abs(sumLo-sqrt(2))> tol && abs(sumHi) > tol)
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










% [EOF] imodwt.m










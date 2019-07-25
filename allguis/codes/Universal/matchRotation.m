function sig = matchRotation(original, s)
% function to support wavsurrogates
% author: Aleksandra Pidde, a.eksandra.pidde@gmail.com
rev = s(end:-1:1);
[sig1, err1] = matchRot(original, s);
[sig2, err2] = matchRot(original, rev);
if err1 < err2
    sig = sig1;
else
    sig = sig2;
end
end

function [sig, err] = matchRot(original, s)
N = length(s);
err = sum((original - s).^2);
sig = s;
[~, dim]  = max(size(s));
for i = 2 : N - 1
   tmp = circshift(s, 1, dim);
   er = sum((original - tmp).^2);
   if er < err
       err = er;
       sig = tmp;
   end
end
end

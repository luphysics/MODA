% update for MODA 19.07.2018
% author: Aleksandra Pidde a.pidde@gmail.com, a.pidde@lancaster.ac.uk

function res = compareMatrix(A, B)
% function helper compating 2 matrix, A and B
nanA = isnan(A); nanB = isnan(B);
res = isequal(nanA, nanB);

infA = isinf(A); infB = isinf(B);
res = res * isequal(infA, infB);

idxA = nanA + infA; idxB = nanB + infB;
res = res * isequal(A(~idxA), B(~idxB));
end

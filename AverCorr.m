function [RlalaM, PlalaM, meanCorr] = AverCorr(A, B)
    %A and B matrxes of individ measurement
    [m,n] = size(A);
    for ii=1:n
        [Rlala,Plala] = corrcoef(A(find(A(:,ii)~=0), ii), B(find(B(:,ii)~=0), ii));
        RlalaM(ii) = Rlala(2,1);
        PlalaM(ii) = Plala(2,1);
    end;
    meanCorr = mean(RlalaM);
end

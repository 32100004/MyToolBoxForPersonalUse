function m=ELM_endmmeber_estimation(Y)

% Input: Y: Hyperspectral Data containing L bands N Pixels
% OutPut: m: Estimated Endmember Numbers
% Reference:  Ben Ismail & Ouiem Bchir
%             Survey on Number of Endmembers Estimation Techniques for Hyperspectral Data Unmixing

[L,N] = size(Y);

% Compute the covariance matrix

K = 1/N*bsxfun(@minus,Y,mean(Y))*transpose(bsxfun(@minus,Y,mean(Y)));
%

% Compute the correlation matrix using Eq (2)

R = 1/N*(Y*Y.');

% - Extract the eigenvalues ofthe covariance matrix T1

eigenvalue_K = sort(eig(K),'descend');

% - Extract the eigenvalues ofthe correlation matrix T2

eigenvalue_R = sort(eig(R),'descend');

% - Compute z = Tl- T2

z = eigenvalue_K - eigenvalue_R;

SigmaArray = sqrt(2/N*( eigenvalue_K.^2 + eigenvalue_R.^2));

G=zeros(1,L);
for i = L : -1 : 1
    LFGVariable = -0.5*(z(i)/SigmaArray(i))^2-log(SigmaArray(i));
    if i==L
        G(i)=LFGVariable;
    else
        G(i) = G(i+1) + LFGVariable;
    end
end

m = find( G == max(G) ) - 1;


end

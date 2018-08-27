function DeltaE=MyEstNoise(Y)

% follow the article "Modified Residual Method for the Estimation of Noise in Hyperspectral Images"   
% slightly different from estNoise of Hysime

        X = Y';
        
        [N,L] = size(X);
        
        % X is assumed to be centered and hence has zero mean for all the bands.
        
        Xc=bsxfun(@minus,X,mean(X));
        
        X=Xc;
        
        DeltaE=zeros(L);
        
        % covariance matrix of centered data C = (1/(N ?1))XT X
        
        C = X'*X/(N-1);
        
        for i = 1 : L
            Xd = X;
            Xd(:,i) = [];
            
            Ct=C;
            Ct(i,:)=[];
            Ct(:,i)=[];
            
            C0 = C(:,i);
            C0(i)=[];
            
            Palpha=Ct\C0;
            
            ri_star = X(:,i) - Xd*Palpha;
            
            DeltaE(i,i)=ri_star'*ri_star/N;
            
        end
        
        
    end
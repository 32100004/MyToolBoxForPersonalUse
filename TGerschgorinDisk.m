function m = TGerschgorinDisk(Y,choice)

% Noise_adjusted TGD method does not test

if nargin ==1
    choice = 'Traditional'
end

% Compute Covariance Matrix
[L,N]=size(Y);

%R = 1/N*bsxfun(@minus,Y,mean(Y))*transpose(bsxfun(@minus,Y,mean(Y)));

switch choice 
    case 'Traditional'
        R=1/N*(Y*Y');
    case 'Noise_adjusted'
        [~,Rn]=estNoise(Y); % Function estNoise from Hysime
        [V_Rn,D_Rn] = eig(Rn);
        Wn=V_Rn*D_Rn^(-0.5);
        R=1/N*(Y*Y');
        R=Wn'*R*Wn;
    otherwise
        error('Methods Must be Traditional or Noise_adjusted!');
end

% Decomposition

C = R(1:end-1,1:end-1);
cl = R(1:end-1,end);
cll = R(end,end);

% SVD decomposition

[Qc,Gamma_c,~]=svd(C);

centre=diag(Gamma_c);

% Construct Q

radius=abs(cl'*Qc);

% Plot Gerschgorin Disk
figure,subplot(211),hold on

for i=1:L-1
    plot_circle(centre(i),0,radius(i));
end

axis equal

subplot(212);tmp_array=zeros(1,L-2);

for i=1:L-2
    tmp_array(i)=radius(i)*centre(i)/radius(i+1)/centre(i+1);
end

plot(tmp_array);

m=find(tmp_array==max(tmp_array));

    function plot_circle(x,y,r)
        mm=linspace(0,2*pi,1000);
        mmx=x+r*sin(mm);
        mmy=y+r*cos(mm);
        plot(mmx,mmy);
    end
end
function [S,par] = STAR(L,B,Phi,M,varargin)
% This script is to solve E/MEG source imaging problem under u-STAR
% framework
% Source Model: S = W*Phi+E
% =================================================================
% [S,par] = ESSTED(L,B,Phi,M,varargin) 
% Input:
%   -L: lead-field matrix
%   -B: E/MEG measurements
%   -Phi: temporal basis functions (TBFs)
%   -M: spatial prior
% Output:
%   -S: reconstructed source
%   -par:
%       par.Phi
%       par.W
%       par.E
%==================================================================
% Author: Feng Zhao
% Data: 2022/5/8
% Reference:
% [1] Bayesian electromagnetic spatio-temporal imaging of extended sources
%     with Markov Random Field and temporal basis expansion
% [2] Probabilistic algorithms for MEG/EEG source reconstruction using 
%     temporal basis functions learned from data.
% [3] u-STAR: a novel framework for spatio-temporal E/MEG source imaging 
%     optimized by microstates

%% Initialize 
[nSensor,nSource] = size(L); % number of channels & source dipoles
nSamp = size(B,2); % snapshots
if norm(B,'fro')<sqrt(numel(B))
    B = B./norm(B,'fro')*sqrt(numel(B));
end
C_noise = eye(nSensor); % noise covariance
F = L/M;
numW = size(Phi,1); % number of TBFs
W = zeros(nSource,numW);
W_t = W;
C_w = cell(numW,1);
C_phi = eye(numW)*1e-6;

lambda = ones(nSource,1)*trace(F*F.')*trace(Phi*Phi.')/(trace(B*B'));
gamma = ones(nSource,1)*trace(F*F.')*trace(Phi*Phi.')/(trace(B*B'));
alpha = ones(numW,1)*trace(F*F.')*trace(Phi*Phi.')/(trace(B*B'));

diagCk = zeros(numW,1);
E = zeros(nSource,nSamp);
dipIdx = 1:nSource;
maxIter =20;
epsilon = 1e-3;
prune = [1e-6,1e-1];
cost = 0;
if nargin>4
   for inar = 1:2:length(varargin)
       Param = lower(varargin{inar});
       Value = varargin{inar+1};
       switch Param
           case 'maxiter'
               maxIter = Value;
           case 'epsilon'
               epsilon = Value;
       
               
       end
   end   
end

for iter = 1:maxIter
    %% check W and TBF
    tbfIdx = abs(1./alpha)>(max(abs(1./alpha))*prune(2));
    if ~prod(tbfIdx)
        Phi = Phi(tbfIdx,:);
        W = W(:,tbfIdx);
        W_t = W_t(:,tbfIdx);
        C_w = C_w(tbfIdx);
        C_phi = C_phi(tbfIdx,tbfIdx);
        alpha = alpha(tbfIdx);
        diagCk = diagCk(tbfIdx);
        numW = sum(tbfIdx);
    end
     
    wIdx = abs(1./gamma(dipIdx))>(max(abs(1./gamma(dipIdx)))*prune(1));
    if ~prod(wIdx)
        dipIdx = dipIdx(wIdx);
    end
    
    %% spatial coefficient
    FAF = F.*repmat(1./gamma',nSensor,1)* F.';
    for k = 1:numW
        tIdx = setdiff(1:numW,k);
        x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
        C_w{k} = FAF+C_noise/x;  
        r_k = (B-L*E)*Phi(k,:).'-L*sum(W(:,tIdx).*repmat(Phi(k,:)*Phi(tIdx,:).'+nSamp*C_phi(k,tIdx),nSource,1),2);
        W_t(:,k) = repmat(1./gamma,1,nSensor).* F.'/C_w{k}*r_k/x;
        diagCk(k) = trace(FAF-FAF/C_w{k}*FAF);
    end
    W = M\W_t;
    %% Phi
    LW = L*W;    
    C_p = (LW.'*LW+diag(alpha+diagCk));
    C_phi = pinv(C_p);
    Phi = C_phi*LW.'*(B-L*E);

    %% residual
    C_e = F.*repmat(1./lambda',nSensor,1)*F.'+C_noise;
    E_t = repmat(1./lambda,1,nSensor).*F.'/C_e*(B-L*W*Phi);
    E = M\E_t;
    
    %% gamma
    mu = zeros(nSource,1);
    for k = 1:numW
        mu = mu+sum((F.'/C_w{k}).*F.',2);
    end
    gamma(dipIdx) = sqrt(mu(dipIdx)./sum(W_t(dipIdx,:).^2,2));
    gamma(gamma>1e16) = 1e16;
    %% alpha
    alpha = 1./(diag(C_phi)+diag(Phi*Phi')/nSamp);
    %% lambda
    beta = sum(F.'/C_e.*F.',2);
    lambda(dipIdx) = sqrt(nSamp*beta(dipIdx)./sum(E_t(dipIdx,:).^2,2));
    lambda(lambda>min(lambda)*1e12) = min(lambda)*1e12;
    %% free energy
    tempk = 0;
    temp = 0;
    FAF = F.*repmat(1./gamma',nSensor,1)* F.';
    for k = 1:numW
        x = (Phi(k,:)*Phi(k,:).')+nSamp*C_phi(k,k);
        temp = temp+W_t(:,k).'.*gamma'*W_t(:,k);
        
        tempk = tempk+log(det((FAF+C_noise/x)))+log(x)*nSensor;
    end  
    cost_old = cost;
    
    cost = -norm(B-L*W*Phi-L*E,'fro')^2+nSamp*trace(LW.'*LW*C_phi)+nSamp*(sum(log(alpha))+log(det(C_phi)))-nSamp*(log(det(C_e)))...
        -temp-trace(E_t.'.*repmat(lambda',nSamp,1)*E_t)-tempk;
    if abs((cost-cost_old)/cost)<epsilon
        break;
    end    
   %% source
    S =E+W*Phi;

end
    par.W = W;
    par.Phi = Phi;
    par.E = E;

end



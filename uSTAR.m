function Sources = uSTAR(L,B,M,varargin)

MSmap = []; %microstate maps
maxIter = 20; % maximum iteration for source imaging update
epsilon = 1e-3; % relative change stop condition
MSopt.minTime = 5; % minimum time length for each microstate
MSopt.polarity = 0; % polarity, usually 0 for resting-state and 1 for ERP
if nargin>4
   for inar = 1:2:length(varargin)
       Param = lower(varargin{inar});
       Value = varargin{inar+1};
       switch Param
           case 'msmap'
               MSmap = Value;
           case 'msopt'
               MSopt = Value;
           case 'maxiter'
               maxIter = Value;
           case 'epsilon'
               epsilon = Value;         
               
       end
   end   
end


if isempty(MSmap)   
    [~,~,Phi] = svd(B);
    Phi = Phi(:,1:2)';
    Sources = STAR(L,B,Phi,M,'maxiter',maxIter,'epsilon',epsilon);
else
    nsamp = size(B,2);
    nsource = size(L,2);
    mslabel = MicroSmooth(B,MSmap,'reject segments',MSopt);
    tmark = [0 diff(mslabel)]~=0;
    clist = find(tmark==1);
    clist = [0 clist nsamp];
    Sources = zeros(nsource,nsamp);
    for numSeg = 1:(sum(tmark)+1)
        MSseries = B(:,(clist(numSeg)+1):clist(numSeg+1));
        [~,~,Dic] = svd(MSseries);
        Dic = Dic(:,1:2)';
        [S_sted] = STAR(L,MSseries,Dic,nbV,'maxIter',maxIter,'epsilon',epsilon);
        Sources(:,(clist(numSeg)+1):clist(numSeg+1)) = S_sted;
    end
end




end
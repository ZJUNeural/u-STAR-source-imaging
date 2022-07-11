%%
load('LeadField.mat');
load('Cortex.mat');
load('Chanlocs.mat');

Gain=bst_gain_orient(LeadField.Gain,LeadField.GridOrient,LeadField.GridAtlas);
[nSensor,nSource] = size(Gain);
remat = eye(nSensor) - ones(nSensor)/nSensor;
Gain = remat*Gain;

fprintf('Load default settings done\n');
%% Generate Simulated EEG Data
SNR = 10; % signal to noise ratio
K = 2; % number of patches
fsrate = 100; % sampling rate
epochT = [-1 1]; % generate ERP
Time = epochT(1):1/fsrate:epochT(2);
StimTime = 0;
S_area = 15; % area of extended source
% ---------------------- Seed Voxel --------------%
ind = randperm(size(Gain,2));
seedvox = ind(1:K); % randomly select the seedvox
% % ------------------ Active Sources -------------%
AreaDef = S_area*1e-4*ones(numel(seedvox),1); % Area of each extended sources (unit: mm^2)
ActiveVoxSeed = num2cell(seedvox);
ActiveVox = [];
[~, VertArea] = tess_area1(Cortex);

for k = 1:numel(seedvox)
    ActiveVoxSeed{k} = PatchGenerate(seedvox(k),Cortex.VertConn,VertArea,AreaDef(k));
    ActiveVox = union(ActiveVoxSeed{k},ActiveVox);
end
Area = sum(VertArea(ActiveVox));

% %------------------ Simulation data ---------------------%
s_real = zeros(nSource,numel(Time));
for k = 1:K
    Basis = erp_generate_signal(1000*(StimTime-epochT(1)),4*randi([90,100]),1*randi([70,100]),fsrate,...
        1000*abs(diff(epochT)));
    Basis = [2*Basis(1)-Basis(2);Basis(:)]';
    s_real(ActiveVoxSeed{k},:) = repmat(Basis,numel(ActiveVoxSeed{k}),1);
end
B = awgn(Gain*s_real,SNR,'measured');
fprintf('Actual SNR is %g\n',20*log10(norm(Gain*s_real,'fro')/norm(B-Gain*s_real,'fro')));
B = B./max(B(:)); % normalize
B = remat*B;

[~,temp] = max(abs(B(:)));
[~,maxT] = ind2sub(size(B),temp);

offset = max(std(B,[],2));
plotdata=bsxfun(@plus,B',(0:nSensor-1)*offset);
plotdata = plotdata';
%% Visualization


% smooth
iVertices=1:size(Cortex.Vertices,1);
% Smoothing factor
SurfSmoothIterations = ceil(300 * 0.3 * length(iVertices) / 100000);
% Calculate smoothed vertices locations
loc_sm=Cortex.Vertices;
loc_sm(iVertices,:) = tess_smooth(Cortex.Vertices(iVertices,:), 0.3, SurfSmoothIterations, Cortex.VertConn(iVertices,iVertices), 1);
% Apply smoothed locations
Cortex.Vertices=loc_sm;

figure;
plot(Time,plotdata,'k','Linewidth',0.8)
set(gca,'YTick',(0:nSensor-1)*offset,'YTickLabel','');
xlabel('Time(sec)');
axis([Time(1) Time(end) -offset-0.1 (nSensor+1)*offset]);
%%
figure('color','k')
hp=patch(gca,'vertices',Cortex.Vertices,'faces',Cortex.Faces,...
    'FaceColor',[0.75 0.75 0.75],'edgecolor','none','facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud');
material dull
% camlight('headlight','infinite');
set(gca,'color','k','Xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0],'CameraPosition', [0 0 0],'CameraViewAngle',6);
view([ -90 90 ])
axis equal
light('Position',[-100,0,-100]*mean(sum(Cortex.Vertices.^2,2)));light('Position',[100,0,100]*mean(sum(Cortex.Vertices.^2,2)));
colorMap = jet;
cmax = max(sum(s_real,2));cmin = 0; 
tlen = size(colorMap,1);
colorMap(1:floor(tlen*0.15),:) = repmat([.75 .75 .75],floor(tlen*0.15),1);
cdata=colorMap(floor(min((sum(s_real,2)-cmin)/(cmax-cmin),1)*(length(colorMap)-1))+1,:);

set(hp,'AlphaDataMapping','none', 'FaceVertexCData',cdata, 'facecolor', 'interp');
% colormap jet;
%% Source Reconstruction
[~,~,Dic] = svd(B);
Dic = Dic(:,1)';
nbV = diag(sum(Cortex.VertConn,2))-0.9*Cortex.VertConn;

[S_sted] = ESSTED(Gain,B,Dic,nbV,'maxIter',20,'wtype','specifi');

%%
s_rec = abs(S_sted);
figure('color','k')
hp=patch(gca,'vertices',Cortex.Vertices,'faces',Cortex.Faces,...
    'FaceColor',[0.75 0.75 0.75],'edgecolor','none','facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud');
material dull
% camlight('headlight','infinite');
set(gca,'color','k','Xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0],'CameraPosition', [0 0 0],'CameraViewAngle',6);
view([ -90 90 ])
axis equal
light('Position',[-100,0,-100]*mean(sum(Cortex.Vertices.^2,2)));light('Position',[100,0,100]*mean(sum(Cortex.Vertices.^2,2)));
colorMap = hot;
cmax = max(sum(s_rec,2));cmin = 0; 
tlen = size(colorMap,1);
colorMap(1:floor(tlen*0.1),:) = repmat([.75 .75 .75],floor(tlen*0.1),1);
cdata=colorMap(floor(min((sum(s_rec,2)-cmin)/(cmax-cmin),1)*(length(colorMap)-1))+1,:);
colormap hot
set(hp,'AlphaDataMapping','none', 'FaceVertexCData',cdata, 'facecolor', 'interp');
%%
figure
PlotSource(s_real(:,maxT),Cortex,'smooth',30,'axes',gca,'clim',[-1 1],'thresh',0.05);
figure
PlotSource(s_rec(:,maxT),Cortex,'smooth',30,'axes',gca,'clim',[-1 1]*max(s_rec(:,maxT)),'thresh',0.05);
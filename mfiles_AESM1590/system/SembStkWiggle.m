function SembStkWiggle(Data,hdr,geo,vmin,vmax,vstep,CMP,dCMP,SMute,VPicks,VA,perc,ExpSmooth)
% SEMBwiggle(Data,TrcHdr,Velocity,CMP,dCMP,VPicks,VA,ExpSmooth) creates an 
% interactive weighted semblance plot of the seismic data 'Data'. The first
% figure contains the semblance analysis while the second shows the picked
% velocity model and associated stack.
%
% INPUTS:
%  Data:         nt x ntr floating point array containing seismic traces.
%                nt is the number of time samples while ntr is the number
%                of traces.
%  hdr:          Numeric  5 x ntr matirx containg trace headers.
%  geo:          Numeric 6 element vector containg geometry information.
%  Velocity:     Vector of velocities for the semblance analysis.
%  CMP:          First CMP to be picked.
%  dCMP:         CMP steping rate when using the left or right keyboard
%                keys. Right increases current CMP by dCMP, while left
%                decreases current CMP by dCMP. When going beyond possible
%                CMPs the current CMP is is clipped to the first or last
%                possible CMP. Up or down keys increase or decrease
%                respectively the current CMP by 10*dCMP.
% OPTIONAL:
%  VPicks=[];    (nt+1) x n matrix containing picked velocities with which
%                to initialize the velocity model. The first row contains
%                the CMP location while the subsequent rows contain the
%                velocity picks for each time sample.
%  VA=true:      Variable area wiggle plot boolean. If true then positive
%                wiggles are filled black.
%                NOTE: For each trace the first and last value are zeroed.
%  ExpSmooth=50: defines an exponentialy decaying smoothing window around
%                zero-lag of exp(-ExpSmooth.*abs(t)).
%
% The first figure is for semblance analysis. On the left is a wiggle plot
% of the data and on the right is a semblance panel of the data for the
% velocities in 'Velocity'. The black line (if present) on the semblance
% plot indicates the expected velocity trend based on velocity picks in
% other CMPS.
%
% Left clicking in the semblance plot will show you the associated Normal
% MoveOut (NMO) as a red line on the wiggle plot. Use this and the
% semblance plot to guide your picking of velocities. Clicking and dragging
% allows you to adjust the location of your velocity pick while
% interactively changing the NMO curve. Right clicking in the semblance
% plot allows you to delete the last velocity pick.
%
% Right clicking in the wiggle plot will cause the wiggle plot to display
% the NMO corrected CMP for your velocity picks for this CMP. You can use
% this to check your picks. Middle clicking in the wiggle plot causes the
% original data to be displayed.
%
% The second figure shows the interpolated velocity model based on your
% picks so far. Below it you will find the associated stack of the data.
% Right clicking on the stack causes the data to be stacked using the
% picked velocity model. You can use this to check your picks.
%
% When you are done picking velocities inside the semblance figure, you can
% navigate to a new CMP by either using the left, up, down, or right
% keyboard keys or by clicking on the desired CMP in either the velocity
% model or the stack.
%
% NOTE: Upon closing this program it will ask you if you would like to save
% your picked velocity model and stack. If you click 'yes' then it will
% save them in either V, Stk, & VPicks or if these are in use in Vn, Stkn
% and  VPicksn, where n is an increasing number so as not to overwrite
% previously exported data. If no velocity model was picked then nothing
% will be returned to the workspace.
%
% Using the Interactive Colorbar
% - Dragging the black bars at the top and bottom of the colorbar
%   changes the clipping level.
% - When right clicking on the colorbar one has the following options:
%   - Gray:              Changes the colormap to grayscale
%   - Inverted Gray:     Changes the colormap to inverted grayscale
%   - Seismic:           Changes the colormap to seismic blue-white-red
%                        This is the default colormap
%   - Colormap Editor:   Opens the matlab colormap editor
%   - Set Limits:        Sets the current colorbar clipping level to the
%                        colorbar limits. This allows for more fine tuned
%                        clipping around the center of the colorbar.
%   - Reset Limits:      Resets colorbar limits to their original values.
%   - Symetric Clipping: Toggles if clipping levels should be symmetric
%                        around the center of the colorbar.
%   - Connect Colorbar:  Connect interactive colorbar with other
%                        interactive colorbars so that they all clip the
%                        same.
%                        NOTE: THIS FUNCTIONALITY SHOULD NOT BE USED FOR
%                        THIS PLOT.
%
% COPYRIGHT:
%    Max Enno Holicki
%    21.09.2016
%    The Netherlands

%% Check Inputs
if(~isnumeric(Data));
	fprintf(2,'SembStkWiggle: Data is not numeric (double, float, integer, etc.).\n');
	return;
elseif(~isreal(Data));
	fprintf(2,'SembStkWiggle: Complex Data is not supported.\n');
	return;
elseif(~ismatrix(Data));
	fprintf(2,'SembStkWiggle: Data is not a matrix.\n');
	return;
end;
if ~isnumeric(hdr)||~ismatrix(hdr)||size(hdr,1)~=5;
    fprintf(2,'SembStkWiggle: hdr is not a numeric matirx of size 5xn.\n');
    fprintf(2,'SembStkWiggle: Is hdr a header matrix?\n');
    return;
end;
if ~isnumeric(geo)||~isvector(geo)||numel(geo)~=7;
    fprintf(2,'SembStkWiggle: geo is not a numeric vector with six elements.\n');
    fprintf(2,'SembStkWiggle: Is geo a geometry matrix?\n');
    return;
end;
if ~isscalar(vmin)||~isnumeric(vmin);
    fprintf(2,'SembStkWiggle: vmin is not a numeric scalar!\n');
    return;
end;
if ~isscalar(vmax)||~isnumeric(vmax);
    fprintf(2,'SembStkWiggle: vax is not a numeric scalar!\n');
    return;
end;
if ~isscalar(vstep)||~isnumeric(vstep);
    fprintf(2,'SembStkWiggle: vstep is not a numeric scalar!\n');
    return;
end;
if(~isnumeric(CMP));
	fprintf(2,'SembStkWiggle: CMP is not numeric (double, float, integer, etc.).\n');
	return;
elseif(~isreal(CMP));
	fprintf(2,'SembStkWiggle: Complex CMP is not supported.\n');
	return;
elseif(~isscalar(CMP));
	fprintf(2,'SembStkWiggle: CMP is not a scalar.\n');
	return;
end;
if(~isnumeric(dCMP));
	fprintf(2,'SembStkWiggle: dCMP is not numeric (double, float, integer, etc.).\n');
	return;
elseif(~isreal(dCMP));
	fprintf(2,'SembStkWiggle: Complex dCMP is not supported.\n');
	return;
elseif(~isscalar(dCMP));
	fprintf(2,'SembStkWiggle: dCMP is not a scalar.\n');
	return;
end;
if(~exist('SMute','var'));SMute=0;
elseif(isempty(SMute));SMute=0;
elseif ~isscalar(SMute)||~isnumeric(SMute);
    fprintf(2,'SembStkWiggle: SMute is not a numeric scalar!\n');
    return;
end;
if(~exist('VA','var'));VA=true;
elseif(isempty(VA));VA=true;
elseif(~islogical(VA)||~isscalar(VA));
	fprintf(2,'SembStkWiggle: VA is not a logical scalar.\n');
	return;
end;
if ~exist('perc','var');perc=0;
elseif(isempty(perc));perc=0;
elseif ~isscalar(perc)||~isnumeric(perc);
    fprintf(2,'SembStkWiggle: perc is not a numeric scalar!\n');
    return;
elseif perc>100||perc<0;
    fprintf(2,'SembStkWiggle: perc should be between 0 and 100!\n');
    return;
end;
if(~exist('ExpSmooth','var'));ExpSmooth=50;
elseif(isempty(ExpSmooth));ExpSmooth=0;
elseif(~isnumeric(ExpSmooth));
	fprintf(2,'SembStkWiggle: ExpSmooth is not numeric (double, float, integer, etc.).\n');
	return;
elseif(~isscalar(ExpSmooth));
	fprintf(2,'SembStkWiggle: ExpSmooth is not a scalar.\n');
	return;
end;
%% Convert hdr to TrcHdr structure
TrcHdr=struct('cdp',num2cell(hdr(5,:)),'offset',num2cell(hdr(2,:)));
%% Generate Time Vector
dt=geo(1)/1000;
t=(0:geo(2)-1)'*dt;
Velocity=permute(vmin:vstep:vmax,[1,3,2]);
%% Sort & Bin CMPs
[Data,TrcHdr]=SeisSort(Data,TrcHdr,{'cdp','offset'},1);
[Data,TrcHdr]=SeisBin(Data,TrcHdr,'cdp');
%% Get CMP Locations
cmp=double(cellfun(@(a)a(1).cdp,TrcHdr));
ns=num2str(numel(num2str(max(cmp))));
CMPproc=false(1,numel(cmp));
%% Set-Up Velocity Model + Stack
V=nan(numel(t),numel(Data)); %Velocity Model
if(~exist('VPicks','var'));
elseif(isempty(VPicks));
elseif(~isnumeric(VPicks));
	fprintf(2,'semb: VPicks is not numeric (double, float, integer, etc.).\n');
	return;
elseif(~ismatrix(VPicks));
	fprintf(2,'semb: VPicks is not a scalar.\n');
	return;
elseif(size(VPicks,1)-1~=numel(t));
	fprintf(2,'semb: Mismatch between size(VPicks,1)-1 & size(Data,1).\n');
	return;
elseif(size(VPicks,2)>size(Data,1));
	fprintf(2,'semb: size(VPicks,2) > size(Data,1).\n');
	return;
else
	if size(VPicks,2)==1;
		% Initialize Picks
		CMPproc(cmp==VPicks(1,1))=true;
		for i=1:size(V,2);V(:,i)=VPicks(2:end,1);end;
	else
		% Initialize Picks
		CMPLoc=zeros(size(VPicks,2),1);
		for i=1:size(VPicks,2);
			CMPLoc(i)=find(cmp==VPicks(1,i));
			CMPproc(CMPLoc(i))=true;
			V(:,CMPLoc(i))=VPicks(2:end,i);
		end;
		% Interpolate
		V(:,1:CMPLoc(1)-1)=repmat(V(:,CMPLoc(1)),[1,CMPLoc(1)-1]);
		for i=2:size(VPicks,2)-1;
			V(:,CMPLoc(i-1)+1:CMPLoc(i)-1)=interp2(CMPLoc([i-1,i]),t,V(:,CMPLoc([i-1,i])),CMPLoc(i-1)+1:CMPLoc(i)-1,t);
		end;
		V(:,CMPLoc(end-1)+1:CMPLoc(end)-1)=interp2(CMPLoc([end-1,end]),t,V(:,CMPLoc([end-1,end])),CMPLoc(end-1)+1:CMPLoc(end)-1,t);
		V(:,CMPLoc(end)+1:end)=repmat(V(:,CMPLoc(end)),[1,size(V,2)-CMPLoc(end)]);
	end;
end;
Stk=zeros(numel(t),numel(Data)); %NMO - Stack
%% Plot 1st CMP
CMP=find(cmp==CMP);
if(isempty(CMP));
	fprintf(2,'SembStkWiggle: CMP does not correspond to a possible cmp location.\n');
	return;
end;
% Extract Offsets
x=double([TrcHdr{CMP}.offset]);
% Compute 1st Semblance
S=semb(Data{CMP},x,t,Velocity,ExpSmooth);
%% Semblance Figure
ddcmp=median(diff([TrcHdr{CMP}.offset]));if isnan(ddcmp);ddcmp=1;end;
if(perc);PERC=prctile(Data{CMP}(:),perc);else PERC=max(abs(Data{CMP}(:)));end;
tmp=Data{CMP}./PERC;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
tmp=bsxfun(@times,tmp,ddcmp);
F(1)=figure('CloseRequestFcn',@closereq,'renderer','painters');
AH(1)=subplot(1,2,1);
%% Plot Variable Area
if(VA); % Check for variable area flag
	tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
	[Face,Vert]=VariableArea(tmp,t,x);
	PH=patch('Faces',Face,'Vertices',Vert,'FaceColor','Black','EdgeColor','none');
    set(PH,'BusyAction','cancel','HitTest','off','Interruptible','off');
	axis ij;
	hold on;
end;
%% Plot Wiggle
tmp=bsxfun(@plus,[TrcHdr{CMP}.offset],tmp);
IH(1)=plot(reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'k');
set(IH(1),'BusyAction','cancel','HitTest','off','Interruptible','off');
if(VA);
	hold off;
else
	axis ij;
end;
axis([min([TrcHdr{CMP}.offset])-ddcmp/2,max([TrcHdr{CMP}.offset])+ddcmp/2,t(1),t(end)]);
xlabel('Offset []');ylabel('Two-Way-Traveltime []');
TH(1)=title(sprintf('CMP Gather: %d',cmp(CMP)));
%% Plot Semblance
AH(2)=subplot(1,2,2);IH(2)=imagesc(permute(Velocity,[1,3,2]),t,S,[0,1]);
hCB(2)=InteractiveColorbar;colormap(jet);
xlabel('NMO Velocity');
TH(2)=title('Semblance Plot');
linkaxes(AH,'y');
% Set UserData, Colormap & Sym
set(AH(1),'UserData',false);set(AH(2),'UserData',true);
% Set Pointer & Window Functions
if(isa(F(1),'matlab.ui.Figure'));
    set(F(1),'Pointer','CrossHair');
else
    set(F(1),'Pointer','Full CrossHair');
end;
set(F(1),'WindowButtonDownFcn',@clickFcn);
set(F(1),'WindowButtonUpFcn',@unclickFcn);
set(F(1),'KeyPressFcn',@keyDownListener);
set(AH(1),'buttondownfcn',@NMOcorrGather);
% Set up cursor lines
hNMO=line(x,NaN(1,numel(x)),'Color','red','LineWidth',1.5,'Parent',AH(1)); %NMO Line Handle
set(hNMO,'BusyAction','cancel','HitTest','off','Interruptible','off');
hVM =line(NaN(3,1),NaN(3,1),'Color','red','LineWidth',1.5,'Parent',AH(2)); %Velocity Line Handle
set(hVM,'BusyAction','cancel','HitTest','off','Interruptible','off');
indord=0; %Vector last clicked velocity profile indeces so that they can be remove
% Set up auxillary lines
hVMo=line(V(:,CMP),t,'Color','white','LineWidth',1.5,'Parent',AH(2)); %OriginalVelocity Line Handle
set(hVMo,'BusyAction','cancel','HitTest','off','Interruptible','off');
%% Plot Velocity Model & NMO Stack
F(2)=figure('CloseRequestFcn',@closereq);
colormap(gray);
AH(3)=subplot(2,1,1);IH(3)=imagesc(cmp,t,V,[min(Velocity),max(Velocity)]);
hCB(3)=InteractiveColorbar;
xlabel('CMP Location []');ylabel('Two-Way-Traveltime []');title('RMS Velocity Model');
AH(4)=subplot(2,1,2);IH(4)=imagesc(cmp,t,Stk,[-1,1]);
hCB(4)=InteractiveColorbar;colormap(gray);hCB(4).Sym=true;
xlabel('CMP Location []');ylabel('Two-Way-Traveltime []');title('Stack');
dCMPh=uicontrol('Style','Edit','String',dCMP,'Position',[0,0,23,23],'Backgroundcolor',[1,1,1],'Foregroundcolor',[0,0,0]);
linkaxes(AH(3:4),'xy');
% Set up auxillary lines
hCMPL(1)=line([cmp(CMP),cmp(CMP)],[t(1)-0.5.*dt,t(end)+0.5.*dt],'Color','red','LineWidth',1.5,'Parent',AH(3)); %Velocity Model CMP-Line
hCMPL(2)=line([cmp(CMP),cmp(CMP)],[t(1)-0.5.*dt,t(end)+0.5.*dt],'Color','red','LineWidth',1.5,'Parent',AH(4)); %Stk Model CMP-Line
set(hCMPL,'BusyAction','cancel','HitTest','off','Interruptible','off');
% Set Window Functions
if(isa(F(1),'matlab.ui.Figure'));
    set(F(1),'Pointer','CrossHair');
else
    set(F(1),'Pointer','Full CrossHair');
end;
set(F(2),'WindowButtonUpFcn',@selectCMP);
drawnow;
%% Get Java Window Handles
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jSemb=getJFrame(F(1));
jStk=getJFrame(F(2));
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%% Window Functions
function UpdateVelModelStk
	tmp=get(hVM,'YData');
	if any(isnan(tmp));return;end
	V(:,CMP)=interp1(tmp,get(hVM,'XData'),t);
	% Process Left
	ind=find(CMPproc(1:CMP-1),1,'last');
	if isempty(ind);
		V(:,1:CMP-1)=repmat(V(:,CMP),1,numel(1:CMP-1));
	else
		V(:,ind+1:CMP-1)=interp2([ind,CMP],t,V(:,[ind,CMP]),ind+1:CMP-1,t);
	end;
	% Process Right
	ind=CMP+find(CMPproc(CMP+1:end),1,'first');
	if isempty(ind);
		V(:,CMP+1:end)=repmat(V(:,CMP),1,numel(CMP+1:numel(cmp)));
	else
		V(:,CMP+1:ind-1)=interp2([CMP,ind],t,V(:,[CMP,ind]),CMP+1:ind-1,t);
	end;
	% Update Velocity Model Image
	set(IH(3),'CData',V);
	% Update Process CMP
	CMPproc(CMP)=true;
end
function UpdatePlot
CMP
	ddcmp=median(diff([TrcHdr{CMP}.offset]));if isnan(ddcmp);ddcmp=1;end;
    if(perc);PERC=prctile(Data{CMP}(:),perc);else PERC=max(abs(Data{CMP}(:)));end;
    tmp=Data{CMP}./PERC;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
	tmp=bsxfun(@times,tmp,ddcmp);
	if(VA);
		tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
		[Face,Vert]=VariableArea(tmp,t,[TrcHdr{CMP}.offset]);
		set(PH,'Faces',Face,'Vertices',Vert);
	end;
	tmp=bsxfun(@plus,[TrcHdr{CMP}.offset],tmp);
	set(IH(1),'XData',reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'YData',reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]));
    axis(AH(1),[min([TrcHdr{CMP}.offset])-ddcmp/2,max([TrcHdr{CMP}.offset])+ddcmp/2,t(1),t(end)]);
	set(IH(2),'CData',semb(Data{CMP},[TrcHdr{CMP}.offset],t,Velocity,ExpSmooth));
	set(AH(1),'XLim',[min([TrcHdr{CMP}.offset])-0.5,max([TrcHdr{CMP}.offset])+0.5]);
	set(TH(1),'String',sprintf(['%s %0',ns,'d'],'CMP Gather:',cmp(CMP)));
	set(hNMO,'XData',[TrcHdr{CMP}.offset],'YData',nan(numel([TrcHdr{CMP}.offset]),1));
	set(hVM,'XData',NaN(3,1),'YData',NaN(3,1));
	set(hVMo,'XData',V(:,CMP));
	set(hCMPL(1),'XData',[cmp(CMP),cmp(CMP)]);
	set(hCMPL(2),'XData',[cmp(CMP),cmp(CMP)]);
end
function keyDownListener(~,event)
	SS=str2double(get(dCMPh,'String'));
	if isnan(SS);return;end;
	switch event.Key
		case 'leftarrow'; %Show previous record
			if CMP>1;
                UpdateVelModelStk;
				if CMP>SS;
					CMP=CMP-SS;
				else
					CMP=1;
				end;
				UpdatePlot;
			end;
		case 'downarrow'; %Take a larger step back
			if CMP>1;
                UpdateVelModelStk;
				SS=10.*SS;
				if CMP>SS;
					CMP=CMP-SS;
				else
					CMP=1;
				end;
				UpdatePlot;
			end;
		case 'rightarrow'; %Show next record
			if CMP<numel(cmp);
                UpdateVelModelStk;
				if CMP+SS<=numel(cmp);
					CMP=CMP+SS;
				else
					CMP=numel(cmp);
				end;
				UpdatePlot;
			end;
		case 'uparrow'; %Take a large step forward
			if CMP<numel(cmp);
                UpdateVelModelStk;
				SS=10.*SS;
				if CMP+SS<=numel(cmp);
					CMP=CMP+SS;
				else
					CMP=numel(cmp);
				end;
				UpdatePlot;
			end;
	end;
end
function clickFcn(varargin)
	if strcmpi(get(gcbf,'SelectionType'),'normal');
		if IH(2)~=gco;
			set(hNMO,'YData',NaN(1,numel([TrcHdr{CMP}.offset])));
		else
			set(gcf,'WindowButtonMotionFcn',@dragFcn);
			dragFcn;
		end;
	end;
end
function dragFcn(varargin)
	% Get mouse location
	pt=get(gca,'CurrentPoint');
	% Update cursor line position
	set(hNMO,'YData',sqrt(pt(3).^2+([TrcHdr{CMP}.offset]/pt(1)).^2));
end
function NMOcorrGather(varargin)
	if(strcmpi(get(gcbf,'SelectionType'),'alt'));
		tmp=get(hVM,'XData');
		if(any(isnan(tmp)));
			if(~all(isnan(tmp)));
				fprintf(2,'SEMB: Warning: NaNs detected in velocity picks.\m')
			end;
			return;
		end;
		Vel=interp1(get(hVM,'YData'),get(hVM,'XData'),t);
		T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,[TrcHdr{CMP}.offset],Vel).^2)); %NMO Times
        tmp=Data{CMP};
		for j=1:numel(TrcHdr{CMP});
			tmp(:,j)=interp1(t,Data{CMP}(:,j),T(:,j),'linear',0); %NMO correct
		end;
        tmp=StretchMute(tmp,[TrcHdr{CMP}.offset],t,Vel,SMute);
        tmp=tmp./PERC;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
		tmp=bsxfun(@times,tmp,ddcmp);
		if(VA);
			tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
			[Face,Vert]=VariableArea(tmp,t,[TrcHdr{CMP}.offset]);
			set(PH,'Faces',Face,'Vertices',Vert);
		end;
		tmp=bsxfun(@plus,[TrcHdr{CMP}.offset],tmp);
		set(IH(1),'XData',reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'YData',reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]));
        set(hNMO,'XData',[TrcHdr{CMP}.offset],'YData',nan(numel([TrcHdr{CMP}.offset]),1));
	elseif(strcmpi(get(gcbf,'SelectionType'),'extend'));
		tmp=Data{CMP}./PERC;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
		tmp=bsxfun(@times,tmp,ddcmp);
		if(VA);
			tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
			[Face,Vert]=VariableArea(tmp,t,[TrcHdr{CMP}.offset]);
			set(PH,'Faces',Face,'Vertices',Vert);
		end;
		tmp=bsxfun(@plus,[TrcHdr{CMP}.offset],tmp);
		set(IH(1),'XData',reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'YData',reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]));
	end;
end
function selectCMP(varargin)
	if strcmpi(get(gcbf,'SelectionType'),'alt');
		if any(gco==IH(4));
			CMPStk;
		end;
	else
		if any(gco==IH([3,4]));
			pt=get(gca,'CurrentPoint');
			UpdateVelModelStk;
			CMP=round(interp1(cmp,1:numel(cmp),pt(1)));
			if CMP<1;CMP=1;elseif CMP>numel(cmp);CMP=numel(cmp);end;
			UpdatePlot;
		end;
	end;
end
function unclickFcn(varargin)
	if IH(2)==gco;
		switch get(gcbf,'SelectionType')
			case 'normal'
				set(gcf,'WindowButtonMotionFcn','');
				% Get mouse location
				pt=get(gca,'CurrentPoint');
				% Update velocity model
				Y=get(hVM,'YData');
				if all(isnan(Y));
					set(hVM,'XData',[pt(1),pt(1),pt(1)],'YData',[t(1),pt(3),t(end)]);
					indord=2;
				else
					Ylind=find(Y>pt(3),1,'first');
					X=get(hVM,'XData');
					if Ylind==2;
						set(hVM,'XData',[pt(1),pt(1),X(Ylind:end)],'YData',[t(1),pt(3),Y(Ylind:end)]);
					elseif Ylind==numel(Y);
						set(hVM,'XData',[X(1:end-1),pt(1),pt(1)],'YData',[Y(1:end-1),pt(3),t(end)]);
					else
						set(hVM,'XData',[X(1:Ylind-1),pt(1),X(Ylind:end)],'YData',[Y(1:Ylind-1),pt(3),Y(Ylind:end)]);
					end;
					indord=[Ylind,indord];
				end;
			case 'alt';
				Y=get(hVM,'YData');
				switch indord(1);
					case 0
					case 2
						if numel(indord)==1;
							set(hVM,'XData',NaN(3,1),'YData',NaN(3,1));
						else
							Y(2)=[];
							X=get(hVM,'XData');X([1,2])=[];
							set(hVM,'XData',[X(1),X],'YData',Y);
						end;
					case numel(Y)-1
						Y(end-1)=[];
						X=get(hVM,'XData');X(end-1)=[];X(end)=X(end-1);
						set(hVM,'XData',X,'YData',Y);
					otherwise
						Y(indord(1))=[];
						X=get(hVM,'XData');X(indord(1))=[];
						set(hVM,'XData',X,'YData',Y);
				end;
				if numel(indord)==1;
					indord=0;
				else
					indord(1)=[];
				end;
		end;
	end;
end
function CMPStk
	% Disable Figures
	jSemb.setEnabled(false);
	jStk.setEnabled(false);
	% Create Wait Bar
	WB=waitbar(0,'CMP Stacking','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
	setappdata(WB,'canceling',0)
	% CMP Stack
	for i=1:numel(cmp);
        off=[TrcHdr{i}.offset];
		if getappdata(WB,'canceling');
			break;
		end;
		T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,[TrcHdr{i}.offset],V(:,i)).^2)); %NMO Times
		Stk(:,i)=0;
		for j=1:numel(TrcHdr{i});
			Stk(:,i)=Stk(:,i)+StretchMute(interp1(t,Data{i}(:,j),T(:,j),'linear',0),off(j),t,V(:,i),SMute); %NMO correct
		end;
		waitbar(i/numel(cmp),WB);
	end;
	% Enable Figures
	jSemb.setEnabled(true);
	jStk.setEnabled(true);
	% Draw Stk
	set(IH(4),'CData',Stk);
	% Delete Waitbar
	delete(WB);
end
function closereq(~,~)
	% Check if user wishes to export velocity model
	selection=questdlg('Export Velocity Model & Stack?','Export Velocity Model & Stack?','Yes','No','Cancel','No'); 
	switch selection, 
		case 'Yes',
			if all(isnan(V));
				fprintf(1,'semb: No velocity picks present. Not exporting velocity model.\n');
			else
				tmp=get(IH(3),'CData');
				if evalin('base','exist(''V'',''var'');');
					i=0;
					while true;
						if ~evalin('base',sprintf('exist(''V%d'',''var'');',i));
							assignin('base',sprintf('V%d',i),tmp);
							fprintf(1,'semb: Saved velocity model as V%d.\n',i);
							break;
						else
							i=i+1;
						end;
					end;
				else
					assignin('base','V',tmp);
					fprintf(1,'semb: Saved velocity model as V.\n');
				end;
				if evalin('base','exist(''VPicks'',''var'');');
					i=0;
					while true;
						if ~evalin('base',sprintf('exist(''VPicks%d'',''var'');',i));
							assignin('base',sprintf('VPicks%d',i),[cmp(CMPproc)';tmp(:,CMPproc)]);
							fprintf(1,'semb: Saved velocity picks as VPicks%d.\n',i);
							break;
						else
							i=i+1;
						end;
					end;
                else
					assignin('base','VPicks',[cmp(CMPproc)';tmp(:,CMPproc)]);
					fprintf(1,'semb: Saved velocity picks as VPicks.\n');
				end;
			end;
			if ~any(Stk);
				fprintf(1,'semb: No NMO-Stack present. Not exporting NMO-Stack.\n');
			else
				if evalin('base','exist(''Stk'',''var'');');
					i=0;
					while true;
						if ~evalin('base',sprintf('exist(''Stk%d'',''var'');',i));
							assignin('base',sprintf('Stk%d',i),get(IH(4),'CData'));
							fprintf(1,'semb: Saved NMO-Stack model as Stk%d.\n',i);
							break;
						else
							i=i+1;
						end;
					end;
				else
					assignin('base','Stk',get(IH(4),'CData'));
					fprintf(1,'semb: Saved NMO-Stack as Stk.\n');
				end;
			end;
			delete(F(1));
			delete(F(2));
		case 'No'
			delete(F(1));
			delete(F(2));
		case 'Cancel'
	  return 
	end;
end
end
%% Helper Function
function [S]=semb(data,x,t,v,ExpSmooth)
	% Helper function
	%% Compute 3D NMO time matrix
	T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,x,v).^2));
	%% Interpolate data to NMO time and stack for each velocity
	tt=exp(-ExpSmooth.*abs(bsxfun(@minus,t,t')));
	S=zeros([numel(t),numel(v)]); %Preallocate semblance
	q=zeros([numel(t),numel(x)]); %Preallocate tempoary container
	b=zeros([numel(t),1]);
	for vi=1:numel(v);
		for j=1:numel(x);
		   q(:,j)=interp1(t,data(:,j),T(:,j,vi),'linear',0); %NMO correct
		end;
		r=sum(q,2);
		C=bsxfun(@rdivide,bsxfun(@times,t.*numel(x)/sum(x.^2),x.^2),T);
		C(isnan(C))=0; %C sometimes becomes NaN due to 0/0. We zero these.
		Crq=sum(bsxfun(@times,tt,sum(bsxfun(@times,r   ,q              ),2)),1)';
		Crr=sum(bsxfun(@times,tt,        numel(x).*r.^2                    ),1)';
		Cqq=sum(bsxfun(@times,tt,sum(                   q.^2            ,2)),1)';
		Brq=sum(bsxfun(@times,tt,sum(bsxfun(@times,r   ,C(:,:,vi).*q   ),2)),1)';
		Brr=sum(bsxfun(@times,tt,sum(bsxfun(@times,r.^2,C(:,:,vi)      ),2)),1)';
		Bqq=sum(bsxfun(@times,tt,sum(                   C(:,:,vi).*q.^2 ,2)),1)';
		A=Crr.*Bqq+Cqq.*Brr;
		Rrq=Crq./(Crq-Brq);
		Rrr=Crr./(Crr-Brr);
		Rqq=Cqq./(Cqq-Bqq);
		ind=(Rrr<Rrq&Rrq<Rqq)|(Rqq<Rrq&Rrq<Rrr);
		b(ind)=(1-(2.*Crq(ind).*Brr(ind).*Bqq(ind)-Brq(ind).*A(ind))./(2.*Brq(ind).*Crr(ind).*Cqq(ind)-Crq(ind).*A(ind))).^-1;
		b(~ind)=Rrq(~ind);
		Ind=b>1|b<0;
		b(Ind)=0;
		Wrq=(1-b).*Crq+b.*Brq;
		Wrr=(1-b).*Crr+b.*Brr;
		Wqq=(1-b).*Cqq+b.*Bqq;
		S(:,vi)=Wrq.^2./(Wrr.*Wqq);
		s=Brq(Ind).^2./(Brr(Ind).*Bqq(Ind));
		IndF=find(Ind);
		IND=S(Ind,vi)>s;
		S(IndF(IND),vi)=s(IND);
	end;
    S=S/max(S(:));
	S(isnan(S))=0;
end
function jframe=getJFrame(hFig)
	jf=get(handle(hFig),'javaframe');
	try
		jframe=jf.fFigureClient.getWindow;
	catch
		try
			jframe=jf.fHG1Client.getWindow;
		catch
			jframe=jf.fHG2Client.getWindow;
		end;
	end;
end
function [Face,Vert]=VariableArea(Data,t,x)
	%% Create Variable Area Patches
	% Create Vertices & Faces
	% Vertices
	AreaFill=cell(2,size(Data,2));
	for i=1:size(Data,2);
		[AreaFill{2,i},AreaFill{1,i}]=ZeroClip(t,Data(:,i));
	end;
	Vert=zeros(sum(cellfun(@numel,{AreaFill{1,:}})),2);
	% x
	j=1;k=0;
	for i=1:size(Data,2);
		temp=AreaFill{1,i};
		k=k+numel(temp);
		Vert(j:k,1)=temp+x(i);
		j=k+1;
	end;
	% t
	j=1;k=0;
	for i=1:size(Data,2);
		temp=AreaFill{2,i};
		k=k+numel(temp);
		Vert(j:k,2)=temp;
		j=k+1;
	end;
	% Faces
	Face=nan(size(Data,2),max(cellfun(@numel,{AreaFill{1,:}}))+size(Data,2));
	j1=1;k1=0;
	for i=1:size(Data,2);
		k=numel(AreaFill{1,i});
		k1=k1+k;
		Face(i,1:k)=j1:k1;
		Face(i,k+1)=j1;
		j1=k1+1;
	end;
	%% ZeroClip through linear interpolation
	function [outX,outY]=ZeroClip(inX,inY)
		tmp=inY>0;
		TMP=abs(diff(tmp));
		outX=zeros(sum(tmp)+sum(TMP),1);
		outY=zeros(sum(tmp)+sum(TMP),1);
		outX(1)=inX(1);
		% First Value
		if(tmp(1));outY(1)=inY(1);else outY(1)=0;end;
		% Remaining Values
		j=2;
		for i=1:numel(inY)-2;
			if(TMP(i));
				outX(j)=inX(i)-inY(i)*(inX(i+1)-inX(i))/(inY(i+1)-inY(i));
				outY(j)=0;j=j+1;
			end;
			if(tmp(i+1));outX(j)=inX(i+1);outY(j)=inY(i+1);j=j+1;end
		end;
		if(TMP(end));
			outX(j)=inX(end-1)-inY(end-1)*(inX(end)-inX(end-1))/(inY(end)-inY(end-1));
			outY(j)=0;j=j+1;
		end;
		outX(j)=inX(end);
		if(tmp(end));outY(j)=inY(end);else	outY(j)=0;end;
	end
end
function data=StretchMute(data,x,t0,v,SMute)
    % Hard Stretch Mute from:
    % Processing of Seismic Reflection Data Using Matlab
    if(SMute>0)data((bsxfun(@rdivide,x,t0.*v).^2)/2>SMute)=0;end;
end
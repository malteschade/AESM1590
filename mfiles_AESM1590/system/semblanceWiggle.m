function [hVM,F,AH,hNMO,hCB]=semblanceWiggle(data,TrcH,geo,vmin,vmax,vstep,SMute,perc,VA,ExpSmooth)
% semblanceWiggle(CMPgather,H_CMPgather,geo,cmin,cmax,velstep,smute,perc,VA,expSmooth)
%
% This function creates an interactive weighted semblance plot of a gather 
% (often a CMP gather) using an interactive colorbar. The gather should
% be in matrix format (no. of time samples  x  no. of traces)
%
% Input:  CMPgather   - seismic gather
%         H_CMPgather - header data of seismic gather (e.g., see segyread)
%         geo         - geometry of seismic data (e.g., see segyread)
%         cmin        - minimum velocity in semblance analysis (m/s)
%         cmax        - maximum velocity in semblance analysis (m/s)
%         velstep     - velocity step between cmin and cmax (m/s)
%                       (e.g.(cmax-cmin)/100)
% Input (optional):
%         smute       - Stretch mute value (default is zero - no mute)
%         perc        - percentile for determining clip of plot 
%                       (default = 0 is no perc; suggested: perc=99)
%         VA          - Variable Area plotting 
%                       (default = true, false is no black-filling)
%         expSmooth   - factor of filter with an exponentially decaying 
%                       smoothing window of exp(-ExpSmooth.*abs(t)). 
%                       By Default this value is 50.
%
% NOTE: Upon closing this program it will ask you if you would like to save
% you picked velocity model. If you click 'yes' then it will save them in
% either t & v, or if these are in use in tn & vn where n is an increasing
% number so as not to overwrite previous picks. If no velocity model was
% picked then t & v will not be returned to the workspace.
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

% Max Holicki - designed for AESB2140
% October 2014
% m.e.holicki@tudelft.nl

% This function implements weighted semblance from:
% Velocity anlaysis using weighted semblance
% Simon Luo & Dave Hale
% Center for Wave Phenomena
% Colorado School of Mines, Golden, CO 80401, USA
% CWP-652

% This software makes use of:
% InteractiveColorbar.m (Modified)
% hey dude
% http://www.mathworks.nl/matlabcentral/fileexchange/45363-interactivecolorbar (13.10.2014)
%
% seismic.m (Modified, incorporated into InteractiveColorbar.m)
% Community Surface Dynamics Modeling System
% https://csdms.colorado.edu/svn/sedflux/sedflux-mfiles/sedflux/seismic.m (13.10.2014)
%% Check Inputs
if ~ismatrix(data)||isscalar(data)||~isnumeric(data);
    fprintf(2,'semblanceWiggle: Data is not a 2D numeric matrix!\n');
    return;
end;
if ~ismatrix(TrcH)||size(TrcH,1)~=5||size(TrcH,2)~=size(data,2)||~isnumeric(TrcH);
    fprintf(2,'semblanceWiggle: TrcH is not a 5 by number of receivers numeric matrix!\n');
    return;
end;
if ~isvector(geo)||numel(geo)~=7||~isnumeric(geo);
    fprintf(2,'semblanceWiggle: geo is not a 7 element geometry vector!\n');
    return;
end;
if ~isscalar(vmin)||~isnumeric(vmin);
    fprintf(2,'semblanceWiggle: vmin is not a numeric scalar!\n');
    return;
end;
if ~isscalar(vmax)||~isnumeric(vmax);
    fprintf(2,'semblanceWiggle: vax is not a numeric scalar!\n');
    return;
end;
if ~isscalar(vstep)||~isnumeric(vstep);
    fprintf(2,'semblanceWiggle: vstep is not a numeric scalar!\n');
    return;
end;
if ~exist('SMute','var');SMute=0;
elseif(isempty(SMute));SMute=0;
elseif ~isscalar(SMute)||~isnumeric(SMute);
    fprintf(2,'semblanceWiggle: SMute is not a numeric scalar!\n');
    return;
end;
if(~exist('VA','var'));VA=true;
elseif(isempty(VA));VA=true;
elseif(~islogical(VA)||~isscalar(VA));
	fprintf(2,'semblanceWiggle: VA is not a logical scalar.\n');
	return;
end;
if ~exist('perc','var');perc=0;
elseif(isempty(perc));perc=0;
elseif ~isscalar(perc)||~isnumeric(perc);
    fprintf(2,'semblanceWiggle: perc is not a numeric scalar!\n');
    return;
elseif perc>100||perc<0;
    fprintf(2,'SembStkWiggle: perc should be between 0 and 100!\n');
    return;
end;
if ~exist('ExpSmooth','var');ExpSmooth=50;
elseif(isempty(ExpSmooth));ExpSmooth=0;
elseif ~isscalar(ExpSmooth)||~isnumeric(ExpSmooth);
    fprintf(2,'semblanceWiggle: ExpSmooth is not a numeric scalar\n');
    return;
end;
%% Create Vectors
x=TrcH(2,:);
t=(0:geo(2)-1)*geo(1)/1000;
v=vmin:vstep:vmax;
%% Convert t,x,v to the 1st, 2nd & 3rd dimension respectively
if iscolumn(x);
    x=x';
end;
if isrow(t);
    t=t';
end;
if isrow(v);
    v=permute(v,[1,3,2]);
else
    v=permute(v,[3,2,1]);
end;
%% Compute 3D NMO time matrix
T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,x.^2,v.^2)));
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
S(isnan(S))=0;
%% Determine Percentile
if(perc);perc=prctile(data(:),perc);else perc=max(abs(data(:)));end;
%% Semblance Figure
ddcmp=median(diff(x));if isnan(ddcmp);ddcmp=1;end;
tmp=data./perc;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
tmp=bsxfun(@times,tmp,ddcmp);
F=figure('CloseRequestFcn',@closereq,'renderer','painters');
AH(1)=subplot(1,2,1);
if(VA); % Check for variable area flag
	tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
	[Face,Vert]=VariableArea(tmp,t,x);
	PH=patch('Faces',Face,'Vertices',Vert,'FaceColor','Black','EdgeColor','none');
    set(PH,'BusyAction','cancel','HitTest','off','Interruptible','off');
	axis ij;
	hold on;
end;
tmp=bsxfun(@plus,x,tmp);
IH(1)=plot(reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'k');
set(IH(1),'BusyAction','cancel','HitTest','off','Interruptible','off');
if(VA);
	hold off;
else
	axis ij;
end;
axis([x(1)-ddcmp/2,x(end)+ddcmp/2,t(1),t(end)]);
xlabel('Offset []');ylabel('Two-Way-Traveltime []');title(sprintf('CMP Gather'));
AH(2)=subplot(1,2,2);I(2)=imagesc(permute(v,[1,3,2]),t,S,[0,1]);
hCB=InteractiveColorbar;
xlabel('NMO Velocity');title('Semblance Plot');
% Set UserData, Colormap & Sym
set(AH(1),'UserData',false);set(AH(2),'UserData',true);
hCB.ColorMap='jet';
% Set Pointer & Window Functions
if(isa(F(1),'matlab.ui.Figure'));
    set(F(1),'Pointer','CrossHair');
else
    set(F(1),'Pointer','Full CrossHair');
end;
set(F,'WindowButtonDownFcn',@clickFcn,'WindowButtonUpFcn',@unclickFcn);
set(AH(1),'buttondownfcn',@NMOcorrGather);
% Set up cursor lines
hNMO=line(x,NaN(1,numel(x)),'Color','red','LineWidth',1.5,'Parent',AH(1)); %NMO Line Handle
hVM =line(NaN(3,1),NaN(3,1),'Color','red','LineWidth',1.5,'Parent',AH(2)); %Velocity Line Handle
indord=0; %Vector last clicked velocity profile indeces so that they can be removed
%% Window Functions
function clickFcn(varargin)
    if strcmpi(get(gcbf,'SelectionType'),'normal');
        if I(2)~=gco;
            set(hNMO,'YData',NaN(1,numel(x)));
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
    set(hNMO,'YData',sqrt(pt(3).^2+(x/pt(1)).^2));
end
function unclickFcn(varargin)
    if I(2)==gco;
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
function NMOcorrGather(varargin)
	if(strcmpi(get(gcbf,'SelectionType'),'alt'));
		tmp=get(hVM,'XData');
		if(any(isnan(tmp)));
			if(~all(isnan(tmp)));
				fprintf(2,'semblanceWiggle: Warning: NaNs detected in velocity picks.\m')
			end;
			return;
		end;
		V=interp1(get(hVM,'YData'),get(hVM,'XData'),t);
		T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,x,V).^2)); %NMO Times
        tmp=data;
		for j=1:numel(x);
			tmp(:,j)=interp1(t,data(:,j),T(:,j),'linear',0); %NMO correct
		end;
        tmp=StretchMute(tmp,x,t,V,SMute);
        tmp=tmp./perc;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
        tmp=bsxfun(@times,tmp,ddcmp);
		if(VA);
			tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
			[Face,Vert]=VariableArea(tmp,t,x);
			set(PH,'Faces',Face,'Vertices',Vert);
		end;
		tmp=bsxfun(@plus,x,tmp);
		set(IH(1),'XData',reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'YData',reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]));
        set(hNMO,'XData',x,'YData',nan(numel(x),1));
	elseif(strcmpi(get(gcbf,'SelectionType'),'extend'));
		tmp=data./perc;tmp(tmp>1)=1;tmp(tmp<-1)=-1;
        tmp=bsxfun(@times,tmp,ddcmp);
		if(VA);
			tmp(1,:)=0;tmp(end,:)=0; %Enforces start and end point are 0.
			[Face,Vert]=VariableArea(tmp,t,x);
			set(PH,'Faces',Face,'Vertices',Vert);
		end;
		tmp=bsxfun(@plus,x,tmp);
		set(IH(1),'XData',reshape([tmp;nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]),'YData',reshape([repmat(t,[1,size(tmp,2)]);nan(1,size(tmp,2))],[numel(tmp)+size(tmp,2),1]));
	end;
end
function closereq(~,~)
    % Check if user wishes to export velocity model
	selection = questdlg('Export Velocity Model?','Export Velocity Model','Yes','No','Cancel','No'); 
	switch selection, 
        case 'Yes',
            if all(isnan(get(hVM,'YData')));
                fprintf(1,'semb: No velocity picks present. Not exporting velocity model.\n');
            else
                if evalin('base','exist(''t'',''var'');');
                    i=0;
                    while true;
                        if ~evalin('base',sprintf('exist(''t%d'',''var'');',i));
                            assignin('base',sprintf('t%d',i),get(hVM,'YData'));
                            fprintf(1,'semb: Saved time picks as t%d.\n',i);
                            break;
                        else
                            i=i+1;
                        end;
                    end;
                else
                    assignin('base','t',get(hVM,'YData'));
                    fprintf(1,'semb: Saved time picks as t.\n');
                end;
                if evalin('base','exist(''v'',''var'');');
                    i=0;
                    while true;
                        if ~evalin('base',sprintf('exist(''v%d'',''var'');',i));
                            assignin('base',sprintf('v%d',i),get(hVM,'XData'));
                            fprintf(1,'semb: Saved velocity picks as v%d.\n',i);
                            break;
                        else
                            i=i+1;
                        end;
                    end;
                else
                    assignin('base','v',get(hVM,'XData'));
                    fprintf(1,'semb: Saved velocity picks as v.\n');
                end;
            end;
            delete(gcf);
        case 'No'
            delete(gcf);
        case 'Cancel'
      return 
	end;
end
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
    if(SMute>0);data((bsxfun(@rdivide,x,t0.*v).^2)/2>SMute)=0;end;
end
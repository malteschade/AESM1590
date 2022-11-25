function sembStk(Data,hdr,geo,Velocity,CMP,V,dCMP,ExpSmooth)
% sembStk(Data,hdr,geo,Velocity,CMP) creates an interactive weighted
% semblance plot of the seismic data Data using a default gray-scale
% interactive colorbar.
%
% INPUTS:
%   Data:       Numeric nt x ntr matrix containing seismic data.
%   hdr:        Numeric  5 x ntr matirx containg trace headers.
%   geo:        Numeric 6 element vector containg geometry information.
%   Velocity:   Numeric n element vector containing NMO velocities.
%   CMP:        Integer index of first CMP location.
%   V:          (DEFAULT:0) Numeric nt x nCMP matrix containing intial
%               velocity model. Useful for updating velocity model.
%   dCMP:       (DEFAULT:10) CMP increment when stepping through gathers.
%   ExpSmooth:  (DEFAULT:50) Exponential semblance smoothing factor.
%
% Figure 1 shows the CMP gather and associated semblance plot. To pick
% velocities left click and drag inside the semblance plot till a good NMO
% velocity has been found. To support picking a red NMO curve will appear
% on top of the CMP gather. To undo the last pick right click in semblance
% plot. As you pick points the red line in the semblance plot will indicate
% the picked velocity model. The blue line (if present) indicates the
% interpolated velocity model based on previous picks. You can use the left
% or right arrow to move between CMP gathers according to the stepping
% increment dCMP.
%
% Figure 2 shows the picked velocity model and the associated NMO stack (if
% data has been stacked). Left clicking inside either plot will change the
% CMP location. Right clicking inside the stack will cause the stack to be
% generated or updated using the given velocity model.
%
% Upon closing either figure one is prompted to save the velocity model and
% stack to the workspace. They will be saved in V, V0, V1, V2 ...
% and Stk, Stk0, Stk1, Stk2. Variables will not be overwritten.
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
% October 2015
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
%
% nestedSortStruct
% Jake Hughey
% http://www.mathworks.com/matlabcentral/fileexchange/28573-nested-sort-of-structure-arrays/

%% Check inputs
if ~isnumeric(Data)||~ismatrix(Data)||numel(Data)<2;
    fprintf(2,'sembStk: Data is not a numeric matirx of size greater than 1.\n');
    return;
end;
if ~isnumeric(hdr)||~ismatrix(hdr)||size(hdr,1)~=5;
    fprintf(2,'sembStk: hdr is not a numeric matirx of size 5xn.\n');
    fprintf(2,'sembStk: Is hdr a header matrix?\n');
    return;
end;
if ~isnumeric(geo)||~isvector(geo)||numel(geo)~=7;
    fprintf(2,'sembStk: geo is not a numeric vector with six elements.\n');
    fprintf(2,'sembStk: Is geo a geometry matrix?\n');
    return;
end;
if ~isnumeric(Velocity)||~isvector(Velocity);
    fprintf(2,'sembStk: velocity is not a numeric vector.\n');
    return;
else
    if iscolumn(Velocity);
        Velocity=permute(Velocity,[3,2,1]);
    else
        Velocity=permute(Velocity,[1,3,2]);
    end;
end;
nCMP=numel(unique(hdr(5,:)));
if ~isnumeric(CMP)||~isscalar(CMP)||mod(CMP,1)~=0||CMP<0||CMP>nCMP;
    fprintf(2,'sembStk: CMP is not a positive integer scalar.\n');
    fprintf(2,'sembStk: Is CMP the index of the 1st CMP to be plotted?.\n');
    return;
end;
if exist('V','var');
    if ~isnumeric(V)||~ismatrix(V)||size(V,1)~=size(Data,1)||size(V,2)~=nCMP;
        fprintf(2,'sembStk: V is not a correct size numeric matrix.\n');
        return;
    end;
else
    V=nan(size(Data,1),nCMP);
end;
if exist('dCMP','var');
    if ~isnumeric(dCMP)||~isscalar(dCMP)||mod(dCMP,1)~=0;
        fprintf(2,'sembStk: dCMP is not a integer scalar.\n');
        return;
    end;
else
    dCMP=10;
end;
if exist('ExpSmooth','var');
    if ~isnumeric(dCMP)||~isscalar(dCMP);
        fprintf(2,'sembStk: ExpSmooth is not a numeric scalar.\n');
        return;
    end;
else
    ExpSmooth=50;
end;
%% Convert hdr to TrcHdr structure
TrcHdr=struct('cdp',num2cell(hdr(5,:)),'offset',num2cell(hdr(2,:)));
%% Generate Time Vector
dt=geo(1)/1000;
t=(0:geo(2)-1)'*dt;
%% Sort & Bin CMPs
[Data,TrcHdr]=SeisSort(Data,TrcHdr,{'cdp','offset'},1);
[Data,TrcHdr]=SeisBin(Data,TrcHdr,'cdp');
%% Get CMP Locations
cmp=cellfun(@(a)a(1).cdp,TrcHdr);
ns=num2str(numel(num2str(max(cmp))));
CMPproc=false(1,numel(cmp));
%% Set-Up Velocity Model + Stack
V=nan(numel(t),numel(Data)); %Velocity Model
Stk=zeros(numel(t),numel(Data)); %NMO - Stack
%% Plot 1st CMP
% Extract Offsets
x=[TrcHdr{CMP}.offset];
% Compute 1st Semblance
S=semb(Data{CMP},x,t,Velocity,ExpSmooth);
% Figure
F(1)=figure('CloseRequestFcn',@closereq);
AH(1)=subplot(1,2,1);IH(1)=imagesc(x,t,Data{CMP},[-1,1].*max(abs(Data{CMP}(:))));
hCB(1)=InteractiveColorbar;
TH(1)=title(sprintf('CMP Gather: %d',CMP));
xlabel('Offset [m]');ylabel('Two Way Traveltime [s]');
AH(2)=subplot(1,2,2);IH(2)=imagesc(permute(Velocity,[1,3,2]),t,S,[0,1]);
hCB(2)=InteractiveColorbar;
xlabel('NMO Velocity');ylabel('Two Way Traveltime [s]');
TH(2)=title('Semblance Plot');
% Set UserData, Colormap & Sym
set(AH(1),'UserData',false);set(AH(2),'UserData',true);
hCB(1).ColorMap='gray';hCB(2).ColorMap='gray';hCB(1).Sym=true;
% Set Pointer & Window Functions
set(F(1),'Pointer','CrossHair');
set(F(1),'WindowButtonDownFcn',@clickFcn);
set(F(1),'WindowButtonUpFcn',@unclickFcn);
set(F(1),'KeyPressFcn',@keyDownListener);
% Set up cursor lines
hNMO=line(x,NaN(1,numel(x)),'Color','red','LineWidth',1.5,'Parent',AH(1)); %NMO Line Handle
hVM =line(NaN(3,1),NaN(3,1),'Color','red','LineWidth',1.5,'Parent',AH(2)); %Velocity Line Handle
indord=0; %Vector last clicked velocity profile indeces so that they can be remove
% Set up auxillary lines
hVMo=line(V(:,CMP),t,'Color','blue','LineWidth',1.5,'Parent',AH(2)); %OriginalVelocity Line Handle
%% Plot Velocity Model & NMO Stack
F(2)=figure('CloseRequestFcn',@closereq);
colormap(gray);
AH(3)=subplot(2,1,1);IH(3)=imagesc(cmp,t,V,[min(Velocity),max(Velocity)]);
xlabel('CMP Location [m]');ylabel('Two Way Traveltime [s]');title('RMS Velocity Model');
AH(4)=subplot(2,1,2);IH(4)=imagesc(cmp,t,Stk);
xlabel('CMP Location [m]');ylabel('Two Way Traveltime [s]');title('Stack');
%dCMPh=uicontrol('Style','Edit','String',dCMP,'Position',[0,0,23,23],...
%    'Backgroundcolor',[1,1,1],'Foregroundcolor',[0,0,0]);
% Set up auxillary lines
hCMPL(1)=line([cmp(CMP),cmp(CMP)],[t(1)-0.5.*dt,t(end)+0.5.*dt],'Color','red','LineWidth',1.5,'Parent',AH(3)); %Velocity Model CMP-Line
hCMPL(2)=line([cmp(CMP),cmp(CMP)],[t(1)-0.5.*dt,t(end)+0.5.*dt],'Color','red','LineWidth',1.5,'Parent',AH(4)); %Stk Model CMP-Line
% Set Window Functions
drawnow;
set(F(2),'WindowButtonUpFcn',@selectCMP);
%% Get Java Window Handles
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%TMP=get(F(1),'JavaFrame');jSemb=TMP.fFigureClient.getWindow;
%TMP=get(F(2),'JavaFrame');jStk=TMP.fFigureClient.getWindow;
%TMP=get(F(1),'JavaFrame');jSemb=TMP.fHG1Client.getWindow;
%TMP=get(F(2),'JavaFrame');jStk=TMP.fHG1Client.getWindow;
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%% Blurring Window

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
        % Update Stk
        T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,[TrcHdr{CMP}.offset],V(:,CMP)).^2)); %NMO Times
        for j=1:numel(TrcHdr{CMP});
            Stk(:,CMP)=Stk(:,CMP)+interp1(t,Data{CMP}(:,j),T(:,j),'linear',0); %NMO correct
        end;
        set(IH(4),'CData',Stk);
        % Update Process CMP
        CMPproc(CMP)=true;
    end
    function UpdatePlot
        set(IH(1),'XData',[TrcHdr{CMP}.offset],'CData',Data{CMP});
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
        tmp=str2double(get(dCMPh,'String'));
        if isnan(tmp);return;end;
        switch event.Key
            case 'leftarrow'; %Show previous record
                if CMP>1;
                    UpdateVelModelStk;
                    if CMP>tmp;
                        CMP=CMP-tmp;
                        UpdatePlot;
                    else
                        CMP=1;
                        UpdatePlot;
                    end;
                end;
            case 'downarrow'; %Take a larger step back
                if CMP>1;
                    tmp=10.*tmp;
                    UpdateVelModelStk;
                    if CMP>tmp;
                        CMP=CMP-tmp;
                        UpdatePlot;
                    else
                        CMP=1;
                        UpdatePlot;
                    end;
                end;
            case 'rightarrow'; %Show next record
                if CMP<numel(cmp);
                    UpdateVelModelStk;
                    if CMP+tmp<=numel(cmp);
                        CMP=CMP+tmp;
                        UpdatePlot;
                    else
                        CMP=numel(cmp);
                        UpdatePlot;
                    end;
                end;
            case 'uparrow'; %Take a large step forward
                if CMP<numel(cmp);
                    tmp=10.*tmp;
                    UpdateVelModelStk;
                    if CMP+tmp<=numel(cmp);
                        CMP=CMP+tmp;
                        UpdatePlot;
                    else
                        CMP=numel(cmp);
                        UpdatePlot;
                    end;
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
        % Update velocity model
        UpdateVelModelStk
        % Disable Figures
%        jSemb.setEnabled(false);
%        jStk.setEnabled(false);
        % Create Wait Bar
        WB=waitbar(0,'CMP Stacking','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(WB,'canceling',0)
        % CMP Stack
        for i=1:numel(cmp);
            if getappdata(WB,'canceling');
                break;
            end;
            T=sqrt(bsxfun(@plus,t.^2,bsxfun(@rdivide,[TrcHdr{i}.offset],V(:,i)).^2)); %NMO Times
            Stk(:,i)=0;
            for j=1:numel(TrcHdr{i});
                Stk(:,i)=Stk(:,i)+interp1(t,Data{i}(:,j),T(:,j),'linear',0); %NMO correct
            end;
            waitbar(i/numel(cmp),WB);
        end;
        % Enable Figures
%        jSemb.setEnabled(true);
%        jStk.setEnabled(true);
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
                if any(isnan(V));
                    fprintf(1,'semb: No velocity picks present. Not exporting velocity model.\n');
                else
                    if evalin('base','exist(''V'',''var'');');
                        i=0;
                        while true;
                            if ~evalin('base',sprintf('exist(''v%d'',''var'');',i));
                                assignin('base',sprintf('v%d',i),get(IH(3),'CData'));
                                fprintf(1,'semb: Saved velocity picks as v%d.\n',i);
                                break;
                            else
                                i=i+1;
                            end;
                        end;
                    else
                        assignin('base','V',get(IH(3),'CData'));
                        fprintf(1,'semb: Saved velocity model as V.\n');
                    end;
                end;
                if any(isnan(Stk));
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
end
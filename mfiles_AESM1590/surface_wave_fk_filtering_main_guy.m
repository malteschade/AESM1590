% Main script to use for F-Kx filtering of surface waves
% based on Assignment 3.1
%
% Input file should be created from SU format by stripping the header
% producing thus a binary file.

clear variables
close all
% clc 

%% Input parameters of the data to be filtered
filename_data = 'refl_3layers_fp50_dx0p5_500_rvz_tapered.bin';
nt = 1001;  % number of time samples
dt = 0.001; % time sampling in seconds
nx = 401;   % number of spatial samples (traces)
dx = 2.5;   % spatial sampling in meters, so distance between two neighboring traces
xfirsttrace = 510; % position of first receiver along the line, for plotting only

fprintf(['   --- \n'])
fprintf(['   Number of time samples is ',num2str(nt,'%5.0f'),' \n'])
fprintf(['   Time sampling is ',num2str(dt,'%5.4f'),' (s) \n'])
fprintf(['   Number of spatial samples (traces) is ',num2str(nx,'%5.0f'),' \n'])
fprintf(['   Spatial sampling is ',num2str(dx,'%5.1f'),' (m) \n'])

% Parameters used in (frequency,x-space) domain
fcut = 200; % frequencies higher than this are expected to be zero and are discarded

% Parameters used in (frequency,horizontal-wavenumber) domain
ntaper     = 10;  % points to taper along wavenumber axis before setting values to zero
taperexcur = 20;  % points to left (or right) of surface-wave energy maximum to start taper
v_sw       = 450; % estimated highest surface-wave velocity

ltxplot = 1; % plot t-x results?  ltxplot=0 is NO; otherwise YES
lfkplot = 1; % plot F-Kx results? lfkplot=0 is NO; otherwise YES
%%

% Allocating memory for data
data = zeros(nt,nx);

% Read binary file into Matlab using function readbinsu.m
data(:,:) = readbinsu(filename_data,nx,nt);

%% Plotting t-x records
if ltxplot ~= 0
    % Potting parameters 
    fs            = 14;  % Fontsize
    lw            = 2;   % Linewidth
    cmap          = gray(256); % colour map for plotting common-source gather
    dispscale     = 0.1; % scaling colorbar to clip data to display weaker arrivals
    dispscalefilt = 0.1; % scaling colorbar to clip data to display filtered result

    % Defining time and space axes
    timeaxis  = linspace(0,nt-1,nt)*dt;
    spaceaxis = linspace(0,nx-1,nx)*dx+xfirsttrace;

    figure;
    imagesc(spaceaxis,timeaxis,data(:,:));
    colormap(cmap)
    colorbar
    caxis([dispscale*min(min(data(:,:))) dispscale*max(max(data(:,:)))])
    xlabel('Horizontal Distance (m)','Fontsize',fs)
    ylabel('Two-way traveltime (s)','Fontsize',fs)
    title('Recorded data (reflections + surface waves)','Fontsize',fs)
    set(gca,'Fontsize',fs)
    set(gca,'LineWidth',lw)
end
%%

% (F,KX)-domain filtering using separate function fk_lintaper_guy:
fprintf(['   --- \n'])
data_swfilt = fk_lintaper_guy(data,dx,dt,nx,nt,fcut,ntaper,taperexcur,v_sw,lfkplot);
fprintf(['   --- \n'])

%% Plotting (F,Kx)-filtered records
if ltxplot ~= 0
    figure;
    imagesc(spaceaxis,timeaxis,(data_swfilt(:,:)));
    colormap(cmap)
    colorbar
    caxis([dispscalefilt*min(min(data_swfilt(:,:))) ...
           dispscalefilt*max(max(data_swfilt(:,:)))])
    xlabel('Horizontal Distance (m)','Fontsize',fs)
    ylabel('Two-way traveltime (s)','Fontsize',fs)
    title('F-Kx (surface-wave) filtered data','Fontsize',fs)
    set(gca,'Fontsize',fs)
    set(gca,'LineWidth',lw)

    fprintf('   Notice F-Kx filtering suppresses surface-wave energy, \n')
    fprintf('                while body-wave energy is preserved. \n')
    fprintf('   Also notice artefacts from the filtering. \n')
end

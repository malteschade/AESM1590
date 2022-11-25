%% 1a
clear;
close all;
clc;

a=50;
b=10;
pow = 1.5;
ra = 3:7;
cl = [10,50];


mymatrix = zeros([b,a]);

mymatrix(:) = repmat((1:a),b,1).*(1:b).';

mymatrix(:) = mymatrix.*(pow.^linspace(0,b-1,b)).';
mymatrix = mymatrix(ra,:);
[m,n] = size (mymatrix);
disp([m,n]);


imagesc(mymatrix);
clim(cl);
colorbar;

%% 1b
clear;
close all;
clc;

path = '/Users/malteschade/Desktop/Seismic Ac/AESM1590/datasets/shotground.segy';
xf = 20:40;
yf = 150:250;

[data, h_csp, geo]  = segyread(path);
data = data(yf, xf);

dt = (yf-1).*geo(1)/1000;   % time vector [s]
dx = (xf-1).*geo(6);        % space vector [m]

ylabel('Time [s]')
xlabel('Space [m]')

plotseis(data,dt, dx);

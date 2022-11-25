function [data_swfilt]=fk_lintaper_guy(data,dx,dt,nx,nt,fcut,ntaper,taperexcur,v_sw,lfkplot)
% Function applying F-Kx filtering to input data containing reflections and surface waves.
% The data are transformed to the F-Kx domain, where the data above a certain velocity line
% is set to zero, with the aim to suppress the surface waves. 
% The velocity line is defined by the estimated fastest surface waves.
% The values are smoothly going to zero using a linear taper.
%    Deyan Draganov, Version November 2020
%    Guy Drijkoningen, adapted August 2021

%% Transforming data to the (frequency,horizontal-wavenumber) domain

fprintf('   FFT from time to frequency domain ... ')
data_f   = fft(data,[],1)*dt;
df       = 1/(nt*dt); % frequency sampling in Hertz
nf       = nt;        % number of frequency samples
freqaxis = linspace(0,nf/2,nf/2+1)*df; %making the frequency axis
fprintf('done\n')

% added by Guy:
if fcut > 1/(2*dt)
    fcut = 1/(2*dt);
    fprintf('   since fcut > f_Nyquist, fcut set to Nyquist = 1//(2*dt)\n')
end
fcutel      = find(freqaxis>=fcut,1,'first');
freqaxiscut = freqaxis(1:fcutel); %the new frequency axis
data_f      = data_f(1:fcutel,:);

fprintf('   IFFT from x-space to horizontal-wavenumber domain ... ')
data_fk = fftshift(nx*ifft(fftshift(data_f,2),[],2)*dx,2);
dk      = 1/(nx*dx); % wavenumber sampling in 1/m
nk      = nx;        % number of wavenumber samples
kaxis   = linspace(-nk/2,nk/2-1,nk)*dk;
fprintf('done\n')

%% Suppressing surface waves via (F,Kx)-domain windowing

% looping over frequencies to find right horizontal wavenumber to taper 
%    to the left or right of it
data_fkfil = data_fk;
% checking which is maximum frequency (element) corresponding to 
%    maximum wavenumber (element) present
fswel_max = find(freqaxis<=(kaxis(nk)*v_sw),1,'last');
% looping only till the frequency that ensures we have a corresponding 
%    pair (frequency,wavenumber)
for freq=1:(min(fswel_max,fcutel)-1)
    % below we calculate a wavenumber value using the current frequency and
    % the given surface-wave velocity
    ksw       = freqaxiscut(freq)/v_sw;
    kswel_pos = find(kaxis>=ksw,1,'first');
    kswel_neg = find(kaxis<=-ksw,1,'last');
    %tapering values along the wavenumber axis starting from kswel_pos 
    %shifted to the left by taperexcur and going to the right
    lengthk  = nk-kswel_pos+taperexcur;
    lintaper = [linspace(1,0,ntaper) zeros(1,lengthk-ntaper)];
    data_fkfil(freq,(kswel_pos-taperexcur+1):end) = ...
       data_fk(freq,(kswel_pos-taperexcur+1):end).*lintaper;
    %tapering values along the wavenumber axis starting from kswel_neg 
    %shifted to the right by taperexcur and going to the left
    lengthk  = kswel_neg+taperexcur;
    lintaper = [zeros(1,lengthk-ntaper) linspace(0,1,ntaper)];
    data_fkfil(freq,1:(kswel_neg+taperexcur)) = ...
       data_fk(freq,1:(kswel_neg+taperexcur)).*lintaper;
end

% Transforming data back to time-space domain
fprintf('   FFT from horizontal-wavenumber to x-space domain ... ')
tempfdata = fftshift(fft(fftshift(data_fkfil,2),[],2)*dk,2);
fprintf('done\n')

% Adding zeroes instead of the removed frequencies to obtain time-domain
%    data as long as before the transformations.
addfreq             = zeros(nf,nx);
addfreq(1:fcutel,:) = tempfdata;

% Transforming data back to the time-space domain
fprintf('   IFFT from frequency to time domain ... ')
data_swfilt = 2*real(nf*ifft(addfreq,[],1)*df);
fprintf('done\n')

%% Some plotting paramaters
if lfkplot ~= 0
    fs = 14; % Fontsize
    lw = 2;  % Linewidth

    figure;
    imagesc(kaxis,freqaxiscut,abs(data_fk(:,:)))
    colorbar
    xlabel('Horizontal Wavenumber (1/m)','Fontsize',fs)
    ylabel('Frequency (Hz)','Fontsize',fs)
    title('Data in F-Kx domain','Fontsize',fs)
    set(gca,'Fontsize',fs)
    set(gca,'LineWidth',lw)

    figure;
    imagesc(kaxis,freqaxiscut,abs(data_fkfil(:,:)))
    colorbar
    xlabel('Horizontal Wavenumber (1/m)','Fontsize',fs)
    ylabel('Frequency (Hz)','Fontsize',fs)
    title('F-Kx filtered data in F-Kx domain','Fontsize',fs)
    set(gca,'Fontsize',fs)
    set(gca,'LineWidth',lw)
end 

end



function [t,v]=getVelocityModel(hVM)
% [t,v]=getVelocityModel(hVM) returns the RMS velocity model picked using
% hVM=semb(data,x,t,v). hVM is the handle to the velocity model in the semb
% plot. This function will return the time picks t and the corresponding
% velocity picks v in row vectors.
%
% NOTE: This function must be used while the semb plot window is still
% open.
t=get(hVM,'YData');
v=get(hVM,'XData');
end
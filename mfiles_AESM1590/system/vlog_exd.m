function vt = vlog_exd(v,t,geo)
% ==== system-script, do not execute manually ====
%
% VLOG_EXD - generates a velocity-time table.
%
% This function generates a 1D vertical velocity vs. time log.
% This version is written to be called from GENERATEVMODEL.

% call for geometry
dt=geo(1);
nt=geo(2);

t2=0:dt:(dt*nt-dt);v2=pwlint(t,v,t2);
vt = v2';

function [t, npq] = NPQformat(mdate, f0, fm, Fo, Fm, Opts)

% Use primary FRRF outputs to compute NPQ secondary parameters 
%
% INPUTS: 
% mdate - sample time points in matlab datenum format
% f0 - minimum Chl fluorescence time series
% fm - maximum Chl fluorescence time series
% Fo - dark regulated minimum Chl fluorescence, scalar
% Fm - dark regulated maximum Chl fluorescence, scalar
% Opts - 'sv', 'nsv', 'qN'
%
% OUTPUTS:
% t - sample time points in minutes
% npq - npq calculated from fluorescence values. NPQ eqn determined by Opts

% Safety Checks:

% Format sampling time array
ts = mdate - mdate(1);
dv = datevec(ts);
t = dv(:,4)*60 + dv(:,5) + dv(:,6)/60;


% Compute NPQ

Fv = Fm - Fo;

npq_sv = (Fm - fm)./fm;
npq_nsv = f0./(fm-f0);
qN = 1 - (fm - Fo)./Fv;

npq_set = find(strcmpi(Opts,[{'nsv'},{'sv'},{'qN'}]));
switch npq_set
    case 1 
        npq = npq_nsv;
    case 2
        npq = npq_sv;
    case 3
        npq = qN;
end



end
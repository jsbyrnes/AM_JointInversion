%This script make a synthetic receiver function filetered at a given band.
%That's it.

clear, clc, close all

model.theta = [ 0 0 0 ];%tilt from verticle of the normal vector
model.phi = [ 0 0 0 ];%strike clockwise from north
model.z = [ 30 200 500 ]; %km
model.vp = [ 7 8 7.5 ]; %km/s
model.A = [ 0 0 0 ]; %no units
model.B = [ 0 0 0 ]; %no units
model.vpvs = [ 1.8 1.8 1.8 ]; %no units
model.C = [ 0 0 0 ]; %no units
model.rho = nafedrake_rho(model.vp); %g/cm^3

baz = 0;

phase = 'SV';

slow = .1;

filter = [0.03 0.1];

dt = .05;

[P_comp, SV_comp, ~] = anirec(phase, dt, slow, baz, model, 1);

[trace_phase, ~] = multitaper2rf(P_comp, SV_comp, dt, 7.5, 400, 2.5, 3, phase, filter);

%compare with a water level deconvolution
trace_waterlevel = water_level_decon(SV_comp, P_comp, 0.01, 7.5, phase, filter, dt);

%compare with the original green's function
nsamples = length(SV_comp);
        
SV_comp = flipud(SV_comp);
P_comp = flipud(-1*P_comp);

[~, arrival_index] = max(SV_comp);
        
SV_comp = SV_comp(arrival_index - 7.5/dt + 1:end);
P_comp = P_comp(arrival_index - 7.5/dt + 1:end);

%how many zeros are needed to get back original length?
new_nsamples = length(SV_comp);
        
SV_comp = [SV_comp; zeros(nsamples - new_nsamples, 1)];
P_comp = [P_comp; zeros(nsamples - new_nsamples, 1)];

%filtered P and normlaize to SV
P_comp_filt = bandpassfilt(P_comp, 0.05, filter(2), filter(1));
SV_comp_filt = bandpassfilt(SV_comp, 0.05, filter(2), filter(1));
amp = max(SV_comp_filt);
P_comp_filt = P_comp_filt/amp;

t = ( 0:dt: (length(trace_waterlevel) - 1)*dt ) - 7.5 ;

figure(2), hold on
plot(t, trace_phase, 'r'); 
plot(t, trace_waterlevel, 'b');
plot(t, P_comp_filt, 'k');
xlabel('Time (seconds)');
xlim([0 50]);
legend('Multitaper receiver function', 'Waterlevel receiver function', 'Filtered Green''s function');

%print(2, 'comparison.png', '-dpng');

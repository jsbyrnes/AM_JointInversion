function [outTraces] = wfResample_jsb(inTraces,Parameters)
%USAGE: [outTraces] = wfResample(inTraces,newSR);
% inTraces, outTraces are irisFetch type structures of arrays
% newSR is the desired sampling rate in sps.
% in outTraces all the .data will have the same number of samples
% the fields sampleCount and sampleRate will be updated accordingly.
% I should add a comment saying the traces were resampled.
% any traces that are shorter than what they are supposed to be, will be
% axed.

%loop through traces and resample as needed

outTraces=inTraces;
for k=1:length(inTraces)

    %outTraces(k).data = interp1(linspace(Parameters.t(1), Parameters.t(end), length(outTraces(k).data)),...
    %    inTraces(k).data, Parameters.t, 'pchip');

    %more formal but double filters the data
    oldSR1 = round(outTraces(k).sampleRate*1000); %keep 3 decimals
    newSR1 = round(Parameters.sample_rate*1000);

    outTraces(k).data        = resample(inTraces(k).data,newSR1,oldSR1);

    outTraces(k).sampleRate  = Parameters.sample_rate;
    outTraces(k).sampleCount = length(outTraces(k).data);

end

%OK, now go back and make sure that the sampleCounts are all the same.
usc=unique([outTraces.sampleCount]);

if length(usc)==1
    %if they are already uniform
    return
elseif length(usc)>2 || (max(usc)-min(usc)) > 1
    %the expected scenario is that there will be no more than two sample
    %counts and they will be X and X+1, if that's not the case, display
    %this message
    warning('something is not as expected, check things carefully')
    keyboard
end

%if you've made it this far then there's more than one sampleCount
%go through and fix that
minSc=min(usc);

for k=1:length(outTraces)
    if outTraces(k).sampleCount > minSc
        %chop off the rest of the samples (which should be only one, but
        %I'm making it general anyways
        outTraces(k).data = outTraces(k).data(1:minSc);
        outTraces(k).sampleCount = length(outTraces(k).data);
    end
end

%todo: I should add a 'comment'
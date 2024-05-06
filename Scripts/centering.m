function metric = centering(ts, fs, tau, ind)

    ts = delay_continuous(ts, fs, tau);

    metric = -ts(ind);% + ts(ind-1) + ts(ind+1);%make it a real spike!

end
function [pick, snr] = sta_lta(data, t0, twin, sta, lta, t)

    dt = t(2) - t(1);

    [~, search_ind1] = min(abs(t - (t0 - twin)));
    [~, search_ind2] = min(abs(t - (t0 + twin)));

    search_ind = search_ind1:search_ind2;

    for k = 1:length(search_ind)

        ltaind = (search_ind(k) - round(lta/dt));

        if any(ltaind<1) || (search_ind(k) + round(sta/dt)) > length(data)

            ratio(k) = 0;

        else

            ratio(k) = mean(data(search_ind(k):(search_ind(k) + round(sta/dt))).^2)...
                /mean(data(ltaind:search_ind(k)).^2);

        end
   
    end

    [snr, pick ] = max(ratio);

    pick = t(search_ind(pick));

end
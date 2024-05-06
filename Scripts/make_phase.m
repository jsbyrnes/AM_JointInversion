function phase = make_phase(t, a, mu, w, n)

    %w = w*(n+1);

    %first, set this to zero mean and 1 width
    t = (t - mu)/w;

    %first make the phase
    phase = exp(-0.5*t.^2);

    if n == 0

        poly = 1;

    elseif n == 1

        poly = -t;

    elseif n == 2

        poly = t.^2 - 1;

    elseif n == 3

        poly = -1*(t.^3 - 3*t);

    elseif n == 4

        poly = t.^4 - 6*t.^2 + 3;

    elseif n == 5

        poly = -1*(t.^5 - 10*t.^3 + 15*t);

    elseif n == 6

        poly = t.^6 - 15*t.^4 + 45*t.^2 - 15;

    elseif n == 7

        poly = -1*(t.^7 - 21*t.^5 + 105*t.^3 - 105*t);

    elseif n == 8
        
        poly = t.^8 - 28*t.^6 + 210*t.^4 - 420*t.^2 + 105;

    elseif n == 9

        poly = -1*(t.^9 - 36*t.^7 + 378*t.^5 - 1260*t.^3 + 945*t);

    elseif n == 10

        poly = t.^10 - 45*t.^8 + 630*t.^6 - 3150*t.^4 + 4725*t.^2 - 945;

    end

    phase = a*poly.*phase;

end

function [ G ] = drive_mtm(C1, C2, C3 , C, E, V, sps)
%DRIVE_MTM Do multi taper decon
%   C1 is the source

    n = length(C1);

    S = C1;
    
    [C1, FFT_C1] = mtmdeconSP(S', C1', C.ntap, C.nmtw, E, V);
    [C2, FFT_C2] = mtmdeconSP(S', C2', C.ntap, C.nmtw, E, V);
    [C3, FFT_C3] = mtmdeconSP(S', C3', C.ntap, C.nmtw, E, V);
    
    %C1 = C2;
    %FFT_C1 = FFT_C2;

    shift = round(C.rf_shift*sps);
    
    C1 = circshift(C1, [shift 0]);
    C2 = circshift(C2, [shift 0]);
    C3 = circshift(C3, [shift 0]);
    
    G.C1 = C1;
    G.FFT_C1 = FFT_C1;
    G.C2 = C2;
    G.FFT_C2 = FFT_C2;
    G.C3 = C3;
    G.FFT_C3 = FFT_C3;

end
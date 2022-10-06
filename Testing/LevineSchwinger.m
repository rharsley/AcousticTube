%%Get the radiation impedance of an unflanged pipe based on Levine
%%and Schwinger
%%@param rho - density of air [kg/m3]
%%@param c - speed of sound in air [m/s]
%%@param a - tube radius [m]
%%@param Nka - number of points to calculate
%%@param fmax - the maximum frequency to calculate up to
%%@return ka - wave number
%%@return ZLS - the radiation impedance
%%@return RLS - the radiation reflection function
%%@Note This is only valid up to ka = 3.82, function will return an 
%% error above this.

function [ka, ZLS, RLS] = LevineSchwinger(rho, c, a, Nka, fmax)
    kmax = 2 * pi * fmax * a / c;
    if kmax > 3.82
        error("Wavenumber greater than 3.82");
    end
    
    ka = linspace(0, kmax, Nka);
    magR = zeros(1, Nka);
    l = zeros(1, Nka);
    deltax = 0.25 * (ka(2) - ka(1));
    xinf = [0.00000001 : 0.001 : 1];
    for nn = 0 : 20%powers of 10
        xinf = [xinf, [1 : 0.001 : 9] * 10^nn];
    end
    
    numinf = log(0.5 ./ (besseli(1, xinf, 1) .*  besselk(1, xinf, 1)));%besseli and besselk normalised but normalisation constant cancels in this case
    numinf(1) = 0;
    
    LSprog = Nka /10;
    LSprog0 = LSprog;
    
    for nn = 1 : Nka
        x = [deltax : deltax : ka(nn) - deltax];
        num = atan2(besselj(1, x), -bessely(1, x));
        den = (x .* sqrt(ka(nn)^2 - x.^2));
        f = num ./ den;
        if(nn != 1)
             magR(nn) = exp(-(2 * ka(nn) / pi) * trapz(x, f));
         else
             magR(1) = 1;
         endif
         
        
        num = log(pi * besselj(1, x) .* sqrt(besselj(1, x).^2 + bessely(1, x).^2));
        den = x .* sqrt(ka(nn)^2 - x.^2);
        f = num ./ den;
        
        deninf = xinf .* sqrt(xinf.^2 + ka(nn)^2);
        finf = numinf ./ deninf;
        if(nn != 1)
            l(nn)= (1 / pi) * (trapz(x, f) + trapz(xinf, finf));
        else
            l(1) = 0.6133;
        endif
        
        if nn > LSprog
            {LSprog / Nka}
            LSprog = LSprog + LSprog0;
        end
    end
    RLS = -magR .* exp(-2i * ka .* l);
    ZLS = rho * c * (1 + RLS) ./ (1 - RLS);
% Defines and normalizes eigenfunction
function f = phibar(t, w, tf, c)

    num = w*cos(w*t)/c + sin(w*t);
    den = 0.5*(1+w^2/c^2)*tf + (w^2/c^2-1)*sin(2*w*tf)/(4*w) + (1-cos(2*w*tf))/(2*c);
    den = sqrt(den);
    
    %f = den./num;
    f = num ./ den;
end

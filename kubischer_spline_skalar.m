function [ a, b, c ] = kubischer_spline_skalar( i, ddot_p, p, h )
% Die Funktion berechnet fuer jedes Intervall i die Paramter der kubischen
% Splines a_i, b_i, und c_i
% i         ... aktuelles Intervall
% ddot_p    ... Vektor 
% p         ... Vektor mit den Stuetzstellen
% h         ... aequidistanter h-Vektor

a = ( ddot_p(i+1) - ddot_p(i) ) / ( 6*h(i) );
b = ddot_p(i) * 1/2.0;
c = ( p(i+1) - p(i) ) / h(i) - h(i) * ( ddot_p(i+1) + 2*ddot_p(i) ) / 6;
end


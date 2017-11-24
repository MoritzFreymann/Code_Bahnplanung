function [ S, dot_S, ddot_S, T ] = p2p_kubisch( W_stuetz, T_ges, delta_T )
% Erzeugt aus N_I Stuetzpunkten in W_stuetz je mit Anfangs- und Endpunkt N_I-1 Trajektorien je in Form eines kubischen Polynoms
% S         := Trajektorie auf Positionsebene
% dot_S     := Trajektorie auf Geschwindigkeitsebene
% ddot_S    := Trajektorie auf Beschleunigungsebene
% T         := Zeitvektor der Trajektorie

% W_stuetz  := Stuetzpunkte
% T_ges     := Dauer der Bewegung/Interpolation
% delta_T   := Taktzeit

% Anzahl der Freiheitsgrade
N_Q       = size( W_stuetz,1 );

% Anzahl der Intervalle
N_I       = size( W_stuetz,2 ) -1;

% Zeitintervall fuer Interpolation
T_I       = T_ges/N_I;        % Zeit fuer ein Teilstueck
t_I       = 0:delta_T:T_I;      % Zeitvektor fuer ein Teilstueck
N_T_I     = length(t_I);        % Anzahl der Zeitpunkte eines Teilstuecks

%% Berechnung der Trajektorie

% Initialisierung fuer Teilstuecke
S_I       = zeros( N_Q, N_T_I-1 );
dot_S_I   = zeros(size(S_I));
ddot_S_I  = zeros(size(S_I));

% Initialisierung fuer Gesamttrajektorie
S         = [];
dot_S     = [];
ddot_S    = [];
T         = [];

%% --- ARBEITSBEREICH: ------------------------------------------------

% Schleife ueber Intervalle
for i=1:N_I
    % Berechnung Koeffizienten
    a = ( -2/T_I^3 ) * ( W_stuetz(:,i+1) - W_stuetz(:,i) );
    b = ( 3/T_I^2 ) *( W_stuetz(:,i+1) - W_stuetz(:,i) );
    d = W_stuetz(:,i);
    
    % Schleife ueber Punkte des Intervalls 
    % (letzter Punkt wird weggelassen, da dieser der erste des naechsten
    % Intervalls ist)
    for k=1:N_T_I-1
        S_I(:,k)       = a * t_I(k)^3 + b * t_I(k)^2 + d;
        dot_S_I(:,k)   = 3*a * t_I(k)^2 + 2*b * t_I(k);
        ddot_S_I(:,k)  = 6*a * t_I(k) + 2*b;
    end
    
    % Bestimme Indizes zur Speicherung des Intervallergebnisses
    first = 1 + (i-1) * (N_T_I-1);
    last = N_T_I-1 + (i-1) * (N_T_I-1);
    
    % Erzeuge Trajektorie 
    S(:,first:last) = S_I;
    dot_S(:,first:last) = dot_S_I;
    ddot_S(:,first:last) = ddot_S_I;
end

% Berechnung letzter Punkt der Trajektorie
S(:,last+1) = W_stuetz(:,N_I+1);
dot_S(:,last+1) = 0;
ddot_S(:,last+1) = 0;

% Bercehnung Zeitvektor
T = 0:delta_T:N_I*t_I(N_T_I);   % Anzahl an Intervallen multipliziert mit der verwendeten Zeit eines Intervalls
%% --- ENDE ARBEITSBEREICH --------------------------------------------
end % function

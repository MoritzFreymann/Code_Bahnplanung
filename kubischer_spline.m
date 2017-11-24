function [ S, dot_S, ddot_S, T ] = kubischer_spline( W_stuetz, T_ges, delta_T )
% Erzeugt aus N_I Stuetzvektoren in W_stuetz eine Trajektorie in Form kubischer Splines
% S         := Trajektorie auf Positionsebene
% dot_S     := Trajektorie auf Geschwindigkeitsebene
% ddot_S    := Trajektorie auf Beschleunigungsebene
% T         := Zeitvektor der Trajektorie

% W_stuetz  := Stuetzpunkte
% T_ges     := Dauer der Bewegung/Interpolation
% delta_T   := Taktzeit

% Erzeugt aus n Stuetzvektoren p_i und einer Zeit T eine Trajektorie aus
% stueckweise zusammengesetzten Polynomen dritten Grades (Spline)


%% --- ARBEITSBEREICH: ------------------------------------------------
%% Dimensionen pruefen/festlegen
%     N_Q = ?;
%     N_I = ?;
%     N_T = ? ;
% Anzahl der Freiheitsgrade
N_Q       = size( W_stuetz,1 );

% Anzahl der Stuetzstellen
N_p = size( W_stuetz,2);
% Anzahl der Intervalle
N_I       = N_p -1;

% Zeitintervall fuer Interpolation
T_I       = T_ges/N_I;          % Zeit fuer ein Teilstueck
t_I       = 0:delta_T:T_I;      % Zeitvektor fuer ein Teilstueck
N_T_I     = length(t_I);        % Anzahl der Zeitpunkte eines Teilstuecks

%T_ges neu setzen, damit T_ges exakt bei N_T erreicht ist
T_ges = N_I*t_I(N_T_I);
%% --- ENDE ARBEITSBEREICH --------------------------------------------

% Aequidistanter h-Vektor
h = T_ges / ( N_I ) * ones( 1, N_I );

% Initialisierung von A und r
A = zeros(N_p);
r = zeros(N_p,N_Q);

%% Berechnung A
% erste Zeile
A(1,1) = 2*h(1);
A(1,2) = h(1);

spalte = 1; % setzte zu beschreibende Spalte auf 1
% Zeilen 2 - n-2
for i=2:N_p-1
    
    A(i,spalte) = h(i-1);               % schreibe erste Spalte
    spalte = spalte+1;                  % eine Spalte weiter
    A(i,spalte) = 2*( h(i) + h(i-1) );  % schreibe zweite Spalte
    spalte = spalte+1;                  % eine Spalte weiter
    A(i,spalte) = h(i);                 % schreibe dritte Spalte
    spalte = spalte-1;                  % eine Spalte zurueck (neue Spalte ist eins links von erster beschribener Spalte der vorherigen Zeile)
    
end
% letzte Zeile
A(N_p,N_p-1) = h(N_p-1);
A(N_p,N_p) = 2*h(N_p-1);

%% Berechnung r
% erste Zeile
r(1,:) = ( 6*( W_stuetz(:,2) - W_stuetz(:,1) )/h(1) )';
% Zeilen 2 - n-2
for i=2:N_p-1
    
    r(i,:) = ( -6*( W_stuetz(:,i) - W_stuetz(:,i-1) )/h(i-1) + 6*( W_stuetz(:,i+1) - W_stuetz(:,i) )/h(i) )';
    
end
% letzte Zeile
r(N_p,:) = ( -6*( W_stuetz(:,N_p) - W_stuetz(:,N_p-1) )/h(N_p-1) )';
%% Berechnung ddot_p
ddot_p = A\r;

%% --- ARBEITSBEREICH: ------------------------------------------------
% Berechne Zeitvektor T fuer die gesamte Zeitspanne
T = 0:delta_T:N_I*t_I(N_T_I); % Anzahl an Intervallen multipliziert mit der verwendeten Zeit eines Intervalls

% Initialisierung der Trajektorie
S = zeros(N_Q, length(T) );
dot_S = zeros( size(S) );
ddot_S = zeros( size(S) );

% Schleife ueber alle Intervalle
for i = 1:N_I
    
    % Berechne Koeffizienten
    a = ( ddot_p(i+1,:) - ddot_p(i,:) )' / ( 6*h(i) );
    b = ddot_p(i,:)' / 2.0;
    c = ( W_stuetz(:,i+1) - W_stuetz(:,i) ) / h(i) - h(i) * ( ddot_p(i+1,:) + 2*ddot_p(i,:) )' / 6;
    d = W_stuetz(:,i);
    
    % Spalte der Gesamt-Trajektorie, in die der erste Punkt  
    % des jeweiligen Intervalls i geschrieben wird.
    % Alle weiteren Punkte des Intervalls sind um k verschoben
    spalte = (i-1)*(N_T_I-1); 
    
    % Berechnung der Trajektorie fuer aktuelles Intervall
    % (letzter Punkt wird weggelassen, da dieser der erste des naechsten
    % Intervalls ist)
    for k=1:N_T_I-1
               
        % Schreibe Werte der Trajektorie vom ersten Punkt des Intervalls
        % bis zum n-1ten Punkt in Gesamt-Trajektorie
        S(:, spalte + k) = a * t_I(k)^3 + b * t_I(k)^2 + c * t_I(k) + d;
        dot_S(:,spalte + k) = 3*a * t_I(k)^2 + 2*b * t_I(k) + c;
        ddot_S(:,spalte + k) = 6*a * t_I(k) + 2*b;
    end
end

% Letzten Punkt der Trajektorie berechnen
S(:,N_I*(N_T_I-1)+1) = a * t_I(N_T_I)^3 + b * t_I(N_T_I)^2 + c * t_I(N_T_I) + d;
dot_S(:,N_I*(N_T_I-1)+1) = 3*a * t_I(N_T_I)^2 + 2*b * t_I(N_T_I) + c;
ddot_S(:,N_I*(N_T_I-1)+1) = 6*a * t_I(N_T_I) + 2*b;

%% --- ENDE ARBEITSBEREICH --------------------------------------------

%% Erzeuge Spline fuer jede Komponente

end


%% Hilfsfunktion

%% --- ARBEITSBEREICH: ------------------------------------------------
%function [ ?, ?, ?] = kubischer_spline_skalar( ? )
%% Variablen fuer Gleichungssystem A* ddot_p = r anlegen
%     A = zeros( ?);
%     r = zeros( ? );

% Eigentliche Spline Koeffizienten
%     a = zeros( ? );
%     b = zeros( ? );
%     c = zeros( ? );
%     d = zeros( ? );

%% Erstelle A-Matrix und r-Vektor


%% Gleichungssystem loesen


%% Koeffizientenberechnung
%     for i = 1:?
%         a(i) = ?;
%         b(i) = ?;
%         c(i) = ?;
%         d(i) = ?;
%     end

%% Spline generieren
%     s = zeros(?);
%     dot_s = zeros(?);
%     ddot_s = zeros(?);

% Schleife ueber alle Intervalle

%% --- ENDE ARBEITSBEREICH --------------------------------------------

%end




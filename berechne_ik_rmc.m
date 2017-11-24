function rob = berechne_ik_rmc(rob,W,K,Option)
    % Berechnung der Inversen Kinematik ueber Resolved Motion Rate Control
    % W     ...positiv definite Wichtungsmatrix
    % K     ...Driftkompensationsmatrix
    % Option...dem Verfahren uebergebene Optionen

    %% --- ARBEITSBEREICH: ------------------------------------------------
    % Ergaenzen Sie bitte hier die fehlenden Stellen
    % ---------------------------------------------------------------------
    % Berechnung der mit W gewichteten Pseudoinversen der Arbeitsraum-Jacobimatrix       
    Jw_pseudo = W\rob.Jw'/( rob.Jw/W * rob.Jw' );
    % J_w_# = W^-1 * J_w^T * (J_w * W^-1 * J_w^T)^-1 ...(Vortrag IK, Folie 21)
    % oder ...(Formel 2.4.57, S.28 Skript)
    
    % Berechnung der Gelenkwinkelgeschwindigkeit ueber Resolved Motion Rate Control
    if strcmp(Option,'drift') == true
        % Aufgabe 2.1 RMC - ohne Driftkompensation
        dot_q_new = Jw_pseudo * rob.dot_w_d;

    elseif strcmp(Option,'driftcomp') == true
        % Aufgabe 2.2 RMC - mit Driftkompensation
        
        % Berechnung 
        e = rob.w_d - rob.w;
        
        % Berechnung effektive Geschwindigkeit
        dot_w_eff = rob.dot_w + K * e;
        % ...(Formel 2.4.37 S.24 Skript)
        
        % Berchnung Gelenkwinkelgeschwindigkeit
        dot_q_new = Jw_pseudo * dot_w_eff;

    else
        % Ungueltige Option
        error('Ungueltige Option gewaehlt!')
    end
    %% --- ENDE ARBEITSBEREICH --------------------------------------------

    % Berechnung der Gelenkwinkelbeschleunigung aus Differenzenquotient
    rob.ddot_q = (dot_q_new-rob.dot_q)/rob.dt;

    % Uebernehmen der berechneten Gelenkwinkelgeschwindigkeit
    rob.dot_q = dot_q_new;

    % Gelenkwinkel ueber explizites Euler-Verfahren berechnen
    rob.q = rob.q+rob.dot_q*rob.dt;
    
end
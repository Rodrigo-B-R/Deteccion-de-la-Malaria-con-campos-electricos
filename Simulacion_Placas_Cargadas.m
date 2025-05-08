function[tray,  X, Y, Ex, Ey, x1, y1, x2, y2] = Simulacion_Placas_Cargadas(L_positiva, L_negativa, distancia, v_positivo, v_negativo, q_sangre, m_sangre, v0_sangre, r0_sangre, h)
    epsilon0 = 8.854e-12;
    N1=500; 
    N2=500;
    
    lambda1 = (epsilon0 * v_positivo) / distancia;
    lambda2 = (epsilon0 * v_negativo) / distancia;
   

    margen_x = 0.1 * distancia;
    margen_y = 0.1 * max(L_positiva, L_negativa);

    % Discretización de placas
    y1 = linspace(0, L_positiva, N1);
    x1 = zeros(1, N1);
    y2 = linspace(0, L_negativa, N2);
    x2 = distancia * ones(1, N2);

    % Mallado del campo
    Xmin = -margen_x; Xmax = distancia + margen_x;
    Ymin = -margen_y; Ymax = max(L_positiva, L_negativa) + margen_y;
    [X, Y] = meshgrid(linspace(Xmin, Xmax, 60), linspace(Ymin, Ymax, 60));
    Ex = zeros(size(X));
    Ey = zeros(size(Y));

    % Calcular campo eléctrico
    for i = 1:N1
        rx = X - x1(i); ry = Y - y1(i);
        r2 = rx.^2 + ry.^2; r = sqrt(r2);
        Ex = Ex + (lambda1 / (2*pi*epsilon0)) * rx ./ (r2 .* r);
        Ey = Ey + (lambda1 / (2*pi*epsilon0)) * ry ./ (r2 .* r);
    end
    for i = 1:N2
        rx = X - x2(i); ry = Y - y2(i);
        r2 = rx.^2 + ry.^2; r = sqrt(r2);
        Ex = Ex + (lambda2 / (2*pi*epsilon0)) * rx ./ (r2 .* r);
        Ey = Ey + (lambda2 / (2*pi*epsilon0)) * ry ./ (r2 .* r);
    end

    r = r0_sangre(:);  
    v = v0_sangre(:);  
    tray = r;
    
    while true
        % Interpolación del campo eléctrico en la posición actual
        Ex_interp = interp2(X, Y, Ex, r(1), r(2), 'linear', 0);
        Ey_interp = interp2(X, Y, Ey, r(1), r(2), 'linear', 0);
        E = [Ex_interp; Ey_interp];  % Campo eléctrico en la posición actual

        % Calcular la fuerza y la aceleración
        F = q_sangre * E;   % fuerza sobre la gota
        a = F / m_sangre;   % aceleración (F = m * a)

        % Actualización de la velocidad y la posición con el método de Euler
        v = v + h * a;      % nueva velocidad
        r = r + h * v;      % nueva posición

        % Guardar la posición en la trayectoria
        tray(:, end+1) = r;

        % Verificar si la gota ha salido de los límites espaciales
        if r(1) < Xmin || r(1) > Xmax || r(2) < Ymin || r(2) > Ymax
            break;  % si se sale del área, detener la simulación
        end
    end

end

    



    





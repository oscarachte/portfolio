% Miniprojet 1 - Equation de diffusion - Methode de Chebyshev
% Oscar Achte, Justin Rozzio, Jules Van Cauwenberghe, Titouan Vaganay

clear; close all; clc;

% --- Parametres ---
alpha_values = [2, 3, 4, 5];
npoints_list = [8, 16, 32, 64, 128];
colors       = lines(length(alpha_values));

errors = zeros(length(alpha_values), length(npoints_list));
hvals  = zeros(1, length(npoints_list));

% --- Boucle principale ---
for ia = 1:length(alpha_values)
    alpha = alpha_values(ia);
    fprintf('--- alpha = %d ---\n', alpha);

    for ip = 1:length(npoints_list)
        npoints = npoints_list(ip);

        % Noeuds de Chebyshev-Gauss-Lobatto sur [0,1]
        xp        = ChebyshevNodes(0, 1, npoints);
        hvals(ip) = 1 / (npoints - 1);

        % Matrices de base et derivees secondes
        [PHI, ~, PHIxx] = ChebyshevInterpolation(npoints, xp);

        % Assemblage du systeme lineaire avec conditions de Dirichlet
        f   = sourceterm(xp, alpha);
        A   = PHIxx;
        rhs = f;
        A(1,:)   = PHI(1,:);   rhs(1)   = 0;           % u(0) = 0
        A(end,:) = PHI(end,:); rhs(end) = (0.5)^alpha;  % u(1) = (1/2)^alpha

        % Resolution : coefficients spectraux
        a = A \ rhs;

        % Erreur L2 sur grille fine (evite le biais aux noeuds CGL)
        xfine            = linspace(0, 1, 1000)';
        [PHI_fine, ~, ~] = ChebyshevInterpolation(npoints, xfine);
        u_h              = PHI_fine * a;
        u_ex             = exactsolution(xfine, alpha);
        dx               = xfine(2) - xfine(1);
        errors(ia, ip)   = sqrt(dx * sum((u_h - u_ex).^2));

        fprintf('  npoints = %3d  ||u - u_h||_L2 = %.4e\n', npoints, errors(ia, ip));
    end
    fprintf('\n');
end

% --- Figure 1 : Convergence L2 (log-log) ---
figure('Name','Convergence L2','Position',[100 100 750 550]);
loglog(hvals, errors(1,:), '-o', 'Color', colors(1,:), ...
       'LineWidth', 1.8, 'MarkerSize', 7, ...
       'DisplayName', sprintf('\\alpha = %d', alpha_values(1)));
hold on; grid on; box on;
for ia = 2:length(alpha_values)
    loglog(hvals, errors(ia,:), '-o', 'Color', colors(ia,:), ...
           'LineWidth', 1.8, 'MarkerSize', 7, ...
           'DisplayName', sprintf('\\alpha = %d', alpha_values(ia)));
end
% Pentes de reference theoriques : ordre r = alpha-1
for ia = 1:length(alpha_values)
    r     = alpha_values(ia) - 1;
    ref_y = errors(ia,1) * (hvals / hvals(1)).^r;
    loglog(hvals, ref_y, '--', 'Color', colors(ia,:)*0.7, ...
           'LineWidth', 1.2, 'DisplayName', sprintf('Pente ref. %d', r));
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
xlabel('h = 1/(p-1)', 'FontSize', 13);
ylabel('||u - u_h||_{L_2}', 'FontSize', 13);
title('Convergence en norme L_2 - Colocation Chebyshev', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 11);
saveas(gcf, 'convergence_L2.png');

% --- Figure 2 : Solution exacte vs approchee pour p=64 ---
figure('Name','Solutions','Position',[100 100 1100 400]);
npoints_plot = 64;
xfine_plot   = linspace(0, 1, 500)';

for ia = 1:length(alpha_values)
    alpha = alpha_values(ia);
    subplot(1, 4, ia);
    hold on; grid on; box on;

    xp = ChebyshevNodes(0, 1, npoints_plot);
    [PHI, ~, PHIxx] = ChebyshevInterpolation(npoints_plot, xp);
    f   = sourceterm(xp, alpha);
    A   = PHIxx; rhs = f;
    A(1,:)   = PHI(1,:);   rhs(1)   = 0;
    A(end,:) = PHI(end,:); rhs(end) = (0.5)^alpha;
    a   = A \ rhs;

    [PHI_p, ~, ~] = ChebyshevInterpolation(npoints_plot, xfine_plot);
    plot(xfine_plot, exactsolution(xfine_plot, alpha), 'k-',  'LineWidth', 2,   'DisplayName', 'Exacte');
    plot(xfine_plot, PHI_p * a,                        'r--', 'LineWidth', 1.8, 'DisplayName', 'Chebyshev');
    xlabel('x', 'FontSize', 11); ylabel('u(x)', 'FontSize', 11);
    title(sprintf('\\alpha = %d, p = %d', alpha, npoints_plot), 'FontSize', 12);
    legend('FontSize', 9, 'Location', 'northwest');
end
sgtitle('Solution exacte vs approchee', 'FontSize', 13);
saveas(gcf, 'solutions_comparison.png');

% --- Ordres de convergence numeriques ---
fprintf('%-8s', 'alpha');
for ip = 2:length(npoints_list)
    fprintf('  p=%d->%d', npoints_list(ip-1), npoints_list(ip));
end
fprintf('\n');
for ia = 1:length(alpha_values)
    fprintf('%-8d', alpha_values(ia));
    for ip = 2:length(npoints_list)
        order = log(errors(ia,ip-1)/errors(ia,ip)) / log(hvals(ip-1)/hvals(ip));
        fprintf('  %7.3f ', order);
    end
    fprintf('\n');
end

% ================================================================
%  FONCTIONS LOCALES
% ================================================================

function f = sourceterm(x, alpha)
% f = u'' : nul sur [0, 1/2], alpha*(alpha-1)*(x-1/2)^(alpha-2) sur (1/2, 1]
    f = zeros(size(x));
    idx = x > 0.5;
    if alpha >= 2
        f(idx) = alpha * (alpha-1) * (x(idx) - 0.5).^(alpha-2);
    end
end

function u = exactsolution(x, alpha)
% u = 0 sur [0, 1/2], (x-1/2)^alpha sur (1/2, 1]
    u = zeros(size(x));
    idx = x > 0.5;
    u(idx) = (x(idx) - 0.5).^alpha;
end
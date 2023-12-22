%%%%%%%%%% LOAD FILES AND SET UP VARIABLES %%%%%%%%%%

stoichiometry_K = readtable('stoichiometry_K.csv');
lbl = stoichiometry_K.Properties.VariableNames;

stoichiometry_K = table2array(stoichiometry_K);
K = stoichiometry_K(:,end);
stoichiometry = stoichiometry_K(:,1:end-1);

mass_conservation_total_masses = readmatrix('mass_conservation_total_masses.csv');
total_masses = mass_conservation_total_masses(:,end);
mass_conservation = mass_conservation_total_masses(:,1:end-1);

%%%%%%%%%% RUN SOLVER %%%%%%%%%%

tic
for i = 1:100
[x, n, delta, error_flag] = equpy(stoichiometry, K, mass_conservation, total_masses, [], 20, 0);
end
time = toc/i*1000;

%%%%%%%%%% PRINT RESULTS AND CREATE FIGURE %%%%%%%%%%

disp(sprintf('execution time --- %f milliseconds ---', time))

close(figure(1))

figure(1)
subplot(1,2,1)
scatter(1:n, delta, 120, 'black', 'filled')
hold on
plot(1:n, delta, 'linewidth', 3, 'color', 'black')
xlabel('steps')
ylabel('||r||')
set(gca, 'box', 'on', 'yscale', 'linear', 'linewidth', 2, 'fontsize', 14)
xlim([1 n])

subplot(1,2,2)
bar(x)
xlabel('species')
ylabel('conc (a.u.)')
xticks(1:length(x))
xticklabels(lbl)
set(gca, 'box', 'on', 'linewidth', 2, 'fontsize', 14)

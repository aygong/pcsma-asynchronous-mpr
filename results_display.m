function results_display(T_ana, T_sim)
% Display the results and the corresponding error
fprintf("-> Analysis = %.4f, Simulation = %.4f\n", T_ana, T_sim);
fprintf("-> Error = %.4f%%\n", abs(T_ana - T_sim) / T_ana * 100);
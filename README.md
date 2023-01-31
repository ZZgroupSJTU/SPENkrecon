# SPENkrecon

simulation case:
    run intR_case.m or nonintR_case.m first (for approximate result);
    and then run optprocess.m (for optimized result);

    in intR_case.m or nonintR_case.m; line 8: RES = 8192;
    its value is better to be bigger, for which the simulation is much closer to the real situation;

    in intR_case.m, set value R(2,3,4...) in line 11 for different cases;
    in nonintR_case.m, we used Numerator and Denumerator to generate nonint acceleration rate;

acq case:
    run run_demo.m

    we provide different objects' cases: ACR, Brain, and Orange;
    in section "load and set" select these included in 'Data' file

    the "power" and "weight" can decisively influence the performance of the optimization result.

    for ACR and Brain, we attached the zoomed in compareing setting in the end. If you want to use it, please keep consistent with the initial data selection


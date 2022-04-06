%% run all methods

parpool('local',13)

clear;
clc;
% Model 1
baseline_fun("model1",2,0.08,"yes","yes","no");
fprintf("Finshed Model 1")
clear;
clc;
% Model 2
baseline_fun("model2",2,0.04,"yes","yes","no");
fprintf("Finshed Model 2")
clear;
clc;
% Model 3
baseline_fun("model3",2,0.02,"yes","yes","no");
fprintf("Finshed Model 3")
clear;
clc;
% Model 4
baseline_fun("model4",2,0.01,"yes","yes","no");
fprintf("Finshed Model 4")
clear;
clc;
% Model 5
baseline_fun("model5",2,0.02,"no","yes","no");
fprintf("Finshed Model 5")
clear;
clc;
% Model 6
baseline_fun("model6",2,0.02,"no","no","no");
fprintf("Finshed Model 6")
clear;
clc;
%Model 7
baseline_fun("model7",1,0.02,"no","yes","no");
fprintf("Finshed Model 7")
clear;
clc;
% Model 8
baseline_fun("model8",1,0.02,"yes","yes","no");
fprintf("Finshed Model 8")
clear;
clc;
% Model 9                                                                       
baseline_fun("model9",3,0.02,"no","yes","no");                                 
fprintf("Finshed Model 9")                                                      
clear;                                                                          
clc; 
% Model 10                                                                      
baseline_fun("model10",3,0.02,"yes","yes","no");                                 
fprintf("Finshed Model 10")                                                      
clear;                                                                          
clc; 
% Model 11                                                                      
baseline_fun("model11",2,0.02,"no","no","yes");                                 
fprintf("Finshed Model 11")                                                      
clear;                                                                          
clc; 



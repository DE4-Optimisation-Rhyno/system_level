%--------------------Subsystem_1 Robot Chassis-----------------------------
%--------------------load in data for Mg Alloy-----------------------------
Sim_data = csvread('Formatted_Data_Norms.csv',1,1,[1,1,40,12]);
Design_variables = (Sim_data(21:30,[1,2,3])); 
Mg_outputs = (Sim_data(21:30,[4,6]));
Mg_data = [Design_variables,Mg_outputs];

%----Calculate Betas Values for mass and FS--------------------------------

Betas_mass = mvregress(Design_variables,Mg_outputs(:,1)); %calculating Beta values with mvregress
Betas_fos = mvregress(Design_variables,Mg_outputs(:,2));

%multi-objective function that combines everything 
fun_1 = @(x1) 0.8*(Betas_mass(1)*x1(1) + Betas_mass(2)*x1(2) + Betas_mass(3)*x1(3)) - 0.2*(Betas_fos(1)*x1(1) + Betas_fos(2)*x1(2) + Betas_fos(3)*x1(3));

%-------------------Subsystem_1 optimisation-----------------------
IG = [max(Design_variables)];
Lb = [300 253.2 1.6]; %applying normalisation to upper and lower bound values 
Ub = [1000 660 30];
A = [0 ,-1, 2];
b = [-250];
[Opti_dims,Opti_vals] = fmincon(fun_1,IG,A,b,Aeq,beq,Lb,Ub,nonlcon,options1);

%--------Using linear regression to find mass again---------------

mass_subsystem_1 = Opti_dims(1)*Betas_mass(1) + Opti_dims(2)*Betas_mass(2) + Opti_dims(3)*Betas_mass(3); %mass of the subsystem

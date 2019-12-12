%% Optimising Design of a Firefighting Robot - System Level Optimisation
clc
clear
tic
%% --------------------Subsystem_1 Robot Chassis-----------------------------
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
Aeq = [];
beq = [];
[Opti_dims,Opti_vals] = fmincon(fun_1,IG,A,b,Aeq,beq,Lb,Ub);

%--------Using linear regression to find mass again---------------

mass_subsystem_1 = Opti_dims(1)*Betas_mass(1) + Opti_dims(2)*Betas_mass(2) + Opti_dims(3)*Betas_mass(3) %mass of the subsystem

%% Subsystem 2 - Propulsion system

%Define lower & upper-bounds
lb = [11,0.6,1];
ub = [300,0.8,3.09];

Ap = [];
bp = [];

%Initial Starting Point, [N2,r]
x0 = [11,0.8,1];

P=500.00;
mu_p=0.71;
Mc=160.000;
g=9.81;
w=[1,0];

%Initial Guess Acceleration + Motor Mass values 
disp(['Initial Acceleration without consideration of constraints: ' num2str(accelws(x0,P,mu_p,Mc,g,w))])
w=[0,1];
disp(['Initial Motor Mass without consideration of constraints: ' num2str(accelws(x0,P,mu_p,Mc,g,w))])

% Optimize for Acceleration
disp('Acceleration Optimized')

%SQP
disp('SQP')
w=[1,0];

options = optimset('Algorithm','sqp');
[x] = fmincon(@(x0) accelws(x0,P,mu_p,Mc,g,w),x0,Ap,bp,Aeq,beq,lb,ub,@(x) conws(x),options);

optAccel=accelws(x,P,mu_p,Mc,g,w);
w=[0,1];
optMotorMass=accelws(x,P,mu_p,Mc,g,w);

%Interior Points
disp('Interior Points')
w=[1,0];
tic
options = optimset('Algorithm','interior-point');
x = fmincon(@(x0) accelws(x0,P,mu_p,Mc,g,w),x0,A,b,Aeq,beq,lb,ub,@(x) conws(x),options);
toc
optAccel=accelws(x,P,mu_p,Mc,g,w);
w=[0,1];
mass_subsystem_2=accelws(x,P,mu_p,Mc,g,w)

%% Subsystem 3 - Battery system

load('BatterySet0812.mat');

PackAh = table2array(BatterySet1112(:,19));
PackCost = table2array(BatterySet1112(:,18));

BattAttainableAh = find(abs(PackAh)>40);
x(BattAttainableAh) = [];

BattAttainablePower = find(abs(MaxPackPower)>250);
x(BattAttainablePower) = [];


[val]=intersect(BattAttainableAh,BattAttainablePower);

AttainableTable = BatterySet1112(val',:);

AttainablePackAh = table2array(AttainableTable(:,19));
AttainablePackCost = table2array(AttainableTable(:,18));
AttainablePackMass = table2array(AttainableTable(:,16));

AttainableBatteryNames = table2array(AttainableTable(:,1));

n2 = AttainableBatteryNames;
allOneString2 = sprintf('%s,' , n2);
allOneString2 = allOneString2(1:end-2); %strip final comma
splitString2 = split(allOneString2,',');
TransposedSplitString2 = splitString2';

%Self Adjudged Pareto Set (Mass vs Cost) - for the Attainable Set

%Figure 1: Max Pack Power vs Pack Mass
figure('Name', 'Figure 7: Attainable Set Mass (Kg) vs Attainable Set Cost (GBP)') %Naming the figure

xpos = [AttainablePackCost];
ypos = [AttainablePackMass];
xlim([130 260]);
ylim([0.5 1.5]);
labels = TransposedSplitString2;
%h = labelpoints(xpos, ypos, labels, 'N', 0.2, 1); 
hold on
sz = 25;
scatter(AttainablePackCost, AttainablePackMass, sz,'filled');

xlabel('Pack Cost (GBP)') %Labelling X-Axis
ylabel('Pack Mass (Kg))') %Labelling Y-Axis
title('Figure 7: Attainable Set Mass (Kg) vs Attainable Set Cost (GBP)') %adds a title to the graph
hold off

% Printing the Values of the Selected Battery 

SelectedBatteryPack = AttainableTable(2,:)
SelectedBatteryMass = table2array(AttainableTable(2,16))

mass_subsystem_3 = 2*SelectedBatteryMass

%% Subsystem 4 - water delivery system
% Parameters
g = 9.81;
mu = 0.00089;
rho = 1000;
H = 12.5; %Average set value set via parametric analysis
ep = 0.0001525;
% Initial Values
x0 = [0.2131,45,750]; %taking the average value informed by lower and upper bounds
% Lower and Upper Variable Bounds from Constraints
lb = [0.0762,10,500];
ub = [0.35,80,1000];
% Linear Constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear Constraints
nonlincon = @wt_nlcon;

% Objective Function in negative null form
objective = @(x) -sin(2.*x(2)).*-mu.*(H-20.*x(1).*sin(x(2))-((x(3).*g)./(2.*pi.*x(1)))-((1./(-1.8.*log10((ep./(7.4.*x(1))).^1.1 + ((241.*x(1))./(25000.*x(3)))))).^2.*mu.*20.*x(1).*x(3).^2)./(pi.^2.*x(1).*3.*g));
rng default
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',x0,'objective',objective,'lb',lb,'ub',ub,'nonlcon',nonlincon)
x_gs = run(gs,problem)

r = x_gs(1); % optimised r value
l = r * 20;
a = x_gs(2); % optimised a value

% Plotting Pareto Front
% Import discrete data for dependent variables from water pump specs
pump_shortlist = csvread('pump_values_shortlist2.csv'); % pump mass, flow rate, pump head
[row,col] = size(pump_shortlist);

% Declare H and f from the shortlist
Fsl = pump_shortlist(:,2);
Hsl = pump_shortlist(:,3);

% Create arrays for storing T and m values
T_calc = zeros(row,1);
m_calc = zeros(row,1);

% For loop to sub in different pump variables to find T and m
for px = 1:row
    T_calc(px) = max_T(ep,mu,g,r,l,a,Fsl(px),Hsl(px));
    m_calc(px) = min_m(Fsl(px));
end

% Finding the pareto front for T against m
plot(-m_calc,T_calc,'o')
hold on
title('Trade-off between spray travel anAd mass using different pumps')
xlabel('Mass (kg)')
ylabel('Spray travel (m)')
set(gcf,'color','w');
hold off

mass_subsystem_4 = 26

%% System Level Optimisation

% Adding masses from all subsystems together
mass_total = mass_subsystem_1 + mass_subsystem_2 + mass_subsystem_3 + mass_subsystem_4
disp(['Final Optimised Robot Mass (kg) = ' num2str(mass_total)]);
toc
%% Subsystem 4 Objective and Non-linear Constraint Functions

%Calculating T - vars r,l,a,F,H
function T = max_T(ep,mu,g,r,l,a,F,H)
X = (ep./(7.4.*r)).^1.1 + ((241.*r)./(25000.*F));
Hf =((1./(-1.8.*log10(X))).^2.*mu.*l.*F.^2)./(pi.^2.*r.*3.*g);
T = sin(2.*a).*-mu.*(H-l.*sin(a)-((F.*g)./(2.*pi.*r))-Hf);
end

%Calculating m - var F
function m = min_m(F)
F_s = F/100;
m = 0.0087.*F_s.^3 -0.1663.*F_s.^2 + 1.9751.*F_s + 24.908;
end

% Non-linear Constraints for water delivery system 
function [c,ceq] = wt_nlcon(x)
    c = [Constraint_g1(x),Constraint_g6(x),Constraint_g11(x)];
    ceq = 0;
end

% Constraint g1: Reynold's number cannot exceed laminar flow
function g1 = Constraint_g1(x)
    % Paramters
    mu = 0.00089;
    rho_water = 1000;
    
    % Variables
    r = x(1);
    F = x(3);
    
    % Constraint equation
    g1 = ((2.*rho_water.*F.*mu.^2)./(pi.*r))-4000;

end

% Constraint g6: Jet reaction force cannot exceed 809N
function g6 = Constraint_g6(x)
    % Paramters
    mu = 0.00089;
    
    % Variables
    r = x(1);
    a = x(2);
    F = x(3);
    
    % Constraint equation
    g6 = ((F.*mu.*cos(a))/(pi.*r.^2))-809;

end

% Constraint g11: Keeping full pipe via ratio between radius and flow rate
function g11 = Constraint_g11(x)
    % Paramters
    rho_water = 1000;

    % Variables
    r = x(1);
    F = x(3);
    
    % Constraint equation
    g11 = (F/r.^0.25)-1.5.*rho_water;

end


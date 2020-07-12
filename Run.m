%%% Script for runing the whole design approach
%%% Identification of parameters of Glucose-Insulin model using both real
%%% and simulated data
% Assumption: we assume the measurement errors are white noise 
% Parameters are [tmax, tp, Ib, ka1, ka2, ka3, sT, sD, sE, tmmax, EGP0, F01, k12, ks, Ip0]

% load data
currentFolder = pwd;
% load(string(currentFolder) + '\Data\SimulationData\DataF.mat')
load(string(currentFolder) + '\Data\ExperimentalData\DataExp.mat')


%removing NAN
NanPos = find(isnan(data.xdata));
if(isempty(NanPos))
    Loc = length(data.ydata);
else
    Loc = NanPos(1)-1;
end
data.tdata = data.tdata(1:Loc);
data.ydata = data.ydata(1:Loc);
data.xdata = data.xdata(1:Loc);


%%% Constant parameters
global w VI V MealTimes MealCHOs BolusIns InitialG

w = 75; % Patient's weight (kg)
MealTimes = [1*10 31*10] ; % Meal time samples
MealCHOs = [25 45]; % Meal carb samples (grams)
BolusIns = [19.8 18.9]; % Bolus values (unit)
InitialG = 11; % initial value of sensor glucose
VI = 120; % Insulin distribution volume (mL/kg)
V = 150; % Glucose distribution volume (mL/kg)
matlab.addons.toolbox.installToolbox('mcmc toolbox.mltbx');


% initial guess for the parameters
Param0 = [55, 9, 0.1, 0.0034, 0.056, 0.024, 0.001841, 0.0005, 0.019, 45, 16.9, 11.1, 0.06, 0.093, 20];


%% Bayesian estimation (MCMC)
% sum squared function
model.ssfun = @ModelSS;

% parameters 
params = {
    {'tmax',    Param0(1),  0, 1000, 55, 50}
    {'tp',      Param0(2),  0, 1000, 9,  10}
    {'Ib',      Param0(3),  0, 100,  0.1, 0.4}
    {'ka1',     Param0(4),  0, 0.1 , 0.0034, 0.0034}
    {'ka2',     Param0(5),  0, 5 , 0.056, 0.056}
    {'ka3',     Param0(6),  0, 2 , 0.024, 0.024}
    {'sT',      Param0(7),  0, 0.1 , 0.001841, 0.001841}
    {'sD',      Param0(8),  0, 0.05 , 0.0005, 0.0005}
    {'sE',      Param0(9),  0, 2 , 0.019, 0.019}
    {'tmmax',   Param0(10), 0, 4500 , 45, 45}
    {'EGP0',    Param0(11), 0, 1700 , 16.9, 16.9}
    {'F01',     Param0(12), 0, 1110 , 11.1, 11.1}
    {'k12',     Param0(13), 0, 6 , 0.06, 0.06}
    {'ks',      Param0(14), 0, 10 , 0.093, 0.093}
    {'Ip0',     Param0(15), 0, 2000 , 20, 20}
    };

%prior information (for residual or noises)
model.S20 = [1];
model.N0  = [4];

% generate an initial "burn-in" chain
options.nsimu = 10000;
options.verbosity = 0;
[results, chain, s2chain]= mcmcrun(model,data,params,options);

% Now re-run starting from the results of the previous run,
% this will take about 20 minutes.
options.nsimu = 50000;
options.verbosity = 0;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);


% % predictive model function
% modelfun = @ModelPred;
modelfun = @(d,th) ModelPred(data.tdata,th,data);

% sample N parameters and calculate predictive plots
nsample = 500;
out = mcmcpred(results,chain,s2chain,data.tdata,modelfun,nsample);

% plot results
dimc = [0.87 0.87 0.87]; % dimmest (lightest) color
lowerBound =[];
medianBound = [];
upperBound = [];
for i=1:length(out.obslims{1,1})
    AA = out.obslims{1,1}{1, i};
    lowerBound(i) = AA(1);
    medianBound(i) = AA(2);
    upperBound(i) = AA(3);
end

medD=[];
for i=1:length(out.predlims{1,1})
    AA = out.predlims{1,1}{1, i};
    medD(i) = AA(2);
end


if(results.nsimu > (options.nsimu)/2) 
    figure
    mcmcplot(chain,[],results,'denspanel',2);
    time=data.tdata;
    figure
    fillyy(time(1:end-1),lowerBound,upperBound,dimc);
    hold on
    plot(time(1:end-1), medianBound)
    plot(time, data.ydata,'*')
end


% figure
% plot(time(1:end-1), medD)
% hold on
% plot(time, data.ydata,'*')
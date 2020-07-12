% The function for calculating the error vector
function Err = ModelEV(params,data)

global w VI V MealTimes MealCHOs BolusIns InitialG

tdata = data.tdata;
ydata = data.ydata;
xdata = data.xdata;

Ip0 = params(15); % initial value of plasma insulin
% calculating other inital values based on initial glucose and initial plasma insulin


% initial values for the ODE
Gs0 = InitialG; 
Q10 = V * InitialG;
Qm10 = 0;
Qm20 = 0;
x10 =  params(7) * Ip0;
x20 =  params(8) * Ip0;
x30 =  params(9) * Ip0;
Q20 = Q10 * x10 /(x20 + params(13));
Qi0 = w*VI*Ip0/10^6;
Qi20 = (params(1)*Qi0/params(2))-params(1)*params(3)*w*VI/10^6;
Qi10 = Qi20;

StatesInit =  [Qi10, Qi20, Qi0, x10, x20, x30, Qm10, Qm20, Q10, Q20, Gs0];

Err = [];

for i=1:length(tdata)-1
    startTime = tdata(i);
    endTime = tdata(i+1);
    % Adding meal carb & bolus
    for j=1:length(MealTimes)
       if(endTime > MealTimes(j) && MealTimes(j) >= startTime)  
           % Adding meal carb (acts like a dirac function in the meal absorption dynamic) 
           StatesInit(7) = StatesInit(7) + (MealCHOs(j) * 5551)/(w);
           % Adding meal bolus (bolus input acts like a dirac function in the insulin absorption dynamic as u(t)) 
           StatesInit(1) = StatesInit(1) + BolusIns(j);
       end
    end
    
    [~,Y] = ode45(@(t,y) Model(y, xdata(i),params), [startTime endTime], StatesInit);

    StatesInit = Y(end,:);
        
    Err = [Err (Y(end,end)-ydata(i))];
end
end
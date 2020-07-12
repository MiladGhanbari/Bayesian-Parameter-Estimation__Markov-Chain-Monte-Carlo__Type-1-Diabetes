%% Glucose-Insulin model (different dynamics implementation)
function StatesDifferential = Model(y, input, Eparam)
    %%% Constant parameters
    global w VI V
    
    %%% Estimated parameters
    tmax = Eparam(1);
    tp = Eparam(2);
    Ib = Eparam(3);
    ka1 = Eparam(4);
    ka2 = Eparam(5);
    ka3 = Eparam(6);
    sT = Eparam(7);
    sD = Eparam(8);
    sE = Eparam(9);
    tmmax = Eparam(10);
    EGP0 = Eparam(11);
    F01 = Eparam(12);
    k12 = Eparam(13);
    ks = Eparam(14);
    
    
   
    
   %%% states vector is assumed as: y = [Qi1, Qi2, Qi, x1, x2, x3, Qm1, Qm2, Q1, Q2, Gs]
   % assigning each state to the corresponding variable
   Qi1 = y(1);  % insulin mass in the first subcutaneous compartments (units)
   Qi2 = y(2);  % insulin mass in the second subcutaneous compartments (units)
   Qi = y(3);   % insulin mass in the plasma (U)
   x1 = y(4);  % the delayed (remote) effects of insulin on glucose distribution (1/min)
   x2 = y(5);  % the delayed (remote) effects of insulin on glucose disposal (1/min)
   x3 = y(6);  % the delayed (remote) effects of insulin on the endogenous glucose production (1/min)
   Qm1 = y(7);  % glucose mass in the first gut compartments (umol/kg)
   Qm2 = y(8);  % glucose mass in the second gut compartments (umol/kg)
   Q1 = y(9);  % glucose mass in the accessible compartments (umol/kg)
   Q2 = y(10);  % glucose mass in the nonaccessible compartments (umol/kg)
   Gs = y(11);   % interstitial glucose concentration (mmol/L)

   % output (differential) vector initialization
   StatesDifferential = zeros(size(y));

   % Insulin absorption dynamics
   dQi1dt = input/60 - Qi1/tmax;
   dQi2dt = (Qi1 - Qi2)/tmax;
   dQidt = - Qi/tp + Qi2/tmax + Ib*w*VI/10^6;

   StatesDifferential(1) = dQi1dt;
   StatesDifferential(2) = dQi2dt;
   StatesDifferential(3) = dQidt;
   
   Ip = (Qi * 10^6)/(w * VI);

   % Insulin action dynamics
   dx1dt = - ka1 * x1 + ka1 * sT  * Ip;
   dx2dt = - ka2 * x2 + ka2 * sD  * Ip;
   dx3dt = - ka3 * x3 + ka3 * sE  * Ip;

   StatesDifferential(4) = dx1dt;
   StatesDifferential(5) = dx2dt;
   StatesDifferential(6) = dx3dt;


   % Meal absorption dynamic 
   dQm1dt = - Qm1/tmmax;
   dQm2dt =  (Qm1 - Qm2)/tmmax;
   Um = Qm2/tmmax;
   
   StatesDifferential(7) = dQm1dt;
   StatesDifferential(8) = dQm2dt;
   
   % Plasma Glucose dynamics
   EGP = EGP0 * (1 - x3);

   if(EGP<0)
       EGP=0;
   end

   dQ1dt = - ((F01*Q1/160)/(1+(Q1/160))) - x1 * Q1 + k12 * Q2 + EGP + Um;
   dQ2dt = x1 * Q1 - (k12 +  x2) * Q2;


   StatesDifferential(9) = dQ1dt;
   StatesDifferential(10) = dQ2dt;
   
   G = Q1/V;

   % Sensor (Interstitial Glucose) dynamics
   dGsdt = ks * (G - Gs);

   StatesDifferential(11) = dGsdt;

end
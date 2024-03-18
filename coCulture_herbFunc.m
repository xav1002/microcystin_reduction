function [t,y] = coCulture_herbFunc(tIn)
    % Initial Conditions
    global y0;
    % iterating over BatchFunction with ode45
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,y] = ode45(@BatchFunction,tIn,y0,opts);

    function dydt = BatchFunction(~,y)
        % 1. dXa/dt
        % 2. dXb/dt
        % 3. dN/dt
        % 4. dP/dt
        % 5. dCO2/dt
        % 6. dMCin/dt
        % 7. dH/dt
        % I

        % constant parameters
        global muMaxA muMaxB KNA KNB KPA KPB KCA KCB KIA KIB kdA kdB YNA YNB YPA YPB ktC Io A SoCO2 YCO2A YCO2B YT KHA KHB YH herbicide;
    
        % Light
        I = (Io/(A*(y(1)+y(2))))*(1-exp(-A*(y(1)+y(2))));

        % monod{Species}{Nutrient}
        monodAN = (y(3)/(y(3)+KNA));
        monodAP = (y(4)/(y(4)+KPA));
        monodAC = (y(5)/(y(5)+KCA));
        monodAI = (I/(I+KIA));
        monodInhibAH = (1/(1+y(7)/KHA));
        monodBN = (y(3)/(y(3)+KNB));
        monodBP = (y(4)/(y(4)+KPB));
        monodBC = (y(5)/(y(5)+KCB));
        monodBI = (I/(I+KIB));
        monodInhibBH = (1/(1+y(7)/KHB));

        n = 1;
        alpha = (y(3)/y(4))^n;
    
        dydt(1) = y(1)*(muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH-kdA);
        dydt(2) = y(2)*(muMaxB*monodBN*monodBP*monodBC*monodBI*monodInhibBH-(kdB));
        dydt(3) = (-YNA*y(1)*(muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH)) + ...
            (-YNB*y(2)*(muMaxB*monodBN*monodBP*monodBC*monodBI*monodInhibBH));
        dydt(4) = (-YPA*y(1)*(muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH)) + ...
            (-YPB*y(2)*(muMaxB*monodBN*monodBP*monodBC*monodBI*monodInhibBH));
        dydt(5) = ktC*(SoCO2-y(5)) + (-YCO2A*y(1)*(muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH)) + ...
            (-YCO2B*y(2)*(muMaxB*monodBN*monodBP*monodBC*monodBI*monodInhibBH));
        if herbicide
            dydt(6) = (YT/alpha)*(y(1)*muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH);
        else
            dydt(6) = (YT*alpha)*(y(1)*muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH);
        end
        dydt(7) = (-YH*y(7))*(y(1)*muMaxA*monodAN*monodAP*monodAC*monodAI*monodInhibAH + ...
            y(2)*muMaxB*monodBN*monodBP*monodBC*monodBI*monodInhibBH);
    
        dydt = dydt';
    end
end
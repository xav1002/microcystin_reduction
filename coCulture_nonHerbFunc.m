% Defining ParamFunc and BatchFunc
% Issues:
% 1. Bad starting parameter values;
% 2. Bad modeling relationships?
% 3. Syntax err
% 4. First smooth data with spline/curve fitting?
function [t,y] = coCulture_nonHerbFunc(tIn)
    % Initial Conditions
    global y0;
    yAdjust = y0(1:6);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,y] = ode45(@BatchFunction,tIn,yAdjust,opts);

    function dydt = BatchFunction(~,y)
        % 1. dXa/dt
        % 2. dXb/dt
        % 3. dN/dt
        % 4. dP/dt
        % 5. dCO2/dt
        % 6. dMCin/dt
        % I

        % constant parameters
        global muMaxA muMaxB KNA KNB KPA KPB KCA KCB KIA KIB kdA kdB YNA YNB YPA YPB ktC Io A SoCO2 YCO2A YCO2B YT kdL;
        
        % variable parameters
    
        % Light
        I = (Io/(A*(y(1)+y(2))))*(1-exp(-A*(y(1)+y(2))));
        % Add co-culture effect for nutrient equations
    
        dydt(1) = muMaxA*y(3)/(y(3)+KNA)*y(4)/(y(4)+KPA)*y(5)/(y(5)+KCA)*I/(I+KIA)*y(1)-kdA*y(1);
        dydt(2) = muMaxB*y(3)/(y(3)+KNB)*y(4)/(y(4)+KPB)*y(5)/(y(5)+KCB)*I/(I+KIB)*y(2)-kdB*y(2);
        dydt(3) = (-YNA*muMaxA*y(3)/(y(3)+KNA)*y(4)/(y(4)+KPA)*y(5)/(y(5)+KCA)*I/(I+KIA)*y(1)) + ...
            (-YNB*muMaxB*y(3)/(y(3)+KNB)*y(4)/(y(4)+KPB)*y(5)/(y(5)+KCB)*I/(I+KIB)*y(2));
        dydt(4) = (-YPA*muMaxA*y(3)/(y(3)+KNA)*y(4)/(y(4)+KPA)*y(5)/(y(5)+KCA)*I/(I+KIA)*y(1)) + ...
            (-YPB*muMaxB*y(3)/(y(3)+KNB)*y(4)/(y(4)+KPB)*y(5)/(y(5)+KCB)*I/(I+KIB)*y(2));
        dydt(5) = ktC*(SoCO2-y(5)) - (YCO2A*muMaxA*y(3)/(y(3)+KNA)*y(4)/(y(4)+KPA)*y(5)/(y(5)+KCA)*I/(I+KIA)*y(1) + ...
            YCO2B*muMaxB*y(3)/(y(3)+KNB)*y(4)/(y(4)+KPB)*y(5)/(y(5)+KCB)*I/(I+KIB)*y(2));
        dydt(6) = YT*y(3)/(y(4))*(muMaxA*y(3)/(y(3)+KNA)*y(4)/(y(4)+KPA)*y(5)/(y(5)+KCA)*I/(I+KIA)*y(1)-kdL*y(2));

        dydt = dydt';
    end
end
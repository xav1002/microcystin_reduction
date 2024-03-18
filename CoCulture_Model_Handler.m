%% Parameter and Data Loading
clear
close all
clc

global y0 temp cultureVol muMaxA muMaxB KNA KNB KPA KPB KCA KCB KIA KIB kdA kdB YNA YNB YPA YPB ktC Io A SoCO2 YCO2A YCO2B YT KHA KHB YH herbicide;

% Load Parameters
paramTable = readtable('CoCulture_Parameters.xlsx');
params = paramTable.Value;

% A - Microcystis aeruginosa, B - Synechococcus elongatus
temp = params(1); % Temperature
cultureVol = params(2); % Culture Volume
muMaxA = params(3); % Max Growth Rate for A
muMaxB = params(4); % Max Growth Rate for B
KNA = params(5); % Nitrogen Half Saturation Constant for A 
KNB = params(6); % Nitrogen Half Saturation Constant for B 
KPA = params(7); % Phosphorus Half Saturation Constant for A 
KPB = params(8); % Phosphorus Half Saturation Constant for B
KCA = params(9); % Carbon Half Saturation Constant for A 
KCB = params(10); % Carbon Half Saturation Constant for B 
KIA = params(11); % Light Half Saturation Constant for A
KIB = params(12); % Light Half Saturation Constant for B
kdA = params(13); % Death Constant for A
kdB = params(14); % Death Constant for B
YNA = params(15); % Nitrogen to Biomass Conversion (gN/gBiomass) for A
YNB = params(16); % Nitrogen to Biomass Conversion (gN/gBiomass) for B
YPA = params(17); % Phosphorus to Biomass Conversion (gP/gBiomass) for A
YPB = params(18); % Phosphorus to Biomass Conversion (gP/gBiomass) for B
ktC = params(19); % CO2 Mass Transfer Coefficient 
Io = params(20); % Incident Light Intensity
A = params(21); % Light Function Constant
SoCO2 = params(22); % Equilibrium Dissolved CO2 Concentration
YCO2A = params(23); % Carbon to Biomass Conversion (gC/gBiomass) for A
YCO2B = params(24); % Carbon to Biomass Conversion (gC/gBiomass) for B
YT = params(25); % Toxin-Cell Mass Fraction (gToxin/gA)
YT_orig = params(25); % Toxin-Cell Mass Fraction (gToxin/gA)
KHA = params(26); % Herbicide Half-Inhibition Constant for A 
KHB = params(27); % Herbicide Half-Inhibition Constant for B 
YH = params(28); % Herbicide Consumption Yield

% Model runtime (0h => 360h)
tSmooth = linspace(0,15*24,100);

% Initial Conditions
% y0(1) = A OD730
% y0(2) = B OD730
% y0(3) = N concentration (g/L)
% y0(4) = P concentration (g/L)
% y0(5) = CO2 concentration (g/L)
% y0(6) = Microcystin (g/L)
% y0(7) = Herbicide (chemical stress) concentration (g/L)
y0 = readtable('CoCulture_Init_Cond.xlsx').Value;
% 0.148 gAlgae/(L*OD730)
y0(1) = y0(1) * 0.148;
y0(2) = y0(2) * 0.148;
y0_2_orig = y0(2);
y0_7_orig = y0(7);

% Co-Culture
% Co-Culture + Herb
% Mono-Culture
% Mono-Culture + Herb
for k=1:1:4
    if k == 1
        y0(2) = y0_2_orig;
        y0(7) = 0;
        YT = YT_orig;
        herbicide = false;
        path = ".\Manuscript_Plots\CoCulture_NoHerb\";
        scenario = "CoCulture_NoHerb_";
    elseif k == 2
        y0(2) = y0_2_orig;
        y0(7) = y0_7_orig;
        YT = YT*80;
        herbicide = true;
        path = ".\Manuscript_Plots\CoCulture_Herb\";
        scenario = "CoCulture_Herb_";
    elseif k == 3
        y0(2) = 0;
        y0(7) = 0;
        YT = YT_orig;
        herbicide = false;
        path = ".\Manuscript_Plots\MonoCulture_NoHerb\";
        scenario = "MonoCulture_NoHerb_";
    elseif k == 4
        y0(2) = 0;
        y0(7) = y0_7_orig;
        YT = YT*80;
        herbicide = true;
        path = ".\Manuscript_Plots\MonoCulture_Herb\";
        scenario = "MonoCulture_Herb_";
    end
    
    % Nutrient condition variation
    nbpts = 5;
    Nvals = linspace(y0(3).*0.01, y0(3), nbpts);
    Pvals = linspace(y0(4).*0.01, y0(4), nbpts);
    FGAVals = zeros(length(Nvals), length(Pvals));
    FGBVals = zeros(length(Nvals), length(Pvals));
    FGMCVals = zeros(length(Nvals), length(Pvals));
    
    % Defining which N,P conditions the time-scale plots will be using
    NCond = 5;
    PCond = 5;
    
    % ode45 loops, with varying N and P initial concentrations
    for l=1:1:length(Nvals)
        for m=1:1:length(Pvals)
            y0(3) = Nvals(l);
            y0(4) = Pvals(m);
            [tCalc,yCalc] = coCulture_herbFunc(tSmooth);
            if l == NCond && m == PCond
                yCalcTS = yCalc;
            end
            % ### increasing P concentration as col# increases
            % ### increasing N concentration as row# increases
            FGAVals(m,l) = yCalc(end,1);
            FGBVals(m,l) = yCalc(end,2);
            FGMCVals(m,l) = yCalc(end,6).*(1E6);
        end
    end
    
    %% Plotting
    fontSize = 13;
    lineWidth = 1.2;
    font = 'Arial';
    fontWeight = 'bold';
    gridAlpha = 0.75;

    % Plotting Cell Concentration Kinetic Model (Co-Culture)
    figure;
    plot(tCalc, yCalcTS(:,1), 'r-', 'LineWidth', lineWidth);
    hold on;
    plot(tCalc, yCalcTS(:,2), 'Color', [0,0.7,0.2], 'LineWidth', lineWidth);
    xlabel("Time (h)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("Cell Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend('{\it Microcystis aeruginosa}', '{\it Synechococcus elongatus}', 'FontSize', 10, 'FontName', font, 'FontWeight', fontWeight);
    if k == 1
        lgd.Location = "east";
    else
        lgd.Location = "northwest";
    end
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    box off;
    ylim([0,1.5]);
    saveas(gcf,path+scenario+"MA_and_2973_Cell_Growth_Kinetics.png");
    hold off;

    % Plotting A Cell Concentration Kinetic Model (Co-Culture, but only A)
    figure;
    plot(tCalc, yCalcTS(:,1), 'r-', 'LineWidth', lineWidth);
    xlabel("Time (h)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("Cell Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend('{\it Microcystis aeruginosa}', 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    if k == 1
        lgd.Location = "east";
    else
        lgd.Location = "northwest";
    end
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    box off;
    ylim([0,1.5]);
    saveas(gcf,path+scenario+"MA_only_Cell_Growth_Kinetics.png");
    hold off;
    
    % Plotting Nutrient Concentration Kinetic Model
    figure;
    hold on;
    plot(tCalc, yCalcTS(:,3),'m-', 'LineWidth', lineWidth);
    plot(tCalc, yCalcTS(:,4),'k-', 'LineWidth', lineWidth);
    % plot(tCalc, yCalcTS(:,5),'Color',[0.8,0.6,0.1], 'LineWidth', lineWidth);
    xlabel("Time (h)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    % lgd = legend('N', 'P', 'CO2', 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend('N', 'P', 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd.Location = "east";
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    box off;
    saveas(gcf,path+scenario+"Nutrient_Kinetics.png");
    hold off;
    
    % Plotting Microcystin Kinetic Model
    % subplot(2,2,3);
    figure;
    hold on;
    plot(tCalc, yCalcTS(:,6).*(1E6), 'b-', 'LineWidth', lineWidth);
    xlabel("Time (h)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("Concentration ({\mug}/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend('Microcystin', 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd.Location = "east";
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    box off;
    ylim([0,10]);
    saveas(gcf,path+scenario+"Toxin_Kinetics.png");
    hold off;
    
    % Plotting Herbicide Kinetic Model
    figure;
    hold on;
    plot(tCalc, yCalcTS(:,7)*(1E3), 'c-', 'LineWidth', lineWidth);
    xlabel("Time (h)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("Concentration (mg/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend('DCMU', 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd.Location = "east";
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    box off;
    saveas(gcf,path+scenario+"Herbicide_Kinetics.png");
    hold off;
    
    %% Curve Fitting and Plotting Final Growth vs Nutrient Concentration    
    interpNbPts = 50;
    qNvals = linspace(Nvals(1),Nvals(end),interpNbPts);
    qPvals = linspace(Pvals(1),Pvals(end),interpNbPts)';
    [qGridNVals,qGridPVals] = meshgrid(qNvals,qPvals);
    
    % 2D interpolation for A final cell concentrations
    interpFGMAGrid = interp2(Nvals, Pvals, FGAVals, qNvals, qPvals, 'makima');
    % 2D interpolation for B final cell concentrations
    interpFG2973Grid = interp2(Nvals, Pvals, FGBVals, qNvals, qPvals, 'makima');
    % 2D interpolation for Microcystin Toxin final concentration
    interpFGMCGrid = interp2(Nvals, Pvals, FGMCVals, qNvals, qPvals, 'makima');

    % plotting final cell concentration of B
    figure;
    Color2(:,:,1) = ones(interpNbPts).*linspace(0,0,interpNbPts); % red
    Color2(:,:,2) = ones(interpNbPts).*linspace(0.5,1,interpNbPts); % green
    Color2(:,:,3) = ones(interpNbPts).*linspace(0.2,0.7,interpNbPts); % blue
    surf(qGridNVals,qGridPVals,interpFG2973Grid,Color2);
    xlabel("N Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("P Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    zlabel("Cell Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    zlim([0,1.5]);
    grid on;
    ax.GridAlpha = gridAlpha;
    saveas(gcf,path+scenario+"2973_Final_Cell_Concentrataion.png");
    
    % plotting final cell concentration of A
    figure;
    Color1(:,:,1) = ones(interpNbPts).*linspace(0.7,1,interpNbPts); % red
    Color1(:,:,2) = ones(interpNbPts).*linspace(0,0,interpNbPts); % green
    Color1(:,:,3) = ones(interpNbPts).*linspace(0,0,interpNbPts); % blue
    surf(qGridNVals,qGridPVals,interpFGMAGrid,Color1);
    xlabel("N Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("P Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    zlabel("Cell Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    zlim([0,1.5]);
    grid on;
    ax.GridAlpha = gridAlpha;
    saveas(gcf,path+scenario+"MA_Final_Cell_Concentration.png");
    
    % plotting final concentration of Microcystin
    figure;
    Color3(:,:,1) = zeros(interpNbPts).*linspace(0,0,interpNbPts); % red
    Color3(:,:,2) = ones(interpNbPts).*linspace(0.2,0.7,interpNbPts); % green
    Color3(:,:,3) = ones(interpNbPts).*linspace(0,1,interpNbPts); % blue
    surf(qGridNVals,qGridPVals,interpFGMCGrid,Color3);
    xlabel("N Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("P Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    zlab = zlabel("Microcystin ({\mug}/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    zlab.Position(2) = zlab.Position(2) + 0.005;
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    zlim([0,40]);
    grid on;
    ax.GridAlpha = gridAlpha;
    saveas(gcf,path+scenario+"Final_Toxin_Concentration.png");
    
    % plotting final A and B cell concentrations on same plot
    figure;
    hold on;
    surf(qGridNVals,qGridPVals,interpFGMAGrid,Color1);
    alpha 0.5;
    surf(qGridNVals,qGridPVals,interpFG2973Grid,Color2);
    alpha 0.5;
    xlabel("N Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    ylabel("P Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    zlabel("Cell Concentration (g/L)", 'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd = legend(["{\it Microcystis aeruginosa}","{\it Synechoccocus elongatus}"], ...
        'FontSize', fontSize, 'FontName', font, 'FontWeight', fontWeight);
    lgd.Location = "northeast";
    ax = gca;
    ax.FontSize = fontSize;
    ax.LineWidth = lineWidth;
    ax.FontName = font;
    ax.TickDir = 'out';
    zlim([0,1.5]);
    view(-45,30);
    grid on;
    ax.GridAlpha = gridAlpha;
    saveas(gcf,path+scenario+"MA_and_2973_Final_Cell_Concentration.png");

    disp("Iteration "+k+" finished")
end
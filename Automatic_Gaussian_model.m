%% Leonard Jacques
% 11/7/2024

% For determining the NLS model kinetics applicable to
% polarization reversal in the ferroelectric switching of any material.

% Refs:  1.  (Supplemental Information) Guido et al. (2024) Adv. Sci., 11(16), 2308797.
%        2.  Jo et al. (2007) Phys. Rev. Lett., 99(26), 267602.

% The first long while loop that spans several hundred lines used the
% arctan approximation for the NLS model which is much faster and helps to
% get parameters close to ideal.
% Then, the real NLS model with the integral is run in iterations with the
% appriximated fitting terms as initial conditions. Since the integral is
% more expensive to runs, it is efficient to start with the approximated
% values.
clc
clear
%% User Inputs
%% Import from File:
xfrac = readtable('Switching_Speed_Measurements\switching exported data\Raw_exported_data\TEK00791.CSV'); % Imports data from a .txt file that is tab delimited performs the function of "text to columns"
nsxfrac = readtable(['Switching_Speed_Measurements\switching exported data\Raw_exported_data\TEK00785.CSV']); % Imports non-switching data, later to be subtracted from the switching data. Variables specific to non-switching data have "ns" in front of it.
d = 236; % sample thickness in nanometers
thickness = d*10^-7; % film thickness (cm)
area = 7.85*10^-5; % sample area in cm^2 7.85*10^-5
tzerosp = 70*0.01; % the fraction of the setpoint voltage to call t0
ErrorStop = 100*0.01; % The percentage of the polarization where to stop calculating for the error in the regression calculations at the end of this script

%% Fitting Parameters
u = 900*10^-9; % half-width, half-maximum of the Lorentzian distribution (in nano-seconds)
s = 30*10^-9; % A normalization constant for fitting the NLS distribution to the experimental data (the middle of the distribution) - the position of this should be at the point of maximum current.

%% imported data management
% for the switching data
xfrac = table2array(xfrac); % converts the imported data type from "table' to "array"
nsxfrac = table2array(nsxfrac); % for the non-switching data
% finding where the function generator output begins
xa = xfrac(:,2); % calling CH1, the function generator output
xa = abs(xa); % making all values positive for referencing applications to follow
xaloc = find(xa==0.12); % 0.10 to 0.24V is the function generator output that signifies the output has begun !!! CHECK THIS VALUE IF THE CODE DOES NOT WORK !
n=isempty(xaloc); % this series of if loops is used to determine where the function generator input begins for the switching data
if n==1
    xaloc = find(xa==0.12);
    n=isempty(xaloc);
end
if n==1
    xaloc = find(xa==0.14);
    n=isempty(xaloc);
end
if n==1
    xaloc = find(xa==0.16);
    n=isempty(xaloc);
end
if n==1
    xaloc = find(xa==0.18);
    n=isempty(xaloc);
end
if n==1
    xaloc = find(xa==0.20);
    n=isempty(xaloc);
end
pos1 = xaloc(1) +100 -2; % holds the value where to start accepting points. +100 is for neglecting the first 100 ns before the amplifier output begins.
% ps 1 also governs the position of t0
% for the non-switchng data (ns)
nsxa = nsxfrac(:,2); % calling CH1, the function generator output
nsxa = abs(nsxa);
nsxaloc = find(nsxa==0.18); % > 0.10 V is the function generator output that signifies the output has begun  !!! CHECK THIS VALUE IF THE CODE DOES NOT WORK !!!
n=isempty(nsxaloc);
if n==1  % this series of if loops is used to determine where the function generator input begins for the non-switching data
    nsxaloc = find(nsxa==0.12);
    n=isempty(nsxaloc);
end
if n==1
    nsxaloc = find(nsxa==0.14);
    n=isempty(nsxaloc);
end
if n==1
    nsxaloc = find(nsxa==0.16);
    n=isempty(nsxaloc);
end
if n==1
    nsxaloc = find(nsxa==0.18);
    n=isempty(nsxaloc);
end
if n==1
    nsxaloc = find(nsxa==0.20);
    n=isempty(nsxaloc);
end
nspos1 = nsxaloc(1) +100 - 2; % holds the value where to start accepting points. +100 is for neglecting the first 100 ns before the amplifier output begins.

%% Adjusting for t0, and choosing the data to simulate
t0adjust = xfrac(pos1:end,3);
maxamp = mode(t0adjust); % assumes the mode (most common value) in the amplifier output is the voltage setpoint
fract0 = tzerosp*maxamp; % calculates the voltage where to start accepting
fract0 = abs(fract0); % for handling negative values
t0adjust = abs(t0adjust); % for handling negative values
t0pos = find(t0adjust>fract0); % finds the position of t0 in the amplifier output greater than the defined percentage of the setpoint voltage where to start the simulation

pos1 = pos1 + t0pos(1); % adjusts the pos1 time for the voltage where t0 is assigned
nspos1 = nspos1 + t0pos(1);
%nspos1=9814;
% shortens the arrays for the function generator (fg), amplifier (amp) and
% TIA (tia) switching curves
fg = xfrac(pos1:end,2);
amp = xfrac(pos1:end,3);
tia = abs(xfrac(pos1:end,4).*1000./22); % TIA output (in mA)
l = length(tia);

%% setting the time at polarization saturation
[maxtiaval,maxTIApos] = max(tia); % maximum current value and position from the TIA
tiafind = tia(maxTIApos:end);
zeroTIApos = find(tiafind < 0.5); % to index the zero positions after the TIA current peak
timelimit = maxTIApos + zeroTIApos(1) + 400; % to identify the location where the polarization is considered to have saturated

% non-switching curve
nsfg = nsxfrac(nspos1:end,2);
nsamp = nsxfrac(nspos1:end,3);
nstia = abs(nsxfrac(nspos1:end,4).*1000./22); % TIA output (in mA)
nsl = length(nstia);

% to set the proper array length of data to calculate
if l <= nsl
    l = l;
elseif l > nsl
    l = nsl;
end
l=9000;% sets the length scale to stop at 9 microseconds
%% Leakge current subtration from the TIA output
NetLeakageCurrent = mean(tia(3000:9000))-mean(nstia(3000:9000)); % calculated the mean of the TIA output beyond the switching transient from 3 to 9 microseconds
tia = tia(1:9000)-nstia(1:9000);
if NetLeakageCurrent > 0
    tia(timelimit+400:end) = tia(timelimit+400:end)-NetLeakageCurrent;
end
%% Imported data mamagement and calculation of experimental volume fraction transformed
% calculates the experimental polarization by subtracting a switching curve
% from the non-switching curve

time = [0.001:0.001:0.001*(l)]; % sets the time in microseconds
time = transpose(time);

i=1;
po = 0;
while i <= l
    polarization(i,1) = po + 10^3*((tia(i,1)-nstia(i,1))*(10^-9))/area; % solves for the polarization by subtracting the switching curve from the non-switching curve
    po = polarization(i);
    i=i+1;
end

expfraction = polarization./polarization(timelimit); % calls the sets the volume fraction of polarization
tcorrection = 0.2; % A correction factor which simultaneously shifts the ext. KAI model and projected nucleation rate forward to a lower time (in order to match the slope) and adjusts the x-axis limit of the plot accordingly.
time2 = [0.001:0.001:l*0.001]-tcorrection; % Sets the array for plotting of the ext. KAI model and the projected nucleation rate peak.

%% Finding maximum value of the current to assign the Avrami constant, n

maxCurrent = max(tia);
maxCurrent = find(tia==maxCurrent); % Assigns the position of the max current in the array
% maxCurrent = find(expfraction>=0.7);
Begin = maxCurrent(1) - 10;
Last = maxCurrent(1) + 10;

timevar = log((10^-9).*time([Begin:Last])); % changes unit of time from nanooseconds to seconds
polarizationvar = log(log(1./(1-expfraction([Begin:Last]))));
Avramin = polyfit(timevar,polarizationvar,1);
lnK = Avramin(2);
Avramin = Avramin(1)

%% Initial fitting of the NLS model
meanFind = find(expfraction>=0.5);
u = meanFind(1)*10^-9;
t = 10^-9; % (beginning of) time in the upper integrand (s)
i=1;
while i<=l
    f = 0.5*(1+erf((t-u)/(s*(2^0.5))));
    Pfraction(i) = f; % sets values in an array for transformation fraction
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end

%% Caussian Cumulative Distribution Function calculation
i = 1;
t = 10^-9; %
while i<=l
    NucGauss = (1/sqrt(2*pi*s^2))*exp(-((t-u)^2)/(2*s^2)); % Gaussian distribution of nucleation
    NucRatePeak(i) = NucGauss; % sets values in an array for the nucleation rate peak
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end
y = transpose(Pfraction); % transposed the Pfraction into columnar data which can be used later for exporting

%% Automatic determination of the fitting parameters t1 and w

% calculates the n value for the KAI model
KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
KAIn = polyfit(timevar,KAIvar,1);

%% Position finding with s
exploc = find(Pfraction>=0.5);
j=1;
while j<80
    Pfloc = find(expfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if KAIn(1)<Avramin
        s = s - (s/30);
    elseif KAIn(1)>Avramin
        s = s + (s/30);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1) = Avramin')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated NLS (Gaussian) model
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = 0.5*(1+erf((t-u)/(s*(2^0.5))));
        Pfraction(i) = f; % sets values in an array for transformation fraction
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
    
    % calculates the n value for the KAI model
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % this will end the while loop when the values begin to toggle
    testKAI(j) = KAIn(1); % stores and array of values of the slopes of the simulated curve
    if j>2
        if testKAI(j)==testKAI(j-2)
            break
        end
    end
    j=j+1;
end
clear testKAI

 % Nucleation rate peak calculation
i = 1;
    t = 10^-9; %
    ActualN = 0; % sets the initial value for ActualNuc, the nucleation rate
    while i<=l % this is the letter "l", not the numbner "1" !!!
        NucGauss = (1/sqrt(2*pi*s^2))*exp(-((t-u)^2)/(2*s^2)); % Gaussian distribution of nucleation
        NucRatePeak(i) = NucGauss; % sets values in an array for the nucleation rate peak
        ActualN(i+1) = ActualN(i) + NucGauss*10^-9; % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor
        ActualNuc(i) = ActualN(i+1); % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor - FOR PLOTTING
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
clear testKAI
% Solving for the regression sum of squares
ErStop = find(expfraction > ErrorStop); % Holds the positions of the polarization fractions over a certain fraction that was transformed. The first value in this array is the position where to stop calcualting the fitting errors.
i=1;
while i <= ErStop(1)
    SSRcal = (Pfraction(i) - expfraction(i))^2; % Regression sum of squares
    SSRhold(i) = SSRcal; % holds the value in a array
    i=i+1;
end
SSR = sum(SSRhold)/ErStop(1); % Regression sum of squares
% Solving the Total Sum of Squares
i=1;
while i <= ErStop(1)
    SSTcal = (expfraction(i) - mean(expfraction))^2; % Regression sum of squares
    SSThold(i) = SSTcal; % holds the value in a array
    i=i+1;
end
SST = sum(SSThold)/ErStop(1); % Total sum of squares
%Rsquare = 1 - SSR/SST; % R^2 value
Rsquare = 1 - SSR; % not exactly R^2 value - because there is only 1 data point for each value on the simulated curve

maxN = max(ActualNuc); % Holds the maximum value for nucleation
%NucPercent = ActualNuc./maxN; % Holds the percentage of nucleation completed
NucPercent = ActualNuc./ActualNuc(timelimit); % Holds the normalized percentage of nucleation completed

time = [0.001:0.001:0.001*(length(Pfraction))]; % sets the time in microseconds
time=transpose(time);
% Plot the extended KAI model vs. Experimental data after automatic fitting
% corrections

figure % creates a new figure
yyaxis left % sets up the left yaxis
p = plot(time,expfraction,time,Pfraction,'k'); % plotted data for the left y-axis
title('NLS (Gaussian dist.) comparison with experimental data') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('fraction of polarization reversal')
ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
xlim([0 1])
p(1).LineWidth = 3; % sets the width of lines in the plot
p(2).LineWidth = 2;
caption = sprintf(' n = %1.2s \n u(mean) = %1.2s sec \n sigma = %1.1s sec \n t_0 voltage s.p. = %1.1s \n SSR = %1.2s',Avramin,u,s,tzerosp,SSR); % text for the caption, limited to one decimal place by %1.1s
text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
%setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(time,NucRatePeak); % plotted data for the right y-axis
ylabel('Distribution Function')
legend('experimental','NLS model','Nuc. rate') % creates a legend for both the left and right y-axes in one box

%maxamp = (maxamp/thickness)*10^-6; % converts "maxamp" from voltage to electric field (MV/cm)
% Calculating the theoretical values for t1 and w
%E0 = 15; % Threshold electric field for pinned domains (MV/cm)
%U = 0.07; % eV/function unit in ZMO (Baska et al., 2024)
%T = 273+19; % Temperature in Kelvin
%kb = 8.617*10^-5; % Boltzmann's constant in eV
%t1theo = 10^(U*E0/kb/T/maxamp)


% Outputs the experimetal v.f., simulated v.f., and N/Ninf nucleation rate
% peak data in an array.
Pfraction = transpose(Pfraction);
NucRatePeak = transpose(NucRatePeak);
NucPercent=transpose(NucPercent);
PrOutput = [time,expfraction,Pfraction,NucRatePeak,NucPercent];
[maxNuc,maxNucPos] = max(NucRatePeak);
amp = abs(amp); % used to calcualte the amplitude of the applied field
AppliedField = (max(amp)/(10^6))/thickness; % The calculated applied field in MV/cm
Properties = [tzerosp,Avramin,SSR]; % a set of output properties
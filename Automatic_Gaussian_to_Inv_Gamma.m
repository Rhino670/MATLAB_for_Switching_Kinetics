%% Leonard Jacques
% 9/12/2024

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

%% TROUBLESHOOTING:
%  the "j" variable inequalities for the while loops in lines 343 and 380
%  may have to be altered slightly to change the number of iterations so
%  the loops do not "run away" and produce obscure values. The values
%  should be whole integers and the number of iterations to performw would
%  likely be between 3 and 7 (j<4 and j<8) depending on the curcumstances,
%  but may be outside those bounds.
clc
clear
%% User Inputs
%% Import from File:
xfrac = readtable('TEK00378.CSV'); % Imports data from a .txt file that is tab delimited performs the function of "text to columns"
nsxfrac = readtable(['TEK00336.CSV']); % Imports non-switching data, later to be subtracted from the switching data. Variables specific to non-switching data have "ns" in front of it.
d = 227; % sample thickness in nanometers
thickness = d*10^-7; % film thickness (cm)
area = 7.85*10^-5; % sample area in cm^2 7.85*10^-5
tzerosp = 70*0.01; % the fraction of the setpoint voltage to call t0
ErrorStop = 100*0.01; % The percentage of the polarization where to stop calculating for the error in the regression calculations at the end of this script

%% For modifying the tail of the switching curve: (1) its location (TLadd), and (2) where to begin leakage current subtraction.
LeakageAdd = -400; % The number of nanoseconds to add to the "timelimit" variable. Value of 2,000 for situations where there is a long tail. Value of 400 or around that value is ok for Gaussian distribution switching in ZMO.
TLadd = 400; % The number of nanoseconds before the timelimit for when to begin leakage current subtraction

%% Initial Fitting Parameters / Matching Parameters
% Gaussian Distribution
u = 900*10^-9; % half-width, half-maximum of the Lorentzian distribution (in nano-seconds)
s = 30*10^-9; % A normalization constant for fitting the NLS distribution to the experimental data (the middle of the distribution) - the position of this should be at the point of maximum current.
% Inverse Gamma Distribution
a = 4; % initial value for alpha (describes the shape) - must be a positive number
B = 4.8*10^-7; % initial value for lambda (describes the scale) - must be a positive value between 0 and infinity
TailMatch = 87*0.01; % the position along the tail where to match the inverse gaussian distribution to the experimental curve
ShapeMatch = 60*0.01; % the position along the curve where to match the inverse gaussian distribution to the experimental curve
% Combined Distribution Distribution
DifferenceCondition = 0.01; % Holds the threshold difference between the experimental curve and the Gaussian fit
GaussPre = 0.5; % Holds the threshold difference between the experimental curve and the Gaussian fit

%% imported data management
% for the switching data
xfrac = table2array(xfrac); % converts the imported data type from "table' to "array"
nsxfrac = table2array(nsxfrac); % for the non-switching data
% finding where the function generator output begins
xa = xfrac(:,2); % calling CH1, the function generator output
xa = abs(xa); % making all values positive for referencing applications to follow
xaloc = find(xa==0.18); % 0.10 to 0.24V is the function generator output that signifies the output has begun !!! CHECK THIS VALUE IF THE CODE DOES NOT WORK !
n=isempty(xaloc); % this series of if loops is used to determine where the function generator input begins for the switching data
if n==1 % this series of if loops is used to determine where the function generator input begins for the switching data
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
if n==1 % this series of if loops is used to determine where the function generator input begins for the non-switching data
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
timelimit = maxTIApos + zeroTIApos(1) + TLadd; % to identify the location where the polarization is considered to have saturated

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

l=9000; % sets the time time limit to 9 microseconds

%% Leakge current subtration from the TIA output
NetLeakageCurrent = mean(tia(3000:9000))-mean(nstia(3000:9000)); % calculated the mean of the TIA output beyond the switching transient from 3 to 9 microseconds
if NetLeakageCurrent > 0
    tia(timelimit+LeakageAdd:end) = tia(timelimit+LeakageAdd:end)-NetLeakageCurrent; % Can MODIFY the term here to modify where the begin the leakage current contributions
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

%% Finding maximum value of the current to assign the Avrami constant, n
maxCurrent = max(tia);
maxCurrent = find(tia==maxCurrent); % Assigns the position of the max current in the array
Begin = maxCurrent(1) - 10;
Last = maxCurrent(1) + 10;

timevar = log((10^-9).*time([Begin:Last])); % changes unit of time from nanooseconds to seconds
polarizationvar = log(log(1./(1-expfraction([Begin:Last]))));
Avramin = polyfit(timevar,polarizationvar,1);
lnK = Avramin(2);
Avramin = Avramin(1)

%% Gaussian Distribution

%% Initial fitting of the Gaussian model
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

%% Distribution Function calculation
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
        NucGauss(i) = (1/sqrt(2*pi*s^2))*exp(-((t-u)^2)/(2*s^2)); % Gaussian distribution of nucleation
        %NucRatePeak(i) = NucGauss; % sets values in an array for the nucleation rate peak
        ActualN(i+1) = ActualN(i) + NucGauss(i)*10^-9; % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor
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
SSRgauss = sum(SSRhold)/ErStop(1); % Regression sum of squares

% Plot the NLS Gaussian model vs. Experimental data after automatic fitting
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
caption = sprintf(' n = %1.2s \n u(mean) = %1.2s sec \n sigma = %1.1s sec \n t_0 voltage s.p. = %1.1s \n SSR = %1.2s',Avramin,u,s,tzerosp,SSRgauss); % text for the caption, limited to one decimal place by %1.1s
text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
%setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(time,NucGauss); % plotted data for the right y-axis
ylabel('Fraction of Nucleation Progress')
%ylim([-0.01 1]) % sets bounds for the right y-axis
legend('experimental','NLS model','Distribution Function') % creates a legend for both the left and right y-axes in one box

% Store the Parameters for the Gaussian Distribution in this order: Polarization Fraction, Distribution Function, Reduced
% sum of squares
PfracGauss = Pfraction; % stores the Gaussian fit Pfraction data
Pfraction = transpose(Pfraction);
NucRatePeak = transpose(NucGauss);
PfracGauss = Pfraction;
GaussStore = [Pfraction,NucRatePeak];

%% Inverse Gamma Distribution

%% Initial fitting of the inv. gamma NLS model
t = 10^-9; % (beginning of) time in the upper integrand (s)
i=1;
while i<=l
    fun = @(t) t.^(a-1).*exp(-t); % partial CDF (cumulative distribution function) of gamma (Erlang) function
    f = integral(fun,B/t,Inf);
    Pfraction(i) = f; % sets values in an array for transformation fraction
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end

%% Probability density function (Nucleation)
i = 1;
t = 10^-9;
while i<=l
    NucInvGamma = (B^a)*t^(-a-1)*exp(-B/t)/gamma(a); % (PDF) Probability distribution function of inverse gamma distribution
    NucRatePeak(i) = NucInvGamma; % sets values in an array for the nucleation rate peak
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end
y = transpose(Pfraction); % transposed the Pfraction into columnar data which can be used later for exporting

%% Automatic determination of the fitting parameters a (alpha) and B (Beta)

% calculates the n value for the KAI model
KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
KAIn = polyfit(timevar,KAIvar,1);


TailLoc = find(expfraction>=TailMatch);
TailLoc = TailLoc(1);
exploc = find(expfraction>=ShapeMatch);

z=1;
while z < 30
    % Position finding with B (scale)
j=1;
while j<4
    if Pfraction(TailLoc)<expfraction(TailLoc)
        B = B + (B/30);
    elseif Pfraction(TailLoc)>expfraction(TailLoc)
        B = B - (B/30);
    elseif Pfraction(TailLoc)==expfraction(TailLoc)
        sprintf('Pfloc(1) = exploc(1)')
    end
    B
    % Cumulative Density Function
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    ActualN = 0; % sets the initial value for ActualNuc, the nucleation rate
    while i<=l % this is the letter "l", not the numbner "1" !!!
        NucInvGamma = (B^a)*t^(-a-1)*exp(-B/t)/gamma(a); % (PDF) Probability distribution function of inverse gamma distribution
        NucRatePeak(i) = NucInvGamma; % sets values in an array for the nucleation rate peak
        ActualN(i+1) = ActualN(i) + NucInvGamma*10^-9; % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor
        ActualNuc(i) = ActualN(i+1); % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor - FOR PLOTTING
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
    maxN = max(ActualNuc); % Holds the maximum value for nucleation
    NucPercent = ActualNuc./ActualNuc(timelimit); % Holds the normalized percentage of nucleation completed
    Pfraction = NucPercent;
    % This breaks the while loop when values begin to toggle around a minimum
    BStore(j) = B;
    if j>2
        if BStore(j)==BStore(j-2)
            break
        end
    end
    j=j+1;
end

%% Shape (body) matching with a (alpha)
exploc = find(expfraction>=ShapeMatch);
j=1;
while j<6
    Pfloc = find(Pfraction>=ShapeMatch); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        a = a - (a/50);
    elseif Pfloc(1)>exploc(1)
        a = a + (a/50);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
    a
    
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    ActualN = 0; % sets the initial value for ActualNuc, the nucleation rate
    while i<=l % this is the letter "l", not the numbner "1" !!!
        NucInvGamma = (B^a)*t^(-a-1)*exp(-B/t)/gamma(a); % (PDF) Probability distribution function of inverse gamma distribution
        NucRatePeak(i) = NucInvGamma; % sets values in an array for the nucleation rate peak
        ActualN(i+1) = ActualN(i) + NucInvGamma*10^-9; % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor
        ActualNuc(i) = ActualN(i+1); % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor - FOR PLOTTING
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
    maxN = max(ActualNuc); % Holds the maximum value for nucleation
    NucPercent = ActualNuc./ActualNuc(timelimit); % Holds the normalized percentage of nucleation completed
    Pfraction = NucPercent;

    % This breaks the while loop when values begin to toggle around a minimum
    aStore(j) = a;
    if j>2
        if aStore(j)==aStore(j-2)
            break
        end
    end
    j=j+1;
    %clear Pfloc
end
z=z+1;
end

% Solving for the regression sum of squares (error)
ErStop = find(expfraction > ErrorStop); % Holds the positions of the polarization fractions over a certain fraction that was transformed. The first value in this array is the position where to stop calcualting the fitting errors.
i=1;
while i <= ErStop(1)
    SSRcal = (Pfraction(i) - expfraction(i))^2; % Regression sum of squares
    SSRhold(i) = SSRcal; % holds the value in a array
    i=i+1;
end
SSR = sum(SSRhold)/ErStop(1); % Regression sum of squares

maxN = max(ActualNuc); % Holds the maximum value for nucleation
%NucPercent = ActualNuc./maxN; % Holds the percentage of nucleation completed
NucPercent = ActualNuc./ActualNuc(timelimit); % Holds the normalized percentage of nucleation completed

time = [0.001:0.001:0.001*(length(Pfraction))]; % sets the time in microseconds
time=transpose(time);
% Plot the Inverse Gamma NLS vs. Experimental data after automatic fitting
figure % creates a new figure
yyaxis left % sets up the left yaxis
p = plot(time,expfraction,time,Pfraction,'k'); % plotted data for the left y-axis
title('NLS (Inverse gamma dist.) comparison with experimental data') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('fraction of polarization reversal')
ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
xlim([0 1])
p(1).LineWidth = 3; % sets the width of lines in the plot
p(2).LineWidth = 2;
caption = sprintf(' n = %1.2s \n a = %1.3s sec \n B = %1.3s sec \n t_0 voltage s.p. = %1.1s \n SSR = %1.3s',Avramin,a,B,tzerosp,SSR); % text for the caption, limited to one decimal place by %1.1s
text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
% setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(time,NucRatePeak); % plotted data for the right y-axis
ylabel('Distribution Function')
legend('experimental','NLS model','Dist. Function') % creates a legend for both the left and right y-axes in one box

%InvGamNuc = transpose(NucRatePeak);
PfracInvGamma = transpose(Pfraction);
InvGammaStore = [PfracInvGamma,NucRatePeak];

%% Combination Algorithm
% Combines the Gaissuan and Inverse Gamma curves into one so the tail of
% the inverse gamma distribution will fit the switching curve

i=1;
j=1;
while i <= l
    if expfraction(i) <= GaussPre % Establishes where the Gaussian distribution closely fits the experimental data
        CombinedPfrac(i) = PfracGauss(i); % Establishes the vector for the combined cumulative distribution function with Gaussian for the body and Inv. Gamma for teh tail
        CombinedNucDist(i) = NucGauss(i); % Establishes the vector for the combined distribution curve with Gaussian for the body and Inv. Gamma for the tail
    elseif expfraction(i) > GaussPre
        GaussDifference = abs(expfraction(i) - PfracGauss(i));
        if GaussDifference <= DifferenceCondition
            CombinedPfrac(i) = PfracGauss(i);
            CombinedNucDist(i) = NucGauss(i);
        elseif GaussDifference > DifferenceCondition
            CombinedPfrac(i) = Pfraction(i);
            CombinedNucDist(i) = NucRatePeak(i); 
            TransitionPfrac(j) = expfraction(i); % Tracks at chich polarization feaction the transition from Gaussian distribution to Inv. Gamma distribution occurs
            o=1;
            j=j+1;
        end
    end
    i=i+1;
end

i=1;
while i <= ErStop(1)
    SSRcal = (CombinedPfrac(i) - expfraction(i))^2; % Regression sum of squares
    SSRhold(i) = SSRcal; % holds the value in a array
    i=i+1;
end
SSR = sum(SSRhold)/ErStop(1); % Regression sum of squares

% Entropy
H = a + log(B*gamma(a)) - ((a+1)*psi(a))

% Time at 
t90 = find(expfraction>=0.9);
t90(1)

TransitionPfrac = TransitionPfrac(1)

% Plot the Combined NLS vs. Experimental data after automatic fitting
figure % creates a new figure
yyaxis left % sets up the left yaxis
p = plot(time,expfraction,time,CombinedPfrac,'k'); % plotted data for the left y-axis
title('Combined NLS comparison with experimental data') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('fraction of polarization reversal')
ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
xlim([0 1.5])
p(1).LineWidth = 3; % sets the width of lines in the plot
p(2).LineWidth = 2;
caption = sprintf(' n = %1.2s \n Trans = %1.2s \n t_0 voltage s.p. = %1.1s \n SSR = %1.3s',Avramin,TransitionPfrac,tzerosp,SSR); % text for the caption, limited to one decimal place by %1.1s
text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
% setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(time,CombinedNucDist); % plotted data for the right y-axis
ylabel('Distribution Function')
legend('experimental','NLS model','Dist. Function') % creates a legend for both the left and right y-axes in one box

CombinedPfrac = transpose(CombinedPfrac);
NucPercent=transpose(NucPercent);
PrOutput = [time,expfraction,-CombinedPfrac,NucRatePeak];
amp = abs(amp); % used to calcualte the amplitude of the applied field
AppliedField = (max(amp)/(10^6))/thickness; % The calculated applied field in MV/cm
Properties = [tzerosp,a,B,Avramin,SSR]; % a set of output properties
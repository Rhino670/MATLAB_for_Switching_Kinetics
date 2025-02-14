%% Leonard Jacques
% November 2024

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
xfrac = readtable('TEK00768.CSV'); % Imports data from a .txt file that is tab delimited performs the function of "text to columns"
nsxfrac = readtable(['TEK00767.CSV']); % Imports non-switching data, later to be subtracted from the switching data. Variables specific to non-switching data have "ns" in front of it.
d = 236; % sample thickness in nanometers
thickness = d*10^-7; % film thickness (cm)
area = 7.85*10^-5; % sample area in cm^2 7.85*10^-5
tzerosp = 70*0.01; % the fraction of the setpoint voltage to call t0
ErrorStop = 100*0.01; % The percentage of the polarization where to stop calculating for the error in the regression calculations at the end of this script

%% Fitting Parameters
w = 0.2; % half-width, half-maximum of the Lorentzian distribution (in nano-seconds)
A = 1.02; % A normalization constant for fitting the NLS distribution to the experimental data
t1 = 900*10^-9; % center of the distribution function

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
l=9000; % sets the time time limit to 9 microseconds

%% Leakge current subtration from the TIA output
NetLeakageCurrent = mean(tia(3000:9000))-mean(nstia(3000:9000)); % calculated the mean of the TIA output beyond the switching transient from 3 to 9 microseconds
if NetLeakageCurrent > 0
    tia(timelimit+800:end) = tia(timelimit+800:end)-NetLeakageCurrent; % Can MODIFY the term here to modify where the begin the leakage current contributions
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
Begin = maxCurrent(1) - 10;
Last = maxCurrent(1) + 10;

timevar = log((10^-9).*time([Begin:Last])); % changes unit of time from nanooseconds to seconds
polarizationvar = log(log(1./(1-expfraction([Begin:Last]))));
Avramin = polyfit(timevar,polarizationvar,1);
lnK = Avramin(2);
Avramin = Avramin(1);

m = Avramin-1 % Avramin-2; % sets the fitting parameter "m" to n-2.
m
%% Initial fitting of the NLS model
t = 10^-9; % (beginning of) time in the upper integrand (s)
i=1;
while i<=l
    f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
    Pfraction(i) = f; % sets values in an array for transformation fraction
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end

%% Nucleation Lorentzian calculation
i = 1;
t = 10^-9; %
while i<=l
    NucLorentz = (A/pi)*(w/((log(t)-log(t1))^2 + w^2)); % Lorentzian distribution of nucleation
    NucRatePeak(i) = NucLorentz; % sets values in an array for the nucleation rate peak
    t=t+(10^-9); % adds 1 nanosecond to the time
    i=i+1;
end
y = transpose(Pfraction); % transposed the Pfraction into columnar data which can be used later for exporting


%figure % for plotting the data used to determine the Avrami n
%plot(timevar,polarizationvar)

%% Plot the extended KAI model vs. Experimental data
%figure % creates a new figure
%yyaxis left % sets up the left yaxis
%p = plot(time,expfraction,(time),Pfraction,'k'); % plotted data for the left y-axis
%title('NLS model comparison with experimental data') % plot title
%xlabel('Time (microseconds)') % generates axis titles
%ylabel('fraction of polarization reversal')
%ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
%p(1).LineWidth = 3; % sets the width of lines in the plot
%p(2).LineWidth = 2;
%caption = sprintf(' n = %1.2s \n w = %1.2s sec \n t1 = %1.1s sec',Avramin,w,t1); % text for the caption, limited to one decimal place by %1.1s
%text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
%setting up the right y-axis
%yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
%plot(time2,NucRatePeak); % plotted data for the right y-axis
%ylabel('Distribution Function')
%legend('experimental','NLS model','Distribution Function') % creates a legend for both the left and right y-axes in one box

%% Automatic determination of the fitting parameters t1 and w

% calculates the n value for the KAI model
KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
KAIn = polyfit(timevar,KAIvar,1);

%% Position finding with t1
exploc = find(expfraction>=0.5);
j=1;
while j<20
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/15);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/15);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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

%% Position finding with t1
j=1;
while j<20
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/50);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/50);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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
%% Position finding with t1
j=1;
while j<20
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/250);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/250);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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
%% Position finding with t1
j=1;
while j<20
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/1000);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/1000);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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
%% Slope matching with w - modifying the fitting parameter, w
j=1;
while j<40 % the maximum number of iterations to run
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    Begin = Pfloc(1) - 10;
    Last = Pfloc(1) + 10;
    % Determinining the slope of the simulated curve in order to compare it to
    % the experimental volume fracton transformed.
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % for modifying the power exponent of a
    if KAIn(1)<Avramin
        w = w - (w/15);
    elseif KAIn(1)>Avramin
        w = w + (w/15);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1)=Avramin')
    end
testm = [KAIn(1),Avramin]
w
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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

%% Slope matching with w - modifying the fitting parameter, w
j=1;
while j<40 % the maximum number of iterations to run
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    Begin = Pfloc(1) - 10;
    Last = Pfloc(1) + 10;
    % Determinining the slope of the simulated curve in order to compare it to
    % the experimental volume fracton transformed.
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % for modifying the power exponent of a
    if KAIn(1)<Avramin
        w = w - (w/50);
    elseif KAIn(1)>Avramin
        w = w + (w/50);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1)=Avramin')
    end
testm = [KAIn(1),Avramin]
w
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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

%% Slope matching with w - modifying the fitting parameter, w
j=1;
while j<40 % the maximum number of iterations to run
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    Begin = Pfloc(1) - 10;
    Last = Pfloc(1) + 10;
    % Determinining the slope of the simulated curve in order to compare it to
    % the experimental volume fracton transformed.
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % for modifying the power exponent of a
    if KAIn(1)<Avramin
        w = w - (w/250);
    elseif KAIn(1)>Avramin
        w = w + (w/250);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1)=Avramin')
    end
testm = [KAIn(1),Avramin]
w
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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

%%% Height matching with A - modifying the fitting parameter, A
%j=1;
%while j<40 % the maximum number of iterations to run
%Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
%Begin = Pfloc(1) - 10;
%    Last = Pfloc(1) + 10;
%    % Determinining the slope of the simulated curve in order to compare it to
%    % the experimental volume fracton transformed.
%    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
%    KAIn = polyfit(timevar,KAIvar,1);
%    % for modifying the power exponent of a
%    if Pfraction(timelimit) < expfraction(timelimit)
%        A = A + 0.005;
%    elseif Pfraction(timelimit) > expfraction(timelimit)
%        A = A - 0.005;
%    elseif Pfraction(timelimit) == expfraction(timelimit)
%        sprintf('Pfraction(timelimit) == expfraction(timelimit)')
%    end
%testm = [KAIn(1),Avramin]
%w
%    % Simulated KAI model
%    % setting up the integral
%    t = 10^-9; % (beginning of) time in the upper integrand (s)
%    i=1;
%    while i<=l
%        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
%        Pfraction(i) = f; % sets values in an array for transformation fraction
%        t=t+(10^-9); % adds 1 nanosecond to the time
%        i=i+1;
%    end
%
%    % calculates the n value for the KAI model
%    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
%    KAIn = polyfit(timevar,KAIvar,1);
%    % this will end the while loop when the values begin to toggle
%    testKAI(j) = Pfraction(timelimit); % stores and array of values of the simulated curve at the "timelimit" where the response is condisered to have saturated
%    if j>2
%        if testKAI(j)==testKAI(j-2)
%            break
%        end
%    end
%    j=j+1;
%end
%clear testKAI

%% Slope matching with w - modifying the fitting parameter, w
j=1;
while j<40 % the maximum number of iterations to run
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    Begin = Pfloc(1) - 10;
    Last = Pfloc(1) + 10;
    % Determinining the slope of the simulated curve in order to compare it to
    % the experimental volume fracton transformed.
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % for modifying the power exponent of a
    if KAIn(1)<Avramin
        w = w - (w/250);
    elseif KAIn(1)>Avramin
        w = w + (w/250);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1)=Avramin')
    end
testm = [KAIn(1),Avramin]
w
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
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

%% Position finding with t1
j=1;
while j<100
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/300);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/100);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        f = (A/pi)*(atan((log(t) - log(t1))/w) + (pi/2));
        Pfraction(i) = f; % sets values in an array for transformation fraction
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
    y = transpose(Pfraction); % transposed the Pfraction into columnar data which can be used later for exporting
    % calculates the n value for the KAI model
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % this will end the while loop when the values begin to toggle
    testKAI(j) = KAIn(1); % stores an array of values of the slopes of the simulated curve
    if j>2
        if testKAI(j)==testKAI(j-2)
            break
        end
    end
    j=j+1;
end

%% reconverting from the approximation with the arctan to the actual NLS
% model using the obtained parameters

% ngrow=Avramin-1; % this sets the exponent "n" in the integral to the Avrami n constant from the experimental curve
ngrow=4.75;

%% Position finding with t1
j=1;
while j<20
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/30);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/30);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        fun = @(logt) (1-exp(-(t./exp(logt)).^ngrow)).*(A./pi).*(w./((logt-log(t1)).^2 + w.^2));
        f = integral(fun,-Inf,Inf);
        Pfraction(i) = f; % sets values in an array for transformation fraction
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
   
    % calculates the n value for the KAI model
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % this will end the while loop when the values begin to toggle
    testKAI(j) = test2peakloc(1); % stores and array of values of the slopes of the simulated curve
    if j>2
        if test2peakloc(1)==testKAI(j-2)
            break
        end
    end
    j=j+1;
end
clear testKAI

%% Position finding with t1
j=1;
while j<5
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    if Pfloc(1)<exploc(1)
        t1 = t1 + (t1/90);
    elseif Pfloc(1)>exploc(1)
        t1 = t1 - (t1/90);
    elseif Pfloc(1)==exploc(1)
        sprintf('Pfloc(1) = exploc(1)')
    end
testm = [KAIn(1),Avramin]
test2peakloc= [Pfloc(1),exploc(1)]
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        fun = @(logt) (1-exp(-(t./exp(logt)).^ngrow)).*(A./pi).*(w./((logt-log(t1)).^2 + w.^2));
        f = integral(fun,-Inf,Inf);
        Pfraction(i) = f; % sets values in an array for transformation fraction
        t=t+(10^-9); % adds 1 nanosecond to the time
        i=i+1;
    end
   
    % calculates the n value for the KAI model
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % this will end the while loop when the values begin to toggle
%    testKAI(j) = test2peakloc(1); % stores and array of values of the slopes of the simulated curve
%    if j>2
%        if test2peakloc(1) == test2peakloc(1)
%            break
%        elseif test2peakloc(1)==testKAI(j-2)
%            break
%        end
%    end
    j=j+1;
end
clear testKAI

%% Slope matching with w - modifying the fitting parameter, w
j=1;
while j<40 % the maximum number of iterations to run
    Pfloc = find(Pfraction>=0.5); % holds the index values where the value of Pfraction, the calculated KAI curve, is greater than 0.5, in order to calculate the n value there.
    Begin = Pfloc(1) - 10;
    Last = Pfloc(1) + 10;
    % Determinining the slope of the simulated curve in order to compare it to
    % the experimental volume fracton transformed.
    KAIvar = log(log(1./(1-Pfraction([Begin:Last]))));
    KAIn = polyfit(timevar,KAIvar,1);
    % for modifying the power exponent of a
    if KAIn(1)<Avramin
        w = w - (w/40);
    elseif KAIn(1)>Avramin
        w = w + (w/40);
    elseif KAIn(1)==Avramin
        sprintf('KAIn(1)=Avramin')
    end
testm = [KAIn(1),Avramin]
w
    % Simulated KAI model
    % setting up the integral
    t = 10^-9; % (beginning of) time in the upper integrand (s)
    i=1;
    while i<=l
        fun = @(logt) (1-exp(-(t./exp(logt)).^ngrow)).*(A./pi).*(w./((logt-log(t1)).^2 + w.^2));
        f = integral(fun,-Inf,Inf);
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
        NucLorentz = (A/pi)*(w/((log(t)-log(t1))^2 + w^2)); % Lorentzian distribution of nucleation
        NucRatePeak(i) = NucLorentz; % sets values in an array for the nucleation rate peak
        ActualN(i+1) = ActualN(i) + NucLorentz*10^-9; % stores the number of actual nuclei which have formed during the polarization reversal process in the capacitor
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
title('Automatic NLS model comparison with experimental data') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('fraction of polarization reversal')
ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
xlim([0 1])
p(1).LineWidth = 3; % sets the width of lines in the plot
p(2).LineWidth = 2;
caption = sprintf(' n = %1.2s \n w = %1.2s sec \n t1 = %1.2s sec \n t_0 voltage s.p. = %1.1s \n A = %1.3s \n SSR = %1.3s',Avramin,w,t1,tzerosp,A,SSR); % text for the caption, limited to one decimal place by %1.1s
text(0.1,0.6,caption); % generates a text box at the coordinates (x,y,caption) according to the left y-axis
%setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(time,NucPercent); % plotted data for the right y-axis
ylabel('Fraction of Nucleation Progress')
ylim([-0.01 1]) % sets bounds for the right y-axis
legend('experimental','NLS model','# Nuclei') % creates a legend for both the left and right y-axes in one box

figure % creates a new figure
yyaxis left % sets up the left yaxis
p = plot(time,expfraction,time,Pfraction,'k'); % plotted data for the left y-axis
title('Automatic NLS model comparison with experimental data') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('fraction of polarization reversal')
ylim([-0.01 1]) % sets bounds for the left y-axis
%xlim([0 (2-tcorrection)])
xlim([0 1])
p(1).LineWidth = 3; % sets the width of lines in the plot
p(2).LineWidth = 2;
caption = sprintf(' n = %1.2s \n w = %1.2s sec \n t1 = %1.2s sec \n t_0 voltage s.p. = %1.1s \n A = %1.3s \n SSR = %1.3s',Avramin,w,t1,tzerosp,A,SSR); % text for the caption, limited to one decimal place by %1.1s
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
Properties = [tzerosp,A,Avramin,SSR]; % a set of output properties
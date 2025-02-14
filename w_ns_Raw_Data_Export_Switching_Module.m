%% Leonard Jacques
%  4/10/2024

%% Switching Speed Oscilloscope Export Manager
% Import data from Excel file containing three channel output from the
% Tektronix Osxilloscope, where:
% CH1: function generator output
% CH2: signal amplifier output (Trex, 50x)
% CH3: TIA (sample) output (divide voltage signal by the gain, 22 ohms, to
% obtain the current in Amperes.
%
% Ouptut: 3 plots, one worksheet with data
% function generator output vs. amplifier output
% amplifier output vs. TIA output (mA)
% amplifier output vs. polarization (uC/cm^2)
% 1 worksheet (columns from left to right):
% fg (V), amp (MV/cm), tia (mA), p (uC/cm^2)


% Results are stored in y. The first column is the applied field, the
% following columns are the dielectric constants, and the following columns
% are the dielectric losses.

% Instructions
% 1) Select the .txt file to import, and store in variable "x".
% 2) Insert the dielectric's thickness, t in nm and and area, a in cm^2

clear
%% User Inputs
%% Import from File:
x = readtable(['TEK00768.CSV']); % Imports data from a .txt file that is tab delimited performs the function of "text to columns"
nsx = readtable(['TEK00767.CSV']); % Imports non-switching data, later to be subtracted from the switching data. Variables specific to non-switching data have "ns" in front of it.
d = 220; % sample thickness in nanometers
a = 2.83*10^-5; % sample area in cm^2
pl = 6; % plot x-axis limit (in microseconds)


%% imported data management
% for the switching data
x = table2array(x); % converts the imported data type from "table' to "array"
nsx = table2array(nsx); % for the non-switching data
% finding where the function generator output begins
xa = x(:,2); % calling CH1, the function generator output
xa = abs(xa); % making all values positive for referencing applications to follow
xaloc = find(xa==0.12); % 0.10 to 0.24V is the function generator output that signifies the output has begun !!! CHECK THIS VALUE IF THE CODE DOES NOT WORK !
n=isempty(xaloc); % this series of if loops is used to determine where the function generator input begins for the switching data
if n==1;
    xaloc = find(xa==0.12);
    n=isempty(xaloc);
end
if n==1;
    xaloc = find(xa==0.14);
    n=isempty(xaloc);
end
if n==1;
    xaloc = find(xa==0.16);
    n=isempty(xaloc);
end
if n==1;
    xaloc = find(xa==0.18);
    n=isempty(xaloc);
end
if n==1;
    xaloc = find(xa==0.20);
    n=isempty(xaloc);
end
pos1 = xaloc(1) -100% +90; % holds the value where to start accepting points. +90 is for neglecting th efirst 90 ns before the amplifier output begins.
% for the non-switchng data (ns)
nsxa = nsx(:,2); % calling CH1, the function generator output
nsxa = abs(nsxa);
nsxaloc = find(nsxa==0.18); % > 0.10 V is the function generator output that signifies the output has begun  !!! CHECK THIS VALUE IF THE CODE DOES NOT WORK !!!
n=isempty(nsxaloc);
if n==1;  % this series of if loops is used to determine where the function generator input begins for the non-switching data
    nsxaloc = find(nsxa==0.12);
    n=isempty(nsxaloc);
end
if n==1;
    nsxaloc = find(nsxa==0.14);
    n=isempty(nsxaloc);
end
if n==1;
    nsxaloc = find(nsxa==0.16);
    n=isempty(nsxaloc);
end
if n==1;
    nsxaloc = find(nsxa==0.18);
    n=isempty(nsxaloc);
end
if n==1;
    nsxaloc = find(nsxa==0.20);
    n=isempty(nsxaloc);
end
nspos1 = nsxaloc(1) -100% +90; % holds the value where to start accepting points. +90 is for neglecting the first 90 ns before the amplifier output begins.

% shortens the arrays for teh function generator (fg), amplifier (amp) and
% TIA (tia)
fg = x(pos1:end,2);
amp = x(pos1:end,3);
tia = abs(x(pos1:end,4).*1000./22); % TIA output (in mA)
l = length(tia);
% for non-switching
nsfg = nsx(nspos1:end,2);
nsamp = nsx(nspos1:end,3);
nstia = abs(nsx(nspos1:end,4).*1000./22); % TIA output (in mA)
nsl = length(nstia);

% to set the proper array length of data to calculate
if l <= nsl
    l = l;
    nstia = nstia(1:l);
elseif l > nsl
    l = nsl; % shortens array lengths for the switching portion
    fg = fg(1:l);
    amp = amp(1:l);
    tia = tia(1:l);
end

% shortens the arrays for teh function generator (fg), amplifier (amp) and
% TIA (tia)
t = [0.001:0.001:0.001*l];
t = transpose(t);

i=1;
po = 0;
while i <= l
    p(i,1) = po + 10^3*((tia(i,1)-nstia(i,1))*(10^-9))/a; % solves for the polarization by subtracting the switching response of the TIA from the non-switching (ns) response of the TIA
    po = p(i);
    tiasum(i) = tia(i)-nstia(i);
    i=i+1;
end
%% function generator output vs. amplifier output
%figure
%yyaxis left % sets up the left yaxis
%plot(t,fg); % plotted data for the left y-axis
%title('Switching Speed Export with Matlab') % plot title
%xlabel('Time (microseconds)') % generates axis titles
%ylabel('Function Generator Output (V)')
%xlim([0 pl])
%setting up the right y-axis
%yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
%plot(t,amp); % plotted data for the right y-axis
%ylabel('Amplifier output (V)')
%legend('Function Generator Output','Amplifier Output') % creates a legend for both the left and right y-axes in one box

amp = amp./(10^6)./(d*10^-7); % converts voltage to applied field

tiadensity=tia./a; %concerts the Tia to current density

set(groot, 'DefaultFigurePosition', [3000, 200, 500, 500])

%% amplifier output vs. TIA output (mA)
figure
yyaxis left % sets up the left yaxis
plot(t,amp); % plotted data for the left y-axis
title('Switching Speed Export with Matlab') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('Amplifier Output (MV/cm)')
xlim([0 pl])
%setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(t,tiasum); % plotted data for the right y-axis
ylabel('TIA Output (mA)')
legend('Amplifier Output','TIA Output') % creates a legend for both the left and right y-axes in one box

%% amplifier output vs. polarization (uC/cm^2)
figure
yyaxis left % sets up the left yaxis
plot(t,amp); % plotted data for the left y-axis
title('Switching Speed Export with Matlab') % plot title
xlabel('Time (microseconds)') % generates axis titles
ylabel('Amplifier Output (MV/cm)')
xlim([0 pl])
%setting up the right y-axis
yyaxis right % sets up the right y axis. syntax below will correspond to the right y-axis.
plot(t,p); % plotted data for the right y-axis
ylabel('Polarization (uC/cm^2)')
legend('Amplifier Output','Polarization') % creates a legend for both the left and right y-axes in one box

% Output fg (V), amp (MV/cm), tia (mA), p (uC/cm^2)
Output = [fg,amp,transpose(tiasum),-p]; % Outputs the specified data in a 4-column array; Can also add time by placing a column "t".
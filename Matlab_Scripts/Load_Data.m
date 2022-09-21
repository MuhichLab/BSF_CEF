function [X,Y,Z,x_exp,Temp,Press,dd,muhg_o] = Load_Data(DFT_file,EXP_file)

%% Load DFT data

dft_data = table2array(readtable(DFT_file));
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X = dft_data(:,2); % mol fract Ba
Y = dft_data(:,1); % delta
Z = dft_data(:,3)/8; % E eV per ABO3


%% Load Experimental data
exp_data = readtable(EXP_file);

Press = table2array(exp_data(:,2));  % partial pressure O2 in bar
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

% %To Remove specific pieces of data - example Temp
% remove = (Temp(:) == 811.15);
% Press(remove)=[];
% Temp(remove) = [];
% x_exp(remove )= [];
% dd(remove) = [];

% %To Remove specific pieces of data - example x
% remove = (x_exp(:) == 0.10);
% Press(remove)=[];
% Temp(remove) = [];
% x_exp(remove )= [];
% dd(remove) = [];

% Calculate O2 chemical potential (convert to eV/K)

Ho = [];
So = [];
for t = 1:length(Temp);
    [H,S] = get_O2_thermo(Temp(t));
    Ho(t) = H;
    So(t) = S;
end
conv = 96.487;  % 1 eV = 96.487 kJ/mol
R_gas = 8.314462/1000/conv; %eV/K
Ho = transpose(Ho);
So = transpose(So);

muhg_o2 = Ho/conv - (Temp).*So/conv + R_gas.*(Temp).*log(Press);

muhg_o = 0.5*muhg_o2;

end
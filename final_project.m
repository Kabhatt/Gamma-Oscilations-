Cm = 1;
Iapp = 2;
Iapp2 = @(t) Iapp;
gL = 0.1;
EL = -65;
ENa = 55;
EK = -90;
gNa =  35;
gK =  9;
gsys= 0.8; % paper uses 0.1 - requires running simulation for longer
Esyn = -75;
theta = 0;
phi = 1;

% Defining the opening rate (alpha) and Closing rate (beta) of k and Na
% channels:
a_m = @(V) -0.1*(V + 35) ./ (exp(-0.1 * (V + 35)) - 1);
b_m = @(V) 4 * exp(-(V + 60) / 18);

a_n = @(V) (-0.01)*(V+34)./(exp(-0.1*(V+34)) - 1);
b_n = @(V) 0.125 * exp(-(V + 44)/80);
m_inf = @(V) a_m(V) ./ (a_m(V) + b_m(V));
n_inf = @(V) a_n(V) ./ (a_n(V) + b_n(V));
a_h = @(V) 0.07 * exp(-(V + 58)/20);
b_h = @(V) 1 ./ (exp(-0.1*(V + 28)) + 1);
sF =@(V) 1./(1+exp(-((V-theta)/2)));
a_s = 12;
b_s = 0.1;

%functions
% eqns = [diff(V,t)== -gNa*m^3*h*(V-ENa)-gK*n^4*(V-EK)+Iapp/Cm];

dt = 0.1;
t = 0:dt:10000;
V1 = zeros(1, length(t));
V2 = zeros(1, length(t));
V1(1) = EL;
V2(1) = EL - 10;

h1 = 1;
n1 = 0;
s1 = 0;
h2 = 1;
n2 = 0;
s2 = 0;

for i = 1:length(t) - 1
    %cell1 k1
    k1v1 = (Iapp ...
        - gL * (V1(i) - EL) ...
        - gNa * m_inf(V1(i))^3 * h1 * (V1(i) - ENa) ...
        - gK * n1^4 * (V1(i) - EK) ...
        - gsys * s2 * (V1(i) - Esyn) ...
        )/Cm;
    k1h1 = phi * (a_h(V1(i)) * (1 - h1) - b_h(V1(i)) * h1);
    k1n1 = phi * (a_n(V1(i)) * (1 - n1) - b_n(V1(i)) * n1);
    k1s1 = sF(V2(i))*(1-s1) - b_s*s1;
    %cell2 k1
    k1v2 = (Iapp2(t(i)) ...
        - gL * (V2(i) - EL) ...
        - gNa * m_inf(V2(i))^3 * h2 * (V2(i) - ENa) ...
        - gK * n2^4 * (V2(i) - EK) ...
        - gsys * s1 * (V2(i) - Esyn) ...
        )/Cm;
    k1h2 = phi * (a_h(V2(i)) * (1 - h1) - b_h(V2(i)) * h2);
    k1n2 = phi * (a_n(V2(i)) * (1 - n1) - b_n(V2(i)) * n2);
    k1s2 = sF(V1(i)) * (1-s2) - b_s*s2;

    % half steps
    av1 = V1(i) + k1v1*dt;
    h_sub1= h1 + k1h1*dt;
    n_sub1= n1 + k1n1*dt;
    s_sub1 = s1 + k1s1*dt;
    %cell2 halfsteps
    aV2 = V2(i) + k1v2*dt;
    h_sub2= h2 + k1h2*dt;
    n_sub2= n2 + k1n2*dt;
    s_sub2 = s2 + k1s2*dt;

    %cell1
    k2v1 = (Iapp ...
        - gL * (av1 - EL) ...
        - gNa * m_inf(av1)^3 * h_sub1 * (av1 - ENa) ...
        - gK * n_sub1^4 * (av1 - EK) ...
        - gsys * s_sub2 * (av1 - Esyn) ...
        )/Cm;
    k2h1 = phi * (a_h(av1) * (1 - h_sub1) - b_h(av1) * h_sub1);
    k2n1 = phi * (a_n(av1) * (1 - n_sub1) - b_n(av1)*n_sub1);
    k2s1 = sF(av1) * (1-s_sub1) - b_s*s_sub1;
    %cell2
    k2V2 = (Iapp2(t(i+1)) ...
        - gL * (aV2 - EL) ...
        - gNa * m_inf(aV2)^3 * h_sub2 * (aV2 - ENa) ...
        - gK * n_sub2^4 * (aV2 - EK) ...
        - gsys * s_sub1 * (aV2 - Esyn) ...
        )/Cm;
    k2h2 = phi * (a_h(aV2) * (1 - h_sub2) - b_h(aV2) * h_sub2);
    k2n2 = phi * (a_n(aV2) * (1 - n_sub2) - b_n(aV2)*n_sub2);
    k2s2 = sF(aV2) * (1-s_sub2) - b_s*s_sub2;

    V1(i+1)= V1(i)+dt*(k1v1+k2v1)/2;
    h1 = h1 + dt*(k1h1+k2h1)/2;
    n1 = n1 + dt*(k1n1+k2n1)/2;
    s1 = s1 + dt* (k1s1 + k2s1) / 2;
    V2(i+1)= V2(i)+dt*(k1v2+k2V2)/2;
    h2 = h2 + dt*(k1h2+k2h2)/2;
    n2 = n2 + dt*(k1n2+k2n2)/2;
    s2 = s2 + dt* (k1s2 + k2s2) / 2;
end

close all;
tiledlayout(3,1);
nexttile
plot(t, V1, 'b');
hold on;
nexttile
plot(t, V2, 'r', linestyle='--');
hold on
nexttile
plot(t, V1, 'b');
hold on;
plot(t, V2, 'r', linestyle='--');
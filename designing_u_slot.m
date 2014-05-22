% Design procedure from: Analysis and Design of Broad-Band Single-Layer
%                        Rectangular U-Slot Microstrip Patch Antennas

%% Initialize
clear;
c = 3*10^8;
%% Step 1: Specify center frequency and upper/lower frequencies
%fres3 = 2.3*10^9;
%fres2 = 2.2*10^9;
%fres4 = 2.4*10^9;
fres3 = 2.15*10^9;
fres2 = 1.8*10^9;
fres4 = 2.5*10^9;

%% Step 2: Select substrate permittivity and thickness
er = 2.2;
% Rule of thumb for T: T >= 0.06*lambda_res3/sqrt(er)
Tmin = 0.06*(c/fres3)/sqrt(er);
% Since we want a thin design, choose the minimum thickness (for now at
% least, we could put another expression here, like T = 2*Tmin or something)

% T = 0.0032; % 3.2 mm
T = Tmin;
%T = 0.00635; % 6.35mm - from paper example 1

%% Step 3: Estimate B + 2*deltaB ~= v_0/(2*sqrt(er)*fres3)
bPlusBFringe = c/2/sqrt(er)/fres3; % assuming v_0 = speed of light

%% Step 4: Calculate A = 1.5 * ( B + 2*deltaB)
A = 1.5*bPlusBFringe;

%% Step 5: Calculate eeff and 2*deltaB
eeff = (er+1)/2 + (er-1)/2*(1 + 12*T/A)^(-1/2);
bFringe = 0.824*T*(eeff+0.3)/(eeff-0.258)*(A/T+0.262)/(A/T+0.813);

%% Step 6: Back calculate the value of B
B = c/2/sqrt(eeff)/fres3 - bFringe;

%% Step 7: Select slot thickness (rule of thumb)
E = c/fres3/60;
F = E;

%% Step 8: Calculate D
%D = c/sqrt(eeff)/fres2 - 2*(B + bFringe - E);
D = c/sqrt(eeff)/fres2 - 2*(bPlusBFringe - E);

%% Step 9: Select C such that: C/A >= 0.3 and C/D >= 0.75
Cmin1 = A*0.3;
Cmin2 = D*0.75;
C = max(Cmin1, Cmin2);
while (true)
    %% Step 10: Effective permittivity and effective length of psuedo-patch 
    eeffpp = (er+1)/2 + (er-1)/2*(1 + 12*T/(D-2*F))^(-1/2);
    psuedoFringe = 0.824*T*(eeffpp+0.3)/(eeffpp-0.258)*((D-2*F)/T+0.262)/((D-2*F)/T+0.813);

    %% Step 11: Calculate H
    H = B-E+psuedoFringe-1/sqrt(eeffpp)*(c/fres4-(2*C+D));

    %% Step 12: Check that it worked out!
    if C + E + H >= B
        fprintf('Done, but C+E+H >= B. Tweak C in step 9 to make H physically realizable.\r\n');
        fprintf('C + E + H = %.3e\r\n', C + E + H);
        fprintf('B = %.3e\r\n', B);
        fprintf('==============\r\n');
        fprintf('Current C = %.5e\r\n', C);
        in = input('What do you want the C to be now? Enter 0 to quit:  ');
        if in == 0
            fprintf('QUITTING\r\n');
            break;
        else
            C = in;
            fprintf('\r\nRunning script with C = %.5e\r\n', C);
        end
    else
        fprintf('Done! Looks good! Result in mm:');
        result.A = A*1000;
        result.B = B*1000;
        result.C = C*1000;
        result.D = D*1000;
        result.E = E*1000;
        result.F = F*1000;
        result.H = H*1000;
        %result.R = R;
        result.T = T*1000;
        %result.offset = offset;
        result.er = er;
        display(result);
        break;
    end
end
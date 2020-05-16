% Design of a Two stage Miller Operational Amplifier
% Script for calculating (W/L) ratios and currents

clear screen
%constant parameters 
vgs = 0.9; %gate source voltage 
vds = 0.6; %drain source voltage >= 0.6V

vtn = 0.5; % threshold voltage for NMOS
vtp = 0.55; % threshold voltage for PMOS
vtn_max = 1.2*vtn;
vtn_min = 0.8*vtn;
vtp_max= 1.2*vtp;
vtp_min = 0.6*vtp; %initially it was vtp_min = 0.8*vtp
cox = 6e-15/1e-6; %oxide layer capacitance = 6fF/um2
lambda_n = 0.21; %channel lenght modulation for NMOS
lambda_p = 0.15; %channel lenght modulation for PMOS
kn = 140e-6; % kn = un*cox
kp = 70e-6; % kp = up*cox
L = 240e-9; %minimum length of the transistor recommended 

VDD = 1.8; %Positive Supply voltage 
VSS = 0; %Negative supply voltage 
SR = 20e6; %Slew rate >= 15V/uS
GBW = 80e6; %Gain bandwidth should be greater than equal to 20MHz
Vin_min = 0.5; %input common mode voltage range (0.5V to 1.3V). Higher Vmin is chosen because of -ve VDS_sat
Vin_max = 1.3; %input common mode voltage range (0.5V to 1.3V)
Vout_min = 0.3; %Output voltage min requirement
Vout_max = 1.5; %Output voltage maximum requirement

%Calculation starts from here...
cl = 6e-12; %load capacitance
cc = 2e-12; %coupling capacitance

if(cl < 2e-12)
    disp('invalid load capacitance value, CL should be greater than or equal to 2pF\n')
    if(cc < (0.22*cl))
        fprintf('\n');
        disp('invalid coupling capacitance value. Cc should be greater than 0.22*CL')
    end 
else 
    %tail current calculation based on slew rate
    i5 = SR*cc; 
    fprintf('\nCurrent through NMOS_5, i5 = %g', i5);
    
    %Aspect ratio calculation for PMOS (W/L)3 = (W/L)4
    s3 = i5/(kp*((VDD-Vin_max-vtp_max+vtn_min)^2));
    fprintf('\nS3 = S4 = %f', s3);
    w3 = s3*L;
    fprintf('\t => W3 = W4 = %g', w3); 
    
    %Verify pole of M3
    cgs3 = 0.67*w3*L*cox; %gate-source capacitance of PMOS_3
    gm3 = sqrt((2*kp*s3*i5)/2); %transconductance for PMOS_3
    
    if((gm3/(2*cgs3))<10*GBW)
        fprintf('\nInvalid pole of M3');
    else 
    %Design for S1 = S2
    gm1 = GBW*2*pi*cc; %transconductance of NMOS_1
    fprintf('\ngm1 = %g', gm1);
    s2 = (gm1^2)/(kn*i5);
    fprintf('\nS1 = S2 = %f', s2);
    w2 = s2*L;
    fprintf('\t => W1 = W2 = %g', w2); 
    
    %Design for (W/L)5
    beta1 = kn*s2; 
    Vds_sat = Vin_min - VSS - sqrt((2*i5)/beta1) - vtn_max; %VDS satuartion voltage for the NMOS_5
    if(Vds_sat <= 100e-3)
        fprintf('\nNOTE: The value of Vds_sat should not be less than 100mV,\nto solve this problem: i5 can be reduced or (W/L)5 incresed');
        fprintf('\nVds_sat for the NMOS_5 = %g', Vds_sat);
        s5 = (2*i5)/(kn*(Vds_sat^2));
        fprintf('\nS5 = %f', s5);
        w5 = s5*L;
        fprintf('\t => W5 = %g', w5);
    else
        fprintf('\nVds_sat for the NMOS_5 = %g', Vds_sat);
        s5 = (2*i5)/(kn*(Vds_sat^2));
        fprintf('\nS5 = %f', s5);
        w5 = s5*L;
        fprintf('\t => W5 = %g', w5);       
    end
    
    %Design for (W/L)6
    gm6 = 2.2*gm1*(cl/cc);
    gm4 = sqrt((2*kp*s3*i5)/2);
    s6 = (gm6/gm4)*s3;
    fprintf('\nS6 = %f', s6);
    w6 = s6*L;
    fprintf('\t => W6 = %g', w6);
    %current passing through PMOS_6
    i6 = (gm6^2)/(2*kp*s6); 
    fprintf('\ni6 = %g', i6);
    
    %Design for (W/L)7
    s7 = (i6/i5)*s5;
    fprintf('\nS7 = %g', s7);
    w7 = s7*L;
    fprintf('\t => W7 = %g', w7);
    
    %Design for (W/L)8
    s8 = s6/(1+(cl/cc));
    w8 = s8*L;
    fprintf('\nW8 = %g', w8);
    fprintf('\r\n');
    end 
    
    %Gain Specifications
    i2 = 0.5*kn*s2*((vgs-vtn)^2)*(1+lambda_n*vds);
    i4 = 0.5*kp*s3*((vgs+vtp)^2)*(1-lambda_p*vds);
    i7 = 0.5*kn*s7*((vgs-vtn)^2)*(1+lambda_n*vds);
    fprintf('i2 = %g', i2);
    fprintf('i4 = %g', i4);
    fprintf('i7 = %g', i7);
    fprintf('gm6 = %g', gm6);
    %gain = 20*log10((2*gm1*gm6)/(i5*(i2+i4)*i6(i6+i7)));
    %fprintf('\nGain = %f', gain);
    %Design for Biasing circuit
    fprintf('Design for Biasing circuit');
end





function detResp = detectorModel(thick)

% Function to compute the spectral response of a detector with a CsI:Tl
% scintillator of a given thickness
%
% Input:
%        thick: Thickness of the scintillator (mm)
%


%e_g2 = 0.834;      % Probability of interaction with K shell - Richard
%w_g2 = 0.870;      % Quantum yield of K fluorescence - Richard
e_g2 = S02_xi(1);
w_g2 = S02_w(1);
%fk   = 0.78;       % Probability of absorption of K photon - Richard

s = ones(150,1); % Spectrum, just ones

en_K = S02_K_edge(1);

[fK f_K_vec] = S02_f_K(s,thick,1);

[g1Av g1] = S01_g1(s,thick,1);

g2_escape = S02_g2_escape(1);

xi_w = zeros(1,150);
xi_w(en_K:150) = e_g2*w_g2*1.2;

[g2aAv eg2a g2a]= S02_g2a(s, thick, g2_escape, 1, e_g2, w_g2);

[g2bAv eg2b g2b]= S02_g2b(s, thick, g2_escape, 1, e_g2, w_g2);

[g2cAv eg2c g2c]= S02_g2c(s, thick, g2_escape, 1, e_g2, w_g2, f_K_vec);

g4 = 0.45;
%detResp = (g1.*((1-xi_w).*g2a + xi_w.*g2b + xi_w.*f_K_vec.*g2c))*g4;
detResp = (g1.*((1-xi_w).*g2a + xi_w.*g2b + xi_w.*f_K_vec.*g2c))*g4;

% 
% for en = 1:150,
%     if (en > en_K),
%         detResp(en) = (g1(en)*(g2a(en)+g2b(en)+g2c(en)))*g4;
%     else
%         detResp(en) = (g1(en)*(g2a(en)))*g4;
%     end
        

end
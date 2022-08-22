%% Emmetteur pi/4-DQPSK/QPSK avec pre-synchronisation temporelle et frequentielle 
close all;
clear;
%% Parametres QPSK
%Te=2*1e-7; % Temps d'echantillonnage des CNA et CAN
Ts=2e-6; % Temps symbole
Fse = 6;
Te = Ts/Fse;
nb = 2; % Nombre de bits/symbole
M=2^nb; % Nombre de symboles de la constellation 
fe=1/Te; % Frequence d'echantillonnage
Ds=1/Ts; % Debit symbole
Db=nb*Ds; % Debit binaire
%Fse=floor(Ts/Te); % Facteur de sous-echantillonnage
f0=1.255e9; % Porteuse pour le calcul du PL
roll_off=0.8;
Nfft=2048; % Nombre de points sur lequel la FFT est effectuee
A=127; % Amplitude
Nb_paquets=48;  % Nb_paquets >=2
Octet_verif = hexToBinaryVector('0xAA');


%% Construction du signal de synchro
preambule = [ones(1,Fse) zeros(1,Fse) ones(1,Fse) zeros(1,2*Fse) ones(1,Fse) zeros(1,Fse) ones(1,Fse) zeros(1,4*Fse-1)];
preambule = repmat(preambule,1,2); % repmat() : repete le motif de la matrice preambule sur une fois en ligne et deux fois en colonne

%% Emetteur
I = imread('background600_400.jpg');
[rows, columns, numberOfColorChannels] = size(I);
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);
R= R(:)';
G= G(:)';
B= B(:)';

sb_Te = [];
for i=1:rows*columns
    sb_Te = [sb_Te R(i) G(i) B(i)];
end


sb_Te_bit_str =dec2bin(sb_Te,8);
sb_Te_bit_str = reshape(sb_Te_bit_str',1,[]);

sb_Te = [];
for j =1:length(sb_Te_bit_str)
    sb_Te(j) = str2num(sb_Te_bit_str(j));
end


Nb = length(sb_Te); % Nombre de bits emis
sb_paquets=reshape(sb_Te,[],Nb_paquets);
sb_paquets = [repmat(Octet_verif.',1,Nb_paquets); sb_paquets];

sl_Te=[];
sl_Te_test=0;
for j=1:(Nb+Nb_paquets*length(Octet_verif))/length(sb_paquets)
    disp(j);
    ssd_Ts = zeros(1,((Nb+Nb_paquets*length(Octet_verif))/2/Nb_paquets)+2);
    %% Modulateur numerique 
        thetaInit = pi/4; % Theta initial pour le mode differentiel
        ssd_Ts(1) = exp(1i*thetaInit); % Premier symbole differentiel
        ss_Ts=zeros(1,(Nb+Nb_paquets*length(Octet_verif))/2/Nb_paquets);
        ss_Ts_1 = exp(1i*thetaInit); % Premier symbole pilote pour correction de phase donc s_b(1:2) = [0 0]
        

        for i=1:(Nb+Nb_paquets*length(Octet_verif))/2/Nb_paquets % bits vers symbole mode diff
            if sb_paquets((i-1)*2+1:i*2,j)==[0;0]
                ss_Ts(i)=exp(1i*pi/4);
            elseif sb_paquets((i-1)*2+1:i*2,j)==[0;1]
                ss_Ts(i)=exp(1i*3*pi/4);
            elseif sb_paquets((i-1)*2+1:i*2,j)==[1;1]
                ss_Ts(i)=exp(1i*5*pi/4);
            elseif sb_paquets((i-1)*2+1:i*2,j)==[1;0]
                ss_Ts(i)=exp(1i*7*pi/4);
            end
            ssd_Ts(i+1)=ssd_Ts(i)*ss_Ts(i); % Creation du symbole differentiel au rang i+1
            
        end
        disp(length(ss_Ts));
        disp(length(ssd_Ts));
        ssd_Ts =[ss_Ts_1,ssd_Ts];
        disp(length(ssd_Ts));
        ss_Te = upsample(ssd_Ts,Fse); % Surechantillonnage de s_sd au rythme Fse 
        disp(length(ss_Te));
        

    %% Filtre de mise en forme
        g=rcosdesign(roll_off,2*M, Fse, 'sqrt');
        sl_Te_conv = conv(ss_Te,g); % Attention lors d'un filtrage numerique, les deux vecteurs a filtrer doivent être sur le même rythme d'echantillonnage !
            if(j==1)
            temp=sl_Te_conv;
            end
        disp(length(sl_Te_conv));
    %% Insertion du preambule et scalling
        sl_Te = [sl_Te,preambule*100,sl_Te_conv.*A/max(abs(sl_Te_conv))]; % Multiplication par A pour pouvoir profiter de l'amplification maximale disponible
        disp(length(sl_Te));
        disp(A/max(abs(sl_Te_conv)));

end
figure, plot(ss_Ts,'*');
title('Constellation I/Q de ss_{Ts}')
xlabel('I')
ylabel ('Q')
%% Adaptation a la radiologicielle
sl_Te=sl_Te(1:end);
yl_I=real(sl_Te);
yl_Q=imag(sl_Te);


yl_Tx=zeros(1,2*length(sl_Te));

yl_Tx(1:2:end)=yl_I;
yl_Tx(2:2:end)=yl_Q;

%% Enrengistrement de la sequence emise 
save("sb.mat","sb_Te");
fidID=fopen('QPSK_Tx.raw','w');
fwrite(fidID,int8(yl_Tx),'int8'); % Sauvegarde des echantillons modules
%fwrite(fidID,int16(yl_Tx),'short')
fclose(fidID);


%% Sauvegarde de données de filtre
csvwrite('filter.txt',g) ;
csvwrite('preambule.txt',preambule) ;
csvwrite('preambule_fliplr.txt',fliplr(preambule)) ;
data_filter = [M Fse roll_off];
csvwrite('Data_filter.txt',data_filter) ;

% 
% % Verification de la nature des signaux 
% % comparaison en DSP theorique et DSP de welch
% figure (1);
% subplot 121
% [DSP_Welch, f] = pwelch(sl_Te(length(preambule):end), ones(1,Nfft), 0, Nfft, fe, 'centered');
% 
% DSP_theo=zeros(size(f));
% for i=1:length(f)
%     if abs(f(i))<(1-roll_off)/(2*Ts)
%         DSP_theo(i)=10*2*A^2/(Ds*Fse);
%     elseif abs(f(i))<(1+roll_off)/(2*Ts)
%         DSP_theo(i)=10*A^2*(1+cos(pi*Ts/roll_off*(abs(f(i))-(1-roll_off)/(2*Ts))))/(Ds*Fse); %s'obtient avec une transformation bilineaire
%     end
% end
% 
% semilogy(f,DSP_Welch,f,DSP_theo);
% grid on;
% title('TF[sl_{Te}(t)] en fonction de la frequence')
% xlabel('Frequence (Hz)')
% legend('Methode de Welch','DSP theorique')
% ylabel ('Amplitude (V²)')
% xlim([-fe/2 fe/2]);
% ylim([10^(-5) 1]);
% 
% subplot 122
scatter(real(sl_Te),...
        imag(sl_Te),'o');
grid on
title('Constellation I/Q de sl_{Te}[n]')
xlabel('I')
ylabel ('Q')

% Diagramme de l'oeil
eyediagram(sl_Te(length(preambule):50000+length(preambule)+1),2*Fse,2*Ts);




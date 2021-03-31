%% loading data 
Dbin=0.05; %0.03 %0.02
Lwin=0.5;  %1.5 %0.5
receiver_z=50;
rec_fov = [20 45 90 180 20 45 90 180 20 45 90 180 20 45 90 180].*pi./180;
num_fov=length(rec_fov);
rec_aperture = [ones(num_fov/4,1).*0.2; ones(num_fov/4,1).*0.4; ones(num_fov/4,1).*0.6; ones(num_fov/4,1).*0.8];
%load('s2.mat','All_Received_Photons');
%All_Received_Photons=All_Received_Photons{1,1};
total_rec_packets=size(All_Received_Photons,1);

num_rx=length(rec_aperture);    %number of reciever
rec_pos = zeros(num_rx,2);  % reciever center

Photon_to_Rec_distance=zeros(total_rec_packets,num_rx);
for j = 1:num_rx             % iterate over every receiver
    rx_x = rec_pos(j,1);
    rx_y = rec_pos(j,2);
    radius = rec_aperture(j)/2;                           % 1/2 diameter of receiver
    cos_rec_fov = cos(rec_fov(j)/2);                    % cos(fov/2) to compare with photon's incident angle
    
    Photon_to_Rec_distance(:,j) = sqrt((rx_x-All_Received_Photons(:,1)).^2 + (rx_y-All_Received_Photons(:,2)).^2);
    Receiver_index = find((Photon_to_Rec_distance(:,j)<=radius) & (All_Received_Photons(:,6) >= cos_rec_fov));
    Reciever_Photons{j}=All_Received_Photons(Receiver_index,1:9);
    %total_time_recived_photon{j}=total_time_photon(Receiver_index,1);
    
    %% Impulse response
    Num_bin_photon {j}= floor((Reciever_Photons{j}(:,9)-receiver_z) / Dbin)+1;       % Partition the distances past the balistic distance into bins of distBinWidth
    Num_bin = round(Lwin/Dbin);
    % Hist_distance= zeros(num_rx,Num_bin);  % Creat an array for the distance bins for each receiver
    
    % if (Num_bin_photon{j} < Num_bin)
    for In_H=1:Num_bin
        index_for_bin=find(Num_bin_photon {j}==In_H);
        Sum_Bin_Weight= sum(Reciever_Photons{j}(index_for_bin,7));
        Hist_distance(In_H,j)=Sum_Bin_Weight;
    end
    
    %          Index_IR=find(Num_bin_photon (:,j)< Num_bin);
    %         aa=Reciever_Photons{j}(Index_IR,7);
    %         Hist_distance(j,Num_bin_photon(Index_IR)+1) = Hist_distance(j,Num_bin_photon(Index_IR)+1)+aa';  % Add the received weight into the appropriate bin
      
end
figure;
rx=4;
n=5e7;
binn=(0:round(Lwin/Dbin)-1);
c=3*1e8;
Tbin=Dbin*1.33/c;
time=binn*Tbin;
time=time';
plot(time,Hist_distance(:,rx)/n);
xlabel('time (second)');
ylabel('normalized received power');
title('Impulse response-harbor water-8m- Divergence 0.5 rad');

%% draw frequency response
fs=1/Tbin;  %zv... sampling frequency
nn=size(Hist_distance(:,rx),1);
FR=fft(Hist_distance(:,rx)/n);
dB_norm_H=db(abs(FR)/max(abs(FR)));

%single sided spectrum...
P1 = dB_norm_H(1:nn/2+1);
P1(2:end-1) = P1(2:end-1);
fre=(0:nn/2)*fs/nn;

figure;
plot(fre,P1);
xlabel('frequency (Hz)');
ylabel('Frequency Response (dB)');
title('Frequency Response-harbor water-8m- Divergence 0.5 rad');

%% draw impulse response







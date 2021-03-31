% MC UWOC simulation Tool : (the modified version of the "PHOTONATOR"
% simulation tool by William Cox)
% This code models optical beam propagation through an underwater channel
% with both absorption and turbidity. 

% The effect of different optical beams are studied.
% The SPF is HG to speed up the simulation time.
%  Writen by Zahra Vali (Ph.D. student of Isfahan University of Technology) and Hamed Noori at July 2016 in Vancouver, UBC.


%% Initial parameters
tic
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));  % to change the seed to have different rand output each time

poolSize=4;
num_photons=1*1e6;
chn_real=1;
attenuationLength=7.55;
c=0.151; % 2.19 harbor  %0.398  coastal %0.151 clear   %0.043 pure
albedo=0.245;   %0.83 harbor   %coastal 0.55  % clear 0.245;   %pure 0.0581
beamDiverg=0;
beamWidth=0.001;

g=0.924;
Dbin = 0.06;       % Bin spacing is 10 mm or ~10GHz bandwidth
Lwin = 12 ;      % Record out to 2 meters from balstic length

% Use dimention limits on the receiver plane. This reduces the size of the dataset
rxXLimMax =32;
rxXLimMin =-32;
rxYLimMax =1;
rxYLimMin =-1;
zLimMin =0;

receiver_z = attenuationLength/c;
b =c * albedo;
a = c-b;

rec_pos = [0,0];
sizeRecPos = size(rec_pos);
num_rx = sizeRecPos(1);

%% CDF calculation
[cdf_scatter_old,angle_old] = generate_scatter_HG(g);
%[cdf_scatter_old,angle_old] = generate_scattercox('measured','petzold_harbor');
%[cdf_scatter_old,angle_old] = generate_scattercox('measured','petzold_clear');
angle=0:pi/100000:pi;
cdf_scatter=interp1(angle_old,cdf_scatter_old,angle);

%% initial photon movement  **************************************
photon = zeros(num_photons,9);
photon(:,7) = ones(num_photons,1);  % set weights to one
photon(:,8) = ones(num_photons,1);  % set all active
%[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6),type]=beamProfile_TEM_lens(num_photons,beamWidth,beamDiverg,'gaussian');
[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6),type]=beamProfile_plane_wave(num_photons,beamWidth);
%[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6)]=beamProfile_TEMnm_nolens(sigma,num_photons,beamWidth,beamDiverg,'gaussian');
%[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6),type]=beamProfile_sph_wave(num_photons,beamWidth,beamDiverg);

%Photon_Initial = photon;

%% open pool
% if isempty(gcp('nocreate'))
%     parpool('local',poolSize)
% end

%% Monte carlo main
for m=1: chn_real
    %tic
   
    [All_Received_Photons{m},Total_Received_Photons_MC{m}]= part2_MC_fnc(g,photon,rxXLimMax,rxXLimMin,rxYLimMax,rxYLimMin, zLimMin,num_photons,c,a,receiver_z,cdf_scatter);
    
    % toc
    
end

num_rx=16;
%autoArrangeFigures()
Hist_distance_total = zeros(num_rx,round(Lwin/Dbin));
Total_Received_Power=zeros(num_rx,chn_real);
Total_Num_Recived_Photon=zeros(num_rx,chn_real);

% obtain data for pdf of recieved photons
for Rec_index=1:num_rx
    for Chan_Index=1:chn_real
        Temp_total = Total_Received_Photons_MC{1,Chan_Index}{1, Rec_index};
        Total_Received_Power(Rec_index,Chan_Index)=sum(Temp_total(:,7));
        Total_Num_Recived_Photon(Rec_index,Chan_Index) =size(Temp_total,1);
        %         aa=Hist_distance{1,Chan_Index}(:, Rec_index);
        %         Hist_distance_total(Rec_index,:) = Hist_distance_total(Rec_index,:) +  aa';
    end
end

%% close pool
% if isempty(gcp('nocreate'))~= 1
%     delete(gcp('nocreate'));
% end

beep
beep

%% save data in file
%dataDir = '/zahra/phdtez/MC-turb/all-code/channel/multiple_scattering';
dataDir='D:\zv\result';

foldername = sprintf('outputData%s',datestr(now,'HH MM_yyyy-mm-dd'));
mkdir(sprintf('%s/%s',dataDir,foldername));
fid = fopen(sprintf('%s/%s/whatsHere.txt',dataDir,foldername),'w');
fprintf(fid,'\n\r Atten.Length:%d \n\r Rx.dist:%d c:%d a:%d b:%d albedo:%d number of photon: %d chn real: %d beamwidth: %d beamDiverg: %d type: %d \n\r',attenuationLength,receiver_z, c,a,b, albedo, num_photons,chn_real,beamWidth, beamDiverg,type);
fclose(fid);
save(sprintf('%s/%s/simResults.mat',dataDir,foldername));

toc

% set up an array for the photons,
% x, y, z, mu_x, mu_y, mu_z, weight, received
% 0, 0, 0, 0   , 0   , 1   , 1     , 1 - initial values
% 1, 2, 3, 4   ,       5   ,     6   ,          7     , 8 - Position
% photon(:,1) == X POSITION
% photon(:,2) == Y POSITION
% photon(:,3) == Z POSITION
% photon(:,4) == ux
% photon(:,5) == uy
% photon(:,6) == uz
% photon(:,7) == WEIGHT
% photon(:,8) == RECEIVED? 1 = No, 0 = Yes (detected), -1 = terminated
% photon(:,9) == total distance

% index from origin, 0 theta (angle between x and z), 0 phi (angle between
% x and y) along x-axis

function [All_Received_Photons,Reciever_Photons] = part2_MC_fnc(g,photon,rxXLimMax,rxXLimMin,rxYLimMax,rxYLimMin, zLimMin,  num_photons,c,a,receiver_z,cdf_scatter)

%% initial conditions ************************************************
totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
total_time = zeros(num_photons,1);   % record total time each photon spent in the channel
each_time = zeros(num_photons,1);   % record  time each photon spent
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
rec_loc = zeros(num_photons,5);     % location of the received photon on the x,y rec. plane     
total_rec_packets = 0; % Total number of packets to cross detector
total_rec_power = 0; % Total power to cross detector
zLimMax = receiver_z;

test1=0;

prob_of_survival=((c-a)/c);       % b/c
max_uz = 1-1e-12;
rouletteConst = 10; % 1/rouletteConst determines the probability of a low power photon surviving
rouletteConst_inv = 1/rouletteConst;
inv_c = 1/c;

if prob_of_survival >= 0.90
    min_power = 1e-4; % minimum power value for photons before they are terminated by rouletting
    %min_power = 0.5
elseif prob_of_survival >= 0.8299
    min_power = 1e-5;
elseif prob_of_survival >= 0.70
    min_power = 1e-5;%1e-5;
else
    min_power = 1e-7;
end

%% initial photon movement **************************************
%[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6)]=beamProfile_TEM_lens(num_photons,beamWidth,beamDiverg,'gaussian');



%% after initial conditions ***********************************
photonsRemaining = num_photons; % count down for each photon received/terminated
%clear theta phi x y


while photonsRemaining> 0 % which is faster? create random values on the fly??  <- check
    
    %     %% trace photon...
    %     clear  trace_photon
    %     %x=2;
    %     trace_photon(:, 1:3)=photon(:, 1:3);
    %     %photon tracing....
    %     %figure
    %
    %     %         plot(trace_photon(3,:) ,trace_photon(1,:) )
    %     %         xlabel('z direction');
    %     %         ylabel('x direction ');
    %     %         %zlabel('z direction');
    %
    %     plot3(trace_photon(:,1) , trace_photon (:,2) , trace_photon (:,3) )
    %     xlabel('x direction');
    %     ylabel('y direction ');
    %     zlabel('z direction');
    %
    %     hold on
    %     pause(0.001)
    
    rand_array = rand(num_photons,3);      % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   % Uniformly distributed over (0,2Pi)
    
    % iterate over every single photon to calculate new position and whether it was received or not.
    for i = 1:num_photons
        
        if (photon(i,8) == 1)   % if the photon is still (1)active, (0)received, (-1)terminated
            
            %% step size calculation ******************
            % find new position increment based on PREVIOUS direction vectors (ux,uy,uz). This
            % takes care of the initial condition problem where photons
            % were scattered immediately on the first iteration.
            r =-(inv_c)*log(rand_array(i,1));     % randomly generate optical path length
%             Temp_index = find(abs(cdf_scatter-rand_array(i,2) )<0.0001); 
%             check_empty=isempty(Temp_index);
%             
%             if check_empty==1
%                 Temp_index = find(abs(cdf_scatter-rand_array(i,2) )<0.001);
%             end
%             theta= cdf_scatter(Temp_index(1));

            %for g!=0.......
            cc=rand_array(i,2);
            us=(1/(2*g)).*(1+ (g^2)- ((1-g^2)./(1+g-2.*g.*cc)).^2 );
            theta=acos(us);  %should be between zero and pi
            phi = rand_array(i,3);  % Phi is uniformly distributed over 0-2pi
            
            
            x_step = r*photon(i,4);  % x step size
            y_step = r*photon(i,5);  % y step size
            z_step= r*photon(i,6);  % z step size
            
            %[r,theta,phi,x_step,y_step,z_step]=step_size(rand_array,i,inv_c,cdf_scatter,photon);
            
            %% if photon pass reciver plane ***********************************
            if ((photon(i,3) + z_step) >= receiver_z)
                
                if photon(i,6) ~= 0  % If the photon has a z-component, mu_z != 0
                    
                    % z distance to receiver plane
                    z_dist_rec_intersection = receiver_z - photon(i,3);     % Zrec-Zphoton
                    % y distance to receiver plane
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);
                    % x distance to receiver plane
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);
                else
                    disp('how can the photon cross the receiver plane when it"s pointing straight up???');
                    
                end
                
                % euclidian distance to the reciever plane
                dist_to_rec = z_dist_rec_intersection / photon(i,6);    % z/mu_z
                
                % Photon exceeds the limits
                if ( (photon(i,1) + x_dist_rec_intersection) > rxXLimMax || (photon(i,1) + x_dist_rec_intersection) < rxXLimMin || (photon(i,2) + y_dist_rec_intersection) > rxYLimMax || (photon(i,2) + y_dist_rec_intersection) < rxYLimMin)
                    photon(i,8) =-1;  % mark as terminated
                    photonsRemaining = photonsRemaining-1;    % decrement outstanding photons
                    continue;
                end
                
                %[photonsRemaining,photon,x_dist_rec_intersection, y_dist_rec_intersection,dist_to_rec,check_out]=photon_passed_rcvplane(xLimMax,xLimMin,yLimMax,yLimMin,zLimMax,Limits_tank,i,photon,receiver_z ,rxXLimMax, rxXLimMin, rxYLimMax, rxYLimMin,photonsRemaining);
                
                photon(i,1) = photon(i,1) + x_dist_rec_intersection;   % x-axis location of reception
                photon(i,2) = photon(i,2) + y_dist_rec_intersection;    % y-axis location of reception
                photon(i,3)=receiver_z;
                photon(i,4) = photon(i,4); % for statistics, should be uniform (mu_x)
                photon(i,5) = photon(i,5); % for statistics, should be uniform (mu_y)
                photon(i,6) = photon(i,6); % incident angle, mu_z
                
                %                 trace_photon(i,1)=rec_loc(i,1);
                %                 trace_photon(i,2)=rec_loc(i,2);
                %                 trace_photon(i,3)=receiver_z;
                
                total_rec_packets = total_rec_packets + 1;
                total_rec_power = total_rec_power + photon(i,7);        % total power at the receiver plane (sum of received photons)
                %rec_dist(i) = totaldist(i)+ dist_to_rec; % individual photon's distance traveled.
                photon(i,8) = 0; % mark photon as received
                photonsRemaining = photonsRemaining-1; % decrement number of outstanding photons
                % update the total distance the photon has traveled
                totaldist(i) = totaldist(i) + dist_to_rec;
                photon(i,9)=totaldist(i);
                
                %% if photon does not pass reciver plane ********************************
            else
                photon(i,1) = photon(i,1) + x_step;         % move to new x position
                photon(i,2) = photon(i,2) + y_step;         % move to new y position
                photon(i,3) = photon(i,3) + z_step;         % move to new z position
                %                 trace_photon(i,1)=photon(i,1);
                %                 trace_photon(i,2)=photon(i,2);
                %                 trace_photon(i,3)=photon(i,3);
                %                 %x=x+1;
                
                % update the total distance the photon has traveled
                totaldist(i) = totaldist(i) + r;
                photon(i,9)=totaldist(i);
                
                %                      if strcmp(Limits_tank,'true')
                %
                %                          % Actual tank dimensions
                %                          xLimMax = 0.5*size_tank; %0.61*sizeMult in meter
                %                          xLimMin = -0.5*size_tank;%-.61*sizeMult;
                %                          yLimMax = 0.5*size_tank;%0.1905*sizeMult;
                %                          yLimMin = -0.5*size_tank;%-.5842*sizeMult;
                %                          zLimMax = receiver_z;
                %                          photon= boundries_reflection(num_photons,n1,n2,Turbulence_check,xLimMax,xLimMin,yLimMax,yLimMin,photon,i,nAir,wallAbsorption);
                %
                %                      end
                
                % If the photon collides with the back wall, terminate
                if (photon(i,3)<zLimMin)
                    photon(i,8)=-1; % mark as terminated
                    photonsRemaining = photonsRemaining-1;    % decrement outstanding photons
                    continue;
                else    % if the photon is still in the boundaries
                    photon(i,7) = photon(i,7)*prob_of_survival;     % reduce weight
                    %% rouleting..............................
                    if  photon(i,7) < min_power
                        if rand() > (rouletteConst_inv) % if the photon is randomly chosen to be terminated photon(i,8) = -1; % -1 indicates a photon was terminated, but not received
                            photon(i,8) = -1;
                            photonsRemaining = photonsRemaining-1;      % decrement outstanding photons
                           continue;
                        else % otherwise the photon gets the energey of terminated photons
                            photon(i,7) = photon(i,7)*rouletteConst;    % shift power of terminated photons to new photon
                            if photon(i,7)==0
                                photon(i,8) = -1;
                               photonsRemaining = photonsRemaining-1;
                               continue;
                            end
                        end
                    end
                    %.............................................
                    
                    if abs(photon(i,6)) > max_uz         % if uz ~ 1
                        photon(i,4) = sign(photon(i,6))*sin(theta)*cos(phi);
                        photon(i,5) = sign(photon(i,6))*sin(theta)*sin(phi);
                        photon(i,6) = sign(photon(i,6))*cos(theta);
                        %disp('mu_z near 1')
                    else
                        sqrt_uz = sqrt(1-photon(i,6)^2);
                        old_ux = photon(i,4);
                        old_uy = photon(i,5);
                        old_uz = photon(i,6);
                        photon(i,4) = (sin(theta)/sqrt_uz)*(old_ux*old_uz*cos(phi)-old_uy*sin(phi)) + old_ux*cos(theta);   % ux
                        photon(i,5) = (sin(theta)/sqrt_uz)*(old_uy*old_uz*cos(phi) + old_ux*sin(phi)) + old_uy*cos(theta);   % uy
                        photon(i,6) = (-sin(theta)*cos(phi))*sqrt_uz + old_uz*cos(theta); % uz
                    end
                    
%                      if photon(i,6)<0
%                         test1=test1+1;
%                         photon(i,8) =-1;  % mark as terminated
%                         photonsRemaining = photonsRemaining-1;    % decrement outstanding photons
%                         continue;
%                      end
                
                    % Normalize the pointing vectors -> ux^2 + uy^2 + uz^2 = 1^2
                    if abs(1-(photon(i,4)^2 + photon(i,5)^2 + photon(i,6)^2)) > 1e-11
                        normLength = sqrt(photon(i,4)^2 + photon(i,5)^2 + photon(i,6)^2);
                        photon(i,4) = photon(i,4) / normLength;
                        photon(i,5) = photon(i,5) / normLength;
                        photon(i,6) = photon(i,6) / normLength;
                        %disp('Vector normalization wrong!');
                    end
                end
                %[photon, totaldist,photonsRemaining,check_out2]=photon_not_passed_rcvplane(num_photons,n1,n2,Turbulence_check,totaldist,photon,x_step,y_step,z_step,r,Limits_tank,size_tank,receiver_z,i,nAir,wallAbsorption,zLimMin,photonsRemaining,prob_of_survival,min_power,rouletteConst_inv,rouletteConst,max_uz,theta,phi );
            end
        end
        
        %         if i<101
        %             Trace_Photon{i}(step_count,1)=photon(i,1);
        %             Trace_Photon{i}(step_count,2)=photon(i,2);
        %             Trace_Photon{i}(step_count,3)=photon(i,3);
        %         end
        %         step_count=step_count+1;
        
    end
    %     %% trace photon...
    %     %clear  trace_photon
    %     %x=2;
    %     trace_photon(:, 1:3)=photon(:, 1:3);
    %     %photon tracing....
    %     %figure
    %
    %     %         plot(trace_photon(3,:) ,trace_photon(1,:) )
    %     %         xlabel('z direction');
    %     %         ylabel('x direction ');
    %     %         %zlabel('z direction');
    %
    %     plot3(trace_photon(:,1) , trace_photon (:,2) , trace_photon (:,3) )
    %     xlabel('x direction');
    %     ylabel('y direction ');
    %     zlabel('z direction');
    %
    %     hold on
    %     pause(0.001)
    
end


% %% output data
%
% rec_loc_final = ones(total_rec_packets,5);
% j = 1;
%
% for i = 1:num_photons % iterate over all photons
%     if (photon(i,8) == 0) % if the photon was received
%         rec_loc_final(j,:)=rec_loc(i,:);       % record the receive location and angles
%         j = j +1; % increment the number of received photons
%     end
% end
%
%
% j = 1;
% total_rec_dist = zeros(total_rec_packets,1);
% rec_weights = zeros(total_rec_packets,1);
% total_rec_power=0;
%
% for i = 1:num_photons
%     if (photon(i,8) == 0)
%         total_rec_dist(j) = rec_dist(i);         % store the distance traveled by each photon
%         rec_weights(j) = photon(i,7); % store the weights of received photons
%         j = j + 1;
%     end
% end
%
%
%
% location_weight(:,1)=rec_loc_final(:,1);  %x of recieved photon
% location_weight(:,2)=rec_loc_final(:,2);  %y of recieved photon
% location_weight(:,3)=rec_weights;   %weigtht of recieved photon
%
% %[ location_weight,rec_weights,total_rec_dist,rec_loc_final]=Output_data(total_rec_packets,num_photons, photon,rec_loc,rec_dist);

%% saving data...............................................
All_Received_Photons = zeros(total_rec_packets,9);
% total_rec_dist = zeros(total_rec_packets,1);
% rec_weights = zeros(total_rec_packets,1);
Recived_indexes=find(photon(:,8) == 0);
All_Received_Photons(:,1:9)=photon(Recived_indexes,1:9);


%% at Reciver .....

rec_fov = [20 45 90 180 20 45 90 180 20 45 90 180 20 45 90 180].*pi./180;
num_fov=length(rec_fov);
rec_aperture = [ones(num_fov/4,1).*0.2; ones(num_fov/4,1).*0.4; ones(num_fov/4,1).*0.6; ones(num_fov/4,1).*0.8];
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
%     Num_bin_photon {j}= floor((Reciever_Photons{j}(:,9)-receiver_z) / Dbin)+1;       % Partition the distances past the balistic distance into bins of distBinWidth
%     Num_bin = round(Lwin/Dbin);
%     % Hist_distance= zeros(num_rx,Num_bin);  % Creat an array for the distance bins for each receiver
%     
%     
%     for In_H=1:Num_bin
%         index_for_bin=find(Num_bin_photon {j}==In_H);
%         Sum_Bin_Weight= sum(Reciever_Photons{j}(index_for_bin,7));
%         Hist_distance(In_H,j)=Sum_Bin_Weight;
%     end
    
  
      
end
end



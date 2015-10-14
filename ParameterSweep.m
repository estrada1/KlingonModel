% filename = 'BerkeleyImpact_TallBody';
% batch = 'Oct 8 Skype'; 
% trial = 'No tail';
% time = char(datetime); 
% mkdir(batch,time)

close all; clear; clc; 

resultsMap = containers.Map; 
m = 0:0.01:.1;
l = 0:0.1:1;
k = 0:0.1:1; 
for mm = m;
    for ll = l; 
       for kk = k
           
           %Denote conditions we're looking at
           disp(['m: ' num2str(mm) 'l: ' num2str(ll) 'k: ' num2str(kk)])
           
           % Run Simulation
           BerkeleyImpact_TallBody(mm, ll, kk, 0.5);
           MGworkspace = importMG('BerkeleyImpact_TallBody', 'thisSimulation'); 

           load(MGworkspace); 

           ReboundMag = sqrt(Fx_rebound.^2 + Fy_rebound.^2);
           MaxCompression = min(Fx_contact);
           
           tMaxCompression = t(find(Fx_contact == MaxCompression,1));
           MaxRebound = max(ReboundMag);
           tMaxRebound = t(find(ReboundMag == MaxRebound,1));
           mass = mm; 
           taillength = ll;
           stiffness = kk; 
           
           WStruct=ws2struct();    
          
           key = ['m' num2str(mm) ' l' num2str(ll) ' k' num2str(kk)];
           resultsMap(key) = WStruct; 
       end
    end
end

keys = resultsMap.keys; 
n = resultsMap.Count;

save('LargerSweep',resultsMap)

stop
%%
% Mass Length Stiffness MaxRebound MaxCompression 
resultsSum = zeros(n,5);

for ii = 1:size(keys,2)
    thisKey = keys{ii}
    thisStruct = resultsMap(thisKey); 
    resultsSum(ii,:) = [thisStruct.mass thisStruct.length ...
        thisStruct.stiffness thisStruct.MaxRebound thisStruct.MaxCompression]; 
end
resultsSum

%% Plot stuff 
l = 0:0.1:.3;
k = 0:0.1:.3; 

constantMass = resultsSum(:,1) == 0.010;
resultsConstantMass = resultsSum(constantMass,:);
resultsConstantMass_Rebound = reshape(resultsConstantMass(:,4),[length(l) length(k)]); 
resultsConstantMass_Compression = reshape(resultsConstantMass(:,5),[length(l) length(k)]); 

figure(1)
surf(l, k, (resultsConstantMass_Rebound-4)/14.3)
xlabel('Tail length [m]')
ylabel('Tail Stiffness [N/rad]')
zlabel('Maximum Stretch in Rebound Spring [m]')

figure(2)
scatter3(resultsSum(:,1),resultsSum(:,2),resultsSum(:,3), (resultsSum(:,4)-4)*40, resultsSum(:,5));
xlabel('Mass')
ylabel('Length')
zlabel('Stiffness')

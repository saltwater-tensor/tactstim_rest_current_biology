% HMM tactstim behavior
source_model = 'PNAI';
trial_type = 'preStim_rest';
downsample_rate = 250;

R1 = [];
R2 = [];

fg = figure(1);
fg1 = figure(2);
fg2 = figure(3);
hold on
colors = jet(size(BEHAVIOR,1));

for b = 1:size(BEHAVIOR,1)
    
    bhvr = BEHAVIOR{b,2};
    response1 = find(ismember({bhvr.Response}, {'Response1'}));
    response2 = find(ismember({bhvr.Response}, {'Response2'}));
    response1_SOA = [bhvr(response1).SOA_ms].';
    response1_SOA = response1_SOA(response1_SOA~=-1);
    response2_SOA = [bhvr(response2).SOA_ms].';
    
    response1_SOA = round(response1_SOA);
    response2_SOA = round(response2_SOA);
    
    
    total_SOA = [response1_SOA;response2_SOA];
    
     [SOA1, SOA2, SOA1_num, SOA2_num] = get_SOA(total_SOA);

     %%
    if  SOA1 > SOA2
        
        aio = SOA1+10;
        bio = SOA1-10;
        
        if bio<=SOA2 && SOA2<=aio
            %Same SOA
            if abs(SOA1-SOA2)<2 %Numerical error
                response1_SOA_crit_A = response1_SOA(response1_SOA == SOA1);
                response1_SOA_crit_B = response1_SOA(response1_SOA == SOA2);
                response1_SOA_crit = [response1_SOA_crit_A;response1_SOA_crit_B];

                response2_SOA_crit_A = response2_SOA(response2_SOA == SOA1);
                response2_SOA_crit_B = response2_SOA(response2_SOA == SOA2);
                response2_SOA_crit = [response2_SOA_crit_A;response2_SOA_crit_B];
                
                if SOA1_num > SOA2_num
                    critical_SOA(b,1) = SOA1;
                    critical_SOA(b,2) = SOA1;
                else
                    critical_SOA(b,1) = SOA2;
                    critical_SOA(b,2) = SOA2;
                end
                


            else
                
                if SOA1_num > SOA2_num
                    critical_SOA(b,1) = SOA1;
                    critical_SOA(b,2) = SOA1;
                    SOA_to_use = SOA1;
                else
                    critical_SOA(b,1) = SOA2;
                    critical_SOA(b,2) = SOA2;
                    SOA_to_use = SOA2;
                end
                
                
                response1_SOA_crit = response1_SOA(response1_SOA == SOA_to_use);
                response2_SOA_crit = response2_SOA(response2_SOA == SOA_to_use);
                critical_SOA(b,1) = SOA_to_use;
                critical_SOA(b,2) = SOA_to_use;
            end
            
        else
            % Two SOAs
                response1_SOA_crit_A = response1_SOA(response1_SOA == SOA1);
                response1_SOA_crit_B = response1_SOA(response1_SOA == SOA2);
                response1_SOA_crit = [response1_SOA_crit_A;response1_SOA_crit_B];
                
                response2_SOA_crit_A = response2_SOA(response2_SOA == SOA1);
                response2_SOA_crit_B = response2_SOA(response2_SOA == SOA2);
                response2_SOA_crit = [response2_SOA_crit_A;response2_SOA_crit_B];
                critical_SOA(b,1) = SOA1;
                critical_SOA(b,2) = SOA2;
        end
        
        
        
      elseif SOA2 > SOA1

            aio = SOA2+10;
            bio = SOA2-10;
        if bio<=SOA1 && SOA1<=aio
            %Same SOA
            if abs(SOA1-SOA2)<2 %Numerical error
                response1_SOA_crit_A = response1_SOA(response1_SOA == SOA1);
                response1_SOA_crit_B = response1_SOA(response1_SOA == SOA2);
                response1_SOA_crit = [response1_SOA_crit_A;response1_SOA_crit_B];

                response2_SOA_crit_A = response2_SOA(response2_SOA == SOA1);
                response2_SOA_crit_B = response2_SOA(response2_SOA == SOA2);
                response2_SOA_crit = [response2_SOA_crit_A;response2_SOA_crit_B];
                
                if SOA1_num > SOA2_num
                    critical_SOA(b,1) = SOA1;
                    critical_SOA(b,2) = SOA1;
                else
                    critical_SOA(b,1) = SOA2;
                    critical_SOA(b,2) = SOA2;
                end                
                

            else
                
                if SOA1_num > SOA2_num
                    critical_SOA(b,1) = SOA1;
                    critical_SOA(b,2) = SOA1;
                    SOA_to_use = SOA1;
                else
                    critical_SOA(b,1) = SOA2;
                    critical_SOA(b,2) = SOA2;
                    SOA_to_use = SOA2;
                end
                
                
                response1_SOA_crit = response1_SOA(response1_SOA == SOA_to_use);
                response2_SOA_crit = response2_SOA(response2_SOA == SOA_to_use);
                critical_SOA(b,1) = SOA_to_use;
                critical_SOA(b,2) = SOA_to_use;
            end
            
        else
            % Two SOAs
                response1_SOA_crit_A = response1_SOA(response1_SOA == SOA1);
                response1_SOA_crit_B = response1_SOA(response1_SOA == SOA2);
                response1_SOA_crit = [response1_SOA_crit_A;response1_SOA_crit_B];

                response2_SOA_crit_A = response2_SOA(response2_SOA == SOA1);
                response2_SOA_crit_B = response2_SOA(response2_SOA == SOA2);
                response2_SOA_crit = [response2_SOA_crit_A;response2_SOA_crit_B];
                critical_SOA(b,1) = SOA1;
                critical_SOA(b,2) = SOA2;
        end            
            
                
                
            
        
            
      end
%%
    
    
    R1_crit_SOA(b) = length(response1_SOA_crit);
    
    R2_crit_SOA(b) = length(response2_SOA_crit);
    
    
    plot_error_bars(mean(response1_SOA),std(response1_SOA),b,colors(b,:),fg1)
    plot_error_bars(mean(response2_SOA),std(response2_SOA),b,colors(b,:),fg2)
    
    R1 = [R1;response1_SOA];
    R2 = [R2;response2_SOA];
    
    
end

figure(4)
hold on
plot(R1_crit_SOA,'LineWidth',3)
plot(R2_crit_SOA,'LineWidth',3)
plot(R2_crit_SOA-R1_crit_SOA,'LineWidth',10)
plot(1:25,repmat(0,25,1),'LineWidth',3)
figure(5)
hold on
p2 = plot(critical_SOA(:,2),'LineWidth',10);
p2.Color(4) = 0.5;
p1 = plot(critical_SOA(:,1),'LineWidth',3,'Marker','*');

save(['critical_SOA_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_'],'critical_SOA','-v7.3')

% 
% 
% figure(6)
% hold on
% scatter(R1_crit_SOA',critical_SOA','filled')
% scatter(R2_crit_SOA',critical_SOA','filled')

function plot_error_bars(data_mean,data_error,positions,colrs,fig_handle)
figure(fig_handle)
hold on
for pos = 1:1:length(data_mean)
    
    e1 = errorbar(positions(pos),data_mean(pos),data_error(pos),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8,'Color',colrs);
    
end

hold off
end

function [SOA1, SOA2, SOA1_num, SOA2_num] = get_SOA(t_SOA)

    
    [n,bin] = hist(t_SOA,unique(t_SOA));
    [~,idx] = sort(-n);
    nums = n(idx);
    SOAs = bin(idx);
    SOA1 = SOAs(1);
    SOA2 = SOAs(2);
    SOA1_num = nums(1);
    SOA2_num = nums(2);

    
    
end
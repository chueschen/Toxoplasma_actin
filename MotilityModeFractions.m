%%% Analyzing frequency of occurence of different Toxoplasma motility modes for different concentrations 
% the actin-stabilizer jasplakinolide %%%

%% Import data

conc = [ 0, 10, 25, 50, 100, 500]; % concentartions of jasplakinolide (jasp) used
% Create an array with... columns = PATCH, helical, circular, twirling, total; rows = each
% concentration of jasp = 0 nM, 10 nM, 25 nM, 50nM, 100nM, 500nM. 

% Manually important counts of motility events for each condition. Some
% day, it would be fantastic to develop an accurate AI-enriched pipeline
% for automatic identification and annotation from movies, but
% unfortunately, we found that manual annotation was more accurate by far
% than other strategies we tried. We blindfolded the data and painstakingly
% went through every movie (blind to the experiment # and condition), 
% annotating all motility events, and we imported those data here.

%% Calculate fractions

% fractions of each motility type
dataFrac = data(1:4,:) ./ data(5,:); % each mode over total gliding events

% also create an array for just unidirectional vs. patch
dataFracPvsU(1,:) = data(1,:) ./ data(5,:);
dataFracPvsU(2,:) = sum(data(2:4,:)) ./ data(5,:);

% Calculate fractions of h, c, p of just unidirectional events
dataFracUni = data(2:4,:) ./ sum(data(2:4,:));

%% Calculate uncertainty

% Calculate weighted standard deviation with bias correction for "reliability weights"
% for example: https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf

dataByExp = cat(3,CH095, CH096, CH097); %concatenate data

% First, for all 4 gliding modes
dataFracSTD = NaN(4,6);
% Keep it simple and just loop through, calculating weight std as defined
% in link above
for i = 1:4
    for j = 1:6
        dataFracSTD(i,j) = sqrt( (sum(dataByExp(5,j,:).*(dataByExp(i,j,:)./dataByExp(5,j,:) - dataFrac(i,j)).^2, 3) )./...
           ( sum(dataByExp(5,j,:),3) - (sum(dataByExp(5,j,:).^2,3)./sum(dataByExp(5,j,:),3)) ) );
    end
end

%% Second, for just patch vs. uni
dataFracPvsUSTD = NaN(2,6);
tempUni1 = dataByExp(2:4,:,:); % creating a 6th row in dataByExp that has a total of uni
tempUni2 = sum(tempUni1,1);
dataByExp(6,:,:) = tempUni2; 
clear tempUni1; clear tempUni2

for i = 1 % patch
    for j = 1:6
        dataFracPvsUSTD(1,j) = sqrt( (sum(dataByExp(5,j,:).*(dataByExp(i,j,:)./dataByExp(5,j,:) - dataFracPvsU(1,j)).^2, 3) )./...
           ( sum(dataByExp(5,j,:),3) - (sum(dataByExp(5,j,:).^2,3)./sum(dataByExp(5,j,:),3)) ) );
    end
end
for i = 6 % uni
    for j = 1:6
        dataFracPvsUSTD(2,j) = sqrt( (sum(dataByExp(5,j,:).*(dataByExp(i,j,:)./dataByExp(5,j,:) - dataFracPvsU(2,j)).^2, 3) )./...
           ( sum(dataByExp(5,j,:),3) - (sum(dataByExp(5,j,:).^2,3)./sum(dataByExp(5,j,:),3)) ) );
    end
end

%% Third, for h vs. c vs. twirling of the unidirectional events
dataFracUniSTD = NaN(3,6);
for i = 2:4
    for j = 1:6
        dataFracUniSTD(i-1,j) = sqrt( (sum(dataByExp(6,j,:).*(dataByExp(i,j,:)./dataByExp(6,j,:) - dataFracUni(i-1,j)).^2, 3) )./...
           ( sum(dataByExp(6,j,:),3) - (sum(dataByExp(6,j,:).^2,3)./sum(dataByExp(6,j,:),3)) ) );
    end
end


%% Stacked bar graph -- just unidirectional vs patch
figure
bar(dataFracPvsU','stacked'); % bar(conc,dataFrac','stacked');
hold on
errorbar(cumsum(dataFracPvsU)',dataFracPvsUSTD','.k');
%set(gca,'xticklabel', {'0 nM', '10 nM', '25 nM', '50 nM', '100 nM', '500 nM'})
set(gca,'xticklabel', conc)
xlabel('jasplakinolide (nM)')
ylabel('fraction of motility events')
legend('patch','unidirectional')

%% Stacked bar graph of fraction h, c, t of just unidirectional cells

figure
bar(dataFracUni','stacked');
hold on
errorbar(cumsum(dataFracUni)', dataFracUniSTD','.k');
set(gca,'xticklabel', conc)
xlabel('jasplakinolide (nM)')
ylabel('fraction of unidirectional motility events')
legend('helical', 'circular', 'twirling')

% Maybe shouldn't include 500 nM condition. Just too few cells.
figure
bar(dataFracUni(:,1:5)','stacked');
hold on
errorbar(cumsum(dataFracUni)', dataFracUniSTD','.k');
set(gca,'xticklabel', conc(1:5))
xlabel('jasplakinolide (nM)')
ylabel('fraction of unidirectional motility events')
legend('helical', 'circular', 'twirling')

%% Plotting -- Stacked bar graph with error bars for all conditions combined
figure
bar(dataFrac','stacked'); % bar(conc,dataFrac','stacked');
hold on
errorbar(cumsum(dataFrac)', dataFracSTD','.k');
set(gca,'xticklabel', conc)
xlabel('jasplakinolide (nM)')
ylabel('fraction of motility events')
legend('patch', 'helical', 'circular', 'twirling')

%% One-way ANOVA % https://www.mathworks.com/help/stats/one-way-anova.html

% Create array with fractions by experiment
dataByExpFrac = dataByExp(1:4,:,:)./dataByExp(5,:,:);
for i = 1 % just for patch fraction
  for j = 1:6
    for k = 1:3
      forANOVApatch(k,j) = dataByExpFrac(i,j,k);
    end
  end
end
[p,tbl,stats] = anova1(forANOVApatch); % this p value is about null hypothesis that ALL groups are same
% Check out multiple comparisons
results = multcompare(stats);

dataByExpFracUni = dataByExp(2:4,:,:)./dataByExp(6,:,:); %6th row is uni totals
for i = 1 % just for helical fraction of uni
  for j = 1:5 % - NOTE: EXCLUDING 500 nM as explained in manuscript. Too few observations 
                % here to say anything useful because at high jasp, almost all cells patch glide.
    for k = 1:3
      forANOVAUni_h(k,j) = dataByExpFracUni(i,j,k);
    end
  end
end
[pUni_h,tblUni_h,statsUni_h] = anova1(forANOVAUni_h);
% Check out multiple comparisons
resultsUni_h = multcompare(statsUni_h);

for i = 2 % just for circular fraction of uni
  for j = 1:5 % 
    for k = 1:3
      forANOVAUni_c(k,j) = dataByExpFracUni(i,j,k);
    end
  end
end
[pUni_c,tblUni_c,statsUni_c] = anova1(forANOVAUni_c);
% Check out multiple comparisons
resultsUni_c = multcompare(statsUni_c);
%
for i = 3 % just for twirling fraction of uni
  for j = 1:5
    for k = 1:3
      forANOVAUni_t(k,j) = dataByExpFracUni(i,j,k);
    end
  end
end
[pUni_t, tblUni_t,statsUni_t] = anova1(forANOVAUni_t);
% Check out multiple comparisons
resultsUni_t = multcompare(statsUni_t);


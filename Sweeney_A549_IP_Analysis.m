% Emmott: A549 IFIT pulldowns

clear
clc

% Experimental design:
%Experiment 1: 15 samples TMTpro:
%3 each of:

%Bead+ Ab WT
%Bead+ Ab KO (control)
%Bead+ Ab WT + MG132 (highest expressing)
%Bead+ Ab KO + MG132 (control)
%Bead Only



% Import data

path = '/Users/ed/Dropbox/Liverpool/Collaborations/Trevor_IPs/'

dat = struct();
dat.prot = readtable([path , 'txt_A549/proteinGroups.txt']);

dat.tmt = readmatrix([path , 'Sweeney_A549TMTperm.csv']);
%% Correct for randomisation
% Reorder TMT corrected RI by randomisation strategy. (Note 15 channels -
% no TMTpro 126C
dat.prot(:,[19:33]) = dat.prot(:,dat.tmt + 18);

%% Data cleanup

% Remove: Rev, Con
dat.prot.PotentialContaminant = categorical(dat.prot.PotentialContaminant);
dat.prot.Reverse              = categorical(dat.prot.Reverse);

dat.prot = dat.prot(dat.prot.PotentialContaminant ~= '+',:);
dat.prot = dat.prot(dat.prot.Reverse ~= '+',:);

% Extract TMT data
dat.mat = table2array(dat.prot(:,19:33));

% Convert 0 to NaN
dat.mat(dat.mat == 0) = NaN;

% Filter on rows with at least 6 valid values
logNaN = sum(isnan(dat.mat),2) >= 7;

dat.mat  = dat.mat(~logNaN , :);
dat.prot = dat.prot(~logNaN , :);
%%
% Median normalise
dat.mat = dat.mat ./ nanmedian(dat.mat);
%%
% Knnimpute
dat.mat = knnimpute(dat.mat);

%%
% Row normalise against mean of the beads only + WT IP samples
%dat.mat = dat.mat ./ mean(dat.mat , 2);
dat.mat = dat.mat ./ mean(dat.mat(:,[13:15]),2);
% Should give a value centred around log2(0).

dat.l2mat = log2(dat.mat);

%% Statistical analysis
% Perform t-test for significance
comparisons = {'wt','ko','wt_mg132','ko_mg132'};
compNum = [1,2,3;4,5,6;7,8,9;10,11,12];

for ii = 1:4
%t-test
[~ , dat.p.(comparisons{ii}) , ~ , ~] = ttest2(log2(dat.mat(:,[compNum(ii,:)])) , log2(dat.mat(:,[13:15])),'Dim',2);

%Correct for multiple hypothesis testing
[dat.fdr.(comparisons{ii}) , dat.q.(comparisons{ii})] = mafdr(dat.p.(comparisons{ii}));

% Calculate mean fold-chage per sample set
dat.mean.(comparisons{ii}) = log2(mean(dat.mat(:,compNum(ii,:)),2) ./ mean(dat.mat(:,[13:15]) , 2));
end
%%
% Two further comparisons:
% WT vs KO
[~ , dat.p.wt_vs_ko , ~ , ~] = ttest2(log2(dat.mat(:,[1:3])) , log2(dat.mat(:,[4:6])),'Dim',2);
[dat.fdr.wt_vs_ko , dat.q.wt_vs_ko] = mafdr(dat.p.wt_vs_ko);
dat.mean.wt_vs_ko = log2(mean(dat.mat(:,[1:3]),2) ./ mean(dat.mat(:,[4:6]) , 2));

% WT MG132 vs KO MG132
[~ , dat.p.wtmg_vs_komg , ~ , ~] = ttest2(log2(dat.mat(:,[7:9])) , log2(dat.mat(:,[10:12])),'Dim',2);
[dat.fdr.wtmg_vs_komg , dat.q.wtmg_vs_komg] = mafdr(dat.p.wtmg_vs_komg);
dat.mean.wtmg_vs_komg = log2(mean(dat.mat(:,[7:9]),2) ./ mean(dat.mat(:,[10:12]) , 2));

comparisons{5} = 'wt_vs_ko';
comparisons{6} = 'wtmg_vs_komg';

%% Plot
compTitles = {'WT vs. Beads-only','KO vs. Beads-only','WT MG132 vs. Beads-only'...
    'KO MG132 vs. Beads-only','WT vs. KO','WT MG132 vs. KO MG132'};

f = figure
t = tiledlayout(2,3)
for ii = 1:numel(comparisons)
    
   % filter results on minimum log2 fold-change and q-value
   dat.pval.(comparisons{ii}) = dat.p.(comparisons{ii}) <= 0.05
   %dat.mn.(comparisons{ii}) = dat.mean.(comparisons{ii}) >= 1  % Technically this is the correct way around. However, the option below is for looking at KO vs. WT.

   dat.mn.(comparisons{ii}) = dat.mean.(comparisons{ii}) >= 1 | dat.mean.(comparisons{ii}) <= -1
   %   dat.mn.(comparisons{ii}) = dat.mean.(comparisons{ii}) <= -1
   dat.both.(comparisons{ii}) = dat.pval.(comparisons{ii}) + dat.mn.(comparisons{ii}) == 2; 
    
   nexttile
   scatter(dat.mean.(comparisons{ii}) , -log10(dat.q.(comparisons{ii})),'filled','k','MarkerFaceAlpha',0.1);
   hold on
   scatter(dat.mean.(comparisons{ii})(dat.both.(comparisons{ii})) , -log10(dat.q.(comparisons{ii})(dat.both.(comparisons{ii}))),'filled','r','MarkerFaceAlpha',0.3);
   
   title(compTitles{ii})
   xlabel('Mean log_2 Sample over Mock')
   ylabel('-log_1_0 Q-value')
   xlim([-6,6]);
   ylim([0,3.5]);
   
hold off
end
title(t,'A549 IFIT IP analysis 3')


%% Figure and Table export
saveas(gcf,[path,'TS_A549_Figure3.pdf']);

for ii = 1:numel(compTitles)
    tab = table()
    tab.protGroups = dat.prot.ProteinIDs(dat.both.(comparisons{ii}));
    tab.mean = dat.mean.(comparisons{ii})(dat.both.(comparisons{ii}));
    tab.q = dat.q.(comparisons{ii})(dat.both.(comparisons{ii}));
    [b,i] = sort(tab.mean,'descend');
    tab = tab(i,:);
    
    writetable(tab,[path,'A549_hits_',comparisons{ii},'.csv']);
    clear tab
end


%%
pcntV = 10; % Goal is to perform PCA on 10% most variable proteins

    %Select calculate variance for each peptide
    dat.matvar  = var(dat.mat , 0 , 2);
    
    % Sort on variance: return index
    [~ , dat.mvIA]  = sort(dat.matvar , 'descend');
    % Grab index for top X% values (defined by pcntV)
    dat.varIa  = dat.mvIA(1:round(numel(dat.mvIA) / pcntV),:);

    % Perform PCA
    [dat.pca,~,~,~,dat.pcaV,~] = pca(log2(dat.mat(dat.varIa , :)))
    
    grpVar = categorical({'WT','WT','WT','KO','KO','KO','WT_MG132','WT_MG132','WT_MG132'...
        'YL_MG132','YL_MG132','YL_MG132','beads','beads','beads'});

    figure
    pcA = 1;
    pcB = 2;
    hgA = gscatter(dat.pca(:,pcA) , dat.pca(:,pcB) , grpVar);
    xlabel(['PC ',num2str(pcA),' (',num2str(round(dat.pcaV(pcA))),'%)']);
    ylabel(['PC ',num2str(pcB),' (',num2str(round(dat.pcaV(pcB))),'%)']);
    title('A549 IFIT1 IPs - PCA')
    
    saveas(gcf,[path,'TS_A549_PCA_Figure.pdf']);
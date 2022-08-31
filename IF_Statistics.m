
% Revise: Figure out how to weigh ring intensity by ring pixel count.


%% Select conditions for data to exclude

exclude_Treatment   = {'iso'};
exclude_Region      = {'Pia'};
exclude_TissueID    = {'010b','014b','018b','021b'};
exclude_Hemisphere  = {};

short_stats_tab = working_tab(...
    ~ismember(working_tab.Treatment,exclude_Treatment)...
    &~ismember(working_tab.Region_name,exclude_Region)...
    &~ismember(working_tab.Tissue_ID,exclude_TissueID)...
    &~ismember(working_tab.Hemisphere,exclude_Hemisphere)...
    ,:);

% long_stats_tab has a row for each image down to the level of 
long_stats_tab = an_tab(...
    ~ismember(an_tab.Treatment,exclude_Treatment)...
    &~ismember(an_tab.Region_name,exclude_Region)...
    &~ismember(an_tab.Tissue_ID,exclude_TissueID)...
    &~ismember(an_tab.Hemisphere,exclude_Hemisphere)...
    ,:);

% Remove rows that correspond to an expansion of rings that went outside
% the image
rows_to_delete = ~long_stats_tab.Rings_exist | long_stats_tab.Rings_exist & long_stats_tab.Rings_n_pix == 0;
long_stats_tab(rows_to_delete,:) = [];

%%
figure("Name","Junk");
cond = working_tab.Hemisphere == "left" &...
    working_tab.Treatment == "ab";
boxchart(working_tab.Region_name(cond),working_tab.Red_intensity_tot(cond),"GroupByColor",working_tab.Group(cond))

%% Comparison Groups
list_regions = unique(working_tab.Region_name);
cond = working_tab.Treatment == "ab";
figure("Name","Groups by region")
X = working_tab.Group;
Y = working_tab.Red_intensity_tot;
for i=1:length(list_regions)
    regcond = working_tab.Region_name == list_regions(i);
    subplot(3,3,i)
    boxchart(X(cond & regcond),Y(cond & regcond),"JitterOutliers","on","MarkerStyle",".")
    p = anova1(Y(cond & regcond),X(cond & regcond),"off");
    title(list_regions(i),["p=",p])
    y = max([quantile(Y(cond & regcond & X == "depression"),[0.25,0.5,0.75]);...
        quantile(Y(cond & regcond & X == "control"),[0.25,0.5,0.75])]);
    ymax = y(2)+(y(3)-y(1))*1.75;
    ylim([0 ymax])
end

%% junk
cond = working_tab.Hemisphere == "left" &...
    working_tab.Treatment == "ab";
figure;
subplot(1,2,1)
gscatter(working_tab.Capillary_area_ratio(cond),working_tab.Red_intensity_tot(cond),working_tab.Group(cond))
subplot(1,2,2)
gscatter(working_tab.Capillary_area_ratio(cond),working_tab.Red_intensity_tot(cond),working_tab.Region_name(cond))


%%
cond = ...
    working_tab.Group == "depression" &...
    working_tab.Region_name == "Subi Ento" &...
    working_tab.Treatment == "ab";

p = anova1(working_tab.Green_intensity_tot(cond),working_tab.Hemisphere(cond));

%% Do repeated measures anova with

% ranova or manova?




%% Get one data point per region per brain, then make a few graphs, excluding redone sections

sst = grpstats(short_stats_tab,{'Group','Patient','Region_name'},{'mean','meanci','sem','std','numel'},"DataVars",{'Nuclei_count','Red_per_nucleus','Capillary_area'});
list_regions = unique(sst.Region_name);
mean_table = grpstats(sst,{'Group','Region_name'},"mean","DataVars",{'mean_Capillary_area','mean_Red_per_nucleus','mean_Nuclei_count'});

figure("Name","Comparison of mean capillary ratio between groups")
X       = sst.Group;
Y       = sst.mean_Capillary_area;
Xmean   = mean_table.Group;
Ymean   = mean_table.mean_mean_Capillary_area;
for i=1:length(list_regions)
    regcond = sst.Region_name == list_regions(i);
    mean_regcond = mean_table.Region_name == list_regions(i);
    subplot(2,4,i)
    scatter(X(regcond),Y(regcond),20,"black","o",'jitter','on')
    [h,p] = ttest2(Y(regcond & X=="depression"),Y(regcond & X=="control"));
    hold on;
    boxchart(Xmean(mean_regcond),Ymean(mean_regcond))
    title(list_regions(i),["p=",p])
    xlabel("Group")
    ylabel("Capillary area ratio, means within patient")
    ylim([0 max(Y)*1.1])
    hold off;
end

figure("Name","Comparison of mean total intensity between groups")
X       = sst.Group;
Y       = sst.mean_Red_per_nucleus;
Xmean   = mean_table.Group;
Ymean   = mean_table.mean_mean_Red_per_nucleus;
for i=1:length(list_regions)
    regcond = sst.Region_name == list_regions(i);
    mean_regcond = mean_table.Region_name == list_regions(i);
    subplot(2,4,i)
    scatter(X(regcond),Y(regcond),20,"black",'jitter','on')
    [h,p] = ttest2(Y(regcond & X=="depression"),Y(regcond & X=="control"));
    hold on;
    boxchart(Xmean(mean_regcond),Ymean(mean_regcond))
    title(list_regions(i),["p=",p])
    xlabel("Group")
    ylabel("Total red intensity, means within patient")
    ylim([0 max(Y)*1.1])
    hold off;
end

figure("Name","Comparison of mean nucleus count between groups")
X       = sst.Group;
Y       = sst.mean_Nuclei_count;
Xmean   = mean_table.Group;
Ymean   = mean_table.mean_mean_Nuclei_count;
for i=1:length(list_regions)
    regcond = sst.Region_name == list_regions(i);
    mean_regcond = mean_table.Region_name == list_regions(i);
    subplot(2,4,i)
    scatter(X(regcond),Y(regcond),20,"black",'jitter','on')
    [h,p] = ttest2(Y(regcond & X=="depression"),Y(regcond & X=="control"));
    hold on;
    boxchart(Xmean(mean_regcond),Ymean(mean_regcond))
    title(list_regions(i),["p=",p])
    xlabel("Group")
    ylabel("Nucleus count, means within patient")
    ylim([0 max(Y)*1.1])
    hold off;
end


%% Intensity-distance statistics by patient

sst = grpstats(long_stats_tab,{'Group','Patient','Region_name','Regions_expansion_radii_microns'},{'mean','meanci','sem','std','numel'},"DataVars",'Rings_int_per_pixel');
list_regions = unique(sst.Region_name);
depression  = sst.Group == "depression";
control     = sst.Group == "control";

figure("Name","Intensity-distance plots by patient")
X = sst.Regions_expansion_radii_microns;
Y = sst.mean_Rings_int_per_pixel;
G = sst.Group;
for i=1:length(list_regions)
    regcond = sst.Region_name == list_regions(i);
    subplot(2,4,i)
    gscatter(X(regcond),Y(regcond),G(regcond),'br','ox',5,'on')
    title(list_regions(i))
    xlabel("Distance from capillaries, µm")
    ylabel("Mean pixel intensity")
    ylim([0 inf])
    legend()
end


%% Intensity-distance statistics by group
% Revise: Should the number reflect mean and std per image or per patient?
% middle_step makes std larger by first making group stats within patient,
% and from the resulting means sst then calculates the means of the
% patient means.
middle_step = grpstats(long_stats_tab,{'Group','Region_name','Patient','Regions_expansion_radii_microns'},{'mean','meanci','sem','std','numel'},"DataVars",'Rings_int_per_pixel');
sst = grpstats(middle_step,{'Group','Region_name','Regions_expansion_radii_microns'},{'mean','meanci','sem','std','numel'},"DataVars",'mean_Rings_int_per_pixel');
list_regions = unique(sst.Region_name);
depression  = sst.Group == "depression";
control     = sst.Group == "control";

figure("Name","Intensity-distance plots by group")
X = sst.Regions_expansion_radii_microns;
Y = sst.mean_mean_Rings_int_per_pixel;
G = sst.Group;
Err = sst.sem_mean_Rings_int_per_pixel;
for i=1:length(list_regions)
    regcond = sst.Region_name == list_regions(i);
    subplot(2,4,i)
%     gscatter(X(regcond),Y(regcond),G(regcond))
    for j = 1:2
        if j==1
            cond = regcond & depression;
            errorbar(X(cond),Y(cond),Err(cond),"LineStyle","--")
            hold on;
        else
            cond = regcond & control;
            errorbar(X(cond),Y(cond),Err(cond),"LineStyle","-")
            hold off;
        end
    end
    title(list_regions(i))
    xlabel("Distance from capillaries, µm")
    ylabel("Mean pixel intensity")
    ylim([0 inf])
    legend("depression","control")
end


%% Try out single graphs

% stats_tab = working_tab(working_tab.Treatment=="ab"&working_tab.Region_name~="Pia"&working_tab.Tissue_ID~="010b"&working_tab.Tissue_ID~="014b"&working_tab.Tissue_ID~="018b"&working_tab.Tissue_ID~="021b",:);
% sst = grpstats(stats_tab,{'Group','Patient','Region_name'},{'mean','meanci','sem','std'},"DataVars",{'Nuclei_count','Red_intensity_tot','Capillary_area_ratio'});
% Grp_Reg = sst.Region_name .*  sst.Group;
% sst.Grp_Reg = Grp_Reg;
% mean_table = grpstats(sst,{'Group','Region_name'},{'mean'},"DataVars",{'mean_Capillary_area_ratio','mean_Red_intensity_tot','mean_Nuclei_count'});
% 
% % 
% % scatter([sst.Group,sst.Region_name],Y)
% errorbar(sst.Grp_Reg,sst.mean_Red_intensity_tot,sst.sem_Red_intensity_tot,'bo',"Marker","_","LineStyle","none","CapSize",3,LineWidth=2)
% % boxchart(sst.Group,sst.mean_Red_intensity_tot)
% 
% boxchart(sst.Grp_Reg,sst.mean_Red_intensity_tot)
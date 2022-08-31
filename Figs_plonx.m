function [an_tab,working_tab] = Figs_plonx(dat_tab,varargin)

%% Resolve optional arguments

for funinputs = 1:nargin-1
    if isstring(varargin{funinputs})
        if varargin{funinputs} == "examineID" && funinputs ~= nargin
            lookupID = varargin{funinputs+1};
        end
    end
end

if ~exist("lookupID","var")
    lookupID = randsample(dat_tab.Tissue_ID,1);
end

%% Get some god damned easy matrices

working_tab = dat_tab;

% Correct irregular tissue names
working_tab.Region_name(working_tab.Region_name == "Interpyramidal") = "Interneuronal";
working_tab.Tissue_ID(working_tab.Tissue_ID == "14b") = "014b";
working_tab.Tissue_ID(working_tab.Tissue_ID == "10b") = "010b";
% Clears unused categories
working_tab.Tissue_ID       = removecats(working_tab.Tissue_ID);
working_tab.Region_name     = removecats(working_tab.Region_name);
working_tab.Project_name    = removecats(working_tab.Project_name);
working_tab.Section_no      = removecats(working_tab.Section_no);
working_tab.Frame_no        = removecats(working_tab.Frame_no);
working_tab.Treatment       = removecats(working_tab.Treatment);
% Gets the number of pixels in each ring by subtracting the number of pixels in the subsequent region from the previous region
% Subtract lower ranked regions around capillaries from the rank one step
% above to achieve rings around the capillaries.
working_ringspix = [working_tab.Regions_pixels(:,1),working_tab.Regions_pixels(:,2:end)-working_tab.Regions_pixels(:,1:end-1)];
working_ringspix(working_ringspix < 0 | isnan(working_ringspix)) = 0;
working_tab.Rings_n_pix = working_ringspix;
% Gets the intensity per pixel in each ring
working_ringsint = [working_tab.Regions_intensity(:,1),working_tab.Regions_intensity(:,2:end)-working_tab.Regions_intensity(:,1:end-1)];
working_ringsint_pp = working_ringsint ./working_ringspix; % Revise: When taking the mean across all samples, whould pixel number be taken into account so that different samples are weighted differently?
working_ringsint_pp(working_ringsint_pp < 0 | isnan(working_ringsint_pp)) = 0;
working_tab.Rings_int_per_pixel = working_ringsint_pp;
% REVISE: Why does the ring intensity act weird when going into high
% expansions?
rings_exist = sum(working_ringsint_pp,2) > 0;
working_tab.Rings_exist = rings_exist;
clear working_ringsint
clear working_ringspix
Capillary_area_ratio = working_tab.Capillary_area ./ (working_tab.Total_pixels ./ (3.9^2));
working_tab = addvars(working_tab,Capillary_area_ratio);
clear Capillary_area_ratio
working_tab.Red_per_nucleus = working_tab.Red_intensity_tot ./ working_tab.Nuclei_count;

%% Load unblinding data
working_tab.Group(:)         = categorical("N/A");
working_tab.Hemisphere(:)    = categorical("N/A");
working_tab.Patient(:)       = categorical("N/A");

unblinding = open("test_unblinding.mat");
unblinding = unblinding.test_unblinding;

list_TissueID = unique(working_tab.Tissue_ID);
for i = 1:length(list_TissueID)
    current_ID = list_TissueID(i);
    working_ID = working_tab.Tissue_ID == list_TissueID(i);
    unblinding_data = find(unblinding.tisID==current_ID);
    working_tab.Group(working_ID) = unblinding.grpDK(unblinding_data);
    working_tab.Hemisphere(working_ID) = unblinding.hemis(unblinding_data);
    working_tab.Patient(working_ID) = unblinding.patients(unblinding_data);
end

working_tab.Group           = removecats(working_tab.Group);
working_tab.Hemisphere      = removecats(working_tab.Hemisphere);
working_tab.Patient         = removecats(working_tab.Patient);


%% Make the anaylsis table, much easier to work with

% addvars(an_tab,zeros(1,size(working_tab,2)),'NewVariableNames',string(working_tab.Properties.VariableNames))

variables = working_tab.Properties.VariableNames;

longest_row = 0;
for i = 1:length(variables)
    if size(working_tab.(variables{i})(1,:),2) > longest_row
        longest_row = size(working_tab.(variables{i})(1,:),2);
    end
end

an_tab = table();
for i = 1:length(variables)
    trans_data = working_tab.(variables{i});
    if size(trans_data,2) == 1
        an_tab.(variables{i}) = repmat(working_tab.(variables{i}),[longest_row,1]);
    else
        an_tab.(variables{i}) = reshape(working_tab.(variables{i}),[numel(working_tab.(variables{i})),1]);
    end
end
an_tab = removevars(an_tab,["name","bytes"]);

%% Make me some pretty figures

ab          = an_tab.Treatment == "ab";
iso         = an_tab.Treatment == "iso";
condition   = an_tab.Region_name    == "Hilus"; ...
%             & an_tab.Section_no     == "2" ...
%             & an_tab.Tissue_ID      == "014" ...
%             & an_tab.Frame_no       == "007";

figure("Name","Intensity by distance from capillaries");
boxchart(an_tab.Regions_expansion_radii_microns(condition & ab),an_tab.Rings_int_per_pixel(condition & ab))
ylim([0 inf])
xlim([0 inf])
hold on;
boxchart(an_tab.Regions_expansion_radii_microns(condition & iso),an_tab.Rings_int_per_pixel(condition & iso))
xlabel("Distance from capillaries, µm")
ylabel("Average intensity of pixels in each ring, 8bit")

figure("Name","Broad Summary");
subplot(2,2,1)
scatter(an_tab.datenum(condition & ab),log10(an_tab.Red_intensity_tot(condition & ab)),'.','blue')
hold on;
scatter(an_tab.datenum(condition & iso),log10(an_tab.Red_intensity_tot(condition & iso)),'.','red')
xlabel("Time of imaging")
ylabel("Total red intensity of images, log10")
subplot(2,2,2)
boxchart(an_tab.Tissue_ID(condition & ab),log10(an_tab.Capillary_area(condition & ab)),"BoxFaceColor","blue","MarkerStyle",".","MarkerColor","blue")
hold on;
boxchart(an_tab.Tissue_ID(condition & iso),log10(an_tab.Capillary_area(condition & iso)),"BoxFaceColor","red","MarkerStyle",".","MarkerColor","red")
xlabel("Tissue ID")
ylabel("Capillary area, log10")

% find(dat_tab.Capillary_area == max(dat_tab.Capillary_area))

%% Make a figure of the specified tissue ID, median intensity by region and treatment

statssummary = grpstats(an_tab,{'Tissue_ID','Treatment','Region_name','Regions_expansion_radii_microns'},{'mean','median','gname'},'DataVars',{'Capillary_area','Rings_int_per_pixel'});
figure;
correctID = categorical(statssummary.gname(:,1)) == lookupID;
for regcond = 1:length(unique(categorical(statssummary.gname(correctID,3))))
    reg_n = unique(categorical(statssummary.gname(correctID,3)));
    reg_n = reg_n(regcond);
    correctreg = categorical(statssummary.gname(:,3)) == reg_n;
    for treat = 1:2
        if treat == 1
            correcttreat  = categorical(statssummary.gname(:,2)) == "ab";
            condition = correcttreat & correctreg & correctID;
            plot(statssummary.Regions_expansion_radii_microns(condition),statssummary.median_Rings_int_per_pixel(condition));
            hold on;
        else
%             correcttreat = categorical(statssummary.gname(:,2)) == "iso";
%             condition = correcttreat & correctreg & correctID;
%             plot(statssummary.Regions_expansion_radii_microns(condition),statssummary.median_Rings_int_per_pixel(condition));
%             hold on;
        end
    end
end
ylim([0 inf])
xlabel("Distance from capillary, µm")
ylabel("Intensity per pixel")
title("Median values for each region",["Tissue ID: ", lookupID])
legend(char(unique(statssummary.gname(correctID,3))))

%% Maybe useful?

% statssummary = grpstats(an_tab,{'Treatment','Region_name'},{'mean','median','gname'},'DataVars',{'Capillary_area','Red_intensity_tot'} );

figure;
    cond = an_tab.Treatment == "ab" & an_tab.Tissue_ID == lookupID;
    grpstats(an_tab.Red_intensity_tot(cond),an_tab.Region_name(cond),0.05);
    hold on;
    cond = an_tab.Treatment == "iso";
    grpstats(an_tab.Red_intensity_tot(cond),an_tab.Region_name(cond),0.05);
legend("ab","iso")
title("Mean red intensity total and 95 % confidence intervals",["Tissue ID: ",lookupID])

%% Summary graphs comparing depressed to control patients

figure("Name","Depression")
    cond = working_tab.Treatment == "ab" & working_tab.Region_name ~= "Pia";
subplot(2,2,1)
boxchart(working_tab.Group,working_tab.Capillary_area_ratio,"JitterOutliers","on","MarkerStyle",".","MarkerSize",3)
ylabel("Ratio of capillary area to total area")
ylim([0,0.1])
subplot(2,2,2)
boxchart(working_tab.Group,working_tab.Red_intensity_tot,"JitterOutliers","on","MarkerStyle",".","MarkerSize",3)
ylabel("Total red intensity")
ylim([0,10^7])
subplot(2,2,3)
boxchart(working_tab.Hemisphere,working_tab.Capillary_area_ratio,"JitterOutliers","on","MarkerStyle",".","MarkerSize",3)
ylabel("Ratio of capillary area to total area")
ylim([0,0.1])
subplot(2,2,4)
boxchart(working_tab.Hemisphere,working_tab.Red_intensity_tot,"JitterOutliers","on","MarkerStyle",".","MarkerSize",3)
ylabel("Total red intensity")
ylim([0,10^7])

%% Intensity-distance plots by group (depression vs control)
% This looks good
figure("Name","Intensity-distance plots for regions and groups, with standard deviation")
statssummary = grpstats(an_tab,{'Group','Treatment','Region_name','Regions_expansion_radii_microns'},{'mean','median','std'},'DataVars',{'Capillary_area','Rings_int_per_pixel'});
ab          = statssummary.Treatment == "ab";
dep         = statssummary.Group == "depression";
control     = statssummary.Group == "control";
list_regions = unique(statssummary.Region_name);
for i=1:length(list_regions)
    for j=1:2
        regcond = statssummary.Region_name == list_regions(i);
        if j ==1
            cond = ab & dep & regcond;
            subplot(3,3,i)
            errorbar(statssummary.Regions_expansion_radii_microns(cond),statssummary.mean_Rings_int_per_pixel(cond),statssummary.std_Rings_int_per_pixel(cond))
            title(list_regions(i))
            hold on;
        else
            cond = ab & control & regcond;
            errorbar(statssummary.Regions_expansion_radii_microns(cond),statssummary.mean_Rings_int_per_pixel(cond),statssummary.std_Rings_int_per_pixel(cond))
            ylim([0 inf])
            xlim([0 inf])
            xlabel("Distance from capillaries, µm")
            ylabel("Mean intensity of pixels in ring, 8bit")
            legend(["depression","control"])
            hold off;
        end
    end
end


%% Make a figure of intensity by distance, per region, per group
% Revise
statssummary = grpstats(an_tab,{'Tissue_ID','Treatment','Region_name','Regions_expansion_radii_microns'},{'mean','median','gname'},'DataVars',{'Capillary_area','Rings_int_per_pixel'});
figure;
correctID = categorical(statssummary.gname(:,1)) == lookupID;
for regcond = 1:length(unique(categorical(statssummary.gname(correctID,3))))
    reg_n = unique(categorical(statssummary.gname(correctID,3)));
    reg_n = reg_n(regcond);
    correctreg = categorical(statssummary.gname(:,3)) == reg_n;
    for treat = 1:2
        if treat == 1
            correcttreat  = categorical(statssummary.gname(:,2)) == "ab";
            condition = correcttreat & correctreg & correctID;
            plot(statssummary.Regions_expansion_radii_microns(condition),statssummary.median_Rings_int_per_pixel(condition));
            hold on;
        else
%             correcttreat = categorical(statssummary.gname(:,2)) == "iso";
%             condition = correcttreat & correctreg & correctID;
%             plot(statssummary.Regions_expansion_radii_microns(condition),statssummary.median_Rings_int_per_pixel(condition));
%             hold on;
        end
    end
end
ylim([0 inf])
xlabel("Distance from capillary, µm")
ylabel("Intensity per pixel")
title("Median values for each region",["Tissue ID: ", lookupID])
legend(char(unique(statssummary.gname(correctID,3))))

%%

% categorical(statssummary.gname(:,1))

% an_tab.bleh = repmat(working_tab.Blue_sat_pc,[2 1])
% 
% 
% for i = 1:length(working_tab.Properties.VariableNames)
%     working_tab.Properties.VariableNames{i}
% end
% 
% an_str = table2struct(an_tab,"ToScalar",true);



% 
% 
% figure("Name","Regions");
% subplot(3,3,1)  
% gscatter(working_tab.Regions_expansion_radii_microns(working_tab.Region_name == "Alveus",:),mean(working_tab.Rings_int_per_pixel(working_tab.Region_name == "Alveus",:),1),working_tab.Treatment(working_tab.Region_name == "Alveus",:))
% 
% [min,max,mean] = grpstats(working_tab.Rings_int_per_pixel,working_tab.Tissue_ID,{'min','max','mean',});
% xes = working_tab.Regions_expansion_radii_microns(1:22,1:18);
% figure;
% scatter(xes',mean')
% 
% rings_int_pp = [reg_int(1,:);reg_int(2:end,:)-reg_int(1:end-1,:)];
% rings_int_pp = rings_int_pp ./ rings_n_pix;
% rings_int_pp(rings_int_pp < 0 | isnan(rings_int_pp)) = nan;
% rings_exp_mic = reg_exp_rad_mic;
% 
% rings_int_pp_ab = rings_int_pp(:,treat == "ab");
%     rings_int_pp_ab = reshape(rings_int_pp_ab,[1,numel(rings_int_pp_ab)]);
% 
% rings_int_pp_iso = rings_int_pp(:,treat == "iso");
%     rings_int_pp_iso = reshape(rings_int_pp_iso,[1,numel(rings_int_pp_iso)]);
% 
% rings_exp_mic_ab = rings_exp_mic(1:numel(rings_int_pp_ab));
% rings_exp_mic_iso = rings_exp_mic(1:numel(rings_int_pp_iso));





% figure;
% gscatter(dat_tab.Regions_expansion_radii_microns,dat_tab.Regions_intensity,string(dat_tab.Region_name))
% 
% errorbar(mean(dat_tab.Regions_expansion_radii_microns,1),mean(dat_tab.Regions_intensity,1),stderror)
% 
% 
% stderror=std(dat_tab.Regions_intensity)/sqrt(length(dat_tab.Regions_intensity))
% 
% 

% an_tab = working_tab;
% for i=1:size(an_tab,2)
%     if size(an_tab(:,i),2) == 1
%         an_tab{:,i} = repmat(an_tab{:,i}(:,i),[1, size(an_tab.Regions_expansion_radii_microns,2)])
%     else
%         an_tab(:,i) = reshape(an_tab{:,i},[numel(an_tab{:,i}), 1]);
%     end
% end
% figure;
% scatter(working_tab.Regions_expansion_radii_microns,working_tab.Regions_intensity)


%     % Treatments
% tab = dat_tab(dat_tab.Treatment == "ab",:);
% tiso = dat_tab(dat_tab.Treatment == "iso",:);
%     % Region names
% tfim = 
% talv = 


% reg_exp_rad_mic = working_tab.Regions_expansion_radii_microns';
% reg_exp_rad_pix = working_tab.Regions_expansion_radii_pixels';
% reg_n_pix = working_tab.Regions_pixels';
% reg_int = working_tab.Regions_intensity';
% cap_area = working_tab.Capillary_area';
% cap_peri = working_tab.Capillary_perimeter';
% n_pix = working_tab.Total_pixels';
% blue_sat = working_tab.Blue_sat_pc';
% green_sat = working_tab.Green_sat_pc';
% red_sat = working_tab.Red_sat_pc';
% proj_name = working_tab.Project_name';
% %     proj_name = repmat(proj_name,
% tis_ID = working_tab.Tissue_ID';
%     
% red_int_tot = working_tab.Red_intensity_tot';
%     
% sec_no = working_tab.Section_no';
%     
% treat = working_tab.Treatment';
%     
% reg_name = working_tab.Region_name';
%     
% frame_no = working_tab.Frame_no';
% %     
% 
% %% Do some figures quick-like
% 
markersize = 1;
redcolor = "red";
greencolor = "green";
bluecolor = "blue";

figure;
subplot(1,3,1)
scatter(an_tab.datenum,an_tab.Red_sat_pc,2)
subplot(1,3,2)
scatter(an_tab.datenum,an_tab.Green_sat_pc,2)
subplot(1,3,3)
scatter(an_tab.datenum,an_tab.Blue_sat_pc,2)

% 
figure("Name","Saturated pixels");
subplot(3,3,1)
boxrst = boxchart(working_tab.Treatment,working_tab.Red_sat_pc,"JitterOutliers","on");
title('Red saturation by treatment')
xlabel("Treatment")
ylabel("% saturated pixels")
boxrst.BoxFaceColor = redcolor;
boxrst.MarkerColor = redcolor;
boxrst.MarkerSize = markersize;
subplot(3,3,2)
boxgst = boxchart(working_tab.Treatment,working_tab.Green_sat_pc,"JitterOutliers","on");
title('Green saturation by treatment')
xlabel("Treatment")
ylabel("% saturated pixels")
boxgst.BoxFaceColor = greencolor;
boxgst.MarkerColor = greencolor;
boxgst.MarkerSize = markersize;
subplot(3,3,3)
boxbst = boxchart(working_tab.Treatment,working_tab.Blue_sat_pc,"JitterOutliers","on");
title('Blue saturation by treatment')
xlabel("Treatment")
ylabel("% saturated pixels")
boxbst.BoxFaceColor = bluecolor;
boxbst.MarkerColor = bluecolor;
boxbst.MarkerSize = markersize;
subplot(3,3,4)
boxrss = boxchart(working_tab.Section_no,working_tab.Red_sat_pc,"JitterOutliers","on");
title('Red saturation by section number')
xlabel("Section number")
ylabel("% saturated pixels")
boxrss.BoxFaceColor = redcolor;
boxrss.MarkerColor = redcolor;
boxrss.MarkerSize = markersize;
subplot(3,3,5)
boxgss = boxchart(working_tab.Section_no,working_tab.Green_sat_pc,"JitterOutliers","on");
title('Green saturation by section number')
xlabel("Section number")
ylabel("% saturated pixels")
boxgss.BoxFaceColor = greencolor;
boxgss.MarkerColor = greencolor;
boxgss.MarkerSize = markersize;
subplot(3,3,6)
boxbss = boxchart(working_tab.Section_no,working_tab.Blue_sat_pc,"JitterOutliers","on");
title('Blue saturation by section number')
xlabel("Section number")
ylabel("% saturated pixels")
boxbss.BoxFaceColor = bluecolor;
boxbss.MarkerColor = bluecolor;
boxbss.MarkerSize = markersize;
subplot(3,3,7)
boxrsr = boxchart(working_tab.Region_name,working_tab.Red_sat_pc,"JitterOutliers","on");
title('Red saturation by region')
xlabel("Region")
ylabel("% saturated pixels")
boxrsr.BoxFaceColor = redcolor;
boxrsr.MarkerColor = redcolor;
boxrsr.MarkerSize = markersize;
subplot(3,3,8)
boxgsr = boxchart(working_tab.Region_name,working_tab.Green_sat_pc,"JitterOutliers","on");
title('Green saturation by region')
xlabel("Region")
ylabel("% saturated pixels")
boxgsr.BoxFaceColor = greencolor;
boxgsr.MarkerColor = greencolor;
boxgsr.MarkerSize = markersize;
subplot(3,3,9)
boxbsr = boxchart(working_tab.Region_name,working_tab.Blue_sat_pc,"JitterOutliers","on");
title('Blue saturation by region')
xlabel("Region")
ylabel("% saturated pixels")
boxbsr.BoxFaceColor = bluecolor;
boxbsr.MarkerColor = bluecolor;
boxbsr.MarkerSize = markersize;

end
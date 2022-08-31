function [image_table,early_breaks] = Pipeline_04(varargin)
% Pipeline_04


% Comments with 'revise' indicate things that should definitely be revisited

% Suggested dil_vec
% dil_vec = [linspace(0,5,11),...
%             linspace(6,10,5),...
%             linspace(12,20,5),...
%             linspace(23,35,5),...
%             linspace(40,70,4),...
%             linspace(80,140,4),...
%             linspace(170,200,2)]

%% Resolving input variables

for funinputs = 1:nargin
    if isstring(varargin{funinputs})
        if varargin{funinputs} == "mf"
            make_figures = varargin{funinputs+1};
        elseif varargin{funinputs} == "dil_vec"
            dil_vec = varargin{funinputs+1};
        elseif varargin{funinputs} == "lookat"
            lookat = varargin{funinputs+1};
        elseif varargin{funinputs} == "folder"
            data_tree = varargin{funinputs+1};
        end
    end
end

%% Hard variables

    % Resolution of images, used in mask calculations where we work with
    % distances in pixels
pixels_per_micron = 3.9;

simple_threshold_limit = 15;

max_time_per_dilation = 10;

%% Setup for dil

% Number of iterations of expansion, i.e. how many concentric disks do
% we want around the capillary
dil_n_radii = 10;
% Step size in microns, this is the requested increase in radius for
% each iteration
dil_step_size = 2;
% Step size in pixels
dil_expansion_step = dil_step_size * pixels_per_micron;

if ~exist("dil_vec","var")
    % Vector of radii of dilation in pixels
    dil_step_mat = round(linspace(0,dil_n_radii-1,dil_n_radii) * dil_expansion_step);
    % Actual radii of dilation in microns
    dil_actual_steps = dil_step_mat / pixels_per_micron;
else
    dil_n_radii = length(dil_vec);
    dil_step_mat = dil_vec * pixels_per_micron;
    dil_actual_steps = dil_step_mat / pixels_per_micron;
end

%% Setup folders

    % Enter path to the data tree
if ~exist("data_tree","var")
    data_tree = 'C:\Users\nikol\Temporary_data\IHC\HCAR1COL4hippoc';
end
cd(data_tree)

if ~exist("Pipeline_output","dir")
mkdir Pipeline_output % Creates a new folder named Pipeline_output in the current path, data_tree
end
output_path = [data_tree, 'Pipeline_output\'];
    % Revise: automatic saving of the output is not implemented. Output
    % must currently be saved manually.

image_structure = dir(strcat('*/*/*/Image.tif')); % Creates a structure with information on all images in the current path
image_table = struct2table(image_structure); % Converts the structure to a table (I am not sure if I need this)

n_rows = size(image_table,1);

early_breaks = zeros(n_rows,3);

%% Create relevant variables in the data table

image_table.Red_intensity_tot(:) = nan;
image_table.Green_intensity_tot(:) = nan;
image_table.Regions_expansion_radii_microns = double.empty(n_rows,0);
image_table.Regions_expansion_radii_pixels = double.empty(n_rows,0);
image_table.Regions_pixels      = double.empty(n_rows,0);
image_table.Regions_intensity   = double.empty(n_rows,0);
image_table.Capillary_area(:)   = nan;
image_table.Capillary_perimeter(:) = nan;
image_table.Nuclei_count(:)     = nan;
image_table.Total_pixels(:)     = nan;
image_table.Blue_sat_pc(:)      = nan;
image_table.Green_sat_pc(:)     = nan;
image_table.Red_sat_pc(:)       = nan;
image_table.Project_name(:)     = categorical("N/A");
image_table.Tissue_ID(:)        = categorical("N/A");
image_table.Section_no(:)       = categorical("N/A");
image_table.Treatment(:)        = categorical("N/A");
image_table.Region_name(:)      = categorical("N/A");
image_table.Frame_no(:)         = categorical("N/A");

%% Loop through all the images and add analysis results to the table

if exist("lookat","var")
    images_to_analyse = lookat; 
else
    images_to_analyse = 1:n_rows;
end

for image_no = images_to_analyse

    fprintf('Working on image %i out of %i. \n',image_no,n_rows)
%% Extracts blue, green, and red images from the Image in the specified folder

working_path_to_image = image_structure(image_no).folder;

raw_image = imread(strcat(working_path_to_image,'/Image.tif'));

blue_band = raw_image(:,:,1:3);
green_band = raw_image(:,:,4:6);
red_band = raw_image(:,:,7:9);
sum_composite = blue_band + green_band + red_band;

total_pixels = numel(blue_band(:,:,1));

% Loads images as I want it.

%% Creating mask from green image

% REVISE: Why not use the entire green band to create the intensity image?
green_gray = im2gray(green_band); % Converts green (capillary) band to grayscale
%green_gray = green_band(:,:,2); % Select only green band of green band
green_gray_foreground = green_gray;
%green_gray_foreground(green_gray_foreground<=simple_threshold_limit) = 0;
greenintensity_total = sum(green_gray_foreground,[1,2]);
simple_threshold_mask = green_gray > simple_threshold_limit;
greenmask = imbinarize(green_gray); % Creates a mask that covers capillaries using otsu's method (k-means clustering, k=2)
greenmask = greenmask & simple_threshold_mask;

% Revise the following. Weird masks may be resolved through regionprops.
% regionprops output EquivDiameter may resolve it. Only accept regions of
% size larger than a 1 micrometer diameter shape.
seclose     = strel('disk',8,0); 
seopen      = strel('disk',5,0); 

% At this point we dilate the mask, and then fill the holes within the mask
% so the interior is robustly filled out. We erode the mask back to its
% original size and the holes remain filled.
% As a last step we open the mask to remove small structures (smaller than
% approximately 1.25 µm).
greenmask_tempdil   = imdilate(greenmask,seclose);
greenmask_tempfill  = imfill(greenmask_tempdil,'holes');
greenmask_filled    = imerode(greenmask_tempfill,seclose);
greenmask_opened    = imopen(greenmask_filled,seopen);

greenmask_final = greenmask_opened;
% Revise dilated mask. The order of operations is very
% important. Lower bound could be set by isotype control
% histograms. Opening of low pixel areas should be checked
% again. Perhaps erythrocytes can be removed from the mask by
% comparing the size of the mask islands to the theoretical size
% of the cell.



% Seems to work fine. Needs testing on ambiguous images. NOTE: image 44 in 
% test SI 001 seems to be good for testing. 
% NOTE: maybe combining all the green images into one and finding the
% appropriate threshold in one go is better. (It was not)

%% Nuclei

% Revise: Do these calculations, please.
% Get number of nuclei, preferably those in focus.
% Do a neighborhood analysis of red intensity?
% Nucleus number can be compared against total intensity later.

% Try smoothing.
% Try multiple thresholds.
% Decolvolution?
% Spørg Camilla.
testfilter = fspecial('disk',4);

% blue_gray = im2gray(blue_band); % Include red and green band
blue_gray = blue_band(:,:,3); % Exclude red and green band

% blue_gray_filled = imfill(blue_gray);
% subplot(2,3,2)
% imshow(blue_gray_filled,[])
nucleus_mask_initial = imbinarize(blue_gray,"adaptive"); % Adaptive or not?
nucleus_mask_initial = nucleus_mask_initial & blue_gray > 10;
nucleus_mask_initial = imopen(nucleus_mask_initial,seopen);

nucleus_mask_dist = bwdist(~nucleus_mask_initial);
nucleus_mask_dist = imfill(nucleus_mask_dist); % Fill or not?
nucleus_mask_dist = imfilter(nucleus_mask_dist,testfilter); % Mean filter ok?

nucleus_water = watershed(-nucleus_mask_dist,8);
nucleus_water(~nucleus_mask_initial) = 0;
rgb = label2rgb(nucleus_water,'jet',[.5 .5 .5]);

nucleus_water_mask = nucleus_water > 0;
nucleus_perimeters = bwperim(nucleus_water_mask);

nucleus_struct = bwconncomp(nucleus_water_mask,4);
nucleus_count = nucleus_struct.NumObjects;


%% Iterative dilation, dil
% 
    % The mask used for analysis
dil_mask = greenmask_final;

dil_regions   =   false([size(dil_mask),dil_n_radii]);
dil_edges     =   false([size(dil_mask),dil_n_radii]);
dil_rings     =   false([size(dil_mask),dil_n_radii]);

dil_dist = bwdist(dil_mask);

    % Making the expanded regions / masks and their perimeters
for i = 1:dil_n_radii
    dil_img_dilation = dil_dist <= dil_step_mat(i); % ALERT! This line takes ages with large structuring elements
    dil_img_perim = bwperim(dil_img_dilation,4);
    dil_regions(:,:,i) = dil_img_dilation;
    dil_edges(:,:,i) = dil_img_perim;
end


dil_rings(:,:,1) = dil_regions(:,:,1);
for i = 1:size(dil_regions,3)-1
    dil_rings(:,:,i+1) = xor(dil_regions(:,:,i+1), dil_regions(:,:,i));
end 

edges_sum = sum(dil_edges,3)>0;



%% Masking red image

% REVISE: Why did I choose to use only the red band of the red band to
% create the intensity image?
red_gray = im2gray(red_band);
% red_gray = red_band(:,:,1); % Exclude blue and green bands

red_gray_foreground = red_gray;
%red_gray_foreground(red_gray_foreground <= simple_threshold_limit) = 0;
redintensity_total_notbg = sum(red_gray_foreground,[1,2]); % REVISE: sort out background pixels (<15)

red_gray_analysis = repmat(red_gray_foreground,1,1,dil_n_radii); % Revise: summing up the intensity inside regions should exclude background pixels

red_gray_regions = red_gray_analysis;
red_gray_regions(~dil_regions) = 0;
region_intensity  = squeeze(sum(red_gray_regions,[1,2]))';

region_areas      = squeeze(sum(dil_regions,[1,2]))';
region_areas_percentage = 100 * region_areas ./ numel(dil_regions(:,:,1));


    
%% Calculate capillary area and perimeter

    % Revise: should capillary area and perimeter be normalised to image
    % area, and should they be multiplied by section thickness to get
    % capillary wall area and capillary volume instead?

capillary_pixels = sum(greenmask_final,[1 2]);
capillary_area = 1/pixels_per_micron * capillary_pixels/pixels_per_micron;
total_area = total_pixels/(pixels_per_micron^2);
capillary_area_percentage = 100 * capillary_area / total_area;
capillary_pixels_percentage = 100 * capillary_pixels / total_pixels;

capillary_perimeter_struct = regionprops(greenmask_final,"Perimeter");
capillary_perimeter_sum = sum([capillary_perimeter_struct(:).Perimeter])/pixels_per_micron;

%% Calculate saturated pixels

    blue_sat = blue_band(:,:,3);
    n_sat_blue = length(find(blue_sat==255));
    green_sat = green_band(:,:,2);
    n_sat_green = length(find(green_sat==255));
    red_sat = red_band(:,:,1);
    n_sat_red = length(find(red_sat==255));
    n_pixels = numel(blue_sat);
    pc_sat_blue     = n_sat_blue*100/n_pixels;
    pc_sat_green    = n_sat_green*100/n_pixels;
    pc_sat_red      = n_sat_red*100/n_pixels;

%% Extract labels

    image_path = image_table.folder{image_no};
%     str_len = strlen(image_path);

    folder_separator_ind = find(image_path == '\');

    label_separators_ind = folder_separator_ind(end-3:end);

    project_name =  image_path(label_separators_ind(1)+1    :label_separators_ind(2)-1);
    section_ID   =  image_path(label_separators_ind(2)+1    :label_separators_ind(3)-1);
    region_name  =  image_path(label_separators_ind(3)+1    :label_separators_ind(4)-1);
    frame_no     =  image_path(label_separators_ind(4)+1    :end);

    section_ID_separator_ind = find(section_ID == ' ');

    tissue_ID   = section_ID(1:section_ID_separator_ind(1)-1);
    section_no  = section_ID(section_ID_separator_ind(1)+1:section_ID_separator_ind(2)-1);
    treatment   = section_ID(section_ID_separator_ind(2)+1:end);



%% Save results in table
% Remember to revise the output parameters

% if isempty(image_table.Regions_expansion_radii_microns)
    image_table.Regions_expansion_radii_microns(image_no,1:dil_n_radii) = dil_actual_steps;
    image_table.Regions_expansion_radii_pixels(image_no,1:dil_n_radii) = dil_step_mat;
    image_table.Regions_pixels(image_no,1:dil_n_radii) = region_areas;
    image_table.Regions_intensity(image_no,1:dil_n_radii) = region_intensity;
% else
%     n_current_values = length(image_table.Regions_expansion_radii_microns(image_no));
%     image_table.Regions_expansion_radii_microns(image_no,n_current_values+1:n_current_values+dil_n_radii+1) = dil_actual_steps;
%     image_table.Regions_expansion_radii_pixels(image_no,n_current_values+1:n_current_values+dil_n_radii+1) = dil_step_mat;
%     image_table.Regions_pixels(image_no,n_current_values+1:n_current_values+dil_n_radii+1) = region_areas;
%     image_table.Regions_intensity(image_no,n_current_values+1:n_current_values+dil_n_radii+1) = region_intensity;
% end

image_table.Red_intensity_tot(image_no) = redintensity_total_notbg;
    %Revise: intensity measures. Check the conversion to grayscale.
image_table.Green_intensity_tot(image_no) = greenintensity_total;
image_table.Capillary_area(image_no)    = capillary_area;
image_table.Capillary_perimeter(image_no) = capillary_perimeter_sum;
image_table.Total_pixels(image_no)      =   total_pixels;
image_table.Nuclei_count(image_no)      =   nucleus_count;
image_table.Blue_sat_pc(image_no)       =   pc_sat_blue;
image_table.Green_sat_pc(image_no)      =   pc_sat_green;
image_table.Red_sat_pc(image_no)        =   pc_sat_red;
image_table.Project_name(image_no,:)    =   categorical(cellstr(project_name));
image_table.Tissue_ID(image_no,:)       =   categorical(cellstr(tissue_ID));
image_table.Section_no(image_no,:)      =   categorical(cellstr(section_no));
image_table.Treatment(image_no,:)       =   categorical(cellstr(treatment));
image_table.Region_name(image_no,:)     =   categorical(cellstr(region_name));
image_table.Frame_no(image_no,:)        =   categorical(cellstr(frame_no));

%% Create figures for diagnostics and illustrations

if exist("make_figures","var")
    if make_figures == "all" || make_figures == "summary"
        figure('Name',['Bands of image ', num2str(image_no), ': ', tissue_ID, ' ', treatment, ' ', region_name]);
        subplot(2,2,1);
        blue_band_marked = blue_band;
        blue_band_marked(repmat(nucleus_perimeters,[1,1,3])) = 255;
        imshow(blue_band_marked);
        title('Blue band with nucleus mask outlined');
        subplot(2,2,2);
        green_band_marked = green_band;
        green_band_marked(cat(3,dil_edges(:,:,1),dil_edges(:,:,1),dil_edges(:,:,1))) = 255; % To display edges of mask
        imshow(green_band_marked);
        title('Green band with capillary mask outlined');
        subplot(2,2,3);
        red_band_marked = red_band;
        red_band_marked(repmat(max(dil_edges(:,:,:),[],3),1,1,3)) = 255;
        imshow(red_band_marked);
        title('Red band with rings outlined');
        subplot(2,2,4)
        imshow(sum_composite)
        title('Sum Composite')
        if make_figures == "all"
            figure;
            subplot(2,3,1)
            imshow(blue_gray,[])
            title("grayscale blue band")
            subplot(2,3,2)
            imshow(nucleus_mask_initial)
            title("initial mask")
            subplot(2,3,3)
            imshow(nucleus_mask_dist,[])
            title("distance transformed mask")
            subplot(2,3,4)
            imshow(rgb)
            title("Segmented nuclei")
            subplot(2,3,5)
            imshow(nucleus_perimeters)
            title("perimeters of segmented nuclei")
            subplot(2,3,6)
            imshow(blue_band_marked)
            title("Original DAPI image with segmentation perimeters indicated")

            figure('Name',['Creating a mask ',num2str(image_no)])
            subplot(2,3,1)
            imshow(green_gray);
            title('Green in grayscale');
            subplot(2,3,2)
            imshow(greenmask);
            title('Initial mask');
            subplot(2,3,3)
            imshow(greenmask_tempdil);
            title('Mask temporarily dilated');
            subplot(2,3,4)
            imshow(greenmask_tempfill);
            title('Mask filled');
            subplot(2,3,5)
            imshow(greenmask_filled);
            title('Mask eroded to original size, with fill');
            subplot(2,3,6)
            imshow(greenmask_opened);
            title('Final mask - opened to remove small structures (less than 5 pixels wide)')

            figure('Name','Iterative expansion of capillaries');
            subplot(2,2,1)
            montage(dil_regions)
            title('Montage of expanding regions')
            subplot(2,2,2)
            montage(dil_rings)
            title('Region rings')
            subplot(2,1,2)
            imshow(edges_sum)
            title('Perimeters of expanding regions')


            red_rings = table();
            red_rings.Intensity_tot = region_intensity;
            red_rings.Expansion_radius = dil_actual_steps;
            red_rings.Ring_pixels = region_areas;
            red_rings.Total_pixels = numel(dil_rings(:,:,1));

            figure('Name','Colocalisation analysis')
            subplot(1,2,1)
            red_marked = red_band;
            red_marked(repmat(edges_sum,[1,1,3])) = 255;
            imshow(red_marked);
            %montage(cat(3,red_gray_regions(:,:,1),red_gray_regions(:,:,2:10)-red_gray_regions(:,:,1:9)),[linspace(0,1,256)',zeros(256,1),zeros(256,1)])
            title("Montage of expanding rings")
            subplot(1,2,2)
            plot(red_rings.Expansion_radius,red_rings.Intensity_tot ./ red_rings.Ring_pixels)
            title('HCAR1 presence near capillaries')
            xlabel('Distance from capillary')
            ylabel('Intensity / pixel')
            ylim([0 inf])

            figure("Name","How much background can we ignore?");
            for i=1:16
                subplot(4,4,i)
                sample_reduction = red_band;
                sample_reduction(repmat(red_gray,[1,1,3]) < i*5-5) = 0;
                imshow(sample_reduction)
                title(i*5-5)
            end
        end
    else
        fprintf(['Oops!' ...
            '\nInput variable "mf" must be one of the following' ...
            '\n\t"all" \t\t- to display detailed figures illustrating many steps of the pipeline' ...
            '\n\t"summary" \t- to display only a figure showing the color bands and their associated mask boundaries' ...
            '\nI will halt the function early for your convenience.'])
        image_table = 0;
        early_breaks = 0;
        return
    end
    pause()
end


end
%% Output string
%     
%     fprintf(['\n%.2f %% of red intensity is contained within %.0f pixels of capillaries.\n' ...
%         '%.2f %% of red intensity is contained outside of that area.\n' ...
%         '%.2f %% of red intensity is contained within a randomly selected group of pixels \n' ...
%         '\t covering the same area as the capillary mask.\n' ...
%         ''], ...
%         percentintensity_near, cap_neighborhood_microns, percentintensity_away,percent_randmask);


% end


%% Revise figures for ppt
% ppt to illustrate pipeline
% cell count vs red int, log scale
% red band with background removed at different thresholds
% differentiate regions on figures
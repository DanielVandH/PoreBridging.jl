%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% DEFINE THE DIRECTORIES, EXTENSION, AND DAYS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all
addpath(genpath(cd));
base_directory = "data";
base_image_directory = "images";
base_boundary_directory = "boundaries";
base_mask_directory = "masks";
base_geometric_directory = base_directory;
save_table_directory = base_directory;
save_matlab_files_directory = base_directory;
geometric_image_directory = fullfile(base_directory, base_image_directory);
geometric_boundary_directory = fullfile(base_directory, base_boundary_directory);
geometric_mask_directory = fullfile(base_directory, base_mask_directory);
geometric_days = [5, 7, 14, 25, 28];
geometric_geo_types = ["Geo", "CRS500", "LWV500", "SWV500"];
ext = ".tif";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% GET ALL THE FILENAMES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
geometric_file_names = [];
geometric_unusable_path = fullfile(geometric_image_directory, "unusable" + ".txt");
geometric_unusable_files = readlines(geometric_unusable_path);
for d = geometric_days
    for geo_type = geometric_geo_types
        r = 1;
        file_path = fullfile(geometric_image_directory, "Day" + num2str(d), geo_type + "_" + num2str(r) + ext);
        file_attr = dir(file_path);
        while ~isempty(file_attr)
            if ~any(fullfile("Day" + num2str(d), geo_type + "_" + num2str(r)) == geometric_unusable_files)
                geometric_file_names = [geometric_file_names; file_path];
            end
            r = r + 1;
            file_path = fullfile(geometric_image_directory, "Day" + num2str(d), geo_type + "_" + num2str(r) + ext);
            file_attr = dir(file_path);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% GET THE SCALE BAR VALUES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_bar_val = [255,255,255];

microns_per_unit_length_geometric = [];
scalebar_lines_geometric = [];
scalebar_lengths_geometric = [];
figure(1)
for file_name = geometric_file_names'
    I = imread(file_name);
    imshow(I, []);
    len = inputdlg("Enter scale bar length:", "Input", [1 35]);
    len = str2double(len{1});
    scalebar_lengths_geometric = [scalebar_lengths_geometric; len];
    cla reset
end
k = 1;
for file_name = geometric_file_names'
    I = imread(file_name);
    len = scalebar_lengths_geometric(k);
    J = I(:, :, 1) == scale_bar_val(1) & ...
        I(:, :, 2) == scale_bar_val(2) & ...
        I(:, :, 3) == scale_bar_val(3);
    regions = regionprops(J);
    [~, max_area_idx] = max(cat(1, regions.Area));
    region_bboxs = cat(1, regions.BoundingBox);
    scalebar_bbox = region_bboxs(max_area_idx, :);
    scalebar_length = scalebar_bbox(3);
    microns_per_unit_length_geometric = [microns_per_unit_length_geometric; len / scalebar_length];
    cla reset
    k = k + 1;
end
close 1

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% GET THE BOUNDARIES OF EACH PORE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pore_boundaries_geometric = {};
pore_boundaries_x_geometric = {};
pore_boundaries_y_geometric = {};
for file_name = geometric_file_names'
    f = figure(1);
    I = imread(file_name);
    imshow(I, []);
    pb_s = {};
    px_s = {};
    py_s = {};
    while true %% Iterate over multiple masks
        [pb,px,py] = roipoly;
        if sum(pb, 'all') > 0 %% Once a polygon has been closed, click the same vertex three times to exit - this creates an empty mask that is detected by
            pb_s = [pb_s; pb];
            px_s = [px_s; px];
            py_s = [py_s; py];
        else
            break
        end
    end
    pore_boundaries_geometric = [pore_boundaries_geometric; {pb_s}];
    pore_boundaries_x_geometric = [pore_boundaries_x_geometric; {px_s}];
    pore_boundaries_y_geometric = [pore_boundaries_y_geometric; {py_s}];
    cla reset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% REMOVE THE BAD MASKS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bad_geometric = [];
k = 1;
figure(1)
for file_name = geometric_file_names'
    kth_pore_boundary = pore_boundaries_geometric{k};
    I = imread(file_name);
    for i = 1:length(kth_pore_boundary)
        imshowpair(I, kth_pore_boundary{i});
        res = inputdlg("Is this mask poor? Enter 0 for poor, 1 for good:" , "Input", [1 35]);
        res = str2double(res{1});
        if res == 0
            bad_geometric = [bad_geometric; k i];
        end
        cla reset
    end
    k = k + 1;
    cla reset
end
for i = 1:size(bad_geometric, 1) % Assumes that no masks in the same image were flagged for removal
    pore_boundaries_geometric{bad_geometric(i, 1)}(bad_geometric(i, 2)) = [];
    pore_boundaries_x_geometric{bad_geometric(i, 1)}(bad_geometric(i, 2)) = [];
    pore_boundaries_y_geometric{bad_geometric(i, 1)}(bad_geometric(i, 2)) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% COMPUTE ALL THE PORE STATISTICS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geometric_areas = {};
geometric_centroids = {};
geometric_perimeters = {};
geometric_smallest_rad = {};
geometric_bbox = {};
for k = 1:length(geometric_file_names)
    areas = [];
    centroids = [];
    perimeters = [];
    radii = [];
    bboxes = [];
    boundaries = pore_boundaries_geometric{k};
    scale = microns_per_unit_length_geometric(k);
    for r = 1:length(boundaries)
        pore = pore_boundaries_geometric{k}{r};
        pore_x = pore_boundaries_x_geometric{k}{r};
        pore_y = pore_boundaries_y_geometric{k}{r};
        props = regionprops(pore, 'Area', 'Perimeter', 'Centroid', 'BoundingBox');
        area = props.Area * scale^2;
        perimeter = props.Perimeter * scale; 
        centroid = props.Centroid;
        bbox = props.BoundingBox;
        [rad, ~] = get_smallest_distance(pore_x, pore_y, centroid);
        rad = rad * scale;
        areas = [areas; area];
        perimeters = [perimeters; perimeter];
        centroids = [centroids; centroid];
        radii = [radii; rad];
        bboxes = [bboxes; bbox];
    end
    geometric_areas = [geometric_areas; {areas}];
    geometric_perimeters = [geometric_perimeters; {perimeters}];
    geometric_centroids = [geometric_centroids; {centroids}];
    geometric_smallest_rad = [geometric_smallest_rad; {radii}];
    geometric_bbox = [geometric_bbox; {bboxes}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% COMPUTE BOUNDARIES OF ALL THE VOIDS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure(1);
set(f, 'KeyPressFcn', @reset_image)
pore_void_boundaries_geometric = {};
pore_void_boundaries_x_geometric = {};
pore_void_boundaries_y_geometric = {};
for k = 1:length(geometric_file_names)
    % Show the image
    file_name = geometric_file_names(k);
    I = imread(file_name);
    imshow(I, []);
    hold on
    boundaries = pore_boundaries_geometric{k};
    % Loop over the boundaries 
    pb_s = {};
    px_s = {};
    py_s = {};
    r = 1;
    while r <= length(boundaries)
        global skip_flag
        skip_flag = false;
        % Plot the boundary 
        pore_x = pore_boundaries_x_geometric{k}{r};
        pore_y = pore_boundaries_y_geometric{k}{r};
        h1 = plot(pore_x, pore_y, 'w-');
        [pb, px, py] = roipoly;
        if ~skip_flag 
            pb_s = [pb_s; pb];
            px_s = [px_s; px];
            py_s = [py_s; py];
            r = r + 1;
        end
        delete(h1)
    end
    pore_void_boundaries_geometric = [pore_void_boundaries_geometric; {pb_s}];
    pore_void_boundaries_x_geometric = [pore_void_boundaries_x_geometric; {px_s}];
    pore_void_boundaries_y_geometric = [pore_void_boundaries_y_geometric; {py_s}];
    cla reset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% GET THE RADIUS LENGTHS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
geometric_radius_line_lengths = {};
geometric_smallest_radius_positions = {};
for k = 1:length(geometric_file_names)
    % Show the image
    file_name = geometric_file_names(k);
    I = imread(file_name);
    imshow(I, []);
    hold on
    boundaries = pore_boundaries_geometric{k};
    % Get the scales 
    scale = microns_per_unit_length_geometric(k);
    % Loop over each mask 
    all_lengths = [];
    all_pos = [];
    for r = 1:length(boundaries)
        % Plot the boundary 
        pore_x = pore_boundaries_x_geometric{k}{r};
        pore_y = pore_boundaries_y_geometric{k}{r};
        pore_void_x = pore_void_boundaries_x_geometric{k}{r};
        pore_void_y = pore_void_boundaries_y_geometric{k}{r};
        cent = geometric_centroids{k}(r, :)';
        h1 = plot(pore_x, pore_y, 'w-');
        h2 = plot(cent(1), cent(2), 'w.', 'MarkerSize', 14);
        h3 = plot(pore_void_x, pore_void_y, 'm-');
        [len, p1p2] = get_smallest_distance(pore_void_x, pore_void_y, cent);
        len = len*scale;
        h4 = plot([cent(1), p1p2(1)], [cent(2), p1p2(2)]);
        pause(0.1);
        all_lengths = [all_lengths; len];
        all_pos = [all_pos; p1p2'];
    end
    geometric_radius_line_lengths = [geometric_radius_line_lengths; {all_lengths}];
    geometric_smallest_radius_positions = [geometric_smallest_radius_positions; {all_pos}];
    cla reset
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% GET THE SUMMARY STATISTICS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geometric_normalised_void_coverage = {};
geometric_normalised_void_perimeter = {};
geometric_normalised_smallest_radius = {};
geometric_void_area = {};
geometric_void_perimeter = {};
geometric_void_smallest_radius = {};
for k = 1:length(geometric_file_names)
    boundaries = pore_boundaries_geometric{k};
    coverages = [];
    perimeters = [];
    radii = [];
    unscaled_coverages = [];
    unscaled_perimeters = [];
    unscaled_radii = [];
    scale = microns_per_unit_length_geometric(k);
    for r = 1:length(boundaries)
        pore_void = pore_void_boundaries_geometric{k}{r};
        pore_void_x = pore_void_boundaries_x_geometric{k}{r};
        pore_void_y = pore_void_boundaries_y_geometric{k}{r};
        pore_void_props = regionprops(pore_void, 'Area', 'Perimeter');
        pore_area = geometric_areas{k}(r);
        pore_perimeter = geometric_perimeters{k}(r);
        pore_rad = geometric_smallest_rad{k}(r);
        pore_void_area = pore_void_props.Area * scale^2;
        pore_void_perimeter = pore_void_props.Perimeter * scale;
        pore_void_rad = geometric_radius_line_lengths{k}(r);
        coverages = [coverages; pore_void_area/pore_area];
        perimeters = [perimeters; pore_void_perimeter/pore_perimeter];
        radii = [radii; pore_void_rad/pore_rad];
        unscaled_coverages = [unscaled_coverages; pore_void_area];
        unscaled_perimeters = [unscaled_perimeters; pore_void_perimeter];
        unscaled_radii = [unscaled_radii; pore_void_rad];
    end
    geometric_normalised_void_coverage = [geometric_normalised_void_coverage; {coverages}];
    geometric_normalised_void_perimeter = [geometric_normalised_void_perimeter; {perimeters}];
    geometric_normalised_smallest_radius = [geometric_normalised_smallest_radius; {radii}];
    geometric_void_area = [geometric_void_area; {unscaled_coverages}];
    geometric_void_perimeter = [geometric_void_perimeter; {unscaled_perimeters}];
    geometric_void_smallest_radius = [geometric_void_smallest_radius; {unscaled_radii}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% PREPARE THE DATASET AND SAVE IMAGES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = [];
Day = [];
geoType = [];
imageIndex = [];
poreIndex = [];
poreArea = [];
porePerimeter = [];
poreRadius = [];
voidArea = [];
voidPerimeter = [];
voidRadius = [];
voidNormalisedCoverage = [];
voidNormalisedPerimeter = [];
voidNormalisedRadius = [];

for k = 1:length(geometric_file_names)
    file_name = geometric_file_names(k);
    % Find the day
    for m = 1:length(geometric_days)
        if contains(file_name, "Day" + geometric_days(m))
            day = geometric_days(m);
            break
        end
    end
    % Find the type
    for m = 1:length(geometric_geo_types)
        if contains(file_name, geometric_geo_types(m))
            geo_type = geometric_geo_types(m);
            break
        end
    end
    % Get the statistics
    all_pore_areas = geometric_areas{k};
    all_pore_perimeters = geometric_perimeters{k};
    all_pore_radii = geometric_smallest_rad{k};
    all_pore_smallest_radius_positions = geometric_smallest_radius_positions{k};
    all_pore_bboxes = geometric_bbox{k};
    all_pore_void_areas = geometric_void_area{k};
    all_pore_void_perimeters = geometric_void_perimeter{k};
    all_pore_void_radii = geometric_void_smallest_radius{k};
    all_normalised_pore_void_areas = geometric_normalised_void_coverage{k};
    all_normalised_pore_void_perimeters = geometric_normalised_void_perimeter{k};
    all_normalised_pore_void_smallest_radius = geometric_normalised_smallest_radius{k};
    all_pore_boundaries = pore_boundaries_geometric{k};
    all_pore_boundaries_x = pore_boundaries_x_geometric{k};
    all_pore_boundaries_y = pore_boundaries_y_geometric{k};
    all_pore_void_boundaries = pore_void_boundaries_geometric{k};
    all_pore_void_boundaries_x = pore_void_boundaries_x_geometric{k};
    all_pore_void_boundaries_y = pore_void_boundaries_y_geometric{k};
    all_pore_centroids = geometric_centroids{k};
    num_pores = length(all_pore_boundaries);
    for r = 1:num_pores
        % Update the table values
        fileName = [fileName; file_name];
        Day = [Day; day];
        geoType = [geoType; geo_type];
        imageIndex = [imageIndex; k];
        poreIndex = [poreIndex; r];
        poreArea = [poreArea; all_pore_areas(r)];
        porePerimeter = [porePerimeter; all_pore_perimeters(r)];
        poreRadius = [poreRadius; all_pore_radii(r)];
        voidArea = [voidArea; all_pore_void_areas(r)];
        voidPerimeter = [voidPerimeter; all_pore_void_perimeters(r)];
        voidRadius = [voidRadius; all_pore_void_radii(r)];
        voidNormalisedCoverage = [voidNormalisedCoverage; all_normalised_pore_void_areas(r)];
        voidNormalisedPerimeter = [voidNormalisedPerimeter; all_normalised_pore_void_perimeters(r)];
        voidNormalisedRadius = [voidNormalisedRadius; all_normalised_pore_void_smallest_radius(r)];
        % Save the mask
        fig = figure;
        I = imread(file_name);
        imshow(I, []);
        hold on
        bbox = all_pore_bboxes(r, :)';
        smallest_rad_pos = all_pore_smallest_radius_positions(r, :)';
        cent = all_pore_centroids(r, :)';
        pore_x = all_pore_boundaries_x{r};
        pore_y = all_pore_boundaries_y{r};
        pore_void_x = all_pore_void_boundaries_x{r};
        pore_void_y = all_pore_void_boundaries_y{r};
        plot(pore_x, pore_y, 'r-', 'LineWidth', 3);
        plot(pore_void_x, pore_void_y, 'm-', 'LineWidth', 3);
        plot([cent(1), smallest_rad_pos(1)], [cent(2), smallest_rad_pos(2)], 'b-', 'LineWidth', 3);
        xlim([bbox(1), bbox(1)+bbox(3)]);
        ylim([bbox(2), bbox(2)+bbox(4)]);
        save_path = fullfile(geometric_mask_directory, geo_type + "DAY" + num2str(day) + "IMAGE" + num2str(k) + "PORE" + num2str(r) + ".png");
        saveas(fig, save_path);
        % Save the pore boundaries
        pore_bnd = [pore_x pore_y];
        pore_void_bnd = [pore_void_x pore_void_y];
        pore_save_path = fullfile(geometric_boundary_directory, "PORE_" + geo_type + "DAY" + num2str(day) + "IMAGE" + num2str(k) + "PORE" + num2str(r) + ".dat");
        pore_void_save_path = fullfile(geometric_boundary_directory, "VOID_" + geo_type + "DAY" + num2str(day) + "IMAGE" + num2str(k) + "PORE" + num2str(r) + ".dat");
        writematrix(pore_bnd, pore_save_path, 'Delimiter', ' ');
        writematrix(pore_void_bnd, pore_void_save_path, 'Delimiter', ' ');
    end
end

close all
dataTable = table(fileName, Day, geoType, imageIndex, ...
    poreIndex, poreArea, porePerimeter, poreRadius, ...
    voidArea, voidPerimeter, voidRadius, ...
    voidNormalisedCoverage, voidNormalisedPerimeter, ...
    voidNormalisedRadius);
dataTable = sortrows(dataTable, 3);
save_path = fullfile(save_table_directory, "dataset.csv");
writetable(dataTable, save_path, 'QuoteStrings', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SAVE MATLAB FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_path = fullfile(save_matlab_files_directory, "files.mat");
save(save_path, ...
    'geometric_file_names', ...
    'geometric_days', ...
    'geometric_geo_types', ...
    'geometric_areas', ...
    'geometric_perimeters', ...
    'geometric_smallest_rad', ...
    'geometric_smallest_radius_positions', ...
    'geometric_bbox', ...
    'geometric_void_area', ...
    'geometric_void_perimeter', ...
    'geometric_void_smallest_radius', ...
    'geometric_normalised_void_coverage', ...
    'geometric_normalised_void_perimeter', ...
    'geometric_normalised_smallest_radius', ...
    'pore_boundaries_geometric', ...
    'pore_boundaries_x_geometric', ...
    'pore_boundaries_y_geometric', ...
    'pore_void_boundaries_geometric', ...
    'pore_void_boundaries_x_geometric', ...
    'pore_void_boundaries_y_geometric', ...
    'geometric_centroids')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% EXTRA FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reset_image(src, event)
    global skip_flag
    if event.Key == 'q'
        skip_flag = true;
    end
end

function [normqp1, p1p2] = smallest_distance(x1, x2, y1, y2, x, y)
    p1 = [x1; y1];
    p2 = [x2; y2];
    q = [x; y];
    qp1 = q - p1;
    p1p2 = p2 - p1;
    t = dot(qp1, p1p2) / dot(p1p2, p1p2);
    ts = max([t; 0]);
    ts = min([ts; 1]);
    p1p2 = p1 + ts * p1p2;
    qp1 = q - p1p2;
    normqp1 = norm(qp1);
end

function [d, p1p2] = get_smallest_distance(x, y, c)
    d = Inf;
    p1p2 = [0, 0];
    for i = 1:(length(x)-1)
        x1 = x(i);
        x2 = x(i+1);
        y1 = y(i);
        y2 = y(i+1);
        [newd, q] = smallest_distance(x1, x2, y1, y2, c(1), c(2));
        if newd < d 
            d = newd;
            p1p2 = q;
        end
    end
end
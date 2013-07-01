% Script for determining colony doubling time from plate reader growth data.
%
% Users of the script should update the following variables:
%     FILENAME: A tab-delimited text file with only the data and a single header
%         row with well names.
%     INTERVAL: Minutes between readings.
%     FLATTEN_FIRST_N_MINUTES: The number of readings to ignore at the
%         beginning of the script.
%
% The output is written to a new file in the same location as the input
% a text file, with extension '.analyzed_growth.csv'. This can be imported into
% Excel as a tab-delimited file.
%
% Written by Jaron Mercer based on original implementation by Harris Wang.
% Updated June 2013 by Gleb Kuznetsov.

clear;
close all;


%%% Set inputs.

% Windows style filename placeholder.
% FILENAME = 'C:\Users\Marc\Documents\0Harvard\Church Rotation\Church Laboratory Notebook\0rEcoli\Plate Reader Data\2013.06.14.C321.Fix.fitness.compare.top.strains_MOD.txt';

% Unix style filename placeholder.
FILENAME = '/home/glebk/Projects/churchlab/fix-recoli-1/experiment_data/2013.06.27.Fix.rE.coli.all.glycerol.replicates_mod.csv';

% Minutes separating each reading.
INTERVAL = 5;

% Sometimes the growth data shows irregular behavior at the beginning.
% We'll copy the value at the following value to all previous values
% to avoid getting an erroneous reading.
% Set this to 0 if you don't want any flattening.
FLATTEN_FIRST_N_MINUTES = 45;


%%% Begin processing

input = importdata(FILENAME, '\t', 1);

% Matrix where rows are consecutive time measurements and each column
% corresponds to a well.
data = input.data;

% Well names.
headers = input.colheaders;

% The number of wells.
num_wells = size(data, 2);

% Flatten the data. See comment for FLATTEN_FIRST_N_MINUTES above.
if FLATTEN_FIRST_N_MINUTES > 0
  flatten_first_n_observations = FLATTEN_FIRST_N_MINUTES / INTERVAL;
  for well = 1:num_wells
    % Copy the value at the nth position to all previous positions.
    data(1:flatten_first_n_observations, well) = ...
        data(flatten_first_n_observations, well);
  end
end


%% Plot data.
whos data;
plot(data);
xlabel('Time (min/5)');
ylabel('OD 600');
title(strcat(headers{1,1}, ' - ', headers{1,size(headers,2)}));


%% Plot log of data.
ln_data = log (data);
ln_data(isinf(ln_data)) = NaN;
ln_data = real(ln_data);
figure (2); plot(ln_data);
xlabel('Time (min/5)');
ylabel('ln(OD 600)');
title(strcat('ln(', headers{1,1}, ' - ', headers{1,size(headers,2)}, ')'));


%%% Find the linear portion of ln(OD 600). Then do linear regression.
%
% We dynamically search for the limits of the linear portion of the log-OD
% plot. We do this by performing a linear regression on a sliding interval
% window, and keeping track of the window with the greatest slope.
%
% This is the main "algorithmic" part of the script and the primary change
% Jaron Mercer made to Harris Wang's version of this script.

doubleTs = zeros(size(ln_data,2),1);
rSqrs = zeros(size(ln_data,2),1);
maxODs = zeros(size(ln_data,2),1);
deltas = zeros(size(ln_data,2),1);
starts = zeros(size(ln_data,2),1);
warnings = zeros(size(ln_data,2),1);

for well = 1:num_wells
    % Initialize state variables.
    maxSlope = 0;
    maxR = 0;
    maxDelta = 0;
    maxStart = 0;

    % Find the greatest slope.
    intervals_since_greatest = 0;
    for delta = 5:7
        for start = 1:(size(ln_data,1)-delta)
            % We expect the greatest interval early on so no need to go more
            % than 2 hours past max.
            if intervals_since_greatest > 24
                break
            end

            x = (linspace(start, start + delta - 1, delta))';
            y = ln_data(start: start + delta - 1, well);
            line = polyfit(x,y,1); % returns 1x2 matrix: [slope, y-intercept]

            if line(1,1) > maxSlope
                maxSlope = line(1,1);
                maxR = corrcoef(x,y);
                maxDelta = delta;
                maxStart = start;
                intervals_since_greatest = 0;
            else
                intervals_since_greatest = intervals_since_greatest + 1;
            end
        end
    end

    % Save output data.
    doubleTs(well, 1) = (log(2) / maxSlope) * INTERVAL;
    rSqrs(well, 1) = (maxR(1, 2)) ^ 2; %save r-squared
    maxODs(well, 1) = max(data(:, well));
    starts(well, 1) = maxStart * INTERVAL;
    deltas(well, 1) = maxDelta;

    % Output a warning for low r-squared values.
    if rSqrs(well, 1) < .99
        disp(strcat('Warning: low confidence on well --', headers(1, well)));
        figure(3); hold  on; plot(ln_data(:, well));
        warnings(well, 1) = 1;
    end
end


%%% Save data to a tab delimited text file.

output_filename = strcat(FILENAME(1:size(FILENAME, 2) - 4), '.analyzed_growth.csv');

output_data = [headers' num2cell(doubleTs) num2cell(rSqrs) num2cell(maxODs) num2cell(starts) num2cell(deltas) num2cell(warnings)];

fid = fopen(output_filename, 'w');

% Write the header row.
fprintf(fid, 'id\twell\tdoubling_time\tr_sqrd\tmax_OD\tstart_time_min\tdelta\twarnings\n');

% Iterate through the well data, writing one row at a time.
for well = 1:num_wells
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n', well, headers{well}, doubleTs(well), rSqrs(well), maxODs(well), starts(well), deltas(well), warnings(well));
end

fclose(fid);

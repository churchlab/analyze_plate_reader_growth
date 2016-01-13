function analyze_growth(filename, opt_interval)
% ANALYZE_GROWTH Compute doubling time given kinetic read time series.
%
%     Args:
%         filename: Full path to kinetic read data. Tab-delimited. First row is
%             well names. Each row is the the value of reads at each time point.
%         opt_interval: Optional. Interval between reads. Defaults to 5 min.
%         opt_flatten_first_n_minutes: Optional. Number of samples to flatten
%             at the beginning of the time series. Helps with issues due to
%             spurious fluctuations at beginning of read. Defaults to 45 min.
%
%     The output is written to a new file in the same location as the input
%     a text file, with extension '.analyzed_growth.csv'. This can be imported
%     into Excel as a tab-delimited file.
%
%     Example usage:
%         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt')
%
%
%     Written by Jaron Mercer based on original implementation by Harris Wang.
%
%     Updates by Gleb Kuznetsov:
%
%         07/2013: Optimization to ignore reads after max is found and other
%             cleanups.
%         06/2015: Convert into function so that users no longer need to edit
%             source code.
%
%     Updates by Tim Wannier:
%
%         01/2016: 
%             - Expanded search window (or delta) to 12 data points for
%               more reliable data.  
%             - Added functionality to blank data with wells that have media 
%               but no actively growing culture.  
%             - Removed arbitrary data flattening and replaced with an argument
%               to keep scanning the linear fit window until a minimum
%               R-squared is reached.

% Close any open windows.
close all;


%%% Parse args.

% Minutes separating each reading.
DEFAULT_INTERVAL = 5;

if exist('opt_interval')
    interval = opt_interval;
else
    interval = DEFAULT_INTERVAL;
end

% Get user input for blank well positions, i.e. those that are averaged and 
% asubtracted from the experimental data.
%
% This section was added by Tim W.

prompt = 'Pleae list the column positions of all sample blanks as a row vector.';
blank_columns = input(prompt);
blank_columns = sort(blank_columns,'descend');

%%% Begin processing

input_data = importdata(filename, '\t', 1);

% Matrix where rows are consecutive time measurements and each column
% corresponds to a well.
data = input_data.data;

% Well names.
headers = input_data.colheaders;

% Rows of data.
num_points = size(data,1);

% Separate blanks from data -- Aded by Tim W.
% The number of blanks.
num_blanks = size(blank_columns,2);

% Initiate blanks matrix.
blanks = zeros(num_points,num_blanks);

% Copy blanks to blanks matrix and remove from data and headers.
for well = 1:num_blanks
    blanks(:,well) = data(:,blank_columns(well));
    data(:,blank_columns(well)) = [];
    headers(:,blank_columns(well)) = [];
end

% Columns of data with no blanks.
num_wells = size(data, 2);

% Average blanks matrix and subtract from data.
blank_avg = mean(blanks.').';
data = data - blank_avg(:,ones(1,num_wells));


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

mus = zeros(size(ln_data,2),1);
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
    % Look at 12 data points (60 minutes) -- Changed by Tim W. from 5 to 7 
    % data points to reduce noise
    for delta = 12
        for start = 1:(size(ln_data,1)-delta)
            % We expect the greatest interval early on so no need to go more
            % than 2 hours past max. 
            % only break if the correlation is good, otherwise search exhausively
            %  -- Condition added by Tim W.
            if intervals_since_greatest > 24 & maxR > .995
                break
            end
            % Only bother fitting if the starting interval absorbance is between
            % 0.02 and 0.7, or within the plate-reader sweet-spot -- Condition added
            % by Tim W.
            if data(start,well) > 0.05 & data(start,well) < 0.7
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
    end

    % Save output data.
    mus(well, 1) = (maxSlope / interval) * 60; % Added µ to the list of outputs -- Tim W.
    doubleTs(well, 1) = (log(2) / maxSlope) * interval;
    rSqrs(well, 1) = (maxR(1, 2)) ^ 2; % save r-squared
    maxODs(well, 1) = max(data(:, well));
    starts(well, 1) = maxStart * interval;
    deltas(well, 1) = maxDelta;

    % Output a warning for low r-squared values.
    if rSqrs(well, 1) < .990
        disp(strcat('Warning: low confidence on well --', headers(1, well)));
        figure(3); hold  on; plot(ln_data(:, well));
        warnings(well, 1) = 1;
    end
end


%%% Save data to a tab delimited text file.

output_filename = strcat(filename(1:size(filename, 2) - 4), '.analyzed_growth.csv');

fid = fopen(output_filename, 'w');

% Write the header row.
fprintf(fid, 'id\twell\tµ_hourly\tdoubling_time_min\tr_sqrd\tmax_OD\tstart_time_min\tdelta\twarnings\n');

% Iterate through the well data, writing one row at a time.
for well = 1:num_wells
  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', well, headers{well}, mus(well), doubleTs(well), rSqrs(well), maxODs(well), starts(well), deltas(well), warnings(well));
end

fclose(fid);

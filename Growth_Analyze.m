%Written by Jaron Mercer.  Stolen in large part from Harris Wang.

%It needs a tab-delimited text file with nothing but the data and a single
%header row that has the well names.  The output file is a tab-delimited
%text file that copies right into Excel.

clear;
close all;

%% load file, set parameters

% Windows style filename placeholder.
% filename = 'C:\Users\Marc\Documents\0Harvard\Church Rotation\Church Laboratory Notebook\0rEcoli\Plate Reader Data\2013.06.14.C321.Fix.fitness.compare.top.strains_MOD.txt';

% Unix style filename placeholder.
filename = '/home/glebk/Projects/churchlab/fix-recoli-1/experiment_data/2013.06.27.Fix.rE.coli.all.glycerol.replicates_mod.csv'

% Minutes separating each reading.
interval = 5;

% Sometimes the growth data shows irregular behavior at the beginning.
% We'll copy the value at the following value to all previous values
% to avoid getting an erroneous reading.
% Set this to 0 if you don't want any flattening.
FLATTEN_FIRST_N_MINUTES = 60;

input = importdata(filename, '\t', 1);
%access by:  input.colheaders{1,k} and input.data{:,k}
data = input.data;
headers = input.colheaders;

% Flatten the data. See comment for FLATTEN_FIRST_N_MINUTES above.
if FLATTEN_FIRST_N_MINUTES > 0
  flatten_first_n_observations = FLATTEN_FIRST_N_MINUTES / interval;
  for well_column = 1:size(data, 2)
    % Copy the value at the nth position to all previous positions.
    data(1:flatten_first_n_observations, well_column) = ...
        data(flatten_first_n_observations, well_column);
  end
end

%% plot data
whos data;
plot (data);

xlabel('Time (min/5)');
ylabel('OD 600');
title(strcat(headers{1,1}, ' - ', headers{1,size(headers,2)}));
%axis([xmin xmax ymin ymax]);


%% plot log of data%
ln_data = log (data);
ln_data(isinf(ln_data)) = NaN;
ln_data = real(ln_data);
figure (2); plot (ln_data);

xlabel('Time (min/5)');
ylabel('ln(OD 600)');
title(strcat('ln(', headers{1,1}, ' - ', headers{1,size(headers,2)}, ')'));
%axis([xmin xmax ymin ymax]);


%% find linear portion of ln(OD 600). then do linear regression.
%This is the tricky part, and the only real change I made to Harris'
%program.  I'm trying to dynaimcally find the limits of the linear portion 
%of each ln(OD 600) plot.
%This does linear regression on various time intervals in the data and
%then takes the highest slope to use in calculating the doubling time.

doubleTs = zeros(size(ln_data,2),1);
rSqrs = zeros(size(ln_data,2),1);
maxODs = zeros(size(ln_data,2),1);
deltas = zeros(size(ln_data,2),1);
starts = zeros(size(ln_data,2),1);
warnings = zeros(size(ln_data,2),1);


for i = 1:size(ln_data,2)
    maxSlope = 0;
    for delta = 5:9
        for j = 1:(size(ln_data,1)-delta)
            x = (linspace(j, j+delta-1, delta))';
            y = ln_data(j:j+delta-1,i);
            line = polyfit(x,y,1); %returns 1x2 matrix: [slope, y-intercept]

            if line(1,1) > maxSlope
                maxSlope = line(1,1);
                maxR = corrcoef(x,y);
                maxDelta = delta;
                maxStart = j;
            end
        end
    end

    doubleTs(i,1) = log(2)/maxSlope*interval; %save doubling time
    rSqrs(i,1) = (maxR(1,2))^2; %save r-squared
    maxODs(i,1) = max(data(:,i));
    starts(i,1) = maxStart;
    deltas(i,1) = maxDelta;

    if rSqrs(i,1) < .99 %output a warning for low r-squared values
        disp(strcat('Warning: low confidence on well --', headers(1,i)));
        figure(3); hold  on; plot(ln_data(:,i));
        %plots ln(OD600) on figure(3) for low r-squared wells
        warnings(i,1) = 1;
    end

end

%% save data to a tab delimited text file: originalfilename_DoubleT.txt
dlmwrite([filename(1:size(filename,2)-4) '_Analyzed_MJL.txt'],[doubleTs rSqrs maxODs starts deltas warnings], '\t')

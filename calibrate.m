% Script to aid in empirical calibration.

% Generate outputs for different mid-log window sizes.
window_range = 15 + (0:13) * 5;
for window = window_range
    output = ['starting_strains_' num2str(window) '.csv'];
    analyze_growth(GROWTH_CURVES, [96], 0, 5, window, true, output)
end

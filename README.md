ANALYZE_GROWTH Compute doubling time given kinetic read time series.

    Args:
        filename: Full path to kinetic read data. Tab-delimited. First row is
            well names. Each row is the the value of reads at each time point.
        opt_blank_wells: Optional array of wells left blank for calibration.
            Provide as integers. E.g. [48, 96] for wells D12 and H12.
        opt_blank_value: Optional value of blank read. If provided,
            opt_blank_wells is ignored.
        opt_interval: Optional. Interval, in minutes, between reads.
            Default 5 min.
        opt_mid_log_interval: Optional. Size of window, in minutes, for which
            we measure linear growth. Default 60 minutes.
        opt_hide_plots: Optional. Boolean. If true, hide plots.
        opt_output_override: Optional. Override the name of the output file.

    The output is written to a new file in the same location as the input
    a text file, with extension '.analyzed_growth.csv'. This can be imported
    into Excel as a tab-delimited file.

    Example usage:

        analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt')

    Blank wells may be provided to be sourced as the average blank read.
    For example, if well H12 (96th well) is blank:

        analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [96])

    Alternatively, one can provide a global blank value to adjust all reads by:

        analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [], 0.09)

    Note for using optional arguments: User must provide values for all
    arguments up to the optional argument you would like to use.
    For example, to set opt_mid_log_interval to 30 min, the command is:

        analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt', [], 0, 5, 30)

    Authors:
        Jaron Mercer - ported to Matlab
        Harris Wang - original draft
        Gleb Kuznetsov (kuznetsov@g.harvard.edu) - primary maintainer
        Tim Wannier

ANALYZE_GROWTH Compute doubling time given kinetic read time series.

     Args:
         filename: Full path to kinetic read data. Tab-delimited. First row is
             well names. Each row is the the value of reads at each time point.
         opt_interval: Optional. Interval between reads. Defaults to 5 min.
         opt_flatten_first_n_minutes: Optional. Number of samples to flatten
             at the beginning of the time series. Helps with issues due to
             spurious fluctuations at beginning of read. Defaults to 45 min.

     The output is written to a new file in the same location as the input
     a text file, with extension '.analyzed_growth.csv'. This can be imported
     into Excel as a tab-delimited file.

     Example usage:
         analyze_growth('/home/glebk/Data/2015_06_10_growth_test.txt')


     Written by Jaron Mercer based on original implementation by Harris Wang.

     Updates by Gleb Kuznetsov:

         07/2013: Optimization to ignore reads after max is found and other
             cleanups.
         06/2015: Convert into function so that users no longer need to edit
             source code.

     Updates by Tim Wannier:

         01/2016: 
             - Expanded search window (or delta) to 12 data points for
               more reliable data.  
             - Added functionality to blank data with wells that have media 
               but no actively growing culture.  Data must be input as a row
               vector of indeces where A1 = 1, A12 = 12, B1 = 13, etc.
             - Removed arbitrary data flattening and replaced with an argument
               to keep scanning the linear fit window until a minimum
               R-squared is reached.

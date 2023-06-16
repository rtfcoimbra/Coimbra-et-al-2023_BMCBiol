# Usage: python3 site_depth_stats_v2.py <input> <output>

import numpy as np
from scipy import stats  # scipy >= v1.5.0
from sys import argv

# open input file (read mode) and output file (write mode)
with open(argv[1], 'r') as fin, open(argv[2], 'w') as fout:

    # create numpy array from input file lines
    depth = np.fromiter((line for line in fin), int)

    # calculate 5th and 95th percentiles
    perc_05 = np.percentile(depth, 5)
    perc_95 = np.percentile(depth, 95)

    # calculate median and median absolute deviation (MAD)
    median = np.median(depth)
    mad = stats.median_abs_deviation(depth)

    # calculate min. and max. depth cutoffs based on the median and MAD
    min_dp = median - (5 * mad)
    max_dp = median + (5 * mad)

    # write output file
    fout.write(f'5th percentile: {perc_05:.0f}\n')
    fout.write(f'95th percentile: {perc_95:.0f}\n')
    fout.write(f'Median: {median:.0f}\n')
    fout.write(f'MAD: {mad:.0f}\n')
    fout.write(f'Minimum depth cutoff (MEDIAN - 5 * MAD): {min_dp:.0f}\n')
    fout.write(f'Maximum depth cutoff (MEDIAN + 5 * MAD): {max_dp:.0f}')
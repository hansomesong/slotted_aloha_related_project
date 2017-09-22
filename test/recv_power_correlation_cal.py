__author__ = 'qsong'


import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr

if __name__ == "__main__":
    csv_file = os.path.join("..", "xx")
    csv_df = pd.read_csv(csv_file, sep=',', header=None, dtype=np.float64)
    print pearsonr(csv_df.values[1:, 0], csv_df.values[1:, 1])
    print pearsonr([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])


__author__ = 'qsong'

import pandas as pd
import numpy as np

csv_file = "/Users/qsong/Downloads/5_probes_to_alexa_top510_all_completed.csv"

df = pd.read_csv(csv_file, sep=";").dropna(axis=1, how='all')

df1 = df.iloc[:, 0:3]

df_min_rtt = df.iloc[:, 3:].apply(lambda x: x.str.split("/", expand=True)[0]).astype(float)

df_min = pd.concat([df1, df_min_rtt], axis=1)


drop_measure_id = set(df_min[df_min_rtt.apply(pd.Series.nunique, axis=1) == 1].iloc[:, 0].values)
print drop_measure_id
print df_min['measurement_id'].isin(drop_measure_id)
df_new = df_min[~df_min['measurement_id'].isin(drop_measure_id)]







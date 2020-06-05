import pandas as pd
import numpy as np
import utility_functions as uf
import utility_functions2 as uf2
'''
# df = pd.DataFrame({'gene': ['A1BG', np.nan, 'P1', 'YYYY'], 'id': [1, 2, 3, 5],
#                    'col3': [11, 12, 13, 14], 'col4': [15, 19, 21, 23]})

df = pd.DataFrame({'gene': ['A1BG', 'PRB1', 'FOXA1', 'TULP2'], 'id': [1, 2, 3, 5],
                   'col3': [11, 12, 13, 14], 'col4': [15, 19, 21, 23]})

df.set_index('gene', inplace=True)
df2 = df.copy()

print(gl)
print(gl2)
# print(type(gl.iloc[1, 1]), type(gl2.iloc[1, 1]))
print(gl.eq(gl2))
'''
'''
df = pd.DataFrame({'gene': ['A1BG', 'PRB1', 'FOXA1', 'TULP2'] * int(10000/4), 'id': np.random.rand(10000),
                   'col3': np.random.rand(10000), 'col4': np.random.rand(10000)})
df.set_index('gene', inplace=True)
df2 = df.copy()

lg = uf.log2(df)
lg2 = uf2.log2(df2)

print(lg)
print('----------')
print(lg2)

assert(df.equals(df2))
'''
'''
df = pd.DataFrame({'gene': range(4), 'id': np.random.randint(100, size=4),
                   'col3': np.random.randint(100, size=4), 'col4': np.random.randint(100, size=4)})
df.set_index('gene', inplace=True)
uf2.zscore(df, axis='row', epsilon=0)
df2 = df.copy()

sm = uf.createStandardizedMatrix(df)
sm2 = uf2.createStandardizedMatrix(df)

# print(sm)
# print(sm2)

assert(sm.equals(sm2))
'''

df = pd.DataFrame({'gene': range(4), 'id': np.random.randint(100, size=4),
                   'col3': np.random.randint(100, size=4), 'col4': np.random.randint(100, size=4)})
uf2.createUpAttributeSetLib(df, '', '')

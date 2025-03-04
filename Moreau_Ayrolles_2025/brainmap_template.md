```python
import pandas as pd
import nibabel as nib
import pathlib as pal
import numpy as np
np.bool = np.bool_
from matplotlib import pyplot as plt
import sys
import os
import seaborn as sns
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
```


```python
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from enigmatoolbox.plotting import plot_subcortical
```


```python
DATA = pd.read_csv("your_file", na_values=["NA"], delimiter=";")

```

## Plot your effect sizes in a brain map


```python
#cortical
coef = DATA['d_val']
coef_HighP_fsa5 = parcel_to_surface(coef, 'aparc_fsa5')
plot_cortical(array_name=coef_HighP_fsa5, 
              filename="submission/Figures/cc_AN_TD_CT.png" , screenshot = "True", 
              surface_name="fsa5", size=(800, 200),interactive=False,
              background=(1, 1, 1), color_range=(-1, 1),cmap='seismic', #cmap='RdYlBu',
               color_bar=False)
```


```python
#subcortical
coef = DATA['Dval_adj']
plot_subcortical(array_name=coef, size=(800, 400), 
                 filename="submission/Figures/cc_AN_TD_ASEG.png" ,  screenshot = "True", 
                 cmap='seismic', color_bar=False, interactive=False, color_range=(-1, 1))
```

## Spin permutations


```python
from enigmatoolbox.permutation_testing import spin_test, shuf_test
```


```python
map1=DATA.asd_meta_d_icv
map2=DATA.Dval
perm_p, perm_null = spin_test(map1, map2, 
                               surface_name='fsa5', 
                               parcellation_name='aparc',
                               type='pearson', 
                               n_rot=5000, 
                               null_dist=True)
```

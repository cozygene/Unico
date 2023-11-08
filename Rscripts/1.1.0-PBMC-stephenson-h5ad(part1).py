#!/usr/bin/env python
# coding: utf-8
#https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10026
# In[1]:


import os
import scanpy as sc
import numpy as np
import pandas as pd

np.random.seed(2023)


# In[2]:


source_file = "../Data/RNA/Simulation-PBMC/covid_portal_210320_with_raw.h5ad"
res_dir = "../Data/RNA/Simulation-PBMC"


# In[3]:


adata = sc.read_h5ad(source_file)


# In[4]:


adata


# In[5]:


# Remove LPS_10hours and LPS_90mins samples
cells = (adata.obs["Status_on_day_collection_summary"] != "LPS_10hours") & (adata.obs["Status_on_day_collection_summary"] != "LPS_90mins")
adata = adata[cells,]


# In[6]:


adata


# In[7]:


# Exclude samples that are resampling of the same patient; arbitrarily take the first sample.
keep_samples = []
patients = np.unique(adata.obs["patient_id"])
for patient in patients:
    keep_samples.append(np.unique(adata.obs["sample_id"][adata.obs["patient_id"] == patient])[0])


# In[8]:


print(len(keep_samples))
keep_samples[0:10]


# In[9]:


adata = adata[np.isin(adata.obs["sample_id"], keep_samples),:]
assert(len(np.unique(adata.obs["patient_id"])) == len(np.unique(adata.obs["sample_id"])))


# In[10]:


adata


# In[11]:


adata.write_h5ad(f"{res_dir}/stephenson.h5ad")



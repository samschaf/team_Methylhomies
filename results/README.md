
Results
========
in this folder we will discuss the imput and outputs of the code as well as interpretations. the breakdown of our codes can be found in the [src folder](https://github.com/STAT540-UBC/team_Methylhomies/tree/master/src/final_codes)

Processing
-------------

Overview of processing steps:

![](/Images/Pipeline_of_Methods.png)

**Aim**

**Imput: **

Dasen data = file: 

QN data = file:

meta data = file

**Output:**
  - files
  
--Dasen filtered data = file: GSE59685_batch_cor.RData

--Dasen filtered + cell corrected data = file:GSE59685_cell_cor.RData

--meta data = file: Meta_brain_cell_batch.RData

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

  - images
  
**conclusion**

Analysis
--------

Sample to sample correlation:

**Aim**
to determine whether the samples correlate well with eachother. If processing normalized the data we should see that samples will correlate with samples within the same cell type and outside of the same cell type.

**inputs**

--Dasen filtered data = file: GSE59685_batch_cor.RData

--Dasen filtered + cell corrected data = file:GSE59685_cell_cor.RData

--meta data = file: Meta_brain_cell_batch.RData

--QN + BMIQ filtered data = file: GSE43414_batch_cor.RData

--QN + BMIQ filtered + cell type corrected data = file: GSE43414_cell_cor.RData

--meta data = file:Meta_batch_cor.RData

**outputs**
- images

![](/Images/Pipeline_of_Methods.png)

**conclusions**




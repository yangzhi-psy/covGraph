
# Code for structural-functional co-variance project
**Author: Zhi Yang Ph.D., Institute of Psychology, Chinese Academy of Sciences**

This repository contains code for implementing a discovery-validation scheme used to learn brain structural-functional co-variance associations from large-samples of neuroimaging data

___

**00.scpt_prepareSingleMetricICA.m**

- script to combine metric maps from all subjects for ICA.
- dependencies
  - subjects_harvard.list
  - subjects_swu.list
  - metric surface maps: 'alff', 'area', 'dc', 'ec', 'falff', 'lgi', 'meancurv', 'reho', 'reh o2', 'sulc', 'thickness', 'volume'
  - fnc_collectVertices.m
  - nii_template.mat
- output
  - combined single-metric surface maps, like: forICA_alff.nii.gz


**01.scpt_performSingleMetricICA.sh**

- script to call melodic to perform ICA.
- dependencies
  - combined single-metric surface maps, like: forICA_alff.nii.gz
  - mask: allOne.nii.gz
- output
  - ICA components stored in alff.ica

**02.scpt_discoverMMCU.m**

- script to discover MMCU from dataset 1.
- dependencies
  - ICA mixing time series in alff.ica
  - gRAICAR (download: https://github.com/yangzhi-psy/gRAICAR)
  - examConnections.m
- output
  - matched components
  - similarity matrices of MMCUs from discovery dataset

**03.scpt_validate.m**

- script to validate MMCUs in dataset 2.
  - the MMCUs found in dataset1 are used as spatial templates
  - spatial regression is used to compute mixing course in dataset2
  - similarity matrices between mixing courses obtained in dataset2 are then computed
  - ICC are computed between the similarity matrics obtained in dataset1 and dataset2.
- dependencies
  - matched components
  - similarity matrices of MMCUs from discovery dataset
  - ICA components from dataset2, stored in alff.ica
  - ICC.m
  - IPN_kendallW.m
- output
  - similarity matrics of MMCUs obtained from dataset2
  - crsST_cc_HVU_SWU.mat

**04.scpt_permuteMMCU.m**
- script for permuting SCMs to generate fake MMCUs and compute their ICCs across  sites
- dependencies
  - ICA components from dataset1 and dataset2, stored in alff.ica
  - ICC.m
- output
  - null distribution of ICCs
  - null_crsST_ICC_permuteMMCU.mat

**05.scpt_overallCov.m**
- script to compute overall co-variance between 12 metrics
  - combined validated MMCUs (ICC>0.6)
- dependencies
  - similarity matrics of MMCUs obtained from dataset1 and dataset2
- output
  - a figure of similarity matrix and a dendrograph.

# mng-rangeland

Is observed change in grassland vegetation index due to rising herding sizes, changing climate, both or neither?

Summary: 

This study explores the relationship between livestock populations and rangeland vegetation conditions using remote sensing at the soum level (the second smallest administrative subdivision of Mongolia), controlling for extreme weather events. We have assembled 35 years of spatially disaggregated data, allowing for a more precise estimation of the correlates of rangeland degradation than has been possible previously. 

In this repo, I uploaded Python script used to download two main datasets: vegetitation indices and weather indicators. The `01-script/gee-python` folder contains jupyter notebook, where the data is accessed and downlaoded, and four python helper modules used in downloading stages. The `environment.yml` file contains important packages and dependencies used in this analysis. 

Simple workflow to reproduce the files downloaded in this study: 

1. Clone this repository
2. Install all the packages and dependencies stored in `environment.yml`. [For more information on creating an environment from an environment.yml file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)
3. Sign-up for Earth-Engine access. [For more information](https://developers.google.com/earth-engine/guides/access)
4. Run jupyter notebook `data_acquisition.ipynb`

The output for vegetation indices will be .csv files for each year by soum regions. The csv files are downloaded and stored in the Google Drive folder named after '--SOUM_LANDSAT_VEG' or '--SOUM_MODIS_VEG' associated with the Google Earth Engine. Similarly, weather outputs will be in the format of csv files for each month by soum regions stored in Google Drive folder named '--SOUM_ERA_WEATHER'. 

If you have any questions, please don't hesitate to contact me at ta346@cornell.edu



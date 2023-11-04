# Disentangling Drivers of Rangeland Degradation: Herd Size versus Climate in Mongolia, 1985-2023

## Summary

This study explores the relationship between livestock herd size and rangeland productivity using remote sensing at the soum level (the second smallest administrative subdivision of Mongolia), controlling for extreme weather events. We have assembled 35 years of spatially disaggregated data, allowing for a more precise estimation of rangeland degradation than has been possible previously.

The current version of draft can be accessed at [here](https://drive.google.com/file/d/1ueFvNf86GHdMPXUL3Mb0OvVT3nVbw2Ji/view)

## Data Acquisition

In this repository, you will find Python scripts used to download main raw datasets used in this study: vegetation indices, weather indicators, plot level vegetation indices. The `01-scripts` directory contains a Jupyter notebook (`data_acquisition.ipynb`) where the data is accessed and downloaded, as well as three Python helper modules used in the data acquisition process. You can set up the required environment using the `environment.yml` file, which includes the necessary packages and dependencies used in this analysis.

## Reproducing the Data

To reproduce the data downloaded in this study, follow the following steps:

1. Clone this repository to your local machine.
2. Install all the packages and dependencies stored in `environment.yml`. You can create an environment from the YAML file using the conda package management system.
3. Sign up for Earth Engine access by visiting the Google Earth Engine's guide [here](https://earthengine.google.com/signup/).
4. Run the Jupyter notebook `data_acquisition.ipynb` to download the data.
5. The output for vegetation indices will be CSV files for each month by soum regions. The CSV files will be downloaded and stored in a Google Drive folder named `--SOUM_LANDSAT_VEG` or `--SOUM_MODIS_VEG` associated with your Google Earth Engine account.
6. Similarly, weather data outputs will be in the format of CSV files for each month by soum regions, and they will be stored in a Google Drive folder named `--SOUM_ERA_WEATHER`.

If you have any questions or need further assistance, please don't hesitate to contact me at ta346@cornell.edu.

# STRAT2-MERGIR-IR
Python script to download and analyze NASA 4-km MERGIR IR from GES DISC OpenDAP Server

**Created by Doug Goetz May 2021**

**Local Directory:**
/User Directory/Mergirdownload.py

**Overview:**
This script performs partial downloads to memory of the MERGIR IR dataset based on an input .csv file that contains time,latitude, and longitude. The input file is intended to contain 30min averaged Stratéole2 EUROS data, but any spatial dataset can be used. The script downloads subsections of the MERGIR data based on the time and a user defined bounding box centered around the gondola position. Each downloaded brightness temperature dataset is saved as .png and overlaid on a map. The IR dataset is also processed to find distance from the gondola position to the nearest brightness temperature pixel below a user defined threshold and the average IR brightness temperature below the gondola within a user defined domain.

The MERGIR


**User defined variables:**

gond_pos_csv: string with full path of input csv file containing gondola time and position

LatBoundarySize: +/-  latitude degrees from rounded gondola postion that defines download bounding box (10 degrees by default)

LonBoundarySize: +/- longitude degrees from rounded gondola postion that defines download bounding box (10 degrees by default)

BTthreshold: the brightness temperature value used to define a convective cell. Defined as 235K in Corcos et al, 2021, JGR Atmospheres

BTdomain: the size domain in km from gondola used to find the average brightness temperature below the gondola (20km has been used with OK results)

**Input file structure:**

the input file should be a .csv file without headers that contains a single column for time since 01/01/2019, gondola latitude in dd.ddddd, and gondola longitude in dd.ddddd

**the Mergir dataset contains imagery at 30 minute and 4km resolution. It is recommended to average the input data to 30 minutes at hh:00 and hh:30. Otherwise the script could take a verrry long time to run.

**Outfile file structure**

An .csv file containing the calculated distance from convection and average brightness temperature below the gondola is created within the local directory. Addtionally, the URL to find/download the imagery for that data point from the OpenDAP server is included

.png files for each unique gondola position and time is created within the directory: /User Directory/MergirImages/



 
 



  



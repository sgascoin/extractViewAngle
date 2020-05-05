# extractViewAngle
[`extractViewAngle.py`](https://github.com/sgascoin/extractViewAngle/blob/master/extractViewAngle.py)  is a script to export points or raster of viewing incidences angles from a Theia L2A product. In raster mode it will write a multiband geotiff file for each band and angle (zenith and azimuth) scaled by 100 and stored as UInt16. In point mode it will write the value of the incidence angles at every point provided in geographic coordinates (no data is -10000).

usage in raster mode: 

`extractViewAngle.py productFolder outputFolder`

usage in point mode (space delimited values): 

`extractViewAngle.py productFolder outputFolder point_table_as_lon_lat.txt`

examples:

`python extractViewAngle.py SENTINEL2A_20180224-103018-463_L2A_T31TGK_C_V2-2 angles`

`python extractViewAngle.py SENTINEL2A_20180224-103018-463_L2A_T31TGK_C_V2-2 angles points.csv`

where `points.csv` contains a list of longitude and latitude separated by a comma:

> 6.16800,44.88358  
> 6.40165,45.03698  
etc.

In raster mode it's really just the reformatting of the metadata file (MDT xml file). In point mode the code additionally checks the detector footprints to solve for ambiguities in areas where the MTD angle grids overlap. No resampling is performed.

It should work for any version of Theia L2A (tested for v1-1 and v2-2). If the detector footprints are not available, in case of ambiguity the code takes the value of the first detector by default. In this configuration, the code should also work with L1C and ESA L2A products as long as the MDT xml file is present at the root of the product directory.

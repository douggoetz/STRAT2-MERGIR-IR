
from pydap.client import open_url
from pydap.cas.urs import setup_session
import pandas as pd
import numpy as np
import datetime
#import pprint
from math import radians, cos, sin, asin, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import imageio
import os

## Credentials for NASA EarthData Account
##follow this link for information on creating an account: https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+For+an+EarthData+Login+Profile
usernm = 'douggoetz'
userpass = 'jdgIl2l2mal!!'

##input csv file that contains gondola time, latitude, and longtitude without headers. Column 0 (time) should be in UTC time since 2019-1-1, latitude (column 1) and longitude (column 2) in decimal degrees
gond_pos_csv = "/Users/jago1642/Documents/Strat2.1/MERGIR/RATS_processed_30min.csv" 
colnames = ['Gondtime', 'Gondlat', 'Gondlon']
gond_data = pd.read_csv(gond_pos_csv, header=None, names = colnames, dtype = np.float64) #make a pandas dataframe that contains the gondola data

##make directory for plotting
cwd = os.getcwd()
if not os.path.isdir(cwd+'/MergirImages'):
    os.mkdir(cwd+'/MergirImages')
    
path = cwd+'/MergirImages'

###global parameters for file processing

## +/- number of lat and lon degrees from rounded gondola postion (dd.dd) to refine Mergir image download for analysis 
LatBoundarySize = 10   # 1° is about 111 km and each pixel is 4km               
LonBoundarySize = 10   # 1° is about 111 km and each pixel is 4km 

## Spatial values from Mergir readme file
Mergirstartlat = -59.9818
deltalat = 0.036383683
Mergirstartlon = -179.9818 
deltalon = 0.036378335

BTthreshold = 235 #brightness temperature for convection from Corcos et al, 2021, JGR Atmospheres
BTdomain = 10  #+/- km from gondola position for averaging domain of BT below gondola


def DatetimeStr(rowpnt):  #this creates a tuples of datetime strings that are part of the NASA GESDISC OpenDAP URL

    Eurostime = gond_data.Gondtime[rowpnt]
    starttime = 1546300800 #UNIX time for 1970-1-1
    dt = datetime.datetime.fromtimestamp(starttime+Eurostime)
    dstruct = dt.timetuple() #makes struct with indexed datetime info
    datestr = str(dstruct[0])+str(dstruct[1]).zfill(2)+str(dstruct[2]).zfill(2)+str(dstruct[3]).zfill(2) #this is the date/day/hour string that is used in the .nc file name
    yearday = str(dstruct[7]) #this is the year day string that is used in the .nc file name

    halfhourflag = 1

    if(dstruct[4]<30): #determine if the first half or second half of the hour should be used
        halfhourflag = 0

    return[datestr, yearday, halfhourflag]

def FindImageBoundaries(rowpnt): #this creates a tuple of index values that are used to download subsections of the full MERGIR image

    ##Round gondola position to 2 decimal places
    Roundlat = round(gond_data.Gondlat[rowpnt], 2) 
    Roundlon = round(gond_data.Gondlon[rowpnt], 2)

    if pd.isna(Roundlat) or pd.isna(Roundlon): ## if either of the gondola values are NAN
    
        return[-1, 0, 0, 0]
        
    else:

        Imagemaxlat = Roundlat + LatBoundarySize
        Imageminlat = Roundlat - LatBoundarySize

        Imagemaxlon = Roundlon + LonBoundarySize
        Imageminlon = Roundlon - LonBoundarySize

        if Imagemaxlat > 60.00 or Imageminlat < -60.00:
            print("refined image boundary is outside of Mergir latitude range")
            print("for gondola position index point" + str(rowpnt))
            exit  ## to do: figure out how to handle edge cases automatically

        if Imagemaxlon > 180.00 or Imageminlon <-180.00:

            print("refined image boundary is outside of Mergir longitude range")
            print("for gondola position index point " + str(rowpnt))
            exit  ## to do: figure out how to handle boundary edges automatically

        MergirMinRow = round((Imageminlat - Mergirstartlat)/deltalat)
        MergirMaxRow = round((Imagemaxlat - Mergirstartlat)/deltalat)
        MergirMinCol = round((Imageminlon - Mergirstartlon)/deltalon)
        MergirMaxCol = round((Imagemaxlon - Mergirstartlon)/deltalon)

        return[MergirMinRow,MergirMaxRow,MergirMinCol,MergirMaxCol]

def DistanceFromBTvalue(bt_arr, poslist, gondlat, gondlon): #this calculates the distace from the gondola to the nearest convection <BTthreshold

    dist_array = np.zeros_like(bt_arr)

    with np.nditer(dist_array, flags=['multi_index'], op_flags=['writeonly']) as it:
        while not it.finished:
            lat = deltalat*(int(poslist[0])+it.multi_index[1])+Mergirstartlat #find the latitude for the given position
            lon = deltalon*(int(poslist[2])+it.multi_index[2])+Mergirstartlon #find the latitude for the given position
            it[0] = haversine(lat,gondlat,lon, gondlon)                    
            is_not_finished = it.iternext()

    locationlist = dist_array[np.where(bt_arr < BTthreshold)]
    
    if locationlist.shape[0] == 0:
        return 0
    else: 
        return np.amin(locationlist)
    
def haversine(lat1,lat2,lon1,lon2): #great circle distance calculator in km

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6378.1 # Radius of earth in kilometers.
    #print (c*r)
    return c * r

def BTBelowGondola(bt_arr, poslist, gondlat, gondlon): #this calculates the average BT from the MERGIR IR image within the BTdomain distance from the gondola position

    dist_array = np.zeros_like(bt_arr)

    with np.nditer(dist_array, flags=['multi_index'], op_flags=['writeonly']) as it:
        while not it.finished:
            lat = deltalat*(int(poslist[0])+it.multi_index[1])+Mergirstartlat #find the latitude for the given position
            lon = deltalon*(int(poslist[2])+it.multi_index[2])+Mergirstartlon #find the latitude for the given position
            it[0] = haversine(lat,gondlat,lon, gondlon)                    
            is_not_finished = it.iternext()

    locationlist = bt_arr[np.where(dist_array < BTdomain)]
    #print(np.nanmean(locationlist))
    return np.nanmean(locationlist)

def MergirImageAnalysis(): # this runs the MERGIR IR download and file processing using pydap

    BTdist2gond = np.zeros_like(gond_data.Gondtime)
    BTbelow = np.zeros_like(gond_data.Gondtime)
    URLnm = np.empty_like(gond_data.Gondtime,dtype=object)

    for i in gond_data.Gondtime.index: #gond_data.shape[0]:

        timelist = DatetimeStr(i)
        boundlist = FindImageBoundaries(i)

        if boundlist[0] == -1:
            print("NAN found in gondola dataset at point " + str(i))

        else:
            datestr = timelist[0]
            year = datestr[0:4]
            short_url = 'https://disc2.gesdisc.eosdis.nasa.gov/opendap/MERGED_IR/GPM_MERGIR.1/'+year+'/'+timelist[1]+'/merg_'+datestr+'_4km-pixel.nc4'
            mergir_url = short_url+'?Tb['+str(timelist[2])+':1:'+str(timelist[2])+']['+str(boundlist[0])+':1:'+str(boundlist[1])+']['+str(boundlist[2])+':1:'+str(boundlist[3])+']'

            
            #example of long url with bounds :https://disc2.gesdisc.eosdis.nasa.gov/opendap/MERGED_IR/GPM_MERGIR.1/2021/326/merg_2021112211_4km-pixel.nc4?Tb[1:1:1][1245:1:1795][6199:1:6749]

            if i == 0:
    
                session = setup_session(usernm, userpass, check_url=short_url) ##user should change login credentials
                dataset = open_url(short_url, session=session)
                Btarr = dataset['Tb'] #calls to brightness temperature matrix found in OpenDAP
                print (mergir_url)
                grid = Btarr[0,boundlist[0]:boundlist[1],boundlist[2]:boundlist[3]] #this downloads the data within new bounding box from OpenDAP to memory
                #pprint.pprint(grid.attributes) ##to print attributes for debugging
                #print(grid) #for debugging
                gridshp = grid.shape 
                #print(grid.data[0,0:gridshp[1],0:gridshp[2]])

                last_timelist = timelist
                last_boundlist = boundlist

            else:

                if last_timelist != timelist or last_boundlist != boundlist: #check to see if new data needs to be downloaded or if the previous grid data will work

                    session = setup_session(usernm, userpass, check_url=short_url) ##user should change login credentials
                    dataset = open_url(short_url, session=session)
                    Btarr = dataset['Tb'] #calls to brightness temperature matrix found in OpenDAP
                    print (mergir_url)
                    grid = Btarr[0,boundlist[0]:boundlist[1],boundlist[2]:boundlist[3]] #this downloads the data within new bounding box from OpenDAP to memory
                    #pprint.pprint(grid.attributes) ##to print attributes for debugging
                    gridshp = grid.shape
                    #print(grid.data[0,0:gridshp[1],0:gridshp[2]])

                    last_timelist = timelist
                    last_boundlist = boundlist
                 
        BTdist2gond[i] = DistanceFromBTvalue(grid, boundlist, gond_data.Gondlat[i], gond_data.Gondlon[i])
        BTbelow[i] = BTBelowGondola(grid, boundlist, gond_data.Gondlat[i], gond_data.Gondlon[i])
        URLnm[i] = mergir_url

        saveimage(grid,boundlist,timelist,i)

    gond_data['GondDist'] = BTdist2gond.tolist()
    gond_data['GondBelow'] = BTbelow.tolist()
    gond_data['FileURL'] = URLnm.tolist()

    gond_data.to_csv('MergirOutput.csv')
    
def saveimage(array,blist,tlist,index): #this function generates plots of the IR imagery from the MERGIR downloads

    lat1 = deltalat*(int(blist[0]))+ Mergirstartlat 
    lat2 = deltalat*(int(blist[1]))+Mergirstartlat
    long1 = deltalon*(int(blist[2]))+Mergirstartlon
    long2 = deltalon*(int(blist[3]))+Mergirstartlon

    arrshape = np.shape(array)
    #print (arrshape)
    latlist = np.empty(arrshape[1],dtype=float)
    lonlist = np.empty(arrshape[2],dtype=float)

    for i in range(len(latlist)):
        latlist[i] = lat1 + (deltalat*i)

    for i in range(len(lonlist)):
        lonlist[i] = long1 +(deltalon*i)   

    m = Basemap(projection='cyl',resolution='l',llcrnrlon =long1,llcrnrlat=lat1,urcrnrlon=long2 ,urcrnrlat=lat2)
    m.drawcoastlines()
    m.drawcountries(linewidth=1.8)
    m.drawmapboundary(fill_color='#99ffff')
    m.drawparallels(np.arange(-60,60,5),labels=[1,1,0,1])  ##configure this to change latitude axis 
    m.drawmeridians(np.arange(-180,180,5),labels=[1,1,0,1]) ##configure this to change longitude axis

    m.plot(gond_data.Gondlon[index],gond_data.Gondlat[index], marker="o", color="green", ls="") #this is the marker for the gondola position

    m.pcolormesh(lonlist, latlist, array[0][:][:],latlon=True, cmap='RdBu_r', shading='auto') ## configure this to change the color scheme of image
    plt.colorbar(label='Bt',shrink = 0.75)
    plt.clim(BTthreshold, 315) ##modify this to adjust z-scaling
    plt.title(tlist[0]+'_'+'Tb['+str(tlist[2])+':1:'+str(tlist[2])+']['+str(blist[0])+':1:'+str(blist[1])+']['+str(blist[2])+':1:'+str(blist[3])+']')
    
    fignm = 'Mergir'+tlist[0]+'_'+'Tb['+str(tlist[2])+':1:'+str(tlist[2])+']['+str(blist[0])+':1:'+str(blist[1])+']['+str(blist[2])+':1:'+str(blist[3])+']'
    plt.savefig(path+'/'+fignm) #Saves each figure as an image
    plt.close()
        
def main():
  # main function that you use to run the program
    MergirImageAnalysis()   
   
if __name__ == '__main__':
    main()


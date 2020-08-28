import glob
from netCDF4 import Dataset
from numpy import *
fmt='https://pmm-gv.gsfc.nasa.gov/pub/gpm-validation/data/gpmgv/orbit_subset/GPM/DPR/2ADPR/V06A/CONUS/%4.4i/%2.2i/%2.2i'
from bs4 import BeautifulSoup

import urllib.request
import os
f2ADPR=glob.glob('2A-DPR/2A-C*')
dFile={}
dOrbs=[]
for f in f2ADPR:
    orb=f.split('.')[5][1:]
    dFile[orb]=f
    dOrbs.append(orb)

zKuL=[[],[],[],[]]
zKaL=[[],[],[],[]]
clutFL=[[],[],[],[]]
reliabFlagNSL=[[],[],[],[]]
reliabFlagMSL=[[],[],[],[]]
pathAttenNSL=[[],[],[],[]]
pathAttenMSL=[[],[],[],[]]
binSfcNSL=[[],[],[],[]]
gvSfcNSL=[[],[],[],[]]
binSTopL=[[],[],[],[]]
binZeroDegL=[[],[],[],[]]
dprWCL=[[],[],[],[]]
def readGRCMB(fname):
    f=Dataset(fname,'r')
    dprRrate=f.variables['precipTotRate_MS'][:,:]
    grRrate=f.variables['GR_RP_rainrate_MS'][:,:]
    gvZ=f.variables['GR_Z_MS'][:,:]
    cmbZ=f.variables['correctedReflectFactor_MS'][:,:]
    lat=f.variables['DPRlatitude_MS'][:] 
    lon=f.variables['DPRlongitude_MS'][:]
    gr_Nw=f.variables['GR_Nw_MS'][:,:]
    zeroDegH=f.variables['zeroDegAltitude_MS'][:]
    bHeight=f.variables['bottomHeight_MS'][:,:]
    tHeight=f.variables['topHeight_MS'][:,:]
    dm=f.variables['GR_Dm_MS'][:,:]
    stormH=f.variables['stormTopAltitude_MS'][:,0]
    pType=f["precipitationType_MS"][:]
    dprWC=f["precipTotWaterCont_MS"][:,:]

    return dprRrate,grRrate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm,stormH,pType,cmbZ, dprWC

fs=glob.glob("data/*K*nc")
s1=[0,0,0,0]
s2=[0,0,0,0]
s1L=[[],[],[],[]]
s2L=[[],[],[],[]]
cmbSfcRainL=[[],[],[],[]]
gvSfcRainL=[[],[],[],[]]
orbL=[]
iM=0
s1M=[0,0,0,0]
s2M=[0,0,0,0]
for l in fs[:]:
    dprRate,grRate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm,stormH,pType,cmbZ,dprWC=readGRCMB(l)
    fh=Dataset(l)
    #stop
    lxss=l.split('.')
    orb=l[3]
    yy=fh['Year_MS'][0]
    mm=fh['Month_MS'][0]
    dd=fh['DayOfMonth_MS'][0]

    pType=(pType/1e7).astype(int)
    a=nonzero(zeroDegH>3750)
    b=nonzero(pType[a]==2)
    c=nonzero(stormH[a][b]<12500)
    iabmax=cmbZ.shape[0]
    ray=fh['rayNum_MS'][:]
    scan=fh['scanNum_MS'][:]
    extL=[]
    if len(c[0])>2:
        idown=0
        for i in a[0][b][c]:
            iab=1
            if bHeight[iab,i]<zeroDegH[i] and iab<iabmax-1:
                iab+=1
            if bHeight[iab,i]<zeroDegH[i] and iab<iabmax-1:
                iab+=1
            if bHeight[iab,i]<zeroDegH[i] and iab<iabmax-1:
                iab+=1
            if dprRate[0,i]>-0.1 and grRate[0,i]>-0.1:
                if stormH[i]<6500:
                    ic=0
                else:
                    if stormH[i]<12000:
                        if cmbZ[iab:,i].max()<35:
                            ic=1
                        else:
                            ic=2
                    else:
                        ic=3
                if ic<4:
                    extL.append((lon[i],lat[i],dprRate[0,i],\
                                     grRate[0,i],dprWC[0,i],ray[i],scan[i],ic))
                s1[ic]+=dprRate[0,i]
                s2[ic]+=grRate[0,i]
                s1L[ic].append(dprRate[0,i])
                s2L[ic].append(grRate[0,i])
        if orb not in orbL:
            orbL.append(orb)
        if orb in dOrbs:
            print(dFile[orb])
            fh=Dataset(dFile[orb])
            dLon=fh['MS/Longitude'][:,:]
            dLat=fh['MS/Latitude'][:,:]
            reliabFlag=fh['MS/SRT/reliabFlag'][:,:]  
            nt=reliabFLag.shape[0]
            for rec in extL:
                ind=argmin((rec[0]-dLon)**2+(rec[1]-dLat)**2)
                i0=int(ind/25)
                j0=ind-i0*25
                if (rec[0]-dLon[i0,j0])**2+(rec[1]-dLat[i0,j0])**2<1e-3:
                    iM+=1
                    if reliabFlag[i0,j0]==1 or reliabFlag[i0,j0]==2 and i0>0 and i0<nt-1:
                        ic=rec[-1]
                        s1M[ic]+=rec[2]
                        s2M[ic]+=rec[3]
                        zKu1=fh['NS/PRE/zFactorMeasured'][i0-1:i0+2,j0+11:j0+14,:]
                        zKa1=fh['NS/PRE/zFactorMeasured'][i0,j0,:]
                        clutF1=fh['NS/PRE/binClutterFreeBottom'][i0,j0+12]
                        reliabFlagNS1=fh['NS/SRT/reliabFlag'][i0,j0+12]
                        reliabFlagMS1=fh['MS/SRT/reliabFlag'][i0,j0]
                        pathAttenNS1=fh['NS/SRT/pathAtten'][i0,j0+12]
                        pathAttenMS1=fh['MS/SRT/pathAtten'][i0,j0]
                        binSfcNS1=fh['NS/PRE/binRealSurface'][i0,j0+12]
                        binSTop1=fh['NS/PRE/binStormTop'][i0,j0+12]
                        binZeroDeg1=fh['NS/VER/binZeroDeg'][i0,j0+12]
                        clutFL[ic].append(clutF1)
                        reliabFlagNSL[ic].append(reliabFlagNS1)
                        reliabFlagMSL[ic].append(reliabFlagMS1)
                        pathAttenNSL[ic].append(pathAttenNS1)
                        pathAttenMSL[ic].append(pathAttenMS1)
                        binSfcNSL[ic].append(binSfcNS1)
                        binSTopL[ic].append(binSTop1)
                        binZeroDegL[ic].append(binZeroDeg1)
                        zKuL[ic].append(zKu1)
                        zKaL[ic].append(zKa1)
                        cmbSfcRainL[ic].append(rec[2])
                        gvSfcRainL[ic].append(rec[3])
                        dprWCL[ic].append(rec[4])
                        
        if len(extL)>1 and 1==0:
            response = urllib.request.urlopen(fmt%(yy,mm,dd))
            html = response.read()
            soup = BeautifulSoup(html)
            text=soup.get_text()
            texts=text.split('\n')
            for fs in texts:
                if fs!='':
                    fss=fs.split()[0]
                    print(fss, orb in fss)
                    if '2A-CS' in fss and orb in fss:
                        cmd='wget -nc '+fmt%(yy,mm,dd)+'/'+fss
                        os.system(cmd)
            #stop


import matplotlib.pyplot as plt
plt.bar(arange(4)*2,s1M)
plt.bar(arange(4)*2+1,s2M)
plt.legend(['CMB','GV-VN'])
plt.xticks(2*arange(4)+0.5,['ETop<6.5','Zmax<35', 'Zmax>35','12<ETop'])
plt.ylabel('Rain Volume')
plt.savefig('CMB_VN_Comparison_CNV_35_MS_RelFact<=2.png')

import pickle

dataOut=[clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,\
    zKuL,zKaL,cmbSfcRainL,gvSfcRainL,dprWCL]

pickle.dump(dataOut,open('dataVN_3D.pklz','wb'))

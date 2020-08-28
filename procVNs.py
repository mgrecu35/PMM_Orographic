import glob
from netCDF4 import Dataset
from numpy import *
def readGRCMB(fname):
    f=Dataset(fname,'r')
    dprRrate=f.variables['precipTotRate_NS'][:,:]
    grRrate=f.variables['GR_RP_rainrate_NS'][:,:]
    gvZ=f.variables['GR_Z_NS'][:,:]
    cmbZ=f.variables['correctedReflectFactor_NS'][:,:]
    lat=f.variables['DPRlatitude_NS'][:] 
    lon=f.variables['DPRlongitude_NS'][:]
    gr_Nw=f.variables['GR_Nw_NS'][:,:]
    zeroDegH=f.variables['zeroDegAltitude_NS'][:]
    bHeight=f.variables['bottomHeight_NS'][:,:]
    tHeight=f.variables['topHeight_NS'][:,:]
    dm=f.variables['GR_Dm_NS'][:,:]
    stormH=f.variables['stormTopAltitude_NS'][:]
    pType=f["precipitationType_NS"][:]

    return dprRrate,grRrate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm,stormH,pType,cmbZ

fs=glob.glob("data/*nc")
s1=[0,0,0,0]
s2=[0,0,0,0]
s1L=[[],[],[],[]]
s2L=[[],[],[],[]]
for l in fs[:]:
    dprRate,grRate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm,stormH,pType,cmbZ=readGRCMB(l)
    pType=(pType/1e7).astype(int)
    a=nonzero(zeroDegH>3750)
    b=nonzero(pType[a]==1)
    c=nonzero(stormH[a][b]<12500)
    iabmax=cmbZ.shape[0]
    if len(c[0])>2:
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
                s1[ic]+=dprRate[0,i]
                s2[ic]+=grRate[0,i]
                s1L[ic].append(dprRate[0,i])
                s2L[ic].append(grRate[0,i])
        #print(s1/s2, l)

import matplotlib.pyplot as plt
plt.bar(arange(4)*2,s1)
plt.bar(arange(4)*2+1,s2)
plt.legend(['CMB','GV-VN'])
plt.xticks(2*arange(4)+0.5,['ETop<6.5','Zmax<35', 'Zmax>35','12<ETop'])
plt.ylabel('Rain Volume')
plt.savefig('CMB_VN_Comparison_CNV_35.png')

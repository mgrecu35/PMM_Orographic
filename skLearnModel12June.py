from numpy import  *
from netCDF4 import Dataset
fname='dBase12June2014.nc'
nc=Dataset(fname,'r')
zProf=nc['zProf'][:,:]/10.
zProf[zProf<0]=0.
woro=-nc['ww'][:,:,:]
zProf=zProf[:,:]
dbz=nc['dbz'][:,:,:]/10.
dbz[dbz<0]=0
dbz0=nc['dbz0'][:,:,:]/10.
dbz0=dbz0[:,:,:]
dbz0[dbz0<0]=0
zsfc=nc['zsfc'][:,:,:]/1e3
zsfc=zsfc[:,:,:]
hagl=nc['hagl'][:,:,:]/1e3
hagl=hagl[:,:,:]

nt,nw,nw=zsfc.shape
zsfc2d=zsfc.reshape(nt,nw*nw)
woro2d=woro.reshape(nt,nw*nw)
#woro2d=hagl.reshape(nt,nw*nw)
from sklearn.cluster import MiniBatchKMeans, KMeans

k_means = MiniBatchKMeans(init='k-means++', n_clusters=49, batch_size=5000, random_state=0)
import pickle
#k_means=pickle.load(open('kmeansZ.pklz','rb'))
#labels_=k_means.predict(zsfc2d)
k_means.fit(zsfc2d)
labels_=k_means.labels_
it=0
import matplotlib.pyplot as plt
fig=plt.figure()
zsfcL=[]
for i in range(7):
    for j in range(7):
        ax=plt.subplot(7,7,it+1)
        a=nonzero(labels_==it)
        zsfcm=zsfc[a[0],:,:].mean(axis=0)
        zsfcL.extend(list(zsfcm.reshape(49)))
        im=plt.pcolormesh(zsfcm,vmax=1.5,cmap='jet',vmin=0.5)
        it+=1
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.savefig('terrain2.png')
fig=plt.figure()
it=0
wmL=[]
for i in range(7):
    for j in range(7):
        ax=plt.subplot(7,7,it+1)
        a=nonzero(labels_==it)
        #print(len(a[0]))
        zm1=10*dbz[a[0],:,:].mean(axis=0)
        zm2=10*dbz0[a[0],:,:].mean(axis=0)
        wm=woro[a[0],:,:].mean(axis=0)
        wmL.extend(list(wm.reshape(49)))
        im=plt.pcolormesh(wm,vmax=0.35,vmin=-0.35,cmap='RdBu_r')
        it+=1
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.savefig('woro2.png')
fig=plt.figure()
it=0
for i in range(7):
    for j in range(7):
        ax=plt.subplot(7,7,it+1)
        a=nonzero(labels_==it)
        #print(len(a[0]))
        zm1=10*dbz[a[0],:,:].mean(axis=0)
        zm2=10*dbz0[a[0],:,:].mean(axis=0)
        zm11=zeros((7,7),float)
        zm12=zeros((7,7),float)
        c=0
        for i1 in a[0]:
            a2=nonzero(dbz[i1,:,:]>0.00)
            if len(a2[0])>1:
                zm11+=10*(dbz[i1,:,:])
                zm12+=10*(dbz0[i1,:,:])
                c+=1
        zm11/=c
        zm12/=c
        wm=woro[a[0],:,:].mean(axis=0)
        im=plt.pcolormesh(zm11,vmax=12,vmin=4,cmap='jet')
        it+=1
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.savefig('averageZ2.png')
fig=plt.figure()
it=0
zratioL=[]
for i in range(7):
    for j in range(7):
        ax=plt.subplot(7,7,it+1)
        a=nonzero(labels_==it)
        #print(len(a[0]))
        zm1=10*dbz[a[0],:,:].mean(axis=0)
        zm2=10*dbz0[a[0],:,:].mean(axis=0)
        zm11=zeros((7,7),float)
        zm12=zeros((7,7),float)
        c=0
        for i1 in a[0]:
            a2=nonzero(dbz[i1,:,:]>0.00)
            if len(a2[0])>1:
                zm11+=10*(dbz[i1,:,:])
                zm12+=10*(dbz0[i1,:,:])
                c+=1
        zm11/=c
        zm12/=c
        wm=woro[a[0],:,:].mean(axis=0)
        zratio=zm12/zm11/0.84
        zratioL.extend(list(zratio.reshape(49)))
        im=plt.pcolormesh(zratio,vmax=1.25,vmin=0.75,cmap='RdBu_r')
        it+=1
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.savefig('zRatio2.png')

fig=plt.figure(figsize=(12,6))
plt.subplot(121)
plt.scatter(zsfcL,zratioL)
plt.xlabel('Surface height')
plt.ylabel('Z enhancement')
plt.subplot(122)
plt.scatter(wmL,zratioL)
plt.xlabel('Orographic velocity')
#plt.ylabel('Z enhancement')
plt.savefig('scatterPlot2.png')

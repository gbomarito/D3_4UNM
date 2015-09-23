import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import random

from Dislocation import Dislocation


class DislocationMngr(object):
    """A manager that takes care of disloaction storage, querying, and steping
    
    Attributes:
        dislocations: list of dislocation objects <list of Dislocaton>
        nu: poisson coefficient <float>
        mu: lame' parameter <float>
        drag: drag coefficient <float>
    """
    
    def __init__(self, n, m,d, dnum=20, b=1.0, area=1.0):
        """Returns an initialized dislocation object"""
        self.nu=n
        self.mu=m
        self.drag=d
        
        self.dislocations=[]
        #self.dislocations.append(Dislocation( np.array((0.,0.,0.)), np.array((1.,0.,0.)), np.array((0.,1.,0.)) ))
        #self.dislocations.append(Dislocation( np.array((1.,0.,0.)), np.array((-1.,0.,0.)), np.array((0.,1.,0.)) ))
        
        for i in range(dnum):
            if random.random()>0.5:
                burgers=np.array((b,0.,0.))
            else:
                burgers=np.array((-b,0.,0.))
            loc=np.array((-area+random.random()*2,-area+random.random()*2,-area+random.random()*2))
            self.dislocations.append(Dislocation( loc, burgers, np.array((0.,1.,0.)) ))
        

    def dd_step(self,dt,sig_ff=None,FE_results=None):
        """Increments dislocations based on a time step (dt)
        optionally finite element results or a far field stress can be added"""
        
        dnum=len(self.dislocations)
        
        #calculate stresses at each dislocation point
        sigma=[np.zeros((3,3))]*dnum
        for i in range(dnum):
            d_i=self.dislocations[i]
        
            #stress cused by all dislocations
            for d_j in self.dislocations:
                sigma[i] += d_j.stress_at_point(d_i.X,self.mu,self.nu)
                
            #far field stress
            if sig_ff is not None:
                sigma[i] += sig_ff
                
            #TODO add FE stress field    
        
        #move dislocations
        for i in range(dnum):
            self.dislocations[i].move(sigma[i],dt,self.drag)


    def stress_at_point(self,Y):
        """calculates the stress at point Y caused by all dislocations"""
        
        sigma=np.zeros((3,3))
        for d_i in self.dislocations:
            sigma += d_i.stress_at_point(Y,self.mu,self.nu)
            
        return sigma


    def disp_at_point(self,Y):
        """calculates the displacement at point Y caused by all dislocations"""
        
        disp=np.zeros(3)
        for d_i in self.dislocations:
            sigma += d_i.disp_at_point(Y,self.mu,self.nu)
            
        return disp


    def distortion_at_point(self,Y):
        """calculates the displacement at point Y caused by all dislocations"""
        
        distortion=np.zeros(3)
        for d_i in self.dislocations:
            distortion += d_i.distortion_at_point(Y,self.mu,self.nu)
            
        return distortion


    def dump(self):
        """dumps the location of all dislocations"""
        
        count=0
        for d_i in self.dislocations:
            print count,": ",d_i.X
            count+=1


    def plot(self,filename):
        """plot of all dislocations is output to a file"""
        
        x=[]
        y=[]
        theta=[]
        
        for d in self.dislocations:
            x.append(d.X[0])
            y.append(d.X[1])
            theta=1
            
        plt.figure()   
        plt.scatter(x,y,marker=r'$\perp$')
        pylab.savefig(filename+".png", bbox_inches='tight')


    def plot_w_stress(self,filename,sig_ff=None,FE_results=None):
        """plot of all dislocations with a stress contour"""
        
        xmax=2.
        ymax=1.
        xmin=-1.
        ymin=-1.
        #find stress field
        delta = 0.05
        sx = np.arange(-1.0, 2.0, delta)
        sy = np.arange(-1.0, 2.0, delta)
        sX, sY = np.meshgrid(sx, sy)
        pts=np.hstack( (sX.reshape(-1,1), sX.reshape(-1,1)) )
        sig_xx=np.zeros(sX.shape)
        sig_yy=np.zeros(sX.shape)
        sig_zz=np.zeros(sX.shape)
        sig_xy=np.zeros(sX.shape)
        sig_xz=np.zeros(sX.shape)
        sig_yz=np.zeros(sX.shape)
        for j in range(len(sx)):
            for i in range(len(sy)):
                loc=np.array((sX[i,j],sY[i,j],0.))
                
                sig=self.stress_at_point(loc) #disloc stress
                
                if sig_ff is not None:          #far field
                    sigma[i] += sig_ff
                
                sig_xx[i,j]=sig[0,0]
                sig_yy[i,j]=sig[1,1]
                sig_zz[i,j]=sig[2,2]
                sig_xy[i,j]=sig[0,1]
                sig_xz[i,j]=sig[0,2]
                sig_yz[i,j]=sig[1,2]
        
        
        plt.figure()        
        plt.contourf(sX,sY,sig_xx, levels=np.linspace(np.min(sig_xx),np.max(sig_xx),100) )
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xx.png", bbox_inches='tight')
        plt.close()
        
        plt.figure()        
        plt.contourf(sX,sY,sig_yy,levels=np.linspace(np.min(sig_yy),np.max(sig_yy),100) )
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_yy.png", bbox_inches='tight')
        plt.close()
        
        plt.figure()        
        plt.contourf(sX,sY,sig_zz,levels=np.linspace(np.min(sig_zz),np.max(sig_zz),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_zz.png", bbox_inches='tight')
        plt.close()
        
        plt.figure()        
        plt.contourf(sX,sY,sig_xy,levels=np.linspace(np.min(sig_xy),np.max(sig_xy),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xy.png", bbox_inches='tight')
        plt.close()
        
        '''
        plt.figure()        
        plt.contourf(sX,sY,sig_xz,levels=np.linspace(np.min(sig_xz),np.max(sig_xz),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xz.png", bbox_inches='tight')
        
        plt.figure()        
        plt.contourf(sX,sY,sig_yz,levels=np.linspace(np.min(sig_yz),np.max(sig_yz),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_yz.png", bbox_inches='tight')
        '''
        
        

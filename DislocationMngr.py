import numpy as np
import matplotlib.pyplot as plt
import pylab

from Dislocation import Dislocation


class DislocationMngr(object):
    """A manager that takes care of disloaction storage, querying, and steping
    
    Attributes:
        dislocations: list of dislocation objects <list of Dislocaton>
        nu: poisson coefficient <float>
        mu: lame' parameter <float>
        drag: drag coefficient <float>
    """
    
    def __init__(self, n, m,d):
        """Returns an initialized dislocation object"""
        self.nu=n
        self.mu=m
        self.drag=d
        
        self.dislocations=[]
        self.dislocations.append(Dislocation( np.array((0.,0.,0.)), np.array((1.,0.,0.)), np.array((1.,0.,0.)) ))
        self.dislocations.append(Dislocation( np.array((1.,0.,0.)), np.array((-1.,0.,0.)), np.array((1.,0.,0.)) ))
        

    def dd_step(self,dt,sig_ff=None,FE_results=None,):
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


    def plot(self,filename):
        """plot of all dislocations is output to a file"""
        
        x=[]
        y=[]
        theta=[]
        
        for d in self.dislocations:
            x.append(d.X[0])
            y.append(d.X[1])
            theta=1
            
            
        plt.scatter(x,y,marker=r'$\perp$')
        pylab.savefig(filename, bbox_inches='tight')
        
        

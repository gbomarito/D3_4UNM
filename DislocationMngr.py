import numpy as np
import Dislocation

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
        
        dislocations=[]
        dislocations.append(Dislocation(np.array((0.,0.,0.)), np.array((1.,0.,0.))), np.array((1.,0.,0.))))
        dislocations.append(Dislocation(np.array((1.,0.,0.)), np.array((-1.,0.,0.))), np.array((1.,0.,0.))))
        

    def dd_step(self,dt,sig_ff=None,FE_results=None,):
        """Increments dislocations based on a time step (dt)
        optionally finite element results or a far field stress can be added"""
        
        dnum=len(dislocations)
        
        #calculate stresses at each dislocation point
        sigma=[np.zeros((3,3))]*dnum
        for i in range(dnum):
            d_i=silocations[i]
        
            #stress cused by all dislocations
            for d_j in disloactions:
                sigma[i] += d_j.stress_at_point(d_i.X)
                
            #far field stress
            if sig_ff is not None:
                sigma[i] += sig_ff
                
            #TODO add FE stress field    
        
        #move dislocations
        for i in range(dnum):
            dislocations[i].move(sigma[i],dt,drag)


    def stress_at_point(self,Y)
        """calculates the stress at point Y caused by all dislocations"""
        
        sigma=np.zeros((3,3))
        for d_i in disloactions:
            sigma += d_i.stress_at_point(Y)
            
        return sigma


    def disp_at_point(self,Y)
        """calculates the displacement at point Y caused by all dislocations"""
        
        disp=np.zeros(3)
        for d_i in disloactions:
            sigma += d_i.disp_at_point(Y)
            
        return disp


    def distortion_at_point(self,Y)
        """calculates the displacement at point Y caused by all dislocations"""
        
        distortion=np.zeros(3)
        for d_i in disloactions:
            distortion += d_i.distortion_at_point(Y)
            
        return distortion
        

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
        dislocations.append(Dislocation(np.array((0,0,0)), np.array((1,0,0))), np.array((1,0,0))))
        

    def dd_step(self,dt,FE_results=None,sig_ff=None):
        """Increments dislocations based on a time step (dt)
        optionally finite element results or a far field stress can be added"""
        
    
    def stress_at_point(self,Y)
        """calculates the stress at point Y caused by all dislocations"""
    
    def disp_at_point(self,Y)
        """calculates the displacement at point Y caused by all dislocations"""
        

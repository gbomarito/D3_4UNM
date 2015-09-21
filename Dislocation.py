import numpy as np

class Dislocation(object):
    """A dislocation
    
    Attributes:
        X: location of dislocation <numpy vector>
        burgers: burgers vector <numpy vector>
        line_vec: line vector <numpy vector>
        slip_plane: slip plane <numpy vector>
    """
    
    def __init__(self, location, bv, sp):
        """Returns an initialized dislocation object"""
        self.X=np.copy(location)
        self.burgers=np.copy(bv)
        self.line_vec=np.array((0,0,1))
        self.slip_plane=np.copy(sp)
        

    def stress_at_point(self,Y,mu,nu):
        """calculates the stress at point Y caused by a dislocation (self) 
        in a material with elastic parameters nu and mu"""
        
        return sigma


    def disp_at_point(self,Y,mu,nu):
        """calculates the displacement at point Y caused by a dislocation (self) 
        in a material with elastic materparameters nu and mu"""
        
        return disp


    def distortion_at_point(self,Y,mu,nu):
        """calculates the distortion at point Y caused by a dislocation (self) 
        in a material with elastic parameters nu and mu"""
        
        return distortion


    def move(self,sigma,dt,drag):
        """moves disloaction self based on the stress (sigma), time step (dt), 
        and drag coeeficient (drag)"""
        
        

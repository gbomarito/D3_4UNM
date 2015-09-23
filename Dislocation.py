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
        self.slip_plane=np.copy(sp)/np.linalg.norm(sp)
        

    def stress_at_point(self,Y,mu,nu):
        """calculates the stress at point Y caused by a dislocation (self) 
        in a material with elastic parameters nu and mu"""
        sigma = self.stress_screw(Y,mu,nu) + self.stress_edge(Y,mu,nu)
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
        """moves dislocation (self) based on the stress (sigma), time step (dt), 
        and drag coeeficient (drag)"""
        
        F = np.cross(np.dot(self.burgers,sigma),self.line_vec)
        self.X += (F - np.dot(F,self.slip_plane)*self.slip_plane)*dt/drag

    def get_velocity(self,sigma,drag)
        """ finds the velocity of dislocation (self) based on the stress (sigma)
        and drag"""
        
	    F = np.cross(np.dot(self.burgers,sigma),self.line_vec)
	    V = (F - np.dot(F,self.slip_plane)*self.slip_plane)/drag
	    
	    return V
        
    def stress_screw(self, Y, mu, nu):
        """finds stress caused by screw component of dislocation (self) at point
        Y in a material with elastic properties mu, nu"""
        
        bscrew = np.dot(self.burgers,self.line_vec)*self.line_vec
	    b = np.linalg.norm(bscrew)

        if np.abs(np.abs(self.line_vec[0])-1)>.0001:
        	yvec=np.cross(self.line_vec,np.array((1.,0.,0.)))
        else:
        	yvec=np.cross(self.line_vec,np.array((0.,1.,0.)))
        
        yvec = yvec/np.linalg.norm(yvec)
        xvec = np.cross(yvec,self.line_vec)
        xvec = xvec/np.linalg.norm(xvec)
        
        g = np.vstack((xvec,yvec,self.line_vec))
        g = np.transpose(g)
        
        r = np.dot(np.transpose(g),Y-self.X)
        
        rn = np.linalg.norm(r)
        
        sig = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        
        if np.abs(rn)>1e-8:
            sig[0,2] = -mu*b*r[1]/(2*np.pi*rn**2)
            sig[2,0] = sig[0,2]
            
            sig[1,2] = -mu*b*r[0]/(2*np.pi*rn**2)
            sig[2,1] = sig[1,2]
        
        return np.dot(np.dot(g,sig),np.transpose(g))

    def stress_edge(self, Y, mu, nu):
        """finds stress caused by edge component of dislocation (self) at point
        Y in a material with elastic properties mu, nu"""
        bedge = self.burgers - np.dot(self.burgers,self.line_vec)*self.line_vec
	    b = np.linalg.norm(bedge)
        
        g = np.vstack((bedge/b,np.cross(self.line_vec,bedge/b),self.line_vec))
        g = np.transpose(g)
        
        r = np.dot(np.transpose(g),Y-self.X)
        
        rn = np.linalg.norm(r)
        
        sig = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        
        if np.abs(rn)>1e-8:
            K = -mu*b/(2*np.pi*(1-nu))
            
            sig[0,0] = K*(r[1]/rn**4)*(r[1]**2 + 3*r[0]**2)
            sig[1,1] = K*(r[1]/rn**4)*(r[1]**2 - r[0]**2)
            sig[0,1] = -K*(r[0]/rn**4)*(r[0]**2 - r[1]**2)
            sig[1,0] = sig[1,2]
            sig[2,2] = K*2*nu*(r[1]**2/rn**2)
        
        return np.dot(np.dot(g,sig),np.transpose(g))

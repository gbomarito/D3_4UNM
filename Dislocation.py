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
        """Returns an initialized dislocation object.  Only dislocations 
           with line vectors in the z direction are possible."""
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
        disp = displacement_edge(Y,mu,nu) + displacement_screw(Y,mu,nu)
        return disp


    def distortion_at_point(self,Y,mu,nu):
        """calculates the distortion at point Y caused by a dislocation (self) 
        in a material with elastic parameters nu and mu"""
        distortion = distortion_screw(Y,mu,nu) + distortion_edge(Y,mu,nu)
        
        return distortion


    def move(self,sigma,dt,drag):
        """moves dislocation (self) based on the stress (sigma), time step (dt), 
        and drag coeeficient (drag).  """
        
        F = np.cross(np.dot(self.burgers,sigma),self.line_vec)
        self.X += (F - np.dot(F,self.slip_plane)*self.slip_plane)*dt/drag

    def move_verlot(self,sigma,dt,drag,V_old,dt_old):
        """moves dislocation (self) based on the stress (sigma), time step (dt), 
        and drag coeeficient (drag).  Verlot integration also requires the velocity and 
           time step from the previous iteration"""
        
        V = self.get_velocity(sigma,drag)
        dX = V*dt + .5*(V-V_old)*dt*dt/dt_old
        self.move_set(dX)

    def move_set(self,dX):
        self.X += dX

    def rotate(self,beta):
        Q = beta
        Q[0,2] = 0.
        Q[1,2] = 0.
        Q[2,2] = 1.
        b = np.linalg.norm(self.burgers)        
        self.burgers = Q*self.burgers
        self.burgers = b*self.burgers/np.linalg.norm(self.burgers)
        self.slip_plane = Q*self.slip_plane
        self.slip_plane = self.slip_plane/np.linalg.norm(self.slip_plane)

    def rotate_and_stretch(self,beta):
        """Proposed duplicate of rotate, but allowing the second order effect of changing the length of Burgers vector"""
        Q = beta
        Q[0,2] = 0.
        Q[1,2] = 0.
        Q[2,2] = 1.
        self.burgers = Q*self.burgers
        self.slip_plane = Q*self.slip_plane
        self.slip_plane = self.slip_plane/np.linalg.norm(self.slip_plane)

    def get_velocity(self,sigma,drag):
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
            C = -mu*b/(2*np.pi*(1-nu))
            
            sig[0,0] = C*(r[1]/rn**4)*(r[1]**2 + 3*r[0]**2)
            sig[1,1] = C*(r[1]/rn**4)*(r[1]**2 - r[0]**2)
            sig[0,1] = -C*(r[0]/rn**4)*(r[0]**2 - r[1]**2)
            sig[1,0] = sig[0,1]
            sig[2,2] = C*2*nu*(r[1]**2/rn**2)
        
        return np.dot(np.dot(g,sig),np.transpose(g))

    def distortion_screw(self, Y, mu, nu):
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
        
        beta = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        if np.abs(rn)>1e-8:
            beta[2,0] = -b*r[1]/(2.*np.pi*rn**2)
            beta[2,1] = b*r[0]/(2.*np.pi*rn**2)

        return np.dot(np.dot(g,beta),np.transpose(g))

    def distortion_edge(self,Y,mu,nu):
        bedge = self.burgers - np.dot(self.burgers,self.line_vec)*self.line_vec
        b = np.linalg.norm(bedge)
        
        g = np.vstack((bedge/b,np.cross(self.line_vec,bedge/b),self.line_vec))
        g = np.transpose(g)
        
        r = np.dot(np.transpose(g),Y-self.X)
        
        rn = np.linalg.norm(r)
        
        beta = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        
        if np.abs(rn)>1e-8:
            C = -b/(4.*np.pi*(1-nu)*rn**2)
            beta[0,0] = C*r[1]*((1.-nu) + 2.*r[0]**2/rn**2)
            beta[0,1] = -C*r[0]*((3.-2.*nu) - 2.*r[1]**2/rn**2)
            beta[1,0] = C*r[0]*((1-nu) + 2.*r[1]**2/rn**2)
            beta[1,1] = C*r[1]*((1.-nu) - 2.*r[0]**2/rn**2)
        
        return np.dot(np.dot(g,beta),np.transpose(g))

    def displacement_screw(self, Y, mu, nu):
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
        
        u = np.array([0.,0.,0.])
        if np.abs(rn)>1e-8:
            u[2] = -b*np.atan2(r[0],r[1])/np.pi/2
        return np.dot(g,u)

    def displacement_edge(self,Y,mu,nu):
        bedge = self.burgers - np.dot(self.burgers,self.line_vec)*self.line_vec
        b = np.linalg.norm(bedge)
        
        g = np.vstack((bedge/b,np.cross(self.line_vec,bedge/b),self.line_vec))
        g = np.transpose(g)
        
        r = np.dot(np.transpose(g),Y-self.X)
        
        rn = np.linalg.norm(r)
        
        beta = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        
        u = np.array([0.,0.,0.])
        if np.abs(rn)>1e-8:
            x = r[0]
            y = r[1]
            u[0] = (.5*b/np.pi)*(atan2(y,x) + x*y/((2. - 2.*nu)*rn**2.))
            u[1] = (-.5*b/np.pi)*(1/(4.-4.*nu))*(np.log(rn**2.) + (x*x - y*y)/rn**2.)
        return np.dot(g,u)


################################################################################
#  STATIC FUNCTIONS
def interaction_energy(A,B, dXa, dXb, mu, nu):
    Ra = B.X - A.X
    R = Ra + dXb - dXa
    Rmag = np.linalg.norm(R)
    Rhat = R/Rmag
    Ramag = np.linalg.norm(Ra)
    Rahat = Ra/Ramag
    ba = A.burgers
    bb = B.burgers
    xi = A.line_vec #Assumes the line vecs are parallel, only valid for 2D
"""
    dE = np.dot(ba,xi)*np.dot(bb,xi)*np.log(Rmag/Ramag)
    dE += (1/(1-nu))*np.dot(np.cross(ba,xi),np.cross(bb,xi))*np.log(Rmag/Ramag)
    dE += (1/(1-nu))*np.dot(np.cross(ba,xi),Rhat)*np.dot(np.cross(bb,xi),Rhat)
    dE += -(1/(1-nu))*np.dot(np.cross(ba,xi),Rahat)*np.dot(np.cross(bb,xi),Rahat)
    dE = -dE*mu/(2*np.pi)
"""
    dE = ba[2]*bb[2]*np.log(Rmag/Ramag)
    dE += (1/(1-nu))*(ba[1]*bb[1] + ba[0]*bb[0])*np.log(Rmag/Ramag)
    dE += (1/(1-nu))*(ba[1]*Rhat[0] - ba[0]*Rhat[1])*(bb[1]*Rhat[0] - bb[0]*Rhat[1])
    dE += -(1/(1-nu))*(ba[1]*Rahat[0] - ba[0]*Rahat[1])*(bb[1]*Rahat[0] - bb[0]*Rahat[1])
    dE = -dE*mu/(2*np.pi)
    
    return dE

def check_for_interaction(A,B, L_glide, L_climb):
    #Checks to see if B is in range of A.  Extinction threshold is an oblate ellipsoid
    r = B.X - A.X
    r_off_plane = np.dot(r,A.slip_plane)*A.slip_plane
    r_in_plane = r - r_off_plane
    d_climb = np.linalg.norm(r_off_plane)
    d_glide = np.linalg.norm(r_in_plane)
    if d_climb**2./Lclimb**2.+ d_glide**2./Lglife**2. < 1.:
        return True
    else:
        return False

def interact(A,B)
    #We shpould consider adding in a dissipation output to aid in energy balance
    #(i.e. calculate the distance they would travel to annihilate)
    products = []
    if abs(np.dot(A.burgers,B.burgers)+1)<.05:
        return products
    else:
        products.append(A)
        products.append(B)
        return products
    

        
        
        
        
        

	

import numpy as np

class Dislocation(object):
    """A dislocation
    
    Attributes:
        X: location of dislocation <numpy vector>
        burgers: burgers vector <numpy vector>
        line_vec: line vector <numpy vector>
        slip_plane: slip plane <numpy vector>
    """
    
    def __init__(self, location, bv, sp, mu,nu):
        """Returns an initialized dislocation object.  Only dislocations 
           with line vectors in the z direction are possible."""
        self.X=np.copy(location)
        self.burgers=np.copy(bv)
        self.line_vec=np.array((0,0,1))
        self.slip_plane=np.copy(sp)/np.linalg.norm(sp)

        self.has_screw_component=abs(bv[2])>np.linalg.norm(bv)/100
        self.has_edge_component=abs(np.linalg.norm(bv[0:1]))>np.linalg.norm(bv)/100
        
        self.bedge = np.array((self.burgers[0],self.burgers[1],0.))
        self.bedge_norm = np.linalg.norm(self.bedge)
        self.edge_const = -mu*self.bedge_norm/(2*np.pi*(1-nu))
        
        self.bscrew = np.array((0.,0.,self.burgers[2]))
        self.bscrew_norm = self.burgers[2] 
        self.screw_const=-mu*self.bscrew_norm/(2*np.pi)
        
        self.update_g_mats()
        

    def stress_at_point(self,Y,mu,nu):
        """calculates the stress at point Y caused by a dislocation (self) 
        in a material with elastic parameters nu and mu"""
        if self.has_edge_component and self.has_screw_component:
            sigma = self.stress_screw(Y,mu,nu) + self.stress_edge(Y,mu,nu)
        elif self.has_edge_component:
            sigma = self.stress_edge(Y,mu,nu)
        elif self.has_screw_component:
            sigma = self.stress_screw(Y,mu,nu)
        else:
            print "***WARNING**** dislocation without screw or edge component!"
        
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
        update_g_mats(self)

    def rotate_and_stretch(self,beta):
        """Proposed duplicate of rotate, but allowing the second order effect of 
        changing the length of Burgers vector"""
        Q = beta
        Q[0,2] = 0.
        Q[1,2] = 0.
        Q[2,2] = 1.
        self.burgers = Q*self.burgers
        self.slip_plane = Q*self.slip_plane
        self.slip_plane = self.slip_plane/np.linalg.norm(self.slip_plane)
        update_g_mats(self)

    def get_velocity(self,sigma,drag):
        """ finds the velocity of dislocation (self) based on the stress (sigma)
        and drag"""
        F = np.cross(np.dot(self.burgers,sigma),self.line_vec)
        V = (F - np.dot(F,self.slip_plane)*self.slip_plane)/drag
        return V
        
    def update_g_mats(self):
        """precaculates the g matrices used in stress calculations"""
        #update the g matrix for screw component
        if np.abs(np.abs(self.line_vec[0])-1)>.0001:
        	yvec=np.cross(self.line_vec,np.array((1.,0.,0.)))
        else:
        	yvec=np.cross(self.line_vec,np.array((0.,1.,0.)))
        
        yvec = yvec/np.linalg.norm(yvec)
        xvec = np.cross(yvec,self.line_vec)
        xvec = xvec/np.linalg.norm(xvec)
        
        self.g_screw = np.vstack((xvec,yvec,self.line_vec))
        self.g_screw = np.transpose(self.g_screw)
        
        #update g matrix for edge component
        self.g_edge = np.vstack((self.bedge/self.bedge_norm,np.cross(self.line_vec,self.bedge/self.bedge_norm),self.line_vec))
        self.g_edge = np.transpose(self.g_edge)
        
    def stress_screw(self, Y, mu, nu):
        """finds stress caused by screw component of dislocation (self) at point
        Y in a material with elastic properties mu, nu"""
        dx = Y-self.X;
        rx = self.g_screw[0,0]*dx[0] + self.g_screw[1,0]*dx[1]
        ry = self.g_screw[0,1]*dx[0] + self.g_screw[1,1]*dx[1] #np.dot(np.transpose(self.g_screw),Y-self.X)
        
        rn2 = rx*rx + ry*ry #np.linalg.norm(r)
        
        if rn2>1e-16:
            sig=np.empty((3,3))
            sig[0,2] = ry*self.screw_const/(rn2)
            sig[2,0] = sig[0,2]
            sig[1,2] = rx*self.screw_const/(rn2)
            sig[2,1] = sig[1,2]
            sig[0,0]=0.
            sig[0,1]=0.
            sig[1,0]=0.
            sig[1,1]=0.
            sig[2,2]=0.
        else:        
            sig = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        
        return np.dot(np.dot(self.g_screw,sig),np.transpose(self.g_screw))  #should be able to do this

    def stress_edge(self, Y, mu, nu):
        """finds stress caused by edge component of dislocation (self) at point
        Y in a material with elastic properties mu, nu"""
        
        dx = Y-self.X;
        rx= self.g_edge[0,0]*dx[0] + self.g_edge[1,0]*dx[1]
        ry= self.g_edge[0,1]*dx[0] + self.g_edge[1,1]*dx[1] #np.dot(np.transpose(self.g_edge),Y-self.X)
        
        rx2=rx*rx
        ry2=ry*ry
        rn2 = rx2 + ry2 #np.linalg.norm(r)
        
        if rn2>1e-16:
            rn4=rn2*rn2
            
            sxx = self.edge_const*(ry/rn4)*(ry2 + 3*rx2)
            syy = self.edge_const*(ry/rn4)*(ry2 - rx2)
            sxy = -self.edge_const*(rx/rn4)*(rx2 - ry2)
            szz = self.edge_const*2*nu*(ry2/(rn2))
            
            tmp1=self.g_edge[0,0]*sxx+self.g_edge[0,1]*sxy
            tmp2=self.g_edge[0,0]*sxy+self.g_edge[0,1]*syy
            
            sig=np.empty((3,3))
            sig[0,0]=self.g_edge[0,0]*tmp1+self.g_edge[0,1]*tmp2
            sig[0,1]=self.g_edge[1,0]*tmp1+self.g_edge[1,1]*tmp2
            sig[1,0]=sig[0,1]
            sig[1,1]=self.g_edge[1,0]*self.g_edge[1,0]*sxx+self.g_edge[1,1]*(2*self.g_edge[1,0]*sxy+self.g_edge[1,1]*syy)
            sig[2,2]=szz
            sig[0,2]=0.
            sig[1,2]=0.
            sig[2,0]=0.
            sig[2,1]=0.
            #np.dot(np.dot(self.g_edge,sig),np.transpose(self.g_edge))
        else:
            sig = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
            
        return sig 

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
    Rmag = np.sqrt(R[0]*R[0] + R[1]*R[1]) #np.linalg.norm(R)
    Rhat = R/Rmag
    Ramag = np.sqrt(Ra[0]*Ra[0] + Ra[1]*Ra[1]) #np.linalg.norm(Ra)
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
    logRRa=np.log(Rmag/Ramag)
    dE = ba[2]*bb[2]*logRRa
    dE += (1/(1-nu))*( (ba[1]*bb[1] + ba[0]*bb[0])*logRRa
                      +(ba[1]*Rhat[0] - ba[0]*Rhat[1])*(bb[1]*Rhat[0] - bb[0]*Rhat[1])
                      -(ba[1]*Rahat[0] - ba[0]*Rahat[1])*(bb[1]*Rahat[0] - bb[0]*Rahat[1]))
    dE *= -mu/(2*np.pi)
    
    return dE

def check_for_interaction(A,B, L_glide, L_climb):
    """Checks to see if B is in range of A.  Extinction threshold is an oblate 
    ellipsoid"""
    r = B.X - A.X
    r_off_plane = (r[0]*A.slip_plane[0] + r[1]*A.slip_plane[1])*A.slip_plane
    r_in_plane = r - r_off_plane
    d_climb2 = (r_off_plane[0]*r_off_plane[0] + r_off_plane[1]*r_off_plane[1]) #np.linalg.norm(r_off_plane)
    d_glide2 = (r_in_plane[0]*r_in_plane[0] + r_in_plane[1]*r_in_plane[1]) #np.linalg.norm(r_in_plane)
    if d_climb2/(L_climb*L_climb)+ d_glide2/(L_glide*L_glide) < 1.:
        return True
    else:
        return False

def interact(A,B):
    #We shpould consider adding in a dissipation output to aid in energy balance
    #(i.e. calculate the distance they would travel to annihilate)
    if abs(np.dot(A.burgers/abs(np.linalg.norm(A.burgers)),B.burgers/abs(np.linalg.norm(B.burgers)))+1)<.05:
        print "detected annihilation"
        return True, None
    else:
        return False, None
    # else: just one product
    #    return True, do_reaction(A,B)
    

        
        
        
        
        

	

import numpy as np

class DislocationSrc(object):

    def __init__(self, location, bv, sp, L_nuc, mu, nu):
        """Returns an initialized dislocation source object.  Only sources 
           with line vectors in the z direction are possible.  Lnuc
            should be much larger than Lannihilation, perhaps by a factor of 2.
            Also, I only trust this with edges so far"""
        self.X=np.copy(location)
        self.burgers=np.copy(bv)
        self.line_vec=np.array((0,0,1))
        self.slip_plane=np.copy(sp)/np.linalg.norm(sp)
        self.tau_nuc = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*L_nuc))
        

    def load(self,sigma,mu,nu, L_glide):
        products = []
        tau = np.dot(self.burgers,np.dot(sigma,self.slip_plane))
        if abs(tau)>self.tau_nuc:
            L = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*tau))
            if L<L_glide:
                print "These dislocations will immediately annihilate!!!"
                #We could add two dipoles?  Output this error and reset the time step?
                #This is unlikely to occur anyway if Lnuc is twice L_glide.
                #We might have no issue if we step again before checking for annihilation
            Xstep = L*self.burgers/np.linalg.norm(self.burgers)/2.
            products.append(Dislocation.Dislocation( self.X + Xstep, self.burgers, self.slip_plane ))
            products.append(Dislocation.Dislocation( self.X - Xstep, -self.burgers, self.slip_plane ))            
        return products

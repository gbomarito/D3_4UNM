import numpy as np

class DislocationSrc(object):

    def __init__(self, location, bv, sp, Lnuc, mu, nu):
        """Returns an initialized dislocation source object.  Only sources 
           with line vectors in the z direction are possible.  Lnuc
            should be much larger than Lannihilation, perhaps by a factor of 2.
            Also, I only trust this with edges so far"""
        self.X=np.copy(location)
        self.burgers=np.copy(bv)
        self.line_vec=np.array((0,0,1))
        self.slip_plane=np.copy(sp)/np.linalg.norm(sp)
        self.taunuc = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*Lnuc))
        

    def load(self,sigma,mu,nu):
        products = []
        tau = np.dot(self.burgers,np.dot(sigma,self.slip_plane))
        if abs(tau)>self.taunuc:
            L = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*tau))
            Xstep = L*self.burgers/np.linalg.norm(self.burgers)/2.
            products.append(Dislocation.Dislocation( self.X + Xstep, self.burgers, self.slip_plane ))
            products.append(Dislocation.Dislocation( self.X - Xstep, -self.burgers, self.slip_plane ))            
        return products

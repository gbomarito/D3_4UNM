import numpy as np
import Dislocation

class DislocationSrc(object):

    def __init__(self, location, bv, sp, L_nuc, t_nuc, mu, nu):
        """Returns an initialized dislocation source object.  Only sources 
           with line vectors in the z direction are possible.  L_nuc
            should be much larger than L_annihilation, perhaps by a factor of 2.
            Also, I only trust this with edges so far"""
        self.X=np.copy(location)
        self.burgers=np.copy(bv)
        self.line_vec=np.array((0,0,1))
        self.slip_plane=np.copy(sp)/np.linalg.norm(sp)
        self.L_nuc = L_nuc
        self.tau_nuc = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*L_nuc)
        self.t_nuc = t_nuc
        self.load_time = 0.
        

    def load(self,sigma,mu,nu, L_glide, dt):
        products = []
        tau = np.dot(self.burgers/np.linalg.norm(self.burgers),np.dot(sigma,self.slip_plane))
        if abs(tau)>self.tau_nuc:
            self.load_time += dt*np.sign(tau)
        else:
            self.load_time = 0
        if abs(self.load_time)>=self.t_nuc:
            L = mu*np.linalg.norm(self.burgers)/(2*np.pi*(1-nu)*tau)
            if abs(L)<L_glide:
                print "These dislocations will immediately annihilate!!!"
                #We could add two dipoles?  Output this error and reset the time step?
                #This is unlikely to occur anyway if Lnuc is twice L_glide.
                #We might have no issue if we step again before checking for annihilation
            """while abs(L)<1.1*L_glide: #A possible fix, simulates multiple nucleation events
                Xstep = Xstep = 1.1*L_glide*self.burgers/np.linalg.norm(self.burgers)/2.
                products.append(Dislocation.Dislocation( self.X + Xstep, self.burgers, self.slip_plane ))
                products.append(Dislocation.Dislocation( self.X - Xstep, -self.burgers, self.slip_plane ))
                tau = tau - mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*1.1*L_glide))
                L = mu*np.linalg.norm(bv)/(2*np.pi*(1-nu)*tau))"""
            Xstep = np.sign(self.load_time)*self.L_nuc*self.burgers/np.linalg.norm(self.burgers)/2.
            products.append(Dislocation.Dislocation( self.X + Xstep, self.burgers, self.slip_plane ))
            products.append(Dislocation.Dislocation( self.X - Xstep, -self.burgers, self.slip_plane ))
            self.load_time = 0           
        return products

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
        """Proposed duplicate of rotate, but allowing the second order effect of 
        changing the length of Burgers vector"""
        Q = beta
        Q[0,2] = 0.
        Q[1,2] = 0.
        Q[2,2] = 1.
        self.burgers = Q*self.burgers
        self.slip_plane = Q*self.slip_plane
        self.slip_plane = self.slip_plane/np.linalg.norm(self.slip_plane)

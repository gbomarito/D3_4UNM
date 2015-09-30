import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import random

import Dislocation
import DislocationSrc


class DislocationMngr(object):
    """A manager that takes care of disloaction storage, querying, and steping
    
    Attributes:
        dislocations: list of dislocation objects <list of Dislocaton>
        nu: poisson coefficient <float>
        mu: lame' parameter <float>
        drag: drag coefficient <float>
    """
    
    def __init__(self, n, m,d, dnum=2, b=1.0, sim_size=1.0, ad_glide=20., ad_climb=5.):
        """Returns an initialized dislocation object"""
        self.nu=n
        self.mu=m
        self.drag=d
        self.sim_size=sim_size
        self.annihilation_dist_glide=ad_glide
        self.annihilation_dist_climb=ad_climb
        
        self.dislocations=[]
        '''
        self.dislocations.append(Dislocation.Dislocation( np.array((-sim_size/2.,0.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((sim_size/2.,sim_size/4.,0.)), np.array((-b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((sim_size*3./4.,0.,0.)), np.array((-b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((sim_size*-3./4.,0.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((sim_size*-5./4.,sim_size/2.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((0.,-sim_size/2.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)) ))
        self.dislocations.append(Dislocation.Dislocation( np.array((0.,sim_size/2.,0.)), np.array((-b,0.,0.)), np.array((0.,1.,0.)) ))'''

        self.sources=[]
        L_nuc = 40.*b
        self.sources.append(DislocationSrc.DislocationSrc( np.array((-sim_size/2.0,-sim_size/2.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)), L_nuc, self.mu, self.nu ))
        self.sources.append(DislocationSrc.DislocationSrc( np.array((sim_size/2.0,sim_size/2.,0.)), np.array((b,0.,0.)), np.array((0.,1.,0.)), L_nuc, self.mu, self.nu ))
        
        
        for i in range(dnum):
            if random.random()>0.5:
                burgers=np.array((b,0.,0.))
            else:
                burgers=np.array((-b,0.,0.))
            loc=np.array((-sim_size+random.random()*2*sim_size,-sim_size+random.random()*2*sim_size,0.0))
            self.dislocations.append(Dislocation.Dislocation( loc, burgers, np.array((0.,1.,0.)) ))
        
        
        
        self.velocities_last=[None]*len(self.dislocations)
        self.dt_last=None
        self.E_bal_skip=False
        self.L_glide = 20.*b
        self.L_climb = 5.*b
        
        

    def dd_step(self,dt,sig_ff=None,FE_results=None, E_TOL=1e-3):
        """Increments dislocations based on a time step (dt)
        optionally finite element results or a far field stress can be added"""
        
        dnum=len(self.dislocations)

        #print sig_ff(np.array([[0.,0.,0.]]))
        
        #calculate stresses at each dislocation point
        sigma=[None]*dnum
        for i in range(dnum):
            d_i=self.dislocations[i]
            sigma[i]=np.zeros((3,3))
        
            #stress caused by all dislocations
            for d_j in self.dislocations:
                if d_i is not d_j:
                    sigma[i] += d_j.stress_at_point(d_i.X,self.mu,self.nu)
                
            #far field stress
            if sig_ff is not None:
                sigma[i] += sig_ff(d_i.X)
                
            #TODO add FE stress field  
            
        #get velocities_next of each dislocation
        velocities_next=[np.zeros(3)]*dnum
        velocities=[np.zeros(3)]*dnum
        for i in range(dnum):
            d_i=self.dislocations[i]
            
        #adaptive timestep loop
        E_balanced=False
        bal_count=0
        while not E_balanced:
            #get current velocities of each dislocation
            for i in range(dnum):
                d_i=self.dislocations[i]
                
                velocities_next[i]=d_i.get_velocity(sigma[i],self.drag)  
                if self.velocities_last[i] is None:
                    velocities[i]=np.copy(velocities_next[i])
                else:
                    velocities[i]=velocities_next[i]+0.5*dt/self.dt_last*(velocities_next[i]-self.velocities_last[i]) #verlot integration
                    #velocities[i]=0.5*(velocities_next[i]+self.velocities_last[i]) #central diff 
            
            #ENERGY CALCS
            #caculate helper for energy dissipation by drag
            v_dot_v = 0.0
            for i in range(dnum):
                v_dot_v += np.dot(velocities[i],velocities[i])

            E_drag=v_dot_v*self.drag*dt
            E_dislocation=0.0
            for i in range(dnum):
                dx_i=velocities[i]*dt
                for j in range(i+1,dnum):
                    dx_j=velocities[j]*dt
                    E_dislocation += Dislocation.interaction_energy(self.dislocations[i],self.dislocations[j],dx_i,dx_j, self.mu,self.nu)
            
            #print "dt: ",dt," drag: ",E_drag," disloc: ",E_dislocation, "ratio: ", E_drag/E_dislocation
            if self.E_bal_skip:
                E_balanced=True
                self.E_bal_skip=False
            else:
                E_balanced = np.abs(E_drag+E_dislocation)/E_drag<E_TOL
            
            if not E_balanced:
                dt=dt/2.0
                
            bal_count+=1;
            if bal_count>100:
                raw_input()

        #update history variables
        self.dt_last=dt
        for i in range(dnum):
            self.velocities_last[i]=np.copy(velocities_next[i])        
    
        #move dislocations
        for i in range(dnum):
            dx_i=velocities[i]*dt
            self.dislocations[i].move_set(dx_i)
            
        #check for reactions
        removal_list=[]
        addition_list=[]
        for i in range(dnum):
            d_i=self.dislocations[i]
            for j in range(dnum):
                if i not in removal_list and j not in removal_list and i is not j:
                    d_j=self.dislocations[j]
                    if Dislocation.check_for_interaction(d_i,d_j,self.annihilation_dist_glide,self.annihilation_dist_climb):
                        doReaction,products=Dislocation.interact(d_i,d_j)
                        if doReaction:
                            removal_list.append(i)
                            removal_list.append(j)
                            if products is not None:
                                for p in products:
                                    addition_list.append(p)
        # remove the reactants
        removal_list.sort()
        for i in reversed(removal_list):
            print "deleteing: ",i
            del self.dislocations[i]
            del self.velocities_last[i]
            self.E_bal_skip=True
            #self.velocities_last=[None]*len(self.dislocations) #reset velocities so energy balance works in next step

        #Nucleate dislocations at sources
        """for i in self.sources:
            thissig = self.stress_at_point(i.X)
            #far field stress
            if sig_ff is not None:
                thissig += sig_ff(i.X)
            res = i.load(thissig,self.mu,self.nu, self.L_glide)
            if len(res) > 0:
                self.E_bal_skip=True
                for j in res:
                    self.dislocations.append(j)
                    self.velocities_last.append(None)#0.0?"""
            
        
        # add reaction products
        #TODO add code to add products of reactions in addition_list
        
        return dt, abs(E_drag+E_dislocation)/E_drag


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
        
        print "d: ",np.linalg.norm(self.dislocations[0].X-self.dislocations[1].X)/np.linalg.norm(self.dislocations[0].burgers)
        #count=0
        #for d_i in self.dislocations:
        #    print count,": ",d_i.X
        #    count+=1


    def plot(self,filename):
        """plot of all dislocations is output to a file"""
        
        x=[]
        y=[]
        theta=[]
        
        h=plt.figure() 
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        for s in self.sources:
            plt.text(s.X[0],s.X[1],'o')
        h.axes[0].set_xticks([])
        h.axes[0].set_yticks([]) 
        plt.axis((-1.5*self.sim_size,1.5*self.sim_size,-1.5*self.sim_size,1.5*self.sim_size))
        pylab.savefig(filename+".png")
        plt.close()


    def plot_w_stress(self,filename,sig_ff=None,FE_results=None):
        """plot of all dislocations with a stress contour"""

        bound = 1e7
        
        #find stress field
        delta = 0.05*self.sim_size
        sx = np.arange(-1.5*self.sim_size, 1.5*self.sim_size, delta)
        sy = np.arange(-1.5*self.sim_size, 1.5*self.sim_size, delta)
        sX, sY = np.meshgrid(sx, sy)
        pts=np.hstack( (sX.reshape(-1,1), sX.reshape(-1,1)) )
        sig_xx=np.zeros(sX.shape)
        sig_yy=np.zeros(sX.shape)
        sig_zz=np.zeros(sX.shape)
        sig_xy=np.zeros(sX.shape)
        sig_xz=np.zeros(sX.shape)
        sig_yz=np.zeros(sX.shape)
        for j in range(len(sx)):
            for i in range(len(sy)):
                loc=np.array((sX[i,j],sY[i,j],0.))
                
                sig=self.stress_at_point(loc) #disloc stress
                
                if sig_ff is not None:          #far field
                    sig += sig_ff(loc)
                
                sig_xx[i,j]=float(sig[0,0])*(abs(sig[0,0])<bound) + np.sign(sig[0,0])*bound*(abs(sig[0,0])>bound)
                sig_yy[i,j]=float(sig[1,1])*(abs(sig[1,1])<bound) + np.sign(sig[1,1])*bound*(abs(sig[1,1])>bound)
                sig_zz[i,j]=float(sig[2,2])*(abs(sig[2,2])<bound) + np.sign(sig[2,2])*bound*(abs(sig[2,2])>bound)
                sig_xy[i,j]=float(sig[0,1])*(abs(sig[0,1])<bound) + np.sign(sig[0,1])*bound*(abs(sig[0,1])>bound)
                sig_xz[i,j]=float(sig[0,2])*(abs(sig[0,2])<bound) + np.sign(sig[0,2])*bound*(abs(sig[0,2])>bound)
                sig_yz[i,j]=float(sig[1,2])*(abs(sig[1,2])<bound) + np.sign(sig[1,2])*bound*(abs(sig[1,2])>bound)
        
        
        h=plt.figure()        
        plt.contourf(sX,sY,sig_xx, levels=np.linspace(-bound,bound,100) )
        h.axes[0].set_xticks([])
        h.axes[0].set_yticks([])
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xx.png")
        plt.close()
        
        h=plt.figure()        
        plt.contourf(sX,sY,sig_yy, levels=np.linspace(-bound,bound,100) )
        h.axes[0].set_xticks([])
        h.axes[0].set_yticks([])
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_yy.png")
        plt.close()
        
        h=plt.figure()        
        plt.contourf(sX,sY,sig_zz, levels=np.linspace(-bound,bound,100) )
        h.axes[0].set_xticks([])
        h.axes[0].set_yticks([])
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_zz.png")
        plt.close()
        
        h=plt.figure()        
        plt.contourf(sX,sY,sig_xy, levels=np.linspace(-bound,bound,100) )
        h.axes[0].set_xticks([])
        h.axes[0].set_yticks([])
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xy.png")
        plt.close()
        
        '''
        plt.figure()        
        plt.contourf(sX,sY,sig_xz,levels=np.linspace(np.min(sig_xz),np.max(sig_xz),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_xz.png", bbox_inches='tight')
        
        plt.figure()        
        plt.contourf(sX,sY,sig_yz,levels=np.linspace(np.min(sig_yz),np.max(sig_yz),100))
        for d in self.dislocations:
            theta=180.0/np.pi * math.atan2(d.burgers[1],d.burgers[0])
            plt.text(d.X[0],d.X[1],r'$\perp$',ha='center',rotation=theta,va='center')
        plt.colorbar()
        pylab.savefig(filename+"_yz.png", bbox_inches='tight')
        '''
        
        

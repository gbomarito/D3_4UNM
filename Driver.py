import numpy as np
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    dt=1e-3
    num_time_steps=200
    
    nu=0.3
    mu=1.
    drag=0.5
    
    #initialize dislocation manager
    dm=DislocationMngr(nu,mu,drag)
    
    #main time loop
    for i in range(num_time_steps):
        dm.dd_step(dt)
        print "\n---:: Time Step ",i," ::---"
        dm.dump()
        if i%10 is 0:
            dm.plot_w_stress("sig{0:03d}".format(i))  
 


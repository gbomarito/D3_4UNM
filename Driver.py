import numpy as np
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    dt=1e-5
    num_time_steps=10
    
    nu=0.3
    mu=1.
    drag=1.
    
    #initialize dislocation manager
    dm=DislocationMngr(nu,mu,drag)
    
    #main time loop
    for i in range(num_time_steps):
        dm.dd_step(dt)
        print "\n---:: Time Step ",i," ::---"
        dm.dump()
        #dm.plot('output.png')


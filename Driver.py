import numpy as np
import DislocationMngr
import Dislocation


################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    dt=1e-5
    num_time_steps=10
    
    nu=
    mu=
    drag=
    
    #initialize dislocation manager
    dm=DisloactionMngr(nu,mu,drag)
    
    #main time loop
    for i in range(num_time_steps):
        dm.dd_step(dt)
        print "\n---:: Time Step ",i," ::---"
        dm.dump

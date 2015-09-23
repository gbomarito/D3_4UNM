import numpy as np
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    dt=5e-5
    num_time_steps=1000
    
    nu=0.3
    mu=27.0
    drag=0.5
    
    #initialize dislocation manager
    burgers=4.05e-10/(2.**.5)
    dm=DislocationMngr(nu,mu,drag,dnum=20, b=burgers, area=burgers*100)
    
    #main time loop
    plot_res=200
    for i in range(num_time_steps):
        dm.dd_step(dt)
        print "\n---:: Time Step ",i," ::---"
        #dm.dump()
        if i%plot_res is 0:
            dm.plot_w_stress("sig{0:04d}".format(i/plot_res))  
 


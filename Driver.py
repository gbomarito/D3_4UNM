import numpy as np
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    suggested_dt=1e-3
    sim_time=1e-1
    
    nu=0.3
    mu=27.0
    drag=1.0
    
    #initialize dislocation manager
    burgers=4.05e-10/(2.**.5)
    dm=DislocationMngr(nu,mu,drag,dnum=15, b=burgers, sim_size=burgers*100)
    
    #main time loop
    plot_res=sim_time/30.0
    t_count=1
    next_plot_time=0.0
    while time<sim_time:
        dt,E_bal = dm.dd_step(suggested_dt)
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "\ttime step: ",dt
        print "\tcurrent time: ",time
        print "\tenergy balance: ",E_bal
        #dm.dump()
        if time> next_plot_time:
            dm.plot_w_stress("sig{0:04d}".format(i/plot_res))  
            next_plot_time+=plot_res
        t_count+=1
 


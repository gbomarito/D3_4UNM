import numpy as np
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    suggested_dt=1e-4
    sim_time=5e-1
    
    nu=0.3
    mu=270000000
    drag=1.0e4
    
    #initialize dislocation manager
    burgers=4.05e-8/(2.**.5)    #
    dm=DislocationMngr(nu,mu,drag,dnum=15, b=burgers, sim_size=burgers*100)
    
    #main time loop
    plot_res=sim_time/100.0
    t_count=1
    next_plot_time=0.0
    plot_num=0
    time=0.
    while time<sim_time:
        dt,E_bal = dm.dd_step(suggested_dt)
        time+=dt
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "\ttime step: ",dt
        print "\tcurrent time: ",time
        print "\tenergy balance: ",E_bal
        dm.dump()
        if time> next_plot_time:
            dm.plot_w_stress("sig{0:04d}".format(plot_num))  
            next_plot_time+=plot_res
            plot_num+=1
        t_count+=1
 


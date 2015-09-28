import numpy as np
import os
from DislocationMngr import DislocationMngr



################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    #os.system("rm sig*.png")
    #os.system("rm D*.png")
    
    max_dt=1e-2
    sim_time=1e-0
    
    nu=0.3
    mu=270000000
    drag=1.0e4
    
    energy_tol=1e-3
    
    #initialize dislocation manager
    burgers=4.05e-8/(2.**.5)    #
    dm=DislocationMngr(nu,mu,drag,dnum=50, b=burgers, sim_size=burgers*100)
    #dm.plot_w_stress("sig{0:04d}".format(0))
    dm.plot("D{0:04d}".format(0))  
    
    #main time loop
    plot_res=sim_time/50.0
    t_count=1
    next_plot_time=plot_res
    plot_num=1
    time=0.
    suggested_dt=max_dt
    while time<sim_time:
    
        #do timestep
        dt,E_bal = dm.dd_step(suggested_dt,E_TOL=energy_tol)
        time+=dt
        
        #adjust timestep based on last timestep
        if dt<suggested_dt:
            suggested_dt=dt
        elif E_bal < energy_tol/10.0:
            suggested_dt=suggested_dt*2 #min(max_dt,suggested_dt*2.0)
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "\ttime step: ",dt
        print "\tcurrent time: ",time
        print "\tenergy balance: ",E_bal
        #dm.dump()
        if time> next_plot_time:
            #dm.plot_w_stress("sig{0:04d}".format(plot_num)) 
            dm.plot("D{0:04d}".format(plot_num))  
            next_plot_time+=plot_res
            plot_num+=1
        t_count+=1
    
    

    #print "CONVERTING GIF"    
    #os.system("convert -delay 10 -loop 0 D*.png animation.gif")
 


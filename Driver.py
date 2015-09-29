import numpy as np
import os
from DislocationMngr import DislocationMngr

def simple_shear():
        sig = sigma=np.zeros((3,3))
        sig[0,1] = 2700000.
        sig[1,0] = sig[0,1]
        return sig

################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    os.system("rm sig*.png")
    os.system("rm D*.png")
    
    
    
    #----:: INPUT ::----
    # simulation params
    starting_dt=1e-2
    sim_time=2e-0
    
    # output params
    plot_res=sim_time/200.
    
    # material params
    nu=0.3
    mu=270000000.
    drag=1.0e4
    burgers=4.05e-8/(2.**.5)
    
    # integration/discretization params
    energy_tol=1e-3
    
    
    
    #----:: ININTIALIZATION ::----
    # initialize dislocation manager
    dm=DislocationMngr(nu,mu,drag,dnum=50, b=burgers, sim_size=burgers*100, ad_climb=5*burgers, ad_glide=20*burgers)
    #dm.plot_w_stress("sig{0:04d}".format(0))
    dm.plot("D{0:04d}".format(0))  
    
    
    
    #----:: MAIN TIME LOOP ::----
    t_count=1
    next_plot_time=plot_res
    plot_num=1
    time=0.
    suggested_dt=starting_dt
    while time<sim_time:
        """lambda X: time*simple_shear()/sim_time"""
        #do timestep
        dt,E_bal = dm.dd_step(suggested_dt,sig_ff= None,E_TOL=energy_tol)
        time+=dt
        
        #adjust timestep based on last timestep
        if dt<suggested_dt:
            suggested_dt=dt
        elif E_bal < energy_tol/10.0:
            suggested_dt=suggested_dt*2
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "\ttime step: ",dt,
        print "\tcurrent time: ",time,
        print "\tenergy balance: ",E_bal
        #dm.dump()
        if time> next_plot_time:
            dm.plot_w_stress("sig{0:04d}".format(plot_num),sig_ff=None) 
            #dm.plot("D{0:04d}".format(plot_num))  
            next_plot_time+=plot_res
            plot_num+=1
        t_count+=1
    
    

    print "CONVERTING GIF"    
    os.system("convert -delay 10 -loop 0 sig*_xy*.png animation.gif")
    os.system("gifview -a animation.gif&")

 


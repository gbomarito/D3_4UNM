import numpy as np
import os
from DislocationMngr import DislocationMngr

################################################################################
def const_shear(loc):
    shear=2e6
    return np.array(((0,shear,0),(shear,0,0),(0,0,0)))

def variable_shear(loc, t):
    shear = 1.8e6#1.3e6 + .5e6*t
    return np.array(((0,shear,0),(shear,0,0),(0,0,0)))


################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    os.system("rm sig*.png")
    os.system("rm D*.png")
    
    
    
    #----:: INPUT ::----
    # simulation params
    starting_dt=1e-9
    sim_time=2.e-8
    
    # output params
    plot_res=sim_time/100.
    
    # material params
    nu=0.3
    mu=27.e10 #dyne/cm^2
    drag=1.e-3 #poise
    burgers=4.05e-8/(2.**.5) #cm
    
    # integration/discretization params
    energy_tol=1.e-3
    
    
    
    #----:: ININTIALIZATION ::----
    # initialize dislocation manager
    dm=DislocationMngr(nu,mu,drag,dnum=20, b=burgers, sim_size=burgers*100, ad_climb=1*burgers, ad_glide=5*burgers)
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
        dt,E_bal = dm.dd_step(suggested_dt,E_TOL=energy_tol, sig_ff=lambda X: variable_shear(X,time/sim_time))
        time+=dt
        
        #adjust timestep based on last timestep
        if dt<suggested_dt:
            suggested_dt=dt
        elif E_bal < energy_tol/10.0:
            suggested_dt=suggested_dt*2
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "time step: ",dt,
        print "\tcurrent time: ",time,
        print "\tenergy balance: ",E_bal
        #dm.dump()
        if time> next_plot_time:
            dm.plot_w_stress("sig{0:04d}".format(plot_num),sig_ff=lambda X: variable_shear(X,time/sim_time)) 
            #dm.plot("D{0:04d}".format(plot_num))  
            next_plot_time+=plot_res
            plot_num+=1
        t_count+=1
    
    

    print "CONVERTING GIF"    
    os.system("convert -delay 10 -loop 0 sig*_xy*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 sig*_xx*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 sig*_xz*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 sig*_yz*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 sig*_zz*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 sig*_yy*.png animation.gif")
    #os.system("convert -delay 10 -loop 0 D*.png animation.gif")    
    os.system("gifview -a animation.gif&")

 


import numpy as np
import os
import optparse
import h5py
from DislocationMngr import DislocationMngr

################################################################################
def variable_shear(loc, t, initial_shear, shear_rate):
    shear = initial_shear+shear_rate*t #1.8e6#1.3e6 + .5e6*t
    return np.array(((0,shear,0),(shear,0,0),(0,0,0)))


################################################################################
if __name__ == "__main__":
    """ Main DDD driver """
    
    
    # Command line arguments:
    parser = optparse.OptionParser()
    parser.add_option("-i", type="string", dest="infile", 
                     help="Path to h5 database file containing the input")
    parser.add_option("-o", type="string", dest="outfile",
                    help="output file name")
    parser.add_option("-a", action="store_true", dest="do_animation")
    parser.add_option("-v", action="store_true", dest="verbose")
    parser.add_option("-V", action="store_true", dest="veryverbose")

    # Default values if not specified:
    parser.set_defaults(infile="input.h5",
                        outfile="out.h5",
                        do_animation=False,
                        verbose=False,
                        veryverbose=False)
    opts, args = parser.parse_args()
    
    
    #make plot folde, and make sure its clean
    if os.path.isdir("plots"):
        for root, dirs, files in os.walk('plots', topdown=False):
            for fname in files:
                name, ext = os.path.splitext(fname)
                ext = ext.lower()
                if ext == ".png":
                    os.remove(root+"/"+fname)
    else:
        os.mkdir("plots")
    
    
    
    #----:: INPUT ::----
    # open input file
    f=h5py.File(opts.infile,"r")
    
    # simulation params
    sim_time=f["Simulation"]["sim_time"][()]
    starting_dt=sim_time/1000;
    energy_tol=f["Simulation"]["energy_tolerance"][()]
    
    # plotting params
    plot_res=sim_time/f["Simulation"]["Plotting"]["plot_num"][()]
    plot_stresses=f["Simulation"]["Plotting"]["plot_stresses"][()]==1
    
    # external stresses
    initial_shear=f["Simulation"]["External_Stresses"]["initial_shear"][()]
    shear_rate=f["Simulation"]["External_Stresses"]["shear_rate"][()]
    
    # close input file
    f.close()
    
    
    #----:: ININTIALIZATION ::----
    # initialize dislocation manager
    dm=DislocationMngr(opts.infile)
    dm.plot("plots/D{0:04d}".format(0))  
    
    
    
    #----:: MAIN TIME LOOP ::----
    t_count=1
    next_plot_time=plot_res
    plot_num=1
    time=0.
    suggested_dt=starting_dt
    while time<sim_time:
        
        #calculate far field stress function
        far_field_stress=lambda X: variable_shear(X,
                                                 time,
                                                 initial_shear,
                                                 shear_rate)
        #do timestep
        dt,E_bal = dm.dd_step(suggested_dt,
                              E_TOL=energy_tol,
                              sig_ff=far_field_stress,
                              verbose=opts.verbose)
        time+=dt
        
        #adjust timestep based on last timestep
        if dt<suggested_dt:
            suggested_dt=dt
        elif E_bal < energy_tol/10.0:
            suggested_dt=suggested_dt*2
        
        #output
        print "\n---:: Time Step ",t_count," ::---"
        print "time step: {0:12.2e}".format(dt),
        print "\tcurrent time: {0:12.5e}".format(time),
        print "\tenergy balance: {0:.6f}".format(E_bal)
        if opts.veryverbose:
            dm.dump()
        if time> next_plot_time:
            if plot_stresses:
                dm.plot_w_stress("plots/sig{0:04d}".format(plot_num),
                                 sig_ff=far_field_stress) 
            dm.plot("plots/D{0:04d}".format(plot_num))  
            next_plot_time+=plot_res
            plot_num+=1
        t_count+=1
    
    

    #do animation if requested
    if opts.do_animation:
        print "CONVERTING GIF"    
        if plot_stresses:
            os.system("convert -delay 10 -loop 0 plots/sig*_xy*.png plots/animation.gif")
            #os.system("convert -delay 10 -loop 0 plots/sig*_xx*.png plots/animation.gif")
            #os.system("convert -delay 10 -loop 0 plots/sig*_xz*.png plots/animation.gif")
            #os.system("convert -delay 10 -loop 0 plots/sig*_yz*.png plots/animation.gif")
            #os.system("convert -delay 10 -loop 0 plots/sig*_zz*.png plots/animation.gif")
            #os.system("convert -delay 10 -loop 0 plots/sig*_yy*.png plots/animation.gif")
        else:
            os.system("convert -delay 10 -loop 0 plots/D*.png plots/animation.gif")    
        os.system("gifview -a plots/animation.gif&")

 


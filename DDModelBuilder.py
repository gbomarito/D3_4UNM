import numpy as np
import h5py
import random
import optparse

################################################################################
def make_model_from_prompt(): 
    """ Create a DD model by getting input from user """
    
    filename = raw_input('Enter a file name: ')
    
    print "\n---Material Properties---"
    nu = float(input('Poisson ratio, nu: '))
    mu = float(input('Shear modulus, mu: '))
    B = float(input('Drag coefficient, B: '))
    
    print "\n---Initial Dislocation Structure---"
    b = float(input('Length of Burgers vector: '))
    sim_size = float(input('Simulation size: '))
    d_rand = int(input('Number of random dislocations: '))
    
    print "\n---Annihilation---"
    L_glide = float(input('Annihilation distance for glide, L_glide: '))
    L_climb = float(input('Annihilation distance for climb, L_climb: '))
    
    print "\n---Sources---"
    s_rand = int(input('Number of randomly located dislocation sources: '))
    L_nuc = float(input('Nucleation distance, L_nuc: '))
    t_nuc = float(input('Nucleation time, t_nuc: '))
    
    
    #random dislocations
    d_bv_list=[]
    d_sp_list=[]
    d_loc_list=[]
    for i in range(d_rand):
        if random.random()>0.5:
            d_bv_list.append(np.array((b,0.,0.)))
        else:
            d_bv_list.append(np.array((-b,0.,0.)))
        d_loc_list.append(np.array((-sim_size+random.random()*2*sim_size,-sim_size+random.random()*2*sim_size,0.0)))
        d_sp_list.append(np.array((0.,1.,0.)))
    d_bv_array=np.array(d_bv_list)
    d_sp_array=np.array(d_sp_list)
    d_loc_array=np.array(d_loc_list)
    
    #random dislocation sources
    s_bv_list=[]
    s_sp_list=[]
    s_loc_list=[]
    s_t_list=[]
    s_L_list=[]
    for i in range(s_rand):
        s_bv_list.append(np.array((b,0.,0.)))
        s_loc_list.append(np.array((-sim_size+random.random()*2*sim_size,-sim_size+random.random()*2*sim_size,0.0)))
        s_sp_list.append(np.array((0.,1.,0.)))
        s_t_list.append(np.array((t_nuc)))
        s_L_list.append(np.array((L_nuc)))
    s_bv_array=np.array(s_bv_list)
    s_sp_array=np.array(s_sp_list)
    s_loc_array=np.array(s_loc_list)
    s_t_array=np.array(s_t_list)
    s_L_array=np.array(s_L_list)
    
    
    print "\n---Simulation---"
    sim_time = float(input('Simulation time: '))
    energy_tol = float(input('Energy balance tolerance: '))
    plot_num = int(input('Number of plots: '))
    plot_stresses = int(input('Do you want plots of stress fields? (1=yes, 0=no): '))
    
    print "\n---External Stress---"
    initial_shear = float(input('Initial shear: '))
    shear_rate = float(input('Shear rate: '))
    
    
    #write file
    f=h5py.File(filename,"w")
    dset=f.create_dataset("Model/nu", data=nu)
    dset=f.create_dataset("Model/mu", data=mu)
    dset=f.create_dataset("Model/B", data=B)
    dset=f.create_dataset("Model/b", data=b)
    dset=f.create_dataset("Model/sim_size", data=sim_size)
    dset=f.create_dataset("Model/L_glide", data=L_glide)
    dset=f.create_dataset("Model/L_climb", data=L_climb)
    dset=f.create_dataset("Model/Dislocations/burgers_vectors", data=d_bv_array)
    dset=f.create_dataset("Model/Dislocations/slip_planes", data=d_sp_array)
    dset=f.create_dataset("Model/Dislocations/locations", data=d_loc_array)
    dset=f.create_dataset("Model/Sources/burgers_vectors", data=s_bv_array)
    dset=f.create_dataset("Model/Sources/slip_planes", data=s_sp_array)
    dset=f.create_dataset("Model/Sources/locations", data=s_loc_array)
    dset=f.create_dataset("Model/Sources/t_nuc", data=s_t_array)
    dset=f.create_dataset("Model/Sources/L_nuc", data=s_L_array)
    
    dset=f.create_dataset("Simulation/sim_time", data=sim_time)
    dset=f.create_dataset("Simulation/energy_tolerance", data=energy_tol)
    dset=f.create_dataset("Simulation/Plotting/plot_num", data=plot_num)
    dset=f.create_dataset("Simulation/Plotting/plot_stresses", data=plot_stresses)
    dset=f.create_dataset("Simulation/External_Stresses/initial_shear", data=initial_shear)
    dset=f.create_dataset("Simulation/External_Stresses/shear_rate", data=shear_rate)
    f.close()
    

################################################################################
def edit_model_from_prompt(): 
    """allows user to edit previously created  model"""
    
        
    filename = raw_input('Enter a file name: ')
    f=h5py.File(filename,"r+")
    still_editing=True
    
    while still_editing:
        print "\n\n---Material Properties---"
        print ' (1)  Poisson ratio, nu: ',f["Model"]["nu"][()]
        print ' (2)  Shear modulus, mu: ',f["Model"]["mu"][()]
        print ' (3)  Drag coefficient, B: ',f["Model"]["B"][()]
        
        print "\n---Dislocation Structure---"
        print ' (4)  Simulation size: ',f["Model"]["sim_size"][()]
        print ' (5)  Dislocations... '
        
        print "\n---Annihilation---"
        print ' (6)  Annihilation distance for glide, L_glide: ',f["Model"]["L_glide"][()]
        print ' (7)  Annihilation distance for climb, L_climb: ',f["Model"]["L_climb"][()]
        
        print "\n---Sources---"
        print ' (8)  Sources... '
        
        print "\n---Simulation---"
        print " (9)  Simulation time: ",f["Simulation"]["sim_time"][()]
        print " (10) Energy Tolerance: ",f["Simulation"]["energy_tolerance"][()]
        print " (11) Number of plots: ",f["Simulation"]["Plotting"]["plot_num"][()]
        print " (12) Stress plot on (1=yes, 0=no): ",f["Simulation"]["Plotting"]["plot_stresses"][()]
        
        print "\n---External Stress---"
        print " (13) Initial shear ",f["Simulation"]["External_Stresses"]["initial_shear"][()]
        print " (14) Shear rate: ",f["Simulation"]["External_Stresses"]["shear_rate"][()]
        
        option=int(input('\nWhich would you like to edit: '))
        
        
        if option==1:
            f["Model"]["nu"][()] = float(input('new Poisson ratio, nu: '))
        elif option==2:
            f["Model"]["mu"][()] = float(input('new Shear modulus, mu: '))
        elif option==3:
            f["Model"]["B"][()] = float(input('new Drag coefficient, B: '))
        elif option==4:
            f["Model"]["sim_size"][()] = float(input('new Simulation size: '))
        elif option==5:
            print "OPTION NOT AVAILABLE"
        elif option==6:
            f["Model"]["L_glide"][()] = float(input('new Annihilation distance for glide, L_glide: '))
        elif option==7:
            f["Model"]["L_climb"][()] = float(input('new Annihilation distance for climb, L_climb: '))
        elif option==8:
            print "OPTION NOT AVAILABLE"
        elif option==9:
            f["Simulation"]["sim_time"][()] = float(input('new Simulation time: '))
        elif option==10:
            f["Simulation"]["energy_tolerance"][()] = float(input('new Energy balance tolerance: '))
        elif option==11:
            f["Simulation"]["Plotting"]["plot_num"][()] = int(input('new Number of plots: '))
        elif option==12:
            f["Simulation"]["Plotting"]["plot_stresses"][()] = int(input('Do you want plots of stress fields? (1=yes, 0=no): '))
        elif option==13:
            f["Simulation"]["External_Stresses"]["initial_shear"][()] = float(input('new Initial shear: '))
        elif option==14:
            f["Simulation"]["External_Stresses"]["shear_rate"][()] = float(input('new Shear rate: '))
        else:
            print "exiting..."
            still_editing=False
            
    f.close()
    

################################################################################
if __name__ == "__main__":
    """ Model builder for dd code """
    
    
    # Command line arguments:
    parser = optparse.OptionParser()
    parser.add_option("-n", action="store_true", dest="new")
    parser.add_option("-e", action="store_true", dest="edit")

    # Default values if not specified:
    parser.set_defaults(new=False,
                        edit=False)
    opts, args = parser.parse_args()
    
    if opts.new:
        make_model_from_prompt()
    elif opts.edit:
        edit_model_from_prompt()

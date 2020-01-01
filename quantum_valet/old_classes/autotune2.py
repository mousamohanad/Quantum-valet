import os
import logging
import numpy as np
from ase.calculators.vasp import Vasp2
from ase.io import read,write
from ase.eos import EquationOfState

class ValetAutotune:
    
    
    def __init__(self,POSCAR_init=None,path_out=None,kpts=None,calculators=None):
        # Initial structure handling
        if POSCAR_init == None:
            print("Please provide an initial structure file.")
            exit()
        elif os.path.isfile(POSCAR_init):
            self.poscar_init = POSCAR_init
            self.system = read(self.poscar_init)
        else:
            logging.error("Cannot locate POSCAR file as {}".format(POSCAR_init))
            self.error_read_log()
            exit()
            
        # Setting up some paths
        if path_out == None:
            print("Please enter a valid out directory.")
            exit()
        else:
            self.path_out = os.path.join("valet_workspace",path_out)
            self.path_log = os.path.join(self.path_out,"log")
            self.poscar_vol = os.path.join(self.path_out,"POSCAR_vol")
        
        # Check if user supplied kpoints
        if kpts != None:
            self.kpts = kpts
        else:
            self.kpts = [4,4,4]
            
        # Check if user supplied calculators
        self.calc_dict = {'encut': Vasp2(xc="PBE",kpts=self.kpts,gamma=True,directory=self.path_out,setups="recommended"),
			  'volume':Vasp2(xc="PBE",kpts=self.kpts,gamma=True,directory=self.path_out,setups="recommended")}
#			  'encut':Vasp2(xc='PBE',kpts=self.kpts,gamma=True,directory=self.path_out),
#                          'volume':Vasp2(xc='PBE',prec='Acc',  # First relaxation for curve
#                                        algo='Fast', 
#                                        icharg=0,  # Calc charge from initial wavefunction
#                                        nelm=100,  # Max electronc SC steps
#                                        ibrion=2,  # Conj. grad. ionic relax
#                                        ediff=1E-4,#0.00005*len(self.system),
#                                        nsw=100,   # Max number ionic steps
#                                        isif=4,    # Position,shape change. Not volume
#                                        lreal='Auto',
#                                        ismear=-5, # Tetrahedon method; requires gamma-centered
#                                        sigma=0.05,
#                                        lwave=False, # No wavecar
#                                        lorbit=11,   # DOS car and lm-decomposed PROCAR,
#                                        kpts=self.kpts, gamma=True,directory=self.path_out)
#                         
#                         }                           
                          #'volume':Vasp2(xc='PBE',kpts=self.kpts,gamma=True,directory=self.path_out)}
                        
        if calculators != None:
            print("Attaching your calculators.")
            if type(calculators) != type({}):
                logging.error("Calculators should be provided in dictionary form. The relevant keys are 'encut', \
                              'volume', ....more to come.")
                self.error_read_log()
                exit()
            else:
                for key in calculators.keys():
                    self.calc_dict[key] = calculators[key]
                    self.calc_dict[key].set(directory=self.path_out) # Make sure output path is correct
        else:
            print("Loading default calculators.")
    
        
        
            
        # Initialize the rest of the attributes
        self.encut = None
        self.band_structure = None
        self.band_gap = None
        
        # Make output directory
        if not os.path.isdir(self.path_log):
            os.system("mkdir -p {}/".format(self.path_log))
        # Start log
        logging.basicConfig(format='%(asctime)s : %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
               filename=os.path.join(self.path_log,"valet.log"), filemode="w",level=logging.INFO)  
      
    def error_read_log(self):
        print("ERROR: See log for error.")
    
    def clean_up(self):
        os.system("cd {}; rm ase* C* D* E* I* K* O* PC* POSCAR POT* vaspr* WA* X*".format(self.path_out))
    
    def set_encut_autotune_params(self,tol=None,max_steps=None,start_encut=None,encut_step=None):
        stop = False
        ens = []
        cnt = 0
        return tol,stop,cnt,max_steps,start_encut,encut_step,ens
    
    def check_stop(self,x,tol):
        perc_diff = lambda x,y: np.abs(np.abs(x-y)/((x+y)/2))*100
        if len(x) < 2:
            return False
        elif perc_diff(x[-1],x[-2]) < tol:
            return True
        else: 
            return False
        
    def do_encut_autotune(self,tol=1E-4,max_steps=32,start_encut=200,encut_step=25,retune=False):
        logging.info("Commencing ENCUT autotune.")
        if (self.encut == None) or (retune):
            sys = read(self.poscar_init)
            calc = Vasp2(xc='PBE',kpts=(4,4,4),directory=self.path_out)
            tol,stop,cnt,max_steps,ENCUT,ENCUT_step,ens = self.set_encut_autotune_params(tol=tol,max_steps=max_steps,
                                                                                         start_encut=start_encut,
                                                                                         encut_step=encut_step)

            print("Autotuning ENCUT...")
            while not stop:
                if 'error' in locals():
                    del error
                if cnt < max_steps:
                    logging.debug("Count {}".format(cnt))
                    calc.set(encut=ENCUT)
                    sys.set_calculator(calc)
                    ens.append(sys.get_potential_energy())
                    stop = self.check_stop(ens,tol)
                    ENCUT += ENCUT_step
                    cnt += 1
                else:
                    logging.debug("Maximum autotune steps reached.")
                    error = "Maximum autotune steps reached."
                    tip = "Tip: Increase either start_encut, encut_step, or max_steps. Decrease tol."
                    stop = True

            try:
                print(error)
                print(tip)
                self.encut = None
            except NameError:
                print("\tENCUT auto-tuned to: {}".format(ENCUT))
                logging.info("Completed ENCUT autotune.")
                self.encut = ENCUT
            print("Done.")  
        else:
            print("ENCUT already autotuned. If you would like to retune, please pass 'retune' as True to 'self.do_encut_autotune'.")
            logging.info("Already autotuned ENCUT.")
            
        self.clean_up()
        
    def set_volume_autotune_params(self,start_scan=None,end_scan=None,num_scans=None,exclude_type=None,only_type=None):
        vols = []
        ens = []
        return start_scan,end_scan,num_scans,exclude_type,only_type,vols,ens

    def do_volume_autotune(self,encut=None,start_scan=0.925,end_scan=1.075,num_scans=20,exclude_type=None,only_type=None):
        if not os.path.isfile(self.poscar_vol):
            encut = encut if encut!=None else self.do_encut_autotune()

            logging.info("Commencing volume scan...")
            print("Determining optimal volume...")
            start_scan,end_scan,num_scans,exclude_type,only_type,vols,ens = self.set_volume_autotune_params(
                start_scan=start_scan,end_scan=end_scan,num_scans=num_scans,exclude_type=exclude_type,
                only_type=only_type)


            # Load POSCAR
            try:
                sys = read(self.poscar_init)
                calc = self.calc_dict['volume']
                calc.set(encut=encut)
                start_cell = sys.get_cell()
            except:
                logging.error("Initial POSCAR file could not be found at {}".format(poscar_init))
                print("Error: See error log for details.")
                return None
            logging.info("Loaded initial POSCAR file.")

            # Do volume scan
            logging.info("Performing volume scan.")
            for x in np.linspace(start_scan,end_scan,num_scans):
                sys.set_cell(x*start_cell,scale_atoms=True)
                sys.set_calculator(calc)
                ens.append(sys.get_potential_energy())
                vols.append(sys.get_volume())
            logging.info("Volume scan complete.")

            # Fit EoS
            logging.info("Fitting equations of state.")
            eos_types = "sjeos taylor murnaghan birch birchmurnaghan pouriertarantola p3 antonschmidt".split()

            if exclude_type != None:
                _ = [eos_types.remove(i) for i in exclude]
            elif only_type != None:
                eos_types = only_type

            
            
            eos_fits = {}
            for typ in eos_types:
                eos = EquationOfState(vols,ens,eos=typ)
                try:
                    v, e, B = eos.fit()
                    eos_fits[typ] = {'volume':v,'energy':e,'buld_modulus':B}
                except:
                    print("Unable to fit type {}.".format(typ))

            # Rescale initial cell to optimized volume
            logging.info("Optimal volume found. Rescaling original cell.")
            vol_avg = np.average([eos_fits[key]['volume'] for key in eos_fits.keys()])
            sys.set_cell(start_cell,scale_atoms=True)

            scale = sys.get_volume()/vol_avg
            sys.set_cell(scale*start_cell,scale_atoms=True)

            # Perform sc-step
            logging.info("Performing self-consistent step.")
            print("Second relaxation")
            #calc.set(isif=3)  # Everything can relax
            sys.set_calculator(calc)
            en = sys.get_potential_energy()
            
            # Set relaxed structure as active structure
            print("Setting relaxed structure to active structure.")
            self.system = sys

            # Save structure
            logging.info("Saving self-consistent calculator.")
            write(self.poscar_vol,sys)

            print("Done.")
            logging.info("Volume scan complete.")
            self.clean_up()
        else:
            sys = read(self.poscar_vol)
            print("Volume optimized system ready.")
        

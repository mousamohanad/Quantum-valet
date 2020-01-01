import os, logging, sys
from ase.io import read,write
import dftvalet
from ase.calculators.vasp import Vasp2
import numpy as np
from ase.eos import EquationOfState


class ValetSystem():

    def __init__(self,xyz_init=None,path_out=None,xyz_format=None):
        if self.set_paths(path_out):
            self.load_structure(xyz_init,xyz_format)
        # Initialize the rest of the attributes
        self.encut = None
        self.band_structure = None
        self.band_gap = None

    # Methods
    
    def set_paths(self,path_out):
        # Setting up some paths
        if path_out == None:
            logging.error("Please enter a valid out directory.")
            return False
        else:
            self.path_out = os.path.join("valet_workspace",path_out)
            self.path_log = os.path.join(self.path_out,"log")
            self.path_xyzs = os.path.join(self.path_out,"structures")
            self.xyzs = {'init':os.path.join(self.path_xyzs,"initial.xyz"),
                         'vol_scan':os.path.join(self.path_xyzs,"volume_scan.xyz"),
                         'opt':os.path.join(self.path_xyzs,"optimized.xyz")}
    
            # Make output directory
            if not os.path.isdir(self.path_log):
                os.system("mkdir -p {}/".format(self.path_log))
            if not os.path.isdir(self.path_xyzs):
                os.system("mkdir -p {}/".format(self.path_xyzs))
                                             
            # Start log
            #handler = logging.StreamHandler(sys.stdout)
            #handler.setLevel(logging.INFO)
            #root = logging.basicConfig(format='%(asctime)s : %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
            #       filename=os.path.join(self.path_log,"valet.log"), filemode="w",level=logging.INFO) 
            #root.addHandler(handler)
            file_handler = logging.FileHandler(os.path.join(self.path_log,"valet.log"))
            stdout_handler = logging.StreamHandler(sys.stdout)
            stdout_handler.setLevel(logging.ERROR)
            handlers = [file_handler, stdout_handler]

            logging.basicConfig(format='%(asctime)s : %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO,handlers=handlers)
            return True
    
    def load_structure(self,xyz_init,xyz_format):
        # Initial structure handling
        if xyz_init == None:
            logging.error("Please provide an initial structure file.")
            return None
        elif os.path.isfile(xyz_init):
            self.xyz_init = xyz_init
            try:
                if xyz_format == None:
                    self.system_init = read(self.xyz_init)
                else:
                    self.system_init = read(self.xyz_int,format=xyz_format)
            except:
                if xyz_format == None:
                    logging.error("Unable to determine file type. Try instantiating system with `xyz_format=type`.")
                    return None
                else:
                    logging.error("Unable to load structure file with format '{}'.".format(xyz_format))
                    return None
        else:
            logging.error("Cannot locate POSCAR file as {}".format(xyz_init))
            return None
        write(self.xyzs['init'],self.system_init)
        print("ValetSystem successfully instantiated.")
            
    #def do_autotune_encut(self):
     #   AutotuneEncut(self)
    #self.do_autotune_encut = AutotuneEncut(self)        

class ValetJob:
    def __init__(self,system,kpts=None):
        if self.attach_system(system):
            self.set_kpts(kpts)

        
    def attach_system(self,system):
        if type(system) != dftvalet.classes.ValetSystem:
            logging.error("Please provide a valid ValetSystem for this job.")
            return False
        else:
            self.valet_system = system
            logging.info("Attached system.")
            print("ValetSystem successfully attached.")
            return True
        
    def clean_up(self):
        os.system("cd {}; rm ase* C* D* E* I* K* O* PC* POSCAR POT* vaspr* WA* X*".format(self.valet_system.path_out))
            
    def set_kpts(self,kpts):
        # Check if user supplied kpoints
        if kpts != None:
            self.kpts = kpts
        else:
            self.kpts = [4,4,4]
        
class AutotuneEncut(ValetJob):
    def __init__(self,system,calc=None,tol=1E-4,max_steps=32,start_encut=200,encut_step=25,retune=False):
        ValetJob.__init__(self,system)
        self.set_calc(calc)
        self.do_encut_autotune(tol=tol,max_steps=max_steps,start_encut=start_encut,encut_step=encut_step,retune=retune)
        
    def set_calc(self,calc=None):
        if calc == None:
            print("Attaching standard calculator.")
            self.valet_system.system_init.set_calculator(Vasp2(xc='PBE',kpts=self.kpts,
                                                               gamma=True,directory=self.valet_system.path_out))
        else:
            print("Attaching custom calculator.")
            self.valet_system.system_init.set_calculator(calc)  
    
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
        if (self.valet_system.encut == None) or (retune):
            sys = self.valet_system.system_init
            calc = sys.calc
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
                self.valet_system.encut = None
            except NameError:
                print("\tENCUT auto-tuned to: {}".format(ENCUT))
                logging.info("Completed ENCUT autotune.")
                self.valet_system.encut = ENCUT
            print("Done.")  
        else:
            print("ENCUT already autotuned. If you would like to retune, please pass 'retune' as True to 'self.do_encut_autotune'.")
            logging.info("Already autotuned ENCUT.")
            
        ValetJob.clean_up(self)
        
class AutotuneVolume(ValetJob):
    def __init__(self,system,calc=None,encut=None,start_scan=0.85,end_scan=1.15,num_scans=15,exclude_type=None,only_type=None):
        ValetJob.__init__(self,system)
        self.set_calc(calc)
        self.encut = (self.valet_system.encut if self.valet_system.encut != None 
                      else AutotuneEncut(system))
        self.do_volume_autotune(encut=self.encut,start_scan=start_scan,end_scan=end_scan,num_scans=num_scans,exclude_type=exclude_type,only_type=only_type)
        
    def set_calc(self,calc=None):
        if calc == None:
            print("Attaching standard calculator.")
            self.valet_system.system_init.set_calculator(Vasp2(xc='PBE',kpts=self.kpts,
                                                               gamma=True,directory=self.valet_system.path_out))
        else:
            print("Attaching custom calculator.")
            self.valet_system.system_init.set_calculator(calc) 

        
     
            
    def set_volume_autotune_params(self,start_scan=None,end_scan=None,num_scans=None,exclude_type=None,
                                   only_type=None):
        vols = []
        ens = []
        return start_scan,end_scan,num_scans,exclude_type,only_type,vols,ens

    def do_volume_autotune(self,encut=None,start_scan=0.85,end_scan=1.15,num_scans=15,exclude_type=None,only_type=None):
        if not os.path.isfile(self.valet_system.xyzs['vol_scan']):
            encut = encut #if encut!=None else self.do_encut_autotune()

            logging.info("Commencing volume scan...")
            print("Determining optimal volume...")
            start_scan,end_scan,num_scans,exclude_type,only_type,vols,ens = self.set_volume_autotune_params(start_scan=start_scan,end_scan=end_scan,num_scans=num_scans,exclude_type=exclude_type,only_type=only_type)


            # Load POSCAR
            #try:
            #    sys = read(self.poscar_init)
            #    calc = self.calc_dict['volume']
            #    calc.set(encut=encut)
            #    start_cell = sys.get_cell()
            #except:
            #    logging.error("Initial POSCAR file could not be found at {}".format(poscar_init))
            #    print("Error: See error log for details.")
            #    return None
            #logging.info("Loaded initial POSCAR file.")

            # Do volume scan
            logging.info("Performing volume scan.")
            sys = self.valet_system.system_init.copy()
            calc = self.valet_system.system_init.calc
            calc.set(encut=encut)
            start_cell=sys.get_cell()
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
                v, e, B = eos.fit()
                eos_fits[typ] = {'volume':v,'energy':e,'buld_modulus':B}

            # Rescale initial cell to optimized volume
            logging.info("Optimal volume found. Rescaling original cell.")
            vol_avg = np.average([eos_fits[key]['volume'] for key in eos_fits.keys()])
            sys.set_cell(start_cell,scale_atoms=True)

            scale = sys.get_volume()/vol_avg
            sys.set_cell(scale*start_cell,scale_atoms=True)

            # Perform sc-step
            logging.info("Performing self-consistent step.")
            sys.set_calculator(calc)
            en = sys.get_potential_energy()

            # Save structure
            logging.info("Saving self-consistent calculator.")
            write(self.valet_system.xyzs['vol_scan'],sys)

            print("Done.")
            logging.info("Volume scan complete.")
            self.clean_up()
        else:
            sys = read(self.valet_system.xyzs['vol_scan'])
            print("Volume optimized system ready.")
        ValetJob.clean_up(self)

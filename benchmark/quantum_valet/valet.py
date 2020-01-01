import os
from ase.calculators.vasp import Vasp2
from ase.io import read
import pickle as pkl
import ase
import numpy as np
from ase.eos import EquationOfState
from ase.dft.bandgap import bandgap

class AutotuneBandgap:
    def __init__(self,system,calc,path_out):
        # Load system and calculator
        if isinstance(system,ase.atoms.Atoms):
            print("System loaded.")
        else:
            print("System must be of type `ase.atoms.Atoms`. Got type `{}`.".format(type(system)))
            return
        if isinstance(calc,ase.calculators.vasp.Vasp2):
            self.system = system
            self.system.set_calculator(calc)
            print("Calculator loaded and attached to system.")
        else:
            print("Calculator must be of type `ase.calculator.vasp.Vasp2`. Got type `{}`.".format(type(calc)))
            return
        self.path_out = os.path.join('.qv_workspace',path_out)
        self.system.calc.set(directory=self.path_out)
        if not os.path.isdir(self.path_out):
            os.system('mkdir -p {}'.format(self.path_out))
            print('Workspace created at {}'.format(self.path_out))

        #self.do_autotune_volume()
        #self.get_bandgap()
        print("Done.")

    def do_autotune_volume(self,start_scan=0.90,end_scan=1.10,num_scan=25,exclude_type=None,only_type=None):

        # Load system and calculator
        system = self.system.copy()
        calc = self.system.calc 
        system.set_calculator(calc)
        # Scale volumes and collect energy & volume measurements
        ens = []
        vols = []
        start_cell = system.get_cell()
        
        print("Performing volume scan.")
        for i,x in enumerate(np.linspace(start_scan,end_scan,num_scan)):
            #print(i)
            scale = x**(1./3.)
            system.set_cell(scale*start_cell,scale_atoms=True)
            vols.append(system.get_volume())
            ens.append(system.get_potential_energy())
        self.volume_vs_energy = [vols,ens]
            
        # Fit to Equations of States
        print("Fitting data to equations of state.")
        eos_types = "sjeos taylor murnaghan birch birchmurnaghan pouriertarantola p3 antonschmidt".split()
        if exclude_type != None:
            _ = [eos_types.remove(i) for i in exclude_type]
        elif only_type != None:
            eos_types = only_type
            
        eos_fits = {}
        errors = []
        
        for typ in eos_types:
            try:
                eos = EquationOfState(vols,ens,eos=typ)
                v, e, B = eos.fit()
                eos_fits[typ] = {'volume':v,'energy':e,'buld_modulus':B}
            except:
                errors.append(typ)

        self.eos_fits = eos_fits
        self.eos_errors = errors if len(errors)>0 else None

        if len(errors) > 0:
            print("WARNING: Was not able to fit equation of state for the following: {}".format(errors))
        # Rescale initial cell to optimized volume
        print("Rescaling with optimized volume.")
        vol_avg = np.average([eos_fits[key]['volume'] for key in eos_fits.keys()])
        system.set_cell(start_cell,scale_atoms=True)
        scale = (vol_avg/system.get_volume())**(1./3.)
        print("Optimized volume: {}".format(vol_avg))
        del self.system
        self.system = system
        self.system.set_cell(scale*start_cell,scale_atoms=True)
        
        print("Performing second relaxation.")
        self.system.get_potential_energy()
        
        #self.clean_up()
        
    def get_bandgap(self):
        print("Calculating bandgap.")
        gap,p1,p2 = bandgap(self.system.calc);
        direct_gap,direct_p1,direct_p2 = bandgap(self.system.calc,direct=True);
        self.bandgap = gap
        self.bandgap_direct = direct_gap

    def clean_up(self):
        os.system("cd {}; rm ase-sort.dat CH* D* E* I* K* OSZ* P* REPORT vasp.out vasprun* X*".format(self.path_out))
    

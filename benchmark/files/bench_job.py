import pickle as pkl
from ase.calculators.vasp import Vasp2
import ase,sys,signal
from quantum_valet.valet import AutotuneBandgap
import numpy as np

i = sys.argv[1]
bench_set = pkl.load(open('bench_set/pickles/bs_part_{}.pkl'.format(i),'rb'))
keys = bench_set.keys()

# Timer
def signal_handler(signum, frame):
    raise Exception("Timed out!")

signal.signal(signal.SIGALRM, signal_handler)
signal.alarm(10800)   # Three hours per job first run

results = []

for j,key in enumerate(keys):
    print(j)
    bs = bench_set[key]
    try:
        # Autotune
        sys = bs['structure_initial']
        setup = bs['potcar']
        calc = Vasp2()
        calc.read_incar('../calculators/INCAR.volume_scan')
        calc.set(xc='pbe',setups=setup,kpts=[5.,5.,5.],gamma=True)
        sys.set_calculator(calc)
        valet = AutotuneBandgap(sys,calc,key)
        try:
            valet.do_autotune_volume()
            error_volume = False
        except:
            error_volume = True
        
        try:
            valet.get_bandgap()
            error_bandgap = False
        except:
            error_bandgap = True
        
        if (not error_volume) and (not error_bandgap):
            result = {'system':valet.system,'bandgap':valet.bandgap,'bandgap_direct':valet.bandgap_direct,'eos_fits':valet.eos_fits,'eos_errors':valet.eos_errors,'volume_vs_energy':valet.volume_vs_energy,'index':key}  
        elif not error_volume:
            result = {'system':valet.system,'bandgap':None,'bandgap_direct':None,'eos_fits':valet.eos_fits,'eos_errors':valet.eos_errors,'volume_vs_energy':valet.volume_vs_energy,'index':key}  
        else:
            result = {'system':None,'bandgap':None,'bandgap_direct':None,'eos_fits':None,'eos_errors':None,'volume_vs_energy':None,'index':key}
        valet.clean_up()
    except Exception:
        # Timeout Error
        result = {'index':key,'timeout':True}

    results.append(result)
   
    with open('cp_bs_part_{}.pkl'.format(i),'wb') as fout:
        pkl.dump(results,fout)

with open('results_bs_part_{}.pkl'.format(i),'wb') as fout:
    pkl.dump(results,fout)

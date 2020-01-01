import pickle as pkl
from ase.calculators.vasp import Vasp2
import ase
from quantum_valet.valet import AutotuneBandgap
import numpy as np

bench_set = pkl.load(open('bench_set/pickles/bs_part_0.pkl','rb'))
key = 'mp-1'
bs = bench_set[key]

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
    results = {'system':valet.system,'bandgap':valet.bandgap,'bandgap_direct':valet.bandgap_direct,'eos_fits':valet.eos_fits,'eos_errors':valet.eos_errors,'volume_vs_energy':valet.volume_vs_energy,'index':key}  
elif not error_volume:
    results = {'system':valet.system,'bandgap':None,'bandgap_direct':None,'eos_fits':valet.eos_fits,'eos_errors':valet.eos_errors,'volume_vs_energy':valet.volume_vs_energy,'index':key}  
else:
    results = {'system':None,'bandgap':None,'bandgap_direct':None,'eos_fits':None,'eos_errors':None,'volume_vs_energy':None,'index':key}

print(results)  

    

 #valet.get_bandgap()
 #
 #
 #bench_sys = bs['structure_final']
 #bench_vol = bench_sys.get_volume()
 #bench_bandgap = bs['bandgap']
 #
 #valet_vol = valet.system.get_volume()
 #valet_sys = valet.system.copy()
 #
 #
 #results = {'bench_bandgap':bench_bandgap,'bench_vol':bench_vol,'bench_sys':bench_sys,'valet_vol':valet_vol,
 #          'valet_sys':valet_sys,'valet_bandgap':valet.bandgap,'valet_bandgap_direct':valet.bandgap_direct}
 #
 #
 #
 #print(results)

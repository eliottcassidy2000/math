"""Independent inherited-Python/raw-cut support check, plus prefix evidence."""
import importlib.util
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited_compiler',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
cores=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]

def require(p,msg):
    if not p:raise RuntimeError(msg)

def main():
    bound=4095
    with tempfile.TemporaryDirectory(prefix='chanpak-support-') as tmp:
        binary=Path(tmp)/'support.exe'
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'support.cpp'),'-o',str(binary)],check=True,capture_output=True,text=True)
        subprocess.run([str(binary),str(HERE/'signatures.txt'),str(bound),str(HERE/'support_12.json')],check=True,capture_output=True,text=True)
    r=json.loads((HERE/'support_12.json').read_text())
    mask=[0]*(bound//2+1);pairs=0
    for ci,core in enumerate(cores):
        for O in range(1<<len(core)):
            for h,A,B in c.generate(c.signature(core,O),bound):
                mask[h//2]|=1<<ci;pairs+=1
    want={str(k):Counter(mask)[k] for k in range(32)}
    require(r['support_mask_counts']==want,'raw Python cuts versus C++ dedup')
    require(r['five_core_missing']==[2*i+1 for i,z in enumerate(mask) if not z],'five misses')
    require(r['small_parent_extra_holes']==[2*i+1 for i,z in enumerate(mask) if not(z&7) and 2*i+1 not in (7,21)],'small parent holes')
    old=json.loads((HERE.parents[1]/'session251_creative_completion_20261111/ternary/explore-results.json').read_text())
    for bits in (20,24):
        large=json.loads((HERE/f'support_{bits}.json').read_text())
        require([n for n in large['small_parent_extra_holes'] if n<=old['bound']]==old['targets_requiring_parent9_in_this_bank'],'inherited prefix')
        require(large['five_core_missing']==[7,21],'new five-core miss')
    return {'status':'PASS','independent_raw_cut_bound':bound,
            'raw_python_parameter_pairs':pairs,'all_32_per_core_support_counts_match':True,
            'inherited_prefix_bound':old['bound'],
            'scope':'Independent complete support through4095; higher runs retain C++ exact inverse and inherited prefix controls, not an independent high-range replay.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))

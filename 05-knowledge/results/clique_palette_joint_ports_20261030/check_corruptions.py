"""Ensure damaged normality certificates are rejected; never modify the real one."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import copy
import importlib.util
import json

HERE=Path(__file__).resolve().parent

def run():
    spec=importlib.util.spec_from_file_location('normality_audit',HERE/'check_normality.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    original=json.loads((HERE/'normality_certificate.json').read_text())
    cases={}
    bad=copy.deepcopy(original); bad['fundamental_points'].pop(); cases['omitted_lattice_class']=bad
    bad=copy.deepcopy(original); bad['simplices'].pop(); cases['missing_simplex']=bad
    bad=copy.deepcopy(original); bad['fundamental_points'][0]['colors'][0]=0
    cases['damaged_palette']=bad
    bad=copy.deepcopy(original); bad['vertex_masks'][0]=1; cases['wrong_support_inventory']=bad
    rejected={}
    with TemporaryDirectory(prefix='palette-certificate-hostiles-') as temp:
        module.HERE=Path(temp)
        for name,cert in cases.items():
            (module.HERE/'normality_certificate.json').write_text(json.dumps(cert),encoding='utf-8')
            try:
                module.geometry()
            except RuntimeError as error:
                rejected[name]=str(error)
            else:
                raise RuntimeError('Damaged certificate accepted: '+name)
    # Recheck the repaired THM-658 coloring itself, without importing its old script.
    offsets=[2,11,0]
    color=lambda x:(offsets[x%3]+7*(x//3))%27
    speeds=list(range(1,12))+[13,24]
    distances=[min((color(x+v)-color(x))%27,(color(x)-color(x+v))%27)
               for x in range(3) for v in speeds]
    if min(distances)<2:
        raise RuntimeError('Circular coloring certificate fails')
    return {'damaged_certificates_rejected':rejected,'native_GW_phase_checks':len(distances),
            'native_GW_minimum_circular_distance':min(distances),
            'scope':'Negative certificate controls and the explicit 27/2 upper-bound witness only.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); args=ap.parse_args()
    text=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')

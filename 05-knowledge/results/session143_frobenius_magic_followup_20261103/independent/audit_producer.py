#!/usr/bin/env python3
"""Replay producer independently and compare its certificate to direct audit."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(c,m):
    if not c:raise ValueError(m)


def main():
    p=argparse.ArgumentParser();p.add_argument('--producer',type=Path,default=Path(__file__).resolve().parent.parent/'producer');args=p.parse_args()
    independent=Path(__file__).resolve().parent
    audit=json.loads(subprocess.check_output([sys.executable,str(independent/'check_independent.py')],text=True))
    producer=args.producer.resolve()
    outputs=[]
    for mode in ([],['-O']):
        data=subprocess.check_output([sys.executable,*mode,str(producer/'check_bridge.py')])
        need(data==(producer/('result.optimized.json' if mode else 'result.json')).read_bytes(),'producer replay bytes differ')
        outputs.append(json.loads(data))
    need(outputs[0]==outputs[1],'producer optimization mismatch')
    prod=outputs[0]
    a={(r['a'],r['b']):(r['t'],r['N1'],r['N2'],r['N3']) for r in audit['models']}
    b={(r['a'],r['b']):(r['t'],r['n5'],r['n25'],r['n125']) for r in prod['models']}
    need(a==b and len(a)==20,'different-field independent point counts mismatch')
    need(prod['magic_square']==audit['magic_grid'] and prod['line_sums']==audit['magic_line_sums'],'magic grid mismatch')
    need(prod['magic_arrangement_count']==audit['magic_arrangement_counts']['N3'],'independent 8! arrangement count mismatch')
    need(prod['class_trace_fibers']==audit['trace_fibers'],'trace fiber mismatch')
    need(prod['twist_orbits']==audit['twist_orbits']==6,'twist orbit mismatch')
    need(prod['tournament_spectrum_fibers']==audit['tournament_control']['spectrum_fiber_partition'],'Leibniz/Newton characteristic polynomial mismatch')
    for mode in ([],['-O']):
        subprocess.check_output([sys.executable,*mode,str(producer/'test_controls.py')])
        subprocess.check_output([sys.executable,*mode,str(producer/'verify_manifest.py')])
    paths=['PROOF.md','check_bridge.py','result.json','result.optimized.json','source_pins.json','test_controls.py','manifest.json']
    pins={p:hashlib.sha256((producer/p).read_bytes()).hexdigest() for p in paths}
    print(json.dumps({'status':'PASS','producer_pins':pins,'pointcount_comparisons':20,'independent_field_modulus':'z^3+2z+1','producer_field_modulus':'z^3+z+1','arrangement_method':'8! center-fixed permutations','tournament_polynomial_paths':['auditor Newton identities','producer Leibniz determinant'],'producer_normal_optimized_replay':'byte-identical to both saved certificates','producer_controls_and_manifest':'passed normal and optimized'},indent=2,sort_keys=True))

if __name__=='__main__':main()

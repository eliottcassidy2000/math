#!/usr/bin/env python3
"""Reproduce the full portable atlas/obstruction package in both Python modes."""
import hashlib,json,subprocess,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def run(script,args=(),optimized=False):
    return subprocess.check_output([sys.executable]+(['-O'] if optimized else [])+[str(ROOT/script)]+list(args))
def main():
    outputs={};sequence=[]
    for script,frozen in [('build_atlas.py','atlas.json'),('check_basis_independent.py','independent.normal.json')]:
        a=run(script);b=run(script,optimized=True)
        need(a==b,script+' mode disagreement');need(a==(ROOT/frozen).read_bytes(),script+' frozen disagreement')
        outputs[script]=hashlib.sha256(a).hexdigest()
    data={}
    for B in range(1,22,2):
        args=['--threshold',str(B)];a=run('forbidden_basis.py',args);b=run('forbidden_basis.py',args,True)
        need(a==b,'basis mode disagreement at '+str(B));obj=json.loads(a);data[B]=obj
        if B in (7,21):need(a==(ROOT/('basis'+str(B)+'.json')).read_bytes(),'frozen basis disagreement')
        sequence.append({'threshold':B,**obj['counts'],'maximum_order':obj['sharp_general_maximum_order']})
    need(sequence==json.loads((ROOT/'basis_sequence.json').read_text()),'sequence disagreement')
    for a,b in ((5,7),(19,21)):
        for key in ('strong_obstructions','nonstrong_obstructions'):
            need(data[a][key]==data[b][key],'equal hereditary classes failed')
    print(json.dumps({'status':'PASS','python_modes':['normal','optimized'],
        'thresholds':list(range(1,22,2)),'all_order_basis_counts':[row['total'] for row in sequence],
        'basis21_classes':data[21]['counts']['total'],'independent_cycle_cover_presentations':512+16384,
        'literal_independent_checks':434,'output_sha256':outputs},indent=2,sort_keys=True))
if __name__=='__main__':main()

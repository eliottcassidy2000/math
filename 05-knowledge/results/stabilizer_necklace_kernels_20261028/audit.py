"""Same-agent alternate paths; deterministic receipts and inherited hostile control."""
from collections import Counter
from pathlib import Path
from math import prod
import hashlib
import importlib.util
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    runs=[]
    for name in ['check_stabilizers','check_kernels','check_cokernels']:
        script=HERE/(name+'.py'); output=script.with_suffix('.json'); baseline=output.read_bytes()
        for flags in [[],['-O']]:
            cp=subprocess.run([sys.executable,'-B',*flags,str(script)],capture_output=True,text=True)
            need(cp.returncode==0,cp.stderr)
            need(output.read_bytes()==baseline,name+' normal/-O exact receipt')
            runs.append({'script':name,'optimized':bool(flags),'checks':json.loads(baseline)['checks']})
    kernels=json.loads((HERE/'check_kernels.json').read_text())['rows']
    coker=json.loads((HERE/'check_cokernels.json').read_text())['rows']
    stabs=json.loads((HERE/'check_stabilizers.json').read_text())['rows']
    for k,row in enumerate(coker):
        need(row['determinant']==kernels[k]['determinant'],'Smith product equals Bareiss/product recurrence')
        need(row['N']-kernels[k]['finite_field_ranks']['19']==row['odd_Smith_factors'],'F19 loss equals odd factors')
        if k<len(stabs): need(row['odd_Smith_factors']==stabs[k]['minimum_C_regular_subset'],'shared depth sum')
    oracle_path=HERE.parent/'relation_tower_frobenius_ports_20261028/check_root_switch.py'
    oracle=load(oracle_path,'hostile_oracle'); hostile=oracle.hostile(); spectra=[]
    for key in ['rows_AAA','rows_BAA']:
        group=oracle.isomorphisms(hostile[key],hostile[key],True); seen=set(); counts=Counter()
        for mask in range(4096):
            if mask in seen: continue
            orbit={oracle.move_mask(mask,g) for g in group}; seen.update(orbit)
            stabilizer=len(group)//len(orbit); counts[stabilizer]+=1
        need(sum(v*len(group)//s for s,v in counts.items())==4096,'hostile labelled mass')
        spectra.append(dict(sorted(counts.items())))
    need([sum(p.values()) for p in spectra]==[176,512],'81/27 orbit census')
    receipt={'audit':'same-agent, alternate exact computations; no independent referee or Lean claim',
             'runs':runs,'cross_checks':'Smith/determinant/rank/minimum regular-set depth sums agree',
             'hostile_81_27':{'automorphism_orders':hostile['automorphism_orders'],
                  'cyclic_triples':hostile['labelled_cyclic_triples'],'Hamilton_paths':hostile['Hamilton_paths'],
                  'binary_orbits_by_stabilizer_order':spectra},'hashes':{}}
    paths=list(HERE.glob('check_*.py'))+list(HERE.glob('check_*.json'))+[oracle_path,
        ROOT/'01-canon/theorems/THM-2172-frobenius-collapse-of-mahler-and-twisted-cyclic-packets.md',
        HERE.parent/'session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md']
    receipt['hashes']={str(f.relative_to(ROOT)):digest(f) for f in paths}
    (HERE/'audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    broken=[]
    for file in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'): continue
            if not (file.parent/target.split('#')[0]).resolve().exists(): broken.append((file.name,target))
    need(not broken,str(broken))
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()

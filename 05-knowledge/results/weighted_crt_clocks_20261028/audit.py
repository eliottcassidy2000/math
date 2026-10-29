"""Same-agent alternate checks and reproducible normal/-O receipts."""
from pathlib import Path
from math import gcd
import hashlib
import importlib.util
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def mod(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def main():
    runs=[]
    for name in ['check_crt','check_weighted','check_affine','check_eta']:
        script=HERE/(name+'.py'); result=script.with_suffix('.json')
        expected=result.read_bytes()
        for flags in [[],['-O']]:
            p=subprocess.run([sys.executable,'-B',*flags,str(script)],capture_output=True,text=True)
            need(p.returncode==0,p.stderr)
            need(result.read_bytes()==expected,name+' deterministic normal/optimized receipt')
            runs.append({'script':name,'optimized':bool(flags),'checks':json.loads(expected)['checks']})
    crt=mod('check_crt'); algebra_checks=0
    for n in range(3,302,2):
        a,b=crt.fib(n-1)-1,crt.fib(n); ln=2*crt.fib(n-1)+b
        need(a*a+a*b-b*b==-ln,'independent norm')
        need(gcd(a,b)==(2 if n%3==0 else 1),'Smith first factor')
        m=ln//4 if n%3==0 else ln
        r=(-a*pow(b,-1,m))%m
        need((r*r-r-1)%m==0 and pow(r,n,m)==1%m,'root and period')
        algebra_checks+=3
    n=55; a,b=crt.fib(n-1)-1,crt.fib(n); ln=2*crt.fib(n-1)+b
    fs=crt.factors(ln); r=(-a*pow(b,-1,ln))%ln
    need(121 in fs,'repeated odd prime factor in period55')
    need(pow(r,55,121)==1 and pow(r,5,121)!=1 and pow(r,11,121)!=1,'period55 modulo121')
    weighted=mod('check_weighted'); g=weighted.permutations_a(2); t=weighted.negation(2)
    c=[p for p in g if all(p[t[i]]==t[p[i]] for i in range(9))]
    # Distinct signed weights, constant on each of the five single-port orbits.
    w=[3,2,-5,7,7,7,-4,-4,-4]
    mean_square=sum(sum(w[i] for i in range(9) if p[i]==i)**2 for p in c)
    need(mean_square==len(c)*((3+2-5+7-4)**2+2*(7-4)**2),'weighted shared-port variance')
    # Pair kernels and joint phase counts checked independently of word decoding.
    for row in json.loads((HERE/'check_crt.json').read_text())['censuses']:
        n=row['n']
        necklaces=sum(sum(gcd(i,d)==1 for i in range(1,d+1))*(2*crt.fib(n//d-1)+crt.fib(n//d))
                      for d in range(1,n+1) if n%d==0)//n
        need(necklaces==row['global_rotation_orbits'],'independent necklace Burnside')
    broken=[]
    for file in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',file.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'): continue
            if not (file.parent/target.split('#')[0]).resolve().exists(): broken.append((file.name,target))
    need(not broken,str(broken))
    dependencies=[HERE.parent/'relation_tower_frobenius_ports_20261028/check_root_switch.py',
                  HERE.parents[2]/'04-computation/a000568_exact_burnside.py',
                  HERE.parent/'lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md',
                  HERE.parent/'lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md',
                  HERE.parent/'cursor_cm_phase_20261028/TRACE_CARRIER.md']
    receipt={'audit':'same-agent; alternate arithmetic paths; no formal or independent referee claim',
             'runs':runs,'odd_ring_algebra_checks':algebra_checks,
             'period55_algebra_only':{'states':ln,'factors':fs,'root':r,'order_mod_121':55},
             'weighted_variance_literal_group':True,'necklace_cross_check':True,'local_links':'all resolve',
             'hashes':{str(p.relative_to(HERE.parents[2])):digest(p) for p in dependencies}}
    for p in sorted(HERE.glob('check_*')):
        if p.suffix in ['.py','.json']: receipt['hashes'][p.name]=digest(p)
    (HERE/'audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()

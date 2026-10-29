"""Independent native controls for the U128 bounded exclusion engine.

The complete horizon-74 run is not repeated or independently audited here.
"""
import argparse, hashlib, json, subprocess
from pathlib import Path
from reverse_reset import native
P=Path(__file__).resolve().parent

def main(binary):
    checks=0
    def need(condition,msg):
        nonlocal checks
        checks+=1
        if not condition:raise ValueError(msg)
    def run(k,b):return json.loads(subprocess.check_output([binary,str(k),str(b)],text=True))
    mins=[None]*17
    # Original source scan, no affine congruence recurrence and no branch pruning.
    for source in range(1,150000,2):
        n=source
        for k in range(1,17):
            h,w=native(n)
            if h<=n:break
            if mins[k] is None:mins[k]=source
            n=h
    controls=[]
    for k in range(1,17):
        need(mins[k] is not None,('finite scan found source',k))
        negative=run(k,mins[k]);positive=run(k,mins[k]+1)
        need(not negative['counterexample'],('exclusive minimum',k))
        need(positive['counterexample'] and positive['source']==mins[k],('inclusive witness',k))
        controls.append(dict(k=k,minimum=mins[k],visited=positive['visited']))
    # Higher-horizon control inherited from the independently audited total-affine search.
    old=json.loads((P.parent/'atomic_tail_reset_20261029/least60.json').read_text())['rows']
    positive=run(55,old[54]['least_source']+1)
    need(positive['source']==old[54]['least_source'],'late branch positive witness')
    n=positive['source'];J=D=0
    words=((1,),(3,1,1,1),(3,3,1,1,1,1,1),(3,1,3,1,1,1,1),(3,1,1,3,1,1,1))
    for letter in positive['word']:
        h,w=native(n);need(h>n and tuple(w)==words[letter],'literal high witness')
        J+=len(w);D+=sum(w);n=h
    need(n==int(positive['endpoint']),'literal high endpoint')
    need(4 in positive['word'],'all five branches remain available')
    receipt=json.loads((P/'bound74.json').read_text())
    return dict(status='PASS',explicit_checks=checks,native_universe='all odd 1<=n<150000, first 16 growing returns',
                small_controls=controls,high_control=positive,high_odd_clock=J,high_division_clock=D,
                bound74=receipt,source_sha256=hashlib.sha256((P/'window74.cpp').read_bytes()).hexdigest(),
                audit_scope='Independent native scan through 16 and audited minimum control at 55; no independent full horizon-74 audit')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('binary');print(json.dumps(main(p.parse_args().binary),indent=2,sort_keys=True))

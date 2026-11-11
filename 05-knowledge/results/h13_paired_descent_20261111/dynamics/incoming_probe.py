"""Exact literal401 repair controls; inherited sources are not re-audited."""
from pathlib import Path
import importlib.util,json

def main():
    p=Path(__file__).resolve().parents[1]/'modular_cycles/probe.py'
    s=importlib.util.spec_from_file_location('native',p)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    records=[]
    for cut in (6,3):
        t=m.c.construct((18,8,3,5,14),cut,0,0)
        trace=[m.state(t,5)]
        for mask in (0,2,32):
            t=m.add(t,mask);trace.append(m.state(t,5))
        records.append({'root':5,'rows':t,'trace':trace})
    expected=[[(31,7,45,13),(31,0,90,13),(95,7,313,41),(401,41,585,41)],
              [(33,7,45,13),(33,0,90,13),(101,7,307,41),(401,41,573,41)]]
    if [r['trace'] for r in records]!=expected:raise RuntimeError(records)
    x,y=(r['rows'] for r in records)
    arcs=[(a,b) for a in range(9) for b in range(a+1,9) if ((x[a]^y[a])>>b)&1]
    if arcs!=[(0,5),(2,5)]:raise RuntimeError(arcs)
    report={'status':'FINITE-EXACT named source constants; all-context consequence proved in INCOMING_REPAIR.md',
            'universe':'two specified initial6-vertex graphs, three specified births, eight rooted states by literal subset DP',
            'records':records,'changed_unordered_edges':arcs,
            'known_boundary':'The same schedule at H13 addressB1 has middlecounts265/265, so no uniform ternary-branch claim.'}
    Path(__file__).with_name('incoming_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('literal401 repair:',[r['trace'][-1] for r in records],'; changed edges:',arcs)

if __name__=='__main__':main()

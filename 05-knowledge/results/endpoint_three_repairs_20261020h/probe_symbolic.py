"""Affine binary compiler retaining a TWO-element exceptional port set.

Infeasibility is exact for each stated order, marks, and arc guards; no
all-order conclusion is inferred from a finite list of such checks.
"""
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import probe as P

OLD_SYMBOLIC=P.OLD_PATH.with_name('probe_symbolic.py')
spec=importlib.util.spec_from_file_location('endpoint_two_affine',OLD_SYMBOLIC)
OLD=importlib.util.module_from_spec(spec)
spec.loader.exec_module(OLD)


def blue_compatible(n,incoming,words,step,extra_arcs=()):
    incoming=set(incoming)
    outgoing={n-1-v for v in incoming}
    P.need(len(incoming)==2 and incoming<=set(range(2,n-1)),'two admissible ports')
    pairs=list(combinations(range(n),2))
    indices={pair:i for i,pair in enumerate(pairs)}
    affine=OLD.Affine()
    def fixed(u,v,bit,reason):
        affine.add(1<<indices[P.edge(u,v)],bit^int(u>v),reason)
    for i in range(n-1):
        fixed(i,i+1,1,f'path:{i}')
    for u,v in pairs:
        image=(n-1-v,n-1-u)
        affine.add((1<<indices[(u,v)])^(1<<indices[image]),0,f'blueP:{u},{v}')
    for v in range(1,n):
        fixed(0,v,int(v not in incoming|{n-1}),f'first:{v}')
    for v in range(n-1):
        fixed(n-1,v,int(v in outgoing|{0}),f'last:{v}')
    for i,(u,v,bit) in enumerate(extra_arcs):
        fixed(u,v,bit,f'extra:{i}:{u},{v}')
    compiled=P.compile_ledger(n,words)
    for u,v,bit in compiled['guards']:
        fixed(u,v,bit,f'direction:{u},{v}')
    if compiled['guard_conflicts']:
        return {'compatible':False,'reason':'internally inconsistent direction guards'}
    prior=set()
    for word in words[:step]:
        prior^=P.edges(word)
    q=words[step]
    for i,j in combinations(range(n),2):
        u,v=q[i],q[j]
        a,b=q[n-1-j],q[n-1-i]
        pair,image=P.edge(u,v),P.edge(a,b)
        rhs=int(u>v)^int(a>b)^int(pair in prior)^int(image in prior)
        affine.add((1<<indices[pair])^(1<<indices[image]),rhs,f'blueQ:{i},{j}')
    if affine.contradiction is not None:
        return {'compatible':False,'contradiction_equations':affine.contradiction}
    assignment=affine.solution()
    rows=[0]*n
    for i,(u,v) in enumerate(pairs):
        rows[u if assignment>>i&1 else v] |= 1<<(v if assignment>>i&1 else u)
    rows=tuple(rows)
    P.source_checks(rows,tuple(sorted(incoming)),tuple(sorted(outgoing)))
    now=rows
    for i,word in enumerate(words):
        P.need(P.path(now,word),'affine ledger directed')
        if i==step:
            P.need(P.blue(now,word),'affine selected step blue')
        now=P.complement(now,word)
    P.need(now==P.complement(rows,tuple(range(n))),'affine literal target')
    return {'compatible':True,'rows':rows,'scores':[row.bit_count() for row in rows]}


def B_words(n,cut=None):
    cut=n-2 if cut is None else cut
    a,b,x,z=0,1,2,n-1
    d,r=tuple(range(3,cut)),tuple(range(cut,n-1))
    P.need(d and r,'nonempty B blocks')
    return ((a,b)+r+(z,x)+d,(b,)+r+(a,z,x)+d,(b,x)+d+r+(a,z))


def main():
    records=[]
    for n in range(8,21):
        for h in sorted({2,4,6,n-6}):
            if h not in range(2,n-2) or h==n-3:
                continue
            for step in range(3):
                answer=blue_compatible(n,{h,n-3},B_words(n),step,[(1,n-2,1)])
                records.append({'n':n,'h':h,'step':step+1,**answer})
                if answer['compatible']:
                    print('BLUE COMPATIBLE',n,h,step+1,answer['rows'],flush=True)
    receipt={'status':'FINITE-EXACT affine feasibility with explicit port and open-middle guards',
             'source_hashes':{'current_symbolic':sha256(Path(__file__).read_bytes()).hexdigest(),
                              'current_probe':sha256(Path(P.__file__).read_bytes()).hexdigest(),
                              'inherited_endpoint_two_symbolic':sha256(OLD_SYMBOLIC.read_bytes()).hexdigest()},
             'records':records}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('checked',len(records),'finite systems',flush=True)


if __name__=='__main__':
    main()

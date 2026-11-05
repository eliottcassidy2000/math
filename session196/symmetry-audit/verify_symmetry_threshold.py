#!/usr/bin/env python3
"""Independent finite controls for the 4151 nonrigid-7p threshold.
No repository imports beyond this audit's standalone primitive verifier.
"""
import itertools
import json
from pathlib import Path
from verify_prime_rigidity import adjacency, hamilton_paths, automorphisms


def strong(out):
    all_vertices=(1<<len(out))-1
    for v in range(len(out)):
        seen=1<<v
        while True:
            new=seen
            for i,row in enumerate(out):
                if seen>>i&1: new|=row
            if new==seen: break
            seen=new
        if seen!=all_vertices: return False
    return True


def literal_paths(out):
    return sum(all(out[a]>>b&1 for a,b in zip(p,p[1:]))
               for p in itertools.permutations(range(len(out))))


def circulant(p,choice):
    steps=[i if choice>>(i-1)&1 else p-i for i in range(1,(p+1)//2)]
    return [sum(1<<((i+s)%p) for s in steps) for i in range(p)]


def inflate(q,r,m):
    blocks=[]; count=0
    for i in range(len(q)):
        b=list(range(count,count+(len(m) if i==r else 1)))
        blocks.append(b); count+=len(b)
    out=[0]*count
    for i,bi in enumerate(blocks):
        for j,bj in enumerate(blocks):
            if q[i]>>j&1:
                for x in bi:
                    for y in bj: out[x]|=1<<y
        if i==r:
            for a,x in enumerate(bi):
                for b,y in enumerate(bi):
                    if m[a]>>b&1: out[x]|=1<<y
    return out


def main():
    result={}
    strong5=[adjacency(5,k) for k in range(1024) if strong(adjacency(5,k))]
    counts={}
    for t in strong5:
        h=hamilton_paths(t)
        assert h==literal_paths(t)
        counts[h]=counts.get(h,0)+1
    assert len(strong5)==544 and min(counts)==9
    result['all_labelled_strong5']={'count':len(strong5),'spectrum_counts':dict(sorted(counts.items())),
        'checks':'All 1024 orientations considered; every strong case checked by subset DP and literal permutations.'}
    strong4=[adjacency(4,k) for k in range(64) if strong(adjacency(4,k))]
    assert len(strong4)==24
    for p in (3,5,7):
        rows=[]
        for choice in range(1<<((p-1)//2)):
            m=circulant(p,choice)
            h=hamilton_paths(m)
            assert h==literal_paths(m)
            c3=hamilton_paths(inflate([2,4,1],0,m))
            c4=sorted({hamilton_paths(inflate(q,r,m)) for q in strong4 for r in range(4)})
            rows.append({'choice':choice,'rows':m,'H':h,'cyclic3_inflation_H':c3,'rooted_strong4_inflation_H':c4})
        result[f'all_circulant_{p}_connection_sets']=rows
    eleven=[hamilton_paths(circulant(11,choice)) for choice in range(32)]
    assert min(eleven)==92411
    result['all_circulant_11_H']=eleven
    two_orbits=[]
    for orientation in (1,2):
        for b in range(8):
            out=[1<<((i+1)%3) for i in range(3)]+[1<<(3+(i+orientation)%3) for i in range(3)]
            for i in range(3):
                for j in range(3):
                    if b>>((j-i)%3)&1: out[i]|=1<<(j+3)
                    else: out[j+3]|=1<<i
            h=hamilton_paths(out)
            assert h==literal_paths(out)
            assert h==(9 if b in (0,7) else 45)
            two_orbits.append({'second_cycle_step':orientation,'cross_mask':b,'H':h})
    result['all_two_C3_orbit_cases']=two_orbits
    witness=[240,225,195,135,142,156,184,512,255,383]
    assert strong(witness)
    h=hamilton_paths(witness); aut=automorphisms(witness)
    assert h==4151 and aut==7
    assert all(593%d for d in range(2,25))
    result['witness']={'rows':witness,'H':h,'automorphism_order':aut,'strong':True,'factorization':[7,593]}
    output=Path(__file__).with_name('symmetry_threshold_exact.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()

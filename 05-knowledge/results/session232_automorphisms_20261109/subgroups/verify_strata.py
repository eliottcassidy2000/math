#!/usr/bin/env python3
from collections import Counter
from itertools import permutations
from pathlib import Path
import hashlib,json
import verify as v

OUT=Path(__file__).resolve().parent

def induced(rows,vertices):
    return tuple(sum(((rows[u]>>w)&1)<<j for j,w in enumerate(vertices)) for u in vertices)

def two_covers(rows):
    n=len(rows); total=0
    for mask in range(1,(1<<n)-1):
        a=[u for u in range(n) if mask>>u&1];b=[u for u in range(n) if not mask>>u&1]
        total+=v.hpaths(induced(rows,a))*v.hpaths(induced(rows,b))
    literal=0
    for p in permutations(range(n)):
        for k in range(1,n):
            if all(rows[p[i]]>>p[i+1]&1 for i in range(n-1) if i!=k-1):literal+=1
    v.require(total==literal,'ordered two-path covers independent check')
    return total

def cyclic_extension(rows):
    n=len(rows);ans=[r|(1<<n) for r in rows]
    ans.extend((1<<(n+1),(1<<n)-1))
    return tuple(ans)

def main():
    controls=[]
    for name,n,steps,h,f2,aut,q in [('R5',5,[1,2],15,90,5,3),('P7',7,[1,2,4],189,1638,21,9)]:
        rows=v.cayley(n,steps)
        v.require(v.hpaths(rows)==h,'base H')
        v.require(two_covers(rows)==f2,'base F2')
        extension=cyclic_extension(rows);he=v.hpaths(extension);ae=len(v.automorphisms(extension))
        v.require(he==3*h+f2,'exact singleton cyclic extension')
        v.require(ae==aut,'extension full automorphism group')
        v.require(he//ae==(27 if n==5 else 105),'module interaction quotient')
        controls.append({'name':name,'rows':rows,'H':h,'F2_nonempty_ordered':f2,'q':q,
            'extension_rows':extension,'extension_H':he,'extension_aut_order':ae,'extension_q':he//ae})
    tests=[]
    for n in range(1,8):
        counts=Counter()
        for rows in v.normalized_census(n):
            g=v.automorphisms(rows);h=v.hpaths(rows)
            for a in v.subgroups(g):
                q=h//len(a)
                if q<23:
                    price=3**(v.valuation(len(a),5)+2*v.valuation(len(a),7))
                    v.require(q==price or q>=3*price,'threefold separation')
                    counts['threefold_separation']+=1
                if len(a)%7==0 and q<23:
                    v.require(q==9 and len(a)==len(g) and v.equality_stack(rows,0,1),'sharp seven stratum')
                    counts['seven']+=1
                if len(a)%5==0 and q<9:
                    v.require(q==3 and len(a)==len(g) and v.equality_stack(rows,1,0),'sharp five stratum')
                    counts['five']+=1
        tests.append({'n':n,'instances':dict(sorted(counts.items()))})
    double5=v.ordinal_sum(v.cayley(5,[1,2]),v.cayley(5,[1,2]))
    dg=v.automorphisms(double5);dh=v.hpaths(double5)
    v.require(len(dg)==25 and dh==225 and dh//len(dg)==9 and v.equality_stack(double5,2,0),'twenty-five stratum')
    data={'status':'FINITE-EXACT','twenty_five_control':{'H':dh,'aut_order':len(dg),'q':dh//len(dg)},'controls':controls,'complete_subgroup_tests':tests,
          'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'shared_verifier_sha256':hashlib.sha256((OUT/'verify.py').read_bytes()).hexdigest()}
    out=json.dumps(data,sort_keys=True,indent=2)+'\n';(OUT/'strata_receipt.json').write_text(out);print(out,end='')
if __name__=='__main__':main()

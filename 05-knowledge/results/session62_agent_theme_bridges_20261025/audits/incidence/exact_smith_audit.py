#!/usr/bin/env python3
"""Exact Smith invariants of the FULL physical matrix, not reduced equations."""
import argparse, json, pathlib
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
from direct_audit import incidence, require
OUT=pathlib.Path(__file__).resolve().parent

def forest(p,L):
    require(1<=L<p, f"forest requires 1<=L<p, got p={p}, L={L}")
    I=set(range(L)); H={j for j in I if (j+9)%p in I}
    G=set(range(p))-I
    adj={v:set() for v in range(p)}
    for j in G:
        a,b=(j-6)%p,(j-3)%p
        adj[a].add(b); adj[b].add(a)
    seen=set(); f=o=0; compatible=True
    for root in range(p):
        if root in seen: continue
        todo=[root]; col={root:0}; seen.add(root)
        while todo:
            a=todo.pop()
            for b in adj[a]:
                if b in col:
                    require(col[a]!=col[b], f"non-bipartite graph at p={p}, L={L}, edge={(a,b)}")
                else: col[b]=1-col[a]; seen.add(b); todo.append(b)
        anchors={color for v,color in col.items() if v not in H}
        if not anchors:
            f+=1
            o |= len(col)%2
        elif len(anchors)>1:
            compatible=False
    return {'free_components':f,'has_odd_free':bool(o),'c_compatible':compatible,
            'generic_dim':L-1+f-o,'char_p_dim':L-1+f-o+int(compatible and o)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=pathlib.Path, default=OUT)
    parser.add_argument("--negative-control", action="store_true", help="Inject one deliberately false forest-dimension claim; the audit must fail.")
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results=[]
    for p in [5,7,11,13,17,19,23,29,31,37,43]:
        for L in range(1,p+1):
            a=incidence(p,L)
            s=smith_normal_form(Matrix(a),domain=ZZ)
            inv=[abs(int(s[j,j])) for j in range(min(s.shape)) if s[j,j]]
            generic=4*L-len(inv)
            torsion=[v for v in inv if v!=1]
            expected=forest(p,L) if L<p else {'generic_dim':2*p-2,'char_p_dim':2*p-1}
            if args.negative_control and (p,L)==(5,1):
                expected['generic_dim']+=1  # Deliberate hostile claim.
            require(generic==expected['generic_dim'], f"full-matrix rank mismatch: p={p}, L={L}, dim={generic}, expected={expected}")
            require(torsion==([p] if expected['char_p_dim']>generic else []), f"full-matrix torsion mismatch: p={p}, L={L}, torsion={torsion}, expected={expected}")
            results.append({'p':p,'L':L,'rank_Q':len(inv),'nonunit_invariants':torsion,'forest':expected})
        print(f'exact full-matrix SNF p={p} complete',flush=True)
    summary={'kind':'FINITE-EXACT over every coefficient field','count':len(results),'physical_primes':[5,7,11,13,17,19,23,29,31,37,43],'L':'every1<=L<=p','errors':0}
    (args.output_dir/'full_smith_results.json').write_text(json.dumps({'summary':summary,'cases':results},indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()

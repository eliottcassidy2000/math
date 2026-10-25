#!/usr/bin/env python3
"""Independent physical-point incidence; no producer imports or reduced rows."""
import argparse, json, pathlib
import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
OUT=pathlib.Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def points(p,L):
    return sorted({((-a-3*b)%p,(j-9*b)%p) for a in (0,1) for b in (0,1) for j in range(L)})

def incidence(p,L):
    pts=points(p,L)
    require(len(pts)==4*L, f"physical point collision: p={p}, L={L}, actual={len(pts)}, expected={4*L}")
    rows=[]
    # A line label is the determinant of its direction and the point.
    # Group actual distinct points, independent of their a,b,j names.
    for dx,dy in ((1,0),(0,1),(3%p,9%p)):
        groups={}
        for k,(x,y) in enumerate(pts):
            groups.setdefault((dx*y-dy*x)%p,[]).append(k)
        for group in groups.values():
            row=[0]*len(pts)
            for k in group: row[k]=1
            rows.append(row)
    return rows

def rank_mod(rows,q):
    a=np.asarray(rows,dtype=np.int64).copy()%q
    nr,nc=a.shape
    r=0
    for c in range(nc):
        nz=np.flatnonzero(a[r:,c])
        if not len(nz): continue
        s=r+int(nz[0]); a[[s,r]]=a[[r,s]]
        a[r]=(a[r]*pow(int(a[r,c]),-1,q))%q
        if r+1<nr:
            a[r+1:]=(a[r+1:]-a[r+1:,c,None]*a[r][None,:])%q
        r+=1
        if r==nr: break
    return r

def predicted(p,L,q):
    return L-1+max(0,L-10)+max(0,L-p+6)+max(0,L-p+3)+(q==p and L>=p-6)

def reduction(p,L):
    I=set(range(L)); H=sorted(set(range(p))-I)
    E=sorted(j for j in I if (j+9)%p in I)
    # Integral reduced equations on r_E plus c. Explicitly retain p*c=0.
    rows=[]
    for j in H:
        rows.append([sum((j-d)%p==k for d in (3,6,9)) for k in E]+[-1])
    rows.extend([[1]*len(E)+[0],[0]*len(E)+[p]])
    return E,rows

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=pathlib.Path, default=OUT)
    parser.add_argument("--negative-control", action="store_true", help="Inject one deliberately false reduced-dimension claim; the audit must fail.")
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tests=[]; small={}; snfdata={}; errors=[]
    ps=[5,7,11,13,17,19,23,29,31,37,43]
    qs=[2,3,5,7,11,13,17,19,23,29,31,37,43,47,101]
    for p in ps:
        small[str(p)]={}
        for L in range(1,p+1):
            rows=incidence(p,L)
            E,b=reduction(p,L)
            d=smith_normal_form(Matrix(b),domain=ZZ)
            inv=[abs(int(d[i,i])) for i in range(min(d.shape)) if d[i,i]]
            snfdata[f'{p},{L}']={'E':len(E),'invariants':inv}
            got={}
            for q in qs:
                dim=4*L-rank_mod(rows,q)
                red=L-1+len(E)+1-sum(v%q!=0 for v in inv)
                if args.negative_control and (p,L,q)==(5,1,2):
                    red+=1  # Deliberate hostile: a false dimension must be rejected.
                got[str(q)]=dim
                if dim!=red: errors.append({'kind':'reduction','p':p,'L':L,'q':q,'direct':dim,'reduced':red})
                if p>=17 and dim!=predicted(p,L,q):
                    errors.append({'kind':'formula','p':p,'L':L,'q':q,'direct':dim,'pred':predicted(p,L,q)})
                tests.append((p,L,q,dim))
            if p<17: small[str(p)][str(L)]=got
        print(f'p={p} complete, tests={len(tests)}, errors={len(errors)}',flush=True)
    summary={'kind':'FINITE-EXACT','physical_primes':ps,'coefficient_characteristics':qs,'L':'every 1<=L<=p','test_count':len(tests),'errors':errors}
    (args.output_dir/'direct_results.json').write_text(json.dumps({'summary':summary,'small':{p:v for p,v in small.items() if int(p)<17},'snf':snfdata},indent=2)+'\n')
    (args.output_dir/'direct_cases.json').write_text(json.dumps(tests)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    require(not errors, f"physical-incidence audit rejected {len(errors)} false claims: {errors}")
if __name__=='__main__': main()

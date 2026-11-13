"""Literal arbitrary-module controls for the protected 462-chain.

No inherited census or generated binary is required.  Run from repo root:
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/module_lift.py
"""
from functools import lru_cache
from pathlib import Path
from math import gcd
import json


RECIPES = (
    ((2,4,1), (0,), (1,0,1,1), 462, -5),
    ((2,36,57,51,35,1), (), (1,0,1), 1386, -13),
    ((2,4,105,115,103,67,3), (), (0,1,1), 4158, -37),
)


def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))


def induced(T, keep):
    return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(keep)) for u in keep)


def opposite(T):
    return tuple(((1<<len(T))-1)^(1<<i)^row for i,row in enumerate(T))


def tables(T):
    n=len(T); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1: continue
            rest=mask^(1<<v)
            if rest:
                dp[mask][v]=sum(dp[rest][u] for u in range(n) if (rest>>u)&1 and (T[u]>>v)&1)
    return [1]+[sum(row) for row in dp[1:]],dp


def profile(T, r):
    H,ends=tables(T)
    start=tables(opposite(T))[1][-1][r]
    q=H[((1<<len(T))-1)^(1<<r)]
    F=sum(ends[mask][r]*H[((1<<len(T))-1)^mask] for mask in range(1<<len(T)))
    return (H[-1],start,q,F)


def birth(T,r):
    n=len(T)
    return tuple(row|((1<<n) if i==r else 0) for i,row in enumerate(T))+((((1<<n)-1)^(1<<r)),)


def build(recipe,M):
    C,left,right,A,B=recipe
    m=len(M); n0=len(C); tail=len(left)+m+len(right); x=n0+tail
    rows=[0]*(x+1)
    for i,row in enumerate(C):
        rows[i]=row|(((1<<tail)-1)<<n0)|((1<<x) if i else 0)
    rows[x]=1
    module=tuple(range(n0+len(left),n0+len(left)+m))
    bits=left+(1,)*m+right
    for i,bit in enumerate(bits):
        v=n0+i
        rows[v]=sum(1<<(n0+j) for j in range(i+1,tail))
        if v in module:
            rows[v]&=~sum(1<<u for u in module)
            rows[v]|=M[v-module[0]]<<module[0]
        if bit:rows[v]|=1<<x
        else:rows[x]|=1<<v
    T=birth(tuple(rows),x)
    return T,module,x,len(T)-1


def run_coeff(T,m,visits):
    """Count reduced exterior words, allowing the marked placeholder r visits."""
    outside=tuple(i for i in range(len(T)) if i!=m)
    full=(1<<len(outside))-1
    @lru_cache(None)
    def rec(mask,last,k):
        if mask==full and k==visits:return 1
        out=0
        if last!=m and k<visits and (last==-1 or (T[last]>>m)&1):
            out+=rec(mask,m,k+1)
        for i,v in enumerate(outside):
            if not (mask>>i)&1 and (last==-1 or (T[last]>>v)&1):
                out+=rec(mask|(1<<i),v,k)
        return out
    return rec(0,-1,0)


def run():
    coefficients=[]
    for recipe in RECIPES:
        T,M,x,y=build(recipe,(0,))
        hcoef=[run_coeff(T,M[0],j) for j in (1,2,3)]
        parent=T[:-1]
        parent=tuple(row&((1<<len(parent))-1) for row in parent)
        qcoef=[run_coeff(parent,M[0],j) for j in (1,2,3)]
        A,B=recipe[-2:]
        assert hcoef==[2*A+B,A,0]
        coefficients.append(dict(A=A,B=B,H=hcoef,deleted_H=qcoef))
    records=[]
    for name,M in [('K1',(0,)),('TT2',tt(2)),('TT3',tt(3)),('C3',(2,4,1))]:
        htab,_=tables(M); h=htab[-1]; full=(1<<len(M))-1
        c=sum(htab[s]*htab[full^s] for s in range(1,full))
        stages=[]
        for j,recipe in enumerate(RECIPES):
            T,block,x,y=build(recipe,M)
            z=profile(T,y); A,B=recipe[-2:]
            q1,q2,_=coefficients[j]['deleted_H']
            expected=((2*A+B)*h+A*c,(q1-1)*h+q2*c,q1*h+q2*c,(2*A+B+1)*h+A*c)
            assert z==expected,(name,j,z,expected)
            assert [i for i,row in enumerate(T) if (row>>y)&1]==[x]
            nr=birth(birth(T,y),len(T))
            # NR retains the first newborn; the second birth has index len(T)+1.
            nr_z=profile(nr,len(T)) if len(nr)<=15 else None
            if nr_z is not None:
                assert nr_z==(3*z[0]+2*h,z[0],z[0],3*z[0]+2*h)
            stages.append(dict(rows=T,module=block,predecessor=x,root=y,state=z,NR_state=nr_z))
        for j in (0,1):
            assert stages[j+1]['state'][0]==3*stages[j]['state'][0]+2*h
            t1,t2=stages[j],stages[j+1]
            assert induced(t1['rows'],tuple(t1['module'])+(t1['predecessor'],t1['root']))==induced(t2['rows'],tuple(t2['module'])+(t2['predecessor'],t2['root']))
        records.append(dict(module=name,module_rows=M,h=h,C2=c,stages=stages))
    primitive=[]
    for name,M,h0 in [('H35',(8,5,9,18,71,31,47),35),('H49',(4,9,18,21,67,31,47),49)]:
        H,_=tables(M); full=(1<<len(M))-1
        h=H[-1]; c=sum(H[s]*H[full^s] for s in range(1,full))
        assert h==h0
        stages=[]
        for j,recipe in enumerate(RECIPES):
            T,block,x,y=build(recipe,M)
            A,B=recipe[-2:];q1,q2,_=coefficients[j]['deleted_H']
            z=((2*A+B)*h+A*c,(q1-1)*h+q2*c,q1*h+q2*c,(2*A+B+1)*h+A*c)
            stages.append(dict(rows=T,module=block,root=y,predecessor=x,state=z,root_content=gcd(*z)))
        assert [z['root_content'] for z in stages]==[7,1,7]
        primitive.append(dict(module=name,module_rows=M,h=h,C2=c,P2=c+2*h,stages=stages,
                              universe='Module invariants independently counted; lifted states use proved run coefficients, not a new large-graph census'))
    result=dict(status='PASS',scope='Three literal inherited recipes; independent reduced-run coefficients through3; four actual modules K1,TT2,TT3,C3; complete rooted subset-DP profiles; direct NR up to order15; two small module certificates for the all-parameter primitive repair',coefficients=coefficients,controls=records,primitive_repairs=primitive)
    Path(__file__).with_name('module_lift_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(coefficients=coefficients,controls=[(r['module'],[s['state'] for s in r['stages']]) for r in records]),indent=2))


if __name__=='__main__':run()

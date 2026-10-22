#!/usr/bin/env python3
"""Exact colored forest-law audits. Python standard library only.

Reproduce: python3 session31-blue-forests/check_blue_forests.py
Universe: tournaments on 0,...,n-1 with every i->i+1 fixed.
Blue means converse followed by label reversal fixes the presentation.
The forest convention is outdegree <=1, directed toward a root.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import json
from probe import edge_orbits, get_lap, forest_coeff

HERE=Path(__file__).resolve().parent

def all_mean(n):
    p=[Q(0)]*n
    for l in range((n-1)//2+1):
        for j in range(n-2*l):
            p[2*l+j]+=Q(comb(n-1-l,l)*comb(n-1-2*l,j)*n**j,4**l*2**j)
    return p

def mask_reflection(mask, edges, n):
    pos={e:k for k,e in enumerate(edges)}
    return sum(((mask>>k)&1)<<pos[(n-1-j,n-1-i)] for k,(i,j) in enumerate(edges))

def energy(L):
    n=len(L)
    return sum((2*L[i][i]-(n-1))**2 for i in range(n))

def sum_law(n,orbs,full):
    totals=[0]*n; esum=0
    for mask in range(1<<len(orbs)):
        L=get_lap(n,orbs,mask); e=energy(L);esum+=e
        if full:
            f=forest_coeff(L)
            assert f[0]==1 and f[1]==n*(n-1)//2
            if n>=3:
                assert 8*f[2]==4*(n*(n-1)//2)**2-n*(n-1)**2-e
            totals=[x+y for x,y in zip(totals,f)]
    return esum,totals

def constraint_sums(n):
    """Independent enumeration of all possible rooted directed forests."""
    edges,orbs=edge_orbits(n); D=len(edges);o=len(orbs)
    oid={e:k for k,orb in enumerate(orbs) for e in orb}
    choices=[[-1]+[j for j in range(n) if j!=i] for i in range(n)]
    blue=[0]*n;all_=[0]*n;accepted=0;conflicts=0
    for parent in product(*choices):
        # Reject any cycle, and reject backwards fixed path arcs.
        if any(p==i-1 for i,p in enumerate(parent) if p>=0):continue
        cyclic=False
        for i in range(n):
            seen=set();j=i
            while j>=0:
                if j in seen:cyclic=True;break
                seen.add(j);j=parent[j]
            if cyclic:break
        if cyclic:continue
        accepted+=1;k=sum(p>=0 for p in parent);req={};r=0;compatible=True
        for i,j in enumerate(parent):
            if j<0 or j==i+1:continue
            r+=1;e=tuple(sorted((i,j))); b=int(i<j);h=oid[e]
            if h in req and req[h]!=b:compatible=False
            req[h]=b
        all_[k]+=2**(D-r)
        if compatible:blue[k]+=2**(o-len(req))
        else:conflicts+=1
    return all_,blue,accepted,conflicts

def bareiss(M):
    """Independent integer determinant control."""
    a=[r[:] for r in M];n=len(a);sgn=1;prev=1
    for k in range(n-1):
        if not a[k][k]:
            pivot=next((i for i in range(k+1,n) if a[i][k]),None)
            if pivot is None:return 0
            a[k],a[pivot]=a[pivot],a[k];sgn=-sgn
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=a[i][j]*p-a[i][k]*a[k][j]
                assert z%prev==0;a[i][j]=z//prev
        for i in range(k+1,n):a[i][k]=0
        prev=p
    return sgn*a[-1][-1]

def main():
    report={'status':'PASS','universe':'Fixed directed path 0->1->...->n-1; all other arcs fair; blue is reflection-fixed','moments':[],'full_polynomials':[],'constraint_enumeration':[],'hostile_controls':[]}
    for n in range(2,10):
        edges,orbs=edge_orbits(n);D=len(edges);o=len(orbs);f=sum(len(z)==1 for z in orbs)
        A=2**D;B=2**o;mu=n*n-3*n+4;delta=n-3 if n%2 else 0
        assert f==(n-1)//2 and o==(D+f)//2==(n-1)**2//4
        e,s=sum_law(n,orbs,n<=8)
        assert e==B*(mu-delta)
        row={'n':n,'D':D,'fixed_edges':f,'orbits':o,'all_tilings':A,'blue_tilings':B,'blue_energy_sum':e,'blue_energy_mean':str(Q(e,B)),'all_energy_mean_formula':mu}
        if n<=7:
            ae,ap=sum_law(n,[[x] for x in edges],n<=6)
            assert ae==A*mu
            row['all_energy_sum_enumerated']=ae
            if n>=4:
                assert Q(ae-e,A-B)==mu+Q(B*delta,A-B)
                row['black_energy_mean']=str(Q(ae-e,A-B))
            if n<=6:assert ap==[int(A*z) for z in all_mean(n)]
        if n>=3:
            row['blue_lines']=B//2
            row['literal_black_lines']=(A-B)//2
            row['fundamental_black_edges']=(A-B)//4
        report['moments'].append(row)
        if n<=8:
            b=[Q(x,B) for x in s];p=all_mean(n);diff=[x-y for x,y in zip(b,p)]
            if n>=3:assert diff[2]==Q(delta,8)
            if n>=4:assert diff[3]==Q(delta*n*(n-3),16)
            report['full_polynomials'].append({'n':n,'blue_total_coefficients':s,'blue_mean_coefficients':list(map(str,b)),'all_mean_coefficients':list(map(str,p)),'blue_minus_all_coefficients':list(map(str,diff))})
            for mask in sorted({0,(1<<o)-1,(1<<o)//3}):
                L=get_lap(n,orbs,mask);c=forest_coeff(L)
                for t in (-2,1,3):
                    M=[[int(i==j)+t*L[i][j] for j in range(n)] for i in range(n)]
                    assert bareiss(M)==sum(x*t**k for k,x in enumerate(c))
                # For a blue tournament, K=2L-nI+J anticommutes with reversal.
                K=[[2*L[i][j]-n*int(i==j)+1 for j in range(n)] for i in range(n)]
                assert all(K[n-1-i][n-1-j]==-K[i][j] for i in range(n) for j in range(n))
        if n<=6:
            ac,bc,count,conflicts=constraint_sums(n)
            assert ac==[int(A*z) for z in all_mean(n)] and bc==s
            report['constraint_enumeration'].append({'n':n,'all_coefficients':ac,'blue_coefficients':bc,'path_compatible_forests':count,'rejected_reflection_conflicts':conflicts})
        if 3<=n<=7:
            full=(1<<D)-1;seen=set();cnt={2:0,4:0}
            for m in range(A):
                if m in seen:continue
                r=mask_reflection(m,edges,n);orb={m,m^full,r,r^full}
                assert len(orb) in (2,4)
                assert (r==m)==(len(orb)==2)
                seen.update(orb);cnt[len(orb)]+=1
            assert cnt[2]==B//2 and cnt[4]==(A-B)//4
            row['klein_orbits_enumerated']=cnt
        print('PASS n=',n,flush=True)
    polys={r['n']:r for r in report['full_polynomials']}
    assert polys[5]['blue_minus_all_coefficients']==['0','0','1/4','5/4','27/16']
    assert polys[6]['blue_minus_all_coefficients']==['0','0','0','0','5/16','15/16']
    report['hostile_controls']=[
      'n=5: replacing correlated blue arcs by their identical marginal means fails already at t^2, difference 1/4.',
      'n=6: equal all/blue score-energy means do not imply equal forest laws; difference (5/16)t^4(1+3t).',
      'n=3: all presentations are blue, so conditional black expectations are undefined.',
      'n=2: the complement involution is identity, so blue-line count B/2 must not be used.'
    ]
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All checks passed. Wrote',HERE/'verification.json')

if __name__=='__main__':main()

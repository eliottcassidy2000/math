#!/usr/bin/env python3
"""Separate inverse-word and coefficient replay for the uniform remainder.

Does not import the positive-table or integer-row producers.
"""
import argparse, hashlib, json, math
from fractions import Fraction as Q
from collections import defaultdict
from pathlib import Path
import numpy as np

def check(ok,m):
    if not ok:raise RuntimeError(m)
def chi(n):return (0,1,-1)[n%3]
def w(n):return (-1,-1,2)[n%3]
def seed(n):
    while n%8==5:n=(n-1)//4
    return n
def initial(N,s):
    ans=0
    while s<=4*N:
        ans+=-3*chi(s) if s<=2*N else w(s);s=4*s+1
    return ans
def words(d,budget):
    if not d:
        yield ()
    else:
        for A in range(1,budget-d+2):
            for tail in words(d-1,budget-A):yield (A,)+tail

def inverse_word(k,word):
    alpha=Q(1);beta=Q(0)
    for A in reversed(word):
        a=1 if A%2 else 2;c=(3+2**a-2**A)//3
        alpha*=Q(2**A,3);beta=(2**A*beta-c)/3
    return 4**k*alpha,4**k*beta+Q(4**k-1,3)

def word_audit(N,d):
    m=(4*N-1).bit_length();K=m+2*d;q=3**d
    direct={s:initial(N,s) for s in range(1,4*N+1,2) if s%8!=5}
    for unused in range(d):
        out=defaultdict(int)
        for s,c in direct.items():
            a=2 if s%8==1 else 1;out[seed((3*s+1)//2**a)]+=(-1)**a*c
        direct=out
    inverse=defaultdict(int);events=defaultdict(lambda:defaultdict(Q));count=0;largest=Q(0)
    top=math.ceil(Q(3,2)**d*(4*N+1))
    for k in range((K-d)//2+1):
        for word in words(d,K-2*k):
            count+=1;alpha,beta=inverse_word(k,word);largest=max(largest,abs(beta/alpha))
            check(alpha.denominator<=q and (q*alpha).denominator==1 and (q*beta).denominator==1,'affine denominator')
            r0=(-int(q*beta)*pow(int(q*alpha),-1,q))%q if q>1 else 0
            sign=(-1)**sum(word)
            for r in (r0,r0+q,r0+2*q):
                s=alpha*r+beta;check(s.denominator==1,'word progression integral')
                outer=sign*w(int(s));inner=sign*(-3*chi(int(s))-w(int(s)))
                events[r][Q(4*N)/alpha]+=outer;events[r][Q(2*N)/alpha]+=inner
                for t in range(r,top+1,3*q):
                    if t<1 or t%2==0 or t%8==5:continue
                    actual=alpha*t+beta
                    if actual<=4*N:inverse[t]+=outer
                    if actual<=2*N:inverse[t]+=inner
    check(count==sum(math.comb(K-2*k,d) for k in range((K-d)//2+1)),'exact binomial word count')
    check(largest<=Q(4,3)*Q(3,2)**d-1,'normalized carry bound for every word')
    check({k:v for k,v in direct.items() if v}=={k:v for k,v in inverse.items() if v},'complete actual inverse and forward rows agree')
    integral=Q(0)
    for r,ev in events.items():
        value=sum(ev.values());before=Q(0)
        for end,delta in sorted(ev.items()):
            integral+=abs(value)*(end-before)/Q(8*q);value-=delta;before=end
        check(value==0,'step function vanishes after the final threshold')
    actual=sum(abs(v) for v in direct.values())
    R=Q(3,2)**d
    outer=sum(abs(v) for t,v in direct.items() if t>R*(2*N+1)-1)
    predicted=0
    lower=Q(2*N+1,2**d);upper=Q(4*N+1,2**d)
    for h in range(math.floor(lower)+1,math.floor(upper)+1):
        t=3**d*h-1
        if t%2 and t%8!=5:predicted+=abs(w(2**d*h-1))
    check(outer==predicted,'exact isolated all-one outer packet')
    check(abs(Q(outer)-Q(N,2**d))<=12,'uniform lower bound for actual mass')
    error=abs(Q(actual)-integral);bound=36*(Q(4,3)*Q(3,2)**d+1)*count
    check(error<=bound,'literal affine plus quadrature bound')
    return dict(N=N,depth=d,K=K,words=count,max_offset=str(largest),actual_l1=actual,isolated_outer_mass=outer,truncated_integral=str(integral),actual_to_truncated_error=str(error))

def band_audit(receipt):
    B=receipt['rows'][0]['bands'];old=None;checks=0
    for d in range(5):
        mod=3**(d+1);rows=[]
        for b in range(B):
            row=[]
            for r in range(mod):
                value=0
                if d==0:
                    x=Q(3,2**b);rr=r
                    while x<=4:
                        value+=-3*chi(rr) if x<=2 else w(rr);x*=4;rr+=1
                else:
                    rho=r
                    for j in range(b//2+1):
                        if rho%3:
                            a=2 if rho%3==1 else 1;bb=b-(a-1)-2*j
                            if bb>=0:value+=(-1)**a*old[bb][((2**a*rho-1)//3)%(mod//3)]
                        rho=(4*rho+1)%mod
                row.append(value);checks+=1
            rows.append(row)
        check([sum(abs(v) for v in row) for row in rows]==receipt['rows'][d]['absolute_band_sums'],'inverse-table absolute bands')
        demod=np.array([[(-1)**(b+d)*v for v in row] for b,row in enumerate(rows)],dtype=np.int64)
        check(hashlib.sha256(demod.tobytes()).hexdigest()==receipt['rows'][d]['coefficient_sha256'],'all signed coefficients independently match')
        old=rows
    old_receipt=json.loads(Path('04-computation/experiments/joint_logistic_current_20261017f/iterated_tariffs.json').read_text())
    for d in range(4):
        exact=Q(old_receipt['depths'][d]['tariff']);row=receipt['rows'][d]
        check(Q(row['lower'])<=exact<=Q(row['upper']),'inherited exact tariff lies in new rational enclosure')
    return dict(independent_coefficient_checks=checks,maximum_depth=4,inherited_tariffs_recovered=4)

def colour_audit():
    r=json.loads(Path('04-computation/experiments/root_positive_transport_20261017g/golden_shared.json').read_text())
    n=r['vertices'];edges=r['edges'];patterns=r['patterns'];witnesses=r['witnesses'];checks=0
    check(len(set(map(tuple,edges)))==len(edges),'edge uniqueness')
    for a,b in edges:check(0<=a<b<n,'simple labelled graph')
    for pat,col in zip(patterns,witnesses):
        check(len(col)==n and list(map(int,col[:8]))==pat,'literal port colours')
        check(all(c in '012345' for c in col),'six colours')
        for a,b in edges:check(col[a]!=col[b],'source edge in boundary witness');checks+=1
    check(len(patterns)==len(witnesses)==570,'complete old boundary list')
    check(all(r['four_colour'][a]!=r['four_colour'][b] for a,b in edges),'four-colour hostile replay')
    return dict(vertices=n,edges=len(edges),mixed_edges=len(r['mixed_edges']),boundary_patterns=len(patterns),colour_edge_checks=checks)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    rec=json.loads(Path('04-computation/experiments/root_positive_transport_20261017g/positive.json').read_text())
    ans=dict(status='FINITE-EXACT independent formulas and witness replay; self-audit, not an independent referee',bands=band_audit(rec),words=[word_audit(N,d) for N in (1,3,16,63,257) for d in range(5)],colour=colour_audit())
    Path(a.output).write_text(json.dumps(ans,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in ans.items() if k!='words'},sort_keys=True));print('word cases',len(ans['words']),'words',sum(v['words'] for v in ans['words']))

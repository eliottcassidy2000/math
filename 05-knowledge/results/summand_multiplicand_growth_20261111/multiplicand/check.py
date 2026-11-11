"""Exact strict-multiplicand birth and inverse controls; no external libraries."""
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path
import hashlib,json,subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
LIMIT=512
SOURCES=[
 '01-canon/theorems/THM-362-natural-operation-graph-shadows.md',
 '01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md',
 '05-knowledge/results/collatz_mod6_20260922_lead_reflection2.md',
 '05-knowledge/results/session259_arithmetic_graph_growth_20261111/multiplicand/RESEARCH.md',
]


def need(x,label):
    if not x:raise RuntimeError(label)


def factor(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out


FAC=[factor(n) for n in range(LIMIT+1)]
DIV=[[] for _ in range(LIMIT+1)]
for d in range(1,LIMIT+1):
    for n in range(d,LIMIT+1,d):DIV[n].append(d)


def mu(n):
    exponents=FAC[n].values()
    return 0 if any(e>1 for e in exponents) else (-1)**len(FAC[n])


def eta(j,n):
    values=FAC[n].values()
    if any(e>=2**(j+1) for e in values):return 0
    return (-1)**sum(e.bit_count() for e in values)


def kernel(d,n):
    power=d;j=0;answer=0;terms=[]
    while power<=n:
        if n%power==0:
            term=eta(j,n//power)
            if term:terms.append((j,term));answer+=term
        power*=power;j+=1
    return answer,terms


def interval_scales(d,n):
    a=FAC[d];e=FAC[n];primes=a.keys()|e.keys()
    lower=max((Fraction(e.get(p,0),a.get(p,0)+2) for p in primes),default=Fraction(0))
    upper=min(Fraction(e.get(p,0),ap) for p,ap in a.items())
    out=[];power=1;j=0
    while power<=upper:
        if lower<power:out.append(j)
        power*=2;j+=1
    return out


def main():
    kernels=[[0]*(LIMIT+1) for _ in range(LIMIT+1)]
    edges=0;births=[];value_counts=Counter();nonzero_scale_counts=Counter()
    for n in range(2,LIMIT+1):
        square=isqrt(n)**2==n
        pairs=[(a,n//a) for a in range(2,isqrt(n)+1) if n%a==0 and a<n//a]
        incoming=sorted({v for pair in pairs for v in pair})
        need(len(incoming)==len(DIV[n])-2-square,'strict birth indegree')
        need(len(incoming)==2*len(pairs),'paired-incidence parity')
        edges+=len(incoming)
        need(edges==sum(len(DIV[x]) for x in range(1,n+1))-2*n-isqrt(n)+2,'cumulative edge count')
        if n<=36:births.append(dict(n=n,prime=len(DIV[n])==2,square=square,factor_packets=pairs,incoming=incoming))
        for d in DIV[n]:
            if d==1:continue
            k,terms=kernel(d,n);kernels[n][d]=k
            need([j for j,t in terms]==interval_scales(d,n),'prime-exponent interval criterion')
            need(len(terms)<=2,'at most two scales')
            if len(terms)==2:
                need(all(a==1 for a in FAC[d].values()),'two scales require squarefree source')
                need(terms[1][0]==terms[0][0]+1,'adjacent scales')
                need(terms[0][1]==(-1)**len(FAC[d])*terms[1][1],'two-scale sign law')
            need(abs(k)<=2,'bounded inverse entries')
            if n==d*d:need(k==1+mu(d),'square diagonal formula')
            value_counts[k]+=1;nonzero_scale_counts[len(terms)]+=1
    # Independent sparse multiplication of the literal incidence matrix by K.
    basis_entries=0
    for n in range(2,LIMIT+1):
        row=[d for d in DIV[n] if d>1 and d*d!=n]
        for b in range(2,n+1):
            value=sum(kernels[d][b] for d in row)
            need(value==int(n==b),'literal matrix times explicit inverse')
            basis_entries+=1
    f=[0]+[0]+[((n*37)%29)-14 for n in range(2,LIMIT+1)]
    g=[0]*(LIMIT+1);online=[0]*(LIMIT+1);closed=[0]*(LIMIT+1)
    for n in range(2,LIMIT+1):
        g[n]=sum(f[d] for d in DIV[n] if d>1 and d*d!=n)
        online[n]=g[n]-sum(online[d] for d in DIV[n] if 1<d<n and d*d!=n)
        closed[n]=sum(kernels[n][d]*g[d] for d in DIV[n] if d>1)
    need(f==online==closed,'signed-weight reconstruction')
    pins={p:hashlib.sha256(subprocess.check_output(['git','show','HEAD:'+p],cwd=ROOT)).hexdigest() for p in SOURCES}
    receipt=dict(status='PASS; finite exact controls for all-order proofs in PROOF.md',
                 limit=LIMIT,universe='strict nonunit graph on2..512; all lower-triangular basis equations; exact integer arithmetic',
                 matrix_entries_checked=basis_entries,source_sha256=pins,births_through36=births,
                 kernel_values_on_divisor_pairs=dict(sorted(value_counts.items())),nonzero_scale_counts=dict(sorted(nonzero_scale_counts.items())),
                 examples={f'{d}->{n}':dict(value=kernel(d,n)[0],terms=kernel(d,n)[1])
                           for d,n in ((2,4),(3,6),(4,16),(6,36),(6,180),(30,900)) if n<=LIMIT},
                 tests=['factor-packet birth law and cumulative edge formula',
                        'closed kernel equals prime-exponent interval/sign classifier',
                        'literal strict incidence matrix times explicit inverse equals identity',
                        'online and closed inverse recover arbitrary signed input',
                        'entrywise bound2 and square-diagonal value1+mu'])
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('births_through36','source_sha256')},indent=2))


if __name__=='__main__':main()

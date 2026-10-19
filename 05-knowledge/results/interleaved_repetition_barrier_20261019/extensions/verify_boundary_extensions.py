#!/usr/bin/env python3
"""Exact controls for separately staged structural extensions; no dependencies."""
from fractions import Fraction
from math import gcd
import hashlib
import json
checks=0

def need(v,label):
    global checks
    checks+=1
    if not v:raise RuntimeError(label)

def branch(bits):
    a,b,M=1,0,1
    for bit in bits:
        if bit=='1':a,b=3*a,3*b+M
        elif bit!='0':raise ValueError(bit)
        M*=2
    return a,b,M

def compose(u,v):
    a,b,M=u;c,d,Q=v
    return c*a,c*b+d*M,Q*M

def digest(n):
    return hashlib.sha256(n.to_bytes((n.bit_length()+7)//8,'big')).hexdigest()

K='110110110110110101011011101010'
bits={'V':4*K+'10','F':'10'}
data={z:branch(w) for z,w in bits.items()}
aV,bV,mV=data['V'];v=Fraction(aV,mV)
need(aV==3**77 and mV==2**122,'V exact prime exponents')
need(1<v and v**9<Fraction(4,3)<v**10,'nine/ten run control')
need(Fraction(2**48,450)>53*10**10+174,'arbitrary-position logarithmic loss absorbed')
need(Fraction(123*5,50)+Fraction(123,61)+498<122*5+1,'complexity inequality starts at m=5')
need(Fraction(123,50)<122,'complexity inequality persists')
endpoint_cases=0
for n in range(1,9):
    for mask in range(1<<n):
        w=format(mask,'0%db'%n)
        a,b,M=branch(w)
        if a>=M:continue
        Delta=M-a
        need(gcd(M,Delta)==1,'endpoint modulus coprime')
        delta=b*pow(M,-1,Delta)%Delta if Delta>1 else 0
        hits=[]
        for d in range(b//M+1):
            eq=(b>M*d and (b-M*d)%Delta==0)
            residue=(d>=delta and (d-delta)%Delta==0 and M*d<b)
            need(eq==residue,'carry-residue iff diagonal source')
            if eq:
                x=(b-M*d)//Delta;z=x
                for bit in w:
                    need(z%2==int(bit),'diagonal equality authenticates parity')
                    z=(3*z+1)//2 if z%2 else z//2
                need(z==x+d,'diagonal endpoint');hits.append((d,x))
            endpoint_cases+=1
        need(bool(hits)==(M*delta<b),'minimal residue iff non-descent exists')
records=[]
for m in (1,10,901,1024):
    a,b,M=1,0,1;word=[];runs=[];run=0;conditional=0;mechanical=[]
    for _ in range(m):
        a,b,M=compose((a,b,M),data['V']);word.append('V');run+=1
        need(a>=M,'expansion before crossing')
        reset=3*a>=4*M;mechanical.append('1' if reset else '0')
        if reset:
            a,b,M=compose((a,b,M),data['F']);word.append('F')
            conditional+=1;runs.append(run);run=0
            need(a>=M,'conditional reset above crossing')
        need(M<=a and 3*a<4*M,'greedy interval invariant')
    if run:runs.append(run)
    a,b,M=compose((a,b,M),data['F']);word.append('F')
    need(a<M,'final strict first crossing')
    need(max(runs)<=10,'all expansion runs bounded')
    need(all(t in (9,10) for t in runs[:-1]),'complete runs nine or ten')
    literal=''.join(bits[z] for z in word);N=len(literal);L=122*m
    need(branch(literal)==(a,b,M),'independent literal coefficients')
    need(6*b<N*M and 6*b<N*a,'strict carry bounds')
    if m>900:need(100*1098<L,'no one-percent pure-expansion power')
    if m>915:need(100*1098<120*m,'rewrite-robust repetition exclusion')
    Delta=M-a;delta=b*pow(M,-1,Delta)%Delta if Delta>1 else 0
    need(M*delta>=b,'no nonnegative endpoint displacement in finite control')
    x=-b*pow(a,-1,M)%M
    need(x>0 and Delta*x>b,'least positive cylinder source descends')
    if m<=10:
        z=x
        for bit in literal:
            need(z%2==int(bit),'literal least-source guard')
            z=(3*z+1)//2 if z%2 else z//2
        need(M*z==a*x+b and z<x,'literal endpoint descent')
    r=(L+99)//100;P=len({literal[i:i+r] for i in range(N-r+1)})
    bound=124*((r+121)//122+3)+2
    need(P<=bound,'literal morphism complexity bound')
    if m>=5:need(P<N-r+1,'factor pigeonhole at one-percent scale')
    mech=''.join(mechanical)
    for k in (1,2,5,10,20):
        if k<=m:need(len({mech[i:i+k] for i in range(m-k+1)})<=k+1,'rotation factor count')
    records.append({'m':m,'N':N,'expansion_bits':L,'explicit_resets':conditional+1,
        'max_expansion_run':max(runs),'r':r,'literal_factor_count':P,
        'complexity_bound':bound,'source_bits':x.bit_length(),'source_sha256':digest(x),
        'delta_bits':delta.bit_length(),'delta_sha256':digest(delta),
        'literal_sha256':hashlib.sha256(literal.encode()).hexdigest(),'literal_replay':m<=10})
trajectory=[27]
for _ in range(180):
    x=trajectory[-1];trajectory.append((3*x+1)//2 if x%2 else x//2)
parity=''.join(str(x%2) for x in trajectory[:-1]);matched=0
for r in (1,2,5,10,20):
    seen={}
    for i in range(len(parity)-r+1):
        w=parity[i:i+r]
        for j in seen.get(w,[]):
            need((trajectory[i]-trajectory[j])%(1<<r)==0,'matched block binary divisibility')
            if trajectory[i]!=trajectory[j]:need(abs(trajectory[i]-trajectory[j])>=1<<r,'distinct-state binary gap')
            else:
                n=min(len(parity)-i,len(parity)-j)
                need(parity[i:i+n]==parity[j:j+n],'exact return preserves aligned continuation')
            matched+=1
        seen.setdefault(w,[]).append(i)
    P=len(seen)
    if P<len(parity)-r+1:need(len({parity[i:i+r] for i in range(P+1)})<=P,'first P+1 starts pigeonhole')
print(json.dumps({'status':'PASS exact controls; universal statements proved in notes',
    'checks':checks,'diagonal_tests':endpoint_cases,'matched_pairs':matched,
    'greedy_controls':records},indent=2,sort_keys=True))

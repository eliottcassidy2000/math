#!/usr/bin/env python3
"""Exact circuit-ratio realization; integer arithmetic and rational sign brackets."""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
from hashlib import sha256
import json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
GATES=0
def need(ok,why):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(why)
def sign(x):return (x>0)-(x<0)
def canonical(x):return json.dumps(x,separators=(',',':'),sort_keys=True).encode()+b'\n'
def rational(x):return str(x)
def evaluate(a,x):
    value=F(0)
    for v in reversed(a):value=value*x+v
    return value
def realize(C,lam=None):
    d=len(C)+2
    prefixes=[F(1)]
    for c in C:prefixes.append(prefixes[-1]*c)
    kappas=[F(comb(d,k-1)*comb(d,k+1),comb(d,k)**2) for k in range(1,d)]
    if lam is None:lam=9*max(k/p for k,p in zip(kappas,prefixes))
    R=[lam*p for p in prefixes];h=[F(1),F(1)]
    for k in range(1,d):h.append(h[k]**2/(R[k-1]*h[k-1]))
    a=[comb(d,k)*h[k] for k in range(d+1)]
    q=[a[k-1]/a[k] for k in range(1,d+1)]
    actual=[h[k]**2/(h[k-1]*h[k+1]) for k in range(1,d)]
    need(actual==R,'literal normalized Newton ratios')
    need([actual[k]/actual[k-1] for k in range(1,d-1)]==list(C),'exact prescribed circuit ratios')
    need(a[0]==1 and a[1]==d,'monic reversal and fixed root-parameter sum')
    need(all(a[k]**2>=9*a[k-1]*a[k+1] for k in range(1,d)),'conservative coefficient separation')
    need(all(q[k]>=9*q[k-1] for k in range(1,d)),'strictly ordered sample addresses')
    samples=[F(0)]+[3*v for v in q]
    values=[evaluate(a,-x) for x in samples]
    need([sign(v) for v in values]==[(-1)**k for k in range(d+1)],'literal alternating rational brackets: d simple negative roots')
    return dict(d=d,C=list(map(str,C)),R=list(map(str,R)),coefficients=list(map(str,a)),samples=list(map(str,samples)),sample_signs=list(map(sign,values)))
def main():
    counts=[];hashes=[];saved=[]
    for d in range(2,9):
        rows=[]
        for word in product((-1,0,1),repeat=d-2):
            C=[F(2)**s for s in word]
            row=realize(C,F(9*2**(d-2)))
            need([sign(F(c)-1) for c in row['C']]==list(word),'every ternary circuit word including exact ties')
            rows.append(row)
            if d in (3,5,8) and (len(set(word))==1 or tuple(word)==tuple((-1,0,1)[i%3] for i in range(d-2))):saved.append(row)
        need(len(rows)==3**(d-2),'complete positional ternary universe')
        counts.append([d,len(rows)]);hashes.append([d,sha256(canonical(rows)).hexdigest()])
    arbitrary=[]
    for C in ([F(1,1000),F(37),F(1),F(11,7)], [F(103,97),F(5,11)], [F(1)]*9):
        arbitrary.append(realize(C))
    # Root parameters (1/2,1,3/2) have fixed sum3 and square sum7/2.
    a=[F(1),F(3),F(11,4),F(3,4)]
    h=[a[k]/comb(3,k) for k in range(4)]
    R1=h[1]**2/h[2]
    need(R1==F(12,11),'fixed two moments fix the first Newton ratio')
    need(R1*F(1,2)<1,'same moment fibre cannot realize C2=1/2 by Newton inequality')
    # Anchored degree-five B has sum13 and sum squares59: R1=338/275.
    need(F(13,5)**2/F(55,10)==F(338,275),'inherited degree-five anchor first ratio')
    need(F(338,275)*F(1,2)<1,'universality does not preserve the anchored second moment')
    cert=dict(status='PROVED analytic realization; FINITE-EXACT controls',counts=counts,complete_bank_sha256=hashes,
      controls=saved,arbitrary_positive_targets=arbitrary,fixed_moment_hostile=dict(root_parameters=['1/2','1','3/2'],R1='12/11',forbidden_C2='1/2'),
      inherited_anchor=dict(sum=13,sum_squares=59,R1='338/275',forbidden_C2='1/2'),
      scope='Arbitrary positive circuit ratios with fixed degree and root sum; root second moment and original factorial/wall coefficients are not preserved.')
    here=Path(__file__).resolve();dest=here.parent
    if here.parent.name=='04-computation':dest=here.parent.parent/'05-knowledge/results'
    data=canonical(cert);(dest/(here.stem+'_certificate.json')).write_bytes(data)
    print('EXACT_CIRCUIT_SURJECTIVITY arbitrary positive target vector; root-parameter sum equals degree')
    print('TERNARY_UNIVERSE',counts,'TOTAL',sum(c for d,c in counts))
    print('ROOT_PROOF literal alternating rational samples under coefficient separation9; no numeric roots')
    print('HOSTILE fixed second moment forces R1; targetC2=1/2 fails Newton in two exact anchored examples')
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('PASS',GATES,'always-active exact gates; actual LF')
if __name__=='__main__':main()

"""Independent closed-product and Sturm audit of exact Newton circuit surjectivity."""
from fractions import Fraction as F
from itertools import product
from math import comb,lcm,gcd
from pathlib import Path
from hashlib import sha256
import argparse
import json
import sys
import sympy as S
sys.stdout.reconfigure(newline="\n")
GATES=0
def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(why)

parser=argparse.ArgumentParser()
here=Path(__file__).resolve().parent
default=here if here.name=="04-computation" else here.parent/"continuing8_20260906_root"
parser.add_argument("--producer-dir",type=Path,default=default)
args=parser.parse_args()
base=args.producer_dir
data=base.parent/"05-knowledge/results" if base.name=="04-computation" else base
stem="continuing8_20260906_newton_universality"
need(sha256((base/(stem+".py")).read_bytes()).hexdigest()=="1fc809f98a154e75dff50bf20c43e8a68418a438a0af9ec08d19377f10e13700","frozen producer source")
saved=json.loads((data/(stem+"_certificate.json")).read_bytes())
need(sha256((data/(stem+"_certificate.json")).read_bytes()).hexdigest()=="8d11f09b961beb042f4146f55459141a2372275995b7b476cc85377633f6c4b0","frozen full certificate")

def sg(x): return (x>0)-(x<0)
def direct_coefficients(C,scale):
    """Solve log h by a triangular product, not the producer's recurrence."""
    d=len(C)+2
    coefficients=[]
    for k in range(d+1):
        h=scale**(-comb(k,2))
        for j,c in enumerate(C,start=2):
            if j<k: h*=c**(-comb(k-j+1,2))
        coefficients.append(comb(d,k)*h)
    return coefficients

def primitive_integer(a):
    D=lcm(*(c.denominator for c in a))
    z=[int(c*D) for c in a]
    g=0
    for c in z: g=gcd(g,c)
    return [c//g for c in z]

sturm_rows=0
x=S.symbols("x")
def referee(C,scale,sturm=False):
    global sturm_rows
    d=len(C)+2
    a=direct_coefficients(C,scale)
    integers=primitive_integer(a)
    need(all(v>0 for v in integers),"primitive positive coefficients")
    actualR=[F(integers[k]**2,integers[k-1]*integers[k+1])
             *F(comb(d,k-1)*comb(d,k+1),comb(d,k)**2) for k in range(1,d)]
    need([actualR[k]/actualR[k-1] for k in range(1,d-1)]==C,"target circuits from primitive integer coefficients")
    need(a[0]==1 and a[1]==d,"exact root-sum gauge")
    q=[F(k,d-k+1)*scale**(k-1)
       * product_value([c**(k-j) for j,c in enumerate(C,start=2) if j<k])
       for k in range(1,d+1)]
    need(q==[a[k-1]/a[k] for k in range(1,d+1)],"independent product formula for root sample addresses")
    need(all(q[k]>=9*q[k-1] for k in range(1,d)),"sample geometric separation")
    samples=[F(0)]+[3*v for v in q]
    signs=[]
    for k,point in enumerate(samples):
        terms=[a[j]*point**j for j in range(d+1)]
        value=sum((-1)**j*term for j,term in enumerate(terms))
        need(terms[k]>sum(terms)-terms[k],"literal strict kth-term domination")
        need(sg(value)==(-1)**k,"alternating root-bracketing signs")
        signs.append(sg(value))
    if sturm:
        p=S.Poly(sum(v*x**k for k,v in enumerate(integers)),x)
        need(p.count_roots(-S.oo,0)==d and S.degree(S.gcd(p,p.diff()))==0,
             "independent exact Sturm count and squarefreeness")
        sturm_rows+=1
    return dict(d=d,C=list(map(str,C)),R=list(map(str,actualR)),coefficients=list(map(str,a)),
                samples=list(map(str,samples)),sample_signs=signs)

def product_value(values):
    result=F(1)
    for v in values: result*=v
    return result

def canonical(obj): return json.dumps(obj,separators=(",",":"),sort_keys=True).encode()+b"\n"
total=0
for d,digest in saved["complete_bank_sha256"]:
    rows=[]
    for word in product((-1,0,1),repeat=d-2):
        C=[F(2)**s for s in word]
        row=referee(C,F(9*2**(d-2)),sturm=d<=5)
        need([sg(c-1) for c in C]==list(word),"all positional ternary words including exact zero signs")
        rows.append(row)
    need(sha256(canonical(rows)).hexdigest()==digest,"complete producer bank independently reconstructed")
    need(len(rows)==3**(d-2),"complete ternary universe")
    total+=len(rows)
    print("TERNARY",d,len(rows),"full digest PASS")

for row in saved["arbitrary_positive_targets"]:
    C=list(map(F,row["C"])); scale=F(row["R"][0])
    need(referee(C,scale,sturm=len(C)<=4)==row,"arbitrary nonternary saved target reconstruction")

# Fresh targets and mean scaling check another normalization independently.
for C in ([F(17,3)], [F(1,7),F(49),F(1,11)], [F(2),F(1),F(1,2),F(1)]):
    d=len(C)+2
    prefix=[F(1)]
    for c in C: prefix.append(prefix[-1]*c)
    kappas=[F(k*(d-k),(k+1)*(d-k+1)) for k in range(1,d)]
    scale=9*max(a/b for a,b in zip(kappas,prefix))
    row=referee(C,scale,sturm=True)
    coeff=list(map(F,row["coefficients"]))
    factor=F(13,2*d)
    scaled=[v*factor**k for k,v in enumerate(coeff)]
    ratios=[scaled[k]**2/scaled[k-1]/scaled[k+1]
            *F(comb(d,k-1)*comb(d,k+1),comb(d,k)**2) for k in range(1,d)]
    need(scaled[1]==F(13,2) and [ratios[k]/ratios[k-1] for k in range(1,d-1)]==C,
         "arbitrary fixed positive root sum keeps every circuit")

# Two exactly fixed first ratios give genuine sidecar-loss hostiles.
need(F(3,3)**2/F(F(11,4),3)==F(12,11),"first elementary moments of the literal three-root hostile")
need(F(12,11)*F(1,2)<1,"Newton inequality excludes desired circuit in fixed second-moment fibre")
need(F(13,5)**2/F(55,10)==F(338,275) and F(338,275)/2<1,"actual anchored degree-five moment firewall")
print("STURM",sturm_rows,"independent exact root counts; no floating roots")
print("TOTAL",total,"complete ternary rows; all rational amplitudes and ties retained")
print("PASS",GATES,"always-active independent exact gates; raw LF")

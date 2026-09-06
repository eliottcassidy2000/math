"""Exact norm-jet bridge controls from two previously audited full certificates.

These are pinned polynomial identity data, not imported producer routines.
The analytic theorem has all heights and 0<=h-r<=6; this bank has h6..10.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
import hashlib,json,sys,argparse
sys.stdout.reconfigure(newline='\n')
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)
def evalpoly(c,x):
    ans=Q(0)
    for v in reversed(c):ans=ans*x+v
    return ans
ap=argparse.ArgumentParser()
here=Path(__file__).resolve().parent
ap.add_argument('--data',type=Path,default=here if here.name=='04-computation' else Path('C:/w/s0905/04-computation'))
args=ap.parse_args()
names=['continuing3_20260906_laurent_endpoints_certificate.json',
       'continuing4_20260906_regular_duality_certificate.json']
pins=['b55bec0ba7f1063396af1a8db9be725f279fd01a403adb94b9b026b9490d8e62',
      '0d5a65f03fc4f4295f3db38bca8609375cfea8805f21499a28e9a5e0d9a1ccd4']
data=[]
for name,pin in zip(names,pins):
    raw=(args.data/name).read_bytes()
    check(hashlib.sha256(raw).hexdigest()==pin,'complete inherited certificate pin')
    data.append(json.loads(raw))
regular={row['H']:row for row in data[1]['rows']}
def regular_norm(H,z):
    if H==0:return Q(1)
    entry=regular[H]['residuals'][-1]
    check(entry['k']==H,'full regular norm rather than trace')
    n=evalpoly(list(map(Q,entry['coefficients'])),z)
    for ell in range(1,H+1):n*=((z+2*ell-1)*(z+2*ell))**ell
    check(n>0,'certified regular complementary norm')
    return n
bank=[]
for family in data[0]['families']:
    h=family['h'];entry=family['characteristic'][-1]
    check(entry['k']==h,'full old norm')
    c=[Q(entry['content'])*a for a in entry['coefficients_ascending']]
    for r in range(max(1,h-6),h+1):
        H=h-r;z=2*r+1
        a0=Q(factorial(2*h+1)*factorial(H),factorial(3*H)*factorial(2*r+1))
        A=Q((-1)**(r-1)*factorial(2*h+1)*factorial(H)*factorial(r-1),factorial(3*h))
        B=Q(2*factorial(2*H)*factorial(2*r),factorial(6*h+3))
        expected=B**r*a0**(2*r+1)*regular_norm(H,z)/(A*factorial(4*h+2)**H)
        # Old coefficient data are ascending in x-1. Expand around x=-r.
        jet=[sum(c[j]*comb(j,k)*Q(-r-1)**(j-k) for j in range(k,len(c))) for k in range(r)]
        for low in jet[:-1]:check(low==0,'all lower jet coefficients vanish')
        check(jet[-1]==expected,'exact leading complementary norm jet')
        check((-1)**(r-1)*expected>0,'nonzero sign with full normalization')
        bank.append((h,r,H))
    print('height',h,'exact diagonals',[H for hh,r,H in bank if hh==h],flush=True)
check(len(bank)==34,'complete h6..10 H0..6 bank')
check(Q(2*factorial(0)*factorial(2),factorial(9))==Q(1,90720),'height1 carry boundary constant')
print('exact jet triples',len(bank))
print('PASS',gates,'always-active exact gates')

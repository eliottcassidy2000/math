"""Exact controls for the 2026-10-13 incoming statement repairs. Run with -O too."""
from math import gcd
from pathlib import Path
import hashlib, json, subprocess

checks=0
def check(p, label):
    global checks
    checks+=1
    if not p: raise RuntimeError(label)

def data(n,k):
    x=n; p,q,c=1,1,0; sig=kap=None; cut=0
    for j in range(1,k+1):
        if x%2:
            x=(3*x+1)//2; p,c=3*p,3*c+q
        else: x//=2
        q*=2
        check(q*x==p*n+c, 'exact affine prefix')
        if x<n and sig is None: sig=j
        if p<q:
            cut=max(cut,(c+q-p-1)//(q-p))
            if kap is None:kap=j
    return sig,kap,cut

check(data(1,2)[:2]==(None,2), 'small exceptional source one')
check(data(5,2)[:2]==(2,2), 'same residue, actual stopping differs')
heights=[]
for k in range(1,10):
    H=max(data(r,k)[2] for r in range(1,2**k+1))
    heights.append(H)
    check(H<3**k, 'simple uniform exceptional height')
    for r in range(1,2**k+1):
        m=r+2**k*((H+2**k)//2**k)
        for n in [m,m+2**k]:
            s,c,_=data(n,k)
            check(n>H and s==c, 'actual and coefficient stopping agree above height')

fixed=[]
for a in range(1,55):
    if gcd(a,55)>1:continue
    powers=[]; z=a
    while z!=1:
        powers.append(z);z=z*a%55
    order=len(powers)+1
    if 54 in powers:
        fixed.append([a,order]);check(order in (2,10), 'only fixed residue degrees two and ten')
    if order==20:check(pow(a,10,55)!=54,'degree twenty does not imply conjugation fixed')
check(pow(2,10,55)==34, 'minimal false degree-twenty example')
# At the ramified primes 5 and 11, the surviving unramified factor excludes conjugation.
check(pow(5,5,11)==1 and all(pow(5,j,11)!=10 for j in range(5)), 'no fixed prime above five')
check(11%5==1, 'no fixed prime above eleven')
check(49**2==2*10**3+4*10**2+1, 'original square curve point')
check(98**2==20**3+4*20**2+4, 'correct coordinate change')
check(196**2!=20**3+4*20**2+16, 'old coordinate change rejected')
gp="""E=ellinit([0,4,0,0,4]); print([E.j,ellglobalred(E)[1],ellrank(E),elltors(E)]);
print(ellminimalmodel(E)[1..5]);
F=ellinit([0,4,0,0,16]);print([F.j,ellglobalred(F)[1]]);
G=ellinit([0,109,0,224,0]);print([G.j,ellglobalred(G)[1],ellrank(G),elltors(G)]); print(ellminimalmodel(G)[1..5]);
"""
r=subprocess.run(['gp','-fq'],input=gp,text=True,capture_output=True,check=True)
check(not r.stderr.strip(),'no hidden GP errors')
lines=r.stdout.strip().splitlines()
check(lines[0]=='[-65536/91, 364, [1, 1, 0, [[0, 2]]], [1, [], []]]','correct square curve rank and conductor')
check(lines[1]=='[0, 1, 0, -5, 7]','correct minimal model')
check(lines[2]=='[-4096/43, 43]','different conductor-forty-three curve')
check(lines[3].startswith('[1408317602329/2153060, 910, [1, 1, 0,'),'fruit rank and conductor survive')
check(lines[4]=='[1, 0, 1, -234, 1352]','fruit minimal model survives')
print(json.dumps(dict(checks=checks,depths=list(range(1,10)),exceptional_heights=heights,
 fixed_residue_classes=fixed,gp_source=gp,gp_output=lines,
 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))

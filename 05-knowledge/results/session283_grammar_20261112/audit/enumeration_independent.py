#!/usr/bin/env python3
"""Independent direct sector reconstruction and exact packet consumption."""
import importlib.util,sys,json,math
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from itertools import combinations
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('producer',ROOT/'enumeration/residual_compiler.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def require(p,s):
    if not p:raise RuntimeError(s)
def parts(n,m=1):
    if n==0:yield []
    for a in range(m,n+1,2):
        for p in parts(n-a,a):yield [a]+p
def direct(n,B):
    cells=defaultdict(lambda:[F(0)]*5);answer=F(0);pc=0
    for p in parts(n):
        pc+=1;counts=Counter(p);den=math.prod(a**c*math.factorial(c) for a,c in counts.items())
        exponent=sum((a-1)//2 for a in p);trueexp=exponent;defects=[]
        for a,b in combinations(p,2):
            g=math.gcd(a,b);trueexp+=g
            if a not in B and b not in B:
                exponent+=min(a,b)
                if min(a,b)>g:defects.append(min(a,b)-g)
            else:exponent+=g
        answer+=F(2**trueexp,den)
        c=[F(1)]
        for d in defects:
            q=1-F(1,2**d);new=[F(0)]*(min(4,len(c))+1)
            for i,x in enumerate(c):
                if i<len(new):new[i]+=x
                if i+1<len(new):new[i+1]+=x*q
            c=new
        delta=min(defects,default=0);key=(len(defects),delta,sum(d>delta for d in defects))
        w=F(2**exponent,den)
        for i,x in enumerate(c):cells[key][i]+=w*x
    return cells,answer,pc
def parsed(packet):
    return {(c['M'],c['delta'],c['R']):list(map(F,c['C'])) for c in packet['cells']}
def ground_exact(packet):
    ans=F(0)
    for c in packet['cells']:
        m,delta,r=c['M'],c['delta'],c['R'];cs=list(map(F,c['C']))
        require(r<=4,'not an exact degree-four packet')
        rho=F(1,2**delta);a=1-1/rho;b=1/rho
        for k in range(min(m,4)+1):
            z=sum(math.comb(m-j,k-j)*a**(k-j)*b**j*cs[j] for j in range(k+1))*rho**m
            ans+=(-1)**k*z
    require(ans.denominator==1,'noninteger exact result')
    return int(ans)
def main():
    small=0;large=[]
    for n in range(25):
        for B in [(),(3,),(3,5,7),(3,7)]:
            d,a,c=direct(n,set(B));p,_=mod.produce(n,B)
            require(d==parsed(p),'small packet');require(a.denominator==1,'Burnside integer');small+=1
    for n in [64,70]:
        B=(3,5,7,9,11);d,a,c=direct(n,set(B));p,_=mod.produce(n,B)
        require(d==parsed(p),'large packet')
        large.append({'n':n,'odd_partitions':c,'exact_bits':int(a).bit_length()})
        if n==64:require(ground_exact(p)==a,'exact consumer n64')
    p150=json.loads((ROOT/'enumeration/packet_n150_B35.json').read_text())
    a150=ground_exact(p150)
    completed=json.loads((ROOT/'enumeration/completion.json').read_text())['rows']
    previous=next(r for r in completed if r['n']==150 and r['backbone_top']==27)
    require(a150==int(previous['recovered_integer']),'independent packet versus modular completion')
    p200,_=mod.produce(200,tuple(range(3,50,2)))
    a200=ground_exact(p200)
    receipt200=next(r for r in completed if r['n']==200 and r['backbone_top']==49)
    require(a200==int(receipt200['recovered_integer']),'n200 separate exact consumer')
    # Rebuild complete finite threshold rule without importing threshold code.
    allpatterns=[p for p in combinations(range(0,23,2),4) if sum(p)<=28]
    hist=Counter();firsts={}
    for a in range(29,6959,2):
        valid=[]
        for offset in allpatterns:
            p=[a+x for x in offset]
            ds=[min(x,y)-math.gcd(x,y) for x,y in combinations(p,2)]
            if min(ds)>0 and ds.count(min(ds))==1:valid.append(sum(p))
        require(bool(valid),'threshold existence')
        first=min(valid);hist[first-4*a]+=1;firsts[a]=first
    require(firsts[37]==168,'B35 threshold')
    # Signed fifth finite difference annihilates factorial moments0..4, changes atom0.
    perturb=[(-1)**j*math.comb(5,j) for j in range(6)]
    for k in range(5): require(sum(perturb[j]*math.comb(j,k) for j in range(k,6))==0,'moment kernel')
    print(json.dumps({'status':'PASS','small_complete_packets':small,'large_complete_packets':large,
        'n150_packet_integer_bits':a150.bit_length(),'n150_matches_independent_modular_completion':True,'n200_independent_consumer_bits':a200.bit_length(),
        'threshold_patterns':len(allpatterns),'threshold_odd_classes':3465,'threshold_histogram':dict(hist)},sort_keys=True,indent=2))
if __name__=='__main__':main()

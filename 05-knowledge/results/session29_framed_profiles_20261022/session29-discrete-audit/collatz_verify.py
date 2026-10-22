"""Exact, fail-closed tests; no convergence oracle and no external dependencies."""
from fractions import Fraction
from itertools import product
import json

def require(p, msg):
    if not p: raise RuntimeError(msg)
def v2(n):
    require(n>0,'positive valuation input')
    return (n & -n).bit_length()-1
def orbit(n,h):
    vals=[]; nodes=[n]
    for _ in range(h):
        z=3*n+1; a=v2(z); n=z>>a
        vals.append(a);nodes.append(n)
    return tuple(vals),nodes
def cylinder(w):
    A=0;C=0
    for a in w:C=3*C+(1<<A);A+=a
    M=1<<(A+1)
    r=((1<<A)-C)*pow(3**len(w),-1,M)%M
    return A,C,r,M
def exponent(r,L):
    require(L>=3 and r%8 in (1,3),'power subgroup guard')
    m=0 if r%8==1 else 1
    for j in range(3,L):
        if pow(3,m,1<<(j+1))!=r%(1<<(j+1)):m+=1<<(j-2)
        require(pow(3,m,1<<(j+1))==r%(1<<(j+1)),'binary lift')
    return m,1<<(L-2)
def compositions(n):
    if n==0:yield ();return
    for a in range(1,n+1):
        for w in compositions(n-a):yield (a,)+w
# Every subgroup at these levels, independently of binary lifting.
subgroups=[]
for L in range(3,15):
    M=1<<L;seq=[];x=1
    for m in range(1<<(L-2)):seq.append(x);x=3*x%M
    require(x==1 and len(set(seq))==len(seq),'exact order')
    require(set(seq)=={r for r in range(M) if r%8 in (1,3)},'subgroup image')
    subgroups.append({'L':L,'order':len(seq)})
words=0; compatible=0; oddparts=0
for A in range(2,12):
    for w in compositions(A):
        a,c,r,M=cylinder(w);words+=1
        for t in (0,1,3):
            n=r+t*M
            require(orbit(n,len(w))[0]==w,'valuation cylinder')
            require(orbit(n,len(w))[1][-1]*2**A==3**len(w)*n+c,'carry and clock')
        k,K=exponent((4*r-1)%(4*M),A+3)
        require(k%2==1 and K==M,'oddpart exponent period')
        u=(3**k+1)//4
        require(orbit(u,len(w))[0]==w,'oddpart universal word realization')
        oddparts+=1
        expected=w[0]==2 or (len(w)>=2 and w[0]==1 and w[1]>=2)
        require((r%8 in (1,3))==expected,'allowed prefix classification')
        if expected:
            compatible+=1;m,P=exponent(r,A+1)
            for k in (0,1):
                n=3**(m+k*P)
                require(orbit(n,len(w))[0]==w,'actual power realization')
# Initial reset and exact descent split, literal powers (not logarithmic floats).
initial=[]
for m in range(1,1001):
    n=3**m;w,nodes=orbit(n,2)
    if m%2==0:
        require(w[0]==2 and nodes[1]<n,'even exponent descent')
    else:
        require(w[0]==1 and w[1]==v2(3**(m+2)+5)-1 and w[1]>=2,'odd exponent reset')
        require((nodes[2]<n)==(m%4==1),'odd exponent descent classification')
    if m<=12:initial.append({'m':m,'word':w,'nodes':nodes})
# Unbounded first-descent delay and exact shortcut/ordinary clocks.
delays=[]
for r in range(19):
    w=(1,2)+(1,)*r;A,C,s,M=cylinder(w)
    require(s==(-13*pow(9,-1,M))%M,'closed hostile cylinder')
    m,P=exponent(s,A+1);n=3**m;actual,nodes=orbit(n,len(w))
    require(actual==w and all(x>n for x in nodes[1:]),'all prefixes grow')
    require(A==r+3 and A+len(w)==2*r+5,'two literal clocks')
    delays.append({'r':r,'least_exponent':m,'period':P,'source_bits':n.bit_length(),'odd_steps':len(w),'shortcut_clock':A,'ordinary_clock':A+len(w)})
# Second exponent distribution via modular counting, not actual large powers.
densities=[]
for b in range(2,13):
    M=1<<(b+3);period=1<<(b+1)
    hits=[m for m in range(1,period,2) if (pow(3,m+2,M)+5)%M==(1<<(b+1))]
    require(len(hits)==1,'one selected upper-bit lift at expanded modulus')
    # Exact b requires 3^(m+2)+5 divisible by 2^(b+1), not 2^(b+2).
    hits2=[m for m in range(1,1<<b,2) if (pow(3,m+2,1<<(b+2))+5)%(1<<(b+2))==(1<<(b+1))]
    require(len(hits2)==1,'one exact exponent class')
    densities.append({'second_valuation':b,'residue':hits2[0],'period':1<<b,'all_exponent_density':str(Fraction(1,1<<b))})
oddpart_delays=[]
for h in range(1,18):
    k,P=exponent((-5)%(1<<(h+3)),h+3)
    u=(3**k+1)//4;w,nodes=orbit(u,h)
    require(w==(1,)*h and all(x>u for x in nodes[1:]),'oddpart first-growth prefix')
    oddpart_delays.append({'h':h,'least_odd_exponent':k,'period':P,'source_bits':u.bit_length(),'shortcut_clock':h,'ordinary_clock':2*h})
guards=[]
for k in (1,3,5,7,9,11):
    u=(3**k+1)//4
    guards.append({'odd_exponent_mod12':k,'source_mod15':u%15,'source_mod35':u%35,'source_mod105':u%105})
    for j in range(1,11):
        require(((3**(k+12*j)+1)//4)%105==u%105,'joint guard period')
# CRT compatibility retains common mod4 rather than inventing six independent choices.
for A in range(2,9):
    for w in compositions(A):
        _,_,r,M=cylinder(w);k,P=exponent((4*r-1)%(4*M),A+3)
        residues={(k+j*P)%12 for j in range(3)}
        require(residues=={s for s in range(1,12,2) if s%4==k%4},'three compatible guard lifts')
# Affine debt composition is not source ownership; only check its actual algebra.
for l1,d1,l2,d2 in [(Fraction(2),Fraction(1),Fraction(3),Fraction(2)),(Fraction(1),Fraction(0),Fraction(4),Fraction(1))]:
    l=l1*l2;d=d1+l1*d2
    require(l-1-d==(l1-1-d1)+l1*(l2-1-d2),'weighted slack cocycle')
print(json.dumps({'status':'PASS','word_clock_cap':11,'words':words,'compatible_words':compatible,'oddpart_universal_words':oddparts,'oddpart_delay_receipts':oddpart_delays,'subgroups':subgroups,'initial_examples':initial,'unbounded_delay_receipts':delays,'second_valuation_residues':densities,'mod15_mod35_register':guards,'warning':'Finite tests verify implementations; all-order claims use accompanying proofs. No Collatz convergence claim.'},indent=2))

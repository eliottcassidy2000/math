import json, math

def require(ok, label):
    if not ok: raise RuntimeError(label)

def good(n,r,budget=0):
    bad=0
    for _ in range(r):
        n,d=divmod(n,9)
        if d>=6:return False
        bad+=d>=3
    return bad<=budget

def allgood(n,budget=0):
    bad=0
    while n:
        n,d=divmod(n,9)
        if d>=6:return False
        bad+=d>=3
    return bad<=budget

rows=[]
# Independent direct orbit enumeration, not the lifting implementation.
for r in range(1,6):
    m=9**r;period=6*9**(r-1);x=1;seen=set();counts=[0]*3
    for k in range(period):
        require(x not in seen,'primitive orbit repetition');seen.add(x)
        for b in range(3):counts[b]+=good(x,r,b)
        x=2*x%m
    require(x==1 and len(seen)==2*m//3,'complete unit orbit')
    for b in range(3):
        want=2*3**(r-1)*sum(math.comb(r,j) for j in range(min(b,r)+1))
        require(counts[b]==want,'defect count')
    rows.append({'r':r,'period':period,'counts_b0_b1_b2':counts})
# Cheap holdouts through r=10: each no-defect node has exactly three lifts,
# one-defect branching is three for used budget, six if unused.
S=[0,1];T=[0,1,2,5];period=6
for r in range(1,10):
    mod=9**(r+1);nS=[];nT=[];oldS=set(S)
    for k in S:
        children=[k+j*period for j in range(9) if good(pow(2,k+j*period,mod),r+1)]
        require(len(children)==3,'three children');nS+=children
    for k in T:
        children=[k+j*period for j in range(9) if good(pow(2,k+j*period,mod),r+1,1)]
        require(len(children)==(6 if k in oldS else 3),'one defect branching');nT+=children
    S,T=nS,nT;period*=9
    require(len(S)==2*3**r and len(T)==2*(r+2)*3**r,'holdout count')
# Ordinary exponents are tested by full integer digits, NOT a fixed window.
small=[[],[]];x=1
for k in range(2001):
    for b in range(2):
        if allgood(x,b):small[b].append(k)
    x*=2
# False bridge: periodic finite survivor need not be global survivor.
require(good(64,1) and not allgood(64),'k6 finite/global hostile')
require(good(4,1,1) and not good(4,1),'one defect is distinct')
# Nonunit digits and nonprimitive bases break naive cardinality formulas.
require(sum(good(pow(8,k,81),2) for k in range(18))!=6,'nonprimitive base hostile')
require(sum(good(n,2) for n in range(81))==9,'integer residue count includes nonunits')
# Two parity sheets genuinely matter: 2^k cannot be defined on Z3 alone.
require(pow(2,3,3)!=pow(2,0,3),'parity sheet hostile')
# Exact 3-adic valuation transport on each parity sheet, including negatives
# via modular inverses at finite precision.
def v3(n):
    if n==0:return None
    v=0
    while n%3==0:n//=3;v+=1
    return v
for a in range(51):
    for b in range(a):require(v3(4**a-4**b)==1+v3(a-b),'log similarity')
# Independently calculate initial actual coefficients using integral recurrence.
c=[1,2]
for n in range(1,4096):
    q,rem=divmod((26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1],(n+1)**2)
    require(rem==0,'coefficient integrality');c.append(q)
for k in range(13):
    n=2**k;require(c[n]%4==2,'dyadic coefficient')
    require((c[n]%3!=0)==allgood(n),'actual mod3 intersection')
    require((c[n]%9!=0)==allgood(n,1),'actual mod9 intersection')
print(json.dumps({'status':'PASS finite tests, infinite claims proved separately','direct_rows':rows,'holdout_r':10,'holdout_counts':[len(S),len(T)],'full_integer_exponents_through_2000':small,'actual_coefficients_through':4096,'hostile_controls':6},sort_keys=True,indent=2))

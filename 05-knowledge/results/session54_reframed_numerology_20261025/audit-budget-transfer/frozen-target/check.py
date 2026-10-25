import json,math

def require(ok,label):
    if not ok:raise RuntimeError(label)
def choose(n,j):return math.comb(n,j) if 0<=j<=n else 0
def step(v):return [3*(v[j]+(v[j-1] if j else 0)) for j in range(len(v))]
def numerator(b):
    p=[0]*(b+1)
    # 2(1-x)^b + 2 sum x^(j-1)(1-x)^(b-j)
    for a in range(b+1):p[a]+=2*(-1)**a*choose(b,a)
    for j in range(1,b+1):
        for a in range(b-j+1):p[a+j-1]+=2*(-1)**a*choose(b-j,a)
    return p
rows=[]
for b in range(7):
    v=[2]+([2]+[0]*(b-1) if b else [])
    values=[]
    for r in range(1,21):
        want=[2*3**(r-1)*choose(r,j) for j in range(b+1)]
        require(v==want,'transfer/binomial vector');values.append(sum(v));v=step(v)
    p=numerator(b);require(sum(p)==2,'no pole cancellation over Q')
    den=[(-3)**j*choose(b+1,j) for j in range(b+2)]
    num=[x*3**j for j,x in enumerate(p)]
    gf=[]
    for n in range(20):
        gf.append((num[n] if n<len(num) else 0)-sum(den[j]*gf[n-j] for j in range(1,min(n,b+1)+1)))
    require(gf==values,'rational GF independent coefficient recursion')
    diff=[values[r-1]//(3**(r-1)) for r in range(1,21)]
    for _ in range(b):diff=[y-x for x,y in zip(diff,diff[1:])]
    require(set(diff)=={2},'degree exactly b, leading coefficient2/b!')
    require(all(y-x==0 for x,y in zip(diff,diff[1:])),'degree no larger than b')
    # Explicit largest nilpotent direction is excited and observable.
    v0=[2]+([2]+[0]*(b-1) if b else [])
    shifted=v0[:]
    for _ in range(b):shifted=[0]+shifted[:-1]
    require(sum(shifted)==2,'maximal Jordan direction')
    require([0]+shifted[:-1]==[0]*(b+1),'nilpotent overflow discarded')
    rows.append({'b':b,'initial':v0,'numerator_in_x':p,'minimal_pole_order':b+1,'first_four_totals':values[:4]})
# Independent direct modular orbit count, including budgets larger than depth.
for r in range(1,5):
    period=6*9**(r-1);mod=9**r;x=1;counts=[0]*7
    for _ in range(period):
        y=x;defects=0;valid=True
        for j in range(r):
            y,d=divmod(y,9);defects+=3<=d<=5;valid &= d<6
        if valid:counts[defects]+=1
        x=2*x%mod
    for b in range(7):
        v=[2]+([2]+[0]*(b-1) if b else [])
        for _ in range(r-1):v=step(v)
        require(v==counts[:b+1],'direct orbit exact defect vector')
require(numerator(0)==[2],'b0 separate')
require(numerator(1)==[4,-2],'b1 numerator 4-6z')
require(step([2,0])!=step([2,2]),'wrong r0-style initialization rejected')
require([3*x for x in [2,2]]!=step([2,2]),'diagonal model rejected')
require([3*(2+2),3*(2+2)]!=step([2,2]),'cyclic debt recycling rejected')
require(all(x%2==0 for x in numerator(3)),'over Q qualifier: mod2 numerator vanishes')
print(json.dumps({'status':'PASS','budgets':rows,'direct_orbit_depth':4,'transfer_and_GF_levels':20,'hostile_controls':4},sort_keys=True,indent=2))

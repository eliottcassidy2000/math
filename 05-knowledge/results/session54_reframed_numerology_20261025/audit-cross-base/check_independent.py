import math,json

def need(ok,label):
    if not ok: raise RuntimeError(label)

def digits(n):
    out=[]
    while n:
        out.append(n%9);n//=9
    return out or [0]

def accept(n,b):
    ds=digits(n)
    return max(ds)<6 and sum(d//3 for d in ds)<=b

def bsgs(g,h,m,order):
    s=math.isqrt(order)+1
    baby={};a=1
    for j in range(s):baby.setdefault(a,j);a=a*g%m
    step=pow(pow(g,s,m),-1,m);a=h
    for i in range(s+1):
        if a in baby:
            k=i*s+baby[a]
            if k<order:return k
        a=a*step%m
    raise RuntimeError('missing discrete log')

rows=[]
for r in range(1,7):
    mod=9**r;period=2*3**(2*r-1)
    need(pow(2,period,mod)==1 and pow(2,period//2,mod)!=1 and pow(2,period//3,mod)!=1,'exact order')
    candidates=[(d,int(d>=3)) for d in (1,2,4,5)]
    for j in range(1,r):
        candidates=[(a+d*9**j,b+int(d>=3)) for a,b in candidates for d in range(6) if b+int(d>=3)<=2]
    exponents={}
    for a,b in candidates:
        k=bsgs(2,a,mod,period)
        need(pow(2,k,mod)==a,'inverse exponent')
        need(k not in exponents,'injectivity');exponents[k]=b
    counts=[sum(b<=budget for b in exponents.values()) for budget in range(3)]
    need(counts==[2*3**(r-1)*sum(math.comb(r,j) for j in range(min(b,r)+1)) for b in range(3)],'counts')
    rows.append([r,period,counts])
# Direct defining binomial sum, then square-root convolution, not recurrence.
c=[]
for n in range(121):
    s=sum(math.comb(n,j)**2*math.comb(n+j,j)*math.comb(2*j,n) for j in range((n+1)//2,n+1))
    if n==0:c.append(1);continue
    q,rem=divmod(s-sum(c[j]*c[n-j] for j in range(1,n)),2)
    need(rem==0,'integral square root');c.append(q)
    need((q%3!=0)==accept(n,0),'mod3 law')
    need((q%9!=0)==accept(n,1),'mod9 law')
    need(q%4==(2 if n&(n-1)==0 else 0),'mod4 law')
need(c[0]==1 and c[1]==2,'constant vs exponent zero')
# Bounded full integer scan is only a countermodel guard.
survivors=[[k for k in range(2001) if accept(1<<k,b)] for b in range(2)]
need(survivors==[[0,1,13],[0,1,2,7,13]],'finite scan')
need(digits(8192)==[2,1,2,2,1],'k13')
need(accept(pow(2,6,9),0) and not accept(2**6,0),'finite not global')
for r in range(1,12):need(pow(2,6*9**(r-1),9**r)==1,'moving witness')
for K in range(1,100001):
    r=1
    while 6*9**(r-1)<K:r+=1
    need((2*3**(r-1))**2<=6*K,'uniform bound constant')
# Negative and positive z on each parity sheet, finite precision.
def valuation(n):
    v=0
    while n and n%3==0:v+=1;n//=3
    return v
for e in (0,1):
    for z in range(-30,31):
        for w in range(-30,z):
            delta=(2**e*(pow(4,z,3**9)-pow(4,w,3**9)))%3**9
            need(valuation(delta)==1+valuation(z-w),'sheet similarity')
need(pow(2,3,3)!=1,'one-sheet obstruction')
need(sum(accept(a,0) for a in range(81))==9,'nonunit trap')
need(sum(accept(pow(8,k,81),0) for k in range(18))!=6,'nonprimitive trap')
need(accept(4,1) and not accept(4,0),'defect budget trap')
# At odd 3-adic precision a new low ternary digit is fixed too;
# ignoring it yields the report's valid (non-sharp) Frostman bound.
for m in range(1,15):
    r=(m+1)//2
    values=[0]
    for j in range(r):values=[a+d*9**j for a in values for d in range(3)]
    bins={}
    for a in values:bins[a%3**m]=bins.get(a%3**m,0)+1
    need(max(bins.values())*3**(m//2)<=len(values),'Frostman cylinder bound')
print(json.dumps({'status':'PASS','independent_digit_enumeration_and_BSGS':rows,'binomial_square_coefficients_through':120,'finite_survivors_through_2000':survivors,'interval_constant_K_through':100000,'negative_and_positive_similarity':True,'Frostman_precision_through':14},indent=2,sort_keys=True))

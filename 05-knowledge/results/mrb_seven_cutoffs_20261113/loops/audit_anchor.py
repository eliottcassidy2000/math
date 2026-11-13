"""Independent audit of sibling anchor classifiers. No sibling imports."""
from collections import Counter,deque
from itertools import product
import json,math

def ck(ok,msg):
    if not ok:raise RuntimeError(msg)

def step(v,c,m=None):
    x,b,d=v
    w=(x+2*b,d-b,x+d) if c=='L' else (2*x+2*b,0,d)
    return tuple(t%m for t in w) if m else w

def phase_direct(v,p,period):
    z=v;mask=0
    for j in range(period):
        if z[1]==0:mask|=1<<(j%2)
        z=step(z,'L',p)
    ck(z==v,'claimed order')
    return mask

def formula7(v):
    x,b,d=v;a=(6*x+3*b+5*d)%7;u=(b-a)%7;w=(5*(d-b)-a)%7;N=(6*u*u+3*u*w+w*w)%7
    return int((a*a-4*N)%7 in {0,3,5,6})+2*int((a*a+4*N)%7 in {0,3,5,6})

def formula11(v):
    x,b,d=v;ell=(2*x+6*b+d)%11
    if ell:return 3
    a=(4*x-b)%11;c=(b-3*x)%11
    if not a and not c:return 3
    if not a or not c:return 0
    z=6*c*pow(a,-1,11)%11
    return 1 if z in {1,3,4,5,9} else 2

def closure(m):
    def enc(v):return v[0]+m*(v[1]+m*v[2])
    def dec(i):x=i%m;i//=m;return x,i%m,i//m
    seed=enc((2%m,0,1%m));seen=bytearray(m**3);seen[seed]=1;q=deque([seed]);count=0
    while q:
        v=dec(q.popleft());count+=1
        for c in 'LK':
            w=enc(step(v,c,m))
            if not seen[w]:seen[w]=1;q.append(w)
    return seen,enc,dec,count

def main():
    phases={};local={}
    for p,period,f in [(7,48,formula7),(11,110,formula11)]:
        by={v:phase_direct(v,p,period) for v in product(range(p),repeat=3)}
        ck(all(f(v)==e for v,e in by.items()),f'classifier{p}')
        phases[p]=by;local[p]=dict(sorted(Counter(by.values()).items()))
    modular={}
    for m in [7,11,14,77]:
        seen,enc,dec,count=closure(m)
        for idx in range(m**3):
            v=dec(idx)
            if m==7:expected=bool(formula7(v))
            elif m==11:expected=bool(formula11(v))
            elif m==14:expected=v[0]%2==0 and v[2]%2==1 and bool(formula7(tuple(t%7 for t in v))&(1<<(v[1]%2)))
            else:expected=bool(formula7(tuple(t%7 for t in v))&formula11(tuple(t%11 for t in v)))
            ck(bool(seen[idx])==expected,f'closure{m}:{v}')
        plane=[(x,0,d) for x in range(m) for d in range(m) if m!=14 or (x%2==0 and d%2==1)]
        ck(all(seen[enc(v)] for v in plane),f'plane reach{m}')
        modular[m]={'reachable':count,'plane':len(plane)}
        if m==77:ck(not seen[enc((5,1,72))],'CRT phantom')
    ready=0;words=0;first=None
    for length in range(12):
        for letters in product('LK',repeat=length):
            word=''.join(letters);v=(2,0,1)
            for c in word:v=step(v,c)
            words+=1;x,b,d=v
            expect=bool(word) and word.endswith('L') and set(word[:-1])<={'K'}
            ck((b==1)==expect,'integral ready '+word)
            ck(x>=2 and d>=1 and d-b>=1,'positive cone')
            if expect:ready+=1
            if b%7==1 and (x+d)%7==0 and first is None:first={'length':length,'word':word,'state':v}
    ck(first['length']==6,'first false-ready length')
    q=(2,0,1)
    for c in 'KLLLKL':q=step(q,c)
    ck(q==(38,15,53),'named false ready')
    # Direct field model: t^2=4t+1 over F7; pair (a,b) means a+b*t.
    def mul(z,w):return ((z[0]*w[0]+z[1]*w[1])%7,(z[0]*w[1]+z[1]*w[0]+4*z[1]*w[1])%7)
    t=(0,1);z=(1,0);orbit=[]
    for _ in range(16):orbit.append(z);z=mul(z,t)
    ck(z==(1,0) and len(set(orbit))==16 and orbit[8]==(6,0),'norm circle order')
    for A,B in product(range(7),repeat=2):
        u=(2*A+4*B)%7;w=(4*A+4*B)%7;N=(A*A+4*A*B-B*B)%7
        ck(N==(6*u*u+3*u*w+w*w)%7,'norm elimination')
    print(json.dumps({'status':'PASS','local_phase_counts':local,'complete_modular_closures':modular,'integral_words':words,'defect_one_words':ready,'first_false_ready':first,'norm_pairs':49,'crt_gcd':math.gcd(48,110),'crt_excluded_product_states':337*1311-429807},indent=2,sort_keys=True))

if __name__=='__main__':main()

"""Exact physical observation charts and nonautonomous scalar compiler."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import json,random
from math import gcd
from itertools import combinations

def check(ok,msg):
    if not ok: raise RuntimeError(msg)
@lru_cache(None)
def g(p):
    if p==0:return 1
    v=p//2
    return g(v) if p%2 else g(v)+v+(v%2)
def params(w):
    p=int(w,2);q=sum((b=='0')<<j for j,b in enumerate(w))
    return g(p),p+1,q+1,1<<len(w),int(q%2==0),int(p%2==0)
def step(x,w):
    A,W,d=x;h,a,b,t,c,f=params(w)
    return h*A+b*W-c*d,a*A+t*W,f*A
def packet(ws):
    x=(1,0,0)
    for w in ws:x=step(x,w)
    return x
def chart(x):
    A,W,d=x;B=step(x,'100')[0];C=step(x,'110')[0]
    check((B-A)%2==0 and (C-A)%2==0,'parity')
    return A,(B-A)//2,(C-A)//2
def unchart(z):
    A,Y,Z=z
    return A,Y-Z,4*A+5*Y-7*Z

def observe_inverse(counts,w1,w2):
    A0,A1,A2=counts;h1,a1,b1,t1,_,_=params(w1);h2,a2,b2,t2,_,_=params(w2)
    numer=A2-h2*A1+(1-b2*a1)*A0
    if numer%(b2*t1):return None
    W=numer//(b2*t1)
    return A0,W,h1*A0+b1*W-A1

def recurrence(Ajm2,Ajm1,Aj,wj,wjp1):
    h,a,b,t,_,_=params(wj);H,aa,B,tt,_,_=params(wjp1)
    numer=(b*H+B*t)*Aj+(b*B*a-b-B*t*h)*Ajm1+B*t*Ajm2
    if numer%b:return None
    return numer//b

def graph(ws):
    rows=[];roots=[]
    for w in ws:
        z=len(rows);roots.append(z);rows.extend([0]*(len(w)+1))
        for i in range(z):
            for j in range(z,len(rows)):rows[i]|=1<<j
        for i,bit in enumerate(w,z+1):
            if bit=='1':rows[z]|=1<<i
            else:rows[i]|=1<<z
            for j in range(i+1,len(rows)):rows[i]|=1<<j
    for i,j in zip(roots,roots[1:]):rows[i]^=1<<j;rows[j]|=1<<i
    return rows,roots[-1]
def literal(ws,full_packet=False):
    rows,r=graph(ws);n=len(rows);full=(1<<n)-1
    for i in range(n):
        check(not rows[i]>>i&1,'loop')
        for j in range(i):check(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,'pair')
    @lru_cache(None)
    def paths(rem,v):
        if not rem:return 1
        ans=0;todo=rem&rows[v]
        while todo:
            bit=todo&-todo;todo-=bit
            ans+=paths(rem^bit,bit.bit_length()-1)
        return ans
    @lru_cache(None)
    def H(s):
        if not s:return 1
        ans=0;todo=s
        while todo:
            bit=todo&-todo;todo-=bit
            ans+=paths(s^bit,bit.bit_length()-1)
        return ans
    a=H(full)
    if not full_packet:return a
    states={(1<<r,r):1};odd=even=0
    for size in range(1,n+1):
        nxt={}
        for (s,v),mult in states.items():
            val=mult*H(full^s)
            if size%2:odd+=val
            else:even+=val
            todo=rows[v]&~s
            while todo:
                bit=todo&-todo;todo-=bit;k=(s|bit,bit.bit_length()-1)
                nxt[k]=nxt.get(k,0)+mult
        states=nxt
    return a,odd+even,odd-even

def mm(A,B):return [[sum(A[i][j]*B[j][k] for j in range(3))for k in range(3)]for i in range(3)]
def matrix(w):
    h,a,b,t,c,f=params(w)
    return [[h,b,-c],[a,t,0],[f,0,0]]
def smith(P):
    d1=gcd(*(x for row in P for x in row));pairs=list(combinations(range(3),2))
    d12=gcd(*(P[i][a]*P[j][b]-P[i][b]*P[j][a] for i,j in pairs for a,b in pairs))
    det=sum(P[0][j]*(P[1][(j+1)%3]*P[2][(j+2)%3]-P[1][(j+2)%3]*P[2][(j+1)%3]) for j in range(3))
    return d1,d12//d1,abs(det)//d12

def main():
    words=[''.join(w) for n in range(1,4) for w in product('01',repeat=n)]
    seeds=[ws for k in (1,2) for ws in product(words,repeat=k)]
    literal_checks=probe_checks=0
    for ws in seeds:
        x=packet(ws);check(x==literal(ws,True),('packet',ws));literal_checks+=1
        z=chart(x);check(unchart(z)==x,('chart',ws))
        for w in ('100','110'):
            check(step(x,w)[0]==literal(ws+(w,)),('probe',ws,w));probe_checks+=1
    active=['1'+''.join(w)+'0' for n in range(2,9) for w in product('01',repeat=n-2)]
    # Exact all integral pairs in this finite family, including the non-Pisot block.
    pair_checks=0
    for i,w1 in enumerate(active):
        for j,w2 in enumerate(active):
            x=packet(seeds[(17*i+j)%len(seeds)]);y=step(x,w1);z=step(y,w2)
            check(observe_inverse((x[0],y[0],z[0]),w1,w2)==x,'three-count inverse')
            pair_checks+=1
    rng=random.Random(2807);recurrence_checks=0;maxbits=0;smith_checks=0
    for trial in range(120):
        ws=[rng.choice(active) for _ in range(100)]
        if trial%3==0:ws[::2]=['10101010']*50
        x=packet(rng.choice(seeds));states=[x]
        for w in ws:states.append(step(states[-1],w))
        A=[z[0] for z in states]
        for j in range(2,len(ws)):
            check(recurrence(A[j-2],A[j-1],A[j],ws[j-1],ws[j])==A[j+1],'recurrence')
            recurrence_checks+=1
        maxbits=max(maxbits,A[-1].bit_length())
        P=[[1,0,0],[0,1,0],[0,0,1]];ell=0
        for k,w in enumerate(ws,1):
            P=mm(matrix(w),P);ell+=len(w)
            d1,d2,d3=smith(P)
            check(d1==1 and d2*d3==1<<ell and d2%(1<<(k-1))==0 and d3%d2==0,'Smith bound')
            check(d2&(d2-1)==0,'Smith dyadic')
            smith_checks+=1
    hostile=[]
    for ws in [('100','0','110'),('100','1','110')]:
        x=packet(ws);counts=[]
        for k in range(3):
            cnt=literal(ws+('10',)*k)
            check(cnt==x[0],'hostile literal');counts.append(cnt)
            x=step(x,'10')
        hostile.append({'words':ws,'packet':packet(ws),'counts':counts})
    check(hostile[0]['counts'][:2]==hostile[1]['counts'][:2],'hostile same pair')
    check(hostile[0]['counts'][2]!=hostile[1]['counts'][2],'hostile difference')
    check(observe_inverse((3,17,1000),'10','10') is None,'invalid lattice guard')
    for A in range(1,20,2):
        for W in range(A,30):
            for d in range(A+1):
                if (W-d)%2==0:check(unchart(chart((A,W,d)))==(A,W,d),'cone chart')
    P=[[1,0,0],[0,1,0],[0,0,1]]
    for k in range(1,41):
        P=mm(matrix('100'),P)
        check(smith(P)==(1,1<<(k-1),1<<(2*k+1)),'sharp Smith family')
    data={'status':'PASS','smith_products':smith_checks,'sharp_smith_depths':40,'literal_root_packets':literal_checks,'literal_equal_order_probes':probe_checks,
          'active_word_lengths':[2,8],'active_words':len(active),'three_count_pairs':pair_checks,
          'random_schedules':120,'blocks_per_schedule':100,'recurrence_instances':recurrence_checks,
          'largest_final_count_bits':maxbits,'hostile':hostile,
          'formal_forbidden_packet':[7,7,1],'formal_probe_counts':[83,69],
          'scope':'Exact identities prove all-word statements; finite tests do not prove Chan-Pak coverage.'}
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()

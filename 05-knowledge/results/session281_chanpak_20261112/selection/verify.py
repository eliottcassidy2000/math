"""Exact certificates; reconstructed after workspace reset. Standard library only."""
from pathlib import Path
from fractions import Fraction as F
from collections import deque
from math import gcd
import json,argparse
HERE=Path(__file__).resolve().parent
CERT=json.loads((HERE/'certificate.json').read_text())
I=((1,0,0),(0,1,0),(0,0,1));ZERO=((0,0,0),)*3

def need(t,msg):
    if not t:raise ValueError(msg)
def mm(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mv(a,v):return tuple(sum(a[i][j]*v[j] for j in range(3)) for i in range(3))
def add(a,b):return tuple(tuple(a[i][j]+b[i][j] for j in range(3)) for i in range(3))
def scale(c,a):return tuple(tuple(c*x for x in row) for row in a)
def mod(a,m):return tuple(tuple(F(x).numerator*pow(F(x).denominator,-1,m)%m for x in row) for row in a)
def power(a,n,m=None):
    out=I
    while n:
        if n&1:out=mm(out,a);out=mod(out,m) if m else out
        a=mm(a,a);a=mod(a,m) if m else a;n//=2
    return out

def inv(a):
    b=[[F(x) for x in a[i]]+[F(x) for x in I[i]] for i in range(3)]
    for j in range(3):
        p=next(i for i in range(j,3) if b[i][j]);b[p],b[j]=b[j],b[p]
        c=b[j][j];b[j]=[x/c for x in b[j]]
        for i in range(3):
            if i!=j:
                c=b[i][j];b[i]=[x-c*y for x,y in zip(b[i],b[j])]
    return tuple(tuple(row[3:]) for row in b)
def e(i,j,c=1):
    a=[list(row) for row in I];a[i][j]+=c;return tuple(map(tuple,a))
def g(p):
    if not p:return 1
    v=p//2;return g(v) if p%2 else g(v)+v+v%2

def block(w):
    p=int(w,2);b=1+sum(1<<i for i,c in enumerate(w) if c=='0')
    return ((g(p),b,-1),(p+1,1<<len(w),0),(1,0,0))
M={w:block(w) for w in CERT['alphabet']}
def word(ws,m=None):
    out=I
    for w in ws:out=mm(M[w],out);out=mod(out,m) if m else out
    return out

def rank3(vectors):
    a=[list(row) for row in zip(*vectors)];r=0
    for j in range(len(vectors)):
        p=next((i for i in range(r,len(a)) if a[i][j]%3),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];c=pow(a[r][j],-1,3);a[r]=[x*c%3 for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                c=a[i][j];a[i]=[(x-c*y)%3 for x,y in zip(a[i],a[r])]
        r+=1
    return r

def all_orbit(m,stop_scalar=False):
    s=(1%m,)*3;prev={s:None};q=deque([s]);vals={s[0]:s}
    expected=set(range(m)) if m%2 else set(range(1,m,2))
    while q and (not stop_scalar or set(vals)!=expected):
        x=q.popleft()
        for w,a in M.items():
            y=tuple(z%m for z in mv(a,x))
            if y not in prev:prev[y]=(x,w);q.append(y);vals.setdefault(y[0],y)
    return prev,vals

def root(k):
    lam=1
    for j in range(1,k):
        f=lambda z:z**3-13*z*z+6*z-8
        if f(lam)%(1<<(j+1)):lam+=1<<j
        need(f(lam)%(1<<(j+1))==0,'lifted root')
    return lam

def literal_packet(words):
    # Inherited Session279 graph builder/subset DP, with the singleton explicit.
    rows=[];roots=[]
    for w in ('',)+tuple(words):
        z=len(rows);roots.append(z);n=len(w)+1;rows.extend([0]*n)
        for a in range(z):
            for b in range(z,z+n):rows[a]|=1<<b
        for i,c in enumerate(w,1):
            if c=='1':rows[z]|=1<<(z+i)
            else:rows[z+i]|=1<<z
            for j in range(i+1,n):rows[z+i]|=1<<(z+j)
    for a,b in zip(roots,roots[1:]):rows[a]^=1<<b;rows[b]|=1<<a
    n=len(rows);full=(1<<n)-1;root=roots[-1]
    dp=[[0]*n for _ in range(1<<n)];rp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    rp[1<<root][root]=1
    for s in range(1,1<<n):
        for v in range(n):
            todo=rows[v]&~s
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                dp[s|bit][j]+=dp[s][v];rp[s|bit][j]+=rp[s][v]
    h=[sum(r) for r in dp];h[0]=1;o=ee=0
    for s in range(1,1<<n):
        term=sum(rp[s])*h[full^s]
        if s.bit_count()%2:o+=term
        else:ee+=term
    return h[full],o+ee,o-ee

def recover_word(prev,x):
    ws=[]
    while prev[x] is not None:x,w=prev[x];ws.append(w)
    return ws[::-1]

def main():
    T,B,C=M['10'],M['100'],M['1000'];s=(1,1,1);s1=mv(T,s)
    floors=[a[0][0]+a[0][1]-1 for a in M.values()];weights=[276,147,173,82,95,82,104]
    need(floors==[5,11,9,23,19,23,17],'growth floors')
    need(sum(weights)==959 and all(u**5*z**4>=1000**5 for u,z in zip(weights,floors)),'4/5 scarcity')
    need(max(a[0][0]+2*a[0][1] for a in M.values())==39,'upper growth')
    U=mm(M['1010'],inv(C));V=mm(M['1100'],inv(C))
    need(U==((1,F(-1,4),F(9,4)),(0,1,2),(0,0,1)),'U')
    need(V==((1,F(-1,8),F(25,8)),(0,1,4),(0,0,1)),'V')
    comm=mm(mm(mm(U,V),inv(U)),inv(V));need(comm==e(0,2,F(-3,4)),'commutator')
    need(mm(power(U,2),inv(V))==((1,F(-3,8),F(19,8)),(0,1,0),(0,0,1)),'12 upper')
    need(mm(U,power(inv(V),2))==((1,0,F(-7,2)),(0,1,-6),(0,0,1)),'23 upper')
    need(mm(mm(inv(T),e(0,2,3)),T)==e(2,0,-3),'31 conjugation')
    need(power(inv(comm),4)==e(0,2,3),'integral transvection')
    for ij,ws in CERT['elementary_mod3'].items():need(word(ws,3)==e(int(ij[0]),int(ij[1])),'base '+ij)
    N=add(e(0,2),scale(-1,I));vectors=[]
    for ws in CERT['conjugate_words']:
        G=word(ws);X=mm(mm(G,N),inv(G));need(mm(X,X)==ZERO,'square zero')
        need(sum(X[i][i] for i in range(3))==0,'trace zero')
        vectors.append(tuple(z for row in mod(X,3) for z in row))
        for j in range(1,5):need(power(add(I,scale(3,X)),3**(j-1))==add(I,scale(3**j,X)),'lift power')
    need(rank3(vectors)==8,'Lie span')
    need(add(add(power(B,3),scale(-13,power(B,2))),add(scale(6,B),scale(-8,I)))==ZERO,'Cayley Hamilton')
    dyadic=[]
    for k in range(1,25):
        m=1<<k;lam=root(k);a=lam-13;c=8*pow(lam,-1,m)%m
        Q=mod(add(add(power(B,2),scale(a,B)),scale(c,I)),m)
        den=(lam*lam+a*lam+c)%m;need(den%2==1,'projector unit')
        P=mod(scale(pow(den,-1,m),Q),m)
        need(mod(mm(P,P),m)==P,'idempotent')
        need(mod(mm(B,P),m)==mod(scale(lam,P),m),'eigenline')
        need(mod(mm(power(B,2*k,m),add(I,scale(-1,P))),m)==ZERO,'transient')
        v=(1,5*pow(lam-8,-1,m)%m,pow(lam,-1,m))
        need(tuple(z%m for z in mv(B,v))==tuple(lam*z%m for z in v),'normalized line')
        amp=[mv(P,x)[0]%m for x in (s,s1)]
        if k>=3:
            need(lam%8==3 and [x%8 for x in amp]==[1,7],'cosets');order=1<<(k-2)
            need(pow(lam,order,m)==1 and pow(lam,order//2,m)!=1,'unit order')
        if k<=12:
            seen=set()
            for x in amp:
                y=x
                for j in range(max(1,m//2)):seen.add(y);y=y*lam%m
            need(seen==set(range(1,m,2)),'dyadic scalar cover')
        dyadic.append({'k':k,'lambda':lam,'amplitudes':amp})
    need(mod(power(B,4),8)==mod(power(B,6),8),'early projector')
    odd=[]
    for m in (3,5,7,9,11,15,21,25,27):
        prev,vals=all_orbit(m)
        expected=sum(gcd(gcd(gcd(a,b),c),m)==1 for a in range(m) for b in range(m) for c in range(m))
        need(len(prev)==expected,'primitive orbit size')
        need(all(gcd(gcd(gcd(x[0],x[1]),x[2]),m)==1 for x in prev),'primitive necessity')
        odd.append({'m':m,'packets':len(prev),'scalar_residues':len(vals)})
    simultaneous=[];examples=[]
    for m in (8,16,32,64,128,24,40,48,72,80,120,144,240,288,360):
        prev,vals=all_orbit(m,True);need(set(vals)==set(range(1,m,2)),'joint scalar cover')
        simultaneous.append({'m':m,'discovered_packets':len(prev),'residues':len(vals)})
        if m in (24,120,360):
            for target in (7,21):
                ws=recover_word(prev,vals[target]);actual=mv(word(ws),s)
                need(actual[0]%m==target and actual[0] not in (7,21),'absolute hostile')
                examples.append({'modulus':m,'target_residue':target,'words':ws,'actual_packet':actual})
    literal=0
    for ws in [()]+[(w,) for w in M]+[('10','100'),('100','10'),('110','100')]:
        need(literal_packet(ws)==mv(word(ws),s),'literal packet');literal+=1
    out={'status':'finite controls for all-modulus proof; no Chan-Pak coverage',
         'scarcity_weights_over_1000':weights,'scarcity_sum_numerator':959,'scarcity_exponent':'4/5',
         'rational_identity_checks':7,'mod3_elementary_words':6,'lie_span_rank':8,
         'dyadic_precisions':dyadic,'odd_packet_orbits':odd,'simultaneous_scalar_controls':simultaneous,
         'hostile_absolute_examples':examples,'literal_controls':literal}
    print(json.dumps(out,indent=2,sort_keys=True))

def residue_witness(m,target):
    need(m>=1,'positive modulus');target%=m;need(m%2 or target%2,'even-modulus target must be odd')
    seed=(1%m,)*3;prev={seed:None};queue=deque([seed]);found=None
    while queue:
        x=queue.popleft()
        if x[0]==target:found=x;break
        for w,a in M.items():
            y=tuple(z%m for z in mv(a,x))
            if y not in prev:prev[y]=(x,w);queue.append(y)
    need(found is not None,'target not reached');ws=recover_word(prev,found);packet=mv(word(ws),(1,1,1))
    need(packet[0]%m==target,'witness residue')
    return {'modulus':m,'residue':target,'words':ws,'actual_packet':packet,
            'tournament_order':1+sum(len(w)+1 for w in ws),'discovered_states':len(prev),
            'scope':'congruence witness, not realization of the integer residue'}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--witness',nargs=2,type=int)
    args=parser.parse_args()
    if args.witness:print(json.dumps(residue_witness(*args.witness),indent=2,sort_keys=True))
    else:main()

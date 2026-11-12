"""Rational triangle-run observability certificates; Python standard library only."""
from fractions import Fraction as F
from itertools import product, combinations
from functools import lru_cache
import json

M=((3,3,-1),(3,4,0),(1,0,0))
RAYS=((1,1,0),(1,1,1),(0,1,0))
I=((1,0,0),(0,1,0),(0,0,1))
PAIRS=tuple(combinations(range(3),2))
T=tuple(tuple(M[i][a]*M[j][b]-M[i][b]*M[j][a] for a,b in PAIRS) for i,j in PAIRS)
def mv(A,v): return tuple(sum(a*b for a,b in zip(row,v)) for row in A)
def mm(A,B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def power(A,n):
    out=I
    while n:
        if n&1: out=mm(out,A)
        A=mm(A,A); n//=2
    return out

def cone(x):
    a,w,d=x
    return 0<=d<=a<=w

def inverse(x):
    a,w,d=x
    if (w-3*d)%4:return None
    parent=(d,(w-3*d)//4,3*d+3*((w-3*d)//4)-a)
    return parent if cone(parent) else None

def geometry(k):
    P=power(M,k); rays=tuple(mv(P,r) for r in RAYS)
    assert k>=1 and all(v[0]>0 for v in rays)
    lo=tuple(min(F(v[j],v[0]) for v in rays) for j in (1,2))
    hi=tuple(max(F(v[j],v[0]) for v in rays) for j in (1,2))
    amin=min(v[0] for v in rays);amax=max(v[0] for v in rays)
    E=max(abs(x[i]*y[j]-x[j]*y[i]) for x,y in combinations(rays,2) for i,j in PAIRS)
    return rays,lo,hi,amin,amax,E

def ceil(q):return -((-q.numerator)//q.denominator)
def recover_packet(n,k):
    rays,lo,hi,*_=geometry(k)
    ranges=[(ceil(n*l),(n*h).numerator//(n*h).denominator) for l,h in zip(lo,hi)]
    if any(l>h for l,h in ranges):return 'reject',None
    if any(l<h for l,h in ranges):return 'ambiguous',None
    packet=(n,ranges[0][0],ranges[1][0])
    x=packet
    for _ in range(k):
        x=inverse(x)
        if x is None:return 'reject',None
    if x[0]<1 or x[0]%2!=1 or (x[1]-x[2])%2:return 'reject',None
    return 'recovered',(packet,x)

@lru_cache(None)
def g(p):
    if not p:return 1
    v=p//2
    return g(v) if p%2 else g(v)+v+v%2

def block(w):
    p=int(w,2);q=sum(1<<j for j,b in enumerate(w) if b=='0')
    return ((g(p),q+1,-int(q%2==0)),(p+1,1<<len(w),0),(int(p%2==0),0,0))
def packet(words):
    x=(1,0,0)
    for w in words:x=mv(block(w),x)
    return x

def rows_for(words):
    rows=[];roots=[]
    for word in words:
        z=len(rows); roots.append(z);n=len(word)+1;rows.extend([0]*n)
        for a in range(z):
            for b in range(z,z+n):rows[a]|=1<<b
        for i,b in enumerate(word,1):
            if b=='1':rows[z]|=1<<(z+i)
            else:rows[z+i]|=1<<z
            for j in range(i+1,n):rows[z+i]|=1<<(z+j)
    for a,b in zip(roots,roots[1:]):rows[a]^=1<<b;rows[b]|=1<<a
    return rows,roots[-1]

def literal_packet(words):
    rows,root=rows_for(words);n=len(rows);full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    rp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    rp[1<<root][root]=1
    for s in range(1,1<<n):
        for v in range(n):
            todo=rows[v]&~s
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                dp[s|bit][j]+=dp[s][v];rp[s|bit][j]+=rp[s][v]
    h=[sum(r) for r in dp];h[0]=1;o=e=0
    for s in range(1,1<<n):
        term=sum(rp[s])*h[full^s]
        if s.bit_count()%2:o+=term
        else:e+=term
    return h[full],o+e,o-e

def main():
    assert T==((3,3,4),(-3,1,0),(-4,0,0))
    rays3=geometry(3)[0]
    lower=min(F(mv(M,r)[0],r[0]) for r in rays3)
    upper=max(F(mv(M,r)[0],r[0]) for r in rays3)
    assert (lower,upper)==(F(97,15),F(797,123))
    norm4=max(sum(abs(x) for x in row) for row in power(T,4))
    assert norm4==992
    rho=norm4*(upper/lower**2)**4
    assert rho<1
    B=96; bounds={}
    for k in range(3,101):
        _,lo,hi,amin,amax,E=geometry(k)
        U=B*F(amax*E,amin**2);bounds[k]=U
        assert all(B*amax*(h-l)<=U for l,h in zip(lo,hi))
        if k>=7:assert U<=rho*bounds[k-4]
    assert all(bounds[k]<1 for k in range(25,29))
    words=[''.join(w) for n in range(1,4) for w in product('01',repeat=n)]
    bank={};literal_count=0
    for length in (1,2):
        for ws in product(words,repeat=length):
            x=packet(ws)
            assert x==literal_packet(ws),(ws,x,literal_packet(ws))
            assert cone(x);bank.setdefault(x,ws);literal_count+=1
    assert len(bank)==84 and max(s[1] for s in bank)==B
    recoveries=0;ambiguous=0; scalar_unique={}
    for k in range(1,41):
        seen={}
        for seed in bank:
            final=mv(power(M,k),seed)
            status,result=recover_packet(final[0],k)
            assert status!='reject'
            if status=='recovered':
                assert result==(final,seed);recoveries+=1
            else:ambiguous+=1
            seen.setdefault(final[0],[]).append(seed)
        scalar_unique[k]=len(seen)
        if k>=25:assert len(seen)==len(bank)
    # Literal counts on both missing-memory controls and their short actual suffixes.
    pair=(('01','00'),('10','0'))
    assert tuple(packet(w) for w in pair)==((9,9,1),(9,9,3))
    pair_rows=[]
    for k in (0,1,2,5,10,25,100):
        states=[mv(power(M,k),packet(w)) for w in pair]
        if k>=1:assert states[0][0]-states[1][0]==2*mv(power(M,k-1),(1,0,0))[0]
        if k<=2:
            for w,x in zip(pair,states):assert literal_packet(w+('10',)*k)==x
        pair_rows.append({'suffix_blocks':k,'counts':[str(x[0]) for x in states],'difference':str(states[0][0]-states[1][0])})
    k=25;outputs={mv(power(M,k),s)[0]:s for s in bank}
    trials={n+d for n in outputs for d in (-2,0,2)}
    accepted=[];rejected=[]
    for n in trials:
        status,result=recover_packet(n,k)
        assert status!='ambiguous'
        ok=status=='recovered' and result[1] in bank
        assert ok==(n in outputs)
        (accepted if ok else rejected).append(n)
    synchronizers=[]
    for n in range(2,13):
        passing=0;bad=[]
        for mid in product('01',repeat=n-2):
            w='1'+''.join(mid)+'0';p=int(w,2)
            q=sum(1<<j for j,b in enumerate(w) if b=='0');t=1<<n;h=g(p)
            coefficient=h*t-(p+1)*(q+1)+1
            fplus=1-(h+t)+coefficient-t
            fminus=-1-(h+t)-coefficient-t
            if fplus<0 and fminus<0:passing+=1
            else:bad.append({'word':w,'f(1)':fplus,'f(-1)':fminus})
            mat=block(w)
            assert sum(mat[i][i] for i in range(3))==h+t
            minors=sum(mat[i][i]*mat[j][j]-mat[i][j]*mat[j][i] for i,j in PAIRS)
            assert minors==coefficient
            determinant=sum(mat[0][j]*(mat[1][(j+1)%3]*mat[2][(j+2)%3]-mat[1][(j+2)%3]*mat[2][(j+1)%3]) for j in range(3))
            assert determinant==t
        if n<=7:assert not bad
        if n==8:assert bad==[{'word':'10101010','f(1)':84,'f(-1)':-1342}]
        synchronizers.append({'word_length':n,'Pisot':passing,'not_Pisot':len(bad),'first_failures':bad[:3]})
    print(json.dumps({
      'physical_seed_words_literal_checked':literal_count,
      'distinct_root_packets':len(bank),'maximum_seed_W':B,
      'uniform_threshold':25,'valid_for':'every integer cone seed with 1<=A and W<=96, at every k>=25',
      'growth_bounds_C3':[str(lower),str(upper)],'wedge_four_step_norm':norm4,'rho_exact':str(rho),
      'four_threshold_bounds':{str(k):{'exact':str(bounds[k]),'decimal_display':float(bounds[k])} for k in range(25,29)},
      'successful_packet_recoveries_k1_to40':recoveries,'sound_ambiguities':ambiguous,
      'scalar_distinct_counts_by_k':scalar_unique,
      'depth25_near_target_trials':len(trials),'accepted_in_seed_bank':len(accepted),'rejected_from_seed_bank':len(rejected),
      'missing_memory_control':pair_rows,
      'active_word_spectral_classification':synchronizers,
      'scope':'Exact suffix decoder and seed-bank selector; not unrestricted tournament-spectrum coverage.'
    },indent=2))
if __name__=='__main__':main()

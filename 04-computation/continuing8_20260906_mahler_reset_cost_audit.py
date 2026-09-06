"""Independent native-first referee for the Mahler fixed-prefix reset theorem.

No producer imports. Addresses are built one native bit at a time, not by
inverting a ternary denominator modulo a binary power. Output is raw LF.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0
def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(why)

here=Path(__file__).resolve().parent
data=here.parent/"05-knowledge/results" if here.name=="04-computation" else here
stem="continuing8_20260906_mahler_reset_cost"
primary=here/(stem+".py")
certificate=data/(stem+"_certificate.json")
need(sha256(primary.read_bytes()).hexdigest()=="77757ed311aac5246cf3a420dea282e84cb7aaaf428abcc928c633e8865f5646","frozen producer source")
need(sha256(certificate.read_bytes()).hexdigest()=="ac3f7f74a5ed80351218669c628a4bf068384f7e8132fa41af0b0828b6646d6c","frozen producer certificate")
saved=json.loads(certificate.read_bytes())

def ceil_orbit(A,n):
    bits=[]
    for _ in range(n):
        bits.append(A&1)
        A=(3*A+1)//2
    return bits,A

def lift(word):
    """Maintain A in [0,2^i), U=T^i(A); choose the next native bit."""
    A,U=0,0
    two,three=1,1
    for c in word:
        native=(c-U)&1
        A+=native*two
        U+=native*three
        need(U%2==c,"native digit makes the desired literal carry")
        U=(3*U+c)//2
        two*=2; three*=3
    return A,U

def oriented_phase(word):
    tail=F(0)
    for c in reversed(word): tail=(F(c)+2*tail)/3
    return tail

def safe_suffixes(word):
    tail=F(0)
    for c in reversed(word):
        tail=(F(c)+2*tail)/3
        if tail>=F(1,2): return False
    return True

def greedy_digits(n):
    num,den=1,1
    digits=[]
    for _ in range(n):
        num*=3; den*=2
        c,num=divmod(num,den)
        digits.append(c)
    return digits

def follower(word):
    d=greedy_digits(len(word)+1)
    q=0
    for i,c in enumerate(word):
        if c>d[q]: return None,i
        q=0 if c<d[q] else q+1
    return q,None

def twoval(v):
    return (v&-v).bit_length()-1

def floor_ratio(N,c):
    k=N.bit_length()-c.bit_length()
    if c*(1<<k)>N: k-=1
    return k

# A complete native-first universe makes the bijection and minimum observable.
native_rows=0; safe_rows=0
for n in range(1,12):
    seen=set()
    for A in range(1<<n):
        w,U=ceil_orbit(A,n)
        need(tuple(w) not in seen,"literal level map injective")
        seen.add(tuple(w)); native_rows+=1
        phase=oriented_phase(w)
        need(F(3**n*A,2**n)+phase*F(3**n,2**n)==U,"dense oriented affine carrier")
        if A and safe_suffixes(w):
            safe_rows+=1
            D,N=phase.numerator,phase.denominator
            j=max(i for i,c in enumerate(w) if c)+1
            need(N==3**j and 0<2*D<N and gcd(D,3)==1,"exact last-one denominator and strict tail")
            b=(N*A+D)//(1<<n)
            need(N*A+D==(1<<n)*b and 0<b<N,"native-first launch formula")
            r=twoval(U); c=b>>r
            need(n+r-A.bit_length()==floor_ratio(N,c),"native clock formula in full native universe")
        need(safe_suffixes(w)==(follower(w)[1] is None),"independent suffix versus follower safety")
    need(len(seen)==1<<n,"complete native level bijection")
print("NATIVE_UNIVERSE",native_rows,"rows through depth11; strict nonzero rows",safe_rows)

cycles=0
for record in saved["prefix_cycles"]:
    word=list(map(int,record["prefix"]))
    m=len(word); phase=oriented_phase(word); D,N=phase.numerator,phase.denominator
    A,U=lift(word)
    period=2*N//3; R=N.bit_length()-1
    need(follower(word)==(0,None),"declared reset endpoint")
    units=[]; histogram=Counter(); active=0
    for offset in range(period):
        n=m+offset
        b=(N*A+D)//(1<<n)
        need(U==3**(n-max(i for i,c in enumerate(word) if c)-1)*b,"native affine exit")
        need(0<A<(1<<n) and (N*A+D)==(1<<n)*b,"least native representative")
        r=twoval(U); c=b>>r; H=A.bit_length(); K=n+r-H
        need(K==floor_ratio(N,c) and 0<=K<=R,"clock and universal bound")
        if H>=m:
            terminal=ceil_orbit(A,H)[1]
            carried,end=ceil_orbit(terminal,K)
            need(carried==[0]*K and end%2==1,"actual terminal runway, not observation runway")
        if n<=14:
            for extra in (1,3,7):
                bits,end=ceil_orbit(A+(1<<n)*extra,n)
                need(bits==word+[0]*offset and end==U+3**n*extra,"all-completion affine controls")
        units.append(b); histogram[K]+=1
        new_bit=U%2
        active+=new_bit
        A+=new_bit*(1<<n)
        U=(3*(U+new_bit*3**n))//2
    need(units==record["unit_cycle"],"every frozen unit-cycle entry by native lifting")
    need(sorted(units)==[x for x in range(1,N) if x%3],"complete unit period")
    need(sorted(histogram.items())==[tuple(v) for v in record["lag_histogram"]],"every saved native lag frequency")
    need(active==period//2,"half native increments on declared phase measure")
    for k in range(R+1):
        # Direct residue enumeration, independent of producer floor-count F(M).
        odds=[c for c in range(1,N,2) if c%3 and (1<<k)*c<N<(1<<(k+1))*c]
        need(histogram[k]==(k+1)*len(odds),"odd-part lag multiplicities")
        need(abs(F(len(odds))-F(N,3*(1<<(k+1))))<=2,"uniform interval discrepancy")
    need(histogram[R]==R+1,"maximum attained at exactly the powers of two")
    tv=(sum(abs(F(histogram[k],period)-F(k+1,1<<(k+2))) for k in range(R+1))+F(R+3,1<<(R+2)))/2
    need(tv==F(record["tv_to_two_head_law"]) and tv<=F((R+2)**2,N),"exact probability normalization and TV bound")
    cycles+=1
    print("CYCLE",record["prefix"],"N",N,"lag",sorted(histogram.items()),"TV",str(tv))

# Fresh phase-distribution controls beyond all saved denominator cycles.
for j in range(1,9):
    N=3**j; R=N.bit_length()-1; P=2*N//3
    hist=Counter(floor_ratio(N,b>>twoval(b)) for b in range(1,N) if b%3)
    need(sum(hist.values())==P,"independent full-residue probability space")
    need(hist[R]==R+1,"all-order sampled sharp lag count")
    tv=(sum(abs(F(hist[k],P)-F(k+1,1<<(k+2))) for k in range(R+1))+F(R+3,1<<(R+2)))/2
    need(tv<=F((R+2)**2,N),"fresh exact fair-coin convergence bound")
need(Counter({0:2,3:4})==Counter(floor_ratio(9,b>>twoval(b)) for b in range(1,9) if b%3),"N9 missing intermediate lags")

family_rows=0
for j in range(1,8):
    N=3**j; P=2*N//3; R=N.bit_length()-1
    for t in (2,3,5):
        E=P*t
        numerator=(1<<(j-1))*((1<<E)-1)
        need(numerator%N==0,"sharp family integral launch")
        A=numerator//N; H=A.bit_length()
        need(H==E+j-1-R and H>=j+2,"sharp family exact native height")
        word,end=ceil_orbit(A,E+j+1)
        need(word==[0]*(j-1)+[1]+[0]*(E-1)+[1,1],"complete independently generated sharp carry word")
        need(follower(word)[1]==H+R+1,"first actual rejection index")
        bits,odd=ceil_orbit(ceil_orbit(A,H)[1],R)
        need(bits==[0]*R and odd==3**(E-1),"maximal postterminal runway equality")
        family_rows+=1
print("SHARP_FAMILY",family_rows,"fresh literal rows at j1..7,t2,3,5")

# Reconstruct the real endpoint independently by repeated rational multiplication.
for record in saved["prefix_cycles"]:
    w=list(map(int,record["prefix"])); phase=oriented_phase(w)
    for extra in (0,1,4,9):
        n=len(w)+extra; A,U=lift(w+[0]*extra); r=twoval(U)
        alpha=F(A)+phase; first=None
        for time in range(n+r+2):
            part=alpha-(alpha.numerator//alpha.denominator)
            if part>=F(1,2): first=(time,part); break
            alpha*=F(3,2)
        need(first==(n+r+1,F(1,2)),"exact first strict real endpoint")
        bits,_=ceil_orbit(A,n+r+1)
        need(follower(bits)[1] is None,"real endpoint alone does not reject ordinary carry prefix")

# Matrix multiplication is a different reconstruction of the inverse observer.
maps=((0,2,3,3),(2,0,1,1))
def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def observer(word):
    M=[[int(i==j) for j in range(4)] for i in range(4)]
    for c in word:
        G=[[int(maps[c][i]==j) for j in range(4)] for i in range(4)]
        M=matmul(G,M)
    return M
constant=[[0,0,0,1] for _ in range(4)]
for word in ("00100","0010011","00100000","00100100"):
    need(observer(list(map(int,word)))==constant,"independent saturated observer matrix")
need(safe_suffixes(list(map(int,"00100"))) and not safe_suffixes(list(map(int,"0010011"))),"finite predicate is not a function of saturated observer")
for text,A,U,N,bit in saved["joint_quotient_hostile"]:
    w,end=ceil_orbit(A,8)
    need(w==list(map(int,text)) and end==U and end%2==bit,"same-clock native hostile")
    need(A.bit_length()==8 and follower(w)==(0,None) and oriented_phase(w).denominator==N,"joint hostile retains clocks and loses arithmetic phase")

# The t>=2 and H>=m qualifications are semantic, not decorative.
need(ceil_orbit(1,4)[0]==[1,0,1,1] and follower([1])==(1,None),"t1 has a terminal match rather than a reset state")
print("QUOTIENT_AND_ENDPOINT controls PASS; no fixed-orbit randomness or infinite survivor inferred")
print("PASS",GATES,"always-active independent exact gates; raw LF")

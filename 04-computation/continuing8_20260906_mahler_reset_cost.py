"""Exact fixed-prefix reset cost, native clock, and oriented real endpoint."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(newline='\n')
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def val(a,p):
    if a==0: raise ValueError('zero valuation')
    k=0
    while a%p==0: a//=p; k+=1
    return k
def orbit(A,n):
    out=[]
    for _ in range(n):
        c=A%2; out.append(c); A=(3*A+c)//2
    return tuple(out),A
def greedy(n):
    z=F(1); d=[]
    for _ in range(n):
        z=3*z/2; a=z.numerator//z.denominator; d.append(a); z-=a
    return d
def follower(word):
    d=greedy(len(word)+1); q=0; resets=0
    for i,c in enumerate(word):
        if c>d[q]: return None,resets,i
        if c<d[q]: q=0; resets+=1
        else: q+=1
    return q,resets,None
def finite_safe(word):
    # All exact oriented truncated suffix values, computed from the right.
    x=F(0)
    for c in reversed(word):
        x=F(2,3)*(c+x)
        if x>=1: return False
    return True
def carrier(word):
    m=len(word)
    C=sum(c*2**i*3**(m-1-i) for i,c in enumerate(word))
    phase=F(C,3**m)
    D,N=phase.numerator,phase.denominator
    residue=(-C*pow(3**m,-1,2**m))%(2**m)
    U=(3**m*residue+C)//2**m
    return C,D,N,residue,U
def cost(word,n):
    C,D,N,r,U=carrier(word)
    b=(D*pow(2,-n,N))%N
    A=(2**n*b-D)//N
    return A,b
def floor_log_ratio(N,c):
    k=0
    while (1<<(k+1))*c<N: k+=1
    return k
def odd_units_up_to(M):
    return (M+1)//2-(M//3+1)//2
def real_first_bad(alpha,maximum):
    for i in range(maximum+1):
        frac=alpha-alpha.numerator//alpha.denominator
        if frac>=F(1,2): return i,frac
        alpha*=F(3,2)
    return None

records=[]
prefixes=['100','00100','00100100']+['0'*(j-1)+'100' for j in range(2,7)]
for text in dict.fromkeys(prefixes):
    w=tuple(map(int,text)); m=len(w)
    C,D,N,r,U=carrier(w); j=max(i for i,c in enumerate(w) if c)+1
    P=2*3**(j-1); R=(N-1).bit_length()-1
    gate(follower(w)[0]==0 and finite_safe(w),'declared prefix is strict and ends at reset state')
    gate(N==3**j and gcd(D,3)==1 and 0<D<N and F(D,N)<F(1,2),'last-carry-one address gives exact reduced real denominator')
    gate(val(U,3)==m-j and 0<U<3**m,'retained affine terminal sidecar valuation')
    units=[]; delays=[]; native_bits=[]; examples=[]
    for n in range(m,m+P):
        A,b=cost(w,n)
        word=w+(0,)*(n-m)
        carried,end=orbit(A,n)
        gate(carried==word and end==3**(n-j)*b,'literal orbit versus complete affine cylinder')
        gate(0<A<2**n and (N*A+D)==2**n*b,'least residue and exact cost equation')
        gate(A==(-D*pow(N,-1,2**n))%(2**n),'independent native modular address')
        if n<=11:
            brute=[B for B in range(1,A+1) if orbit(B,n)[0]==word]
            gate(brute==[A],'literal least positive launch in bounded control')
        for k in (0,1,2,5):
            carried2,end2=orbit(A+2**n*k,n)
            gate(carried2==word and end2==3**(n-j)*b+3**n*k,'all-completion affine cone controls')
        rr=val(b,2); odd=b>>rr; h=A.bit_length()
        gate(n+rr-h==floor_log_ratio(N,odd),'exact termination-to-odd-state clock')
        gate(n+rr-h<=R and rr<=R,'sharp denominator runway ceiling')
        nxt,_=cost(w,n+1)
        gate(nxt==A+2**n*(b%2),'actual next native activity bit')
        later,later_b=cost(w,n+P)
        gate(later_b==b,'one full native phase period')
        later_h=later.bit_length()
        gate(later_h>=m and n+P+rr-later_h==n+rr-h,'full period gives the same actual postterminal reset count')
        gate(2**n<=N*A+D and n<=(N*A+D).bit_length()-1,'exact all-launch height ceiling')
        if b==1: gate(n==(N*A+D).bit_length()-1,'sharp lower launch cost attained')
        continued,oddstate=orbit(end,rr)
        gate(continued==(0,)*rr and oddstate%2==1,'exact additional zero carries at observation')
        if h>=m:
            terminal=orbit(A,h)[1]
            gap=n+rr-h
            z,at_odd=orbit(terminal,gap)
            gate(z==(0,)*gap and at_odd%2==1,'actual post-native-termination zero reset runway')
        if n<=45 or n in (m,m+P-1):
            bad=real_first_bad(F(A)+F(D,N),n+rr+1)
            gate(bad==(n+rr+1,F(1,2)),'literal original phase hits exact first half endpoint')
        units.append(b); delays.append(n+rr-h); native_bits.append(b%2)
        if n<m+3 or b in (1,1<<R):
            examples.append([n,str(A),b,h,rr,n+rr-h])
    gate(sorted(units)==[b for b in range(1,N) if b%3],'exact complete unit cost spectrum and period')
    gate(sum(native_bits)==P//2,'half native ones per exact phase cycle')
    gate(max(delays)==R,'sharp observed clock-delay bound in a complete period')
    histogram=Counter(delays)
    expected={k:(k+1)*(odd_units_up_to((N-1)//2**k)-odd_units_up_to(N//2**(k+1))) for k in range(R+1)}
    gate(all(histogram[k]==expected[k] for k in range(R+1)),'exact clock-lag frequency from primitive odd parts')
    gate(histogram[R]==R+1,'exact multiplicity of maximal postterminal runway')
    tv=(sum(abs(F(histogram[k],P)-F(k+1,2**(k+2))) for k in range(R+1))+F(R+3,2**(R+2)))/2
    gate(tv<=F((R+2)**2,N),'exact total variation control toward the two-head fair-coin failure law')
    print('PREFIX',text,'D/N',str(F(D,N)),'period',P,'max native-to-odd lag',R,'native ones',sum(native_bits))
    records.append({'prefix':text,'D':D,'N':N,'period':P,'lag_maximum':R,
                    'unit_cycle':units,'lag_histogram':sorted(histogram.items()),'tv_to_two_head_law':str(tv),'examples':examples})

# An explicit all-parameter sharp family; the proof supplies the universal identities.
family=[]
for j in range(1,8):
    N=3**j; P=2*3**(j-1); R=(N-1).bit_length()-1
    for t in (2,3,4):
        E=P*t; n=E+j-1-R
        A=2**(j-1)*(2**E-1)//N
        gate(2**(j-1)*(2**E-1)%N==0,'sharp family is an ordinary integer')
        gate(A.bit_length()==n and n>=j+2,'sharp family exact actual terminal depth')
        wanted=(0,)*(j-1)+(1,)+(0,)*(E-1)+(1,1)
        c,last=orbit(A,len(wanted))
        gate(c==wanted,'literal entire sharp-family carry word through first rejection')
        # After each initial reset, only the short states0,1,2 occur until the final11.
        q=0; first_bad=None
        for i,bit in enumerate(c):
            digit=(1,0,1)[q]
            if bit>digit: first_bad=i; break
            q=0 if bit<digit else q+1
        gate(first_bad==n+R+1,'exact first follower rejection, not merely a rational pseudo-orbit')
        z,u=orbit(orbit(A,n)[1],R)
        gate(z==(0,)*R and u==3**(E-1),'sharp unbounded postterminal reset runway')
        gate(N>2**R and N<2**(R+1),'sharp integer logarithmic sidecar')
        if j<=4:
            gate(real_first_bad(F(A)+F(2**(j-1),N),n+R+1)==(n+R+1,F(1,2)),'sharp family real boundary independently followed')
        family.append([j,t,E,n,R,str(A),first_bad])
    print('SHARP FAMILY j',j,'N',N,'postterminal reset count',R,'first rejection depth h+R+1')

# Literal right-composition of the inherited Rule30 inverse-boundary generators.
maps=((0,2,3,3),(2,0,1,1))
def monoid(word):
    f=tuple(range(4))
    for c in word: f=tuple(f[x] for x in maps[c])
    return f
safe=tuple(map(int,'00100')); unsafe=safe+(1,1)
gate(monoid(safe)==monoid(unsafe)==(3,3,3,3),'actual Rule30 saturated state aliases safe and unsafe Mahler words')
gate(finite_safe(safe) and not finite_safe(unsafe),'strict Mahler finite predicate does not factor through spatial monoid')
joint=[]
for text,A,U,nextbit in (('00100000',180,4617,1),('00100100',148,3798,0)):
    w=tuple(map(int,text)); C,D,N,r,u=carrier(w)
    gate((r,u)==(A,U) and orbit(A,8)==(w,U),'literal same-depth terminal joint hostile')
    gate(A.bit_length()==8 and follower(w)[0]==0 and monoid(w)==(3,3,3,3),'same native terminal depth/follower/reset monoid')
    gate(U%2==nextbit,'opposite next reset versus match')
    joint.append([text,A,U,N,nextbit])
gate(joint[0][3:]==[27,1] and joint[1][3:]==[729,0],'lost reduced denominator and immediate event both differ')

# Inherited ordinary-terminal hostile retained as a separate control.
for A,w,rejection in ((8,'0001',4),(13,'1001',10)):
    c,u=orbit(A,12)
    gate(c[:4]==tuple(map(int,w)) and follower(c)[2]==rejection,'inherited A8/A13 future-decision hostile')
gate(follower(tuple(map(int,'1011')))[2]==3,'native-one terminal unsafe control')
print('JOINT QUOTIENT HOSTILE',joint)
print('SCOPE finite zero-tail cylinders and exact ordinary launch families; no Mahler infinite survivor or physical Rule30 conclusion')
certificate={'prefix_cycles':records,'sharp_family_controls':family,'joint_quotient_hostile':joint,
             'rule30_constant':[3,3,3,3],'safe_word':'00100','unsafe_word':'0010011'}
blob=(json.dumps(certificate,separators=(',',':'),sort_keys=True)+'\n').encode()
here=Path(__file__).resolve().parent
out=here.parent/'05-knowledge/results' if here.name=='04-computation' else here
(out/(Path(__file__).stem+'_certificate.json')).write_bytes(blob)
print('CERTIFICATE SHA256',sha256(blob).hexdigest())
print('PASS',G,'always-active exact gates; raw LF')

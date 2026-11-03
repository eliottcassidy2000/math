"""Actual return matrices, exact phase-preserving legalizers, and the 49/50/51 completion."""
from collections import deque,Counter
from itertools import product
from pathlib import Path
from math import gcd
import json
from observer import HERE,need

M3=(61,-2,9,0);M5=(421,-2,33,0);ID=(1,0,0,1)
def mul(A,B,m):
    a,b,c,d=A;e,f,g,h=B
    return ((a*e+b*g)%m,(a*f+b*h)%m,(c*e+d*g)%m,(c*f+d*h)%m)
def power(A,n,m):
    out=ID
    while n:
        if n&1:out=mul(A,out,m)
        A=mul(A,A,m);n//=2
    return out
def mv(A,v,m):return ((A[0]*v[0]+A[1]*v[1])%m,(A[2]*v[0]+A[3]*v[1])%m)
def canonical(A):
    f=next(v for v in A if v);a=pow(f,-1,7);return tuple(a*x%7 for x in A)

# F49=F7[theta], theta^2=4theta+2. Integers0,...,48 encode a+7b.
def fadd(x,y):return ((x%7+y%7)%7)+7*((x//7+y//7)%7)
def neg(x):return (-x%7)+7*(-(x//7)%7)
def fmul(x,y):
    a,b=x%7,x//7;c,d=y%7,y//7
    return (a*c+2*b*d)%7+7*((a*d+b*c+4*b*d)%7)
def fpow(x,n):
    out=1
    while n:
        if n&1:out=fmul(out,x)
        x=fmul(x,x);n//=2
    return out
def divide(x,y):need(y!=0,'field denominator');return fmul(x,fpow(y,47))
def act(A,r):
    if r==50:return 50
    if r==49:num=A[3];den=A[1]
    else:num=fadd(A[2],fmul(A[3],r));den=fadd(A[0],fmul(A[1],r))
    return divide(num,den) if den else 49 if num else 50
def alpha(v):
    b,t=v;return (5*b+t)%7+7*(t%7)
def pair_point(v,w):
    x=v[0]%7+7*(w[0]%7);y=v[1]%7+7*(w[1]%7)
    return divide(y,x) if x else 49 if y else 50
def ratio11(v):
    b,t=v;return t*pow(b,-1,11)%11 if b%11 else 11
def r3(r):return act11(r,6,9)
def act11(r,a,c):return 0 if r==11 else c*pow((a-2*r)%11,-1,11)%11 if (a-2*r)%11 else 11

def main():
    A=tuple(x%7 for x in M3);B=tuple(x%7 for x in M5)
    words={ID:''};queue=[ID]
    for T in queue:
        for letter,G in (('3',A),('5',B)):
            U=mul(G,T,7)
            if U not in words:words[U]=words[T]+letter;queue.append(U)
    expected={T for T in product(range(7),repeat=4) if (T[0]*T[3]-T[1]*T[2])%7}
    need(set(words)==expected and len(words)==2016,'full GL2(F7) return group')
    targets=[(1,1,0,1),(1,0,1,1),(3,0,0,1)]
    proof_words=[]
    for target in targets:
        word=words[target];out=ID
        for letter in word:out=mul(A if letter=='3' else B,out,7)
        need(out==target,'explicit elementary generator word');proof_words.append({'target':target,'word':word})
    group=sorted({canonical(T) for T in words});need(len(group)==336,'projective group order')
    rational=set(range(7))|{49};nonrational=set(range(50))-rational
    need({act(T,49) for T in group}==rational and {act(T,7) for T in group}==nonrational,'eight and forty-two projective orbits')
    signatures0={};signatures1={}
    for r in range(51):
        coarse=tuple(act(T,r) in (49,50) for T in group)
        fine=tuple((act(T,r) in (49,50),act(T,r) in (7,50)) for T in group)
        signatures0.setdefault(coarse,[]).append(r);signatures1.setdefault(fine,[]).append(r)
    need(sorted(map(len,signatures0.values()))==[1]*9+[42] and len(signatures1)==51,'minimal 10 versus 51 zero-test carriers')
    theta=7;beta=fadd(theta,2)
    logarithms={fpow(beta,j):j for j in range(48)}
    need(len(logarithms)==48 and fpow(beta,8)==3 and fpow(beta,7)==fadd(1,neg(beta)),'primitive beta norm and Frobenius')
    for v in product(range(7),repeat=2):need(alpha(mv(B,v,7))==fmul(beta,alpha(v)),'actual C5 phase encoding')
    need(power(A,3,7)==ID and power(A,1,7)!=ID,'uniform phase-neutral C3 clock')
    # Common guard clock for one, two or three actual source alternatives.
    controller=(0,1,3,4,5,6,7,11);legalizers=[]
    for precision in range(1,5):
        modulus=7**precision;period=3*7**(precision-1)
        need(power(M3,period,modulus)==ID,'lifted integral phase period')
        if precision>1:need(power(M3,period//7,modulus)!=ID,'exact lifted period')
        need(power(M5,48*7**(precision-1),modulus)==ID,'lifted C5 phase clock')
        if precision>1:need(power(M5,48*7**(precision-2),modulus)!=ID,'exact C5 lift')
        checked=0;worst=0
        for arity in range(1,4):
            for states in product(controller,repeat=arity):
                success=[]
                for j in range(arity+1):
                    phases=list(states)
                    for _ in range((j*period)%4):phases=[r3(r) for r in phases]
                    if all(r!=7 for r in phases):success.append(j)
                need(success,'shared guarded phase-preserving C5 lift')
                worst=max(worst,min(success));checked+=1
        legalizers.append({'precision':precision,'modulus':modulus,'phase_period_C3':period,'C5_period':48*7**(precision-1),'tuples_checked':checked,'maximum_extra_blocks':worst})
    hostile=(0,7,3,11)
    need(7 in hostile,'four-phase hostile includes forbidden state')
    phases=list(hostile)
    for j in range(4):
        need(7 in phases,'four alternatives cover every guard clock');phases=[r3(r) for r in phases]
    phases=[]
    for k in range(1,31):
        b,t=(1659*2**k+93,553*2**k-45) if k%2==0 else (1401*2**k+157,467*2**k+35)
        v=(b,t);e=logarithms[alpha(v)]
        need(gcd(b,t)==1 and (3*b-2*t)%11!=0,'selected return seed')
        if k<=6:phases.append({'k_mod6':k%6,'pair_mod7':[b%7,t%7],'alpha':alpha(v),'beta_log':e,'Frobenius_clock_shift':6*e%48})
    # Actual independent pair from k2 and k3; the shared C5 entry keeps it independent.
    v=(6729,2167);w=(11365,3771);v=mv(M5,v,7);w=mv(M5,w,7)
    need((v[0]*w[1]-w[0]*v[1])%7!=0,'actual independent source pair')
    point=pair_point(v,w)
    need({act(T,point) for T in group}==nonrational,'all forty-two directions from the supplied pair')
    # Frobenius commutes with every scalar-extension return matrix.
    for T in (A,B):
        for r in range(50):
            fr=fpow(r,7) if r<49 else 49
            ar=act(T,r);far=fpow(ar,7) if ar<49 else 49
            need(act(T,fr)==far,'Frobenius-compatible projective action')
    lifts=[]
    for matrix,period in ((M3,3),(M5,48)):
        residue=power(matrix,period,49)
        need(all((a-b)%7==0 for a,b in zip(residue,ID)),'first integral lift')
        leading=tuple(((a-b)//7)%7 for a,b in zip(residue,ID))
        need(any(leading),'nonzero first lift implies exact all-order odd-prime clock')
        lifts.append({'period_mod7':period,'leading_mod7':leading})
    need(mul(tuple(lifts[0]['leading_mod7']),tuple(lifts[0]['leading_mod7']),7)==ID,'involutory first C3 lift')
    out={'status':'PROVED_SCOPED','M3':M3,'M5':M5,'GL2_size':2016,'PGL2_size':336,'integral_lift_certificates':lifts,
         'elementary_generator_words':proof_words,'GL2_words':[{'matrix':T,'word':word} for T,word in sorted(words.items())],
         'coarse_classes':list(signatures0.values()),'fine_class_count':len(signatures1),'theta':theta,'beta':beta,
         'legalizers':legalizers,'four_source_hostile':hostile,'H497_seed_phases':phases,
         'actual_paired_seed_mod7':[v,w],'actual_paired_seed_projective_point':point}
    (HERE/'field_clock.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: GL2(F7), 48-step sourced clock, shared legalizers, and exact 10-to-51 selector refinement')
if __name__=='__main__':main()

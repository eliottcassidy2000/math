"""Exact debt-cone controls; no universal entry or first-exit claim."""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import json

def need(ok, message):
    if not ok: raise RuntimeError(message)

def coefficients(w):
    A,d,B,C,D=1,1,0,0,0
    for a in w:
        A,B,C,D,d=(3*A,3*B+d,3*C+((d<<a) if a==1 else 0),
                    3*D+((1<<a)-4)*d if a>=3 else 3*D,d<<a)
    need(A+B-d==C-D,'shifted decomposition')
    return A,d,B,C,D

def direct(w):
    L,P,D=F(1),F(0),F(0)
    for a in reversed(w):
        P+=(a==1)*L
        if a>=3:D+=(1-F(4,1<<a))*L
        L*=F(3,1<<a)
    return L,P,D

def exit_word(n):
    x=n; w=[]
    for _ in range(10000):
        y=3*x+1; a=(y&-y).bit_length()-1; x=y>>a; w.append(a)
        if x<n:return w,x
    raise RuntimeError(('unresolved source under explicit cap',n))

def replay(n,w):
    x=n; vals=[x]
    for a in w:
        y=3*x+1; need((y&-y).bit_length()-1==a,'exact integer clock')
        x=y>>a; vals.append(x)
    return vals

def main():
    checks=0; bank=[]
    for m in range(1,8):
        for w in product(range(1,5),repeat=m):
            A,d,B,C,D=coefficients(w); L,P,E=direct(w)
            need((L,P,E)==(F(A,d),F(C,d),F(D,d)),'independent suffix sums')
            resets=[i for i,a in enumerate(w) if a>=3]
            if resets:
                need(E.denominator==1<<(sum(w[resets[0]:])-2),
                     'exact reduced debt denominator retains first reset clock')
            else:need(E==0,'no reset has zero debt')
            checks+=1
            if A>=d and D<=A-d and D and len(bank)<32:bank.append(w)
    for u in bank:
        for v in bank:
            A,d,B,C,D=coefficients(u+v)
            need(D<=A-d,'cone closed under changing generators')
    base=(1,1,1,1,3); A,d,_,_,D=coefficients(base)
    need((A,d,D)==(243,128,64),'block constants')
    repetitions=[]
    for k in (1,2,3,4,8,16,32,64):
        # The smallest terminal valuation making the whole word contractive.
        a=max(2,(3**(5*k+1)).bit_length()-7*k)
        w=base*k+(a,); A,d,B,C,D=coefficients(w)
        need(A<d,'terminal contraction')
        r=((d-B)*pow(A,-1,2*d))%(2*d)
        n=r+2*d
        if (d-A)*n<=B:n+=2*d*((B-(d-A)*n)//((d-A)*2*d)+1)
        vals=replay(n,w); h=vals[-1]
        need(min(vals[:-1])==n and h<n,'first descent in displayed family')
        need(F(D,d)<F(179,115)<2<=n-h,'uniform 179/115 debt bound')
        need(n-h-F(D,d)>F(51,115),'depth-independent payment margin')
        need(A*(n-1)+C<d*(n-1),'consequence: positive charge paid')
        repetitions.append({'k':k,'odd_steps':len(w),'K':sum(w),'last_a':a,
                            'n':n,'h':h,'D':str(F(D,d))})
    stats=Counter(); hostile=[]
    for n in range(3,100002,2):
        w,h=exit_word(n); A,d,B,C,D=coefficients(w)
        need(F(D,d)<n-h,'finite positive-charge control')
        L,s,_,_,E=coefficients(w[:-1])
        admitted=E<=L-s
        stats['sources']+=1;stats['cone_admitted' if admitted else 'cone_unresolved']+=1
        if sum(a>=3 for a in w)>=3:stats['three_or_more_resets']+=1;stats['three_reset_cone_admitted']+=admitted
        if not admitted and len(hostile)<8:
            hostile.append({'n':n,'h':h,'word':w,'prefix_multiplier':str(F(L,s)),
                            'prefix_debt':str(F(E,s)),'D':str(F(D,d))})
    bad=(1,1,1,3); A,d,_,_,D=coefficients(bad)
    need(F(D,A-d)==F(32,17)>1,'expansion alone is not enough')
    print(json.dumps({'status':'FINITE-EXACT controls of the proved conditional cone',
        'word_controls':checks,'exact_debt_denominator_controls':checks,
        'changing_generator_pairs':len(bank)**2,
        'base_kappa':'64/115','payment_margin_lower_bound':'51/115',
        'positive_cylinder_families':repetitions,
        'source_control':dict(stats),'unresolved_cone_hostiles':hostile,
        'expanding_block_hostile':{'word':bad,'kappa':'32/17'}},indent=2))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact finite controls for RESET_DEBT_COMPOSITION.md; no external data."""
from fractions import Fraction as F
from itertools import product
import json

checks = 0

def need(test, label):
    global checks
    checks += 1
    if not test:
        raise ValueError(label)

def valuation(n):
    if n <= 0:
        raise ValueError('valuation requires positive input')
    return (n & -n).bit_length() - 1

def coeff(word):
    A = 0
    C = 0
    S = 1
    for a in word:
        C = 3*C + 2**A
        A += a
        S *= 3
    M = 2**A
    r = ((M-C)*pow(S, -1, 2*M)) % (2*M)
    return S, M, C, r

def law(n, word):
    for a in word:
        z = 3*n+1
        if valuation(z) != a:
            return False
        n = z // 2**a
    return True

def replay(n, word, stopped=False):
    start = n
    literal = [n]
    odd = [n]
    for a in word:
        need(not stopped or n > 1, 'ROOT cannot be continued')
        need(n % 2 == 1, 'odd input')
        z = 3*n + 1
        need(valuation(z) == a, 'exact odd valuation')
        for _ in range(a):
            need(z % 2 == 0, 'literal halving exact')
            z //= 2
            literal.append(z)
        n = z
        odd.append(n)
    first = next((i for i, v in enumerate(literal) if v < start), None)
    return odd, literal, first

def compose(p, q):
    l1, d1 = p
    l2, d2 = q
    return l1*l2, d1+l1*d2

# General affine controls, including signed debts and multipliers below one.
parameters = [(l, d) for l in [F(1,2), F(8,9), F(1), F(4,3), F(2)]
              for d in [F(-1), F(0), F(1,3), F(5,9), F(1)]]
for p, q in product(parameters, repeat=2):
    l1,d1=p; l2,d2=q
    L,D=compose(p,q)
    e1=l1-1-d1; e2=l2-1-d2
    need(L-1-D == e1+l1*e2, 'weighted slack')
    for x in [F(0),F(1),F(3,2),F(11)]:
        need(((x+d1)/l1+d2)/l2 == (x+D)/L, 'chronological affine composition')
    need((1+D)/L - (1+d1)*(1+d2)/L == d2*e1/L,
         'defect is not multiplicative')
    if d1>=0 and d2>=0 and e1>=0 and e2>=0:
        need(D>=0 and D<=L-1, 'nonnegative cone closure')
    if l1>1 and l2>1:
        weighted=((l1-1)*d1/(l1-1)+l1*(l2-1)*d2/(l2-1))/(L-1)
        need(D/(L-1)==weighted, 'normalized debt barycentre')

# Independent exact-cylinder check: enumerate every positive odd residue in
# the full modulus, not only the residue claimed by the affine formula.
word_count=0
residues_checked=0
for k in range(1,5):
    for word in product(range(1,4),repeat=k):
        S,M,C,r=coeff(word)
        good=[]
        for n in range(1,2*M,2):
            residues_checked += 1
            if law(n,word):
                good.append(n)
        need(good==[r], 'full exact valuation cylinder')
        for t in range(3):
            n=r+2*M*t
            odd,literal,first=replay(n,word)
            need(odd[-1]*M==S*n+C, 'independent endpoint')
            need(len(literal)-1==sum(word), 'literal clock is valuation sum')
            z0=(S*r+C)//M
            need(odd[-1]==z0+2*S*t, 'inverse image progression')
            need(F(M*odd[-1]-C,S)==n, 'inverse only on true image')
        word_count += 1

# Infinite first-reset-bank formula, finite hostile controls through h=100.
first_reset_examples=[]
for h in range(101):
    a=2
    while 2**(h+a)<2*3**(h+1)-2**(h+1):
        a+=1
    word=(1,)*h+(a,)
    S,M,C,r=coeff(word)
    need((S,M,C)==(3**(h+1),2**(h+a),3**(h+1)-2**(h+1)),
         'first-reset explicit formula')
    need(F(S,M)<=F(3,4), 'uniform first-reset slope')
    L=F(M,S); D=F(C,S)
    need(D<=L-1, 'first-reset bank budget')
    n=r if r>1 else r+2*M
    odd,literal,first=replay(n,word,stopped=True)
    need(odd[-1]<n, 'first-reset descends')
    N=3**(h+1)*(n+1)-2**(h+1)
    j=next(j for j in range(1,a+1) if N<2**(h+j)*n)
    need(first==h+j, 'exact first-exit clock')
    if h<8:
        first_reset_examples.append({'h':h,'a':a,'residue':r,
                                     'modulus':2*M,'first_exit':first})

bank={'A':(1,2,3),'B':(1,3),'E':(1,1,3)}
expected={'A':(27,64,23,11,5),'B':(9,16,5,19,4),'E':(27,32,19,55,5)}
for name,w in bank.items():
    S,M,C,r=coeff(w)
    need((S,M,C,r)==expected[name][:4], 'concrete block data')
    for t in range(10):
        n=r+2*M*t
        odd,literal,first=replay(n,w,stopped=True)
        need(first==expected[name][4], 'uniform concrete first exit')
        need(odd[-1]<n, 'guarded concrete descent')
        if name=='A':
            need(odd[2]-n==2+16*t, 'first reset in A really grows')
    need(C<=4*(M-S), 'common barrier four')
    need(F(S,M)<=F(27,32), 'common slope bound')
    if name in ('A','B'):
        need(C<=M-S, 'unit bank')
    else:
        need(C>M-S, 'unit budget not necessary')

# Hostile reset: actual lawful orbit, not just an affine formula.
S,M,C,r=coeff((1,2))
need((S,M,C,r)==(9,8,5,11), 'growing reset formula')
need(replay(11,(1,2),True)[0]==[11,17,13], 'growing reset witness')
need(F(M,S)<1 and F(C,S)>0, 'positive debt obstructs unit cone')

# Every short word in three macros: use exact combined residue and replay all
# three lifts. Do not assume separate source cylinders may be chosen freely.
schedules_checked=0
max_clock=0
for k in range(1,7):
    for schedule in product(bank,repeat=k):
        w=sum((bank[s] for s in schedule),())
        S,M,C,r=coeff(w)
        L=F(1);D=F(0); E=F(0)
        for name in schedule:
            s,m,c,_=coeff(bank[name])
            l=F(m,s);d=F(c,4*s)
            E+=L*(l-1-d)
            L,D=compose((L,D),(l,d))
        need(L==F(M,S) and D==F(C,4*S), 'composition matches entire word')
        need(L-1-D==E and E>=0, 'all-length proof finite controls')
        for t in range(3):
            n=r+2*M*t
            odd,literal,first=replay(n,w,True)
            need(odd[-1]<n, 'whole schedule descends')
            need(F(odd[-1]-4)<=F(n-4)/L, 'barrier excess bound')
        schedules_checked += 1
        max_clock=max(max_clock,sum(w))

# Much longer words, including nonconstant order.
for schedule in [('A',)*200, ('B',)*200, ('E',)*200,
                 tuple('ABE'[(j*j+3*j+1)%3] for j in range(257))]:
    w=sum((bank[s] for s in schedule),())
    S,M,C,r=coeff(w)
    need(C<=4*(M-S), 'long schedule budget')
    odd,literal,first=replay(r,w,True)
    need(odd[-1]<r, 'long schedule physical descent')
    schedules_checked += 1
    max_clock=max(max_clock,sum(w))

# Finite-versus-infinite quantifier controls, including a rational 2-adic
# fixed point which is not a positive integer.
for k in range(1,101):
    n=1+2*4**k
    odd,_,_=replay(n,(2,)*k,True)
    for j in range(k+1):
        need(odd[j]==1+2*3**j*4**(k-j), 'all-two arbitrary finite depth')
    S,M,C,r=coeff(bank['A']*k)
    need(37*C==23*(64**k-27**k), 'repeated A intercept')
    need((37*r-23)%M==0, 'repeated A fixed 2-adic address')
    need((37*r-23)%(2*M)==0, 'full repeated A odd-terminal address')
need(F(27*23,37)+23==F(64*23,37), 'A rational fixed point')
need(23%37!=0, 'fixed point is not an ordinary integer')

print(json.dumps({'status':'PASS', 'scope':'Finite exact controls; all-length claims proved in companion note',
                  'checks':checks, 'independent_valuation_words':word_count,
                  'odd_residues_enumerated':residues_checked,
                  'macro_schedules_checked':schedules_checked,
                  'maximum_literal_clock':max_clock,
                  'first_reset_examples':first_reset_examples,
                  'hostiles':{'growing_first_reset':[11,17,13],
                              'unit_budget_not_necessary':[55,83,125,47],
                              'nested_A_limit':'23/37, not a positive integer'},
                  'global_collatz_conclusion':False},sort_keys=True,indent=2))

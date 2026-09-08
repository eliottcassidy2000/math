"""Independent exhaustive deletion and exact integral audit; no producer import."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
from fractions import Fraction as Q
import json, sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'05-knowledge/results' if HERE.name=='04-computation' else HERE
STEM=Path(__file__).stem
gates=0
def need(test,label):
    global gates
    gates+=1
    if not test: raise RuntimeError(label)
def board_universe(n):
    pairs=list(combinations(range(n),2))
    # Independent row-product universe; no recursive deficit pruning.
    for rows in product(pairs,repeat=n):
        degree=Counter(c for pair in rows for c in pair)
        if all(degree[c]==2 for c in range(n)): yield rows
def evaluate(rows):
    pts=[(r,c) for r,pair in enumerate(rows) for c in pair]
    m=len(pts); fam=[]
    for slope in (1,-1,2):
        lines={}
        for i,(r,c) in enumerate(pts):
            key=c-slope*r
            lines[key]=lines.get(key,0)|(1<<i)
        fam.append(list(lines.values()))
    active2=[mask for line in fam[:2] for mask in line if mask.bit_count()>2]
    active3=active2+[mask for mask in fam[2] if mask.bit_count()>2]
    def exact(active):
        for k in range(m+1):
            for indices in combinations(range(m),k):
                deletion=sum(1<<i for i in indices)
                if all((mask & ~deletion).bit_count()<=2 for mask in active): return k
        raise RuntimeError('empty survivor must be feasible')
    t2,t3=exact(active2),exact(active3)
    safe=(1<<m)-1
    for mask in active2: safe &= ~mask
    bonus=sum(max(0,(mask&safe).bit_count()-2) for mask in fam[2])
    singleton=[mask for line in fam[:2] for mask in line if mask.bit_count()==1]
    isolated=sum(mask.bit_count()==3 and all(sum(bool(one&(1<<i)) for one in singleton)==2
                 for i in range(m) if mask&(1<<i)) for mask in fam[2])
    def dual(lines):
        best=0
        for choice in range(1<<len(lines)):
            union=0
            for j,mask in enumerate(lines):
                if choice>>j&1: union |= mask
            best=max(best,union.bit_count()-2*choice.bit_count())
        return best
    b2,b3=dual(active2),dual(active3)
    need(t2==b2,'two-direction exact Boolean dual')
    need(t3>=b3>=t2+bonus>=t2+isolated,'independent third-direction chain')
    return (t2,t3,bonus,isolated,b3),active3,pts
census={}
for n in range(2,6):
    values=[evaluate(rows)[0] for rows in board_universe(n)]
    need(len(values)=={2:1,3:6,4:90,5:2040}[n],'complete independent board universe')
    need(all(v[3]==0 for v in values),'isolated event absent before six')
    census[n]={'boards':len(values),'strict_gain':sum(v[1]>v[0] for v in values)}
need(census[5]['strict_gain']==142,'exact strict-gain census')
named={
 'isolated':((0,1),(1,3),(0,5),(2,3),(2,4),(4,5)),
 'safe_two':((0,1),(2,3),(0,4),(1,2),(3,4)),
 'triangle':((0,3),(1,2),(0,4),(3,4),(1,2)),
 'naive':((0,1),(0,1),(2,3),(2,4),(3,4)),
 'dual_gain':((0,1),(0,3),(2,4),(3,4),(1,2))}
checks={key:evaluate(rows) for key,rows in named.items()}
need(checks['isolated'][0][:4]==(3,4,1,1),'minimal isolated witness')
need(checks['safe_two'][0][:3]==(0,2,2),'two disjoint extra deletions')
need(checks['dual_gain'][0][0:3]==(1,2,0) and checks['dual_gain'][0][4]==2,'Boolean dual improves safe subset')
values,lines,pts=checks['triangle']
need(values[1]==2 and values[4]==1 and len(lines)==3,'three-direction integral gap')
weights=[Q(1,2) if i in (0,5,6) else Q(0) for i in range(10)]
need(sum(weights)==Q(3,2) and all(sum(weights[i] for i in range(10) if mask>>i&1)>=1 for mask in lines),'fractional primal cost')
need(all(sum(Q(1,2) for mask in lines if mask>>i&1)<=1 for i in range(10)),'fractional dual same cost')
need(abs(s.Matrix([[int(mask>>i&1) for i in (0,5,6)] for mask in lines]).det())==2,'nonunimodular actual cell minor')
naive,_,pts=checks['naive']; occ=Counter(c-2*r for r,c in pts)
need(naive[1]==4<naive[0]+sum(max(0,v-2) for v in occ.values()),'overlapping deletion cost cannot simply be added')

# Independent geometry: integrate absolute values after all breakpoints are split.
q,r=s.symbols('q r',real=True)
limits=[(-2,-1,-q/2,1),(-1,s.Rational(-1,2),-q/2,(1-q)/2),
        (s.Rational(-1,2),0,-q/2,(1-q)/2),(0,1,0,(1-q)/2)]
Ls=[]; Js=[]
for qa,qb,lo,hi in limits:
    qm=s.Rational(qa+qb,2)
    breaks=[lo,hi]
    for point in [-q,(1-q)/3]:
        if bool((point-lo).subs(q,qm)>0) and bool((hi-point).subs(q,qm)>0): breaks.append(point)
    breaks=sorted(set(breaks),key=lambda z:float(s.sympify(z).subs(q,qm)))
    integral=0
    for a,b in zip(breaks,breaks[1:]):
        mid=s.sympify((a+b)/2)
        e1=r+q; e2=3*r+q-1
        sign1=s.sign(e1.subs(r,mid).subs(q,qm)); sign2=s.sign(e2.subs(r,mid).subs(q,qm))
        integral+=s.integrate(2-sign1*e1-sign2*e2,(r,a,b))
    Ls.append(s.expand(hi-lo));Js.append(s.factor(integral))
expected=[(q+2)*(q+5)/6,-(q-1)*(q+2)/3,-(q-1)*(q+2)/3,(q-4)*(q-1)/6]
for J,want in zip(Js,expected):need(s.simplify(J-want)==0,'exact seven-line geometry')
I3=sum(s.integrate(L**3,(q,qa,qb)) for L,(qa,qb,_,_) in zip(Ls,limits))
I4=sum(s.integrate(L**4,(q,qa,qb)) for L,(qa,qb,_,_) in zip(Ls,limits))
I2J=sum(s.integrate(L**2*J,(q,qa,qb)) for L,J,(qa,qb,_,_) in zip(Ls,Js,limits))
need((I3,I4,I2J)==(s.Rational(3,16),s.Rational(7,80),s.Rational(187,720)),'integral constants')
cost=s.simplify(2*(I4+3*I2J)/I3)
need(cost==s.Rational(416,45) and I3/6==s.Rational(1,32),'uniform Jensen mean cost and triple density')

cert={'status':'FINITE-EXACT independent controls; all-size proof in audit report',
 'universe':'all 2137 simple two-regular boards n=2..5 via full row product, plus five named controls',
 'census':census,'named_values':{k:list(v[0]) for k,v in checks.items()},
 'value_order':['tau2','tau3','J3','I3','beta3'],
 'integrals':list(map(str,(I3,I4,I2J,cost))),
 'bonus':'exp(-416/45)/4','always_active_gates':gates}
(OUT/(STEM+'_certificate.json')).write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
print('INDEPENDENT THREE-DIRECTION AUDIT: complete 2137 boards and five named controls')
print('CHAIN tau3 >= beta3 >= tau2+J3 >= tau2+I3; triangle 1 < 3/2 < 2')
print('GEOMETRY integral L^3=3/16, L^4=7/80, L^2J=187/720; cost=416/45')
print('UNIFORM BONUS exp(-416/45)/4; original-cell repair and column-transposition concentration')
print('Always-active exact gates:',gates)

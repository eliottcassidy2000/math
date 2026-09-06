#!/usr/bin/env python3
"""Exact norm-five continuum profiles and tiny inherited finite bases.

Uses SymPy for exact symbolic identities/inequalities, and a pinned literal-
audited H63 table as data. No repository mathematics is imported.
"""
import argparse
import csv
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from hashlib import sha256
from pathlib import Path
import sys
import sympy as sp

CHECKS=0
R=F(3,14)
t=sp.symbols('t',real=True)


def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(label)


def nonnegative(expr,lo,hi,label):
    bad=sp.reduce_inequalities([sp.Lt(sp.factor(expr),0),t>lo,t<hi],t)
    need(bad.as_set()==sp.EmptySet,label)


def sectors():
    return [
        ('I',sp.Rational(1,2)-t,(2,2,-1),sp.Rational(1,2),0,sp.Rational(1,4),
         [(2*t+7)/8,1-t/4,2*t*t-t+1]),
        ('II',1-t/2,(1,2,-2),(2+t)/4,0,sp.Rational(2,3),
         [(9*t+14)/16,1-t/16,(t*t-t+2)/2]),
        ('III',2-2*t,(2,1,-2),(2-t)/2,sp.Rational(1,2),sp.Rational(2,3),
         [(t+7)/8,(16-9*t)/8,2*t*t-3*t+2]),
        ('IV',t+sp.Rational(1,2),(2,-2,1),(1+t)/2,0,sp.Rational(1,2),
         [(16*t+7)/(8*(2*t+1)),sp.Integer(1),
          1+t-t*(4*t-1)**2/(4*(2*t+1))])]


def symbolic_profiles():
    signature=[]
    for name,b,v,K,lo,hi,J in sectors():
        w=(t,b,sp.Integer(1))
        need(sp.simplify(sum(v[i]*w[i] for i in range(3)))==0,('signed sector identity',name))
        for i in range(3):
            p=w[(i+1)%3]*w[(i+2)%3]
            s=w[(i+1)%3]+w[(i+2)%3]
            m=abs(v[i])
            nonnegative(s/m-K,lo,hi,('complete raw cutoff',name,i))
            g=sp.factor(m*K-s+2*p)
            formula=sp.factor(2*K-g*g/(2*m*p))
            need(sp.simplify(formula-J[i])==0,('exact cap-integral formula',name,i))
            if name=='IV' and i==2:
                nonnegative(-g,lo,sp.Rational(1,4),'IV third projection stays capped on initial range')
                nonnegative(g,sp.Rational(1,4),hi,'IV third cap activates only after quarter')
                need(sp.simplify(J[i].subs(t,sp.Rational(1,4))-(2*K).subs(t,sp.Rational(1,4)))==0,
                     'IV cap-switch values match')
            else:
                nonnegative(g,lo,hi,('cap-switch before full support endpoint',name,i))
        if name=='I':
            pieces=[(lo,sp.Rational(1,8),0),(sp.Rational(1,8),hi,2)]
            peak=sp.Rational(29,32)
            need(sp.factor(J[2]-J[0])==((2*t-1)*(8*t-1))/8,'I selector crossover factor')
        elif name=='II':
            pieces=[(lo,sp.Rational(1,8),0),(sp.Rational(1,8),hi,2)]
            peak=sp.Rational(121,128)
            need(sp.factor(J[2]-J[0])==((t-2)*(8*t-1))/16,'II selector crossover factor')
        elif name=='III':
            pieces=[(lo,sp.Rational(9,16),0),(sp.Rational(9,16),hi,2)]
            peak=sp.Rational(121,128)
            need(sp.factor(J[2]-J[0])==((t-1)*(16*t-9))/8,'III selector crossover factor')
        else:
            pieces=[(lo,sp.Rational(1,4),0),(sp.Rational(1,4),hi,0)]
            peak=sp.Rational(15,16)
        for left,right,pivot in pieces:
            actual=list(J)
            if name=='IV' and right<=sp.Rational(1,4):
                actual[2]=2*K
            for i in range(3):
                nonnegative(actual[i]-actual[pivot],left,right,('exact selector on piece',name,left,right,i))
            nonnegative(peak-actual[pivot],left,right,('exact continuum peak bound',name,left,right))
        displayed=list(map(str,J))
        if name=='IV':
            displayed[2]='1+t-t*max(4*t-1,0)^2/(4*(2*t+1))'
        signature.append((name,tuple(displayed),str(peak)))
    need(F(3,49)*F(121,128)==F(363,6272),'uniform selected continuum bound')
    return signature


def physical_integral_controls():
    # Direct all-roof integration over normalized z, unrelated to selecting J_i.
    examples=[(F(1,16),F(7,16),F(1)),(F(1,8),F(15,16),F(1)),
              (F(5,8),F(11,16),F(1)),(F(9,16),F(7,8),F(1)),
              (F(1,8),F(5,8),F(1)),(F(3,8),F(7,8),F(1))]
    for w in examples:
        cases=classify(w)
        need(bool(cases),('rational physical integral sector',w))
        for name,v in cases:
            cap=(F(2),F(0))
            lines=[cap]
            limits=[]
            for i in range(3):
                p=w[(i+1)%3]*w[(i+2)%3]
                s=w[(i+1)%3]+w[(i+2)%3]
                lines.append((s/p,-F(abs(v[i]),1)/p))
                limits.append(s/abs(v[i]))
            K=min(limits)
            cuts={F(0),K}
            for (a,b),(c,d) in combinations(lines,2):
                if b!=d:
                    x=(c-a)/(b-d)
                    if 0<x<K:
                        cuts.add(x)
            area=F(0)
            cs=sorted(cuts)
            for left,right in zip(cs,cs[1:]):
                middle=(left+right)/2
                A,B=min(lines,key=lambda L:L[0]+L[1]*middle)
                area+=A*(right-left)+B*(right*right-left*left)/2
            need(area==F(7,8),('independent full physical envelope integral',w,name))
    need(F(3,49)*F(7,8)==F(3,56),'physical continuum from saturated plane area')


def classify(w):
    a,b,c=w
    result=[]
    if c==2*(a+b):result.append(('I',(2,2,-1)))
    if 2*c==a+2*b:result.append(('II',(1,2,-2)))
    if 2*c==2*a+b:result.append(('III',(2,1,-2)))
    if 2*b==2*a+c:result.append(('IV',(2,-2,1)))
    return result


def exact_row(w,v):
    a,b,c=w
    q=F(3,7*c)
    K=min(F(3*(sum(w)-w[i]),14*abs(v[i])) for i in range(3))
    end=(K.numerator-1)//K.denominator
    E=[F(0),F(0),F(0)]
    mass=F(0)
    support=[]
    for k in range(1,end+1):
        if k%3==0:
            continue
        terms=[min(q,F(3*(sum(w)-w[i])-14*k*abs(v[i]),14*w[(i+1)%3]*w[(i+2)%3]))
               for i in range(3)]
        need(min(terms)>0,('full live multiplier',w,k))
        for i in range(3):E[i]+=2*terms[i]
        mass+=2*min(terms)
        support.extend((tuple(k*x for x in v),tuple(-k*x for x in v)))
    return tuple(E),mass,len(support)


def finite_bases(path):
    raw=Path(path).read_bytes()
    need(sha256(raw).hexdigest()=='c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa',
         'frozen independently literal-audited H63 data')
    with Path(path).open(newline='',encoding='utf-8') as f:
        table={tuple(int(row[k]) for k in ('a','b','c')):row for row in csv.DictReader(f,delimiter='\t')}
    need(len(table)==10074,'full table universe size')
    universe=[w for w in combinations([n for n in range(1,51) if n%3],3)
              if gcd(gcd(w[0],w[1]),w[2])==1 and classify(w)]
    need(len(universe)==174,'independently generated complete norm-five H50 universe')
    network_eq=[];physical_eq=[];head27=0;head45=0;sector_counts={}
    for w in universe:
        cases=classify(w)
        need(len(cases)==1,('unique sector under primitive ternary-unit filters',w))
        name,v=cases[0]
        sector_counts[name]=sector_counts.get(name,0)+1
        E,mass,N=exact_row(w,v)
        row=table[w];D=int(row['denominator'])
        need(E==tuple(F(int(row['E'+str(i)+'_numerator']),D) for i in range(3)),('every projection agrees with literal-audited table',w))
        need(mass==F(int(row['mass_numerator']),D) and N==int(row['raw_carriers']),('physical mass and complete multiplier count',w))
        need(min(E)<=F(46,665),('sharp selected-network finite head',w))
        need(mass<=F(51,770),('sharp physical finite head',w))
        if min(E)==F(46,665):network_eq.append(w)
        if mass==F(51,770):physical_eq.append(w)
        head27+=w[2]<=27
        head45+=w[2]<=45
    need(head27==48,'complete original H27 base')
    need(network_eq==[(2,19,20)] and physical_eq==[(1,11,20)],'sharp finite equality boundaries')
    E,mass,N=exact_row((10,11,16),(1,2,-2))
    need(E==(F(17,176),F(9,140),F(3,55)) and mass==F(331,6160), 'switched-roof hostile retained')
    need(min(E)-mass==F(1,1232),'no false fixed-projection identity')
    return len(universe),head27,head45,sector_counts,network_eq,physical_eq


def speed_vertices(v):
    vertices = set()
    for free in range(3):
        if v[free] == 0:
            continue
        other = [j for j in range(3) if j != free]
        for values in product((F(0),F(1)),repeat=2):
            w = [F(0)]*3
            for j,x in zip(other,values):
                w[j] = x
            w[free] = -sum(v[j]*w[j] for j in other)/v[free]
            if 0 <= w[free] <= 1:
                vertices.add(tuple(w))
    need(bool(vertices), ("nonempty normalized speed polytope",v))
    return sorted(vertices)


def fill_value(items, amount, reverse):
    total = F(0)
    for ratio, capacity in sorted(items, reverse=reverse):
        take = min(capacity,amount)
        total += take*ratio
        amount -= take
    need(amount == 0, "continuous resource is filled exactly")
    return total


def slice_width(v,w,delta):
    """LP width via fractional knapsack; no error polygon vertices or clipping."""
    pivot = next(j for j in range(3) if v[j])
    j,k = (pivot+1)%3,(pivot+2)%3
    coefficients = [F(0)]*3
    coefficients[j],coefficients[k] = -F(w[k],v[pivot]),F(w[j],v[pivot])
    items,free_width = [],F(0)
    for p,ell in zip(v,coefficients):
        if p:
            # Change e to sign(p)e, then shift [-R,R] to [0,2R].
            # The constant objective shift cancels between max and min.
            slope = ell*(1 if p>0 else -1)/abs(p)
            items.append((slope,2*R*abs(p)))
        else:
            free_width += 2*R*abs(ell)
    amount = F(delta)+R*sum(map(abs,v))
    need(0 <= amount <= 2*R*sum(map(abs,v)), "error slice resource range")
    return free_width+fill_value(items,amount,True)-fill_value(items,amount,False)


def pattern_slope(p):
    S = sum(p)
    defects = tuple(d for d in range(-3*S,3*S+1)
                    if 14*abs(d)<3*S and ((d%3 == 0) == all(x%3 for x in p)))
    unit = all(x%3 for x in p)
    rho = F(2 if unit else 1,3)
    intercept = F(4,3)*len(defects) if unit else F(len(defects))
    maximum = F(0)
    for negative in range(3):
        if p[negative] == 0:
            continue
        v = tuple(-x if j == negative else x for j,x in enumerate(p))
        for w in speed_vertices(v):
            value = rho*sum((slice_width(v,w,d) for d in defects),F(0))
            maximum = max(maximum,value)
    return defects,maximum,intercept



def cutoff_controls():
    A=F(363,6272)
    need(F(4,7)/(F(11,140)-A)==F(17920,649)<28,'requested strict norm-five 11/140 cutoff')
    need(F(4,7)/(F(46,665)-A)==F(340480,6731)<51,'sharp selected-network tail cutoff')
    need(F(4,7)/(F(51,770)-F(3,56))==F(1760,39)<46,'sharp physical tail cutoff')
    need(F(46,665)<F(11,140) and F(51,770)<F(11,140),'both sharp sector ceilings are strictly below nonadditive target')


def nonadditive_cutoff_certificate():
    # Reuses the already independently audited 747-pattern width algorithm,
    # but recomputes every record and checks its frozen complete digest.
    patterns=sorted({tuple(sorted(p)) for p in product(range(19),repeat=3)
                     if sum(x!=0 for x in p)>=2
                     and gcd(gcd(p[0],p[1]),p[2])==1
                     and sum(x%3==0 for x in p)<=1
                     and tuple(sorted(p)) not in ((0,1,1),(1,1,1),(1,1,2),(1,2,2))})
    need(len(patterns)==747,'complete inherited parity-free coefficient box')
    digest=sha256();thresholds=[]
    for p in patterns:
        ds,slope,B=pattern_slope(p)
        if not ds:
            need((p,slope,B)==((0,1,2),F(0),F(0)),'empty defect list discharged as N=0')
        need(slope<=F(15,98),('inherited coefficient slope',p))
        threshold=B/(F(11,60)-slope)
        need(threshold<=F(35280,199),('exact individual nonadditive count cutoff',p))
        thresholds.append((threshold,p,slope,B))
        digest.update(repr((p,ds,slope,B)).encode())
    need(digest.hexdigest()=='cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c',
         'every exact coefficient record matches the audited fourth-round box')
    top=max(thresholds)
    need(top==(F(35280,199),(7,13,18),F(17,147),F(12)),'complete worst individual box cutoff')
    need(top[0]<178,'entire small-coefficient box count-safe from178')
    alpha=F(6,49)+F(4,7*19)
    delta=F(11,60)-alpha
    A=F(2,7)/delta;B=F(4,3)/delta
    need((alpha,delta,A,B)==(F(142,931),F(1721,55860),F(15960,1721),F(74480,1721)),
         'large-coefficient count constants')
    need(53*A+B==F(920360,1721)<535,'S<=53 large-coefficient tail')
    g54=F(3*54*54,16)-54*A-B
    need(g54==F(18547,6884)>0 and F(3*54,8)>A,'S>=54 quadratic tail and monotonicity')
    return len(patterns),top,alpha,delta,A,B,g54,digest.hexdigest()


def old_ceiling_subset_controls(path):
    additive=[];norm4=[]
    with Path(path).open(newline='',encoding='utf-8') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            a,b,c=(int(row[k]) for k in ('a','b','c'))
            D=int(row['denominator'])
            mass=F(int(row['mass_numerator']),D)
            selected=min(F(int(row['E'+str(i)+'_numerator']),D) for i in range(3))
            data=((a,b,c),mass,selected)
            if c<=61 and a+b==c:additive.append(data)
            if c<=34 and (c==2*a+b or c==a+2*b or 2*b==a+c):norm4.append(data)
    need(len(additive)==146 and len(norm4)==88,'complete inherited old-ceiling subuniverses')
    need([row for row in additive if min(row[1:])<=F(6,77)]==[((1,4,5),F(1,28),F(1,28))],
         'only additive old-ceiling exception in complete H61')
    need([row for row in norm4 if max(row[1:])>F(6,77)]==[((2,11,20),F(11,140),F(11,140))],
         'only norm-four old-ceiling exceedance in complete H34')
    need([row for row in norm4 if F(6,77) in row[1:]]==[((1,5,11),F(6,77),F(6,77))],
         'only norm-four old-ceiling equality in complete H34')
    need(F(3,49)+F(4,7*35)==F(19,245)<F(6,77),'strict norm-four old-ceiling tail')
    need(F(9,98)-F(6,7*62)==F(237,3038)>F(6,77),'strict additive old-ceiling lower tail')
    return len(additive),len(norm4),F(19,245),F(237,3038)


def main():
    sys.stdout.reconfigure(newline='\n')
    parser=argparse.ArgumentParser()
    parser.add_argument('--head-table',required=True)
    args=parser.parse_args()
    profiles=symbolic_profiles()
    physical_integral_controls()
    head=finite_bases(args.head_table)
    cutoff_controls()
    nonadditive=nonadditive_cutoff_certificate()
    old_ceiling=old_ceiling_subset_controls(args.head_table)
    print('scope=primitive sorted distinct positive ternary-unit triples with a signed norm-five(1,2,2) relation')
    print('exact_continuum_profiles='+str(profiles))
    print('uniform_selected_continuum=363/6272; physical_continuum=3/56')
    print('complete_inherited_bases='+str(head))
    print('selected_network_canonical_reconstruction=46/665, equality iff(2,19,20); strict tail c>=51; matches incomingTHM4441')
    print('physical_canonical_reconstruction=51/770, equality iff(1,11,20); strict tail c>=46; matches incomingTHM4441')
    print('requested_11/140_bound=strict at all heights; continuum tail c>=28 plus48-row H27')
    print('fixed_roof_hostile=(10,11,16); minE-physical=1/1232')
    print('general_nonadditive_cutoff_certificate='+str(nonadditive))
    print('general_nonadditive_scope=analytic c>=535 reduction; fullH535 native result remains a separate premise')
    print('old_ceiling_combined_corollary_controls='+str(old_ceiling))
    print('optimization_live_gates='+str(CHECKS))


if __name__=='__main__':
    main()

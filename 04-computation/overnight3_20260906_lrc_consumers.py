"""Exact probes of the two inputs to the scale-three Haar consumer.

Independent literal phases; positive distinct integral speeds, ternary-unit
tails, no parity assumption. All-height proof is in the companion report.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, lcm, prod
import json
from pathlib import Path
import sys
import sympy as s

sys.stdout.reconfigure(newline="\n")
CHECKS=0


def need(ok, label):
    global CHECKS
    CHECKS+=1
    if not ok:
        raise RuntimeError(label)


def sheet(w, j, total):
    unit = total // (42 * w)
    residue = (-w * j) % 3
    intervals = []
    for k in range(w + 1):
        center = 14 * (residue + 3 * k) * unit
        left, right = max(0, center - 3 * unit), min(total, center + 3 * unit)
        if left < right:
            intervals.append((left, right))
    return intervals


def intersect(lists, closed=False):
    current = list(lists[0])
    for other in lists[1:]:
        answer = []
        i = j = 0
        while i < len(current) and j < len(other):
            left = max(current[i][0], other[j][0])
            right = min(current[i][1], other[j][1])
            if left < right or (closed and left == right):
                answer.append((left, right))
            if current[i][1] < other[j][1]:
                i += 1
            elif current[i][1] > other[j][1]:
                j += 1
            else:
                i += 1
                j += 1
        current = answer
    return current


def tail_mass(w):
    total = 42 * prod(w)
    deck = [[sheet(v, j, total) for j in range(3)] for v in w]
    length = sum(b-a for perm in permutations(range(3))
                 for a, b in intersect([deck[i][perm[i]] for i in range(3)]))
    return F(length, total)


def body_components(speeds):
    denominator = 14 * lcm(*speeds)
    safe = [(0, denominator)]
    for w in speeds:
        unit = denominator // (14*w)
        good = [((14*k+1)*unit, (14*k+13)*unit) for k in range(w)]
        safe = intersect([safe, good], closed=True)
    return denominator, safe


def body_mass(speeds):
    denominator, safe = body_components(speeds)
    return F(sum(b-a for a,b in safe), denominator)


def adjacent_formula(b):
    """Standalone physical integration for (1,b,b+1), b=1 mod3."""
    value = F(0)
    for k in range(1, (3*(b+1)-1)//14 + 1):
        if k % 3 == 0:
            continue
        t = F(k,3)
        if t <= F(b,14):
            value += (F(2*b+1,14)-t)/(b*(b+1))
        else:
            value += (F(b+1,14)-t)/b
    return 6*value


def additive_tent(a,b):
    c = a+b
    value = F(0)
    for k in range(1,(3*c-1)//14+1):
        if k%3 == 0:
            continue
        if 14*k <= 3*a:
            value += F(3,7*c)
        elif 14*k <= 3*b:
            value += F(3*(a+2*b)-14*k,14*b*c)
        else:
            value += F(3*c-14*k,14*a*b)
    return 2*value


def adjacent_polynomials():
    h=s.symbols('h',integer=True,nonnegative=True)
    result=[]
    for r in range(1,42,3):
        b=42*h+r
        K=9*h+(3*r)//14
        J=3*h+((3*r)//14)//3
        count=K-J
        sumk=(K*(K+1)-3*J*(J+1))/2
        numerator=s.expand((6*b+3)*count-14*sumk)
        for k0 in range((3*r)//14+1,(3*r+2)//14+1):
            if k0%3:
                numerator+=(b+1)*(3*b+3-14*(9*h+k0))
        numerator=s.expand(numerator)
        lower=s.expand(77*numerator-42*b*(b+1))
        upper=s.expand(42*b*(b+1)-55*numerator)
        # b>=7 requires h>=1 for residues1 and4; all others allow h>=0.
        shifted=s.Poly(lower.subs(h,h+1) if r in (1,4) else lower,h)
        need(all(v>0 for v in shifted.all_coeffs()),'strict all-height lower bound')
        up=s.Poly(upper,h)
        need(all(v>=0 for v in up.all_coeffs()),'all-height upper bound')
        need((up.eval(0)==0)==(r==10),'unique upper-bound equality')
        need(s.limit(numerator/(7*b*(b+1)),h,s.oo)==s.Rational(9,98),'uniform residue limit')
        for j in range(1 if r==1 else 0,6):
            need(F(numerator.subs(h,j))/(7*(42*j+r)*(42*j+r+1))==adjacent_formula(42*j+r),
                 'quasipolynomial matches physical finite sum')
        row={'residue':r,'numerator':str(numerator),'lower_difference':str(lower),'upper_difference':str(upper)}
        result.append(row)
        print('ADJACENT_RESIDUE',r,'NUMERATOR',numerator)
    return result


def main():
    controls={(1,5,7):F(8,245),(1,5,11):F(6,77),(2,5,7):F(22,245),
              (1,7,8):F(31,392),(1,10,11):F(6,55),(2,11,20):F(11,140),
              (1,2,4):F(0),(1,4,5):F(1,28)}
    for w,wanted in controls.items():
        need(tail_mass(w)==wanted,'literal physical control')
        print('TAIL_CONTROL',w,wanted)
    units = [n for n in range(1, 61) if n % 3]
    rows=[]
    bad = []
    for w in combinations(units, 3):
        if gcd(*w) != 1:
            continue
        mass = tail_mass(w)
        rows.append((mass,w))
        if mass > F(6,77):
            bad.append((w,mass))
    best=max(mass for mass,w in rows)
    ties=[w for mass,w in rows if mass==best]
    nonadditive=[(w,mass) for w,mass in bad if w[0]+w[1]!=w[2]]
    need(len(rows)==8664 and len(bad)==136,'complete primitive head')
    need(best==F(6,55) and ties==[(1,10,11)],'finite maximum')
    need(min((w[2],w) for w,mass in bad)==(7,(2,5,7)),'minimal parity hostile')
    need(nonadditive==[((2,11,20),F(11,140))],'nonadditive repair hostile')
    print('FINITE_HEAD',len(rows),'PRIMITIVE_TRIPLES_HEIGHT60 MAX',best,'TIES',ties,'VIOLATIONS',len(bad))
    print('NONADDITIVE_VIOLATIONS',nonadditive)
    count=0
    for a in range(1,30):
        for b in range(a+1,60):
            if a%3 and b%3 and (a+b)%3 and gcd(a,b)==1:
                need(tail_mass((a,b,a+b))==additive_tent(a,b),'independent additive-tent identity control')
                count+=1
    for b in range(4,1001,3):
        need(adjacent_formula(b)==additive_tent(1,b),'two analytic formulas agree')
    print('ADDITIVE_TENT_LITERAL_CONTROLS',count,'ADJACENT_CONTROLS',333)
    polynomials=adjacent_polynomials()
    C=(1,2,3,5,7,8,9,11,12,13)
    wanted=(F(1217,8820),F(17231,194040),F(16277,194040),F(14249,252252))
    for height in range(10,14):
        body_rows=[(body_mass(c),c) for c in combinations(range(1,height+1),10)]
        best_body=min(body_rows)
        fails=sum(mass<F(6,77) for mass,c in body_rows)
        need(best_body[0]==wanted[height-10],'finite body minimum')
        need(fails== (12 if height==13 else 0),'first failing body height')
        print('TEN_BODY_HEIGHT',height,'COUNT',len(body_rows),'MINIMUM',best_body,'BELOW_GATE',fails)
    denominator,components=body_components(C)
    isolated=[F(a,denominator) for a,b in components if a==b]
    need(isolated==[F(i,14) for i in (1,3,5,9,11,13)],'weak endpoint sidecar')
    need(body_mass(C)==F(14249,252252)<F(6,77),'recovered floor obstruction')
    full=tuple(3*c for c in C)+(1,5,11)
    need(body_mass(full)==F(25331,756756),'full-row exact safety measure')
    need(all(min((F(w,14)%1),1-(F(w,14)%1))>=F(1,14) for w in full),'known fixed-clock safety')
    print('RECOVERED_THM530_BODY',C,'MASS',body_mass(C),'ISOLATED',[str(x) for x in isolated])
    print('FULL_ROW_ODD_TAILS', (1,5,11),'MASS',body_mass(full),'KNOWN_SAFE_PHASE','1/14')
    certificate={'scope':'FINITE-EXACT tail head, recovered finite body extremizer, PROVED adjacent all-height formula',
                 'tail_height':60,'tail_universe':len(rows),'tail_violations':[[list(w),str(v)] for w,v in bad],
                 'adjacent_polynomials':polynomials,
                 'body_components':[[str(F(a,denominator)),str(F(b,denominator))] for a,b in components]}
    path=Path(__file__).resolve().parents[1]/'05-knowledge/results/overnight3_20260906_lrc_consumers_certificates.json'
    path.write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS exact parity and body-floor obstructions; sharp adjacent-family law')
    print('CHECKS',CHECKS)


if __name__ == '__main__':
    main()

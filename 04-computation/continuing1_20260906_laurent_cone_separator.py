"""Actual A2B3 obstruction to all positive Laurent-weighted coupled windows.

Standalone exact arithmetic, literal charge counts, original phase scaling.
No mathematical producer imports, floating roots, LP calls, or disabled gates.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def trim(a):
    a = list(map(F, a))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(*rows):
    return trim([sum((a[i] if i < len(a) else 0) for a in rows)
                 for i in range(max(map(len, rows)))])


def scale(a, c):
    return trim([c*v for v in a])


def cv(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def rem(a, p):
    a = trim(a)
    while len(a) >= len(p) and any(a):
        c, shift = a[-1]/p[-1], len(a)-len(p)
        for i, v in enumerate(p):
            a[i+shift] -= c*v
        a = trim(a)
    return a


def at(a, x):
    out = F(0)
    for c in reversed(a):
        out = out*x+c
    return out


def interval(a, I):
    lo = hi = F(0)
    for c in reversed(a):
        products = [v*x for v in (lo, hi) for x in I]
        lo, hi = min(products)+c, max(products)+c
    return lo, hi


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def ff(n, r):
    return factorial(n)//factorial(n-r) if 0 <= r <= n else 0


def coeff(a, i):
    return a[i] if 0 <= i < len(a) else F(0)


def divide_power(a, d):
    need(d >= 0 and all(c == 0 for c in a[:d]), 'exact positive phase-power removal')
    return trim(a[d:] or [0])


def carrier(kernel, g):
    # kernel[i] multiplies u^(2i) s^(degree(kernel)-i).
    q = len(kernel)-1
    return [trim([sum(choose(g, j-2*i)*kernel[i] if q-i == e else 0
                       for i in range(q+1)) for e in range(q+1)])
            for j in range(g+2*q+1)]


def u_product_coefficient(A, B, degree):
    return add(*[cv(A[j], B[degree-j]) for j in range(len(A))
                 if 0 <= degree-j < len(B)])


def actual(h):
    x, q, k, g = 1, h+1, 2*h+1, 3*h+2
    support = (-6*h-3, 1, 3*h+3)
    beta = [(-1)**(q-i)*choose(q+2*i, q-i) for i in range(q+1)]
    C = [(-1)**(q-1-i)*choose(q+1+2*i, q-1-i) for i in range(q)]
    D = [(-1)**(q-1-i)*choose(q+2*i, q-1-i) for i in range(q)]
    H, HC, HD = [carrier(a, g) for a in (beta, C, D)]
    n, power = g+2*q, 2*x-1
    first = [F((-1)**j*factorial(g), factorial(x+j)*factorial(3*h-3*j)*factorial(1+2*j))
             for j in range(h+1)]
    P = scale(first, 1/first[-1])
    target = [F((-1)**(j % 2)*factorial(2*g), factorial(2*x+j)*factorial(6*h-3*j)*factorial(2+2*j))
              for j in range(-1, 2*h+1)]
    need(H[k] == scale([0]*x+first, (-1)**x), 'original first-row zero with its exact phase power')
    need(n-k == 3*q and k < n-k, 'actual paired-grid size k')
    need(gcd(g, -support[0]) == 1, 'genuine first possible support return mass')

    # Literal charge enumeration, independent of the binomial-path compiler.
    for mass in range(1, 2*g+1):
        rows = {}
        for a in range(mass+1):
            for b in range(mass-a+1):
                c = mass-a-b
                if sum(i*v for i,v in zip((a,b,c),support)):
                    continue
                need(mass in (g,2*g), 'literal return-mass firewall')
                z = mass//g
                need((c-z) % 2 == 0, 'literal channel congruence')
                j = (c-z)//2
                rows[j] = factorial(mass)//(factorial(a)*factorial(b)*factorial(c))
        if mass == g:
            need(rows == {j:int((-1)**j*first[j]) for j in range(h+1)}, 'literal first multinomial row')
        elif mass == 2*g:
            need(rows == {j:int((-1)**(j % 2)*target[j+1]) for j in range(-1,2*h+1)}, 'literal doubled multinomial row including lower carry')
        else:
            need(not rows, 'no intervening support returns')

    windows = []
    pairs = [cv(H[k-d],H[k+d]) for d in range(1,k+1)]
    central = cv(H[k],H[k])
    for r in range(k):
        raw = add(*[scale(cv(H[j],H[2*k-j]),ff(j,r)*ff(2*k-j,r))
                    for j in range(r,2*k-r+1)])
        paired = add(scale(central,ff(k,r)**2),
                     *[scale(pairs[d-1],2*ff(k-d,r)*ff(k+d,r)) for d in range(1,k+1)])
        need(raw == paired, 'literal derivative square versus opposite pairs')
        windows.append(divide_power(raw,power))
    skip_product = divide_power(u_product_coefficient(HC,HD,2*k-2),power-1)
    need(target == add(windows[0],scale(skip_product,-2)), 'actual midpoint response with the indispensable factor s')
    need(rem(divide_power(central,power),P) == [0], 'central product vanishes in the true first-row quotient')

    # Independent triangular inversion of every allowed two-sided kernel.
    N = n-k
    for r in range(k):
        for z in range(N):
            values = [F(ff(k-d,r)*ff(k+d,r)*ff(N-d,z)*ff(N+d,z)) for d in range(1,k+1)]
            cone = []
            for t in range(k):
                d = k-t
                basis = lambda i: ff(k-d,i)*ff(k+d,i)
                cone.append((values[d-1]-sum(c*basis(i) for i,c in enumerate(cone)))/basis(t))
            need(all(c >= 0 for c in cone) and any(cone), 'complete two-sided cone reduction')
            for d in range(1,k+1):
                need(values[d-1] == sum(c*ff(k-d,t)*ff(k+d,t) for t,c in enumerate(cone)), 'every paired kernel reconstructed')
    return dict(h=h,x=x,q=q,k=k,g=g,n=n,support=support,P=P,first=first,
                H=H,HC=HC,HD=HD,beta=beta,C=C,D=D,windows=windows,target=target,
                residues=[rem(w,P) for w in windows],skip_product=skip_product)


def functional(a, y):
    return sum(c*coeff(a,i) for i,c in enumerate(y))


def isolate(p, seed, negative_polynomials):
    a,b = map(F,seed)
    need(at(p,a)*at(p,b)<0, 'initial disjoint positive root bracket')
    for step in range(160):
        if all(interval(row,(a,b))[1]<0 for row in negative_polynomials):
            break
        m=(a+b)/2
        value=at(p,m)
        need(value!=0, 'nonrational midpoint in retained isolation')
        if at(p,a)*value<0:
            b=m
        else:
            a=m
    else:
        raise ArithmeticError('root-sign isolation budget exhausted')
    need(at(p,a)*at(p,b)<0, 'retained exact root bracket')
    for row in negative_polynomials:
        need(interval(row,(a,b))[1]<0, 'uniform rational sign enclosure')
    return (a,b)


def shifted_remainder(a,p,d):
    base=[0,1] if d>=0 else scale(p[1:],-1/p[0])
    out=rem(a,p)
    for _ in range(abs(d)):
        out=rem(cv(out,base),p)
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate')
    args=parser.parse_args()
    cases=[actual(h) for h in (1,2,3)]
    certificate={}
    one,two,three=cases
    p=one['P']; ratio=F(6305,1714)
    need(rem(add(one['target'],scale(one['windows'][0],-ratio)),p)==[0], 'degree-one first-row positive constant control')
    print('POSITIVE h=1 x=1 support=(-9,1,6) T=(6305/1714)J0 modulo P')

    p=two['P'];y=(-1376,-139)
    need(p == [1,-10,1], 'actual quadratic first-row normalization')
    values=[functional(row,y) for row in two['residues']]
    target_value=functional(rem(two['target'],p),y)
    need(values == [872320,19900736,277245248,1956810240,4509872640], 'quadratic separator positive on every generator')
    need(target_value == -126496, 'quadratic actual response outside the complete constant cone')
    multiplier=[F(277182676568,191390856673),F(75869733936,191390856673)]
    need(all(c>0 for c in multiplier), 'positive linear phase repair')
    need(rem(add(two['target'],scale(cv(multiplier,two['windows'][0]),-1)),p)==[0], 'exact repaired actual identity')
    print('QUADRATIC support=(-15,1,9) P=s^2-10s+1 dual=',y,'window_values=',values,'target=',target_value)
    print('REPAIR multiplier=',','.join(map(str,multiplier)))
    certificate['quadratic']={'P':p,'residues':two['residues'],'T':rem(two['target'],p),'dual':y,'values':values,'target':target_value,'repair':multiplier}

    p=three['P'];y=(-32559418467680575845,-813366686033280291,-19065454274144095)
    need(p == [F(-1,3),14,-28,1], 'actual cubic first-row normalization')
    target=rem(three['target'],p)
    target_value=functional(target,y)
    need(target_value == -8569669792006914054574, 'cubic actual target strictly separated')
    moments=[]
    for row in three['residues']:
        e={d:functional(shifted_remainder(row,p,d),y) for d in (-1,0,1,2)}
        need(e[0]>=0 and e[1]>=0, 'two central nonnegative dual moments')
        need(e[-1]>3*e[0], 'strict backward convexity gate')
        need(e[2]>e[1], 'strict forward convexity gate')
        moments.append(e)
        # A separate quotient recurrence checks both directions at integer d.
        sequence={0:e[0],1:e[1],2:e[2]}
        for d in range(3,17):
            sequence[d]=28*sequence[d-1]-14*sequence[d-2]+sequence[d-3]/3
        for d in range(-1,-17,-1):
            sequence[d]=3*(sequence[d+3]-28*sequence[d+2]+14*sequence[d+1])
        for d in range(-16,17):
            direct=functional(shifted_remainder(row,p,d),y)
            need(direct == sequence[d], 'independent forward/backward cubic recurrence')
            need(direct>=0, 'signed Laurent-power hostile controls')
    need(moments[0][0]==0 and moments[6][1]==0, 'two defining dual boundary rays')
    N=[y[2]-28*y[1]+14*y[0],y[1]-28*y[0],y[0]]
    seeds=[(F(1,50),F(1,30)),(F(1,3),F(1)),(F(27),F(28))]
    signed_multiplier=[F(142115562175391338833022911,115962939903341750549938130),
                       F(137850584919079100136401223,115962939903341750549938130),
                       F(-825111792094668079242879,23192587980668350109987626)]
    need(signed_multiplier[0]>0 and signed_multiplier[1]>0 and signed_multiplier[2]<0,
         'rootwise-positive repair requires a negative coefficient in this representative')
    need(rem(add(three['target'],scale(cv(signed_multiplier,three['windows'][0]),-1)),p)==[0],
         'actual signed-coefficient multiplier identity')
    intervals=[isolate(p,I,[N,target,*three['residues'],scale(signed_multiplier,-1)]) for I in seeds]
    need(intervals[0][1]<intervals[1][0]<intervals[1][1]<intervals[2][0], 'three disjoint roots exhaust degree three')
    need(F(1,3)<intervals[1][0]<intervals[1][1]<1, 'middle-root bounds used in all-integer proof')
    # The actual doubled row is negative here; cone nonmembership is not noncancellation failure.
    need(functional(rem(add(three['target'],scale(three['windows'][0],-1)),p),y)==target_value,
         'same separator for the actual beta-skip payment')
    print('CUBIC support=(-21,1,12) first_mass=11 second_mass=22 P=s^3-28s^2+14s-1/3')
    print('DUAL',','.join(map(str,y)),'TARGET',target_value)
    print('LAGRANGE_NUMERATOR',','.join(map(str,N)))
    for r,e in enumerate(moments):
        print('FOUR_LEVELS r=',r,' values=',','.join(str(e[d]) for d in (-1,0,1,2)),
              'backward_slack=',e[-1]-3*e[0],'forward_slack=',e[2]-e[1])
    for i,I in enumerate(intervals):
        print('ROOT',i+1,'interval=',','.join(map(str,I)),'all_W_and_T_negative=True N_negative=True')
    certificate['cubic']={'P':p,'residues':three['residues'],'T':target,'dual':y,'target':target_value,
                           'moments':moments,'N':N,'intervals':intervals,'signed_multiplier':signed_multiplier,
                           'full_target':three['target'],'full_first':three['first'],
                           'beta':three['beta'],'C':three['C'],'D':three['D']}
    print('ROOTWISE_POSITIVE_REPAIR coefficient_signs=+,+,- multiplier=',','.join(map(str,signed_multiplier)))

    # Natural mixed/Euler terms: retain their type, not just their dual sign.
    q,g,k=three['q'],three['g'],three['k']
    A=carrier(three['beta'],g-1)
    C=carrier(three['C'],g-1)
    B=[[F(0)]]+[scale(row,F(2,3)) for row in C]
    length=max(len(A),len(B))
    lower=[add(A[i] if i<len(A) else [0],B[i] if i<len(B) else [0]) for i in range(length)]
    need(scale(three['H'][k],k)==scale(lower[k-1],g), 'inherited same-zero lowering with literal original coefficient')
    lower_parts=[u_product_coefficient(A,A,2*k-2),scale(u_product_coefficient(A,B,2*k-2),2),u_product_coefficient(B,B,2*k-2)]
    lower_parts=[divide_power(row,1) for row in lower_parts]
    lower_square=divide_power(u_product_coefficient(lower,lower,2*k-2),1)
    need(add(*lower_parts)==lower_square, 'complete lowered square retains its mixed term')
    lower_values=[functional(rem(row,p),y) for row in [*lower_parts,lower_square]]
    need(lower_values[0]<0 and lower_values[1]<0 and lower_values[2]>0 and lower_values[3]>0,
         'natural cross-wall terms lose the common-zero square predicate')
    need(rem(A[k-1],p)!=[0], 'the isolated binomial-divided carrier does not preserve the selected zero')
    print('LOWERED_PARTS dual_A2,dual_2AB,dual_B2,dual_L2=',','.join(map(str,lower_values)))
    certificate['lowered_parts']={'residues':[rem(row,p) for row in lower_parts], 'values':lower_values,
                                  'selected_A_remainder':rem(A[k-1],p)}
    # Cheap adversarial sign control for the invalid selected-zero deletion.
    # A=(1-2u)^2, H=(1+u)A has H_2=0 but [u^2]A^2=24>0.
    bad=[1,-4,4];badH=cv([1,1],bad)
    need(badH[2]==0 and coeff(cv(bad,bad),2)==24, 'division-by-binomial square-sign hostile')
    print('HOSTILE H=(1+u)(1-2u)^2 H2=0 but divided_carrier_square_coefficient=24')
    print('PASS',GATES,'always-active gates; all-integer Laurent separation is analytic, not inferred from the finite recurrence bank')
    if args.certificate:
        def convert(value):
            if isinstance(value,F):return str(value)
            if isinstance(value,dict):return {str(k):convert(v) for k,v in value.items()}
            if isinstance(value,(list,tuple)):return [convert(v) for v in value]
            return value
        Path(args.certificate).write_bytes((json.dumps(convert(certificate),sort_keys=True,indent=2)+'\n').encode('utf-8'))


if __name__=='__main__':
    main()

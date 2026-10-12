"""Exact controls for fixed-source charge, rhombus orders, and scale bills.

Run from the repo root. Python, SymPy, PySAT, and PARI/GP are required.
No assertion disappears under -O; no timeout is classified as a theorem.
The companion notes prove the infinite statements and separate PARI's
unconditional descent bounds from the finite enumerations.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, gcd, isqrt
import json
import hashlib
import subprocess
import sympy as s
from pysat.solvers import Solver
import source_context_localization_20261011 as field

checks = 0


def check(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)


def v2(n):
    return (n & -n).bit_length()-1


def orbit(n, cap=1000):
    states, word = [n], []
    while n != 1:
        if len(word) == cap:
            raise RuntimeError('orbit cap is pending, not a zero score')
        a = v2(3*n+1)
        n = (3*n+1) >> a
        states.append(n)
        word.append(a)
    return states, word


def charge(word):
    return sum(a % 2 for a in word)


def normalized(states, word):
    if not word:
        return Q(1)
    K = sum((a-1)//2 for a in word)
    return Q(3*states[0]+1, 4)*Q(1, 4**K)*Q(3, 4)**(len(word)-1)


def prefix_at_precision(n, d):
    """The U(1)=1 convention is for shadows, never a strict ROOT receipt."""
    word, A, P, B, x = [], 0, 1, 0, n
    while True:
        a = v2(3*x+1)
        if A+a >= d:
            break
        B, P = 3*B+2**A, 3*P
        A += a
        x = (3*x+1) >> a
        word.append(a)
    return word, A, P, B, x


def terminal_phase(y, t):
    """Unique k mod 3^t with 4^k = 3y+1 mod 3^(t+1)."""
    k = 0
    for j in range(1, t+1):
        modulus = 3**(j+1)
        candidates = [k+c*3**(j-1) for c in range(3)
                      if pow(4, k+c*3**(j-1), modulus) == (3*y+1) % modulus]
        check(len(candidates) == 1, 'unique ternary terminal lift')
        k = candidates[0]
    return k


def shadow(n, d):
    word, A, P, B, y = prefix_at_precision(n, d)
    phase = terminal_phase(y, len(word))
    k = phase
    threshold = max(y.bit_length()+2, (d-A+1)//2)
    if k < threshold:
        k += ((threshold-k+P-1)//P)*P
    h = (4**k-1)//3
    numerator = 2**A*h-B
    check(numerator % P == 0, 'terminal source is integral')
    m = numerator//P
    check(m > n and (m-n) % (2**d) == 0, 'marked source precision retained')
    states, actual = orbit(m, len(word)+2)
    check(actual == word+[2*k], 'strict first-hit shadow controller')
    check(charge(actual) == charge(word), 'terminal change adds no odd charge')
    z = normalized(states, actual)
    check(Q(3,4)*Q(3,8)**charge(word) < z <= Q(1,2**charge(word)),
          'shadow normalized score interval')
    return {'n': n, 'd': d, 'odd_charge': charge(word), 'phase': phase,
            'period': P, 'k': k, 'source_bits': m.bit_length()}


def odd_budget_deadline(n, budget):
    # t even letters give (n-1)*2^budget*(3/4)^t < 2.
    t = 0
    while (n-1)*2**budget*3**t >= 2*4**t:
        t += 1
    return budget+t


def collatz_controls():
    rows = []
    for n in range(1, 10002, 2):
        states, word = orbit(n)
        o = charge(word)
        z = normalized(states, word)
        check(Q(3,4)*Q(3,8)**o < z <= Q(1,2**o), 'all declared ROOT scores')
        check(len(word) <= odd_budget_deadline(n,o), 'charge-only deadline')
        # Independently check the multiplicative carry identity.
        factor = Q(1,2**o)
        for y in states[1:-1]:
            check(y >= 5 and y % 3 != 0, 'interior states avoid 3')
            factor *= Q(3*y,3*y+1)
        check(factor == z, 'independent normalized product')
        if n in (1,3,7,27,223,4591):
            rows.append({'n':n, 'odd_charge':o, 'odd_steps':len(word),
                         'deadline':odd_budget_deadline(n,o), 'Z':str(z)})
    shadows = [shadow(n,d) for n in range(1,64,2) for d in range(1,8)]
    for n in (1,3,7,27,223):
        ss, w = orbit(n)
        terminal_precision = sum(w)+1
        wd, _, _, _, _ = prefix_at_precision(n,terminal_precision)
        check(charge(wd) == charge(w), 'staircase stabilizes on proved source')
        if n>1:
            for budget in range(7):
                T=odd_budget_deadline(n,budget)
                H=1+2**budget*(n-1)
                D=1+(3*H+1).bit_length()*T
                used=0; A=0; stop=None
                for j,a in enumerate(w):
                    used+=a%2; A+=a
                    if used>budget or ss[j+1]==1:
                        stop=j+1; break
                check(stop is not None and stop<=T and A<D,
                      'finite source-specific event and precision invoice')
            budget=charge(w)
            lower=Q(3,3*n+1)*Q(3,8)**budget
            m=(n-1)//2
            nu=lambda a: Q(1,(a+1)**6)-Q(1,(a+2)**6)
            depth=1
            while True:
                q=(3*n+1)//2*(4**depth-1)//3
                error=Q(9,448*q**8)
                if error<=nu(m)*lower/2: break
                depth*=2
            for dd in (depth,2*depth,4*depth):
                q=(3*n+1)//2*(4**dd-1)//3
                coeff=nu(m)+Q(1,4**dd)*nu(m+q)
                check(lower-Q(9,448*q**8)/coeff >= lower/2,
                      'paid original-source floor survives nested refinements')
    # H(z)=3z/(3z+4) supplies an exact telescoping comparison.
    z = s.symbols('z')
    H = lambda t: 3*t/(3*t+4)
    check(s.cancel(H(z)/H(4*z/3)-(3*z+3)/(3*z+4)) == 0,
          'even-block telescoping identity')
    return {'universe':'odd n <= 10001; shadows odd n <= 63, precisions 1..7',
            'examples':rows, 'shadow_count':len(shadows),
            'shadow_examples':[r for r in shadows if r['n']==27],
            'universal_bounded_charge':'OPEN'}


def mod5(q):
    q=Q(q)
    check(q.denominator % 5 != 0, 'order reduction has no denominator 5')
    return q.numerator*pow(q.denominator,-1,5) % 5


def ring_residue(point):
    # alpha -> 1, sqrt(-55) -> 0.
    return sum(mod5(x) for x in point[:4]) % 5


def geometry_controls():
    f=field
    eta=(Q(-3,8),Q(0),Q(0),Q(0),Q(1,8),Q(0),Q(0),Q(0))
    R=[f.ZERO,f.neg(eta),f.scale(f.ONE,Q(3,4)),f.neg(f.conj(eta))]
    P=[f.ZERO]
    for i in range(4):
        P.append(f.add(P[-1],f.AP[i]))
    colour=(0,1,0,1,2)
    pts=P
    rows=[]
    for h in (1,0,2):
        rotation=f.power(eta,h)
        pts=list(dict.fromkeys(f.add(a,f.mul(rotation,b)) for a,b in product(pts,R)))
        residues=[ring_residue(p) for p in pts]
        edges=[]
        adj=[set() for _ in pts]
        for i,j in combinations(range(len(pts)),2):
            if f.norm(f.add(pts[i],f.neg(pts[j]))) == f.ONE:
                edges.append((i,j)); adj[i].add(j); adj[j].add(i)
                check((residues[i]-residues[j])%5 in (1,4), 'every actual unit difference')
                check(colour[residues[i]] != colour[residues[j]], 'explicit order colouring')
        boundary={pts.index(p) for p in P}
        check(len(boundary)==5, 'original pentagon coordinates retained')
        for i in range(5):
            check(pts.index(P[(i+1)%5]) in adj[pts.index(P[i])], 'original odd cycle retained')
        remain=set(range(len(pts))); interior=remain-boundary
        while True:
            peel=[i for i in interior if len(adj[i]&remain) <= 5]
            if not peel: break
            interior.difference_update(peel); remain.difference_update(peel)
        row={'added_rotation_power':h,'vertices':len(pts),'edges':len(edges),
             'relative_6_core':len(interior),'chi':3}
        if h==1:
            row['points']=[[str(x) for x in p] for p in pts]
            row['all_unit_edges']=edges
            row['colouring']=[colour[r] for r in residues]
        rows.append(row)
    check([(r['vertices'],r['edges'],r['relative_6_core']) for r in rows] ==
          [(18,35,0),(72,212,0),(216,780,211)], 'frozen overlap incidence census')
    # beta=sqrt(-11)=sqrt(-55)/sqrt(5), sqrt(5)=1+2(alpha+alpha^-1).
    sqrt5=f.add(f.ONE,f.scale(f.add(f.ALPHA,f.conj(f.ALPHA)),2))
    w=(Q(0),)*4+(Q(1),)+(Q(0),)*3
    beta=f.scale(f.mul(w,sqrt5),Q(1,5))
    check(f.mul(beta,beta)==f.scale(f.ONE,-11),'quadratic subfield coordinate')
    rho=f.scale(f.add(sqrt5,beta),Q(1,4))
    check(f.mul(rho,rho)==eta and f.norm(rho)==f.ONE,'half-angle saturates the order')
    escape=f.scale(f.add(f.scale(f.ONE,-7),f.scale(beta,4)),Q(1,15))
    check(f.norm(escape)==f.ONE,'five-pole escape direction')
    walk=[f.ZERO,f.ONE,f.add(f.ONE,escape),f.scale(f.ONE,Q(1,15))]
    for i,j in combinations(range(4),2):
        check((f.norm(f.add(walk[i],f.neg(walk[j])))==f.ONE)==(j==i+1),
              'complete induced rational-fifteenth path')
    # Norm discriminants for rational golden-chord spacings.
    D=s.symbols('D'); A=4*D-1-D**2; B=2*D-1
    check(s.expand(A*A+A*B-B*B-(D**2-7*D+1)*(D**2-3*D+1))==0,
          'golden seam norm factorization')
    rational_hits=[]
    for q in range(1,121):
        for p in range(1,2*q):
            if gcd(p,q)!=1: continue
            d=Q(p,q); D=d*d; A=4*D-1-D*D; B=2*D-1; N=A*A+A*B-B*B
            for twist in (1,5):
                z=N/twist
                if z>=0 and isqrt(z.numerator)**2==z.numerator and isqrt(z.denominator)**2==z.denominator:
                    rational_hits.append([str(d),twist])
    check(rational_hits==[['1',5]],'declared rational spacing census')
    # Recover the stronger full-field screen routed by THM-4558. The small
    # colouring was found afresh; verifying it requires neither SAT nor an
    # imported data file. Here i^2=-1 in F_(19^2).
    C19=[4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0]
    directions=[(x,y) for x,y in product(range(19),repeat=2) if (x*x+y*y)%19==1]
    differences=sorted({(x+y)%19 for x,y in directions})
    check(differences==[1,2,4,6,9,10,13,15,17,18], 'complete projected circle')
    for a,(x,y) in product(range(19),directions):
        check(C19[a]!=C19[(a+x+y)%19], 'all projected norm-one edges coloured')
    check(s.legendre_symbol(-11,19)==-1 and 19%5==4, 'Frobenius is conjugation in degree-eight field')
    return {'overlap_iterations':rows,'half_angle':'(sqrt(5)+sqrt(-11))/4',
            'escape_unit':'(-7+4sqrt(-11))/15', 'rational_spacing_scan':rational_hits,
            'field19_projection_colour':C19,
            'six_colour_context':'OPEN in Q(zeta55); impossible in Q(zeta5,sqrt(-11))'}


def elliptic_controls():
    program='''for(i=1,4,d=[1,5,55,11][i];E=ellinit([0,-6*d,0,5*d^2,0]);print([d,ellrank(E),elltors(E)]));\n'''
    run=subprocess.run(['gp','-q','-f'],input=program,text=True,capture_output=True,check=True)
    lines=run.stdout.strip().splitlines()
    check(len(lines)==4 and not run.stderr.strip(),'PARI descent completed')
    expected=['[1, [0, 0, 0, []], [4, [2, 2], [[5, 0], [0, 0]]]]',
              '[5, [0, 0, 0, []], [4, [2, 2], [[25, 0], [0, 0]]]]',
              '[55, [0, 0, 2, []], [4, [2, 2], [[275, 0], [0, 0]]]]',
              '[11, [1, 1, 0, [[180, 1950]]], [4, [2, 2], [[55, 0], [0, 0]]]]']
    check(lines==expected,'unconditional rank bounds and torsion, no analytic rank')
    check(1950**2==180*(180-11)*(180-55),'positive rank-one point')
    check(Q(180,605)==Q(6,11)**2,'rank-one auxiliary coordinate')
    check(Q(6,11)**2+Q(4,55)==Q(224,605),'remaining rhombus lift obligation')
    check(isqrt(224)**2!=224 and isqrt(605)**2!=605,'quotient point fails rational-square lift')
    return {'curves':'E_d: Y^2=X(X-d)(X-5d)','descent_output':lines,
            'positive_point_lifts_to_rhombus':False}


def full_field_controls():
    """Genuine degree-40 assemblies, exact ring filters, verified colour lifts."""
    z=s.symbols('z'); poly=s.Poly(s.cyclotomic_poly(55,z),z)
    relation=[(i,int(poly.nth(i))) for i in range(40) if poly.nth(i)]
    def reduce(a):
        a=list(a)+[0]*max(0,40-len(a))
        for d in range(len(a)-1,39,-1):
            if a[d]:
                for i,v in relation: a[d-40+i]-=a[d]*v
        return tuple(a[:40])
    zero=(0,)*40; one=(1,)+(0,)*39
    powers=[reduce([0]*k+[1]) for k in range(55)]
    def add(a,b): return tuple(x+y for x,y in zip(a,b))
    def scale(a,n): return tuple(n*x for x in a)
    def mul(a,b):
        c=[0]*79
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):
                    if y: c[i+j]+=x*y
        return reduce(c)
    def conj(a):
        v=zero
        for i,c in enumerate(a):
            if c: v=add(v,scale(powers[-i%55],c))
        return v
    alpha=powers[11]
    sqrt5=add(one,scale(add(alpha,conj(alpha)),2))
    beta=zero
    for j in range(1,11):
        beta=add(beta,scale(powers[5*j],int(s.legendre_symbol(j,11))))
    w=mul(sqrt5,beta)
    check(mul(w,w)==scale(one,-55),'degree-40 Gauss square root')
    gamma=powers[5]
    quintic=add(add(scale(powers[25],2),mul(add(one,scale(beta,-1)),powers[20])),
                add(scale(powers[15],-2),scale(powers[10],2)))
    quintic=add(quintic,add(scale(mul(add(one,beta),gamma),-1),scale(one,-2)))
    check(quintic==zero,'explicit quintic terminal for genuine field escape')
    check(s.totient(55)==5*8,'odd-degree circle-intersection barrier')
    R=[zero,add(scale(one,3),scale(w,-1)),scale(one,6),add(scale(one,3),w)]
    pent=[zero]; eleven=[zero]
    for j in range(4): pent.append(add(pent[-1],powers[11*j]))
    for j in range(10): eleven.append(add(eleven[-1],powers[5*j]))
    base=list(dict.fromkeys(add(a,b) for a,b in product(pent,eleven)))
    plain=list(dict.fromkeys(add(scale(a,8),b) for a,b in product(base,R)))
    eta=add(scale(one,-3),w)  # numerator; denominator 8
    rotated0=[mul(eta,b) for b in R]
    rotated1=[mul(mul(eta,powers[5]),b) for b in R]
    seam0=list(dict.fromkeys(add(scale(a,64),b) for a,b in product(base,rotated0)))
    seam1=list(dict.fromkeys(add(scale(a,64),b) for a,b in product(base,rotated1)))
    two=list(dict.fromkeys(add(scale(a,8),b) for a,b in product(plain,rotated1)))
    unions=[[],[]]
    for r in (1,3,4,5,9):
        # This is an automorphism over K: alpha and the quadratic Gauss sum
        # beta are fixed, while gamma goes to gamma^r.
        e=next(e for e in range(1,55) if e%5==1 and e%11==r)
        image_beta=zero
        for j in range(1,11):
            image_beta=add(image_beta,scale(powers[5*j*e%55],int(s.legendre_symbol(j,11))))
        check(image_beta==beta,'Galois controller fixes original rhombus field')
        pg=[zero]
        for j in range(10): pg.append(add(pg[-1],powers[5*r*j%55]))
        bb=list(dict.fromkeys(add(a,b) for a,b in product(pent,pg)))
        rr=[mul(mul(eta,powers[5*r%55]),b) for b in R]
        unions[0].extend(add(scale(a,64),b) for a,b in product(bb,rr))
        pp=list(dict.fromkeys(add(scale(a,8),b) for a,b in product(bb,R)))
        unions[1].extend(add(scale(a,8),b) for a,b in product(pp,rr))
    unions=[list(dict.fromkeys(points)) for points in unions]
    prime=331; root=pow(int(s.primitive_root(prime)),6,prime)
    check(pow(root,55,prime)==1 and pow(root,11,prime)!=1 and pow(root,5,prime)!=1,
          'exact order 55 in filter field')
    check(int(poly.eval(root))%prime==0,'filter evaluation is a ring homomorphism')
    rp=[pow(root,k,prime) for k in range(40)]; rc=[pow(root,-k,prime) for k in range(40)]
    rows=[]
    for name,pts,den in [('P5+P11',base,1),('base+R',plain,8),
                         ('base+etaR',seam0,64),('base+eta*zeta11*R',seam1,64),
                         ('base+R+eta*zeta11*R',two,64),
                         ('Galois union of single-rhombus seams',unions[0],64),
                         ('Galois union of double-rhombus seams',unions[1],64)]:
        vals=[(sum(x*y for x,y in zip(a,rp))%prime,
               sum(x*y for x,y in zip(a,rc))%prime) for a in pts]
        edges=[]; survivors=0
        for i,j in combinations(range(len(pts)),2):
            possible=(vals[i][0]-vals[j][0])*(vals[i][1]-vals[j][1])%prime == den*den%prime
            if possible or name=='P5+P11':
                diff=add(pts[i],scale(pts[j],-1))
                actual=mul(diff,conj(diff))==scale(one,den*den)
                if not possible: check(not actual,'independent complete small filter control')
                if actual: edges.append((i,j))
            survivors+=possible
        clauses=[[3*i+c+1 for c in range(3)] for i in range(len(pts))]
        clauses += [[-(3*i+c+1),-(3*j+c+1)] for i,j in edges for c in range(3)]
        with Solver(name='g3',bootstrap_with=clauses) as solver:
            check(solver.solve(),'finite full-field three-colouring exists')
            model=set(solver.get_model())
            colours=[next(c for c in range(3) if 3*i+c+1 in model) for i in range(len(pts))]
        for i,j in edges: check(colours[i]!=colours[j],'returned colour lift on every unit edge')
        ports=[pts.index(scale(v,den)) for v in pent]
        for i in range(5): check(tuple(sorted((ports[i],ports[(i+1)%5]))) in edges,'pinned C5 lower bound')
        if name in ('base+R','base+R+eta*zeta11*R','Galois union of double-rhombus seams'):
            check(all(scale(v,den//8) in pts for v in R),'original four rhombus ports retained')
        digest=hashlib.sha256(json.dumps([pts,edges],separators=(',',':')).encode()).hexdigest()
        rows.append({'name':name,'vertices':len(pts),'edges':len(edges),
                     'exact_norm_tests_after_filter':int(survivors),'chi':3,
                     'geometry_and_edges_sha256':digest,'colouring':colours})
    check([(r['vertices'],r['edges']) for r in rows]==[(54,108),(216,648),(188,550),(206,613),(824,3276),(994,3063),(3976,16228)],
          'frozen genuine full-field universe')
    return rows


def scales_controls():
    x,q=1,0
    pell=[]
    for k in range(1,13):
        x,q=3*x+8*q,x+3*q
        check(x*x-8*q*q==1,'inherited square-triangular Pell orbit')
        check(((x-1)//2)*((x+1)//2)==2*q*q,'square triangular identity')
        D=1+63*q*q; A=1-47*q*q; B=2*x*q
        check(A*A+55*B*B==D*D,'Pell-generated norm-one direction')
        pell.append({'k':k,'x':x,'q':q,'denominator_before_reduction':D})
    m,c,v=s.symbols('m c v')
    M=[s.integrate(s.integrate(v**6*(m+s.Rational(2,3)+v*(c-s.Rational(2,3)))**j,
                            (v,0,1)),(c,1,2)) for j in range(3)]
    coeff=[s.expand(28*(-1)**j*comb(7+j,j)*M[j]) for j in range(3)]
    check(coeff[0]==4 and coeff[1]==-32*m-s.Rational(134,3),'eighth/ninth-power coefficients')
    # Exact pointwise Taylor bounds imply the infinite alias bounds by summation.
    def kernel(x): return 2*(Q(1,(x+1)**6)-Q(1,(x+2)**6))/(3*x+2)
    for mm,qq,k in product((0,1,13,111),(1,2,5,19),(2,3,7,31)):
        actual=kernel(mm+k*qq)
        terms=[Q(coeff[j].subs(m,mm))/(k*qq)**(8+j) for j in range(3)]
        check(terms[0]+terms[1] <= actual <= terms[0]+terms[1]+terms[2],
              'alternating exact capacity bounds')
        check(256*kernel(mm+2*k*qq)-actual > 0,'positive two-scale capacity difference')
    # A genuine arbitrary atomic law defeats importing capacity monotonicity.
    # Its sole mass is at the odd alias m+q, absent from the 2q cylinder.
    check(256*0-1 < 0,'positive capacity is not positive actual-measure readout')
    return {'pell_first_rows':pell[:5], 'capacity_coefficients':[str(a) for a in coeff],
            'odd_zeta_scale_limit':'(16*m+67/3)*(zeta(9)-1)',
            'new_irrationality_or_ROOT_positivity':'NOT ESTABLISHED'}


if __name__=='__main__':
    out={'collatz':collatz_controls(),'geometry':geometry_controls(),'full_field':full_field_controls(),
         'elliptic':elliptic_controls(),'scales':scales_controls()}
    out['checks']=checks
    print(json.dumps(out,sort_keys=True,indent=2))

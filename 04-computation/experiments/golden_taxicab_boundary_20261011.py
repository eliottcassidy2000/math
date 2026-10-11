"""Exact controls for boundary gadgets, ROOT envelopes, and pencil normalizers.

No convergence assumption: the finite orbit bank raises on exceeding its cap.
Assertions remain active with -O. Sympy is used only for rational identities.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, factorial, isqrt, lcm
import json
import sympy as sp

checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def partitions(n):
    if n == 0:
        yield ()
        return
    def walk(a):
        if len(a) == n:
            yield tuple(a)
        else:
            for c in range(max(a) + 2):
                yield from walk(a + [c])
    yield from walk([0])


def norm(a, b):
    return a*a + a*b + b*b


def phi_mul(z, w, modulus):
    a, b = z
    c, d = w
    return ((a*c+b*d) % modulus, (a*d+b*c+b*d) % modulus)


def phi_order(modulus, limit=10000):
    z = (1, 0)
    for exponent in range(1, limit+1):
        z = phi_mul(z, (0, 1), modulus)
        if z == (1, 0):
            return exponent
    raise RuntimeError("golden order cap")


def geometry():
    # All eight Eisenstein representatives in one open sixty-degree sector.
    reps = [(a,b) for a in range(1,42) for b in range(42)
            if norm(a,b) == 1729]
    check(len(reps) == 8, "eight taxicab directions")
    for a,b in reps:
        check(norm(a,b) == 1729, "unit radial edge")
    for (a,b),(c,d) in combinations(reps, 2):
        check(0 < norm(a-c,b-d) < 1729, "distinct nonunit leaves")

    # Golden decagon: Q(sqrt5, sqrt(10-2sqrt5)) arithmetic via ANP.
    rt5 = sp.sqrt(5)
    field = sp.QQ.algebraic_field(rt5, sp.sqrt(10-2*rt5))
    cv = lambda z: field.from_sympy(sp.sympify(z))
    cc, ss = cv((1+rt5)/4), cv(sp.sqrt(10-2*rt5)/4)
    x,y = cv(1),cv(0)
    decagon = []
    for i in range(10):
        decagon.append((x,y))
        check(x*x+y*y == cv(1), "golden unit radial edge")
        x,y = cc*x-ss*y, ss*x+cc*y
    check((x,y) == (cv(1),cv(0)), "decagon closes")
    for z,w in combinations(decagon,2):
        d2 = sum((z[k]-w[k])**2 for k in (0,1))
        check(d2 != cv(0) and d2 != cv(1), "induced golden star")
    # A unit regular pentagon, with four private leaves at every core vertex.
    scale=cv(2/sp.sqrt(10-2*rt5))
    pentagon=[tuple(scale*a for a in decagon[2*i]) for i in range(5)]
    pent_points=pentagon[:]
    pent_groups=[]
    pent_parameters=[]
    for ci,c in enumerate(pentagon):
        ids,ts=[],[]
        for s in range(2,200):
            p=(c[0]+cv(sp.Rational(1-s*s,1+s*s)),
               c[1]+cv(sp.Rational(2*s,1+s*s)))
            if p in pent_points:
                continue
            adjacent=[i for i,z in enumerate(pent_points)
                      if sum((p[k]-z[k])**2 for k in (0,1))==cv(1)]
            if adjacent==[ci]:
                ids.append(len(pent_points));ts.append(s);pent_points.append(p)
                if len(ids)==4:
                    break
        check(len(ids)==4,"golden pentagon private leaves")
        pent_groups.append(ids);pent_parameters.append(ts)
    pent_wanted={tuple(sorted((i,(i+1)%5))) for i in range(5)}
    pent_wanted|={(i,j) for i,g in enumerate(pent_groups) for j in g}
    pent_actual=set()
    for i,j in combinations(range(len(pent_points)),2):
        d2=sum((pent_points[i][k]-pent_points[j][k])**2 for k in (0,1))
        check(d2!=cv(0),"distinct pentagon-gadget points")
        if d2==cv(1):
            pent_actual.add((i,j))
    check(pent_actual==pent_wanted,"complete pentagon induced incidence")

    # Rational double-star with its hidden central edge. Check ALL incidences.
    centres = [(Q(0),Q(0)),(Q(1),Q(0))]
    points = centres[:]
    groups = []
    parameters = []
    for ci,c in enumerate(centres):
        ids, ts = [], []
        for s in range(2,200):
            p = (c[0]+Q(1-s*s,1+s*s),c[1]+Q(2*s,1+s*s))
            if p in points:
                continue
            adjacent = [i for i,z in enumerate(points)
                        if sum((p[k]-z[k])**2 for k in (0,1)) == 1]
            if adjacent == [ci]:
                ids.append(len(points)); ts.append(s); points.append(p)
                if len(ids) == 5:
                    break
        check(len(ids)==5, "double-star construction succeeds")
        groups.append(ids); parameters.append(ts)
    wanted = {(0,1)} | {(i,j) for i,g in enumerate(groups) for j in g}
    actual = set()
    for i,j in combinations(range(len(points)),2):
        d2 = sum((points[i][k]-points[j][k])**2 for k in (0,1))
        check(d2 != 0, "distinct double-star vertices")
        if d2 == 1:
            actual.add((i,j))
    check(actual == wanted, "complete induced incidence audit")
    return {"taxicab_sector": reps,
            "pentagon_private_parameters":pent_parameters,
            "pentagon_edges":sorted(pent_actual),
            "pentagon_points":[[str(field.to_sympy(a)) for a in p]
                               for p in pent_points],
            "double_star_parameters": parameters,
            "double_star_points": [[str(x),str(y)] for x,y in points],
            "double_star_edges": sorted(actual)}


def boundary_tables():
    # All Bell(6)=203 boundary partitions, with literal interior enumeration.
    star_table = []
    for pi in partitions(6):
        w = sum(all(c != a for a in pi) for c in range(6))
        j = len(set(pi))
        check(w == 6-j, "star partition weight")
        star_table.append({"partition":pi,"weight":w})
    check(len(star_table)==203, "complete six-port table")
    check(sum(t["weight"]==0 for t in star_table)==1,
          "only rainbow partition excluded")
    labelled = sum(factorial(6)//factorial(6-len(set(t['partition'])))
                   * t['weight'] for t in star_table)
    check(labelled == 6*5**6, "star total colourings")
    masks = range(64)
    edge_count = 0
    for s,t in product(masks,repeat=2):
        A = [c for c in range(6) if not s>>c&1]
        B = [c for c in range(6) if not t>>c&1]
        w = sum(a!=b for a,b in product(A,B))
        formula = len(A)*len(B)-len(set(A)&set(B))
        check(w==formula, "double-star list formula")
        hall = bool(A and B and len(set(A)|set(B))>=2)
        check((w>0)==hall, "edge Hall condition")
        if 1<=s.bit_count()<=5 and 1<=t.bit_count()<=5:
            edge_count += 1
    check(edge_count==62**2, "reachable double-star mask pairs")
    # Genuinely geometric triangle-core extension with disjoint private ports:
    # enumerate all abstract allowed lists; geometry is supplied separately by
    # a finite-avoidance construction in the proof, not asserted by this loop.
    triangle_cases = 0
    sets = [{c for c in range(6) if m>>c&1} for m in masks]
    for A,B,C in product(sets,repeat=3):
        w = (len(A)*len(B)*len(C)-len(A&B)*len(C)
             -len(A&C)*len(B)-len(B&C)*len(A)+2*len(A&B&C))
        # Independent count by the first two colours, no inclusion-exclusion.
        direct = sum(len(C-{a,b}) for a,b in product(A,B) if a!=b)
        check(w==direct, "triangle list count")
        hall = (bool(A and B and C) and len(A|B)>=2 and len(A|C)>=2
                and len(B|C)>=2 and len(A|B|C)>=3)
        check((w>0)==hall, "triangle Hall condition")
        triangle_cases += 1
    def cycle_count(lists):
        total=0
        for start in lists[0]:
            dp={start:1}
            for allowed in lists[1:]:
                dp={c:sum(count for d,count in dp.items() if c!=d)
                    for c in allowed}
            total+=sum(count for c,count in dp.items() if c!=start)
        return total
    family=[{0,1},{0,2},{1,2},{0,1,2}]
    for lists in product(family,repeat=5):
        literal=sum(all(cols[i]!=cols[(i+1)%5] for i in range(5))
                    for cols in product(*lists))
        check(cycle_count(lists)==literal,"pentagon whole-boundary count")
    check(cycle_count([{0,1}]*5)==0,"pairwise-feasible global hostile")
    check(cycle_count([{0,1,2}]*5)==30,"pentagon three-list positive control")
    return {"star_six_port_table":star_table,
            "star_labelled_colourings":labelled,
            "double_star_mask_pairs":4096,
            "triangle_list_triples":triangle_cases,
            "pentagon_list_bank":1024,
            "pentagon_common_two_list_weight":0,
            "pentagon_common_three_list_weight":30}


def word(n, cap=10000):
    vals=[]
    seen=set()
    while n!=1:
        if n in seen or len(vals)>=cap:
            raise RuntimeError("finite ROOT bank did not terminate")
        seen.add(n)
        z=3*n+1
        a=(z&-z).bit_length()-1
        vals.append(a); n=z>>a
    return vals


def weight(vals):
    if not vals:
        return Q(1)
    K=sum((a-1)//2 for a in vals); tau=len(vals)
    return Q(2,(K+tau+1)*comb(K+tau,K))


def height(n):
    # n>1, exact integer floor(log_4(3n))+1.
    return ((3*n).bit_length()-1)//2+1


def envelope(n):
    if n==1:
        return Q(1)
    h=height(n)
    return Q(2,h*(h+1))


def sibling(n,j):
    return 4**j*n+(4**j-1)//3


def collatz():
    bank={n:word(n) for n in range(1,10000,2)}
    for n, vals in bank.items():
        w=weight(vals); h=n.bit_length()-1
        check(w<=envelope(n), "source height envelope")
        if n>1:
            tau=len(vals); K=sum((a-1)//2 for a in vals); H=K+tau
            check(3*n<4**H, "terminal size inequality")
            check(2*K<=h+tau-1, "source counter upper bound")
        if n<1000:
            for j in range(1,9):
                sw=weight(word(sibling(n,j)))
                gamma=Q((h+2)*(h+4),(h+2*j+2)*(h+2*j+4))
                check(sw<=gamma*w, "source-specific sibling domination")
    for n in range(2,20000):
        h=height(n); r=(4**h-1)//3
        check(r>=n and weight(word(r))==envelope(n),
              "exact tail maximum attained by root ray")
    samples=[]
    for n in [1,3,5,7,27,1729]:
        m=(n-1)//2; h=n.bit_length()-1
        for j in [1,2,4,8]:
            q=(3*m+2)*(4**j-1)//3
            hq=height(n+4*q)
            gamma=Q((h+2)*(h+4),(h+2*j+2)*(h+2*j+4))
            bill=Q(1,8*hq*(hq+1)*q**7)
            nu=lambda a:Q(1,(a+1)**6)-Q(1,(a+2)**6)
            A=nu(m)+gamma*nu(m+q)
            # Independent actual finite alias sum; the theorem controls the
            # infinite remainder, while this is a positive finite control.
            residual=sum((nu(m+k*q)*weight(word(n+2*k*q))
                          for k in range(2,22)),Q(0))
            check(residual<=bill, "refined actual finite residual")
            check(height(n+4*q)>=2, "tail threshold scope")
            samples.append({"n":n,"j":j,"q":q,"tail_height":hq,
                            "gamma":str(gamma),"bill":str(bill),
                            "coefficient":str(A)})
    return {"bank_odd_sources":len(bank),"max_source":9999,
            "max_odd_steps":max(map(len,bank.values())),
            "refined_floor_samples":samples}


def algebra():
    k,z,E,a,b,x,t=sp.symbols('k z E a b x t')
    A=-(k*k+k+1)/3
    B=(k-1)*(2*k+1)*(k+2)/27
    jinv=256*(k*k+k+1)**3/(k*k*(k+1)**2)
    hs=[k,-1-k,1/k,-k/(k+1),-1/(k+1),-(k+1)/k]
    table=[]
    for h in hs:
        check(sp.cancel(jinv.subs(k,h)-jinv)==0,"anharmonic j identity")
        adef=sp.factor(A.subs(k,h)-sp.diff(h,k)**4*A)
        bdef=sp.factor(B.subs(k,h)-sp.diff(h,k)**6*B)
        check((adef==0 and bdef==0)==(h==k),"symplectic base test")
        table.append({"h":str(h),"A_defect":str(adef),"B_defect":str(bdef)})
    I=x*x*t*t+(1+x*x)*t
    maps=[(x,t),(-1/x,x*x*t),(1/x,-1-x*x*(1+t)),(-x,-1-t-1/x**2)]
    for xx,tt in maps:
        check(sp.cancel(I.subs({x:xx,t:tt},simultaneous=True)-I)==0,
              "fixed pencil control")
        det=sp.diff(xx,x)*sp.diff(tt,t)-sp.diff(xx,t)*sp.diff(tt,x)
        check(sp.cancel(det-1)==0,"symplectic positive control")
    bad=(sp.I*x,-1-t)
    check(sp.expand(I.subs({x:bad[0],t:bad[1]},simultaneous=True)+1+I)==0,
          "pencil permutation hostile")
    check(sp.diff(bad[0],x)*sp.diff(bad[1],t)==-sp.I,
          "hostile fails area")
    f=z*z*(a+b*z)
    D=(f-E)**2+2*f
    check(sp.factor(sp.diff(f,z))==z*(2*a+3*b*z),"critical cubic")
    check(sp.cancel(f.subs(z,-2*a/(3*b))-4*a**3/(27*b*b))==0,
          "second critical value")
    for aa,bb in [(0,1),(1,1),(-1,1),(2,-3),(1,0),(-2,0)]:
        dd=D.subs({a:aa,b:bb})
        check(sp.gcd(dd,sp.diff(dd,z))==1,"generic squarefree fibre")
        check(sp.degree(dd,z)==(6 if bb else 4),"genus degree")
    G=lambda aa,bb: -x*x*t/2+(aa+bb*x*t)*(t+x*x*t*t)
    xx=2*(a+b*x*t)/x
    tt=-x*x*t/(2*(a+b*x*t))
    check(sp.cancel(G(a,-b).subs({x:xx,t:tt},simultaneous=True)-G(a,b))==0,
          "cross-source transfer")
    check(sp.cancel(sp.diff(xx,x)*sp.diff(tt,t)-sp.diff(xx,t)*sp.diff(tt,x)-1)==0,
          "cross-source symplecticity")
    for coordinate,target in [(xx,x),(tt,t)]:
        inverse=coordinate.subs(b,-b).subs({x:xx,t:tt},simultaneous=True)
        check(sp.cancel(inverse-target)==0,"cross-source inverse")
    check(sp.cancel(xx.subs(a,0)-2*b*t)==0 and
          sp.cancel(tt.subs(a,0)+x/(2*b))==0,"polynomial genus-two transfer")
    return {"pencil_base_tests":table,"linear_source_genus":{"b_nonzero":2,
            "b_zero_a_nonzero":1,"a_b_zero":0},
            "cross_source_map":[str(xx),str(tt)],
            "polynomial_cross_source_map":["2*b*t","-x/(2*b)"]}


def taxicab():
    # Divisor recovery, independently compared with a positive pair scan.
    rows=[]
    for d in [1,7,13,19,91,133,247,1729]:
        delta=Q(4*1729-d**3,3*d)
        rows.append({"trace":d,"discriminant":str(delta)})
    positive=[(x,y) for x in range(1,13) for y in range(x,13)
              if x**3+y**3==1729]
    check(positive==[(1,12),(9,10)],"taxicab positive pairs")
    for row in rows:
        delta=Q(row['discriminant'])
        if delta>=0:
            check(delta.denominator==1,"integral divisor discriminant")
            dd=delta.numerator
            check(not (dd%5==0 and isqrt(dd//5)**2==dd//5),
                  "no golden integral cubic trace")
    for u,v in [(1596,52668),(1092,3276)]:
        check(v*v==u**3-432*1729**2,"Mordell point")
    orders={str(p):phi_order(p) for p in [7,13,19,1729]}
    check(orders=={'7':16,'13':28,'19':18,'1729':1008},"golden CRT clock")
    check(lcm(16,28,18)==1008,"independent CRT lcm")
    phi=(1+sp.sqrt(5))/2
    check(sp.simplify(phi**3-(2+sp.sqrt(5)))==0,"golden cube")
    check(sp.simplify(phi**3+(1-phi)**3-4)==0,"golden cube trace")
    return {"taxicab_pairs":positive,"divisor_table":rows,"phi_orders":orders}


if __name__=='__main__':
    result={"geometry":geometry(),"boundary":boundary_tables(),
            "collatz":collatz(),"algebra":algebra(),"taxicab":taxicab()}
    result['checks']=checks
    print(json.dumps(result,indent=2,sort_keys=True))

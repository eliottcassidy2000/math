"""Exact image recurrence and precision obstruction for a native periodic clock."""
from collections import Counter
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
OTHER = HERE.parent/'overnight_rank_return_20261029'
sys.path.insert(0,str(OTHER))
from paley_pencil import diagonal, local_layers

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def mm(a,b,q=None):
    out = [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
    return out if q is None else [[x % q for x in row] for row in out]


def power(a,r,q):
    out = eye(len(a))
    while r:
        if r & 1:
            out = mm(out,a,q)
        a = mm(a,a,q)
        r //= 2
    return out


def val(a,p,cap=None):
    if a == 0:
        if cap is None:
            raise ValueError('unexpected exact zero')
        return cap
    v = 0
    while a % p == 0 and (cap is None or v < cap):
        a //= p
        v += 1
    return v


def eval_poly(cs,x,q):
    ans = 0
    for a in reversed(cs):
        ans = (ans*x+a) % q
    return ans


def image_log(h,r):
    return 12*h+5*max(h-2,0)+max(h-2-2*r,0)


def preperiod(h):
    return (max(h-4,0)+1)//2


def run():
    src = json.loads((OTHER/'periodic_defect.json').read_text())
    lift = json.loads((HERE/'periodic_frobenius_lift.json').read_text())
    P = src['recurrence_low_first']
    n, precision = 18, 12
    q = 3**precision
    B = src['transfer']
    Bt = list(map(list,zip(*B)))
    O = src['observable']
    C = [[int(j == i+1) for j in range(n)] for i in range(n-1)]+[[-a for a in P[:-1]]]
    need(len(P) == n+1 and P[-1] == 1, 'actual monic recurrence')
    need(mm(C,O) == mm(O,Bt), 'literal integral source/clock intertwiner')
    ds = diagonal(O)
    vals = [val(d,3) for d in ds]
    need(vals == [0]*12+[2]*6, 'integral observer lattice precision profile')
    need(val(P[0],3) == 2 and P[1] % 3 != 0, 'unique simple contracting Hensel root of valuation two')
    rho, roots = 0, [0]
    for h in range(1,precision):
        candidates = [rho+t*3**h for t in range(3) if eval_poly(P,rho+t*3**h,3**(h+1)) == 0]
        need(len(candidates) == 1, 'unique compatible Hensel digit')
        rho = candidates[0]
        roots.append(rho)
    need(val(rho,3) == 2, 'exact contracting root valuation')
    G = [0]*n
    G[-1] = 1
    for i in range(n-1,0,-1):
        G[i-1] = (P[i]+rho*G[i]) % q
    need((P[0]+rho*G[0]) % q == 0, 'Hensel polynomial division remainder')
    need(G[0] % 3 != 0, 'complementary recurrence factor has unit constant term')
    gr = eval_poly(G,rho,q)
    need(gr % 3 != 0, 'coprime factors provide integral spectral projector')
    I, E = eye(n), [[0]*n for _ in range(n)]
    for a in reversed(G):
        E = [[(x+a*I[i][j]) % q for j,x in enumerate(row)] for i,row in enumerate(mm(C,E,q))]
    E = [[x*pow(gr,-1,q) % q for x in row] for row in E]
    need(mm(E,E,q) == E, 'Hensel projector idempotent')
    need(mm(C,E,q) == [[rho*x % q for x in row] for row in E], 'contracting line carries multiplication by rho')
    need(E[0][0] % 3 == 1, 'projector has a primitive first-coordinate eigenline')
    EO = mm(E,O,q)
    need(min(val(x,3,precision) for row in EO for x in row) == 2, 'source lattice meets contracting line at valuation exactly two')
    need(all(EO[i][j] == pow(rho,i,q)*EO[0][j] % q for i in range(n) for j in range(35)), 'literal contracting eigenvector chart')
    initial = [[x % q] for x in src['values'][:n]]
    need(mm(O,[[x] for x in src['initial']]) == [[x] for x in src['values'][:n]], 'actual initial boundary weights retained')
    slow_initial = mm(E,initial,q)
    need(val(slow_initial[0][0],3,precision) == 4, 'actual initial trajectory has contracting amplitude of valuation four')
    U = [[(a-b) % q for a,b in zip(row,erow)] for row,erow in zip(I,E)]
    A = power(C,726,q)
    delta = [[(a-b) % q for a,b in zip(row,irow)] for row,irow in zip(A,I)]
    delta_unit = mm(delta,U,q)
    stable_return = mm(delta_unit,initial,q)
    need(min(val(x,3,precision) for row in delta_unit for x in row) == 1,
         '726th clock power is identity plus exactly three on the stable sector')
    need(min(val(x,3,precision) for row in stable_return for x in row) == 1,
         'actual stable trajectory has no extra first-return divisibility')
    checks = []
    for r in (0,1,2,3,4,7,11,13,487):
        shifted = mm(O,power(Bt,r,q),q)
        layers = local_layers(shifted,3,precision)
        v = [i for i,c in enumerate(layers) for _ in range(c)]
        for h in range(1,precision+1):
            observed = sum(max(h-e,0) for e in v)
            need(observed == image_log(h,r), 'independent local elimination agrees with all-order image formula')
            checks.append({'precision':h,'steps':r,'image_log_order':observed})
    exact_smith = []
    shifted = O
    for r in range(4):
        smith = diagonal(shifted)
        v = [val(a,3) for a in smith]
        need(v == [0]*12+[2]*5+[2+2*r], 'independent integer Smith at initial clock powers')
        exact_smith.append({'steps':r,'ternary_valuations':v})
        shifted = mm(C,shifted)
    # Verify the actual trajectory's transient through each tested precision
    # on the contracting component; the complementary factor is a unit.
    for h in range(1,precision+1):
        m = 3**h
        amp = slow_initial[0][0] % m
        t = 0
        while amp:
            amp = amp*rho % m
            t += 1
        need(t == preperiod(h), 'exact actual-source transient length')
    need(image_log(5,1) == 76 and image_log(5,487) == 75, 'mod243 image-cardinality conjugacy obstruction')
    need(all(image_log(h,1) == image_log(h,487) for h in range(1,5)), 'sharp first precision for this image-cardinality test')
    need(preperiod(7) == 2 and (preperiod(7)+486)//487 == 1, 'mod2187 actual-orbit preperiod obstruction')
    need(gcd(487,726) == 1 and 487 % 3 != 0, 'acceleration is a unit on every eventual cycle')
    orbit_lifts = []
    for h in range(1,7):
        tail, period = preperiod(h), 726*3**(h-1)
        def canonical(r):
            return r if r < tail else tail+(r-tail) % period
        images = [canonical(487*r) for r in range(tail+period)]
        need(len(set(images)) == len(images), 'source-fixed index map bijective through precision six')
        need(all(images[canonical(r+1)] == canonical(images[r]+487) for r in range(tail+period)),
             'literal recurrence-orbit conjugacy through precision six')
        orbit_lifts.append({'precision':h,'orbit_size':tail+period,'bijection':True})
    m, period = 3**7, 726*3**6
    x1 = mm(C,initial,m)
    xlate = mm(power(C,period+1,m),initial,m)
    need(x1 != xlate and x1[0] != xlate[0], 'actual scalar distinguishes transient point from its later cycle counterpart')
    need(mm(power(C,487,m),initial,m) == mm(power(C,487*(period+1),m),initial,m),
         'literal collision under the forced source-fixed Frobenius index map')
    # Full literal orbit audit at the first four precision levels. No
    # exponent fit: the all-height statement uses the binomial lift proof.
    literal_orbits = []
    for h in range(1,5):
        m = 3**h
        coeffs = [a % m for a in P[:-1]]
        state = tuple(v % m for v in src['values'][:n])
        seen, t = {}, 0
        while state not in seen:
            seen[state] = t
            state = state[1:]+((-sum(a*b for a,b in zip(coeffs,state)) % m),)
            t += 1
        need(seen[state] == preperiod(h) and t-seen[state] == 726*3**(h-1),
             'independent complete actual recurrence-state orbit')
        literal_orbits.append({'precision':h,'preperiod':seen[state],'period':t-seen[state]})

    Phi = lift['frobenius_in_observation_coordinates']
    columns = [tuple(x % 3 for x in col) for col in zip(*O[:12])]
    images = [tuple(row[0] for row in mm(Phi,[[x] for x in col],3)) for col in columns]
    failures = [i for i,v in enumerate(images) if v not in columns]
    need(len(set(columns)) == 22 and len(failures) == 27, 'full boundary-response permutation obstruction')
    need(Counter(images) != Counter(columns), 'no permutation of all35 boundary states induces the Frobenius lift')
    T = lift['clock']
    need(mm(lift['frobenius'],T,3) != mm(T,lift['frobenius'],3), 'no pure bulk word can equal the lifted Frobenius')
    return {'status':'PROVED all-height image/initial-transient formulas with FINITE-EXACT lattice certificates',
            'source_sha256_canonical':hashlib.sha256(json.dumps(src,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'precision':precision,'recurrence':P,'clock':C,'observation_smith':ds,
            'hensel_root_residues':roots,'contracting_root_modulus':q,'contracting_root':rho,
            'unit_factor_low_first_modulus':G,'contracting_projector':E,
            'contracting_source_first_row':EO[0],'contracting_initial_amplitude':slow_initial[0][0],
            'source_line_valuation':2,'initial_line_valuation':4,
            'image_log_formula':'12*h + 5*max(h-2,0) + max(h-2-2*r,0)',
            'actual_initial_preperiod_formula':'ceil(max(h-4,0)/2)',
            'actual_eventual_period_formula':'726*3^(h-1)',
            'stable_first_return_vector':stable_return,'literal_initial_orbit_checks':literal_orbits,
            'source_fixed_orbit_lifts':orbit_lifts,
            'first_forced_orbit_collision':{'precision':7,'r':1,'s':period+1,
                'b_r':x1[0][0],'b_s':xlate[0][0],'equal_after_multiplying_indices_by':487},
            'exact_integer_smith_checks':exact_smith,'modular_checks':checks,
            'boundary_permutation':{'distinct_columns':22,'outside_image_indices':failures,
                'witness_index':failures[0],'source_column':columns[failures[0]],'moved_column':images[failures[0]]},
            'gates':GATES}


if __name__ == '__main__':
    out = run()
    (HERE/'periodic_precision_obstruction.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('gates','image_log_formula','actual_initial_preperiod_formula','source_line_valuation','initial_line_valuation')}))

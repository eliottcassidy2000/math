"""Independent literal-charge / symbolic-window / Bernstein referee."""
from pathlib import Path
from fractions import Fraction
from math import comb, factorial, gcd
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
STEM = "continuing1_20260906_laurent_cone_separator"
PINS = {
    ".py": "5bd98bc1fe072fcfc88b34af46b2310e1a4e184dcefdff6c5a508a820d82d00b",
    ".out": "207b86469f1b5e8b5bb18a5cd1cd2a931c4c0c596a5fe66ba354923fe3069fee",
    "_certificate.json": "4b1ee5770b484e4164e692fbf2934f4099800b0b85d379ce01d8afef71040cc0",
}
gates = 0

def check(ok, label):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(label)

for ending, pin in PINS.items():
    path = HERE / (STEM + ending)
    if not path.exists() and ending == ".out":
        path = HERE.parent / "05-knowledge/results" / (STEM + ending)
    check(hashlib.sha256(path.read_bytes()).hexdigest() == pin, "frozen " + ending)
cert = json.loads((HERE / (STEM + "_certificate.json")).read_text())
s, u, v, z = S.symbols("s u v z")

def poly(x):
    return S.Poly(S.expand(x), s, domain=S.QQ)

def rem(x, p):
    return poly(x).rem(p).as_expr()

def coeffs(x, size):
    q = poly(x)
    return [q.nth(j) for j in range(size)]

def listed(xs):
    return [S.Rational(x) for x in xs]

def ell(x, p, y):
    return sum(a*b for a, b in zip(coeffs(rem(x, p), p.degree()), y))

def falling(n, r):
    return factorial(n) // factorial(n-r) if n >= r >= 0 else 0

def sq_coeff(coefs, degree):
    return S.expand(sum(coefs[i] * coefs[degree-i]
                        for i in range(len(coefs))
                        if 0 <= degree-i < len(coefs)))

def product_coeff(a, b, degree):
    return S.expand(sum(a[i] * b[degree-i] for i in range(len(a))
                        if 0 <= degree-i < len(b)))

def ucoeffs(x):
    q = S.Poly(S.expand(x), u)
    return [q.nth(i) for i in range(q.degree()+1)]

def bernstein_sign(x, left, right, sign):
    """A different exact enclosure: convex hull of Bernstein coefficients."""
    pz = S.Poly(poly(x).as_expr().subs(s, left+(right-left)*z), z)
    n = pz.degree()
    bs = [sum(pz.nth(j)*S.Rational(comb(k,j),comb(n,j))
              for j in range(k+1)) for k in range(n+1)]
    check(all(sign*b > 0 for b in bs), "strict Bernstein enclosure")
    return bs

data = {}
literal_rows = 0
window_count = 0
for h in (1, 2, 3):
    g, q, k = 3*h+2, h+1, 2*h+1
    charges = (-6*h-3, 1, 3*h+3)
    check(gcd(g, 6*h+3) == 1, "return clock primitive")
    check(all(x % g == 1 for x in charges), "all charges one modulo g")
    fibers = {}
    # Literal charges and multinomial weights; no prescribed return parameter j.
    for mass in range(1, 2*g+1):
        row = []
        for na in range(mass+1):
            for nb in range(mass-na+1):
                nc = mass-na-nb
                if na*charges[0]+nb*charges[1]+nc*charges[2] == 0:
                    weight = factorial(mass)//(factorial(na)*factorial(nb)*factorial(nc))
                    row.append((na, nb, nc, weight))
        fibers[mass] = row
        check(bool(row) == (mass % g == 0), "exact first-return clock")
        literal_rows += len(row)
    first = 0
    target = 0
    for na, nb, nc, weight in fibers[g]:
        j = na-1
        check((na, nb, nc) == (1+j, 3*h-3*j, 1+2*j), "first literal counts")
        first += weight*(-s)**j
    for na, nb, nc, weight in fibers[2*g]:
        j = na-2
        check((na, nb, nc) == (2+j, 6*h-3*j, 2+2*j), "doubled literal counts")
        target += weight*(-1)**j*s**(j+1) if j >= 0 else -S.Integer(weight)
    first, target = S.expand(first), S.expand(target)
    p = poly(first).monic()
    # Build the COMPLETE path coefficients by counting binary words, not binom(q+2i,...).
    def paths(east, north):
        if east < 0 or north < 0:
            return 0
        dp = [[0]*(north+1) for _ in range(east+1)]
        dp[0][0] = 1
        for x in range(east+1):
            for y in range(north+1):
                if x or y:
                    dp[x][y] = (dp[x-1][y] if x else 0)+(dp[x][y-1] if y else 0)
        return dp[east][north]
    KB = KC = KD = 0
    for j in range(-1, h+1):
        # Reverse the full Laurent beta index j into i=h-j.
        i = h-j
        KB += (-1)**(q-i)*paths(3*h-3*j, 1+j)*s**(q-i)*u**(2*i)
    for j in range(-1, h):
        i = h-1-j
        KC += (-1)**(q-1-i)*paths(3*h-1-3*j,1+j)*s**(q-1-i)*u**(2*i)
        KD += (-1)**(q-1-i)*paths(3*h-2-3*j,1+j)*s**(q-1-i)*u**(2*i)
    path_beta = S.Poly(KB.subs(s,1),u)
    beta_v = S.Poly(sum(path_beta.nth(2*i)*v**i for i in range(q+1)),v)
    check(beta_v.count_roots(0,S.oo)==q,"alternating displayed beta has positive roots")
    if h == 3:
        check([beta_v.nth(i) for i in range(q+1)]==listed(cert["cubic"]["beta"]),"complete original beta kernel")
        for kernel,key in ((KC,"C"),(KD,"D")):
            raw = S.Poly(kernel.subs(s,1),u)
            check([raw.nth(2*i) for i in range(q)]==listed(cert["cubic"][key]),"complete original contiguous kernel")
    H = S.expand((1+u)**g*KB)
    HC = S.expand((1+u)**g*KC)
    HD = S.expand((1+u)**g*KD)
    hs, cs, ds = map(ucoeffs, (H, HC, HD))
    check(hs[k] == S.expand(-s*first), "exact selected coefficient")
    check(S.expand(sq_coeff(hs,2*k)/s-2*product_coeff(cs,ds,2*k-2)-target) == 0,
          "literal carried midpoint, correct relative factor")
    check(poly(target).degree() == 2*h+1 and poly(target).nth(0) < 0, "retained inverse carry")
    rays = []
    for r in range(k):
        # Symbolic differentiated polynomial, then ordinary coefficient extraction.
        deriv = ucoeffs(S.diff(H,u,r))
        jr = S.cancel(sq_coeff(deriv,2*k-2*r)/s)
        check(S.denom(jr) == 1 or not S.denom(jr).has(s), "exact polynomial division")
        rays.append(S.expand(jr))
    n = len(hs)-1
    # Independently solve the paired-grid interpolation; no constructive cone recurrence.
    grid = S.Matrix([[S.prod((k-i)**2-d*d for i in range(r)) for r in range(k)]
                     for d in range(1,k+1)])
    check(grid.det() != 0, "paired grid independent")
    inverse = grid.inv()
    for r in range(k):
        for right in range(n-k):
            w = S.Matrix([S.prod((k-i)**2-d*d for i in range(r)) *
                          S.prod((n-k-i)**2-d*d for i in range(right))
                          for d in range(1,k+1)])
            coords = inverse*w
            check(all(c >= 0 for c in coords) and any(c > 0 for c in coords), "window cone coordinates")
            # Literal reverse / derivative / reverse / derivative, independent of grid.
            rev = S.expand(sum(a*u**(n-i) for i,a in enumerate(hs)))
            dr = S.diff(rev,u,right)
            back = S.Poly(dr,u)
            right_poly = S.expand(sum(back.nth(i)*u**(n-right-i) for i in range(n-right+1)))
            win = ucoeffs(S.diff(right_poly,u,r))
            response = S.cancel(sq_coeff(win,2*k-2*r)/s)
            check(rem(response-sum(c*ray for c,ray in zip(coords,rays)),p) == 0,
                  "literal two-sided window in quotient cone")
            window_count += 1
    data[h] = (p, first, target, rays, H, HC, HD)
    print("LITERAL", h, charges, "return masses", g, 2*g, "rows", len(fibers[g]), len(fibers[2*g]),
          "rays", k, "windows", k*(n-k))

p, first, target, rays, *_ = data[1]
check(rem(target-S.Rational(6305,1714)*rays[0],p) == 0,"linear positive repair")
p, first, target, rays, *_ = data[2]
quad = cert["quadratic"]
y = listed(quad["dual"])
check(coeffs(p.as_expr(),3) == listed(quad["P"]), "quadratic first")
check(coeffs(rem(target,p),2) == listed(quad["T"]), "quadratic target")
check([ell(j,p,y) for j in rays] == listed(quad["values"]), "quadratic dual rays")
check(ell(target,p,y) == S.Rational(quad["target"]) < 0,"quadratic separation")
repair = sum(c*s**i for i,c in enumerate(listed(quad["repair"])))
check(all(c>0 for c in listed(quad["repair"])) and rem(target-repair*rays[0],p)==0,"quadratic positive phase repair")
for j, recorded in zip(rays,quad["residues"]):
    check(coeffs(rem(j,p),2)==listed(recorded),"quadratic literal ray")

p, first, target, rays, H, HC, HD = data[3]
c = cert["cubic"]
y = listed(c["dual"])
check(coeffs(p.as_expr(),4)==listed(c["P"]),"cubic first monic")
check(coeffs(first,4)==listed(c["full_first"]),"cubic first original scale")
check(coeffs(target,8)==listed(c["full_target"]),"cubic doubled original scale")
check(coeffs(rem(target,p),3)==listed(c["T"]),"cubic target remainder")
check(ell(target,p,y)==S.Rational(c["target"])<0,"cubic target separated")
check(ell(rays[0],p,y)==0 and ell(s*rays[6],p,y)==0,"both exposed boundary rays")
sinverse = S.invert(poly(s),p).as_expr()
check(rem(s*sinverse,p)==1,"inverse phase exists")
for r,j in enumerate(rays):
    check(coeffs(rem(j,p),3)==listed(c["residues"][r]),"cubic literal ray")
    levels = {-1:ell(j*sinverse,p,y),0:ell(j,p,y),1:ell(s*j,p,y),2:ell(s*s*j,p,y)}
    for d,e in levels.items():
        check(e==S.Rational(c["moments"][r][str(d)]),"exact four levels")
    check(levels[0]>=0 and levels[1]>=0,"central values nonnegative")
    check(levels[-1]>3*levels[0] and levels[2]>levels[1],"outer differences point outward")
    # Fresh quotient-ring matrix-power path, positive and inverse powers.
    mult = S.Matrix.hstack(*[S.Matrix(coeffs(rem(s**(i+1),p),3)) for i in range(3)])
    rvec = S.Matrix(coeffs(rem(j,p),3))
    dual = S.Matrix([y])
    vals = {d:(dual*(mult**d)*rvec)[0] for d in range(-9,10)}
    for d in range(-9,7):
        check(vals[d+3]==28*vals[d+2]-14*vals[d+1]+vals[d]/3,"companion recurrence")
    check(all(e>=0 for e in vals.values()),"independent bounded companion powers")
    print("FOUR_LEVELS", r, *(str(levels[d]) for d in (-1,0,1,2)))

N = S.expand(y[0]*(s*s-28*s+14)+y[1]*(s-28)+y[2])
check(coeffs(N,3)==listed(c["N"]),"Lagrange dual numerator")
multiplier = sum(x*s**i for i,x in enumerate(listed(c["signed_multiplier"])))
check(rem(target-multiplier*rays[0],p)==0,"rootwise-positive signed repair exact")
check(poly(multiplier).nth(2)<0,"repair lies outside coefficientwise positive class")
intervals = [tuple(map(S.Rational,lr)) for lr in c["intervals"]]
check(all(0 < a < b for a,b in intervals),"positive ordered intervals")
check(intervals[0][1] < intervals[1][0] < intervals[1][1] < intervals[2][0],"disjoint exhaustive cubic intervals")
check(S.Rational(1,3)<intervals[1][0]<intervals[1][1]<1,"middle-root numerical-free bounds")
for index,(left,right) in enumerate(intervals):
    check(p.eval(left)*p.eval(right)<0,"rational endpoint sign change")
    # Degree three and three disjoint sign changes imply exactly one simple root in each.
    bernstein_sign(p.diff().as_expr(),left,right,1 if index!=1 else -1)
    bernstein_sign(N,left,right,-1)
    for j in rays:
        bernstein_sign(rem(j,p),left,right,-1)
    bernstein_sign(rem(target,p),left,right,-1)
    bernstein_sign(multiplier,left,right,1)
print("SPECTRAL roots positive/simple; beta in (1/3,1); dual signs (-,+,-); all seven rays negative")
print("ALL-INTEGER theorem: normalized spectral sequence has strict positive second difference; four levels confine its integer minimum to {0,1}")

# Original zero-preserving lowering and its failed deletion, without producer formulas.
g, k = 11, 7
A = S.cancel(H/(1+u))
Blow = S.Rational(2,3)*u*S.cancel(HC/(1+u))
low = S.expand(A+Blow)
ac,bc,lc = map(ucoeffs,(A,Blow,low))
check(S.expand(k*S.Poly(H,u).nth(k)-g*S.Poly(low,u).nth(k-1))==0,"lowered selected-zero identity")
parts = [S.cancel(sq_coeff(ac,2*k-2)/s),S.cancel(2*product_coeff(ac,bc,2*k-2)/s),S.cancel(sq_coeff(bc,2*k-2)/s)]
partvals = [ell(x,p,y) for x in parts]
check(partvals+[sum(partvals)]==listed(cert["lowered_parts"]["values"]),"mixed pieces and complete lowering")
for x,recorded in zip(parts,cert["lowered_parts"]["residues"]):
    check(coeffs(rem(x,p),3)==listed(recorded),"lowered original polynomial remainder")
selected_a = rem(S.Poly(A,u).nth(k-1),p)
check(coeffs(selected_a,3)==listed(cert["lowered_parts"]["selected_A_remainder"]) and selected_a!=0,"deleting lowering piece loses zero")
hostile = (1-2*u)**2
check(S.expand((1+u)*hostile).coeff(u,2)==0 and S.expand(hostile**2).coeff(u,2)==24,"deleted-piece sign hostile")

print("SOURCE universe h=1,2,3; literal return rows",literal_rows,"coupled windows",window_count)
print("PASS",gates,"always-active gates; no producer import; no numeric root approximations")

"""Independent native-wall/grid referee; no producer imports."""
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
BASE = Path(__file__).resolve().parent
STEM = "continuing5_20260906_lrc_extremal_location"
PROFILE = "overnight12_20260906_lrc_decoder_descent_inherited_profiles.json"
GATES = 0


def gate(ok, name):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(name)


def locate(name):
    for path in (BASE/name, BASE.parent/"05-knowledge"/"results"/name,
                 Path("C:/w/s0905/04-computation")/name):
        if path.exists():
            return path
    raise FileNotFoundError(name)


for suffix, pin in [
    (".py", "6526dd4306c31861792c5e938b37e498675cb1f39f9a05b8158d2f0940d44659"),
    (".out", "954c713447f9ff4ba162f0318195d12b5a983c6e646ca54821864a661a625ed8"),
]:
    gate(hashlib.sha256(locate(STEM+suffix).read_bytes()).hexdigest() == pin, "frozen producer "+suffix)
raw = locate(PROFILE).read_bytes()
gate(hashlib.sha256(raw).hexdigest() == "935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f", "inherited profile pin")
profiles = {int(k): {(int(c),tuple(w)) for c,w in obj["profiles"]}
            for k,obj in json.loads(raw)["levels"].items()}


def frac(x):
    return x-x.numerator//x.denominator


def ceil(x):
    return -((-x.numerator)//x.denominator)


def danger(x, speed):
    r = (speed*x.numerator) % x.denominator
    return 14*min(r,x.denominator-r) < x.denominator


def physical_walls(speed):
    return {frac(F(14*k+sign,14*speed)) for k in range(speed) for sign in (-1,1)}


def grid_hits(n, alpha, speeds=(3,308)):
    """Literal strict predicates, integer arithmetic at each actual grid point."""
    denominator = n*alpha.denominator
    bits = 0
    value = alpha.numerator
    for j in range(n):
        ok = True
        for u in speeds:
            residue = u*value % denominator
            if 14*min(residue,denominator-residue) >= denominator:
                ok = False
                break
        if ok:
            bits |= 1 << j
        value += alpha.denominator
    return bits


p,q,t,e = 3,308,16704,4
n = t//e
word = (12,16,72,58,64,9,9)
gate(n == 4176 and gcd(p,q) == 1 and p+q == 311, "displayed selected pair")
gate(all(311%d for d in range(2,18)) and 311%3 == 2 and 311 <= 356, "literal strict inert atlas entry")
gate(all(t%d == 0 for d in word) and reduce(gcd,word) == 1, "sheet divisibilities and primitive word")
gate((e*gcd(n,p),e*gcd(n,q)) == word[:2], "forced margins 12,16")
subsets = 0
for k in range(1,7):
    for chosen in combinations(range(7),k):
        c = reduce(gcd,(word[i] for i in chosen))
        complement = tuple(sorted(gcd(c,word[j]) for j in range(7) if j not in chosen))
        gate((c,complement) in profiles[7-k], "every projected proper-subset word")
        subsets += 1
gate(subsets == 126, "complete word scope")
costs = tuple(d*ceil(F(t,7*d)) for d in word)
excess = sum(costs)-t
gate(excess == 188, "literal word excess")

# Independent spatial construction: sweep all raw danger boundaries, rather
# than intersecting pairs of intervals. Test the literal predicates on cells.
spatial = sorted(physical_walls(p) | physical_walls(q))
gate(len(spatial) == 622, "complete raw spatial boundary set")
components = []
for i,lo in enumerate(spatial):
    hi = spatial[i+1] if i+1 < len(spatial) else spatial[0]+1
    mid = (lo+hi)/2
    if danger(mid,p) and danger(mid,q):
        gate(not (danger(lo,p) and danger(lo,q)), "open component excludes lower wall")
        gate(not (danger(hi,p) and danger(hi,q)), "open component excludes upper wall")
        components.append((lo,hi))
gate(len(components) == 45, "45 exact open intersection components")
lengths = Counter(hi-lo for lo,hi in components)
gate(lengths == {F(1,7*q):43,F(1,14*q):2}, "43 full clipped widths and two half widths")
measure = sum((hi-lo for lo,hi in components),F(0))
gate(measure == F(1,49), "exact measure")
separate = sum(ceil(n*(hi-lo))-1 for lo,hi in components)
gate(separate == 43 and e*separate == 172, "separate component lower credit")
projected = {frac(n*x) for component in components for x in component}
gate(len(projected) == 90, "90 intersection endpoint residues")

# The independent uniform certificate uses the LARGER raw phase-wall set
# of the two danger predicates, not the producer's 90 intersection walls.
native_walls = sorted({frac(n*x) for x in spatial})
gate(len(native_walls) == 156, "156 native-predicate phase walls")
gate(projected <= set(native_walls), "intersection changes among native walls")
cache = {}
def state(alpha):
    if alpha not in cache:
        cache[alpha] = grid_hits(n,alpha)
    return cache[alpha]

wall_counts = []
cell_counts = []
active = set()
for i,w in enumerate(native_walls):
    prev = native_walls[i-1] if i else native_walls[-1]-1
    nxt = native_walls[i+1] if i+1 < len(native_walls) else native_walls[0]+1
    left,right = (prev+w)/2,(w+nxt)/2
    L,W,R = state(left),state(w),state(right)
    gate(W & ~(L|R) == 0, "strict wall creates no spurious dangerous point")
    if L != W or R != W:
        active.add(w)
    wall_counts.append(W.bit_count())
    cell_counts.append(R.bit_count())
gate(active == projected, "literal point changes recover exactly the 90 active walls")
gate(min(wall_counts+cell_counts) == 84 and max(wall_counts+cell_counts) == 87,
     "exact uniform min84 max87 on every translate")
gate(state(F(53,539)).bit_count() == 84, "quoted minimum attainment")
gate(state(F(74,77)).bit_count() == 87, "quoted maximum attainment")
gate(state(F(0)).bit_count() >= 84, "origin translate retained")

# Separate path for representative phases: direct membership in
# the spatial cells, independent of the native modular danger predicates.
for alpha in (F(0),F(53,539),F(74,77),F(1,2),F(1,7)):
    count = 0
    for j in range(n):
        xx = frac((alpha+j)/n)
        count += any(lo < xx < hi or lo < xx+1 < hi for lo,hi in components)
    gate(count == state(alpha).bit_count(), "literal spatial membership versus modular grid")

# Exact finite controls of the e-fold lifting identity; d=e*k with (k,n)=1.
lift_checks = 0
for k in (1,5,n+1,2*n-1):
    d = e*k
    gate(gcd(d,t) == e and gcd(k,n) == 1, "actual sheet multiplicity")
    gate((gcd(t,d*p),gcd(t,d*q)) == (12,16), "physical selected pair retains forced margins")
    for alpha in (F(0),F(1,3),F(53,539),F(74,77)):
        actual = grid_hits(t,alpha,(d*p,d*q)).bit_count()
        expected = e*grid_hits(n,k*alpha).bit_count()
        gate(actual == expected and actual >= 336, "literal full-grid e-fold overlap count")
        lift_checks += 1
gate(gcd(12,t) != e, "wrong assumed cofactor gcd is excluded")
credit = e*min(wall_counts+cell_counts)
gate(credit == 336 and credit-excess == 148, "conditional weak-safe lift count")

print("Displayed tuple: t16704, e4, primitive pair(3,308), word(12,16,72,58,64,9,9).")
print("126 projected words PASS; literal marginal costs "+str(costs)+"; E="+str(excess)+".")
print("Spatial native-wall sweep: 622 boundaries,45 components,measure1/49; separate credit43*4=172.")
print("Independent phase universe:156 raw walls +156 open cells; exactly90 active intersection walls.")
print("Wall count histogram: "+str(dict(sorted(Counter(wall_counts).items()))))
print("Open-cell count histogram: "+str(dict(sorted(Counter(cell_counts).items()))))
print("Uniform count [84,87]; min at53/539; max at74/77; actual pair credit336>E188.")
print("16 literal full-grid lift controls PASS; conditional weak-safe lower count148.")
print("Only the displayed selected-pair/word control is eliminated; clock16704 is not globally excluded.")
print("PASS: "+str(GATES)+" always-active exact gates.")

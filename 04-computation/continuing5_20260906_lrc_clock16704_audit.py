"""Independent one-clock referee: native walls, knapsack, full word masks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, comb
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
BASE = Path(__file__).resolve().parent
STEM = "continuing5_20260906_lrc_clock16704"
GATES = 0
T = 16704


def gate(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def locate(name):
    for p in (BASE/name, BASE.parent/"05-knowledge"/"results"/name,
              Path("C:/w/s0905/04-computation")/name,
              Path("C:/w/s0905/05-knowledge/results")/name):
        if p.exists():
            return p
    raise FileNotFoundError(name)


def digest(x):
    return hashlib.sha256(json.dumps(x,separators=(",",":")).encode()).hexdigest()


def ceiling(a,b):
    return (a+b-1)//b


def frac(x):
    return x-x.numerator//x.denominator


profile_raw = locate("overnight12_20260906_lrc_decoder_descent_inherited_profiles.json").read_bytes()
gate(hashlib.sha256(profile_raw).hexdigest() == "935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f", "profile supplier pin")
profile = json.loads(profile_raw)
profiles = {int(k): {(c,tuple(word)) for c,word in obj["profiles"]}
            for k,obj in profile["levels"].items()}
domain = tuple(d for d in profile["levels"]["6"]["gcds"] if T%d == 0)
gate(len(domain) == 18, "complete allowed divisor domain")
caps = (90,30,9,4,2,1,1)
capacity = {d:sum(d<=cap for cap in caps) for d in domain}
native_cost = {d:d*ceiling(T,7*d) for d in domain}
maxima = [json.loads(line[8:]) for line in locate("third_20260906_grid_full_words.out").read_text().splitlines()
          if line.startswith("MAXIMUM ")]
gate(digest([[t,E] for t,E,*_ in maxima]) == "ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca", "inherited maximum-table semantics")
M = next(E for t,E,*_ in maxima if t == T)
gate(M == 188, "inherited unconditional full-word maximum")
gate(hashlib.sha256(locate(STEM+".py").read_bytes()).hexdigest() == "2d12351426bdf67c0801bdc46b16c50e35f6f6600b49e8e5323d4cc975bc2d3a", "frozen producer source pin")
producer_raw = locate(STEM+".out").read_bytes()
gate(hashlib.sha256(producer_raw).hexdigest() == "6d6acf82d798074bd5fac208f79c6b91bf7b9df9f425f1204d9e56324f5cac76", "frozen producer output pin")
producer_lines = producer_raw.decode().splitlines()
produced_candidates = json.loads(next(line[len("CANDIDATES "):] for line in producer_lines if line.startswith("CANDIDATES ")))
produced_prefix = json.loads(next(line[len("LOCATED_PREFIX "):] for line in producer_lines if line.startswith("LOCATED_PREFIX ")))


def pair_cost(a,b):
    """Five-slot bounded knapsack in native ceiling costs, not a sorted bag."""
    used = Counter((a,b))
    if any(d not in capacity or count > capacity[d] for d,count in used.items()):
        return None
    dp = [0]+[None]*5
    for d in domain:
        available = capacity[d]-used[d]
        new = [None]*6
        for occupied,value in enumerate(dp):
            if value is None:
                continue
            for count in range(min(available,5-occupied)+1):
                slot = occupied+count
                val = value+count*native_cost[d]
                if new[slot] is None or val > new[slot]:
                    new[slot] = val
        dp = new
    return native_cost[a]+native_cost[b]+dp[5]-T


def allowed_sum(s):
    # Enumerate prime factors independently from a precomputed small-prime list.
    primes = [p for p in range(2,s+1) if all(p%d for d in range(2,int(p**0.5)+1))]
    for p in primes:
        exponent = 0
        while s%p == 0:
            exponent += 1
            s //= p
        if exponent and (p%3 != 2 or exponent > 2):
            return False
    return s == 1


def spatial(p,q):
    """All raw integer-coordinate walls and a literal midpoint predicate."""
    L = 14*p*q
    cuts = {0,L}
    for v in (p,q):
        width = L//(14*v)
        for k in range(v):
            for sign in (-1,1):
                cuts.add((14*k+sign)*width % L)
    ordered = sorted(cuts)
    intervals = []
    den = 2*L
    for lo,hi in zip(ordered,ordered[1:]):
        num = lo+hi
        if all(14*min(v*num % den,(-v*num) % den) < den for v in (p,q)):
            intervals.append((lo,hi))
    gate(intervals[0][0] == 0 and intervals[-1][1] == L, "literal origin component")
    intervals = [(intervals[-1][0]-L,intervals[0][1])] + intervals[1:-1]
    return L,intervals


atlas = []
geometry = {}
for total in range(3,357):
    if not allowed_sum(total):
        continue
    for p in range(1,(total+1)//2):
        q = total-p
        if p < q and gcd(p,q) == 1:
            atlas.append((p,q))
            geometry[p,q] = spatial(p,q)
gate(len(atlas) == 5855 and len({p+q for p,q in atlas}) == 94, "complete strict atlas")
divisors = [e for e in range(1,7) if T%e == 0]
gate(divisors == [1,2,3,4,6], "all inherited sheet divisors")
rows = []
rejects = Counter()
pair_cache = {}
for e in divisors:
    n = T//e
    for p,q in atlas:
        L,I = geometry[p,q]
        separate = e*sum(ceiling(n*(b-a),L)-1 for a,b in I)
        if separate > M:
            rejects["component_cost"] += 1
            continue
        a,b = e*gcd(n,p),e*gcd(n,q)
        key = tuple(sorted((a,b)))
        if key not in pair_cache:
            pair_cache[key] = pair_cost(a,b)
        Epair = pair_cache[key]
        if Epair is None:
            rejects["forced_margin"] += 1
            continue
        if separate > min(M,Epair):
            rejects["paired_cost"] += 1
            continue
        rows.append([e,p,q,a,b,separate,Epair])
rows.sort()
gate(len(rows) == 988 and rows == produced_candidates, "all988 retained rows match independently")
gate(rejects == {"component_cost":28060,"forced_margin":227}, "complete exclusion accounting")
gate(len(rows)+sum(rejects.values()) == 5855*5, "no omitted pair-sheet input")


def literal_grid(n,alpha,p,q):
    den = n*alpha.denominator
    num = alpha.numerator
    count = 0
    for _ in range(n):
        count += all(14*min(v*num % den,(-v*num) % den) < den for v in (p,q))
        num += alpha.denominator
    return count


def event_sweep(n,p,q):
    """Enter/exit sweep; no component floor/ceiling point-count formula."""
    L,I = geometry[p,q]
    enters,exits = Counter(),Counter()
    for a,b in I:
        enters[n*a % L] += 1
        exits[n*b % L] += 1
    walls = sorted(enters.keys() | exits.keys())
    # The cell immediately preceding the first wall, in an unwrapped period.
    alpha = F(walls[-1]+walls[0]+L,2*L)
    current = literal_grid(n,alpha,p,q)
    initial = current
    values = [current]
    for w in walls:
        at_wall = current-exits[w]
        current = at_wall+enters[w]
        values.extend((at_wall,current))
    gate(current == initial, "complete enter-exit period")
    return min(values),max(values),len(walls)


gate(len(produced_prefix) == 107, "declared bounded located prefix")
for index,record in enumerate(produced_prefix):
    row = rows[index]
    gate(record[:7] == row, "prefix ordering agrees with full candidate set")
    e,p,q,*_ = row
    lo,hi,walls = event_sweep(T//e,p,q)
    gate((e*lo,lo,hi,walls) == (record[7],record[8],record[10],record[12]), "independent complete located profile")
    gate(literal_grid(T//e,F(record[9]),p,q) == lo, "literal minimum owner")
    gate(literal_grid(T//e,F(record[11]),p,q) == hi, "literal maximum owner")
    gate((e*lo <= min(M,row[6])) == (index == 106), "first located survivor only at index106")

# Selected hostile gets an additional complete RAW native-wall partition.
row = rows[106]
gate(row == [4,23,323,4,4,180,207], "selected old-test survivor")
e,p,q = row[:3]
n = T//e
native = set()
for v in (p,q):
    for k in range(v):
        for sign in (-1,1):
            native.add(frac(F(n*(14*k+sign),14*v)))
native = sorted(native)
values = []
for i,w in enumerate(native):
    nxt = native[i+1] if i+1 < len(native) else native[0]+1
    values.append(literal_grid(n,w,p,q))
    values.append(literal_grid(n,(w+nxt)/2,p,q))
gate(len(native) == 692 and min(values) == 45 and max(values) == 92, "all692 native walls and692 cells")
gate(literal_grid(n,F(0),p,q) == 45, "all separate minima attained at origin")
gate(e*min(values) == 180, "location alone leaves the pair")

# Complete forced-word conditioning, using a subset-gcd bitmask recurrence.
def valid_word(word):
    gs = [0]*128
    for mask in range(1,128):
        bit = mask & -mask
        gs[mask] = gcd(gs[mask^bit],word[bit.bit_length()-1])
    if gs[127] != 1:
        return False
    for mask in range(1,127):
        c = gs[mask]
        remaining = tuple(sorted(gcd(c,word[j]) for j in range(7) if not(mask >> j & 1)))
        if (c,remaining) not in profiles[7-mask.bit_count()]:
            return False
    return True


tested = valid = 0
best = -1
owners = []
cost_hist = Counter()
for tail in combinations_with_replacement(domain,5):
    tested += 1
    word = (4,4)+tail
    if not valid_word(word):
        continue
    valid += 1
    E = sum(native_cost[d] for d in word)-T
    cost_hist[E] += 1
    if E > best:
        best,owners = E,[word]
    elif E == best:
        owners.append(word)
gate(tested == comb(22,5) == 26334 and valid == 2422, "complete conditioned multiset universe")
gate(best == 134 and (4,4,3,9,36,58,64) in owners, "full-word conditional maximum and owner")
gate(not any(E >= 180 for E in cost_hist), "no compatible full word attains necessary cost180")
gate(180-best == 46, "conditional actual weak-safe surplus")

print("Independent atlas5855 x sheetdivisors5: candidates988, complete row-by-row agreement.")
print("Candidates SHA256 "+digest(rows))
print("Prefix107: native enter/exit sweeps and literal extremum owners agree; first survivor index106.")
print("Selected pair(e,p,q)=(4,23,323), forced(4,4): sep180, located180, relaxed Epair207, inherited M188.")
print("Selected native partition:692 walls+692 cells, min45 at0, max92; location alone supplies no gain.")
print("Conditioning:26334 five-slot completions,2422 valid full words,max134 at(4,4,3,9,36,58,64).")
print("Conditional credit180 exceeds134 by46; no actual or full-conditioned survivor is claimed.")
print("Candidates107..987 located profiles not audited; whole-clock elimination remains unproved here.")
print("PASS: "+str(GATES)+" always-active exact gates.")

"""Independent complete-clock referee: native walls, knapsack, subset gcds."""
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
COMPLETE = STEM+"_complete"
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
for suffix,expected in ((".py","2c1e624743c6357257558ed8b8f291667577aad02e36a567b71f4d3ee1eec05f"),
                        (".out","33820ee3775baaedb0b8516b1a837b184b6ffdd279266614e479079f215e49fa"),
                        ("_certificate.json","8ab0b9fe8791707ac780f32030103eeec71de6830613619c8333672563a15345")):
    gate(hashlib.sha256(locate(COMPLETE+suffix).read_bytes()).hexdigest()==expected,"complete producer frozen "+suffix)
certificate=json.loads(locate(COMPLETE+"_certificate.json").read_bytes())


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
gate(rows==certificate["candidates"],"complete certificate retains exactly the independently generated universe")
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



full_words=[]
conditional={tuple(key):None for key in combinations_with_replacement(domain,2)}
visits=0
for word in combinations_with_replacement(domain,7):
    visits+=1
    if not valid_word(word):
        continue
    E=sum(native_cost[d] for d in word)-T
    full_words.append([list(word),E])
    pairs={tuple(sorted((word[i],word[j]))) for i in range(7) for j in range(i+1,7)}
    for key in pairs:
        old=conditional[key]
        if old is None or E>old[0]:
            conditional[key]=[E,list(word)]
gate(visits==346104 and len(full_words)==19073, "complete full-word universe")
gate(sum(x is not None for x in conditional.values())==161,"all realized forced pairs")
gate(max(E for w,E in full_words)==188,"independent unconditional optimum")
gate(digest([w for w,E in full_words])==certificate["valid_word_sha256"]=="e585b7c1aa6f5bfd55abc2d3caaf2b245789a9ace580d8abd46a629d0fa49eef","every valid sorted word: complete semantic digest")
gate(len(certificate["conditional_maxima"])==171,"complete forced-pair table domain")
gate([[a,b] for a,b,E,w in certificate["conditional_maxima"]]==[list(k) for k in conditional],"all distinct and repeated pair keys retained")
for a,b,E,owner in certificate["conditional_maxima"]:
    own=conditional[a,b]
    if own is None:
        gate(E is None and owner is None,"absent pair is absent from every full word")
        continue
    gate(E==own[0],"independent complete maximum for each forced pair")
    gate(owner is not None and valid_word(owner),"producer attaining word passes independent full profile test")
    need=Counter((a,b))
    have=Counter(owner)
    gate(all(have[d]>=count for d,count in need.items()),"attaining word contains both distinguished occurrences")
    gate(sum(native_cost[d] for d in owner)-T==E,"literal ceiling cost attains conditional maximum")
    gate(E<=pair_cost(a,b) and E<=188,"conditional maximum refines both relaxations")
count_sep=count_loc=0
lowest=None
decisions=[]
located_by_row={tuple(x["candidate"]):x for x in certificate["located_profiles"]}
gate(len(located_by_row)==153==len(certificate["located_profiles"]),"complete unique located-profile certificate")
for row in rows:
    e,p,q,dp,dq,Csep,Epair=row
    bound=conditional[tuple(sorted((dp,dq)))]
    if bound is None or Csep>bound[0]:
        count_sep+=1
        gate(bound is not None,"all actual old candidates have a compatible projected word")
        decisions.append(row+[bound[0],"separate",Csep,Csep-bound[0]])
        continue
    lo,hi,walls=event_sweep(T//e,p,q)
    surplus=e*lo-bound[0]
    gate(surplus>0,"each remaining candidate closes by location")
    item_cert=located_by_row[tuple(row)]
    gate((lo,hi,walls)==(item_cert["minimum"],item_cert["maximum"],item_cert["walls"]),"full native event profile agrees")
    gate((bound[0],e*lo,surplus)==(item_cert["conditional"],item_cert["credit"],item_cert["surplus"]),"same-pair maximum and location credit")
    gate(literal_grid(T//e,F(item_cert["min_phase"]),p,q)==lo,"literal strict minimum owner")
    gate(literal_grid(T//e,F(item_cert["max_phase"]),p,q)==hi,"literal strict maximum owner")
    decisions.append(row+[bound[0],"located",e*lo,surplus])
    item=(surplus,row,bound[0],lo,hi,walls)
    if lowest is None or item<lowest:
        lowest=item
    count_loc+=1
gate((count_sep,count_loc)==(835,153),"complete decision split")
gate(decisions==certificate["decisions"],"all988 complete decisions compared entry by entry")
gate(lowest==(18,[6,132,221,72,6,174,226],180,33,68,86),"exact worst located surplus and controls")

# Strong independent native control: use every projected wall of either danger
# predicate, including walls inactive in their intersection, at the worst row.
e,p,q=6,132,221
n=T//e
native=set()
for v in (p,q):
    for k in range(v):
        for sign in (-1,1):
            native.add(frac(F(n*(14*k+sign),14*v)))
native=sorted(native)
wall_counts=[]
cell_counts=[]
for i,w in enumerate(native):
    nxt=native[i+1] if i+1<len(native) else native[0]+1
    wall_counts.append(literal_grid(n,w,p,q))
    cell_counts.append(literal_grid(n,(w+nxt)/2,p,q))
gate(min(wall_counts+cell_counts)==33 and max(wall_counts+cell_counts)==68,"all raw native walls and cells: worst row")
gate(e*min(wall_counts+cell_counts)-conditional[6,72][0]==18,"strict positive worst located surplus")

# Existing hostile boundaries are retained in the complete decision table.
lookup={tuple(r[:3]):r for r in decisions}
gate(lookup[4,3,308][5:]==[172,219,188,"located",336,148],"location-only inherited boundary")
gate(lookup[4,23,323][5:]==[180,207,134,"separate",180,46],"owner-only inherited boundary")

# Physical multiplicity follows from a cyclic permutation, tested literally at
# both the strict boundary alpha=0 and non-wall rational translates. Let d=e*k
# with gcd(k,n)=1. Each normalized residue is taken exactly e times.
for e,p,q in ((4,3,308),(4,23,323),(6,132,221)):
    n=T//e
    for k in (1,n+1,2*n-1):
        gate(gcd(k,n)==1,"unit lift control")
        for alpha in (F(0),F(1,3),F(17,101)):
            den=n*alpha.denominator
            actual=0
            for j in range(T):
                num=alpha.numerator+k*j*alpha.denominator
                actual+=all(14*min(v*num%den,(-v*num)%den)<den for v in (p,q))
            gate(actual==e*literal_grid(n,alpha,p,q),"all t-grid lifts equal e times the n-grid count")

# No other clock is recomputed. Remove exactly t=16704 from the inherited full
# array, whose independent certificate and current truth surface are retained.
old_sets=[json.loads(line[10:]) for line in locate("third_20260906_grid_refined.out").read_text().splitlines() if line.startswith("SCALE_SET ")]
old=next(x["survivors"] for x in old_sets if x["name"]=="full_words_coupled")
gate(len(old)==8202 and old==sorted(set(old)),"complete sorted inherited clock set")
gate(digest(old)=="4f29481d984ead40d0144556ce1c45dce210e30b964bb65835a7904ca6353e59"==certificate["old_scale_sha256"],"inherited clock-set pin")
gate([t for t in old if t>14904]==[T],"only16704 lies above14904")
new=[t for t in old if t!=T]
gate(len(new)==8201 and max(new)==14904,"precise inherited-set deletion and new maximum")
gate(new==certificate["new_scales"] and digest(new)==certificate["new_scale_sha256"]=="ddbe0e091d36d54c8f6a7c8ea631bbf363d799b9adaa2ec4fe4cf56250d11a76","entire resulting clock array agrees")
print("FULL",visits,len(full_words),"pairs",sum(x is not None for x in conditional.values()),"wordhash",digest(full_words))
print("DECISIONS",count_sep,count_loc,"WORST",lowest)
print("Every one of171 forced-pair maxima/owners and988 decisions independently matches the frozen certificate.")
print("All153 located profiles: independent raw-interval event sweep plus literal extremum owners.")
print("Worst-row native partition",len(native),"walls and",len(native),"cells; min33,max68; strict endpoint conventions retained.")
print("Full t-grid multiplicity:27 literal phase/unit controls; both previous coordinate hostiles retained.")
print("Inherited8202-clock set minus16704:8201 clocks,maximum14904,newSHA",digest(new))
print("Scoped physical conclusion: clock16704 is excluded in the inherited actual decoder domain; LRC14 remains open.")
print("PASS:",GATES,"always-active exact gates.")

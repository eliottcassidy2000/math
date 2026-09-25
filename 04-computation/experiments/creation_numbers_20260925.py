"""Exact fruit input audit and finite signed Collatz certificates.

Standard library only. --verify reads the saved output and replays ordinary
steps without invoking the trajectory generator or elliptic-curve operations.
All checks remain enabled under python -O.
"""
import argparse
import importlib.util
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "05-knowledge/results/creation_numbers_20260925.out"
RAW = [
    "15,447,680,210,874,616,644,195,131,501,991,983,748,566,432,566,956,543,170,002,663,489,825,320,203,527,7999",
    "4,373,612,677,928,697,257,861,252,602,371,390,152,816,537,558,161,613,618,621,437,993,378,423,467,772,036",
    "368,751,317,941,299,998,271,978,115,652,254,748,254,929,799,689,719,709,962,831,374,716,372,246,340,555,799",
]


def check(ok, label="check failed"):
    if not ok:
        raise RuntimeError(label)


def v2(n):
    check(n > 0)
    return (n & -n).bit_length() - 1


def fruit(t):
    a, b, c = t
    return Fraction(a, b+c) + Fraction(b, a+c) + Fraction(c, a+b)


def cubic(t):
    a, b, c = t
    s = a+b+c
    return s**3 - 6*s*(a*b+a*c+b*c) + 7*a*b*c


def decimal(q):
    with localcontext() as ctx:
        ctx.prec = 40
        return str(Decimal(q.numerator)/Decimal(q.denominator))


def trace(n, sigma, label, limit=100000):
    initial = v2(n)
    x = n >> initial
    seen, states, exponents = {}, [], []
    ordinary_clock = initial
    first4 = None
    first_below = None
    peak = n
    while x not in seen and len(states) < limit:
        seen[x] = len(states)
        states.append(x)
        z = 3*x+sigma
        peak = max(peak, z)
        k = v2(z)
        # Ordinary time: one numerator step, then k halvings.
        if sigma == 1 and first4 is None and k >= 2 and z == 2**k:
            first4 = ordinary_clock + 1 + k-2
        ordinary_clock += 1+k
        x = z >> k
        exponents.append(k)
        if first_below is None and x < n:
            first_below = [len(exponents), ordinary_clock, x]
    check(x in seen, f"bound exhausted for {label}; no certificate produced")
    mu = seen[x]
    cycle = states[mu:]
    result = {
        "label": label, "start": n, "sigma": sigma,
        "initial_halvings": initial, "exponents": exponents,
        "preperiod": mu, "cycle": cycle, "basin_minimum": min(cycle),
        "ordinary_steps_to_cycle": initial+sum(1+k for k in exponents[:mu]),
        "ordinary_steps_to_4": first4,
        "first_below_start": first_below, "peak_ordinary": peak,
    }
    return result, states


def verify_certificate(c):
    """Independent replay: ordinary integer division, explicit every-step parity."""
    n, sigma = c["start"], c["sigma"]
    check(n > 0 and sigma in (-1, 1))
    x = n
    for _ in range(c["initial_halvings"]):
        check(x % 2 == 0)
        x //= 2
    check(x % 2 == 1)
    states = []
    seen = set()
    clock = c["initial_halvings"]
    first4 = 0 if n == 4 else None
    peak = n
    for k in c["exponents"]:
        check(isinstance(k, int) and k >= 1 and x > 0 and x % 2 == 1)
        check(x not in seen, "premature repeat")
        states.append(x)
        seen.add(x)
        x = 3*x+sigma
        clock += 1
        peak = max(peak, x)
        if x == 4 and first4 is None:
            first4 = clock
        for _ in range(k):
            check(x % 2 == 0, "false valuation certificate")
            x //= 2
            clock += 1
            if x == 4 and first4 is None:
                first4 = clock
        check(x % 2 == 1, "non-exact valuation certificate")
    mu = c["preperiod"]
    check(0 <= mu < len(states) and states[mu] == x)
    check(states[mu:] == c["cycle"])
    check(min(c["cycle"]) == c["basin_minimum"])
    check(c["ordinary_steps_to_cycle"] == c["initial_halvings"]+sum(1+k for k in c["exponents"][:mu]))
    check(c["peak_ordinary"] == peak)
    if sigma == 1:
        check(c["cycle"] == [1] and first4 == c["ordinary_steps_to_4"])
    return states


def verify_all(data):
    check(data["raw_input_strings"] == RAW)
    a,c,d=[int(s.replace(",","")) for s in RAW]
    check((d-9)%10 == 0)
    b=(d-9)//10
    check(data["literal_acd"] == [a,c,d])
    check(data["repaired_abc"] == [a,b,c])
    expected={(label,sigma):n for label,n in [("a",a),("b",b),("c",c),("literal_d",d)] for sigma in (-1,1)}
    for record in data["hostile_fruit_triples"]:
        for j,n in enumerate(record["triple"]):
            expected[f"fruit_{record['m']}G_{record['k']}T_coord{j}",-1]=n
    recovered={}
    for cert in data["certificates"]:
        key=cert["label"],cert["sigma"]
        check(key not in recovered and expected[key] == cert["start"])
        recovered[cert["label"],cert["sigma"]]=verify_certificate(cert)
    check(set(recovered)==set(expected) and len(recovered)==14)
    for join in data["common_future_joins"]:
        left=recovered[join["left"],join["sigma"]]
        right=recovered[join["right"],join["sigma"]]
        positions={x:i for i,x in enumerate(right)}
        first=next(([i,positions[x],x] for i,x in enumerate(left) if x in positions),None)
        stored=join["first_left_common_odd_state"]
        check(first == (list(stored) if stored is not None else None))
    for rec in data["hostile_fruit_triples"]:
        t = rec["triple"]
        check(all(x > 0 for x in t) and gcd(gcd(*t[:2]), t[2]) == 1)
        check(cubic(t) == 0 and fruit(t) == 4)
    check(cubic(data["repaired_abc"]) == 0)
    check(cubic(data["literal_acd"]) != 0)
    check(data["literal_F4"] == cubic(data["literal_acd"]))
    check(data["literal_fruit_sum"] == str(fruit(data["literal_acd"])))


def build():
    a, c, d = [int(s.replace(",", "")) for s in RAW]
    check((d-9) % 10 == 0)
    b = (d-9)//10
    check(a == 154476802108746166441951315019919837485664325669565431700026634898253202035277999)
    check(b == 36875131794129999827197811565225474825492979968971970996283137471637224634055579)
    check(c == 4373612677928697257861252602371390152816537558161613618621437993378423467772036)
    check(cubic((a,b,c)) == 0 and fruit((a,b,c)) == 4)
    check(cubic((a,c,d)) != 0)
    for t in [(a,b,c), (a,c,d), (1,1,1), (8,2,1)]:
        x,y,z=t
        check((fruit(t)-4)*(x+y)*(y+z)*(z+x) == cubic(t))
    check(fruit((1,1,1)) == Fraction(3,2))
    r, t = Fraction(a,b+c), Fraction(b*c,(b+c)**2)
    check(t == (r**3-3*r*r-3*r+1)/(6-r))
    check(r > 2 and (r-2)**2 > 3)
    check(4*r > 7 and (4*r-7)**2 < 65)
    certs, state_map, features = [], {}, {}
    for label,n in [("a",a),("b",b),("c",c),("literal_d",d)]:
        odd = n >> v2(n)
        features[label] = {"digits":len(str(n)), "v2":v2(n),
            "odd_part":odd, "odd_mod256":odd%256,
            "v2_odd_plus1":v2(odd+1), "v2_odd_minus1":v2(odd-1),
            "minus_germ_rank":v2(3*odd-1), "plus_germ_rank":v2(3*odd+1)}
        for sigma in (-1,1):
            cert, states = trace(n,sigma,label)
            run=v2(odd+sigma)-1
            check(cert["exponents"][:run]==[1]*run)
            check(cert["exponents"][run]>=2)
            features[label]["initial_one_halving_run_"+str(sigma)]=run
            certs.append(cert)
            state_map[label,sigma] = states
    joins = []
    for sigma in (-1,1):
        for left,right in combinations(features,2):
            ls,rs=state_map[left,sigma],state_map[right,sigma]
            positions={x:i for i,x in enumerate(rs)}
            pair=next(((i,positions[x],x) for i,x in enumerate(ls) if x in positions),None)
            joins.append({"sigma":sigma,"left":left,"right":right,
                "first_left_common_odd_state":pair})
    # Reuse the inherited exact group law; its main/GP work never runs on import.
    path=ROOT/"04-computation/experiments/collatz_mod6_20260921_fruit_rank_positive_multiples.py"
    spec=importlib.util.spec_from_file_location("inherited_fruit",path)
    inherited=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inherited)
    G=(Fraction(-4),Fraction(28)); T=(Fraction(56),Fraction(728))
    hostile=[]
    for m,k in [(13,1),(17,0)]:
        p=inherited.add(inherited.mul(m,G),inherited.mul(k,T))
        triple=inherited.fruit_triple(p)
        check(inherited.on_curve(p) and all(x>0 for x in triple))
        check(cubic(triple)==0 and fruit(triple)==4)
        basins=[]
        for j,n in enumerate(triple):
            cert,_=trace(n,-1,f"fruit_{m}G_{k}T_coord{j}")
            certs.append(cert);basins.append(cert["basin_minimum"])
        check(len(set(basins))<3)
        hostile.append({"m":m,"k":k,"triple":triple,
            "digits":[len(str(n)) for n in triple],"minus_basins":basins})
    data={"status":"FINITE-EXACT certificates; elementary proofs in note",
        "raw_input_strings":RAW,"literal_acd":[a,c,d],"repaired_abc":[a,b,c],
        "repair":"literal d=10b+9; b=(d-9)/10; a,c unchanged",
        "literal_F4":cubic((a,c,d)),"literal_fruit_sum":str(fruit((a,c,d))),
        "literal_fruit_decimal":decimal(fruit((a,c,d))),
        "ratios":{"b/c":decimal(Fraction(b,c)),"a/b":decimal(Fraction(a,b)),
            "a/c":decimal(Fraction(a,c)),"d/c":decimal(Fraction(d,c)),
            "a/(b+c)":decimal(r),"bc/(b+c)^2":decimal(t)},
        "features":features,"common_future_joins":joins,
        "hostile_fruit_triples":hostile,"certificates":certs,
        "universe":"Four literal/repaired distinct integers on both sheets; two inherited positive fruit triples on minus; 100000 odd-step bound each. All completed."}
    verify_all(data)
    # Hostile to a vacuous or disabled verifier: alter the very first valuation.
    bad=dict(certs[0]);bad["exponents"]=list(bad["exponents"])
    bad["exponents"][0]+=1
    caught=False
    try: verify_certificate(bad)
    except RuntimeError: caught=True
    check(caught,"corrupted certificate was accepted")
    OUT.write_text(json.dumps(data,separators=(",",":"))+"\n",encoding="utf-8")
    print("PASS: literal and repaired arithmetic, ratios, 14 trajectory certificates, and corruption control")
    print("ratios",data["ratios"])
    print("literal fruit sum",data["literal_fruit_decimal"])
    for cc in certs:
        print(cc["label"],cc["sigma"],"initial/odd/ordinary",cc["initial_halvings"],cc["preperiod"],cc["ordinary_steps_to_cycle"],"basin",cc["basin_minimum"],"ordinary_to4",cc["ordinary_steps_to_4"])
    print("joins",joins)
    print("output",OUT)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    if args.verify:
        data=json.loads(OUT.read_text(encoding="utf-8"))
        verify_all(data)
        print("PASS: standalone ordinary-step verifier",len(data["certificates"]),"certificates")
    else:
        build()

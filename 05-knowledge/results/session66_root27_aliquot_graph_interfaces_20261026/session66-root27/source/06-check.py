"""Exact tree/carry, dual-section and native-guard controls; standard library."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def val(n, p=2):
    require(n > 0, "positive valuation")
    r = 0
    while n % p == 0:
        n //= p
        r += 1
    return r


def carry(word):
    total = 0
    for i, bit in enumerate(word):
        total = 3*total + int(bit)*(1 << i)
    return total


def qvalue(z, word):
    return 3**len(word)*z + carry(word)


def residue(word):
    modulus = 1 << len(word)
    return (-carry(word)*pow(3**len(word), -1, modulus)) % modulus


def bits(z, n):
    out = []
    for _ in range(n):
        bit = z % 2
        out.append(str(bit))
        z = (3*z+bit)//2
    return "".join(out), z


def forest_ok(word, q):
    for i, bit in enumerate(word):
        q += 3*int(bit)-1
        if i < len(word)-1 and q <= 0:
            return False
    return q == 0


@lru_cache(None)
def trees(e):
    if e == 0:
        return ("0",)
    out = []
    for a in range(e):
        for b in range(e-a):
            c = e-1-a-b
            out.extend("1"+u+v+w
                       for u in trees(a) for v in trees(b) for w in trees(c))
    return tuple(out)


@lru_cache(None)
def forests(q, e):
    if q == 0:
        return ("",) if e == 0 else ()
    return tuple(u+v for a in range(e+1)
                 for u in trees(a) for v in forests(q-1, e-a))


def native_until(z, a, limit):
    word = []
    executed = 0
    for _ in range(limit+1):
        c = z % 2
        word.append(str(c))
        if a == 0 and c == 0:
            return {"halted": True, "z": z, "a": a,
                    "executed": executed, "word": "".join(word)}
        if executed == limit:
            break
        require(a >= 0, "native negative credit")
        a += 3*c-1
        z = (3*z+c)//2
        executed += 1
    return {"halted": False, "z": z, "a": a,
            "executed": executed, "word": "".join(word)}


def v_return(z):
    ell = val(3*z+1)
    return 3**(ell-1)*((3*z+1) >> ell), ell


def v_inverse(y):
    r = val(y, 3)
    m = y//3**r
    numerator = (1 << (r+1))*m-1
    if numerator % 3:
        return None
    z = numerator//3
    require(0 < z < y and z % 2 == 1, "dual inverse boundary")
    require(v_return(z) == (y, r+1), "dual inverse replay")
    return z, r+1


def main():
    stats = Counter()
    tree_counts, forest_counts = [], []
    universes = []
    for e in range(6):
        words = trees(e)
        expected = comb(3*e, e)//(2*e+1)
        require(len(words) == len(set(words)) == expected, "complete tree count")
        tree_counts.append(expected)
        universes.extend((1, w) for w in words)
    for e in range(4):
        words = forests(6, e)
        expected = 6*comb(3*e+6, e)//(3*e+6)
        require(len(words) == len(set(words)) == expected, "complete forest count")
        forest_counts.append(expected)
        universes.extend((6, w) for w in words)

    # Each full language word gets its own literal positive-integer lift.
    for q, word in universes:
        require(forest_ok(word, q), "forest prefix guards")
        r = residue(word)
        z = r + (1 << len(word))
        actual, end = bits(z, len(word))
        require(actual == word and qvalue(z, word) == (1 << len(word))*end,
                "exact cylinder and affine endpoint")
        run = native_until(z, q-1, len(word))
        require(run["halted"] and run["word"] == word
                and run["executed"] == len(word)-1,
                "last symbol is blocked, not executed")
        seed = 7 if q == 1 else 39
        require(qvalue(seed, word) % (1 << len(word)) != 0,
                "bounded marked-source tree census")
        for cut in range(len(word)+1):
            u, v = word[:cut], word[cut:]
            require(carry(word) == 3**len(v)*carry(u)+(1 << len(u))*carry(v),
                    "ordered affine concatenation")
            stats["concatenation_cuts"] += 1
        rotation = word[1:]+word[0]
        require(2*carry(rotation) == 3*carry(word)
                -int(word[0])*(3**len(word)-(1 << len(word))),
                "rotation carries source information")
        stats["literal_forest_lifts"] += 1

    # Every binary word through length8: finite-cylinder converse independent
    # of tree syntax, including incompatible source residues.
    for n in range(1, 9):
        for tup in product("01", repeat=n):
            word = "".join(tup)
            r = residue(word)
            require(bits(r, n)[0] == word, "binary cylinder positive control")
            require(bits(r ^ (1 << (n-1)), n)[0] != word,
                    "last-bit cylinder hostile")
            stats["binary_cylinders"] += 1

    # Leaf replacement controls on all trees through4 internal nodes and
    # all six-tree forests through2 internal nodes, at four literal sources.
    leaf_words = [w for e in range(5) for w in trees(e)]
    leaf_words += [w for e in range(3) for w in forests(6, e)]
    for word in leaf_words:
        for p, bit in enumerate(word):
            if bit != "0":
                continue
            suffix = word[p+1:]
            new = word[:p]+"1000"+suffix
            odd = 3**(len(suffix)+3)-38*carry(suffix)
            require(odd % 2 == 1, "leaf correction is odd")
            for z in (3, 7, 15, 39):
                before, after = qvalue(z, word), qvalue(z, new)
                require(after == 27*before+(1 << p)*odd, "leaf carry identity")
                h, new_h = val(before), val(after)
                if p != h:
                    require(new_h == min(p, h), "only active leaf improves precision")
                else:
                    require(new_h >= h+1, "active leaf improves at least one bit")
                stats["leaf_expansion_controls"] += 1

    # Dual cross-section and its guard: all odd sources1..1999, q1..8.
    image_bank = {}
    for z in range(1, 2000, 2):
        y, ell = v_return(z)
        require(y > z and y not in image_bank, "dual growth and injectivity")
        image_bank[y] = (z, ell)
        literal, target = bits(z, ell)
        require(literal == "1"+"0"*(ell-1) and target == y,
                "dual literal return")
        for q in range(1, 9):
            run = native_until(z, q-1, ell)
            if ell <= q+2:
                require(not run["halted"] and run["executed"] == ell
                        and (run["z"], run["a"]) == (y, q+2-ell),
                        "dual full return guard")
            else:
                require(run["halted"] and run["executed"] == q+2,
                        "dual blocked before formal return")
            stats["dual_guard_controls"] += 1
    for y in range(1, 2000, 2):
        require(v_inverse(y) == image_bank.get(y), "complete dual inverse universe")
    require(v_inverse(7) is None, "dual root7")
    require(v_return(5)[0] == 27 and v_return(7)[0] == 11,
            "pointwise growth is not order preservation")

    # Exact warmup of both sections, with the initial credit0 allowed at odd z.
    z, a = 7, 0
    warmup = [(z, a)]
    for c in (1, 1, 1, 0):
        require(z % 2 == c, "blank prefix")
        a += 3*c-1
        z = (3*z+c)//2
        warmup.append((z, a))
    require(warmup == [(7,0),(11,2),(17,4),(26,6),(39,5)], "section matching")

    # A bounded frontier is checked against V, not used as a larger orbit census.
    word, z, q, h, total_ell = "0", 7, 1, 0, 0
    for step in range(32):
        require(word[h:] == "0"*q and len(word) == h+q, "live frontier suffix")
        require(qvalue(7, word) == (1 << h)*3**q*z, "live frontier carry")
        require(val(qvalue(7, word)) == h and word[h] == "0", "active leaf")
        word = word[:h]+"1000"+word[h+1:]
        target, ell = v_return(z)
        new_h = val(qvalue(7, word))
        require(new_h == h+ell, "frontier precision is dual return length")
        q += 3-ell
        require(q >= 1, "bounded frontier remains open")
        z, h = target, new_h
        total_ell += ell
        require(q == 1+3*(step+1)-total_ell, "zero-offset dual ballot ledger")
        stats["frontier_return_controls"] += 1

    # Pumped families: all repeats have a proof; bounded repeats cross-check it.
    fuel_histogram = Counter()
    pump_count = 0
    for e in range(1, 6):
        for tree in trees(e):
            block = tree[:-1]
            L, C = len(block), carry(block)
            A, B = 3**L, 1 << L
            d7 = val((A-B)*7+C)
            fuel_histogram[d7] += 1
            require(d7 % L != 0, "entire marked-source pump family excluded")
            for source in (2, 5, 7, 15, 39):
                fuel = val((A-B)*source+C)
                for repeats in range(5):
                    candidate = block*repeats+"0"
                    require(forest_ok(candidate, 1), "pump is a tree")
                    realized = qvalue(source, candidate) % (1 << len(candidate)) == 0
                    require(realized == (fuel == L*repeats), "exact pumped-family iff")
                    stats["pump_cylinder_controls"] += 1
            pump_count += 1
    require(pump_count == 344 and dict(fuel_histogram) ==
            {1:183, 2:111, 3:10, 4:31, 5:8, 6:1}, "finite pump bank")
    for repeats in range(16):
        word = "100"*repeats+"0"
        require(carry(word) == 27*(27**repeats-8**repeats)//19,
                "denominator19 rightcomb carry")
    require(val(19*7+9) == val(19*15+9) == 1, "both sources exclude every rightcomb")
    require([residue(w) for w in ("100000000", "010000000")] == [341,398],
            "same forest counts different arithmetic sources")

    hostile = native_until(15, 0, 30)
    require(hostile == {"halted": True, "z": 78894, "a": 0, "executed": 21,
                       "word": "1111000010000010100000"}, "noncomb halt hostile")
    require(forest_ok(hostile["word"], 1) and hostile["word"].count("1") == 7
            and residue(hostile["word"]) == 15, "noncomb arithmetic certificate")

    # Candidate charge is refuted on the actual root27 ray, within its old prefix.
    x, a = 27, 6
    for index in range(20):
        k = val(3*x-1)
        d, e = val(x-1), val(x-1, 3)
        y = 3**(k-1)*((3*x-1) >> k)
        new_a = a+2*k-3
        if k == 1:
            require((val(y-1), val(y-1, 3)) == (d-1,e+1)
                    and new_a+val(y-1,3) == a+e, "two-prime run invariant")
        else:
            require(val(y-1,3) == 0 and new_a-(a+e) == 2*k-3-e,
                    "reset charge loss")
        if index == 19:
            charge_hostile = {"index": index, "x": x, "a": a, "k": k,
                              "e": e, "before": a+e, "after": new_a}
            require(charge_hostile == {"index":19,"x":57971215,"a":21,"k":2,
                                       "e":4,"before":25,"after":22},
                    "actual marked-ray charge hostile")
        x, a = y, new_a

    return {"status":"PASS", "scope":"ELEMENTARY_PROOFS_AND_BOUNDED_EXACT_CONTROLS",
            "tree_counts_e0_to5":tree_counts, "six_forest_counts_e0_to3":forest_counts,
            "controls":dict(sorted(stats.items())), "pump_family_count":pump_count,
            "pump_fuel_histogram":dict(sorted(fuel_histogram.items())),
            "blank_warmup_z_a":warmup, "noncomb_halt":hostile,
            "actual_root27_charge_failure":charge_hostile,
            "open":"Universal root7 tree-source exclusion / root27 debt inequality"}


if __name__ == "__main__":
    print(json.dumps(main(), sort_keys=True, indent=2))

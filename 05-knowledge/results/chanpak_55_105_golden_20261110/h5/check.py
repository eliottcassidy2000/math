#!/usr/bin/env python3
"""Fixed-core two-flank attachment compiler; bounded exact inverse controls."""
import hashlib
import itertools
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    q, bit = divmod(n, 2)
    return g(q) if bit else g(q) + q + (q & 1)


def word(n, width=None):
    return [(n >> i) & 1 for i in reversed(range(n.bit_length() if width is None else width))]


def induced(rows, keep):
    return tuple(sum(((rows[i] >> j) & 1) << k for k, j in enumerate(keep)) for i in keep)


@lru_cache(None)
def hp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, count in enumerate(dp[mask]):
            if not count:
                continue
            todo = rows[i] & ~mask
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def literal_hp(rows):
    return sum(all((rows[p[i]] >> p[i+1]) & 1 for i in range(len(p)-1))
               for p in itertools.permutations(range(len(rows))))


def tournament(rows):
    return all(not (row >> i) & 1 and row < 1 << len(rows) for i, row in enumerate(rows)) and all(
        ((rows[i] >> j) & 1)+((rows[j] >> i) & 1) == 1
        for i in range(len(rows)) for j in range(i+1,len(rows)))


def attach(rows, bits):
    n = len(rows)
    return tuple(row | ((1-bits[i]) << n) for i, row in enumerate(rows)) + (
        sum(bit << i for i, bit in enumerate(bits)),)


def graph(core, a, b, outmask, widths=None):
    left = word(a, None if widths is None else widths[0])
    right = [1-bit for bit in reversed(word(b, None if widths is None else widths[1]))]
    n, shift = len(left)+len(core)+len(right), len(left)
    rows = [sum(1 << j for j in range(i+1, n)) for i in range(n)]
    for i, row in enumerate(core):
        rows[shift+i] = (row << shift) | sum(1 << j for j in range(shift+len(core), n))
    bits = left + [(outmask >> i) & 1 for i in range(len(core))] + right
    return attach(tuple(rows), bits)


@lru_cache(None)
def coefficients(core):
    n, h = len(core), hp(core)
    start, end = [0]*n, [0]*n
    sp, ep, sm, em = [0]*n, [0]*n, [0]*n, [0]*n
    wplus = wminus = 0
    for length in range(1, n+1):
        for p in itertools.permutations(range(n), length):
            if not all((core[p[i]] >> p[i+1]) & 1 for i in range(length-1)):
                continue
            rest = [i for i in range(n) if i not in p]
            weight = hp(induced(core, rest))
            sign = (-1)**length
            sp[p[0]] += weight
            ep[p[-1]] += weight
            sm[p[0]] += sign*weight
            em[p[-1]] += sign*weight
            wplus += weight
            wminus += sign*weight
            if length == n:
                start[p[0]] += 1
                end[p[-1]] += 1
    require(wminus == -h and sm == [-v for v in end] and em == [-v for v in start],
            "signed path reciprocity")
    q = h+wplus
    table = []
    for outmask in range(1 << n):
        o = [i for i in range(n) if (outmask >> i) & 1]
        inside = [i for i in range(n) if i not in o]
        t = hp(attach(core, [(outmask >> i) & 1 for i in range(n)]))
        table.append({"outmask": outmask, "e": sum(ep[i] for i in inside),
                      "s": sum(sp[i] for i in o), "left_parity": sum(start[i] for i in inside),
                      "right_parity": sum(end[i] for i in o), "t": t})
    return {"h": h, "q": q, "start": start, "end": end,
            "cover_start": sp, "cover_end": ep, "cuts": table}


def value(cf, cut, a, b):
    return (cf["q"]*a*b + cut["e"]*a + cut["s"]*b
            + cf["h"]*(g(a)+g(b)-2) + cut["t"]
            + cut["left_parity"]*(a & 1) + cut["right_parity"]*(b & 1))


def support(core, limit):
    """Complete all-order numerical support <= limit, with one reconstructible witness."""
    cf = coefficients(core)
    answer = {}
    candidates = 0
    def record(cut, a, b):
        nonlocal candidates
        candidates += 1
        n = value(cf, cut, a, b)
        if n <= limit:
            candidate = (cut["outmask"], a, b)
            if n not in answer or (a.bit_length()+b.bit_length(), candidate) < (
                    answer[n][1].bit_length()+answer[n][2].bit_length(), answer[n]):
                answer[n] = candidate
    for cut in cf["cuts"]:
        e, s, q, t = cut["e"], cut["s"], cf["q"], cut["t"]
        if e and s:
            for a in range(max(0, (limit-t)//e+1)):
                for b in range(max(0, (limit-t-e*a)//(q*a+s)+1)):
                    record(cut, a, b)
        elif e:  # O empty; A=0 reduces to h*g(B).
            record(cut, 0, 0)
            for r in range(1, limit//cf["h"]):
                record(cut, 0, 2*r)
            for a in range(1, max(0, (limit-t)//e+1)):
                for b in range(max(0, (limit-t-e*a)//(q*a)+1)):
                    record(cut, a, b)
        else:  # O full; converse of preceding inverse.
            record(cut, 0, 0)
            for r in range(1, limit//cf["h"]):
                record(cut, 2*r, 0)
            for b in range(1, max(0, (limit-t)//s+1)):
                for a in range(max(0, (limit-t-s*b)//(q*b)+1)):
                    record(cut, a, b)
    return answer, candidates


def canonical_rows(rows):
    return min(induced(rows, p) for p in itertools.permutations(range(len(rows))))


def main():
    cores = {"H1": (0,), "H3": (2,4,1), "H5": (4,1,10,3)}
    limit = 1001
    controls = Counter()
    for name, core in cores.items():
        require(tournament(core), "valid core tournament")
        cf = coefficients(core)
        if len(core) > 1:
            # Recovered Session197 ONE_CORE_BINARY_INVERSE equation(8):
            # stationary attachment beats every core vertex except named port j.
            for j in range(len(core)):
                cut = cf["cuts"][((1 << len(core))-1) ^ (1 << j)]
                h, q = cf["h"], cf["q"]
                start, end = cf["start"][j], cf["end"][j]
                cover_start, cover_end = cf["cover_start"][j], cf["cover_end"][j]
                require((cut["e"],cut["s"],cut["left_parity"],cut["right_parity"],cut["t"])
                        == (cover_end,q-h-cover_start,start,h-end,h-start+cover_end),
                        "inherited stationary named-port coefficient map")
                controls["stationary_named_core_ports"] += 1
                for a in range(8):
                    for b in range(8):
                        inherited = (h*(g(a)+g(b)-1)+(cover_end+q*b)*(a+1)
                                     -start*((a+1)&1)-b*(cover_start+h)+(b&1)*(h-end))
                        require(value(cf,cut,a,b) == inherited,
                                "inherited stationary binary formula restriction")
                        controls["stationary_parameter_checks"] += 1
        for total in range(4):
            for aa in range(total+1):
                bb = total-aa
                for a in range(1 << aa):
                    for b in range(1 << bb):
                        for cut in cf["cuts"]:
                            rows = graph(core, a, b, cut["outmask"], (aa,bb))
                            require(tournament(rows), "valid constructed tournament")
                            require(hp(rows) == value(cf,cut,a,b), "graph compiler")
                            if len(rows) <= 6:
                                require(literal_hp(rows) == hp(rows), "literal compiler")
                                controls["literal_permutation_checks"] += 1
                            controls["subset_dp_checks"] += 1
    # Independent tiny labelled order4 universe: exactly one strong H=5 class.
    classes = set()
    pairs = list(itertools.combinations(range(4),2))
    for bits in range(64):
        rows = [0]*4
        for j,(u,v) in enumerate(pairs):
            if (bits >> j) & 1:
                u,v=v,u
            rows[u] |= 1 << v
        if hp(tuple(rows)) == 5:
            classes.add(canonical_rows(tuple(rows)))
        controls["order4_labelled"] += 1
    require(len(classes) == 1 and canonical_rows(cores["H5"]) in classes, "unique H5 Q4")
    images, candidate_counts = {}, {}
    for name, core in cores.items():
        images[name], candidate_counts[name] = support(core,limit)
        # Verify every selected inverse witness by an independent subset DP.
        for target,(outmask,a,b) in images[name].items():
            rows = graph(core,a,b,outmask)
            require(hp(rows) == target, "inverse graph witness")
            controls["inverse_witnesses"] += 1
    low = set(images["H1"]) | set(images["H3"])
    combined = low | set(images["H5"])
    misses = [n for n in range(1,limit+1,2) if n not in combined]
    require([n for n in range(1,56,2) if n not in low] == [7,21,55], "prior55 control")
    require(misses == [7,21,173,335,663], "complete bounded bank inverse")
    require(105 not in images["H5"], "no H5 parent of105 at any order")
    samples = {}
    for n in (55,105):
        samples[n] = {}
        for name, core in cores.items():
            if n in images[name]:
                outmask,a,b = images[name][n]
                rows = graph(core,a,b,outmask)
                deletions = [hp(induced(rows,[j for j in range(len(rows)) if j != i]))
                             for i in range(len(rows))]
                samples[n][name] = {"outmask":outmask,"A":a,"B":b,
                                    "rows":rows,"deletions":deletions}
    # Uniform Q4 cut: direct check that treating it as a scalar module loses reentry.
    rows = graph(cores["H5"],1,1,0)
    quotient = graph((0,),1,1,0)
    require(hp(rows) != 5*hp(quotient), "uniform scalar collapse is false")
    core9 = (2,4,9,19,7)
    witness173 = graph(core9,1,1,2)
    require(tournament(core9) and tournament(witness173), "valid173 construction")
    require(hp(core9) == literal_hp(core9) == 9, "independent H9 core")
    require(hp(witness173) == literal_hp(witness173) == 173, "independent173 witness")
    require(value(coefficients(core9),coefficients(core9)["cuts"][2],1,1) == 173,
            "generic compiler H9 control")
    deletions173 = []
    for i in range(len(witness173)):
        sub = induced(witness173,[j for j in range(len(witness173)) if j != i])
        require(hp(sub) == literal_hp(sub), "independent173 deletion")
        deletions173.append(hp(sub))
    require(deletions173 == [27,71,51,27,51,71,77,9], "sharp173 deletion floor")
    payload = {"status":"FINITE-EXACT: bounded target inverse, not a size-cutoff search",
               "limit":limit,"coefficients":coefficients(cores["H5"]),
               "controls":dict(controls),"inverse_candidate_counts":candidate_counts,
               "image_sizes":{k:len(v) for k,v in images.items()},
               "misses_H1_H3_H5":misses,
               "new_values_from_H5":[n for n in sorted(images["H5"]) if n not in low],
               "samples":samples,"uniform_reentry":{"rows":rows,"H":hp(rows),
                   "quotient_rows":quotient,"quotient_H":hp(quotient)},
               "sharp173":{"core":core9,"outmask":2,"A":1,"B":1,
                   "rows":witness173,"H":173,"deletions":deletions173},
               "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(payload,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()

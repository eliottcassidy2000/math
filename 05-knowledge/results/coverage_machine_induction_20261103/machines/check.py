#!/usr/bin/env python3
"""One fixed FRACTRAN reader -> named two-port births -> actual strong sources.

No external inputs. Normal and -O runs must print byte-identical JSON.
"""
import json
from fractions import Fraction
from math import prod


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes(n):
    out = []
    x = 2
    while len(out) < n:
        if all(x % p for p in out if p * p <= x):
            out.append(x)
        x += 1
    return out


REG = dict(zip(("R", "P", "H", "Q", "S"), (2, 3, 5, 7, 11)))
STATES = ["main", "restore", "div", "divback"]
STATES += [f"{b}{phase}" for b in (0, 1)
           for phase in ("r", "rb", "p", "pb", "s", "sb")]
STATES += ["halt"]
CONTROL = dict(zip(STATES, primes(22)[5:]))


def program():
    # (source, required/consumed register units, target, produced units)
    rules = [
        ("main", {"R": 1}, "restore", {}),
        ("main", {}, "halt", {}),
        ("restore", {}, "div", {"R": 1}),
        ("div", {"R": 2}, "divback", {"Q": 1}),
        ("div", {"R": 1}, "1r", {}),
        ("div", {}, "0r", {}),
        ("divback", {}, "div", {}),
    ]
    for b in (0, 1):
        rules += [
            (f"{b}r", {"Q": 1}, f"{b}rb", {"R": 1}),
            (f"{b}rb", {}, f"{b}r", {}),
            (f"{b}r", {}, f"{b}p", {}),
            (f"{b}p", {"P": 1}, f"{b}pb", {"S": 4, "H": 3-b}),
            (f"{b}pb", {}, f"{b}p", {}),
            (f"{b}p", {}, f"{b}s", {}),
            (f"{b}s", {"S": 1}, f"{b}sb", {"P": 1}),
            (f"{b}sb", {}, f"{b}s", {}),
            (f"{b}s", {}, "main", {"P": 2*b, "H": 2*b}),
        ]
    return rules


RULES = program()
FRACTIONS = [Fraction(CONTROL[t] * prod(REG[k]**v for k, v in inc.items()),
                      CONTROL[s] * prod(REG[k]**v for k, v in dec.items()))
             for s, dec, t, inc in RULES]


def register_run(n):
    state = "main"
    r = dict(R=n, P=2, H=3, Q=0, S=0)
    steps, emitted, returns = 0, [], [(n, 2, 3)]
    while state != "halt":
        for s, dec, t, inc in RULES:
            if state == s and all(r[k] >= v for k, v in dec.items()):
                for k, v in dec.items():
                    r[k] -= v
                for k, v in inc.items():
                    r[k] += v
                if s == "div" and t in ("0r", "1r"):
                    emitted.append(int(t[0]))
                if t == "main":
                    returns.append((r["R"], r["P"], r["H"]))
                    need(r["Q"] == r["S"] == 0, "scratch at return")
                state = t
                steps += 1
                break
        else:
            raise RuntimeError("nonhalting state has no rule")
        need(steps < 2000000, "declared simulation cap")
    return r, steps, emitted, returns


def literal_run(n):
    x = CONTROL["main"] * 2**n * 3**2 * 5**3
    steps = 0
    while True:
        for f in FRACTIONS:
            if x % f.denominator == 0:
                x = x // f.denominator * f.numerator
                steps += 1
                break
        else:
            return x, steps
        need(steps < 2000000, "literal cap")


def ca_word(n):
    # Radius-one local rule (data, marker)'_i=(data_(i+1),marker_i).
    # Blank is -1; a finite binary input occupies cells 0,...,l-1.
    bits = [int(x) for x in bin(n)[:1:-1]] if n else []
    cells = {i: (b, int(i == 0)) for i, b in enumerate(bits)}
    cells.setdefault(0, (-1, 1))
    out = []
    for _ in range(len(bits) + 1):
        data, marker = cells[0]
        need(marker == 1, "stationary CA output port")
        if data == -1:
            break
        out.append(data)
        lo, hi = min(cells), max(cells)
        cells = {i: (cells.get(i+1, (-1, 0))[0],
                     cells.get(i, (-1, 0))[1])
                 for i in range(lo-1, hi+1)}
        cells = {i: v for i, v in cells.items() if v != (-1, 0)}
    need(out == bits, "radius-one CA reading")
    return out


def tm_word(bits):
    # A two-tape finite-state Turing transducer; input is read-only.
    # I0,I1 write header10. Read writes the current bit; Zero writes0.
    state, head, output, steps = "I0", 0, [], 0
    while state != "halt":
        if state == "I0":
            output.append(1)
            state = "I1"
        elif state == "I1":
            output.append(0)
            state = "Read"
        elif state == "Read":
            if head == len(bits):
                state = "halt"
            else:
                output.append(bits[head])
                state = "Zero"
        else:
            need(state == "Zero", "TM control")
            output.append(0)
            head += 1
            state = "Read"
        steps += 1
    need(steps == 2*len(bits)+3, "TM output clock")
    return output


def build(bits):
    # One apex 0; base orientation is increasing labels. Seed is word10.
    w = [1, 0]
    rows = [2, 4, 1]
    ports = (0, 2)
    births = []
    for b in bits:
        for c in (b, 0):
            n = len(rows)
            previous = ports
            for i in range(1, n):
                rows[i] |= 1 << n
            if c:
                rows[0] |= 1 << n
                rows.append(0)
            else:
                rows.append(1)
            ports = (0, n)
            births.append(dict(vertex=n, apex_beats_new=c, old_ports=list(previous),
                               new_ports=list(ports)))
            w.append(c)
    return rows, ports, w, births


def path_count(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            value = dp[m][v]
            todo = rows[v] & ~m
            while value and todo:
                bit = todo & -todo
                u = bit.bit_length() - 1
                dp[m | bit][u] += value
                todo ^= bit
    return sum(dp[-1])


def strong(rows):
    n = len(rows)
    for v in range(n):
        reach = 1 << v
        while True:
            new = reach
            for u in range(n):
                if reach >> u & 1:
                    new |= rows[u]
            if reach == new:
                break
            reach = new
        if reach != (1 << n)-1:
            return False
    return True


def decode(rows, apex=0):
    base = [v for v in range(len(rows)) if v != apex]
    # Transitive base has pairwise distinct outdegrees 0,...,m-1.
    mask = sum(1 << v for v in base)
    base.sort(key=lambda v: -(rows[v] & mask).bit_count())
    need(all(rows[a] >> b & 1 for i, a in enumerate(base) for b in base[i+1:]),
         "named apex deletion must be transitive")
    w = [rows[apex] >> v & 1 for v in base]
    need(w[:2] == [1, 0] and len(w) % 2 == 0, "code header")
    need(all(w[i] == 0 for i in range(3, len(w), 2)), "code separators")
    bits = w[2::2]
    need(not bits or bits[-1] == 1, "canonical positive binary input")
    return sum(b << i for i, b in enumerate(bits))


def chart(rows):
    n = len(rows)
    need(all(rows[i] >> (i+1) & 1 for i in range(n-1)), "fixed Hamilton path")
    pairs = [(i, j) for i in range(n) for j in range(i+2, n)]
    d = {(i, j): int(not (rows[i] >> j & 1)) for i, j in pairs}
    defect = {(i, j): d[i, j] ^ d[n-1-j, n-1-i] for i, j in pairs}
    return d, defect


def complement(rows):
    n = len(rows)
    return [sum(1 << j for j in range(n) if j != i and
                (bool(rows[i] >> j & 1) if abs(i-j) == 1
                 else not bool(rows[i] >> j & 1))) for i in range(n)]


def relabel(rows, perm):
    out = [0]*len(rows)
    for i, row in enumerate(rows):
        for j in range(len(rows)):
            if row >> j & 1:
                out[perm[i]] |= 1 << perm[j]
    return out


def main():
    need(len(FRACTIONS) == 25, "fixed instruction count")
    records = []
    for n in range(32):
        r, steps, bits, returns = register_run(n)
        need(ca_word(n) == bits, "FRACTRAN/CA emitted word")
        p, h = 2, 3
        q = n
        for j, b in enumerate(bits):
            q //= 2
            p, h = 4*p+2*b, h+(3-b)*p+2*b
            need(returns[j+1] == (q, p, h), "complete return clock")
        need(r == dict(R=0, P=p, H=h, Q=0, S=0), "final registers")
        rows, ports, w, births = build(bits)
        need(tm_word(bits) == w, "Turing transducer / graph word")
        need(strong(rows), "strong return")
        need(decode(rows) == n, "integer decoding")
        perm = list(reversed(range(len(rows))))
        need(decode(relabel(rows, perm), perm[0]) == n, "marked isomorphism decoding")
        d, defect = chart(rows)
        need(sum(defect.values()) == 2*(2*len(bits)-n.bit_count()), "exact reflection defect")
        need(not any(defect.values()) if n == 0 else any(defect.values()), "line colour")
        need(complement(complement(rows)) == rows, "native complement involution")
        if n:
            need(strong(complement(rows)), "strong native return")
        rec = dict(input=n, bits=bits, prefix=p, H=h, order=len(rows),
                   steps=steps, blue=not any(defect.values()),
                   reflection_defect_weight=sum(defect.values()),
                   rows=rows, ports=list(ports), births=births)
        if n < 16:
            value, literal_steps = literal_run(n)
            need(value == CONTROL["halt"]*3**p*5**h, "literal FRACTRAN output")
            need(literal_steps == steps, "literal step alignment")
            need(path_count(rows) == h, "independent subset path count")
            rec["literal_and_DP"] = True
            rec["native_partner_H"] = path_count(complement(rows))
        records.append(rec)
    # Hostile: noncanonical words 0 and 1 produce distinct marked graphs, same H9.
    z, _, _, _ = build([0])
    o, _, _, _ = build([1])
    need(z != o and path_count(z) == path_count(o) == 9, "H forgets input address")
    # Guard ordering matters: at R=2, the two-unit division rule must beat one-unit.
    div = [rule for rule in RULES if rule[0] == "div"]
    need(div[0][1] == {"R": 2} and div[1][1] == {"R": 1}, "priority is semantic")
    print(json.dumps(dict(
        status="FINITE-EXACT controls for proved all-input transducer",
        scope=dict(register_CA_TM_inputs="0..31", literal_and_DP_inputs="0..15",
                   peak_DP_order=11, no_spectrum_coverage_claim=True),
        register_primes=REG, control_primes=CONTROL,
        fractions=[f"{f.numerator}/{f.denominator}" for f in FRACTIONS],
        records=records, hostile_same_H9_words=["1000", "1010"]),
        sort_keys=True, indent=2))


if __name__ == "__main__":
    main()

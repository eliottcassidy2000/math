"""Exact -17 return cylinders and their extension of the inherited bank.

Run from any directory. All checks remain active under Python -O. The only
imported research code regenerates the explicitly bounded 171-row old bank;
its rows are independently compared with the retained certificate table.
"""
from fractions import Fraction
from pathlib import Path
import json
import re
import runpy

ROOT = Path(__file__).resolve().parents[2]
WORD = (1, 1, 1, 2, 1, 1, 4)
P, Q, CENTER = 2187, 2048, 17


def need(condition, message):
    if not condition:
        raise ValueError(message)


def v2(n):
    need(n > 0, "positive valuation argument")
    return (n & -n).bit_length() - 1


def step(n):
    need(n > 0 and n % 2 == 1, "positive odd source")
    z = 3*n + 1
    a = v2(z)
    return z >> a, a


def compose_word(word):
    """Generic composition, independent of the fixed-point closed form."""
    multiplier, denominator, carry = 1, 1, 0
    for a in word:
        multiplier, denominator, carry = (
            3*multiplier, denominator*2**a, 3*carry+denominator)
    return multiplier, denominator, carry


def follows(n, word):
    for a in word:
        n, actual = step(n)
        if actual != a:
            return False
    return True


def cylinder(m):
    need(m >= 1, "positive block count")
    t = 1
    while 2**t*(Q**m-CENTER) <= P**m-CENTER:
        t += 1
    beta = CENTER*pow(P**m, -1, 2**t) % 2**t
    need(beta > 0 and beta % 2 == 1, "positive odd phase")
    return dict(m=m, t=t, beta=beta, K=11*m+t,
                residue=beta*Q**m-CENTER)


def intersects(a, b):
    return (a[0]-b[0]) % 2**min(a[1], b[1]) == 0


def minimal_cylinders(rows):
    kept = []
    for r, k in sorted(set(rows), key=lambda x: (x[1], x[0])):
        if not any(k >= h and r % 2**h == s for s, h in kept):
            kept.append((r, k))
    return kept


def mass(rows):
    return sum((Fraction(1, 2**k) for _, k in rows), Fraction())


def bank():
    inherited = runpy.run_path(str(ROOT / "04-computation/experiments/entry_20260927_recursive.py"))
    entries = inherited["bank"]()
    retained = {}
    for line in (ROOT / "05-knowledge/results/reset_20260926_swaplift.out").read_text(encoding="utf-8").splitlines():
        if re.fullmatch(r"\d+(?: \d+){8}", line):
            q, j, a, prev, b, clipped, optimized, r, k = map(int, line.split())
            retained[q] = (j, a, prev, b, optimized, r, k)
    need(len(entries) == len(retained) == 171, "inherited row universe")
    for q, row in entries.items():
        need(tuple(row[key] for key in ("J", "A", "P", "b", "R", "residue", "K")) == retained[q],
             "regenerated bank differs from saved table")
    rows = minimal_cylinders((row["residue"], row["K"]) for row in entries.values())
    need(len(rows) == 65, "minimal bank universe")
    need(max(row["J"]+1 for row in entries.values()) == 41, "old selected horizon")
    need(mass(rows) == Fraction(6985206796614369409, 2**65), "old density")
    return entries, rows


def main():
    entries, old = bank()
    need(compose_word(WORD) == (P, Q, 2363), "base transfer")
    negative = [-17]
    for a in WORD:
        z = 3*negative[-1]+1
        need(v2(abs(z)) == a, "negative-cycle valuation")
        negative.append(z//2**a)
    need(negative == [-17, -25, -37, -55, -41, -61, -91, -17], "cycle")

    generated = 0
    for k in range(1, 31):
        transfer = compose_word(WORD*k)
        need(transfer == (P**k, Q**k, CENTER*(P**k-Q**k)), "repeat transfer")
        for a in (2, 4, 6, 10, 18):
            source = a*Q**k-CENTER
            n = source
            for expected in WORD*k:
                n, actual = step(n)
                need(actual == expected and n > source, "generated repeat/threshold")
            need(n == a*P**k-CENTER, "repeat endpoint")
            generated += 1
    accepted = rejected = 0
    for source in range(1, 8192, 2):
        for k in range(1, 4):
            guard = v2(source+CENTER) >= 11*k+1
            need(follows(source, WORD*k) == guard, "guard iff actual word")
            accepted += guard
            rejected += not guard

    family_checks = 0
    rows = [cylinder(m) for m in range(1, 51)]
    for row in rows:
        m, t, beta = row["m"], row["t"], row["beta"]
        for lift in (0, 1, 7, 29):
            b = beta+lift*2**t
            source = b*Q**m-CENTER
            need(source % 2**row["K"] == row["residue"], "cylinder membership")
            n = source
            actual_word = []
            for i in range(1, 7*m+1):
                n, a = step(n)
                actual_word.append(a)
                need((n < source) == (i == 7*m), "exact first descent")
            z = b*P**m-CENTER
            expected_word = WORD*(m-1)+WORD[:-1]+(4+v2(z),)
            need(tuple(actual_word) == expected_word, "full terminal word")
            need(n == z//2**v2(z), "closed terminal endpoint")
            mul, den, carry = compose_word(expected_word)
            need(mul*source+carry == den*n, "independent affine composition")
            family_checks += 1

    # One finite exact bank disjointness certificate applies to every m>=2.
    parent = ((-CENTER) % 4096, 12)
    need(all(not intersects(parent, (row["residue"], row["K"]))
             for row in entries.values()), "all raw bank rows miss parent cylinder")
    need(all(not intersects(parent, row) for row in old), "minimal bank misses parent")
    overlap = []
    for row in rows:
        own = (row["residue"], row["K"])
        inter = sum((Fraction(1, 2**max(own[1], r[1]))
                     for r in old if intersects(own, r)), Fraction())
        need(inter == (Fraction(1, 4096) if row["m"] == 1 else 0), "precise old overlap")
        overlap.append(inter)
        need(row["residue"] % 8 == 7, "minus17 source residue")
    # The entire inherited minus5 family has residue -5 mod8=3.
    need(not intersects((7, 3), (3, 3)), "disjoint from all minus5 cylinders")
    need(all(not intersects((a["residue"], a["K"]), (b["residue"], b["K"]))
             for i, a in enumerate(rows) for b in rows[i+1:]), "pairwise new disjointness")
    first_twenty = [(row["residue"], row["K"]) for row in rows[:20]]
    added = mass(first_twenty)-Fraction(1, 4096)
    union = minimal_cylinders(old+first_twenty)
    need(mass(union) == mass(old)+added, "independent exact union density")
    tail = Fraction(1, (P-1)*P**20)

    # Fixed hostile: the preceding budget cannot guarantee this selected return.
    m, b = 11, 1
    need(cylinder(m)["t"] == 2 and v2(b*P**m-CENTER) == 1, "budget boundary")
    source = b*Q**m-CENTER
    n = source
    for _ in range(7*m):
        n, _ = step(n)
        need(n > source, "weakened-budget hostile has no earlier descent")
    need(n == (b*P**m-CENTER)//2 and n >= source, "hostile endpoint")
    example = 2**22-CENTER
    trace = [example]
    for _ in range(14):
        trace.append(step(trace[-1])[0])

    report = dict(
        status="PROVED scoped families; FINITE-EXACT controls; Collatz OPEN",
        inherited_bank=dict(raw_rows=171, minimal_rows=65, max_odd_steps=41,
                            density=str(mass(old)), saved_table_agreement=True),
        cycle=negative, word=list(WORD), transfer=[P, Q, 2363],
        generated_repeat_controls=dict(k="1..30", coefficients=[2,4,6,10,18], count=generated),
        guard_controls=dict(odd_sources="1..8191", k="1..3", accepted=accepted, rejected=rejected),
        first_descent_controls=dict(m="1..50", coefficient_lifts=[0,1,7,29], count=family_checks),
        first_sixteen_cylinders=rows[:16],
        coverage=dict(old_contains_exactly_m=[1], all_m_ge_2_new=True,
                      raw_bank_misses_parent={"residue":parent[0], "K":parent[1]},
                      disjoint_from_all_inherited_minus5=True,
                      added_density_m2_to20=str(added), added_density_decimal=float(added),
                      augmented_old_density_m1_to20=str(mass(union)),
                      remaining_density_strictly_less_than=str(tail)),
        hostile=dict(m=m, required_t=2, weakened_t=1, b=b, source=source,
                     selected_step=7*m, endpoint=n),
        first_added_trace=trace,
        checks_active_under_optimized_python=True)
    print(json.dumps(report, indent=2))
    print("PASS: all checks active under -O")


if __name__ == "__main__":
    main()

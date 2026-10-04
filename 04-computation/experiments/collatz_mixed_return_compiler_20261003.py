"""Compile expanding heads plus negative cycles into guarded return cylinders.

Exact integer arithmetic, an explicit first-descent predicate, and an ordinary
parity polynomial are retained together. All checks survive Python -O.
"""
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path
import json
import re
import runpy

ROOT = Path(__file__).resolve().parents[2]
HEAD = (1, 2)
CYCLE = (1, 1, 1, 2, 1, 1, 4)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(n > 0, "positive valuation input")
    return (n & -n).bit_length()-1


def odd_step(n):
    need(n > 0 and n % 2 == 1, "positive odd source")
    z = 3*n+1
    a = v2(z)
    return z >> a, a


def word_summary(word):
    A = S = 0
    for a in word:
        need(a >= 1, "positive valuations")
        S = 3*S+2**A
        A += a
    return len(word), A, S


def all_prefixes_expand(word):
    A = 0
    for i, a in enumerate(word, 1):
        A += a
        if 3**i <= 2**A:
            return False
    return True


def cycle_data(word):
    need(bool(word) and all_prefixes_expand(word), "expanding cycle prefixes required")
    q, B, S = word_summary(word)
    P, Q = 3**q, 2**B
    c, rem = divmod(S, P-Q)
    need(rem == 0 and c > 0 and c % 2 == 1, "positive odd integer cycle center")
    n = -c
    for a in word:
        z = 3*n+1
        need(z < 0 and v2(-z) == a, "actual negative-cycle valuation")
        n = z//2**a
    need(n == -c, "negative cycle closes")
    return q, B, P, Q, c


def compile_family(head, cycle, m):
    """A sufficient compiler. Rejection means unsupported by this scheme."""
    need(m >= 1 and all_prefixes_expand(head), "positive count and expanding head")
    p, A, S = word_summary(head)
    q, B, P, Q, c = cycle_data(cycle)
    H, M = 2**A, 3**p
    C = H*c+S
    g = gcd(C, M)
    M0, C0 = M//g, C//g
    need(H*Q**m > C0 and g*Q**m > c, "uniform positive-source conditions")
    t = 1
    while True:
        D = 2**t*H*Q**m-M*P**m
        debt = 2**t*C0-M0*c
        if D > max(0, debt):
            break
        t += 1
    beta = c*pow(g*P**m, -1, 2**t) % 2**t
    residue3 = C0*pow(H*Q**m, -1, M0) % M0 if M0 > 1 else 0
    u0 = beta
    if M0 > 1:
        u0 += 2**t*((residue3-beta)*pow(2**t, -1, M0) % M0)
    K = A+B*m+t
    numerator = H*u0*Q**m-C0
    need(numerator % M0 == 0, "compiled integrality")
    n0 = numerator//M0
    need(0 < n0 < 2**K and n0 % 2 == 1, "canonical positive odd cylinder")
    return dict(head=list(head), cycle=list(cycle), m=m, p=p, A=A, S=S,
                q=q, B=B, P=P, Q=Q, center=c, forced_ternary=g,
                denominator=M0, offset=C0, t=t, parameter=u0,
                parameter_modulus=M0*2**t, residue=n0, K=K,
                first_descent=p+q*m)


def instantiate(row, lift):
    need(lift >= 0, "nonnegative coefficient lift")
    u = row["parameter"]+lift*row["parameter_modulus"]
    n = (2**row["A"]*u*row["Q"]**row["m"]-row["offset"])//row["denominator"]
    z = row["forced_ternary"]*u*row["P"]**row["m"]-row["center"]
    return u, n, z


def parity_positions(head, cycle, m):
    """Exponents in H(z), and nominal ordinary time; final extra zeros omitted."""
    positions = []
    time = 0
    for a in tuple(head)+tuple(cycle)*m:
        positions.append(time)
        time += a+1
    return positions, time


def beta_mul(pair):
    """Multiply A+B*phi by beta=phi^-1=phi-1."""
    a, b = pair
    return b-a, a


def golden_fraction(bits):
    """Exact Horner reading sum bits[j]*beta^(j+1), pair A+B*phi."""
    value = (0, 0)
    for bit in reversed(bits):
        value = beta_mul((value[0]+bit, value[1]))
    return value


def verify_instance(row, lift):
    u, source, z = instantiate(row, lift)
    need(source == row["residue"]+lift*2**row["K"], "source cylinder step")
    tau = v2(z)
    need(tau >= row["t"], "terminal budget")
    word = tuple(row["head"])+tuple(row["cycle"])*(row["m"]-1)
    word += tuple(row["cycle"][:-1])+(row["cycle"][-1]+tau,)
    n = source
    for i, expected in enumerate(word, 1):
        n, actual = odd_step(n)
        need(actual == expected, "actual valuation word")
        need((n < source) == (i == len(word)), "original-source exact first descent")
    need(n == z//2**tau, "terminal closed form")
    p, A, S = word_summary(word)
    need(3**p*source+S == 2**A*n, "independent composed transfer")
    positions, nominal_time = parity_positions(row["head"], row["cycle"], row["m"])
    ordinary_time = nominal_time+tau
    bits = [int(j in positions) for j in range(ordinary_time)]
    # Independent ordinary-map replay verifies every polynomial position.
    raw = source
    for bit in bits:
        need(raw % 2 == bit, "ordinary parity position")
        raw = 3*raw+1 if bit else raw//2
    need(raw == n, "ordinary endpoint")
    # Independent power-sum and Horner evaluation in Z[phi].
    power, value = (1, 0), (0, 0)
    selected = set(positions)
    for j in range(ordinary_time):
        power = beta_mul(power)
        if j in selected:
            value = (value[0]+power[0], value[1]+power[1])
    need(value == golden_fraction(bits), "golden polynomial evaluation")
    return dict(source=source, endpoint=n, ordinary_time=ordinary_time,
                extra_terminal_zeros=tau, golden_constant=list(value),
                golden_tail_multiplier=list(power))


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


def old_bank():
    inherited = runpy.run_path(str(ROOT / "04-computation/experiments/entry_20260927_recursive.py"))
    entries = inherited["bank"]()
    retained = {}
    for line in (ROOT / "05-knowledge/results/reset_20260926_swaplift.out").read_text(encoding="utf-8").splitlines():
        if re.fullmatch(r"\d+(?: \d+){8}", line):
            q, j, a, prev, b, clipped, optimized, r, k = map(int, line.split())
            retained[q] = (j, a, prev, b, optimized, r, k)
    need(len(entries) == len(retained) == 171, "old bank row count")
    for q, row in entries.items():
        need(tuple(row[key] for key in ("J", "A", "P", "b", "R", "residue", "K")) == retained[q], "old saved-table match")
    rows = minimal_cylinders((r["residue"], r["K"]) for r in entries.values())
    need(len(rows) == 65, "minimal old bank count")
    return entries, rows


def bounded_head_ranking(old):
    """Exact new mass against entire inherited families, not just samples.

    Candidate first-descent times are <=25. Inherited cylinders whose exact
    first-descent times exceed 25 cannot intersect any candidate cylinder.
    """
    heads = [head for p in range(1, 5)
             for head in product(range(1, 5), repeat=p)
             if all_prefixes_expand(head)]
    candidates = [(head, [compile_family(head, CYCLE, m) for m in range(1, 4)])
                  for head in heads]
    last_time = max(row["first_descent"] for _, rows in candidates for row in rows)
    baseline = list(old)
    included = []
    for head, cycle in (((), (1, 2)), ((), CYCLE), (HEAD, CYCLE)):
        retained = []
        for m in range(1, last_time//len(cycle)+1):
            row = compile_family(head, cycle, m)
            if row["first_descent"] <= last_time:
                baseline.append((row["residue"], row["K"]))
                retained.append(m)
        included.append(dict(head=list(head), cycle=list(cycle), m=retained))
    baseline = minimal_cylinders(baseline)
    ranked = []
    for head, rows in candidates:
        added = Fraction()
        for row in rows:
            candidate = (row["residue"], row["K"])
            if any(b[1] <= candidate[1] and intersects(candidate, b) for b in baseline):
                continue
            added += Fraction(1, 2**candidate[1])-sum(
                (Fraction(1, 2**b[1]) for b in baseline if intersects(candidate, b)), Fraction())
        ranked.append((added, head, rows))
    ranked.sort(key=lambda entry: (-entry[0], entry[1]))
    need(len(heads) == 13 and last_time == 25, "bounded head search universe")
    need(ranked[0][1] == (1,), "ranked best head")
    need(sum(extra > 0 for extra, _, _ in ranked) == 12, "positive head count")
    return dict(head_lengths="1..4", letters="1..4", cycle=list(CYCLE), m="1..3",
                expanding_heads=len(heads), candidate_cylinders=3*len(heads),
                last_first_descent=last_time, baseline_cylinders=len(baseline),
                inherited_family_truncations_justified_by_first_descent=included,
                ranking=[dict(head=list(head), new_density=str(extra),
                              cylinders=[{k:r[k] for k in ("m", "residue", "K", "t")} for r in rows])
                         for extra, head, rows in ranked])


def main():
    golden = runpy.run_path(str(ROOT / "04-computation/experiments/golden_digit_carry_20261003.py"))
    golden_checks = 0

    def check_external_reader(row, result):
        positions, _ = parity_positions(row["head"], row["cycle"], row["m"])
        constant = golden["evaluate_laurent"]({-j-1: 1 for j in positions})
        tail = golden["phi_power"](-result["ordinary_time"])
        need(list(constant) == result["golden_constant"], "independent Laurent reader constant")
        need(list(tail) == result["golden_tail_multiplier"], "independent Laurent reader time")

    rows = [compile_family(HEAD, CYCLE, m) for m in range(1, 51)]
    checks = 0
    for row in rows:
        need(row["forced_ternary"] == 3 and row["denominator"] == 3 and row["offset"] == 47, "mixed reduced constants")
        for lift in (0, 1, 7, 29):
            result = verify_instance(row, lift)
            check_external_reader(row, result)
            golden_checks += 1
            need(v2(result["source"]+5) == 5, "outside old minus5 family")
            need(v2(3*result["source"]+47) == 11*row["m"]+3, "source-decodable block count")
            checks += 1
    generic_checks = 0
    generic_outside_domain = []
    for head in ((), (1,), (1, 1), (1, 2), (1, 1, 2)):
        for cycle in ((1,), (1, 2), CYCLE):
            for m in range(1, 9):
                try:
                    row = compile_family(head, cycle, m)
                except ValueError as exc:
                    need(str(exc) == "uniform positive-source conditions", "unexpected generic rejection")
                    generic_outside_domain.append(dict(head=list(head), cycle=list(cycle), m=m))
                    continue
                for lift in (0, 2):
                    result = verify_instance(row, lift)
                    check_external_reader(row, result)
                    golden_checks += 1
                    generic_checks += 1
    rejected = []
    for head, cycle, reason in (((2,), CYCLE, "nonexpanding head"),
                                 (HEAD, (1, 1, 2), "noninteger cycle center"),
                                 (HEAD, (2,), "nonexpanding cycle")):
        try:
            compile_family(head, cycle, 1)
        except ValueError:
            rejected.append(reason)
        else:
            raise ValueError("hostile compiler input was accepted")
    entries, old = old_bank()
    ranking = bounded_head_ranking(old)
    parent = (667, 11)
    need(all(not intersects(parent, (r["residue"], r["K"])) for r in entries.values()), "raw bank misses mixed parent")
    need(all(not intersects(parent, r) for r in old), "minimal bank misses mixed parent")
    for row in rows:
        need(row["residue"] % 2048 == 667, "all-height parent class")
        need(row["residue"] % 8 == 3, "outside minus17 source residue")
    need(all(not intersects((a["residue"], a["K"]), (b["residue"], b["K"]))
             for i, a in enumerate(rows) for b in rows[i+1:]), "mixed classes disjoint")
    first_twenty = [(r["residue"], r["K"]) for r in rows[:20]]
    added = mass(first_twenty)
    combined = minimal_cylinders(old+first_twenty)
    need(mass(combined) == mass(old)+added, "independent union density")
    tail = Fraction(1, 9*(2187-1)*2187**20)

    # Retain only the best new head from the bounded ranking.
    best_rows = [compile_family((1,), CYCLE, m) for m in range(1, 51)]
    best_parent = (2719, 12)
    need(all(not intersects(best_parent, (r["residue"], r["K"])) for r in entries.values()), "raw bank misses winning parent")
    for row in best_rows:
        need(row["residue"] % 4096 == 2719, "winning all-height parent")
        need(row["residue"] % 8 == 7, "winning head outside minus5 and mixed")
        for lift in (0, 1, 7, 29):
            result = verify_instance(row, lift)
            check_external_reader(row, result)
            golden_checks += 1
            need(v2(result["source"]+17) == 4, "winning head outside minus17")
            need(v2(3*result["source"]+35) == 11*row["m"]+1, "winning decoder")
    need(all(not intersects((a["residue"], a["K"]), (b["residue"], b["K"]))
             for i, a in enumerate(best_rows) for b in best_rows[i+1:]), "winning classes disjoint")
    best_first_twenty = [(r["residue"], r["K"]) for r in best_rows[:20]]
    best_added = mass(best_first_twenty)
    need(mass(minimal_cylinders(combined+best_first_twenty)) == mass(combined)+best_added, "winning union mass")
    best_hostile_m, best_hostile_u = 5, 5
    best_hostile_n = (2*best_hostile_u*2048**best_hostile_m-35)//3
    best_hostile_z = best_hostile_u*2187**best_hostile_m-17
    need(best_rows[4]["t"] == 2 and v2(best_hostile_z) == 1, "winning weak budget")
    best_hostile_end = best_hostile_n
    for _ in range(7*best_hostile_m+1):
        best_hostile_end, _ = odd_step(best_hostile_end)
        need(best_hostile_end > best_hostile_n, "winning weak budget fails")
    need(best_hostile_end == best_hostile_z//2, "winning hostile endpoint")

    # Omitting the extra division when m=9 admits a wholly rising selected word.
    m, parameter = 9, 11
    need(rows[m-1]["t"] == 2 and parameter*2048**m % 3 == 1, "weak-budget guard")
    z = 3*parameter*2187**m-17
    need(v2(z) == 1, "weak-budget valuation")
    source = (8*parameter*2048**m-47)//3
    endpoint = source
    for _ in range(7*m+2):
        endpoint, _ = odd_step(endpoint)
        need(endpoint > source, "weak budget fails before and at selected return")
    need(endpoint == z//2, "hostile endpoint identity")
    example = verify_instance(rows[0], 0)
    trace = [example["source"]]
    for _ in range(9):
        trace.append(odd_step(trace[-1])[0])
    positions, nominal_time = parity_positions(HEAD, CYCLE, 1)
    report = dict(status="PROVED compiler on stated input domain; FINITE-EXACT checks; Collatz OPEN",
                  independent_golden_reader=dict(script="golden_digit_carry_20261003.py", count=golden_checks),
                  mixed_instances=dict(m="1..50", lifts=[0,1,7,29], count=checks),
                  generic_instances=dict(heads=[[],[1],[1,1],[1,2],[1,1,2]],
                                         cycles=[[1],[1,2],list(CYCLE)], m="1..8", lifts=[0,2], count=generic_checks,
                                         outside_stated_sufficient_domain=generic_outside_domain),
                  compiler_rejection_controls=rejected,
                  bounded_head_ranking=ranking,
                  retained_winning_head=dict(head=[1], instances=200, m="1..50", lifts=[0,1,7,29],
                                            missed_parent={"residue":2719,"K":12},
                                            disjoint_from_minus5_minus17_and_mixed=True,
                                            density_m1_to20=str(best_added), decimal=float(best_added),
                                            tail_strictly_less_than=str(Fraction(1,3*(2187-1)*2187**20)),
                                            first_twelve=[{k:r[k] for k in ("m","t","parameter","residue","K","first_descent")} for r in best_rows[:12]],
                                            example=verify_instance(best_rows[0],0),
                                            weak_budget_hostile=dict(m=best_hostile_m, parameter=best_hostile_u,
                                                                    required_t=2, used_t=1, source=best_hostile_n,
                                                                    selected_step=36, endpoint=best_hostile_end)),
                  first_twelve_cylinders=[{k:r[k] for k in ("m","t","parameter","parameter_modulus","residue","K","first_descent")} for r in rows[:12]],
                  coverage=dict(raw_bank_rows=171, minimal_bank_rows=65, saved_table_agreement=True,
                                missed_parent={"residue":667,"K":11},
                                disjoint_from_minus5_and_minus17=True,
                                added_density_m1_to20=str(added), decimal=float(added),
                                augmented_old_bank_density=str(mass(combined)),
                                tail_strictly_less_than=str(tail)),
                  weak_budget_hostile=dict(m=m, parameter=parameter, required_t=2, used_t=1,
                                          source=source, selected_step=7*m+2, endpoint=endpoint),
                  example=dict(**example, trace=trace, parity_positions=positions,
                               nominal_ordinary_time=nominal_time,
                               certificate="Theta(source)=golden_constant+golden_tail_multiplier*Theta(endpoint)"))
    print(json.dumps(report, indent=2))
    print("PASS: all checks active under -O")


if __name__ == "__main__":
    main()

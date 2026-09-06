"""Consumer of the separately audited full-word maximum-cost table.

Reuses the pinned, independently audited baseline pair compiler. The separate
independent refined audit imports neither producer and uses literal geometry.
"""
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
import importlib.util
import json
import sys

sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE_SHA = "1c4070e17aaf1825d07899f8b2e056d0e2f0b05f224e90997f732e079814eb3e"
BASE_SET_SHA = "a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58"
MAXIMA_SHA = "ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca"
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(message)


def canonical(x):
    return json.dumps(x, separators=(",", ":"))


def digest(x):
    return sha256(canonical(x).encode()).hexdigest()


def load():
    source = ROOT / "04-computation/third_20260906_grid_bootstrap.py"
    need(sha256(source.read_bytes()).hexdigest() == BASE_SOURCE_SHA, "pinned audited baseline implementation")
    spec = importlib.util.spec_from_file_location("refined_grid_baseline", source)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    baseline = None
    for line in (ROOT / "05-knowledge/results/third_20260906_grid_bootstrap.out").read_text().splitlines():
        if line.startswith("SCALE_SET "):
            row = json.loads(line[10:])
            if row["name"] == "profile6_coupled":
                baseline = row["survivors"]
    need(baseline is not None and len(baseline) == 8301 and digest(baseline) == BASE_SET_SHA,
         "entire audited input clock set")
    table_path = ROOT / "05-knowledge/results/third_20260906_grid_full_words.out"
    lines = table_path.read_text().splitlines()
    rows = [json.loads(line[8:]) for line in lines if line.startswith("MAXIMUM ")]
    need([row[0] for row in rows] == baseline, "exact clock alignment, no omitted input scales")
    need(digest([[t, E] for t, E, *_ in rows]) == MAXIMA_SHA, "pinned separately audited exact maxima")
    need("MAXIMA_SHA256 " + MAXIMA_SHA in lines and any(x.startswith("PASS gates=") for x in lines),
         "supplier replay completed successfully")
    raw = base.PROFILE_PATH.read_bytes()
    need(sha256(raw).hexdigest() == base.PROFILE_SHA, "same inherited full profile supplier")
    data = json.loads(raw)
    profiles = {int(k): {(c, tuple(w)) for c, w in level["profiles"]}
                for k, level in data["levels"].items()}
    return base, baseline, rows, data["levels"]["6"]["gcds"], profiles


def literal_owner(t, E, owner, profiles):
    need(len(owner) == 7 and gcd(*owner) == 1, "complete primitive maximizing state word")
    need(all(t % d == 0 for d in owner), "every owner state divides the declared clock")
    need(sum(d * (((t + 7 * d - 1) // (7 * d))) for d in owner) - t == E,
         "literal owner incidence cost attains supplied maximum")
    for k in range(1, 7):
        for I in combinations(range(7), k):
            c = gcd(*(owner[i] for i in I))
            word = tuple(sorted(gcd(c, owner[j]) for j in range(7) if j not in I))
            need((c, word) in profiles[7 - k], "literal full complement-word membership")


def main():
    base, baseline, table, allowed, profiles = load()
    pairs, literal = base.atlas()
    aggregate_set, component_set, coupled_set = [], [], []
    witnesses = {}
    full_cost = {}
    for t, E, owner, prefix_visits in table:
        literal_owner(t, E, owner, profiles)
        full_cost[t] = E
        avail, costs, entries, old_E = base.bag(t, allowed)
        need(0 <= E <= old_E, "full-word maximum refines the independent-value cost relaxation")
        es = [e for e in range(6, 0, -1) if t % e == 0 and base.aggregate(t, e) <= E]
        if not es:
            continue
        aggregate_set.append(t)
        hit_component = False
        hit_coupled = False
        cache = {}
        for e in es:
            for pair in pairs:
                A = base.pair_aggregate(t, e, pair)
                if A > E:
                    continue
                C = base.component(t, e, pair)
                need(C >= max(0, A), "individual nonnegative component credit dominates merged credit")
                if C > E:
                    continue
                hit_component = True
                p, q = pair[:2]
                dp, dq = e * gcd(t // e, p), e * gcd(t // e, q)
                key = tuple(sorted((dp, dq)))
                if key not in cache:
                    cache[key] = base.forced_excess(dp, dq, avail, costs, entries)
                paired_E = cache[key]
                if paired_E is not None and C <= min(E, paired_E):
                    hit_coupled = True
                    witnesses[t] = [E, paired_E, e, p, q, dp, dq, C]
                    break
            if hit_coupled:
                break
        if hit_component:
            component_set.append(t)
        if hit_coupled:
            coupled_set.append(t)
    sets = {"full_words_aggregate": aggregate_set, "full_words_components": component_set,
            "full_words_coupled": coupled_set}
    expected = {"full_words_aggregate": (8288, 21600), "full_words_components": (8202, 16704),
                "full_words_coupled": (8202, 16704)}
    for name, values in sets.items():
        need((len(values), max(values)) == expected[name], "complete refined finite census: " + name)
        need(values == sorted(set(values)) and set(values) <= set(baseline), "full input-set refinement")
    need(component_set == coupled_set, "additional relaxed pair cap deletes no further full-cost clocks")
    excluded = sorted(set(baseline) - set(coupled_set))
    need(len(excluded) == 99 and [t for t in coupled_set if t > 14904] == [16704],
         "all99 new exclusions and the isolated upper survivor")
    # This is a positive compatibility witness in the declared profile quotient.
    t, e, p, q = 16704, 4, 3, 308
    owner = [12, 16, 72, 58, 64, 9, 9]
    literal_owner(t, 188, owner, profiles)
    need(full_cost[t] == 188, "boundary owner attains the global full-word upper bound")
    need((e * gcd(t // e, p), e * gcd(t // e, q)) == (12, 16), "boundary pair is compatible with its owner")
    pair = next(row for row in pairs if row[:2] == (p, q))
    C = base.component(t, e, pair)
    need(C == 172 < 188, "fully paired boundary owner survives the sufficient overlap count")
    print("FINITE-EXACT: full-word cost consumer leaves8202 scales,maximum16704; no physical realizability or LRC(14) closure.")
    print("BASELINE_SET_SHA256", BASE_SET_SHA)
    print("MAXIMA_SHA256", MAXIMA_SHA)
    print("BOUNDARY", canonical([t, e, p, q, [12, 16], 188, C, owner]))
    print("EXCLUDED_FROM_BASELINE", canonical(excluded))
    print("LAST20_RELAXED_WITNESSES", canonical([[t, witnesses[t]] for t in coupled_set[-20:]]))
    for name, values in sets.items():
        print("SCALE_SET", canonical({"name": name, "count": len(values), "maximum": max(values),
                                      "sha256": digest(values), "survivors": values}))
    print("Exact gates:", GATES + base.GATES)
    print("PASS")


if __name__ == "__main__":
    main()

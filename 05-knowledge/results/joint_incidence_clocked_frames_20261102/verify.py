"""Exact, bounded certificates; no assertions disabled by python -O.

Run from any directory. Writes verification.json beside this file.
The continuum proofs are in STABILITY_42.md, not inferred from grids.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial
from pathlib import Path
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def q(x):
    return x * (1 - x)


def triangle_checks():
    grid = [F(i, 20) for i in range(21)]
    count = 0
    for p, a, b in product(grid, repeat=3):
        tau = p*a*(1-b) + (1-p)*(1-a)*b
        z = min(p, 1-p, F(7, 20))
        w = q(p)
        require(tau >= z*(q(a)+q(b))+q(a)*q(b), "triangle surplus")
        require(z-w >= w*w, "polynomial surplus")
        count += 1
    p = a = F(7, 8)
    b = F(1)
    tau = p*a*(1-b)+(1-p)*(1-a)*b
    wedge = q(p)*q(a)+q(p)*q(b)+q(a)*q(b)
    require(tau/wedge == F(64, 49) < F(4, 3), "uniform hostile")
    # Coefficient identity in (3), checked independently on rational r,s.
    for a, b in product(grid, repeat=2):
        p = F(7, 20)
        r, s = a-b, a+b-1
        lhs = 16*(p*a*(1-b)+(1-p)*(1-a)*b-p*(q(a)+q(b))-q(a)*q(b))
        rhs = F(1, 5)+F(4, 5)*s*s+F(44, 5)*r*r-(s*s-r*r)**2-F(12, 5)*r
        require(lhs == rhs, "triangle polynomial identity")
    # Local dual on an independent grid including unbounded-S representatives.
    dual_count = 0
    for s in (F(1, 6), F(4, 23), F(2, 11), F(4, 21), F(1, 4), F(1, 3)):
        for i in range(25):
            r = s*F(i, 24)
            w = r*(1-r)
            for j in range(81):
                mass = F(j, 16)
                lhs = max(w*mass, w*mass/2+mass*mass/4)+r*r*mass
                rhs = s*mass-(s/4)*(4*s-4*w)
                require(lhs >= rhs, "local incidence dual")
                dual_count += 1
    capacity_count = 0
    for p in [F(i, 40) for i in range(41)]:
        w = q(p)
        rho = min(p, 1-p, F(7, 20))-w
        for k in [F(i, 64) for i in range(17)]:
            for mass in [F(i, 16) for i in range(65)]:
                lhs = max(w*mass, w*mass/2+mass*mass/4)+rho*mass
                require(lhs >= k*mass-k*k*(1-4*w), "full quarter capacity potential")
                capacity_count += 1
    return {"triangle_grid": count, "identity_grid": len(grid)**2,
            "local_dual_grid": dual_count, "full_quarter_capacity_grid": capacity_count,
            "uniform_4_over_3_hostile_ratio": str(tau/wedge)}


def global_checks():
    """All labelled five-vertex weights in {0,1/2,1}, scaled by 320."""
    n = 5
    edges = list(combinations(range(n), 2))
    index = {e: i for i, e in enumerate(edges)}
    triangles = [(index[i, j], index[j, k], index[i, k])
                 for i, j, k in combinations(range(n), 3)]
    wedges = [(i, j) for i, e in enumerate(edges) for j, f in enumerate(edges)
              if i < j and len(set(e) & set(f)) == 1]
    count = 0
    for bits in product(range(3), repeat=len(edges)):
        uncertain = {i for i, b in enumerate(bits) if b == 1}
        remaining = set(uncertain)
        selected = []
        grouped_320 = 0
        while remaining:
            i = min(remaining)  # Every nonzero weight equals 1/4.
            incident = {j for j in remaining if set(edges[i]) & set(edges[j])}
            mass = len(incident)-1
            grouped_320 += max(20*mass, 10*mass+5*mass*mass)+8*mass
            selected.append(i)
            remaining -= incident
        tau_320 = 40*sum(bits[i]*bits[j]*(2-bits[k])
                        +(2-bits[i])*(2-bits[j])*bits[k] for i, j, k in triangles)
        wedge_320 = 20*sum(i in uncertain and j in uncertain for i, j in wedges)
        surplus_320 = 8*sum(j != i and j in uncertain and bool(set(edges[i]) & set(edges[j]))
                           for i in selected for j in range(len(edges)))
        require(tau_320 >= wedge_320+surplus_320 >= grouped_320, "joint global budget")
        count += 1
    return {"all_five_vertex_ternary_weight_tournaments": count}


def endpoint_checks():
    c0 = F(1243, 5000)
    taylor = sum(F(26, 5)**j/factorial(j) for j in range(17))
    require(taylor > F(1250, 7), "initial exponential bound")
    require(c0/F(13, 10) > F(4, 21), "small variance exclusion")
    chain = [c0, F(1, 4), F(51, 200), F(27, 100), F(8, 25),
             F(12, 25), F(3, 4), F(9, 10), F(49, 50)]
    records = []
    for n in range(42, 48):
        r = n//2
        kappa = F(4, r)*(1-F(3, 4)**r-(F(1, 4)**r if n % 2 == 0 else 0))
        bootstrap_slacks = []
        for c, next_c in zip(chain, chain[1:]):
            mass = c/kappa-n*kappa/8
            require(0 < mass < F(r, 4), "mass bound domain")
            slack = 1-(1-mass/r)**r-F(1, 4)**r-next_c
            require(slack > 0, f"bootstrap at n={n}, c={c}")
            bootstrap_slacks.append(str(slack))
        c, s = chain[-1], F(4, r)
        t = s*(1-s)
        mass = c/kappa-n*kappa/8
        heavy_bound = (4*mass-4*t*r)/(1-4*t)
        require(heavy_bound > 9, "ten heavy edges")
        Q = r+F(33, 2)+F(8, r)
        tails = [Q*F(3, 4)**(r-L)/L for L in range(1, r-9)]
        require(max(tails) < c-4*s, "tail payment")
        require(F(1, 6) <= s < F(7, 20), "dual range")
        require((s/4)*r == 1, "exact final cancellation")
        records.append({"n": n, "kappa": str(kappa), "heavy_lower_bound": str(heavy_bound),
                        "max_tail": str(max(tails)), "tail_budget": str(c-4*s),
                        "bootstrap_slacks": bootstrap_slacks})
    failures = []
    for n in (40, 41):
        r = n//2
        kappa = F(4, r)*(1-F(3, 4)**r-(F(1, 4)**r if n % 2 == 0 else 0))
        mass = c0/kappa-n*kappa/8
        first_step = 1-(1-mass/r)**r-F(1, 4)**r
        require(first_step < c0, "expected scalar bootstrap hostile")
        failures.append(n)
    return {"initial_taylor_slack": str(taylor-F(1250, 7)),
            "bootstrap_chain": [str(c) for c in chain], "orders": records,
            "certificate_only_failures": failures}


def main():
    result = {"status": "FINITE-EXACT; continuum claims proved in STABILITY_42.md",
              "triangle": triangle_checks(), "global": global_checks(),
              "endpoints": endpoint_checks()}
    path = Path(__file__).with_name("verification.json")
    path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"triangle": result["triangle"], "global": result["global"],
                      "certified_orders": [r["n"] for r in result["endpoints"]["orders"]]}, sort_keys=True))


if __name__ == "__main__":
    main()

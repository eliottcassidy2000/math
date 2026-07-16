#!/usr/bin/env python3
"""
THE GENERAL-CLUSTER RESONANCE LAW — the multi-comb assembly
boxeph-2026-07-16-S26 (executes THM-886's named next step)

For ANY cluster E (section s, swing set R_s, P = 7 lcm E), decompose the signed
endpoint sum by owners: S(n) = sum_e S_e(n).

 (I)   PER-OWNER COMB IDENTITY (exact, every owner, every integer a):
         S_e(e·a) = nu_e(a) := sum_{c mod 7} N_c^{(e)} e(ac/7),
       7-periodic tooth amplitudes; N_c^(e) = signed count of e-owned active
       endpoints with boundary index j == c (mod 7).  [phase collapse:
       e(ea·j/(7e)) = e(aj/7).]
 (II)  PER-OWNER OFF-COMB PROFILE:
         |S_e(m)| <= min(M_e, R_e / sin(pi ||m/e||))   (m not on the comb),
       R_e = total run count of the per-(class, sign) active index sets in Z_e
       (each run is a geometric block at frequency m/e).
 (III) MULTI-COMB PROFILE: |S(m)| <= sum_e profile_e(m) pointwise on Z_P.
 (IV)  COMB TRANSVERSALITY: combs eZ and e'Z meet exactly on lcm(e,e')Z; CRT
       controls compound near-hits.
 (V)   THE ASSEMBLED LAW: with |khat_P(n)| <= 1/(2 P^2 sin^2(pi n/P)),
         Q_s(w) <= 2 pi^2 sum_{l != 0 mod P} khat-bound(l) [sum_e profile_e(lw)]^2,
       and the resonance classifier is the DIOPHANTINE position of w against
       the cluster's own comb system: w is resonant for owner e at depth l iff
       ||l w / e|| is small (in Z_P units) — i.e. Q_s(w) is governed by how
       LONELY w is against E itself at mode scale.  The recursion is structural.

Battery: family sanity + multi-large-owner + balanced (compact-core-like)
+ near-AP clusters; exactness referees + full-Z_P profile scans + the law's
cover on resonant/compound/random w.

Pure Python. Reuses S25's endpoint machinery.
"""

import sys
from math import gcd, lcm, pi, sin
import cmath

sys.path.insert(0, '04-computation')
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints, Qs_bilinear


def owner_data(E, s):
    """Return P, M, and per-owner dicts: positions (as P-integers), signs,
    N_c 7-vectors, run counts R_e, M_e."""
    P = 7 * lcm(*E)
    pts = endpoints(E, s)
    data = {}
    for p, sg, o in pts:
        d = data.setdefault(o, {"pos": [], "sgn": [], "N": [0] * 7, "runs": None})
        j = p * 7 * o
        assert j.denominator == 1
        j = int(j)
        d["pos"].append(int(p * P))
        d["sgn"].append(sg)
        d["N"][j % 7] += sg
        d.setdefault("idx", {}).setdefault((j % 7, sg), set()).add(j // 7)
    for e, d in data.items():
        R = 0
        for (c, sg), st in d.get("idx", {}).items():
            R += sum(1 for m in st if (m - 1) % e not in st)
        d["runs"] = R
        d["M"] = len(d["pos"])
    return P, len(pts), data


def S_owner(d, n, P):
    return sum(sg * cmath.exp(2j * pi * (n % P) * q / P)
               for sg, q in zip(d["sgn"], d["pos"]))


def nu_hat(d, a):
    return sum(d["N"][c] * cmath.exp(2j * pi * a * c / 7) for c in range(7))


def profile_e(e, d, m, P):
    """the (I)/(II) upper profile for owner e at frequency m."""
    r = (m * e) % P
    if r % P == 0 or (m % P) * e % P == 0:
        pass
    # distance of m/e to the nearest integer, in the Z_P sense:
    # m/e mod 1 with m in Z_P: ||m/e|| = dist(m mod e-ish...) — compute via
    # the exact rational m/e mod 1:
    x = (m % P) / e % 1.0
    dist = min(x, 1 - x)
    if dist < 1e-12:               # on the comb
        a = round((m % P) / e)
        return min(d["M"], abs(nu_hat(d, a)) + 1e-9)
    return min(d["M"], d["runs"] / sin(pi * dist))


def referee_cluster(E, s, name, full_scan=True, wbattery=()):
    P, M, data = owner_data(E, s)
    if M == 0:
        print(f"  [{name}] R_s empty at section {s}; skipped")
        return
    owners = sorted(data)
    print(f"  [{name}] s={s} P={P} M={M}; owners: " +
          ", ".join(f"{e}(M={data[e]['M']},R={data[e]['runs']})" for e in owners))

    # (I) exact comb identities
    ok1 = True
    for e in owners:
        for a in (1, 2, 3, 5, 7, 11):
            lhs = S_owner(data[e], e * a, P)
            if abs(lhs - nu_hat(data[e], a)) > 1e-8:
                ok1 = False
    print(f"      (I) per-owner comb identity S_e(ea) = nu_e(a): {ok1}")

    # (II)+(III) profile scan
    if full_scan and P <= 20000:
        pos_all = [q for e in owners for q in data[e]["pos"]]
        sgn_all = [sg for e in owners for sg in data[e]["sgn"]]
        viol2 = viol3 = 0
        worst3 = 0.0
        for m in range(1, P):
            tot = 0.0
            Sm = 0j
            okm = True
            for e in owners:
                se = S_owner(data[e], m, P)
                pe = profile_e(e, data[e], m, P)
                if abs(se) > pe + 1e-6:
                    viol2 += 1
                    okm = False
                tot += pe
                Sm += se
            if abs(Sm) > tot + 1e-6:
                viol3 += 1
            worst3 = max(worst3, abs(Sm) / max(tot, 1e-12))
        print(f"      (II) per-owner profile: violations {viol2}; "
              f"(III) summed profile: violations {viol3}, worst ratio {worst3:.3f}")
    else:
        print(f"      (II/III) P = {P} too large for full scan; combs + samples only")
        viol = 0
        for e in owners:
            for a in range(1, 30):
                if abs(S_owner(data[e], e * a, P)) > min(
                        data[e]["M"], abs(nu_hat(data[e], a))) + 1e-6:
                    viol += 1
        print(f"      comb-point checks: violations {viol}")

    # (V) assembled law on a w battery
    for w in wbattery:
        q = Qs_bilinear(E, s, w)
        # bound: C1*M low-band absorber + sum over l of khat-bound * profile^2
        tot = 12.0 * M
        L = min(P // 2, 4000)
        for l in range(1, L + 1):
            m = (l * w) % P
            if m == 0:
                continue
            prof = sum(profile_e(e, data[e], m, P) for e in owners)
            nt = min(l, P - l)
            tot += 2 * pi * pi * (1 / (2 * P * P * sin(pi * l / P) ** 2)) * prof * prof
        # resonance score: the Diophantine position of w against the combs
        score = 0.0
        for l in range(1, 31):
            for e in owners:
                x = ((l * w) % P) / e % 1.0
                dd = min(x, 1 - x)
                amp = (abs(nu_hat(data[e], round((l * w % P) / e))) if dd < 1e-9
                       else data[e]["runs"] / sin(pi * dd))
                score = max(score, min(data[e]["M"], amp) ** 2 / (8 * l * l))
        print(f"      w={w:>6} gcd={gcd(w,P):>3}: Q_s/M = {q/M:8.2f}; "
              f"law bound/M = {tot/M:9.1f} (cover {q <= tot + 1e-6}); "
              f"resonance score/M = {score/M:7.2f}")


if __name__ == "__main__":
    print("=" * 76)
    print("THE GENERAL-CLUSTER RESONANCE LAW -- multi-comb referees")

    # A: family sanity (t = 60: P = 420)
    referee_cluster([1, 2, 3, 4, 5, 6, 60], 0, "family {1..6,60}",
                    wbattery=(11, 59, 61, 120 + 7, 187))

    # B: two mid owners
    referee_cluster([1, 2, 3, 4, 5, 36, 60], 0, "two-owner {1..5,36,60}",
                    wbattery=(11, 37, 71, 36 * 3 + 1, 60 * 2 + 5))

    # C: two large owners sharing a factor
    referee_cluster([1, 2, 3, 4, 5, 56, 84], 0, "two-large {1..5,56,84}",
                    wbattery=(11, 57, 85, 56 * 2 + 3, 84 * 3 + 1))

    # D: balanced compact-core-like (no small runners)
    for s in range(7):
        referee_cluster([12, 15, 20, 21, 28, 30, 35], s,
                        "balanced {12,15,20,21,28,30,35}",
                        wbattery=(11, 421, 12 * 5 + 1) if s == 0 else ())

    # E: near-AP mid cluster
    referee_cluster([8, 9, 10, 12, 14, 15, 18], 0, "near-AP {8,9,10,12,14,15,18}",
                    wbattery=(11, 9 * 7 + 1, 2521))

    print()
    print("done")

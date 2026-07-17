#!/usr/bin/env python3
"""
THE VALUATION-CONDUCTOR PAIRING LAW (boxeph-2026-07-17-S62)
Owner directive: prove the 7-part concentration law (LEM-032(E)'s named open).

THE THEOREM (LEM-033).  Fix a class g | P, q = P/g, and a prime p with
alpha = v_p(q) >= 1.  Expand the cross-spectrum on the class through the
cross-pair sum X(gu) = sum_{ordered cross pairs (k,k')} eps_k eps_{k'}
e(u t/q), t = (p_k - p_{k'}) mod q, and grade the pairs by j = v_p(t)
(capped at alpha; t = 0 has j = alpha).  For a character chi mod q with
beta = v_p(cond chi):

  (i)  beta >= 1:  X-hat_g(chi) sees EXACTLY the grade j = alpha - beta.
  (ii) beta  = 0:  X-hat_g(chi) sees EXACTLY the grades j in {alpha-1, alpha}.

Proof (both directions, 3 lines each):
  (kernel averaging)  chi is trivial on K_gamma = {1 + (q/p^gamma) h} for
  every gamma <= alpha - beta (v_p of the kernel modulus >= 1 keeps K_gamma
  a clean group of units); averaging e(u t/q) over K_gamma multiplies it by
  [p^gamma | t]; so chi sees only j >= alpha - beta.
  (depth kill)  a pair of grade j has u_p-dependence of additive depth
  d = alpha - j in the p-component U_{p^alpha}; the component character has
  conductor p^beta; the component sum vanishes unless d = beta exactly
  (depth < beta: chi_p nontrivial on the stabilizing subgroup 1 + p^d Z;
  depth > beta: chi_p constant on cosets where the additive character is a
  nontrivial average -- for beta = 0 the depth-1 term survives with the
  Ramanujan weight, depths >= 2 vanish by mu(p^d) = 0).  So j <= alpha-beta
  for beta >= 1, and j >= alpha - 1 for beta = 0.  Intersect with (i).  QED

CONSEQUENCE (the concentration law): the conductor profile of the frame
spectrum IS the joint p-adic valuation profile of the signed cross-difference
multiset.  Top census masses sit at conductors with FULL 7-part because the
grade-0 (7-coprime) cross differences dominate: the owner-type/section-class
selection rules (Part 3) say which pairs those are.

Parts:
  1. per-prime marginal pairing referee -- ALL characters of U_q, all primes
     p | q, all grade cells: incompatible cells vanish, partition reassembles
     (main clusters at g = 1, plus deeper classes, plus 2 stress geometries).
  2. joint multi-prime cells on a character subsample (the full conductor <->
     joint-valuation pairing).
  3. owner-type / section-class selection rules at p = 7 (elementary
     valuations): MIX pairs (7-owner x 7-free) have grade >= 1 iff the
     7-owner endpoint's section-boundary class c' == 0; FF pairs are deep
     (grade >= v_7-part of both positions); MM grade >= 1 is the twisted
     class-resonance c e' == c' e (mod 7).  Verified against the data.
  4. census shares: c-hat mass by beta_7 (deficient fraction per cluster)
     next to the graded pair energies (Parseval per grade).
  5. cross-pair Gauss expansion spot referee: phi(q) X-hat_g(chi) =
     sum_pairs eps eps' G_q(chi-bar, t), G the generalized Gauss sum; plus
     one full closed-form assembly c-hat(chi) = sum_g (2/g^2) L_{P/g}(2,chi)
     X-hat_g(chi) with X-hat_g from compatible cells only.

Pure Python.  Reuses S61 (UGroup, cluster_spectra, cross_all, chat_direct,
L-values) and S26 (owner_data).
"""

import sys
from math import gcd, pi, sin
import cmath

sys.path.insert(0, '04-computation')
from lrc14_parity_lvalue_boxeph_S61 import (UGroup, cluster_spectra, cross_all,
                                            chat_direct, factorize, divisors,
                                            L2_primitive, euler_restrict)
from lrc14_general_resonance_law_boxeph_S26 import owner_data


def vp(n, p, cap):
    if n == 0:
        return cap
    v = 0
    while n % p == 0 and v < cap:
        n //= p
        v += 1
    return v


def cross_pairs(E, s):
    """ordered cross-owner pairs: (ss, delta mod P, e, e', c, c') with
    c = section-boundary class (j mod 7) of the first endpoint, c' of the
    second."""
    P, M, data = owner_data(E, s)
    owners = sorted(data)
    pts = []
    for e in owners:
        d = data[e]
        for sg, q in zip(d["sgn"], d["pos"]):
            j = q * 7 * e // P          # boundary index: q = j P/(7e)
            assert j * P == q * 7 * e
            pts.append((e, sg, q, j % 7))
    pairs = []
    for (e, sg, q, c) in pts:
        for (ee, sg2, q2, c2) in pts:
            if e == ee:
                continue
            pairs.append((sg * sg2, (q - q2) % P, e, ee, c, c2))
    return P, data, owners, pairs


def pairing_referee(name, P, X, pairs, g, subsample=40, quiet=False):
    """Part 1+2 on class g: marginal cells all chars, joint cells subsample."""
    q = P // g
    G = UGroup(q)
    root = [cmath.exp(2j * pi * t / q) for t in range(q)]
    al = {p: vp(q, p, 99) for p in (2, 3, 5, 7) if q % p == 0}
    primes = sorted(al)
    # grade every pair; build per-pair phase vectors ONCE, add into cells
    U = G.units
    nU = len(U)
    marg = {p: [[0j] * nU for _ in range(al[p] + 1)] for p in primes}
    joint = {}
    npairs = 0
    for (ss, dP, e, ee, c, cc) in pairs:
        t = dP % q
        npairs += 1
        vec = [ss * root[(u * t) % q] for u in U]
        js = tuple(vp(t, p, al[p]) for p in primes)
        for pi_, p in enumerate(primes):
            cell = marg[p][js[pi_]]
            for i in range(nU):
                cell[i] += vec[i]
        if js not in joint:
            joint[js] = [0j] * nU
        cj = joint[js]
        for i in range(nU):
            cj[i] += vec[i]
    tol = 1e-6 * max(npairs, 1)
    # direct X-hat for every character + marginal cell checks
    worst_inc = 0.0
    worst_part = 0.0
    worst_joint = 0.0
    nchecks = 0
    chi_cache = {}
    for ci, js_char in enumerate(G.chars()):
        chibar = [G.chi(js_char, u).conjugate() for u in U]
        direct = sum(X[(g * u) % P] * chibar[i] for i, u in enumerate(U)) / nU
        f = G.conductor(js_char)
        betas = {p: vp(f, p, 99) for p in primes}
        tot_cells = 0j
        for p in primes:
            a, b = al[p], betas[p]
            comp = {a - b} if b >= 1 else ({a - 1, a} if a >= 1 else {0})
            for j in range(a + 1):
                T = sum(marg[p][j][i] * chibar[i] for i in range(nU)) / nU
                if j not in comp:
                    worst_inc = max(worst_inc, abs(T))
                    nchecks += 1
            if p == primes[0]:
                tot_cells = sum(
                    sum(marg[p][j][i] * chibar[i] for i in range(nU)) / nU
                    for j in range(a + 1))
        worst_part = max(worst_part, abs(tot_cells - direct))
        # joint cells on a subsample
        if ci % max(1, len(G.chars()) // subsample) == 0:
            pred = 0j
            for js, cj in joint.items():
                ok = True
                for pi_, p in enumerate(primes):
                    a, b = al[p], betas[p]
                    comp = {a - b} if b >= 1 else ({a - 1, a} if a >= 1
                                                   else {0})
                    if js[pi_] not in comp:
                        ok = False
                        break
                T = sum(cj[i] * chibar[i] for i in range(nU)) / nU
                if ok:
                    pred += T
                else:
                    worst_joint = max(worst_joint, abs(T))
            worst_part = max(worst_part, abs(pred - direct))
    if not quiet:
        print(f"  [{name} g={g}] q={q} phi={nU} pairs={npairs} "
              f"primes {al}: incompatible-cell max = {worst_inc:.2e} "
              f"({nchecks} cells x chars); partition/joint max err = "
              f"{worst_part:.2e}; joint incompatible max = {worst_joint:.2e}")
    assert worst_inc < tol and worst_part < tol and worst_joint < tol
    return marg, al, G


def selection_rules(name, P, pairs, alpha7):
    """Part 3: owner-type x grade table at p=7 + THE UNIVERSAL SHALLOW-
    RESONANCE RULE: t mod 7 = c a_e - c' a_{e'} with a_e = (P/(7e)) mod 7,
    hence grade_7 >= 1  iff  c a_e == c' a_{e'} (mod 7).  Corollaries:
    MIX (one a = 0): grade >= 1 iff the shallow side's class c' == 0;
    FF-deep (both a = 0): grade >= 1 always."""
    from collections import defaultdict
    tab = defaultdict(lambda: [0, 0.0])          # (type, j7) -> [count, signed]
    rule_viol = 0
    for (ss, dP, e, ee, c, cc) in pairs:
        j7 = vp(dP % P, 7, alpha7)
        a1 = (P // (7 * e)) % 7
        a2 = (P // (7 * ee)) % 7
        t1 = "7" if e % 7 == 0 else "F"
        t2 = "7" if ee % 7 == 0 else "F"
        typ = "MIX" if t1 != t2 else ("MM" if t1 == "7" else "FF")
        tab[(typ, j7)][0] += 1
        tab[(typ, j7)][1] += ss
        if (j7 >= 1) != ((c * a1 - cc * a2) % 7 == 0):
            rule_viol += 1
    print(f"  [{name}] owner-type x grade-7 table (count, signed sum):")
    for typ in ("MIX", "FF", "MM"):
        row = "    " + typ + ":  " + "  ".join(
            f"j={j}: {tab[(typ, j)][0]}/{tab[(typ, j)][1]:+.0f}"
            for j in range(alpha7 + 1) if (typ, j) in tab)
        if any((typ, j) in tab for j in range(alpha7 + 1)):
            print(row)
    print(f"    UNIVERSAL RULE (grade>=1 iff c a_e == c' a_e' mod 7, "
          f"a_e = (P/7e) mod 7): violations {rule_viol}")
    assert rule_viol == 0


def census_shares(name, P, X, W, marg7, alpha7):
    """Part 4: c-hat mass by beta_7 + graded pair energies."""
    G = UGroup(P)
    cross = cross_all(P, X, W, G.units)
    triv = tuple(0 for _ in G.orders)
    mass_by_beta = {}
    tot = 0.0
    for js in G.chars():
        if js == triv:
            continue
        m = abs(chat_direct(G, js, cross)) ** 2
        b = vp(G.conductor(js), 7, 99)
        mass_by_beta[b] = mass_by_beta.get(b, 0.0) + m
        tot += m
    nU = len(UGroup(P).units)
    print(f"  [{name}] c-hat mass by beta_7 (conductor 7-level): " + ", ".join(
        f"beta={b}: {mass_by_beta.get(b, 0.0) / tot * 100:.1f}%"
        for b in range(alpha7 + 1)))
    defic = sum(v for b, v in mass_by_beta.items() if b < alpha7) / tot
    print(f"      deficient-conductor fraction (beta_7 < {alpha7}) = "
          f"{defic * 100:.1f}%")
    Uq = UGroup(P).units
    en = []
    for j in range(alpha7 + 1):
        en.append(sum(abs(z) ** 2 for z in marg7[j]) / len(Uq))
    s = sum(en)
    print("      graded pair energy (class g=1, Parseval/phi): " + ", ".join(
        f"j={j}: {e / s * 100:.1f}%" for j, e in enumerate(en)))
    return defic


def gauss_spot(name, E, s, nchi=2):
    """Part 5: phi(q) X-hat = sum eps eps' G_q(chi-bar,t) at g=1, plus one
    full closed-form assembly on the smallest cluster."""
    P, data, owners, pairs = cross_pairs(E, s)
    _, _, _, X, W = cluster_spectra(E, s)
    G = UGroup(P)
    root = [cmath.exp(2j * pi * t / P) for t in range(P)]
    evens = [js for js in G.chars()
             if G.parity(js) == 1 and G.conductor(js) > 1]
    evens.sort(key=lambda js: (G.conductor(js), G.order_of(js)))
    sel = [evens[0], evens[len(evens) // 2]][:nchi]
    for js in sel:
        chibar = {u: G.chi(js, u).conjugate() for u in G.units}
        direct = sum(X[u] * chibar[u] for u in G.units)
        Gcache = {}
        tot = 0j
        for (ss, dP, e, ee, c, cc) in pairs:
            t = dP % P
            if t not in Gcache:
                Gcache[t] = sum(chibar[u] * root[(u * t) % P]
                                for u in G.units)
            tot += ss * Gcache[t]
        err = abs(direct - tot) / (1 + abs(direct))
        print(f"  [{name}] Gauss expansion at cond={G.conductor(js)} "
              f"ord={G.order_of(js)}: phi X-hat direct = "
              f"{direct / len(G.units):+.4f}; pair-Gauss sum = "
              f"{tot / len(G.units):+.4f} (rel err {err:.1e})")
        assert err < 1e-9
    # one full assembly: c-hat(chi) via compatible cells only, weights = L-values
    js = sel[0]
    f = G.conductor(js)
    Lstar = L2_primitive(G, js, f)
    cross = cross_all(P, X, W, G.units)
    direct = chat_direct(G, js, cross)
    tot = 0j
    for g in divisors(P):
        if g == P:
            continue
        q = P // g
        if q % f != 0:
            continue
        Gq_al = {p: vp(q, p, 99) for p in (2, 3, 5, 7) if q % p == 0}
        betas = {p: vp(f, p, 99) for p in Gq_al}
        Uq = [u for u in range(1, q + 1) if gcd(u, q) == 1]
        rootq = [cmath.exp(2j * pi * t / q) for t in range(q)]
        xh = 0j
        for (ss, dP, e, ee, c, cc) in pairs:
            t = dP % q
            ok = True
            for p, a in Gq_al.items():
                b = betas[p]
                j = vp(t, p, a)
                comp = {a - b} if b >= 1 else ({a - 1, a} if a >= 1 else {0})
                if j not in comp:
                    ok = False
                    break
            if not ok:
                continue
            for u in Uq:
                xh += ss * rootq[(u * t) % q] * G.chistar(js, u, f).conjugate()
        xh /= len(Uq)
        tot += (2.0 / (g * g)) * euler_restrict(G, js, f, q, Lstar) * xh
    err = abs(direct - tot) / (1 + abs(direct))
    print(f"  [{name}] FULL closed form (L-values x compatible cells only): "
          f"c-hat direct = {direct:+.4f}; assembled = {tot:+.4f} "
          f"(rel err {err:.1e})")
    assert err < 1e-7


if __name__ == "__main__":
    print("THE VALUATION-CONDUCTOR PAIRING LAW -- LEM-033 referee (boxeph S62)")
    print("=" * 78)
    print("PART 1+2 -- the pairing (incompatible cells vanish; partition exact)")
    mains = [([12, 15, 20, 21, 28, 30, 35], 0, "balanced"),
             ([1, 2, 3, 4, 5, 36, 60], 0, "two-owner"),
             ([1, 2, 3, 4, 5, 6, 60], 0, "family")]
    kept = {}
    for E, s, name in mains:
        P, data, owners, pairs = cross_pairs(E, s)
        _, _, _, X, W = cluster_spectra(E, s)
        # sanity: pair sum rebuilds X (20 samples)
        for m in range(1, P, max(1, P // 20)):
            xm = sum(ss * cmath.exp(2j * pi * m * dP / P)
                     for (ss, dP, e, ee, c, cc) in pairs)
            assert abs(xm.imag) < 1e-6 and abs(xm.real - X[m]) < 1e-6
        marg, al, Gq = pairing_referee(name, P, X, pairs, 1)
        kept[name] = (P, X, W, pairs, marg, al)
        # a deeper class with alpha_7 >= 1
        for g in (5, 3, 2):
            if (P // g) % 7 == 0 and P % g == 0:
                pairing_referee(name, P, X, pairs, g)
                break
    print()
    print("PART 1s -- stress geometries (odd-q classes; fresh 7-owner mixes)")
    for E, s, g, name in [([8, 9, 10, 12, 14, 15, 18], 0, 8, "near-AP"),
                          ([1, 2, 3, 4, 5, 56, 84], 0, 8, "two-large")]:
        P, data, owners, pairs = cross_pairs(E, s)
        _, _, _, X, W = cluster_spectra(E, s)
        pairing_referee(name, P, X, pairs, g)
    print()
    print("PART 3 -- owner-type / section-class selection rules at p = 7")
    for E, s, name in mains + [([8, 9, 10, 12, 14, 15, 18], 0, "near-AP"),
                               ([1, 2, 3, 4, 5, 56, 84], 0, "two-large")]:
        P, data, owners, pairs = cross_pairs(E, s)
        selection_rules(name, P, pairs, vp(P, 7, 99))
    print()
    print("PART 4 -- census shares by beta_7 vs graded pair energy")
    for name in kept:
        P, X, W, pairs, marg, al = kept[name]
        census_shares(name, P, X, W, marg[7], al[7])
    print()
    print("PART 5 -- Gauss expansion + full closed-form assembly (family)")
    gauss_spot("family", [1, 2, 3, 4, 5, 6, 60], 0)
    print("=" * 78)
    print("done")

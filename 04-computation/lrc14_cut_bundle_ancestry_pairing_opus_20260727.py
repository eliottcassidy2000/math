#!/usr/bin/env python3
"""Hostile-audited exact EXTERNAL CONTROL for the THM-2508 Section 6
bilinear pairing.

This file was first promoted as a physical pairing on one THM-2471/2478
Boolean ancestry cospan.  That interpretation is false.  Its two inputs are
independently constructed data:

* a synthetic weighted root-correlation control copied from THM-2471's exact
  companion (not an exhibited physical row of THM-2471); and
* THM-2506's canonical two-row defect, which belongs to the already-empty
  high-septimal branch.

There is no common base point, atom product, owner label, or physical morphism
joining them.  THM-2478 (21a)-(21b) licenses equality of a printed future
collision-colour label with a fixed old target-component label.  It does NOT
identify that temporal character with THM-2506's static cut-row character.
Thus the summation label ``alpha`` below is an external algebraic alignment,
not a typed-sheet bridge.

All arithmetic is exact: integers, Fractions, and reduction modulo Phi_91
(resp. Phi_13) in Z[zeta_91].  Embeddings: zeta_13 = zeta_91^7 and
xi_7 = zeta_91^13.  Q(zeta_91) is a field, so a reduced vector is zero iff
the cyclotomic element is zero.

PART 1  Synthetic THM-2471 companion control: service 169*I_r, twelve
        nonzero root colours, signed ledger, energy floor, and separately
        calculated owner/deep margins.
PART 2  THM-2508 cut bank on THM-2506's canonical two-row defect: all 5,184
        primitive modes nonzero and factorized; 864 zero-cut-character and
        504 invariant-readout hostile controls.
PART 3  Typed-sheet obstruction: temporal root u and static cut row h are
        different sheets.  The cut torsor is not transported by the graft.
PART 4  External pairing
          P_{tau,a}(beta) = sum_{alpha in F_13^*}
                            Chat(-alpha) Psi_{tau,a}(alpha,beta).
        The displayed controls give 432/432 nonzero placements in both leg
        conventions.  K_{0,beta}=0 excludes the trivial Fourier character
        alpha=0; it says nothing about the physical coordinate value C(0).
        A Boolean hostile C=delta_2 has full nontrivial root-colour support,
        while the defect has full primitive-mode support, yet 108/432
        pairings vanish in each convention (36/36 at dual tau=4).  Full
        factor support therefore does not prevent cross-alpha cancellation.
        A separate direct evaluator checks the displayed 432/432 result at
        five split primes in addition to the Q(zeta_91) proof.
PART 5  Separate finite-grid graft, target-neutrality, and owner-loop controls.
        These margins are not factors of the PART 4 external product.
PART 6  Exact covariance, kappa dressing, and both torsor collapses for the
        external product.  This is equivariant algebraic data, not a physical
        current.
PART 7  Separate forward-graft/rebase boundary replay.

Conventions:
 C1  The numeric root object is THM-2471's companion
     ``weighted_root_control`` reconstructed verbatim: a synthetic exact
     positive control, not a canon-realized physical row.
 C2  The ledger anchors J=Chat_unnorm/169; fixed nonzero scaling does not
     affect nonvanishing.
 C3  Cut alpha is numerically aligned with the root-control colour alpha.
     No theorem identifies their typed sheets.  Both Chat(-alpha) (dual) and
     Chat(+alpha) (untwisted) are reported.
 C4  The separate finite-grid graft uses a scalar rational old current; it
     does not turn PART 4 into a physical product.
"""

from fractions import Fraction as F

P13 = 13
P7 = 7
N91 = 91
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------- cyclotomics

def poly_div_exact(numerator, denominator):
    remainder = list(numerator)
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator):
        coefficient = remainder[-1]
        shift = len(remainder) - len(denominator)
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            remainder[index + shift] -= coefficient * value
        while remainder and remainder[-1] == 0:
            remainder.pop()
    require(not remainder, "exact cyclotomic division")
    return quotient


def cyclotomic_polynomial(order, cache):
    if order in cache:
        return cache[order]
    polynomial = [-1] + [0] * (order - 1) + [1]
    for divisor in range(1, order):
        if order % divisor == 0:
            polynomial = poly_div_exact(polynomial, cyclotomic_polynomial(divisor, cache))
    cache[order] = polynomial
    return polynomial


PHI91 = cyclotomic_polynomial(N91, {1: [-1, 1]})
DEG91 = len(PHI91) - 1
require(DEG91 == 72, "Phi_91 degree")


def reduce91(raw):
    """Reduce a length<=91 coefficient list modulo Phi_91 -> 72-tuple."""
    coefficients = list(raw) + [0] * (N91 - len(raw))
    for exponent in range(N91 - 1, DEG91 - 1, -1):
        leading = coefficients[exponent]
        if leading:
            shift = exponent - DEG91
            for index, value in enumerate(PHI91):
                coefficients[index + shift] -= leading * value
    require(all(value == 0 for value in coefficients[DEG91:]), "remainder degree")
    return tuple(coefficients[:DEG91])


ZERO72 = tuple(0 for _ in range(DEG91))


def dict_to_reduced(poly):
    raw = [0] * N91
    for exponent, value in poly.items():
        raw[exponent % N91] += value
    return reduce91(raw)


def dmul(left, right):
    answer = {}
    for le, lv in left.items():
        for re, rv in right.items():
            e = (le + re) % N91
            answer[e] = answer.get(e, 0) + lv * rv
    return {e: v for e, v in answer.items() if v}


def z13(exponent):
    """zeta_13^exponent as a zeta_91 monomial dict."""
    return {(7 * (exponent % P13)) % N91: 1}


def x7(exponent):
    """xi_7^exponent as a zeta_91 monomial dict."""
    return {(13 * (exponent % P7)) % N91: 1}


def phi13_reduce(values, colour):
    """Reduce sum_s values[s] zeta_13^(-colour*s) modulo Phi_13 -> 12-tuple."""
    raw = [F(0) for _ in range(P13)]
    for s, value in enumerate(values):
        raw[(-colour * s) % P13] += value
    return tuple(raw[j] - raw[12] for j in range(12))


# ------------------------------------------------- PART 1: the THM-2471 stalk

def build_stalk():
    """Synthetic weighted control from THM-2471's exact companion.

    Its 4 strata, root u in F_13, 5 Boolean sheets, and 7 atom labels verify
    the displayed stalk identities, but no theorem identifies this numeric
    object with a physical THM-2471 row.
    """
    strata, sheets, atoms = 4, 5, 7
    fa = [[[0] * sheets for _ in range(P13)] for _ in range(strata)]
    ea = [[[0] * sheets for _ in range(P13)] for _ in range(strata)]
    f_atom = [[[(2 * y + u + a) % atoms for a in range(sheets)]
               for u in range(P13)] for y in range(strata)]
    e_atom = [[[(3 * y + 2 * u + a) % atoms for a in range(sheets)]
               for u in range(P13)] for y in range(strata)]

    for y in range(strata):
        f_roots = ((2 * y) % P13, (2 * y + 4) % P13)
        e_roots = ((2 * y + 1) % P13, (2 * y + 7) % P13)
        for a in range(1 + y % 4):
            fa[y][f_roots[0]][a] = 1
        for a in range(1 + (y + 2) % 4):
            fa[y][f_roots[1]][a] = 1
        for a in range(1 + (y + 1) % 4):
            ea[y][e_roots[0]][a] = 1
        for a in range(1 + (y + 3) % 4):
            ea[y][e_roots[1]][a] = 1

    u_w = [[F(sum(fa[y][r]), sheets) for r in range(P13)] for y in range(strata)]
    v_w = [[F(sum(ea[y][r]), sheets) for r in range(P13)] for y in range(strata)]
    require(any(x not in (0, 1) for row in u_w + v_w for x in row),
            "stalk weights must be genuinely non-Boolean marginals")
    require(all(u_w[y][r] * v_w[y][r] == 0 for y in range(strata) for r in range(P13)),
            "same-root disjointness A(y,u)F(y,u)=0 (THM-2471 eq (7))")

    corr_stratum = [[sum(u_w[y][(r + s) % P13] * v_w[y][r] for r in range(P13))
                     for s in range(P13)] for y in range(strata)]
    correlation = [sum(corr_stratum[y][s] for y in range(strata)) / strata
                   for s in range(P13)]
    service = sum(correlation)
    collision = service / (P13 * P13)   # I_r; service = 169*I_r (eq (7))
    require(correlation[0] == 0, "C(0)=0 same-root disjointness coordinate (eq (14))")
    require(service > 0, "positive service M=169*I_r (eq (7))")
    require(service == sum((sum(u_w[y]) * sum(v_w[y]) for y in range(strata)),
                           F(0)) / strata, "service identity (eq (6)-(7))")

    # Twelve nonzero colours (eq (9)) and the signed ledger (eq (10)).
    chat = {k: phi13_reduce(correlation, k) for k in range(P13)}
    for k in range(1, P13):
        require(any(chat[k]), f"root colour {k} vanished (eq (9))")
    require(chat[0] == tuple([service] + [F(0)] * 11), "Chat(0)=169*I_r")
    ledger = [sum(chat[k][j] for k in range(1, P13)) for j in range(12)]
    require(ledger == [-service] + [F(0)] * 11,
            "signed ledger sum_{k!=0} J(k) = -I_r (eq (10), J=Chat/169)")

    # Energy floor (eq (11)) through exact Parseval in Q.
    chat_energy_all = sum(x * x for x in correlation) / P13
    j_energy = (chat_energy_all - (service / P13) ** 2) / (P13 * P13)
    require(j_energy >= collision * collision / 12, "energy floor (eq (11))")

    # Arrival-refined atom matrix (eq (29)-(31)): rich support.
    arrival = [[F(0)] * 7 for _ in range(7)]
    for y in range(strata):
        a_m = [F(0)] * 7
        e_m = [F(0)] * 7
        for root in range(P13):
            for sheet in range(sheets):
                if fa[y][root][sheet]:
                    a_m[f_atom[y][root][sheet]] += F(1, sheets)
                if ea[y][root][sheet]:
                    e_m[e_atom[y][root][sheet]] += F(1, sheets)
        for w in range(7):
            for z in range(7):
                arrival[w][z] += a_m[w] * e_m[z] / strata
    require(sum(sum(row) for row in arrival) == service, "arrival matrix total")
    arrival_support = sum(1 for row in arrival for x in row if x)

    # Source-refined owner loop (eq (34)-(37)): e*P_O = e forces every
    # occupied source sheet onto the unique owner atom O, and the
    # source-refined current inherits O; the matrix is one loop.
    OWNER = 0
    source = [[F(0)] * 7 for _ in range(7)]
    for y in range(strata):
        a_mass = sum(sum(F(x, sheets) for x in fa[y][r]) for r in range(P13))
        e_mass = sum(sum(F(x, sheets) for x in ea[y][r]) for r in range(P13))
        source[OWNER][OWNER] += a_mass * e_mass / strata
    require(source[OWNER][OWNER] == service and
            all(source[w][z] == 0 for w in range(7) for z in range(7)
                if (w, z) != (OWNER, OWNER)),
            "source-refined owner loop M_{O,O}=169*I_r (eq (37))")

    # Deep-root sidecar (eq (44)-(47)) for the realized stalk: d=13 with
    # strict-profile deep speed C=2*13^5: delta=13, d_0=1, h=C/d=2*13^4,
    # 13|h => the deep phase is base-only; a mod d_0 is trivially retained.
    d_norm = 13
    deep_speed = 2 * 13 ** 5
    delta = 13   # gcd(deep_speed, d_norm)
    require(deep_speed % d_norm == 0 and delta == d_norm, "gcd(C,d)")
    d0 = d_norm // delta
    require(d0 == 1 and (deep_speed // d_norm) % 13 == 0,
            "d_0=1 and 13|h: base-only deep descent (eq (44)-(46))")

    return {
        "strata": strata, "sheets": sheets,
        "corr_stratum": corr_stratum, "correlation": correlation,
        "service": service, "collision": collision,
        "chat": chat, "j_energy": j_energy,
        "arrival_support": arrival_support, "owner": OWNER,
        "owner_loop": source[OWNER][OWNER],
        "d_norm": d_norm, "clock": 1, "d0": d0,
    }


# ----------------------------------------------- PART 2: the THM-2508 cut bank

def canonical_defect():
    """THM-2506 eq (23)-(25) canonical two-row defect of the explicit
    guard/AP/blocker cover (centres 84,71,33,47,58; blocker source 5)."""
    d = [[0] * P7 for _ in range(P13)]
    d[0][5] = 1
    d[0][3] = -1
    d[1][5] = 1
    d[1][4] = -1
    require(all(sum(row) == 0 for row in d), "row-zero law (THM-2508 eq (2))")
    require(sum(abs(x) for row in d for x in row) == 4, "L1=4 <= 18")
    return d


def radon(defect, tau, a, c):
    """THM-2508 eq (3): R_{tau,a,c}(v) = sum_r d(v - tau*rep(ar+c), r)."""
    return [sum(defect[(v - tau * ((a * r + c) % P7)) % P13][r]
                for r in range(P7)) for v in range(P13)]


def psi_dict(defect, tau, a, alpha, beta):
    """THM-2508 eq (13): Psi_{tau,a}(alpha,beta) as a zeta_91 dict."""
    raw = {}
    for c in range(P7):
        xexp = (13 * ((-beta * c) % P7)) % N91
        row = radon(defect, tau, a, c)
        for v in range(P13):
            if row[v]:
                e = (xexp + 7 * ((-alpha * v) % P13)) % N91
                raw[e] = raw.get(e, 0) + row[v]
    return raw


def dtilde_dict(defect, alpha, gamma):
    """THM-2506 eq (4) mixed transform as a zeta_91 dict."""
    raw = {}
    for h in range(P13):
        for r in range(P7):
            if defect[h][r]:
                e = (7 * ((-alpha * h) % P13) + 13 * ((-gamma * r) % P7)) % N91
                raw[e] = raw.get(e, 0) + defect[h][r]
    return raw


def kernel_dict(u, beta):
    """THM-2508 eq (15): K_{u,beta} = sum_{s=0}^{6} zeta^{-us} xi^{-beta s}."""
    raw = {}
    for s in range(P7):
        e = (7 * ((-u * s) % P13) + 13 * ((-beta * s) % P7)) % N91
        raw[e] = raw.get(e, 0) + 1
    return raw


def translate_defect(defect, A, H, B, C):
    """THM-2508 eq (6) pullback d^g(h,r) = d(A^{-1}(h-H), B^{-1}(r-C))."""
    Ainv = pow(A, -1, P13)
    Binv = pow(B, -1, P7)
    return [[defect[(Ainv * (h - H)) % P13][(Binv * (r - C)) % P7]
             for r in range(P7)] for h in range(P13)]


def cut_bank_controls(defect):
    # Positive control: all 5,184 primitive cut coefficients nonzero and
    # factorized as K_{alpha tau,beta} * dtilde(alpha,-beta a) (eq (14)).
    kd_cache = {}
    checked = 0
    for alpha in range(1, P13):
        for beta in range(1, P7):
            for tau in range(1, P13):
                for a in range(1, P7):
                    key = (alpha, beta, tau, a)
                    product = dmul(kernel_dict((alpha * tau) % P13, beta),
                                   dtilde_dict(defect, alpha, (-beta * a) % P7))
                    kd_cache[key] = product
                    direct = dict_to_reduced(psi_dict(defect, tau, a, alpha, beta))
                    expected = dict_to_reduced(product)
                    require(direct == expected, f"factorization failed {key}")
                    require(direct != ZERO72, f"primitive cut coefficient zero {key}")
                    checked += 1
    require(checked == 5184, "5,184-mode census")

    # Hostile control 1: beta=0 recreates the pure F_13 pushforward no-go
    # (eq (20)): all 864 zero-cut-character coefficients vanish.
    zero_checks = 0
    for alpha in range(1, P13):
        for tau in range(1, P13):
            for a in range(1, P7):
                require(dict_to_reduced(psi_dict(defect, tau, a, alpha, 0)) == ZERO72,
                        "zero cut character survived")
                zero_checks += 1
    require(zero_checks == 864, "zero-cut census")

    # Hostile control 2: the gauge-invariant linear readout (eq (21)):
    # sum_v R_{tau,a,c}(v) = 0 on every one of the 504 components.
    invariant_checks = 0
    for tau in range(1, P13):
        for a in range(1, P7):
            for c in range(P7):
                require(sum(radon(defect, tau, a, c)) == 0,
                        "invariant linear readout survived")
                invariant_checks += 1
    require(invariant_checks == 504, "component census")

    # Trivial-character exclusion: K_{0,beta}=0 for beta!=0, so the alpha=0
    # Fourier mode never enters a beta!=0 pairing.  This is a statement about
    # a character index, not the separate physical coordinate value C(0).
    for beta in range(1, P7):
        require(dict_to_reduced(kernel_dict(0, beta)) == ZERO72,
                "K_{0,beta} must vanish for beta!=0")
        for tau in range(1, P13):
            for a in range(1, P7):
                require(dict_to_reduced(psi_dict(defect, tau, a, 0, beta)) == ZERO72,
                        "alpha=0 cut coefficient must vanish")
    return kd_cache, checked, zero_checks, invariant_checks


# --------------------------------------------------- PART 3: typed-sheet audit

def typed_sheet_audit(stalk):
    d_norm, clock = stalk["d_norm"], stalk["clock"]
    R_clock = 13 ** clock
    require(d_norm % P7 != 0 and R_clock % P7 != 0,
            "stalk sheets a in Z/d, b in Z/13^K are 13-primary: no F_7 sheet")
    orbit = sorted({pow(13, L, P7) for L in range(1, 20)})
    require(orbit == [1, 6], "graft multiplier 13^L mod 7 reaches only {+1,-1}")
    require(set(range(1, P7)) - set(orbit) == {2, 3, 4, 5},
            "the full F_7^* cut-torsor action is NOT transported by the graft")
    return orbit


# ----------------------------------------------------- PART 4: the pairing

def chat_scaled_dict(stalk, colour):
    """100*Chat(colour) as an integer zeta_91 dict (stalk denominators 100)."""
    raw = {}
    for s in range(P13):
        value = stalk["correlation"][s] * 100
        require(value.denominator == 1, "correlation scaling")
        if value:
            e = (7 * ((-colour * s) % P13)) % N91
            raw[e] = raw.get(e, 0) + int(value)
    return raw


def pairing_raw(stalk, defect, tau, a, beta, twist, kd_cache=None):
    """100 * sum_{alpha!=0} Chat(twist*alpha) Psi_{tau,a}(alpha,beta),
    as an unreduced zeta_91 dict.  twist=-1 is the dual (conjugate) leg."""
    total = {}
    for alpha in range(1, P13):
        if kd_cache is not None and beta != 0:
            psi = kd_cache[(alpha, beta, tau, a)]
        else:
            psi = psi_dict(defect, tau, a, alpha, beta)
        term = dmul(chat_scaled_dict(stalk, (twist * alpha) % P13), psi)
        for e, v in term.items():
            total[e] = total.get(e, 0) + v
    return total


def pairing_census(stalk, defect, kd_cache):
    results = {}
    zero_pairs = []
    for twist, label in ((-1, "dual"), (1, "untwisted")):
        nonzero = 0
        per_beta = {beta: 0 for beta in range(1, P7)}
        witness = None
        for tau in range(1, P13):
            for a in range(1, P7):
                for beta in range(1, P7):
                    reduced = dict_to_reduced(
                        pairing_raw(stalk, defect, tau, a, beta, twist, kd_cache))
                    if reduced != ZERO72:
                        nonzero += 1
                        per_beta[beta] += 1
                        if witness is None:
                            witness = (tau, a, beta, reduced)
                    else:
                        zero_pairs.append((label, tau, a, beta))
        results[label] = (nonzero, per_beta, witness)

    # beta=0 hostile: the paired current with the cut label summed away
    # vanishes identically (each Psi(alpha,0)=0 by eq (20)).
    for tau in (1, 5, 12):
        for a in (1, 3, 6):
            require(dict_to_reduced(
                pairing_raw(stalk, defect, tau, a, 0, -1)) == ZERO72,
                "beta=0 pairing must vanish")

    # alpha=0 exclusion at the paired level: adding the trivial-character
    # term Chat(0)*Psi(0,beta) changes nothing.  Do not identify this index
    # with the physical coordinate condition C(0)=0.
    for beta in range(1, P7):
        extra = dmul(chat_scaled_dict(stalk, 0), psi_dict(defect, 1, 1, 0, beta))
        require(dict_to_reduced(extra) == ZERO72,
                "trivial Fourier character must contribute exactly zero")

    # Consistency: factorized pairing equals direct-Psi pairing at witness.
    direct = {}
    for alpha in range(1, P13):
        term = dmul(chat_scaled_dict(stalk, (-alpha) % P13),
                    psi_dict(defect, 1, 1, alpha, 1))
        for e, v in term.items():
            direct[e] = direct.get(e, 0) + v
    require(dict_to_reduced(direct) ==
            dict_to_reduced(pairing_raw(stalk, defect, 1, 1, 1, -1, kd_cache)),
            "direct/factorized pairing consistency")

    # Diagonal factors are individually nonzero in this external alignment.
    # The delta_2 hostile below shows that their sum can nevertheless cancel;
    # THM-2478 (21) does not type the cut character as a collision character.
    for alpha in range(1, P13):
        term = dict_to_reduced(dmul(chat_scaled_dict(stalk, (-alpha) % P13),
                                    kd_cache[(alpha, 1, 1, 1)]))
        require(term != ZERO72, "diagonal colour component vanished")
    return results, zero_pairs


def delta2_hostile(defect):
    """Boolean full-support factors whose external alpha-sum cancels.

    Take one predecessor packet at root 0 and one disjoint arrival packet at
    root 2.  Its correlation is C=delta_2, so every nontrivial root Fourier
    colour is a monomial and hence nonzero.  For the dual convention,

      sum_{alpha != 0} zeta_13^(2 alpha) Psi(alpha,beta)
        = 13 sum_c R_{tau,a,c}(2) xi_7^(-beta c),

    because every R component has total sum zero.  At tau=4 every one of the
    42 values R_{4,a,c}(2) is zero, so all 36 (a,beta) placements vanish.
    """
    arrival = [F(0)] * P13
    predecessor = [F(0)] * P13
    arrival[2] = F(1)
    predecessor[0] = F(1)
    require(all(arrival[r] * predecessor[r] == 0 for r in range(P13)),
            "delta_2 packets must be same-root disjoint")
    correlation = [sum(arrival[(r + shift) % P13] * predecessor[r]
                       for r in range(P13)) for shift in range(P13)]
    require(correlation == [F(int(s == 2)) for s in range(P13)],
            "Boolean packets must give C=delta_2")
    require(correlation[0] == 0 and sum(correlation) > 0,
            "delta_2 must have C(0)=0 and positive service")
    for alpha in range(1, P13):
        require(any(phi13_reduce(correlation, alpha)),
                f"delta_2 root colour vanished at alpha={alpha}")

    hostile_stalk = {"correlation": correlation}
    zeros = {}
    zero_sets = {}
    for twist, label in ((-1, "dual"), (1, "untwisted")):
        placements = []
        for tau in range(1, P13):
            for a in range(1, P7):
                for beta in range(1, P7):
                    value = dict_to_reduced(
                        pairing_raw(hostile_stalk, defect, tau, a, beta, twist))
                    if value == ZERO72:
                        placements.append((tau, a, beta))
        zeros[label] = len(placements)
        zero_sets[label] = placements
        require(len(placements) == 108,
                f"delta_2 {label} exact zero census changed")

    zero_tau = {}
    for label, placements in zero_sets.items():
        zero_tau[label] = {
            tau: sum(t == tau for t, _a, _beta in placements)
            for tau in range(1, P13)
            if any(t == tau for t, _a, _beta in placements)
        }
    require(zero_tau["dual"] == {4: 36, 6: 36, 12: 36},
            "delta_2 dual zero-tau classification changed")
    require(zero_tau["untwisted"] == {1: 36, 3: 36, 7: 36},
            "delta_2 untwisted zero-tau classification changed")

    tau4 = [(tau, a, beta) for tau, a, beta in zero_sets["dual"] if tau == 4]
    require(len(tau4) == 36, "delta_2 dual tau=4 must vanish for all a,beta")
    radon_teeth = 0
    for a in range(1, P7):
        for c in range(P7):
            require(radon(defect, 4, a, c)[2] == 0,
                    f"delta_2 cancellation mechanism failed at a={a},c={c}")
            radon_teeth += 1
    require(radon_teeth == 42, "delta_2 Radon-tooth census")
    return zeros, zero_tau, tau4[0], radon_teeth


def prime_factors(number):
    factors = []
    candidate = 2
    while candidate * candidate <= number:
        if number % candidate == 0:
            factors.append(candidate)
            while number % candidate == 0:
                number //= candidate
        candidate += 1
    if number > 1:
        factors.append(number)
    return factors


def is_prime(number):
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def primitive_root(prime):
    factors = prime_factors(prime - 1)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // q, prime) != 1 for q in factors):
            return candidate
    raise RuntimeError(f"no primitive root modulo {prime}")


def independent_five_prime_audit(stalk, defect):
    """Independent direct Radon/DFT evaluation at five p=1 (mod 91).

    This path does not use Phi_91 reduction, psi_dict, dtilde_dict,
    kernel_dict, pairing_raw, or kd_cache.  It is only an independent check;
    exact nonvanishing is supplied by the Q(zeta_91) computation above.
    """
    primes = (547, 911, 1093, 2003, 2549)
    direct_radon = {}
    for tau in range(1, P13):
        for a in range(1, P7):
            for c in range(P7):
                direct_radon[(tau, a, c)] = tuple(
                    sum(defect[(v - tau * ((a * r + c) % P7)) % P13][r]
                        for r in range(P7))
                    for v in range(P13))

    total_nonzero = 0
    total_beta0_zero = 0
    per_prime = {}
    for prime in primes:
        require(is_prime(prime) and (prime - 1) % N91 == 0,
                f"invalid split prime {prime}")
        generator = primitive_root(prime)
        root91 = pow(generator, (prime - 1) // N91, prime)
        require(pow(root91, N91, prime) == 1 and
                pow(root91, N91 // P7, prime) != 1 and
                pow(root91, N91 // P13, prime) != 1,
                f"root does not have exact order 91 modulo {prime}")
        root13 = pow(root91, P7, prime)
        root7 = pow(root91, P13, prime)

        correlation100 = []
        for value in stalk["correlation"]:
            scaled = value * 100
            require(scaled.denominator == 1, "five-prime correlation scaling")
            correlation100.append(int(scaled) % prime)
        chat = {}
        for colour in range(P13):
            chat[colour] = sum(
                correlation100[s] * pow(root13, (-colour * s) % P13, prime)
                for s in range(P13)) % prime

        psi = {}
        for tau in range(1, P13):
            for a in range(1, P7):
                for alpha in range(1, P13):
                    for beta in range(P7):
                        value = 0
                        for c in range(P7):
                            cut_phase = pow(root7, (-beta * c) % P7, prime)
                            row = direct_radon[(tau, a, c)]
                            value += cut_phase * sum(
                                row[v] * pow(root13, (-alpha * v) % P13, prime)
                                for v in range(P13))
                        psi[(tau, a, alpha, beta)] = value % prime

        prime_nonzero = 0
        prime_beta0_zero = 0
        for twist in (-1, 1):
            for tau in range(1, P13):
                for a in range(1, P7):
                    for beta in range(1, P7):
                        value = sum(
                            chat[(twist * alpha) % P13] *
                            psi[(tau, a, alpha, beta)]
                            for alpha in range(1, P13)) % prime
                        require(value != 0,
                                f"five-prime external pairing zero mod {prime} at "
                                f"{twist,tau,a,beta}")
                        prime_nonzero += 1
                    value0 = sum(
                        chat[(twist * alpha) % P13] * psi[(tau, a, alpha, 0)]
                        for alpha in range(1, P13)) % prime
                    require(value0 == 0,
                            f"five-prime beta=0 hostile failed mod {prime}")
                    prime_beta0_zero += 1
        require(prime_nonzero == 864 and prime_beta0_zero == 144,
                f"five-prime census failed mod {prime}")
        per_prime[prime] = (prime_nonzero, prime_beta0_zero)
        total_nonzero += prime_nonzero
        total_beta0_zero += prime_beta0_zero
    require(total_nonzero == 4320 and total_beta0_zero == 720,
            "five-prime aggregate census")
    return primes, per_prime, total_nonzero, total_beta0_zero


# --------------------------------------------- PART 5: forward graft, exact grid

DRIFT = [2, -1, 3, 0, -2]   # old signed current; integral q = 2/5 != 0 (C4)


def graft_colour_grid(stalk):
    """THM-2478 eq (17)/(21): delayed collision colours g_k(d*13^L y)
    against an old current b(y), exact midpoint integrals on a common
    refinement grid; all twelve colours nonzero at the first lawful
    delays and the BV error decreases."""
    strata = stalk["strata"]
    corr_stratum = stalk["corr_stratum"]
    corr_total = [sum(corr_stratum[y][s] for y in range(strata))
                  for s in range(P13)]
    q_old = sum(map(F, DRIFT)) / len(DRIFT)
    require(q_old != 0, "old current must be charged")

    limit = {k: tuple(q_old / strata * x for x in phi13_reduce(corr_total, k))
             for k in range(1, P13)}
    errors = {}
    nonzero_colours = {}
    for L in (2, 3):   # lawful delays: L >= K+1 with K=1
        m = stalk["d_norm"] * 13 ** L
        cells = 20 * m
        mass = [0] * strata
        for j in range(cells):
            two_j = 2 * j + 1
            mass[(strata * m * two_j) // (2 * cells) % strata] += \
                DRIFT[(len(DRIFT) * two_j) // (2 * cells) % len(DRIFT)]
        mass_frac = [F(x, cells) for x in mass]
        coupled = [sum(mass_frac[y] * corr_stratum[y][s] for y in range(strata))
                   for s in range(P13)]
        err_total = F(0)
        count = 0
        for k in range(1, P13):
            vec = phi13_reduce(coupled, k)
            if any(vec):
                count += 1
            err_total += sum(abs(a - b) for a, b in zip(vec, limit[k]))
        nonzero_colours[L] = count
        errors[L] = err_total
    require(nonzero_colours[2] == 12 and nonzero_colours[3] == 12,
            "all twelve colours must survive at every lawful delay")
    require(errors[3] < errors[2], "BV error must decrease (eq (20))")
    return errors


def graft_neutrality(stalk):
    """THM-2478 eq (19a): at L=K+1 every stalk leg multiplier (13^L for the
    current/word and bare-source legs, 13^{L-K} for the ancestry leg) makes
    every old /13 target co-shift integral; at L=K the ancestry leg fails."""
    K = stalk["clock"]
    L = K + 1
    for multiplier in (13 ** L, 13 ** (L - K), 13 ** L):
        for theta in range(P13):
            require((multiplier * theta) % P13 == 0,
                    "full-stalk target neutrality at L=K+1")
    require(any((13 ** 0 * theta) % P13 != 0 for theta in range(1, P13)),
            "L=K ancestry-leg neutrality must FAIL (sharp boundary)")


def graft_owner_loop(stalk):
    """Owner-type retention: the grafted service is still one owner loop.
    Gate G(z)=E(z)Q(T^K z) delayed by N_L; the source-refined atom algebra
    is untouched (e*P_O=e persists), so the grafted matrix has the single
    entry (O,O) with positive mass."""
    owner_cells = [1, 1, 0, 0, 0]
    word_cells = [0, 1, 1, 0, 1]
    K, L = stalk["clock"], 2
    strata = stalk["strata"]
    sdens25 = [int(sum(stalk["corr_stratum"][y][s] for s in range(P13)) * 25)
               for y in range(strata)]
    cells = 20 * 13 ** (L - 1 + K)
    total = 0
    for j in range(cells):
        two_j = 2 * j + 1
        scale = 13 ** (L - 1)
        gate = (owner_cells[(5 * scale * two_j) // (2 * cells) % 5]
                * word_cells[(5 * scale * (13 ** K) * two_j) // (2 * cells) % 5])
        if gate:
            total += sdens25[(strata * two_j) // (2 * cells) % strata]
    grafted = F(total, 25 * cells)
    require(grafted > 0, "grafted owner-loop service must stay positive")
    matrix = {(stalk["owner"], stalk["owner"]): grafted}
    require(len(matrix) == 1, "grafted matrix is one owner loop")
    return grafted


# --------------------------------------------------- PART 6: kappa/carry audit

def kappa_audit(stalk, defect, kd_cache):
    tau, a, alpha, beta = 1, 1, 1, 1

    # THM-2508 eq (19) under all pure CRT translations (A=1,B=1):
    # Psi^{d^g} = zeta^{-alpha H} xi^{beta a C} Psi^d.
    base_psi = psi_dict(defect, tau, a, alpha, beta)
    for H in range(P13):
        for C in range(P7):
            moved = translate_defect(defect, 1, H, 1, C)
            left = dict_to_reduced(psi_dict(moved, tau, a, alpha, beta))
            phase = dmul(z13(-alpha * H), x7(beta * a * C))
            require(left == dict_to_reduced(dmul(phase, base_psi)),
                    f"translation covariance failed H={H} C={C}")

    # One full gauge g=(A,H,B,C)=(2,3,2,4):
    # Psi^{d^g}_{tau,a}(alpha,beta)
    #   = zeta^{-alpha H} xi^{beta a C} Psi^d_{A^{-1}tau, aB}(A alpha, beta).
    A, H, B, C = 2, 3, 2, 4
    moved = translate_defect(defect, A, H, B, C)
    left = dict_to_reduced(psi_dict(moved, tau, a, alpha, beta))
    right = dmul(dmul(z13(-alpha * H), x7(beta * a * C)),
                 psi_dict(defect, (pow(A, -1, P13) * tau) % P13,
                          (a * B) % P7, (A * alpha) % P13, beta))
    require(left == dict_to_reduced(right), "full gauge covariance (19) failed")

    # The pairing is NOT invariant: some H changes it; the exact
    # kappa-dressing identity P_H = sum_alpha zeta^{-alpha H} term_alpha holds.
    base_pairing = dict_to_reduced(pairing_raw(stalk, defect, tau, a, beta, -1, kd_cache))
    term = {}
    for al in range(1, P13):
        term[al] = dmul(chat_scaled_dict(stalk, (-al) % P13),
                        kd_cache[(al, beta, tau, a)])
    changed = 0
    collapse_H = {}
    for H in range(P13):
        moved = translate_defect(defect, 1, H, 1, 0)
        pH_raw = pairing_raw(stalk, moved, tau, a, beta, -1)
        pH = dict_to_reduced(pH_raw)
        dressed = {}
        for al in range(1, P13):
            piece = dmul(z13(-al * H), term[al])
            for e, v in piece.items():
                dressed[e] = dressed.get(e, 0) + v
        require(pH == dict_to_reduced(dressed), f"kappa dressing failed H={H}")
        if pH != base_pairing:
            changed += 1
        for e, v in pH_raw.items():
            collapse_H[e] = collapse_H.get(e, 0) + v
    require(changed > 0, "pairing must move under the kappa carry")
    require(dict_to_reduced(collapse_H) == ZERO72,
            "H-torsor (kappa) collapse of the pairing must vanish")

    # Septimal side: pure C-translation multiplies the pairing by the global
    # phase xi^{beta a C}; averaging over the C-torsor kills it.
    collapse_C = {}
    for C in range(P7):
        moved = translate_defect(defect, 1, 0, 1, C)
        pC_raw = pairing_raw(stalk, moved, tau, a, beta, -1)
        phase = x7(beta * a * C)
        expected = {}
        base_raw = pairing_raw(stalk, defect, tau, a, beta, -1, kd_cache)
        for e, v in dmul(phase, base_raw).items():
            expected[e] = expected.get(e, 0) + v
        require(dict_to_reduced(pC_raw) == dict_to_reduced(expected),
                f"cut-phase covariance failed C={C}")
        for e, v in pC_raw.items():
            collapse_C[e] = collapse_C.get(e, 0) + v
    require(dict_to_reduced(collapse_C) == ZERO72,
            "C-torsor (cut) collapse of the pairing must vanish")
    return changed


# ------------------------------------------------ PART 7: rebase boundary replay

def rebase_boundary():
    """THM-2478 eq (23)-(26): forward graft keeps the deep probe sheet-free;
    rebasing at L > lambda makes a mod 13^{L-lambda} essential (root-uniform
    danger/safe hostile at the strict 1/14 radius)."""
    deep = 2 * 13 ** 5
    z = F(4, 31)
    shallow_delay = 4
    phases = {(F(deep * (z + a), 13 ** shallow_delay) % 1)
              for a in range(13 ** shallow_delay)}
    require(len(phases) == 1, "deep phase descends before its valuation")

    lost_delay, valuation = 6, 5
    lost_modulus = 13 ** (lost_delay - valuation)
    unit_inverse = pow(deep // 13 ** valuation, -1, lost_modulus)
    for root in range(P13):
        a0 = (root * unit_inverse) % lost_modulus
        a1 = ((root + 1) * unit_inverse) % lost_modulus
        phase0 = (F(deep * a0, 13 ** lost_delay) - F(root, P13)) % 1
        phase1 = (F(deep * a1, 13 ** lost_delay) - F(root, P13)) % 1
        require(phase0 == 0 and phase1 == F(1, P13),
                "essential ancestry residue (eq (25)-(26))")
        require(phase0 < F(1, 14) and phase1 > F(1, 14),
                "danger/safe sheet hostile at radius 1/14")
    return lost_modulus


# --------------------------------------------------------------------- driver

def fmt72(vec, limit=8):
    entries = [(i, x) for i, x in enumerate(vec) if x]
    head = ",".join(f"z91^{i}:{x}" for i, x in entries[:limit])
    tail = ",..." if len(entries) > limit else ""
    return f"{{{head}{tail}}}; nonzero_coeffs={len(entries)}/72"


def main():
    print("LRC14 CUT-BUNDLE x ROOT-CORRELATION PAIRING -- HOSTILE-AUDITED EXACT CONTROL")
    print("status: FINITE-EXACT EXTERNAL CONTROL; not a physical/common-ancestry current")
    print("objects: synthetic weighted control from THM-2471's companion EXTERNALLY"
          " multiplied by THM-2508's cut bundle on THM-2506's high-branch defect")

    stalk = build_stalk()
    print("[P1] synthetic root control: service=169*I_r with I_r="
          f"{stalk['collision']}; service={stalk['service']}; C(0)=0; "
          f"colours_nonzero=12/12; signed_ledger=sum_(k!=0)J(k)=-I_r=PASS")
    print(f"[P1] energy_floor J>=I^2/12: {stalk['j_energy']}>="
          f"{stalk['collision'] ** 2 / 12}: PASS; "
          f"arrival_matrix_support={stalk['arrival_support']}; "
          f"separate source-margin owner loop M_OO={stalk['owner_loop']}=169*I_r")
    print(f"[P1] separate deep margin: d={stalk['d_norm']}, C_deep=2*13^5, d_0="
          f"{stalk['d0']}; neither owner nor deep labels occur in the P4 product")

    defect = canonical_defect()
    kd_cache, modes, zeros, invariants = cut_bank_controls(defect)
    print(f"[P2] cut bank on THM-2506 canonical defect: primitive_modes={modes}"
          " nonzero+factorized (eq (14)-(16)): PASS  [positive control]")
    print(f"[P2] hostile controls: beta=0 vanishings={zeros}/864 PASS; "
          f"invariant linear readouts sum_v R=0: {invariants}/504 PASS; "
          "K_(0,beta)=0 trivial-Fourier-character exclusion: PASS")

    orbit = typed_sheet_audit(stalk)
    print("[P3] typed sheets: stalk root u (temporal collision fibre) and cut-"
          "bundle h (static quotient-stalk row, THM-2506 (28)) are DIFFERENT"
          " typed sheets (THM-2506 Sec.6; THM-2508 Sec.6 scope)")
    print("[P3] NO LICENSED CUT/COLLISION BRIDGE: THM-2478 (21a)-(21b) aligns a"
          " fixed old target-component label with a future collision label only;"
          " it does not type static cut alpha as temporal collision k.  The P4"
          " alpha-alignment is external.  No septimal identification exists:"
          f" stalk sheet moduli (13,13) are 13-primary, graft orbit 13^L mod 7={orbit}"
          " (parity subgroup, not F_7^*); cut phase rides the static leg, parity-pinned")

    results, zero_pairs = pairing_census(stalk, defect, kd_cache)
    for label in ("dual", "untwisted"):
        nonzero, per_beta, witness = results[label]
        beta_line = ",".join(f"beta={b}:{per_beta[b]}/72" for b in range(1, P7))
        print(f"[P4] pairing[{label}]: nonzero {nonzero}/432 over"
              f" (tau,a,beta) in F_13^* x F_7^* x F_7^*; per cut character: {beta_line}")
        if witness is not None:
            tau, a, beta, vec = witness
            print(f"[P4] pairing[{label}] witness (tau,a,beta)=({tau},{a},{beta}):"
                  f" 100*P={fmt72(vec)}")
    if zero_pairs:
        print(f"[P4] EXACT ZERO PLACEMENTS: {zero_pairs}")
    else:
        print("[P4] displayed external control: no beta!=0 placement vanishes;"
              " beta=0 pairing=0 PASS; trivial alpha=0 character contributes 0"
              " (this is NOT the physical coordinate statement C(0)=0);"
              " all twelve diagonal colour components nonzero PASS")

    hostile_zeros, hostile_tau, hostile_first, hostile_teeth = delta2_hostile(defect)
    print("[P4-hostile] Boolean C=delta_2: same-root disjoint, positive service,"
          " C(0)=0, root colours nonzero 12/12; canonical defect primitive modes"
          " nonzero 5184/5184")
    print(f"[P4-hostile] nevertheless external zeros dual={hostile_zeros['dual']}/432,"
          f" untwisted={hostile_zeros['untwisted']}/432; dual tau=4 gives 36/36"
          f" zeros, first={hostile_first}; mechanism R_(4,a,c)(2)=0 for"
          f" {hostile_teeth}/42 teeth.  Full factor support != nonzero alpha-sum.")
    print(f"[P4-hostile] exact zero-tau classes: dual={hostile_tau['dual']};"
          f" untwisted={hostile_tau['untwisted']}")

    primes, prime_census, prime_nonzero, prime_beta0 = \
        independent_five_prime_audit(stalk, defect)
    prime_line = ",".join(
        f"p={p}:{prime_census[p][0]}/864" for p in primes)
    print(f"[P4-independent] direct finite-field Radon/DFT path ({prime_line});"
          f" aggregate beta!=0 nonzero={prime_nonzero}/4320;"
          f" beta=0 zeros={prime_beta0}/720 PASS")

    errors = graft_colour_grid(stalk)
    graft_neutrality(stalk)
    grafted = graft_owner_loop(stalk)
    print(f"[P5] SEPARATE synthetic grid graft (K=1): twelve colours nonzero"
          " at L=2 and L=3;"
          f" BV error L=2:{errors[2]} > L=3:{errors[3]} (decreasing) PASS")
    print("[P5] full-stalk target neutrality at L=K+1 (multipliers 13^L,13^(L-K))"
          " PASS; sharp failure at L=K PASS;"
          f" separate grafted owner margin M_OO={grafted}>0 PASS; these checks"
          " do not license or enrich the P4 external product")

    changed = kappa_audit(stalk, defect, kd_cache)
    print(f"[P6] kappa/carry audit: covariance (19) under all 91 CRT translations"
          " + one full gauge PASS; pairing moves under kappa"
          f" ({changed}/12 nonzero carries change it); exact kappa-dressing identity PASS")
    print("[P6] invariant collapses: sum over H-torsor (kappa) = 0 PASS;"
          " sum over C-torsor (cut) = 0 PASS -- the invariant readout is dead"
          " (THM-2508 (20)-(21) at the paired level); only the external"
          " equivariant pairing survives; no physical sidecar-retention"
          " conclusion follows")

    lost_modulus = rebase_boundary()
    print("[P7] rebase boundary: forward graft keeps the old deep probe at x"
          " (sheet-free; d_0=1); rebasing at L=6>lambda=5 breaks the ancestry"
          f" residue a mod {lost_modulus} (root-uniform danger/safe hostile) --"
          " this is a separate graft control, not ancestry for P4")

    print("conventions: C1 synthetic weighted control (not a realized row);"
          " C2 ledger-anchored J=Chat/169; C3 untyped external alpha alignment,"
          " dual and untwisted variants; C4 separate scalar grid graft")

    dual_nonzero = results["dual"][0]
    untw_nonzero = results["untwisted"][0]
    if dual_nonzero == 432 and untw_nonzero == 432 and not zero_pairs:
        print("VERDICT: FINITE-EXACT EXTERNAL CONTROL.  The displayed unrelated"
              " coefficient arrays have dual and untwisted pairing nonzero"
              " 432/432 in Q(zeta_91), with five-prime independent agreement.")
        print("NON-UNIVERSALITY: C=delta_2 has full nontrivial root-colour support"
              " and the defect has full primitive support, but 108/432 pairings"
              " vanish in each convention.  Multiplying unrelated nonzero"
              " coefficients and then summing does not produce a current.")
        print("NO CONSEQUENCE: no common ancestry cospan, physical current, target"
              " charge, owner retention, typed cut/collision bridge, live-row"
              " transplant, or LRC(14) branch exclusion has been proved.")
        print("BRANCH-TRANSPLANT QUANTIFIER: for every live low-septimal row R in"
              " the 165-row bank, construct one positive-measure Boolean base"
              " Omega_R carrying collision, owner, word, deep and cut data from"
              " the SAME atoms, plus a lawful typed intertwiner; then choose a"
              " row-dependent but omega-independent (tau_R,a_R,beta_R) and prove"
              " integral_{Omega_R} sum_{alpha!=0} Gamma_R(omega,-alpha)"
              " Psi_{d_R(omega),tau_R,a_R}(alpha,beta_R) dmu != 0.  The product"
              " must be formed before marginalization/DFT/integration; pointwise"
              " factor support is insufficient because alpha and omega can cancel.")
    else:
        print("VERDICT: NEGATIVE/PARTIAL -- exact zero placements listed above;"
              " see zero set for the breaking sidecar.")
    print(f"ALL EXACT CHECKS PASSED; exact_require_checks={CHECKS}")


if __name__ == "__main__":
    main()

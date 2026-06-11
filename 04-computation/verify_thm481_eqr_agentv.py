#!/usr/bin/env python3
"""
verify_thm481_eqr_agentv.py -- ADVERSARIAL audit of the THM-481 eQR proof sketch.

PROOF SKETCH under attack (p prime, p = 7 mod 8):
  H = order-(p+1) skew-Hadamard bordered from the Paley tournament:
      H[inf,inf]=+1, H[inf,j]=+1, H[i,inf]=-1, H[i,i]=+1, H[i,j]=chi(j-i).
  (a) row inf of (J-H)/2 is zero; row i has support (i+QNR) u {inf}.
  (b) |QNR|=(p-1)/2 odd (p=3 mod 4)  =>  inf-bit of each generator = parity of
      its core  =>  C(H) := rowspan_{GF(2)}(J-H)/2 = parity-extension of the
      cyclic code <theta>, theta(x) = sum_{s in QNR} x^s.
  (c) p = -1 mod 8 => theta is the idempotent of one of the two QR codes of
      dimension (p+1)/2  =>  C(H) = extended QR code of length p+1, Type II.

AUDIT POINTS (per the task):
  (1) sign convention: chi(j-i) = -1 iff j-i in QNR; repo Paley convention is
      i->j iff chi(j-i)=+1, so the core of (J-H)/2 is the IN-neighborhood
      (loss) matrix A^T, with row supports i+QNR.  Checked structurally.
  (2) where '2 is a QR mod p iff p = +-1 mod 8' is needed: Frobenius squaring
      sends theta to sum x^{2s}; theta idempotent  <=>  2*QNR = QNR  <=>
      2 in QR.  Checked, plus NEGATIVE CONTROL p=11 (p=3 mod 8: 2 in QNR,
      theta NOT idempotent, <theta> = full space, C(H) = [12,11] even-weight
      code -- NOT a QR code).  Second control p=17 (p=1 mod 8: chi(-1)=+1, H
      as defined is NOT Hadamard at all, and |QNR| even kills the parity step).
  (3) idempotent claim at p=7 and p=23: theta^2 = theta mod (x^p-1) over
      GF(2); dim <theta> = (p+1)/2; and WHICH code: build GF(2^m), pick a
      primitive p-th root alpha, form q(x)=prod_{r in QR}(x-alpha^r) and
      n(x)=prod_{s in QNR}(x-alpha^s), check both lie in GF(2)[x], and that
      <theta_QNR> = <n(x)> = N (zeros at NONRESIDUE powers of alpha), while
      <theta_QR> = <q(x)> = Q.  (Q,N labels are alpha-dependent; swapping
      alpha -> alpha^s, s in QNR, swaps them.  The invariant statement:
      theta_QNR generates the QR code whose zero set is {alpha^s: s in QNR}
      for the SAME alpha used to define q(x).)
  (4) parity-extension convention: literature 'extended QR code' appends an
      overall-parity coordinate x_inf = sum_i x_i (binary case; MWS Ch.16 S.4
      use x_inf = -gamma*sum x_i, gamma=1 over GF(2)).  Checked: C(H) is
      exactly {(c, wt(c) mod 2) : c in <theta>}, every codeword has even
      total weight, and C(H) restricted-to-core (puncture at inf) = <theta>.

MAIN equality + Type II checks at p=7 ([8,4,4] = e8-hat, wd {0:1,4:14,8:1})
and p=23 ([24,12,8] = extended Golay, wd {0:1,8:759,12:2576,16:759,24:1}).
Also: ext(N) != ext(Q) as subspaces, but j -> s*j (s in QNR, inf fixed) maps
one onto the other (so 'THE extended QR code' is well-defined only up to
equivalence -- the sketch should say WHICH: it is ext(N), the QNR one).

Everything built from scratch; no project imports.
Author: agent verifier, 2026-06-11.
"""

import sys

FAILURES = []


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(f"{name}: {detail}")


# ---------------------------------------------------------------- GF(2)[x]
def pmod(a, b):
    db = b.bit_length() - 1
    while a and a.bit_length() - 1 >= db:
        a ^= b << (a.bit_length() - 1 - db)
    return a


def pmul(a, b, mod=None):
    r = 0
    sh = 0
    while b:
        if b & 1:
            r ^= a << sh
        b >>= 1
        sh += 1
    if mod is not None:
        r = pmod(r, mod)
    return r


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def rref(vectors):
    """Canonical reduced-echelon basis of the GF(2) span (bitmask vectors)."""
    basis = []
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    for i in range(len(basis)):
        for j in range(len(basis)):
            if i != j:
                lead = 1 << (basis[j].bit_length() - 1)
                if basis[i] & lead:
                    basis[i] ^= basis[j]
    return sorted(basis, reverse=True)


def span_words(basis):
    words = [0]
    for b in basis:
        words += [w ^ b for w in words]
    return words


def wdist(words):
    d = {}
    for w in words:
        wt = bin(w).count("1")
        d[wt] = d.get(wt, 0) + 1
    return dict(sorted(d.items()))


def in_span(v, basis):
    for b in basis:
        v = min(v, v ^ b)
    return v == 0


# ------------------------------------------------- irreducible factor finder
def irreducible_factors(poly):
    """Brute-force factorization of a small GF(2)[x] polynomial into
    irreducibles (trial division by increasing degree).  Fine for p <= 23."""
    factors = []
    d = 1
    while poly.bit_length() - 1 >= 2 * d:
        found = True
        while found:
            found = False
            for c in range(1 << (d - 1)):
                cand = (1 << d) | (c << 1) | 1  # monic, constant term 1
                if d == 1:
                    cand = 0b11  # only x+1 has nonzero constant term
                if pmod(poly, cand) == 0:
                    factors.append(cand)
                    # divide out
                    q, r = 0, poly
                    db = cand.bit_length() - 1
                    while r and r.bit_length() - 1 >= db:
                        sh = r.bit_length() - 1 - db
                        q |= 1 << sh
                        r ^= cand << sh
                    assert r == 0
                    poly = q
                    found = True
                    break
        d += 1
    if poly.bit_length() - 1 >= 1:
        factors.append(poly)
    return factors


def poly_str(g):
    return " + ".join(
        f"x^{k}" if k > 1 else ("x" if k == 1 else "1")
        for k in range(g.bit_length() - 1, -1, -1)
        if (g >> k) & 1
    )


# ============================================================ per-prime audit
def audit(p):
    print("\n" + "=" * 72)
    print(f"AUDIT p = {p}   (p mod 8 = {p % 8}, p mod 4 = {p % 4})")
    print("=" * 72)
    INF = p
    NB = p + 1
    XP1 = (1 << p) | 1  # x^p + 1 over GF(2)

    QR = sorted({(i * i) % p for i in range(1, p)})
    QNR = sorted(set(range(1, p)) - set(QR))
    print(f"QR  = {QR}")
    print(f"QNR = {QNR}")

    # ---- audit point (1): chi sign convention vs repo Paley tournament
    def chi(a):
        a %= p
        assert a != 0
        return 1 if a in QR else -1

    check("(1) chi(j-i) = -1  <=>  j-i in QNR (definition of chi)",
          all((chi(d) == -1) == (d in QNR) for d in range(1, p)))
    skew_ok = (p % 4 == 3)
    check(f"chi(-1) = {chi(p-1)}; skewness needs chi(-1)=-1, i.e. p=3 mod 4 "
          f"(here p%4={p%4})", (chi(p - 1) == -1) == skew_ok)

    # build H per the sketch
    H = [[0] * NB for _ in range(NB)]
    H[INF][INF] = 1
    for j in range(p):
        H[INF][j] = 1
        H[j][INF] = -1
    for i in range(p):
        for j in range(p):
            H[i][j] = 1 if i == j else chi(j - i)

    is_skew = all(H[i][j] + H[j][i] == (2 if i == j else 0)
                  for i in range(NB) for j in range(NB))
    is_had = all(sum(H[i][k] * H[j][k] for k in range(NB))
                 == (NB if i == j else 0)
                 for i in range(NB) for j in range(NB))
    if p % 4 == 3:
        check("H + H^T = 2I (skew)", is_skew)
        check(f"H H^T = {NB} I (Hadamard)", is_had)
    else:
        check("CONTROL p=1 mod 4: H as defined is NOT skew and NOT Hadamard",
              (not is_skew) and (not is_had),
              f"skew={is_skew}, hadamard={is_had}")

    # repo tournament convention: i->j iff chi(j-i)=+1  (a038375_solver.c etc)
    # => core of H = I + A - A^T with A = Paley adjacency; core of (J-H)/2 = A^T
    A = [[1 if (i != j and (j - i) % p in QR) else 0 for j in range(p)]
         for i in range(p)]
    core_match = all(H[i][j] == (1 if i == j else A[i][j] - A[j][i])
                     for i in range(p) for j in range(p))
    if p % 4 == 3:
        check("(1) core(H) = I + A - A^T for repo Paley tournament A[i][j]=1 "
              "iff (j-i) in QR", core_match)
    else:
        check("CONTROL p=1 mod 4: A is NOT a tournament (chi(-1)=+1, arcs "
              "doubled) so core(H) != I + A - A^T", not core_match)

    # rows of (J-H)/2 as bitmasks
    def row_mask(row):
        m = 0
        for j in range(NB):
            v = (1 - row[j]) // 2
            assert v in (0, 1)
            if v:
                m |= 1 << j
        return m

    rows = [row_mask(H[i]) for i in range(NB)]
    check("(a) row inf of (J-H)/2 is zero", rows[INF] == 0)
    supp_ok = all(
        rows[i] == (sum(1 << ((i + s) % p) for s in QNR) | (1 << INF))
        for i in range(p))
    check("(a) row i of (J-H)/2 has support (i+QNR) u {inf}", supp_ok)
    print("    note: with repo convention i->j iff chi(j-i)=+1, the core of "
          "(J-H)/2 is A^T (the LOSS matrix), support of row i = in-nbrs of i.")

    # ---- audit point (b): parity of generators
    wt_core = (p - 1) // 2
    odd_core = wt_core % 2 == 1
    check(f"(b) |QNR| = (p-1)/2 = {wt_core} is "
          f"{'ODD' if odd_core else 'EVEN'}; odd <=> p=3 mod 4 (here p%4={p%4})",
          odd_core == (p % 4 == 3))
    gen_parity_ok = all(
        ((rows[i] >> INF) & 1) == (bin(rows[i] & ((1 << p) - 1)).count("1") % 2)
        for i in range(p))
    if p % 4 == 3:
        check("(b) inf-bit of every generator = parity of its core", gen_parity_ok)
    else:
        check("CONTROL p=1 mod 4: inf-bit != core parity (|QNR| even) -- the "
              "parity-extension step of the sketch FAILS here", not gen_parity_ok)
        print("    (stopping this control here: premise 'H skew-Hadamard' "
              "already failed)")
        return

    # C(H)
    CH = rref(rows)
    print(f"dim C(H) = {len(CH)}")

    # ---- audit point (2): theta idempotency and the role of 2 in QR
    theta = 0
    for s in QNR:
        theta |= 1 << s
    theta_R = 0
    for r in QR:
        theta_R |= 1 << r

    two_qr = 2 in QR
    check(f"(2) 2 in QR is {two_qr}; expected (p=+-1 mod 8) = "
          f"{p % 8 in (1, 7)}", two_qr == (p % 8 in (1, 7)))
    doubled = sorted((2 * s) % p for s in QNR)
    check(f"(2) 2*QNR == QNR  <=>  2 in QR  (here 2*QNR = "
          f"{'QNR' if doubled == QNR else 'QR' if doubled == QR else '??'})",
          (doubled == QNR) == two_qr)
    idem = pmul(theta, theta, XP1) == theta
    check(f"(3) theta^2 = theta mod (x^p - 1) over GF(2): {idem}; must equal "
          f"'2 in QR' = {two_qr}", idem == two_qr)
    if not idem:
        check("CONTROL: theta^2 = theta_QR instead (Frobenius sends QNR-sum "
              "to QR-sum when 2 in QNR)", pmul(theta, theta, XP1) == theta_R)

    # cyclic code generated by theta = span of its p cyclic shifts
    def cshift(v, s):
        return ((v << s) | (v >> (p - s))) & ((1 << p) - 1) if s else v

    Ntheta = rref([cshift(theta, s) for s in range(p)])
    dim_target = (p + 1) // 2
    if idem:
        check(f"(3) dim <theta> = (p+1)/2 = {dim_target}",
              len(Ntheta) == dim_target, f"dim = {len(Ntheta)}")
    else:
        check(f"CONTROL: dim <theta> = {len(Ntheta)} != (p+1)/2 = {dim_target}"
              " (no QR code without 2 in QR)", len(Ntheta) != dim_target)

    # ---- parity-extension equality:  C(H) =? {(c, wt(c) mod 2): c in <theta>}
    ext_basis = rref([b | ((bin(b).count('1') % 2) << INF) for b in Ntheta])
    check("(b) C(H) = parity-extension of <theta>  [exact subspace equality]",
          CH == ext_basis)
    all_even = all(bin(w).count("1") % 2 == 0 for w in span_words(CH))
    check("(4) every codeword of C(H) has EVEN total weight (inf coordinate "
          "is an overall parity check -- matches the literature's extension)",
          all_even)
    punctured = rref([b & ((1 << p) - 1) for b in CH])
    check("(4) puncturing C(H) at inf gives exactly <theta>", punctured == Ntheta)

    if not idem:
        wd = wdist(span_words(CH))
        print(f"CONTROL p={p}: weight distribution of C(H): {wd}")
        check(f"CONTROL p={p} (p=3 mod 8): C(H) is the FULL even-weight "
              f"[{NB},{NB-1}] code, NOT an extended QR code",
              len(CH) == NB - 1 and all_even)
        return

    # ---- audit point (3) continued: WHICH QR code does theta generate?
    # Build GF(2^m) from an irreducible factor of (x^p-1)/(x+1); alpha = x.
    cof = XP1
    # divide out (x+1)
    q, r = 0, cof
    while r and r.bit_length() >= 2:
        sh = r.bit_length() - 2
        q |= 1 << sh
        r ^= 0b11 << sh
    assert r == 0
    facs = irreducible_factors(q)
    degs = sorted(f.bit_length() - 1 for f in facs)
    print(f"(x^p-1)/(x+1) factors into irreducibles of degrees {degs}")
    f1 = facs[0]
    m = f1.bit_length() - 1
    print(f"GF(2^{m}) = GF(2)[y]/({poly_str(f1)}),  alpha = y")

    def fmul(a, b):
        return pmul(a, b, f1)

    apow = [1] * p
    for k in range(1, p):
        apow[k] = fmul(apow[k - 1], 2)  # alpha = y = bitmask 2
    check("alpha is a primitive p-th root of unity (alpha^p = 1, alpha != 1)",
          fmul(apow[p - 1], 2) == 1 and apow[1] != 1)

    def evaluate(polymask, point):
        """Evaluate GF(2)[x] poly at a GF(2^m) point (Horner)."""
        acc = 0
        for k in range(polymask.bit_length() - 1, -1, -1):
            acc = fmul(acc, point)
            if (polymask >> k) & 1:
                acc ^= 1
        return acc

    # q(x) = prod_{r in QR} (x - alpha^r),  n(x) = prod_{s in QNR}(x - alpha^s)
    def root_product(exps):
        poly = [1]  # coefficient list over GF(2^m), index = degree
        for e in exps:
            root = apow[e]
            new = [0] * (len(poly) + 1)
            for d, c in enumerate(poly):
                new[d + 1] ^= c            # x * c x^d
                new[d] ^= fmul(c, root)    # (+ root) * c x^d  (char 2)
            poly = new
        return poly

    qcoef = root_product(QR)
    ncoef = root_product(QNR)
    q_in_gf2 = all(c in (0, 1) for c in qcoef)
    n_in_gf2 = all(c in (0, 1) for c in ncoef)
    check("(2)+(3) q(x)=prod_{QR}(x-alpha^r) has GF(2) coefficients "
          "(needs QR closed under doubling, i.e. 2 in QR)", q_in_gf2)
    check("(2)+(3) n(x)=prod_{QNR}(x-alpha^s) has GF(2) coefficients", n_in_gf2)
    qmask = sum((c << d) for d, c in enumerate(qcoef))
    nmask = sum((c << d) for d, c in enumerate(ncoef))
    print(f"q(x) = {poly_str(qmask)}")
    print(f"n(x) = {poly_str(nmask)}")
    check("(x+1) q(x) n(x) = x^p - 1", pmul(0b11, pmul(qmask, nmask)) == XP1)

    # zero sets of theta and theta_R among alpha-powers
    zset_theta = sorted(e for e in range(p) if evaluate(theta, apow[e]) == 0)
    zset_thR = sorted(e for e in range(p) if evaluate(theta_R, apow[e]) == 0)
    print(f"zeros of theta_QNR among alpha^e: e in {zset_theta}")
    print(f"zeros of theta_QR  among alpha^e: e in {zset_thR}")
    check("(3) WHICH CODE: theta_QNR vanishes exactly on the NONRESIDUE "
          "powers of alpha  =>  <theta_QNR> = N = <n(x)>, dim (p+1)/2",
          zset_theta == QNR)
    check("(3) control: theta_QR vanishes exactly on the RESIDUE powers  =>  "
          "<theta_QR> = Q = <q(x)>", zset_thR == QR)
    check("(3) divisibility cross-check: n(x) | theta_QNR and q(x) !| theta_QNR",
          pmod(theta, nmask) == 0 and pmod(theta_R, nmask) != 0
          and pmod(theta, qmask) != 0)
    # Gauss periods over GF(2): eta_R + eta_N = 1 (sum over all nonzero = 1)
    eta_R = 0
    for r in QR:
        eta_R ^= apow[r]
    eta_N = 0
    for s in QNR:
        eta_N ^= apow[s]
    print(f"Gauss periods: eta_R = {eta_R}, eta_N = {eta_N} "
          f"(in GF(2^{m}); eta_R + eta_N must be 1)")
    check("eta_R + eta_N = 1 and {eta_R, eta_N} = {0, 1} (p = -1 mod 8: "
          "periods are RATIONAL, i.e. in GF(2) -- this is '2 in QR' again)",
          (eta_R ^ eta_N) == 1 and {eta_R, eta_N} == {0, 1})

    # idempotent-of-N check in the strict sense: theta in <n>, theta generates
    Ncode = rref([nmask << i for i in range(p - (nmask.bit_length() - 1))])
    Qcode = rref([qmask << i for i in range(p - (qmask.bit_length() - 1))])
    check("(3) <theta_QNR> = <n(x)> as subspaces (theta is THE idempotent "
          "generator of N)", Ntheta == Ncode)
    check("dim N = dim Q = (p+1)/2",
          len(Ncode) == dim_target and len(Qcode) == dim_target)

    # ---- main claim: C(H) = ext(N); and ext(N) vs ext(Q)
    def parity_extend(basis):
        return rref([b | ((bin(b).count('1') % 2) << INF) for b in basis])

    eN = parity_extend(Ncode)
    eQ = parity_extend(Qcode)
    check("MAIN: C(H) = extended QR code ext(N) (the QNR/theta one) "
          "[exact subspace equality]", CH == eN)
    check("ext(N) != ext(Q) as subspaces (so 'THE extended QR code' is only "
          "defined up to equivalence; the sketch should name ext(N))",
          eN != eQ)
    print(f"    C(H) == ext(Q)? {CH == eQ}")
    # multiplier equivalence j -> s*j (s in QNR), inf fixed, maps ext(N)->ext(Q)
    s0 = QNR[0]

    def mult_perm(v):
        w = v & (1 << INF)
        for j in range(p):
            if (v >> j) & 1:
                w |= 1 << ((s0 * j) % p)
        return w

    check(f"equivalence: coordinate map j -> {s0}*j (inf fixed) maps ext(N) "
          "onto ext(Q)", rref([mult_perm(v) for v in eN]) == eQ)

    # ---- Type II checks
    words = span_words(CH)
    wd = wdist(words)
    print(f"weight distribution of C(H): {wd}")
    check(f"dim C(H) = (p+1)/2 = {dim_target}", len(CH) == dim_target)
    check("self-dual: dim = n/2 and all basis inner products even",
          len(CH) == NB // 2 and all(bin(a & b).count("1") % 2 == 0
                                     for a in CH for b in CH))
    check("doubly even (all weights = 0 mod 4) => Type II for p = -1 mod 8",
          all(w % 4 == 0 for w in wd))
    if p == 7:
        check("p=7: C(H) is the [8,4,4] code e8-hat, wd {0:1,4:14,8:1}",
              wd == {0: 1, 4: 14, 8: 1})
    if p == 23:
        check("p=23: C(H) is [24,12,8] (extended Golay; unique by Pless 1968)",
              wd == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1})


# ================================================================== run all
print("ADVERSARIAL AUDIT of THM-481 eQR proof sketch")
print("main cases p = 7, 23 (p = 7 mod 8); controls p = 11 (3 mod 8: kills "
      "idempotency), p = 17 (1 mod 8: kills skewness AND parity step)")

for p in (7, 23, 11, 17):
    audit(p)

print("""
========================================================================
ATTRIBUTION (web check, 2026-06-11 -- the identity is NOT new):
  Jon-Lark Kim, Patrick Sole, "Skew Hadamard designs and their codes",
  Designs, Codes and Cryptography 49 (2008) 135-145,
  DOI 10.1007/s10623-008-9173-y (open access).
    Proposition 2: for a skew Hadamard design A of parameters
      (4n-1, 2n-1, n-1) over GF(p) with p | n, the codes C(A) (= span of A
      extended by an ALL-ONE column) and D(A) (= span of I+A extended and
      augmented) are SELF-DUAL.
    Corollary 1: over GF(2) these codes are TYPE II; moreover
      <(I+A^T)~> is the even part of C(A) and C(A) = <(I+A^T)~, 1>.
    Proposition 3: l := 4n-1 prime, p | n prime, p a QR mod l, H the Paley
      Hadamard matrix of order 4n, A its associated incidence matrix:
      "C(A)^ is the extended quadratic residue code Q^ over GF(p)."
      For p = 2: 2|n and 2 in QR(l)  <=>  l = 7 mod 8 -- exactly the
      sketch's claim (their proof = the same idempotent route).
    Their idempotent reference: MacWilliams & Sloane, The Theory of
      Error-Correcting Codes, Chap. 16, Theorem 4 (idempotent of the QR
      code over GF(p); the Gauss-sum normalization theta=+1 fixes the
      QR-vs-QNR labeling -- 'conventions matter' resolved there).
    Their order-8 case (S 4.1): C(A) = binary Hamming [8,4,4]; order-48
      (S 4.10): QR codes over GF(2) and GF(3) since 2,3 in QR(47).
  Earlier program (codes of Hadamard designs, p-rank):
    E. F. Assmus Jr., J. D. Key, "Hadamard matrices and their designs: a
    coding-theoretic approach", Trans. AMS 330 (1992) 269-293; and
    "Designs and their Codes", CUP 1992 (S 2.10 covers QR codes) --
    cited by Kim-Sole as [1,2]; the (p+1)/2 rank statement they credit to
    T. S. Michael.  The Assmus-Key 'update' survey contains no Paley-QR
    statement (checked); the sharp statement is Kim-Sole Prop. 3.
  Type II / self-dual iff p = -1 mod 8: classical, e.g. MacWilliams-Sloane
    Ch. 16 / Huffman-Pless Fundamentals Ch. 6 (for p = +1 mod 8 the
    extended QR code is only formally self-dual).
========================================================================
""")

print("=" * 72)
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S):")
    for f in FAILURES:
        print("  -", f)
    sys.exit(1)
else:
    print("RESULT: ALL CHECKS PASSED -- sketch confirmed at p=7, 23 with the")
    print("discrepancies/clarifications noted in the printout (which-code,")
    print("loss-matrix orientation, controls showing where each congruence")
    print("hypothesis is needed).")
    sys.exit(0)

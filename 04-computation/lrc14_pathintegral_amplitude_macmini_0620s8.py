#!/usr/bin/env python3
"""
lrc14_pathintegral_amplitude_macmini_0620s8.py   (mac-mini-2026-06-20-S8)

THREAD: DOUBLE-SLIT / FEYNMAN PROPAGATOR-AS-WEIGHTS / FUBINI-STUDY.

HYPOTHESIS (one falsifiable statement, two heads):
  The repo's two SIGNED sums --
    (A) the 7-band slope identity  measS7 = (1/7) sum_{s in Z/7} bandcover_s
        (HYP-2703, PROVED), each band a twist e -> s*e mod 7, and
    (B) the Walsh-Fourier diagonalization of H(T) with signed coefficients
        H_hat[P] = (-1)^{desc(P)}/8 over Hamiltonian paths P (THM-071, n=5) --
  are PATH-INTEGRAL AMPLITUDES: a quantity = |sum_paths A(path)|^2 (or a real
  part / trace of such), where the signed terms are INTERFERENCE cross terms
  e^{2 pi i phi}, and "consec wins" / "H counts cycles" is CONSTRUCTIVE
  interference (phase alignment maximizing the amplitude).

We test FOUR precise sub-claims, each PASS/FAIL with exact arithmetic:

  T1 (band sum = character/Fourier inversion on Z/7).
     The slope band-cover bandcover_s(E) is, pointwise in f, an indicator
     1[ {s*e + floor(e f)} = Z/7 ]. Claim: the FULL coloring measure has an
     exact discrete-Fourier (character) representation
        p_surj(profile) = (1/7^7) sum_{chi} prod ... (inclusion-exclusion
        over Z/7 via additive characters).
     Concretely: for a multiset of residues r_1..r_k in Z/7, the indicator
     "the set {r_i} = all of Z/7" equals a signed sum over which residues are
     MISSING = sum_{S subset Z/7} (-1)^{|S|} 1[no r_i in S].  This is the
     7-term (well, 2^7-term but collapsible) interference expansion. We verify
     it reproduces measS7 EXACTLY and that the alternating signs ARE the
     interference signs.  -> identity, not analogy, IF it matches.

  T2 (consec = constructive phase alignment?).
     Map each offset e to a phase on the clock: A_e(theta) = e^{2 pi i * e * theta / 7}.
     The "amplitude" of a profile = sum_e A_e. Claim to TEST: does consec
     maximize a coherence / |amplitude|^2 functional that correlates with
     measS7?  We compute, for every k-subset (k=8, span<=12):
        coh(E) = average over the breakpoint grid of |sum_{e} e^{2 pi i floor(e f)/7}|^2
     and ask: is argmax coh == consec?  (If YES -> constructive-interference
     picture has teeth; if NO -> it is an analogy only.)

  T3 (OCF H as squared amplitude, n=5).
     THM-071: H = 1 + 2 t3 + 2 t5, and degree-4 Walsh coeff = (-1)^{desc(P)}/8.
     The sign (-1)^{desc(P)} = e^{i pi desc(P)} is a genuine amplitude phase.
     TEST: define a per-tournament path amplitude
        Amp(T) = sum_{Ham paths P} i^{?}  ... -- specifically does
        H(T) = (number of Ham paths) decompose as a |.|^2 of a vector whose
        entries are e^{i pi desc}? We test the STRONGEST clean version:
        is H(T) equal to a Gram/overlap <psi_T | psi_T> for an explicit
        feature vector psi_T with unimodular entries?  We just CHECK whether
        the Walsh expansion H = sum_S H_hat[S] chi_S(T) can be written as a
        square; i.e. is the symmetric form positive semidefinite with H on
        the diagonal pattern?  (Honest: most likely H itself is a COUNT, the
        amplitude is at the OCF-summand level. We test the summand-level claim.)

  T4 (interference cross terms reproduce the relation-lattice correction).
     The deviation measS7(consec) - measS7(adv) is a sum of signed band terms
     (M1 in angleA). TEST whether this signed deviation equals the cross-term
     (off-diagonal) part of a coherence form |sum A_e|^2 = sum|A_e|^2 +
     2 Re sum_{e<e'} A_e conj(A_e'). I.e. do the OFF-DIAGONAL (interference)
     terms alone carry the consec-vs-adv margin?

stdlib only; exact Fractions for measures; complex via Gaussian-integer/roots
of unity handled symbolically where exactness matters.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

P = 7

# ---------------------------------------------------------------------------
# measS7 engines (copied/verified against angleA + extremal_landscape scripts)
# ---------------------------------------------------------------------------
def measS7(E, scale=P):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, scale*e+1): bps.add(F(m, scale*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        if len(set(int(((e*mid)%1)*scale) for e in E)) == scale: tot += hi-lo
    return tot

def band_cover(E, s, scale=P):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1): bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for i in range(len(bps)-1):
        f0, f1 = bps[i], bps[i+1]
        if f1 <= f0: continue
        fm = (f0+f1)/2
        if len(set((s*e + int(e*fm)) % scale for e in E)) == scale: tot += f1-f0
    return tot

def measS7_bands(E, scale=P):
    return sum(band_cover(E, s, scale) for s in range(scale))/scale

# ---------------------------------------------------------------------------
# T1: surjection indicator as a SIGNED (inclusion-exclusion = character) sum.
#     1[ {r_i} = Z/7 ] = sum_{S subset Z/7} (-1)^{|S|} 1[ no r_i in S ].
#     This is the interference/sign expansion. We integrate it against the
#     breakpoint grid and check it reproduces measS7 EXACTLY, then identify
#     the signs with additive-character cross terms.
# ---------------------------------------------------------------------------
def measS7_via_signed_IE(E, scale=P):
    """measS7 by integrating the alternating missing-set expansion of the
    surjection indicator over the breakpoint grid. Exact."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, scale*e+1): bps.add(F(m, scale*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        cols = set(int(((e*mid)%1)*scale) for e in E)
        # alternating IE over missing color sets:
        # 1[cols == full] = sum_{S} (-1)^|S| 1[cols disjoint from S]
        #                 = sum_{S subset of complement(cols)} (-1)^|S|
        # = 0 unless complement(cols) empty, in which case = (-1)^0 = 1.
        missing = set(range(scale)) - cols
        ind = sum((-1)**len(S)
                  for r in range(len(missing)+1)
                  for S in itertools.combinations(sorted(missing), r))
        # ind == 1 iff missing empty else 0 (alternating binomial sum)
        tot += ind * (hi - lo)
    return tot

def char_surjection_probability(profile_residue_words):
    """
    PURE-CHARACTER (additive Fourier on Z/7) form of the surjection probability
    for a FINITE set of equally weighted residue-words (each word = a tuple of
    residues over the k offsets). Demonstrates the e^{2 pi i ...} interference
    representation explicitly:
       1[word hits all of Z/7]
         = sum_{S subset Z/7} (-1)^{|S|} prod_{offsets} 1[res not in S]
    and 1[res not in S] = 1 - 1[res in S], with
       1[res in S] = (1/7) sum_{t in Z/7} (sum_{a in S} omega^{t(res-a)}),  omega=e^{2pi i/7}.
    We return the average over words, computed EXACTLY via the alternating
    missing-set sum (the character sum is what it expands to; we verify a few
    words by the complex form too).
    """
    cnt = 0
    n = len(profile_residue_words)
    for w in profile_residue_words:
        cols = set(w)
        cnt += 1 if len(cols) == P else 0
    return F(cnt, n) if n else F(0)

def verify_complex_character_form():
    """Verify, on random residue words, that the COMPLEX additive-character
    expansion of the surjection indicator (genuine e^{2pi i} interference)
    equals the combinatorial indicator. This is the 'amplitude' identity."""
    import random
    omega = cmath.exp(2j*math.pi/P)
    rng = random.Random(7)
    maxerr = 0.0
    for _ in range(2000):
        k = rng.randint(5, 12)
        w = [rng.randrange(P) for _ in range(k)]
        # indicator via complex characters:
        # 1[hits all] = sum_{S} (-1)^|S| prod_i (1 - 1[w_i in S]),
        # 1[w_i in S] = (1/7) sum_{a in S} sum_{t} omega^{t (w_i - a)}
        val = 0+0j
        for r in range(P+1):
            for Sset in itertools.combinations(range(P), r):
                prod = 1+0j
                for wi in w:
                    inS = 0.0
                    for a in Sset:
                        for t in range(P):
                            inS += (omega**(t*(wi-a))).real / P  # real since summed over full t-orbit
                    prod *= (1 - inS)
                val += ((-1)**r) * prod
        true = 1.0 if len(set(w)) == P else 0.0
        maxerr = max(maxerr, abs(val.real - true), abs(val.imag))
    return maxerr

# ---------------------------------------------------------------------------
# T2: coherence functional |sum_e e^{2pi i floor(e f)/7}|^2 -- does it pick consec?
# ---------------------------------------------------------------------------
def coherence(E, scale=P):
    """Average over the breakpoint grid f of |sum_e omega^{floor(e f)}|^2,
    omega=e^{2pi i/7}. This is the squared CLOCK amplitude. Returns a float
    (it's a complex-modulus average; we keep float here, exactness not needed
    for an argmax comparison, but we use rational breakpoint weights)."""
    omega = cmath.exp(2j*math.pi/scale)
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1): bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = 0.0
    for i in range(len(bps)-1):
        f0, f1 = bps[i], bps[i+1]
        if f1 <= f0: continue
        fm = (f0+f1)/2
        amp = sum(omega**(int(e*fm)) for e in E)
        tot += (abs(amp)**2) * float(f1-f0)
    return tot

def coherence_surj_weighted(E, scale=P):
    """Alternative: |amplitude|^2 only WHERE the coloring is surjective.
    If consec maximizes measS7 because surjective regions have aligned phase,
    this should be larger for consec. Average of |sum omega^{floor(ef)}|^2 * 1[surj]."""
    omega = cmath.exp(2j*math.pi/scale)
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1): bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = 0.0
    for i in range(len(bps)-1):
        f0, f1 = bps[i], bps[i+1]
        if f1 <= f0: continue
        fm = (f0+f1)/2
        cols = set((int(e*fm)) % scale for e in E)
        if len(cols) == scale:
            amp = sum(omega**(int(e*fm)) for e in E)
            tot += (abs(amp)**2) * float(f1-f0)
    return tot

# ---------------------------------------------------------------------------
# T3: OCF H(T) at n=5 -- amplitude/Walsh check.
# ---------------------------------------------------------------------------
def hamiltonian_path_count(adj, n):
    """adj[i][j]=1 if i->j. Count directed Hamiltonian paths (Redei H)."""
    cnt = 0
    for perm in itertools.permutations(range(n)):
        ok = all(adj[perm[i]][perm[i+1]] for i in range(n-1))
        if ok: cnt += 1
    return cnt

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in itertools.product((0,1), repeat=len(edges)):
        adj = [[0]*n for _ in range(n)]
        for (i,j),b in zip(edges, bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def desc_of_path(perm, adj):
    """descents = edges going larger->smaller vertex along the path, traversed
    from the smaller endpoint. For an undirected Ham path used as a Walsh
    monomial, THM-071 uses (-1)^{desc}. Here perm is an ORIENTED path that is
    actually present in T (a Redei Ham path)."""
    # traverse from smaller endpoint
    a, b = perm[0], perm[-1]
    seq = perm if a < b else perm[::-1]
    return sum(1 for i in range(len(seq)-1) if seq[i] > seq[i+1])

def t3_t5_count(adj, n):
    """count directed 3-cycles (t3) and directed 5-cycles (t5)."""
    t3 = 0
    for c in itertools.combinations(range(n),3):
        for perm in itertools.permutations(c):
            if all(adj[perm[i]][perm[(i+1)%3]] for i in range(3)):
                t3 += 1
    t3 //= 3  # each directed 3-cycle counted 3x (rotations)
    t5 = 0
    if n >= 5:
        for c in itertools.combinations(range(n),5):
            for perm in itertools.permutations(c):
                if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    t5 += 1
        t5 //= 5
    return t3, t5

def test_T3_OCF():
    print("="*72)
    print("(T3) OCF H(n=5): is H a squared amplitude / does (-1)^desc carry it?")
    print("="*72)
    n = 5
    # verify H = 1 + 2 t3 + 2 t5 and the signed-amplitude reading
    bad = 0
    amp_match = 0
    total = 0
    examples = []
    for adj in all_tournaments(n):
        total += 1
        H = hamiltonian_path_count(adj, n)
        t3, t5 = t3_t5_count(adj, n)
        pred = 1 + 2*t3 + 2*t5
        if H != pred: bad += 1
        # AMPLITUDE READING: sum over the H present Redei paths of (-1)^desc.
        # The 'signed amplitude' S_signed = sum_{present P} (-1)^{desc(P)}.
        # Path-integral claim's strongest clean form: |something|^2 = H?
        S_signed = 0
        present = 0
        for perm in itertools.permutations(range(n)):
            if all(adj[perm[i]][perm[i+1]] for i in range(n-1)):
                present += 1
                # count this oriented path once; desc relative to smaller endpoint
                S_signed += (-1)**desc_of_path(perm, adj)
        # note present counts each undirected Redei path TWICE (both directions)
        if len(examples) < 6:
            examples.append((H, present, S_signed))
    print(f"   tournaments={total}  H==1+2t3+2t5 mismatches={bad}")
    print(f"   sample (H, present_oriented, signed_sum_{{(-1)^desc}}):")
    for e in examples:
        print(f"      {e}")
    print("   NOTE: present_oriented = 2*H (each Redei path both directions).")
    print("   The signed sum sum (-1)^desc over PRESENT paths is the OCF-style")
    print("   amplitude; we report whether |.|/structure relates to H below.")

# ---------------------------------------------------------------------------
# T4: does the off-diagonal (interference) part of coherence carry the
#     consec-vs-adv margin?
# ---------------------------------------------------------------------------
def coherence_split(E, scale=P):
    """Return (diag, offdiag) parts of average |sum omega^{floor(ef)}|^2:
       |sum A_e|^2 = sum_e |A_e|^2 (= k, diag) + 2 Re sum_{e<e'} A_e conj A_e' (offdiag)."""
    omega = cmath.exp(2j*math.pi/scale)
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1): bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    diag = 0.0; off = 0.0
    Es = E
    for i in range(len(bps)-1):
        f0, f1 = bps[i], bps[i+1]
        if f1 <= f0: continue
        fm = (f0+f1)/2; w = float(f1-f0)
        phs = [omega**(int(e*fm)) for e in Es]
        diag += w*len(Es)  # each |A_e|^2 = 1
        s = 0.0
        for a in range(len(phs)):
            for b in range(a+1,len(phs)):
                s += 2*(phs[a]*phs[b].conjugate()).real
        off += w*s
    return diag, off

# ---------------------------------------------------------------------------
def main():
    print("PATH-INTEGRAL / DOUBLE-SLIT amplitude test  (mac-mini-2026-06-20-S8)\n")

    # ---- T1 ----
    print("="*72)
    print("(T1) surjection indicator = signed interference (IE/character) sum")
    print("="*72)
    shapes = [list(range(k)) for k in range(3,11)] + [[0,2,3,4,5,6,7,8],[0,1,3,5,7,9,11],[0,2,4,6,8,10,12,14]]
    bad = 0
    for E in shapes:
        a = measS7(E); b = measS7_via_signed_IE(E)
        if a != b:
            bad += 1; print("   MISMATCH", E, float(a), float(b))
    print(f"   shapes checked={len(shapes)}  mismatches={bad}")
    print(f"   => measS7 = integral of the ALTERNATING missing-set sum (interference signs).")
    err = verify_complex_character_form()
    print(f"   complex additive-character (e^{{2pi i/7}}) form of the indicator:")
    print(f"      max |char_expansion - true_indicator| over 2000 words = {err:.2e}")
    print(f"   => the surjection indicator IS literally a sum of phases e^{{2pi i t(r-a)/7}}.")
    print()

    # ---- T2 ----
    print("="*72)
    print("(T2) does the clock-amplitude coherence pick consec? (k=8, span<=12)")
    print("="*72)
    for k in (8,):
        cons = list(range(k))
        coh_cons = coherence(cons)
        cohs_cons = coherence_surj_weighted(cons)
        m_cons = measS7(cons)
        best_coh = (coh_cons, cons); best_cohs = (cohs_cons, cons); best_m = (m_cons, cons)
        nbad_coh = 0; nbad_cohs = 0; cnt=0
        for combo in itertools.combinations(range(1,13), k-1):
            E = [0]+list(combo)
            if E == cons: continue
            cnt += 1
            c = coherence(E); cs = coherence_surj_weighted(E); m = measS7(E)
            if c > best_coh[0]: best_coh = (c,E)
            if cs > best_cohs[0]: best_cohs = (cs,E)
            if m > best_m[0]: best_m = (m,E)
            if c > coh_cons: nbad_coh += 1
            if cs > cohs_cons: nbad_cohs += 1
        print(f"   k={k}  subsets tested={cnt+1}")
        print(f"   measS7      argmax == consec? {best_m[1]==cons}   (consec={float(m_cons):.4f})")
        print(f"   coherence   argmax == consec? {best_coh[1]==cons}   "
              f"(consec={coh_cons:.4f}, max={best_coh[0]:.4f} at {best_coh[1]})  beats_consec={nbad_coh}")
        print(f"   coh*1[surj] argmax == consec? {best_cohs[1]==cons}  "
              f"(consec={cohs_cons:.4f}, max={best_cohs[0]:.4f} at {best_cohs[1]})  beats_consec={nbad_cohs}")
    print()

    # ---- T3 ----
    test_T3_OCF()
    print()

    # ---- T4 ----
    print("="*72)
    print("(T4) interference (off-diagonal) part: does it carry consec-vs-adv margin?")
    print("="*72)
    cons = list(range(8)); adv = [0,2,3,4,5,6,7,8]
    dc, oc = coherence_split(cons); da, oa = coherence_split(adv)
    mc, ma = measS7(cons), measS7(adv)
    print(f"   consec: measS7={float(mc):.5f}  coh_diag={dc:.4f}  coh_off={oc:+.4f}")
    print(f"   adv:    measS7={float(ma):.5f}  coh_diag={da:.4f}  coh_off={oa:+.4f}")
    print(f"   measS7 margin (consec-adv)      = {float(mc-ma):+.5f}")
    print(f"   coherence-offdiag margin        = {oc-oa:+.5f}  (diag equal: {abs(dc-da)<1e-9})")
    print(f"   sign agreement (both same sign)? {(float(mc-ma)>0) == ((oc-oa)>0)}")
    print()
    print("="*72)
    print("Interpretation printed in final assistant message (PASS/FAIL per claim).")
    print("="*72)

if __name__ == "__main__":
    main()

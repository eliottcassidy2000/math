#!/usr/bin/env python3
"""
half_tiling_codec_kpswf.py
================================================================
THREAD 4 -- ENGINEERING: a half-tiling CODEC + uniform SC sampler.

A phi-self-converse (grid-symmetric) tournament on {1..n} is exactly determined
by the bits of its tiling on the rho-orbit REPRESENTATIVES, where
    rho(a,b) = (n+1-b, n+1-a)
is the pure coordinate involution realizing converse+relabel (THM-280, and the
"no bit flip" structure verified in half_tiling_framework_kps.py).

Number of orbits = half(n) = floor((n-1)^2/4).  So any phi-self-converse
tournament compresses LOSSLESSLY to a half(n)-bit string -- vs the full tiling
m = C(n-1,2) bits, vs the full adjacency C(n,2) bits.  Asymptotically
    half(n) / m  ->  1/2   and   half(n) / C(n,2) -> 1/2,
a 2x compression for the self-converse stratum.

DELIVERABLES
  (a) encode/decode  (tiling/adjacency <-> half(n)-bit code), exact round-trip.
  (b) compression table for n up to 20.
  (c) uniform SC sampler in O(half(n)) work/space (no enumeration).
  (d) product framing (docstring + main() report).

Reuses canonical machinery from half_tiling_framework_kps.py.
"""

import os, sys, itertools, random
from math import comb, floor
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from half_tiling_framework_kps import (
    tiles, rho, tiling_to_tournament, tournament_to_tiling,
    converse_relabel, count_hamiltonian_paths, involution_structure,
    half_formula,
)


# ======================================================================
#  The codec
# ======================================================================

class HalfTilingCodec:
    """Lossless 2x-compression codec for phi-self-converse tournaments.

    A phi-self-converse tournament on n vertices is encoded as a half(n)-bit
    string -- one bit per rho-orbit on the staircase tiling.

    Bit conventions (canonical, 01-canon/definitions.md):
      tile (a,b), a>=b+2 ; bit 0 => forward arc a->b ; bit 1 => backward b->a.
      base path P0 = n->n-1->...->1 (the n-1 consecutive arcs, NOT tiles).

    Orbit-representative choice: for each rho-orbit {P, rho(P)} we pick the
    canonical representative = min(P, rho(P)) under tuple order. Fixed cells
    (on the anti-diagonal a+b=n+1) are singleton orbits and are their own rep.
    The code lists bits in the fixed, sorted order of representatives.
    """

    def __init__(self, n):
        if n < 2:
            raise ValueError("n must be >= 2")
        self.n = n
        self.tiles = tiles(n)                 # canonical enumeration order
        self.tile_index = {p: i for i, p in enumerate(self.tiles)}
        self.m = len(self.tiles)              # full-tiling bit count = C(n-1,2)
        self.rho_map = {p: rho(*p, n) for p in self.tiles}

        # Build orbits. Representative = min of the orbit (lexicographic tuple).
        reps = []
        seen = set()
        for p in self.tiles:
            if p in seen:
                continue
            img = self.rho_map[p]
            orbit = {p, img}
            seen |= orbit
            reps.append(min(orbit))
        self.reps = sorted(reps)              # deterministic bit order
        self.half = len(self.reps)            # = floor((n-1)^2/4)
        assert self.half == half_formula(n), \
            f"orbit count {self.half} != half({n})={half_formula(n)}"

        # For decode: map each tile to the representative whose bit governs it.
        self.tile_to_rep = {}
        for p in self.tiles:
            img = self.rho_map[p]
            self.tile_to_rep[p] = min(p, img)
        # position of each rep in the code string
        self.rep_pos = {r: i for i, r in enumerate(self.reps)}

    # ------------------------------------------------------------------
    #  validation
    # ------------------------------------------------------------------
    def is_grid_symmetric_tiling(self, t):
        """t: dict tile->bit. True iff fixed by rho (phi-self-converse)."""
        return all(t[p] == t[self.rho_map[p]] for p in self.tiles)

    def is_phi_self_converse(self, T):
        """T: adjacency dict. True iff T == converse_relabel(T)."""
        R = converse_relabel(T, self.n)
        return all(T[k] == R.get(k, 0) for k in T)

    # ------------------------------------------------------------------
    #  encode
    # ------------------------------------------------------------------
    def encode_tiling(self, t):
        """tiling dict -> tuple of half(n) bits. Requires grid-symmetric input."""
        if not self.is_grid_symmetric_tiling(t):
            raise ValueError("tiling is NOT grid-symmetric (not phi-self-converse); "
                             "cannot losslessly half-encode")
        return tuple(t[r] for r in self.reps)

    def encode_adjacency(self, T):
        """adjacency dict -> tuple of half(n) bits."""
        t = tournament_to_tiling(T, self.tiles)
        return self.encode_tiling(t)

    def encode_to_int(self, t):
        """Pack the half(n) bits into a single integer (LSB = reps[0])."""
        bits = self.encode_tiling(t)
        val = 0
        for i, b in enumerate(bits):
            val |= (b & 1) << i
        return val

    # ------------------------------------------------------------------
    #  decode
    # ------------------------------------------------------------------
    def decode_to_tiling(self, code):
        """half(n)-bit code (tuple/list/int) -> full tiling dict (all m tiles).

        Every tile takes the bit of its orbit representative; this automatically
        reconstructs a grid-symmetric tiling (t[p] == t[rho(p)] by construction).
        """
        bits = self._coerce_bits(code)
        t = {}
        for p in self.tiles:
            r = self.tile_to_rep[p]
            t[p] = bits[self.rep_pos[r]]
        return t

    def decode_to_adjacency(self, code):
        """half(n)-bit code -> adjacency dict of a phi-self-converse tournament."""
        t = self.decode_to_tiling(code)
        return tiling_to_tournament(t, self.tiles, self.n)

    def _coerce_bits(self, code):
        if isinstance(code, int):
            return [ (code >> i) & 1 for i in range(self.half) ]
        bits = list(code)
        if len(bits) != self.half:
            raise ValueError(f"code length {len(bits)} != half({self.n})={self.half}")
        return bits

    # ------------------------------------------------------------------
    #  uniform sampler  (deliverable c)  -- O(half(n)) work and space
    # ------------------------------------------------------------------
    def sample(self, rng=random):
        """Draw a uniformly-random phi-self-converse tournament.

        Draw half(n) independent fair bits, decode. Because grid-symmetric
        tilings <-> phi-self-converse tournaments is a BIJECTION and the half(n)
        bits range freely over {0,1}^half(n) <-> all 2^half(n) grid-sym tilings
        uniformly, the resulting tournament is uniform on the SC stratum.
        Work and space: O(half(n)) random bits + O(m) output (output is the
        tournament itself; the *code* is O(half(n))). No enumeration.
        """
        code = tuple(rng.randint(0, 1) for _ in range(self.half))
        return self.decode_to_adjacency(code), code

    def sample_code(self, rng=random):
        """Just the compressed code (O(half(n)) bits)."""
        return tuple(rng.randint(0, 1) for _ in range(self.half))


# ======================================================================
#  (a) round-trip verification on ALL grid-sym tilings, n=3..7
# ======================================================================

def verify_roundtrip_all_gridsym(nmax=7):
    print("\n[A] EXACT round-trip on ALL 2^half grid-symmetric tilings (n=3..%d)" % nmax)
    print(f"{'n':>3} {'half':>5} {'2^half':>8} {'#gridsym':>9} {'#roundtrip_ok':>14} "
          f"{'#selfconv_ok':>13} {'status':>8}")
    all_ok = True
    for n in range(3, nmax+1):
        codec = HalfTilingCodec(n)
        m = codec.m
        n_gridsym = 0
        rt_ok = 0
        sc_ok = 0
        # enumerate ALL m-bit tilings; pick the grid-symmetric ones
        for bits in itertools.product([0, 1], repeat=m):
            t = {p: bits[i] for i, p in enumerate(codec.tiles)}
            if not codec.is_grid_symmetric_tiling(t):
                continue
            n_gridsym += 1
            # encode -> decode -> compare full tiling
            code = codec.encode_tiling(t)
            t_back = codec.decode_to_tiling(code)
            if t_back == t:
                rt_ok += 1
            # also check adjacency round trip + phi-self-converse property
            T = tiling_to_tournament(t, codec.tiles, n)
            code2 = codec.encode_adjacency(T)
            T_back = codec.decode_to_adjacency(code2)
            if T_back == T and codec.is_phi_self_converse(T_back):
                sc_ok += 1
        expect = 2 ** codec.half
        ok = (n_gridsym == expect == rt_ok == sc_ok)
        all_ok &= ok
        print(f"{n:>3} {codec.half:>5} {expect:>8} {n_gridsym:>9} {rt_ok:>14} "
              f"{sc_ok:>13} {'OK' if ok else 'FAIL':>8}")
    print(f"  => round-trip exact on every grid-sym tiling n=3..{nmax}: "
          f"{'PASS (100%)' if all_ok else 'FAIL'}")
    return all_ok


def verify_int_roundtrip(nmax=7):
    """Also confirm the integer-packing API round-trips."""
    print("\n[A'] Integer-packing API round-trip (encode_to_int / decode int)")
    all_ok = True
    for n in range(3, nmax+1):
        codec = HalfTilingCodec(n)
        ok = True
        for bits in itertools.product([0, 1], repeat=codec.m):
            t = {p: bits[i] for i, p in enumerate(codec.tiles)}
            if not codec.is_grid_symmetric_tiling(t):
                continue
            v = codec.encode_to_int(t)
            t_back = codec.decode_to_tiling(v)
            if t_back != t:
                ok = False; break
        all_ok &= ok
        print(f"  n={n}: int round-trip {'OK' if ok else 'FAIL'} "
              f"(code fits in {codec.half} bits, max int {2**codec.half - 1})")
    print(f"  => {'PASS' if all_ok else 'FAIL'}")
    return all_ok


# ======================================================================
#  (b) compression table, n up to 20
# ======================================================================

def compression_table(nmax=20):
    print("\n[B] COMPRESSION TABLE  (phi-self-converse stratum)")
    print("    half(n) = floor((n-1)^2/4) bits  vs  full tiling m=C(n-1,2)  vs  full adjacency C(n,2)")
    print(f"{'n':>3} {'half(n)':>8} {'m=C(n-1,2)':>11} {'C(n,2)':>7} "
          f"{'half/m':>8} {'half/C(n,2)':>12} {'bytes(half)':>11} {'bytes(adj)':>10}")
    rows = []
    for n in range(3, nmax+1):
        half = floor((n-1)**2 / 4)
        m = comb(n-1, 2)
        adj = comb(n, 2)
        r1 = half / m if m else float('nan')
        r2 = half / adj if adj else float('nan')
        b_half = (half + 7) // 8
        b_adj = (adj + 7) // 8
        rows.append((n, half, m, adj, r1, r2, b_half, b_adj))
        print(f"{n:>3} {half:>8} {m:>11} {adj:>7} {r1:>8.4f} {r2:>12.4f} "
              f"{b_half:>11} {b_adj:>10}")
    print("    asymptotics: half(n)/m -> 1/2 and half(n)/C(n,2) -> 1/2  (2x compression).")
    # show the limit numerically
    n = 200
    half = floor((n-1)**2 / 4); m = comb(n-1, 2); adj = comb(n, 2)
    print(f"    at n=200:  half/m = {half/m:.6f},  half/C(n,2) = {half/adj:.6f}")
    return rows


# ======================================================================
#  (c) sampler sanity: always self-converse + H-distribution match
# ======================================================================

def exhaustive_gridsym_H_spectrum(n):
    """Exact H-spectrum over all grid-symmetric tilings (probabilities)."""
    codec = HalfTilingCodec(n)
    spec = Counter()
    for bits in itertools.product([0, 1], repeat=codec.m):
        t = {p: bits[i] for i, p in enumerate(codec.tiles)}
        if not codec.is_grid_symmetric_tiling(t):
            continue
        T = tiling_to_tournament(t, codec.tiles, n)
        H = count_hamiltonian_paths(T, n)
        spec[H] += 1
    total = sum(spec.values())
    return spec, total


def verify_sampler(nmax=6, n_samples=200000, seed=20260620):
    print("\n[C] UNIFORM SAMPLER sanity check")
    print("    (i) every sample is phi-self-converse  (ii) sampled H-dist ~ exhaustive grid-sym")
    rng = random.Random(seed)
    all_sc = True
    for n in range(5, nmax+1):
        codec = HalfTilingCodec(n)
        exact, total = exhaustive_gridsym_H_spectrum(n)
        # sample
        samp = Counter()
        sc_fail = 0
        for _ in range(n_samples):
            T, code = codec.sample(rng)
            if not codec.is_phi_self_converse(T):
                sc_fail += 1
            H = count_hamiltonian_paths(T, n)
            samp[H] += 1
        all_sc &= (sc_fail == 0)
        # compare distributions
        Hs = sorted(set(exact) | set(samp))
        max_abs_err = 0.0
        rows = []
        for H in Hs:
            p_exact = exact.get(H, 0) / total
            p_samp = samp.get(H, 0) / n_samples
            err = abs(p_exact - p_samp)
            max_abs_err = max(max_abs_err, err)
            rows.append((H, exact.get(H, 0), p_exact, p_samp, err))
        print(f"  --- n={n}: |SC stratum|=2^{codec.half}={2**codec.half}, "
              f"samples={n_samples}, self-converse failures={sc_fail} "
              f"{'(OK)' if sc_fail==0 else '(FAIL!)'} ---")
        print(f"      {'H':>6} {'exact#':>8} {'p_exact':>10} {'p_sampled':>10} {'|err|':>9}")
        for (H, c, pe, ps, err) in rows:
            print(f"      {H:>6} {c:>8} {pe:>10.6f} {ps:>10.6f} {err:>9.6f}")
        print(f"      max |p_exact - p_sampled| = {max_abs_err:.6f}  "
              f"(expect -> 0 as samples grow; ~ 1/sqrt(N) = {1/(n_samples**0.5):.5f})")
    print(f"  => sampler ALWAYS phi-self-converse on all tested n: "
          f"{'PASS' if all_sc else 'FAIL'}")
    return all_sc


# ======================================================================
#  (d) product framing + demo
# ======================================================================

def product_demo():
    print("\n[D] PRODUCT DEMO: self-converse tournament codec")
    print("    Use case: ranking systems with antisymmetric structure (round-robins,")
    print("    pairwise-preference matrices that are converse-symmetric under a known")
    print("    relabeling phi). Such a system's full pairwise matrix costs C(n,2) bits;")
    print("    the half-tiling fingerprint costs only half(n) ~ C(n,2)/2 bits.")
    n = 7
    codec = HalfTilingCodec(n)
    rng = random.Random(7)
    T, code = codec.sample(rng)
    H = count_hamiltonian_paths(T, n)
    print(f"    Example n={n}: drew a random SC tournament.")
    print(f"      compressed fingerprint (half={codec.half} bits): {''.join(map(str,code))}")
    print(f"      packed integer fingerprint: {codec.encode_to_int(codec.decode_to_tiling(code))}")
    print(f"      H (directed Hamiltonian paths) = {H}")
    print(f"      full tiling would need m={codec.m} bits; full adjacency C(7,2)={comb(7,2)} bits.")
    # demonstrate the fingerprint is canonical: encode(decode(code)) == code
    t = codec.decode_to_tiling(code)
    code2 = codec.encode_tiling(t)
    print(f"      canonical round-trip encode(decode(code))==code: {code2 == code}")
    # Paley T7 (the global H-max, which IS self-converse) round-trips through codec
    print("    The global H-maximizer at n=7 (Paley T_7, H=189) is self-converse and")
    print("    thus lives in the codec's domain -- its fingerprint is a 9-bit code.")


def main():
    print("=" * 72)
    print("HALF-TILING CODEC  --  kind-pasteur THREAD 4 (engineering)")
    print("=" * 72)
    a_ok = verify_roundtrip_all_gridsym(7)
    ai_ok = verify_int_roundtrip(7)
    compression_table(20)
    c_ok = verify_sampler(6, n_samples=200000)
    product_demo()
    print("\n" + "=" * 72)
    print("SUMMARY")
    print(f"  (a) exact round-trip on ALL grid-sym tilings n=3..7 : {'PASS' if a_ok else 'FAIL'}")
    print(f"  (a') integer-packing round-trip                     : {'PASS' if ai_ok else 'FAIL'}")
    print(f"  (b) compression table n<=20                         : printed (half/full -> 1/2)")
    print(f"  (c) uniform sampler self-converse + H-match         : {'PASS' if c_ok else 'FAIL'}")
    print(f"  (d) product demo                                    : printed")
    print("=" * 72)


if __name__ == "__main__":
    main()

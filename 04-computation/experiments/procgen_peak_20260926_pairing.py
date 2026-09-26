#!/usr/bin/env python3
"""procgen_peak_20260926_pairing -- T5 of the peak lane: can the Collatz pairing family (THM-4470/4475)
catch at the peak?

Pairing family: pairs {2i-1, 2i}, one bit per pair.  bit 0 (Collatz): 2i-1 -> 3i-1, 2i -> i.
bit 1 (flip): 2i-1 -> i-1, 2i -> 3i.  P_L: every n >= 3 falls below itself within L steps.
Greedy constructions (processing n = 3, 4, ..., N in order, as in THM-4475):
  the default path of n uses frozen bits and reads free bits as 0; if it descends within L, the pairs of
  the points used (times 0..d-1) are frozen; otherwise n is rescued by flipping ONE free pair and the new
  path is frozen.  Rescue orders:
    'afb'  : THM-4475's rule (A = own pair, F = first w < 2n on the path that is the image of an up-move
             and w = 3 mod 4, B = pair of T(n)); order A,F,B if n = 3 mod 8, else F,B,A
    'peak' : the odd path points v = 2 (mod 3) (times 0..L-1) in DECREASING order of value (catch as high as
             possible); v = 2 (mod 3) isolates the partner v+1 = 0 (mod 3) exactly as in THM-4475 (no odd
             Collatz preimage, so only its own halving chain enters it); partner-safe: a flip at v = 2i-1 is
             accepted only if n descends AND the partner 2i (now sent up to 3i) descends below itself under
             the current bits; both certificates are frozen.  If no such flip works, THM-4475's options
             A (own pair), F, B are tried.
  If no single flip works (rare: partners 2i of earlier flips, whose own pair is already flipped), a
  depth-4 search over flips at free odd path points (earliest first) is used; a source that still fails is
  counted as stuck (the runner checks there are none).
Densities: flipped pair indices i <= N/2 (final: a flip made while processing n has index > n/2), and
i <= N/4 (THM-4475's audit convention)."""


class Bits:
    """0 free, 1 frozen at Collatz, 2 frozen flipped; bytearray up to cap, dict above."""
    __slots__ = ("cap", "arr", "over")

    def __init__(self, cap):
        self.cap = cap
        self.arr = bytearray(cap)
        self.over = {}


def run_pairing(L, N, mode, cap):
    B = Bits(cap)
    arr, over, CAP = B.arr, B.over, cap

    def getb(i):
        return arr[i] if i < CAP else over.get(i, 0)

    def setb(i, v):
        if i < CAP:
            arr[i] = v
        else:
            over[i] = v

    def path(n):
        pts = [n]
        x = n
        for _ in range(L):
            i = (x + 1) >> 1
            b = arr[i] if i < CAP else over.get(i, 0)
            if x & 1:
                x = (x - 1) >> 1 if b == 2 else (3 * x + 1) >> 1
            else:
                x = (3 * x) >> 1 if b == 2 else x >> 1
            if x < n:
                return pts, True
            pts.append(x)
        return pts, False

    def freeze(pts):
        for x in pts:
            i = (x + 1) >> 1
            if (arr[i] if i < CAP else over.get(i, 0)) == 0:
                setb(i, 1)

    def try_flip(n, v):
        i = (v + 1) >> 1
        if getb(i) != 0:
            return False
        setb(i, 2)
        p2, ok2 = path(n)
        if ok2:
            freeze(p2)
            return True
        setb(i, 0)
        return False

    def try_flip_safe(n, v):
        """flip the pair of the odd point v = 2i-1 only if n then descends AND the partner 2i = v+1 (pushed up
        to 3i) descends below itself under the current bits; freeze both certificates."""
        i = (v + 1) >> 1
        if getb(i) != 0:
            return False
        setb(i, 2)
        p2, ok2 = path(n)
        if ok2:
            p3, ok3 = path(v + 1)
            if ok3:
                freeze(p2)
                freeze(p3)
                return True
        setb(i, 0)
        return False

    def T(x):
        return x >> 1 if not (x & 1) else (3 * x + 1) >> 1

    st = {"rescued": 0, "single": 0, "multi": 0, "stuck": 0, "flips": 0, "at_own": 0, "rank_sum": 0}
    for n in range(3, N + 1):
        pts, ok = path(n)
        if ok:
            freeze(pts)
            continue
        st["rescued"] += 1
        cand = pts[:L]
        done = False
        if mode == "afb":
            def optA():
                return (n & 1) == 1 and try_flip(n, n)

            def optF():
                for k in range(1, len(pts)):
                    v, w = pts[k - 1], pts[k]
                    if (v & 1) and getb((v + 1) >> 1) != 2 and w == (3 * v + 1) >> 1 and w % 4 == 3 and w < 2 * n:
                        return try_flip(n, w)
                return False

            def optB():
                t = T(n)
                return (t & 1) == 1 and try_flip(n, t)

            order = (optA, optF, optB) if n % 8 == 3 else (optF, optB, optA)
            for fn in order:
                if fn():
                    done = True
                    break
            if done:
                st["single"] += 1
                st["flips"] += 1
        else:
            assert mode == "peak"
            # partner isolation (as in THM-4475): only odd points v = 2 (mod 3), so the partner v+1 = 0 (mod 3)
            # has no odd Collatz preimage and is entered only from its own halving chain
            odd = sorted({v for v in cand if (v & 1) and v % 3 == 2}, reverse=True)
            for rank, v in enumerate(odd):
                if try_flip_safe(n, v):
                    done = True
                    st["single"] += 1
                    st["flips"] += 1
                    st["rank_sum"] += rank
                    if v == n:
                        st["at_own"] += 1
                    break
            if not done:
                # fall back to THM-4475's source-level options A, F, B
                t = T(n)
                for v in ([n] if (n & 1) else []) + [w for w in pts[1:] if (w & 1) and w % 4 == 3 and w < 2 * n][:1] + ([t] if t & 1 else []):
                    if try_flip(n, v):
                        done = True
                        st["single"] += 1
                        st["flips"] += 1
                        st["fallback_afb"] = st.get("fallback_afb", 0) + 1
                        break
        if done:
            continue
        # fallback (rare; happens for partners 2i of earlier flips, whose own pair is already flipped):
        # depth-limited search over flips at free odd points of the current path, earliest first
        chosen = []

        def dfs(depth):
            p2, ok2 = path(n)
            if ok2:
                return True
            if depth == 0:
                return False
            seen = set()
            for v in p2[:L]:
                if not (v & 1) or v in seen:
                    continue
                seen.add(v)
                i = (v + 1) >> 1
                if getb(i) != 0:
                    continue
                setb(i, 2)
                chosen.append(i)
                if dfs(depth - 1):
                    return True
                chosen.pop()
                setb(i, 0)
            return False

        if dfs(4):
            p2, ok2 = path(n)
            freeze(p2)
            st["multi"] += 1
            st["flips"] += len(chosen)
        else:
            st["stuck"] += 1
    fails = 0
    for n in range(3, N + 1):
        if not path(n)[1]:
            fails += 1
    h, q4 = N // 2, N // 4
    nf2 = sum(1 for i in range(1, h + 1) if getb(i) == 2)
    nf4 = sum(1 for i in range(1, q4 + 1) if getb(i) == 2)
    st.update({"fails": fails, "dens_half": nf2 / h, "dens_quarter": nf4 / q4, "overflow_entries": len(over)})
    return st

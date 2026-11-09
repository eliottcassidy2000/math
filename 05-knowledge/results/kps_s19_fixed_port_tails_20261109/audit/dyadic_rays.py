"""S19 audit: THM-5085 (3) pure-dyadic statements, base 2, cores of order <= 7.
(a) single-bridge sink rays with power-of-2 leading coefficient that land just below a power of two
    (value h0 + 2^s R0 = 2^K - c, s >= 1, R0 a power of 2, small c) - tests the README sentence
    'the values just below powers of 2 always need non-dyadic leading parts';
(b) bridge-appending rays: after a word ending in B the state is (h, R, -hu) and appending B^m gives
    h_m = [h - (R - hu)] + 2^m (R - hu); look for pure-dyadic R - hu with offset h - (R - hu) <= 0 or in {7, 21}.
"""
import sys
cores = [tuple(map(int, l.split())) for l in open(sys.argv[1])]
def is_pow2(v): return v > 0 and v & (v - 1) == 0
ex = []
for h0, R0, D0, hu in cores:
    if not is_pow2(R0): continue
    for s in range(1, 25):
        v = h0 + (1 << s) * R0
        K = v.bit_length()
        c = (1 << K) - v
        if 1 <= c <= 3 and v < 10**6:
            ex.append((v, c, h0, R0, D0, hu, s))
ex.sort()
print("(a) single-bridge pure-dyadic sink-ray values of the form 2^K - c, c <= 3 (value, c, h0, R0, D0, hu, s):")
for e in ex[:12]: print("   ", e)
print("    total found below 10^6:", len(ex))
# (b) bridge-appending rays from states reached by words of length <= 6 ending in B
found = []
for h0, R0, D0, hu in cores:
    from itertools import product
    for L in range(1, 7):
        for w in product('SB', repeat=L):
            if w[-1] != 'B': continue
            h, R, D = h0, R0, D0
            for c in w:
                if c == 'S': R, D = 2 * R, 0
                else: h, R, D = h + (R + D), 2 * R - hu, -hu   # x = 2: B adds x(R+D)/2 = R + D
            lead = R - hu
            if is_pow2(lead):
                off = h - lead
                if off <= 0 or off in (7, 21):
                    found.append((off, lead, h0, R0, D0, hu, ''.join(w)))
found = sorted(set(found))
print("(b) bridge-appending rays h_m = off + 2^m * lead with lead a power of 2 and off <= 0 or off in {7, 21}:")
for f in found[:12]: print("   ", f)
print("    total (core, word) pairs found:", len(found))
# (c) PROOF.md Cor. (c): sink-ray offsets o(I) = h0 - hu Q(I) (x = 2) equal to 7 or 21, with |I| = 2, I = {1, 1+d}
def Q2(d): return 2 if d == 1 else 2 ** (d - 1)
hits = []
for h0, R0, D0, hu in cores:
    for d in range(1, 12):
        o = h0 - hu * Q2(d)
        if o in (7, 21): hits.append((o, h0, R0, D0, hu, d))
hits.sort()
print("(c) sink rays with offset 7 or 21 (offset, h0, R0, D0, hu, d for I = {1, 1+d}); values o + 2^s R0 (1 + 2^d):")
for t in hits[:6]: print("   ", t)
print("    total:", len(hits))

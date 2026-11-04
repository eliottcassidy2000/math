"""NOTE-level probe of the packet's Section 6 'Sharpness' remark: the prime-eight quotient Q8 is coarse-blind
mod 4 only for the NONCONSTANT coarse part (equivalently in S188's equal-profile setting).  With arbitrary
same-size target interiors (this packet's setting) the constant term H(W)-H(W') can be 2 mod 4."""
from aud_lib import compose, converse, ham_start_end, run_jobs, transitive_from_order
Q8 = [238, 4, 240, 118, 195, 146, 34, 74]
S, E = ham_start_end(Q8)
print('Q8 S =', S, ' E =', E, ' p =', [(s + e) % 2 for s, e in zip(S, E)])
TT2 = transitive_from_order([0, 1], 2)
TT3 = transitive_from_order([0, 1, 2], 3)
C3 = [2, 4, 1]
for label, src0, tgt0 in (('equal interiors (TT3 / TT3)', TT3, TT3), ('TT3 source / C3 target', TT3, C3)):
    W, _ = compose(Q8, [src0] + [TT2] * 7)
    Wp, _ = compose(converse(Q8), [tgt0] + [TT2] * 7)
    r = run_jobs([('H', len(W), W, None, [(0, 0)]), ('H', len(Wp), Wp, None, [(0, 0)])])
    print('%-28s H(W)=%d H(W\')=%d  difference mod 4 = %d' % (label, r[0][0], r[1][0], (r[0][0] - r[1][0]) % 4))

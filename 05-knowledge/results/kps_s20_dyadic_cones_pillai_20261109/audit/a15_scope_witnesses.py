# a15: explicit witnesses that dyadic cone steps can re-land on a dyadic port OUTSIDE the attached gadget
# (so THM-5109 (E2)/(E3) need the qualifier "new port inside the attached gadget" everywhere they are summarised).
# Each witness: base (T, v) with R_v(T) = 2^j, cone step = one-arc join with TT_m at entry p >= 2 (value h + 2^(p-1) R),
# and a vertex w of T (not in the gadget) with R_w(N) a power of two.  All counts by the audit's own a_paths.exe.
from a_lib import *
def show(name, T, v, m, p):
    rT = run([T])[0]; h, e, R = rT['H'], rT['e'][v], rT['R'][v]
    assert ispow2(R)
    N = one_arc_join(T, v, TT(m), p - 1); r = run([N])[0]
    assert r['H'] == h + 2 ** (p - 1) * R
    hits = [(w, r['e'][w], r['R'][w]) for w in range(len(T)) if ispow2(r['R'][w])]
    rows = [sum(1 << j for j in range(len(N)) if N[i][j]) for i in range(len(N))]
    print('%s: base triple %s at v=%d (|T|=%d); TT_%d entry p=%d -> H(N)=%d = %d + 2^%d; old-vertex dyadic ports (w, e_w, R_w): %s; '
          'N out-masks %s' % (name, (h, e, R), v, len(T), m, p, r['H'], h, (p - 1) + R.bit_length() - 1, hits, rows))
# (1) 1 -> 5: T = TT_2 rooted at its sink (triple (1,1,1)); the old source of T becomes a dyadic port of value 5
show('W1', TT(2), 1, 3, 3)
# (2) (1,0,4) -> 9: T = TT_3 rooted at its source; the dyadic port lands on an old vertex
show('W2', TT(3), 0, 2, 2)
# (3) (33,0,128) -> 289: base from the audit census (first realisation of the triple in a03_triples_le8.txt)
for l in open('a03_triples_le8.txt'):
    n, c, v, h, e, R = map(int, l.split())
    if (h, e, R) == (33, 0, 128):
        show('W3', from_code(c, n), v, 2, 2); break

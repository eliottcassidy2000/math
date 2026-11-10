# a14: Session236 offsets/PROOF.md section 2 (audited 2026-11-09 16:19, before S20) already gives dyadic-port rays:
# "T rows [12,29,24,16,1], port 1: H=5, R=32 ... H(J_m)=2^(m+5)+5"; "T rows [14,28,24,16,1], port 2: H=9, R=16 ...
# H(J_m)=2^(m+4)+9".  Re-verify the seeds' triples and the rays with the audit's own counter (rows = out-neighbour masks).
from a_lib import *
def from_rows(rows):
    n = len(rows); A = empty(n)
    for i, m in enumerate(rows):
        for j in range(n):
            if m >> j & 1: A[i][j] = 1
    for i in range(n):
        for j in range(i + 1, n): assert A[i][j] + A[j][i] == 1
    return A
for rows, port in (([12, 29, 24, 16, 1], 1), ([14, 28, 24, 16, 1], 2)):
    T = from_rows(rows); r = run([T])[0]
    h, e, R = r['H'], r['e'][port], r['R'][port]
    rays = []
    for m in range(1, 5):
        # J_m: m sinks appended (fixed port), then a singleton B with the single arc B -> port reversed (one-arc join)
        Tm = ordinal_sum(T, TT(m))
        N = one_arc_join(Tm, port, TT(1), 0)
        rays.append(run([N])[0]['H'])
    print('Session236 seed rows', rows, 'port', port, ': (h,e,R) =', (h, e, R), '; H(J_m), m=1..4 =', rays,
          '; equals h + 2^m R:', rays == [h + 2 ** m * R for m in range(1, 5)])

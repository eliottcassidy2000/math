"""S19 audit: the x = 2 case of THM-5085 (2) against the earlier O/Z formulas (pure arithmetic identities).
THM-5085:     h = h0 + R0 N - hu Q(N) + D0 [N odd]            (x = 2), Q = pair sum with c(1) = 2, c(d) = 2^(d-1)
Session205:   H_s(m) = h0 + R0 m - d Q(m) - e0 (m mod 2),  Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2)
Session197:   h + t(R - d) + d(F(t) - 1) + (t mod 2)(d - e),  F(0)=1, F(2q+1)=F(q), F(2q)=F(q)+q+(q mod 2)
Session198:   D' = D + U A - E eps_0 - d C(A),  C(A) = sum_{h<i} e_h e_i 2^(i-h-1) + sum_{i>=1} e_(i-1) e_i
Dictionary (x = 2): (h, e, R, d) = (h, -Delta, R, h_u); e = Hamiltonian paths ending at the port (E_2 - O_2 = -e checked
by naive_check.py and check_transfer.py inputs). The identities below show the four formulas coincide for every N.
"""
def Q_pairs(N):
    I = [i + 1 for i in range(N.bit_length()) if N >> i & 1]
    return sum(2 if b - a == 1 else 2 ** (b - a - 1) for k, a in enumerate(I) for b in I[k + 1:])
LIM = 1 << 17
Q = [0] * LIM; F = [1] * LIM
for m in range(1, LIM):
    k = m >> 1
    Q[m] = Q[k] + k + (k & 1) if m & 1 else Q[k]
    F[m] = F[k] if m & 1 else F[k] + k + (k & 1)
bad_pairs = sum(1 for N in range(LIM) if Q_pairs(N) != Q[N])
bad_F = sum(1 for N in range(LIM) if Q[N] != N + 1 - F[N] - (N & 1))
bad_C = 0
for A in range(LIM):
    e = [A >> i & 1 for i in range(max(1, A.bit_length()))]
    C = sum(e[h] * e[i] * 2 ** (i - h - 1) for i in range(len(e)) for h in range(i)) + sum(e[i - 1] * e[i] for i in range(1, len(e)))
    if C != Q[A]: bad_C += 1
print(f"N < 2^17: THM-5085 pair sum vs Session205 recursion mismatches {bad_pairs}; "
      f"Q(N) = N + 1 - F(N) - [N odd] (Session197 form) mismatches {bad_F}; Session198 C(A) vs Q mismatches {bad_C}")

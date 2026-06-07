#!/usr/bin/env python3
"""
THM-438 ADDENDUM-10 (monad-explorer-2026-06-07, 12th session).

Two durable findings about the cycle-rank triangle t(k,m) of even-series patterns
(S_k = sum_m (-1)^m t(k,m) = (-1)^k C_k):

 (A) THE DIAGONAL t(k,k) = A088368(k) IS A NAMED, CLOSED-FORM OBJECT:
     A088368 = "number of partitions of [n] into sets of NONCROSSING LISTS" (Callan,
     arXiv:0711.4841).  G.f. A(x) satisfies A(x/F(x)) = F(x), F(x)=sum_{n>=0} n! x^n.
     Established asymptotic (Kotesovec, OEIS):  a(n) ~ e * n!.

 (B) CORRECTION OF ADD-9 point (6).  ADD-9 flagged "A088368(m) ~ e*m!" as FALSE and
     proposed "A088368(m) ~ m!(m+2)/2 instead", from the ratios at m<=7.  That is the
     ERROR: the ratio a(m)/m! OVERSHOOTS e, peaks at m=8 (~4.36), then DESCENDS slowly
     back toward e.  ADD-9 only sampled the rising side (m<=7).  Also a transcription
     slip: A088368(7) = 21477, not "22417".  See MISTAKE-063.

 (C) The Catalan collapse runs between TWO CATALOGUED ENDPOINTS; the path is structureless:
     diagonal  t(k,k)=A088368(k) ~ e*k!  (wild, A088368) -->  signed row sum (-1)^k C_k
     (tame, Catalan).  Every intermediate sequence is OEIS-NEGATIVE (recorded so not retried).

 (D) Proven row-7 entries (columns m=1..4 are fully determined by their proved column GFs;
     diagonal from A088368), plus the signed-sum constraint linking the two unknowns.
"""
from math import comb, factorial as fact, e

# ---- the known triangle (k<=6), from paley_starstar_triangle_*_monad ----
T = {1:[1], 2:[1,3], 3:[1,9,13], 4:[1,18,72,69],
     5:[1,30,230,580,421], 6:[1,45,560,2626,4845,2867]}

A088368 = [1,1,3,13,69,421,2867,21477,175769,1567273,15213955,160727997,
           1846282381,23013527421,310284575683,4506744095141,70199956070705,
           1167389338452753,20636801363971139,386304535988493101,7630926750477398037]

print("="*72)
print("(A)+(C) diagonal = A088368 ; signed row sum = (-1)^k C_k  [both CATALOGUED]")
print("="*72)
for k,row in T.items():
    signed = sum((-1)**m * row[m-1] for m in range(1,k+1))
    Ck = comb(2*k,k)//(k+1)
    print(f"  k={k}: diag t(k,k)={row[-1]:>6} =A088368({k})={A088368[k]:>6} | "
          f"signed={signed:>5}  (-1)^k C_k={(-1)**k*Ck:>5}  match={signed==(-1)**k*Ck}")

print()
print("="*72)
print("(B) ADD-9 pt(6) REFUTED: a(n)/n! overshoots e, peaks at n=8, descends to e")
print("="*72)
print(f"  e = {e:.6f}")
peak_n, peak_v = 0, 0.0
for n in range(2,len(A088368)):
    r = A088368[n]/fact(n)
    if r>peak_v: peak_v, peak_n = r, n
    add9_pred = (n+2)/2.0   # ADD-9's wrong "m!(m+2)/2" => ratio (m+2)/2
    print(f"  n={n:2d}  a(n)/n!={r:7.4f}   [ADD-9 (n+2)/2={add9_pred:5.1f}]")
print(f"  --> ratio PEAKS at n={peak_n} ({peak_v:.4f}); then strictly DECREASES toward e.")
print(f"  --> ADD-9's '(n+2)/2' fits only the rising side n<=7 and diverges after the peak.")

print()
print("="*72)
print("(D) PROVEN ROW 7 (columns m=1..4 from proved column GFs; diag from A088368)")
print("="*72)
# proved falling-factorial column forms t(k,m)=(k)_m h_m(k), VERIFIED m<=4
def falling(k,m):
    p=1
    for j in range(m): p*=(k-j)
    return p
h2 = lambda k: 3/2
h3 = lambda k: (5*k-2)/6
h4 = lambda k: (181*k*k-219*k+50)/720
def t(k,m):
    if m==1: return 1
    if m==2: return round(falling(k,2)*h2(k))
    if m==3: return round(falling(k,3)*h3(k))
    if m==4: return round(falling(k,4)*h4(k))
    return None
# sanity vs known triangle
ok = all(t(k,m)==T[k][m-1] for k in T for m in range(1,min(k,4)+1))
print(f"  column forms reproduce known triangle (m<=4): {ok}")
row7 = {m:t(7,m) for m in (1,2,3,4)}
row7[7] = A088368[7]
print(f"  t(7,1..4) = {row7[1]}, {row7[2]}, {row7[3]}, {row7[4]}   t(7,7)=A088368(7)={row7[7]}")
# signed-sum constraint:  sum_m (-1)^m t(7,m) = (-1)^7 C_7 = -429
C7 = comb(14,7)//8
known = -row7[1]+row7[2]-row7[3]+row7[4]-row7[7]   # m=1,2,3,4,7 (m=5 minus, m=6 plus)
# full: -t1+t2-t3+t4 -t5 +t6 -t7 = -429  => t6 - t5 = -429 - (known)
gap = -C7 - known
print(f"  C_7={C7};  signed-sum identity =>  t(7,6) - t(7,5) = {gap}")
print(f"  (one linear relation between the two still-uncomputed entries t(7,5),t(7,6))")

print()
print("="*72)
print("(E) the bridge polynomial h_m(k): two ends -> {e, 0}")
print("="*72)
print("  WILD end  h_m(m)   = A088368(m)/m!  -> e        (Kotesovec)")
print("  TAME end  h_m(-1)  = (2^m-1)/((-1)^m m!) -> 0    (super-exponential)")
for m in range(2,9):
    wild = A088368[m]/fact(m)
    tame = (2**m-1)/((-1)**m * fact(m))
    print(f"  m={m}: h_m(m)={wild:7.4f}   h_m(-1)={tame:+.5f}")
print("  => h_m interpolates the bridge from 0 (tame/Mersenne) to e (wild/noncrossing-lists).")

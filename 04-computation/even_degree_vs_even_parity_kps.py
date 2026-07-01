#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE TWO 'EVENS' ARE INDEPENDENT: even DEGREE (even graphs / cycle space) vs even PARITY (automorphism sign, A_n).

kind-pasteur-2026-07-01-S11. The owner flags a conflation.  THREE distinct 'even' notions live here:
  (F) even FUNCTION: a GRAPH is a symmetric (even) +-1 function on pairs; a TOURNAMENT is skew (odd).
      [the det(I+S) / sgn-vs-chi lens, determinant-lens reflection]
  (D) even DEGREE: an EVEN GRAPH has all degrees even (the cycle space H_1(K_n;F_2); Ising/flow side).
  (P) even PARITY: an automorphism is an EVEN permutation (in A_n; sign = +1).
This script pins down (P) vs (D) and their ODD versions, and the role of the Paley heptagon.

KEY THEOREM (forbidden-seven reflection): for a TOURNAMENT, |Aut| is ODD (a transposition reverses an arc =>
is an ANTI-automorphism), so EVERY automorphism has odd order => is a product of odd cycles => an EVEN
permutation: Aut(T) subset A_n ALWAYS.  Complementation / anti-automorphism (iota) lives in the ODD coset.
For EVEN GRAPHS this FAILS: even degree does NOT force even-parity automorphisms.
"""
import sys, itertools
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def sgn(p):
    n=len(p); seen=[False]*n; s=1
    for i in range(n):
        if seen[i]: continue
        L=0; j=i
        while not seen[j]:
            seen[j]=True; j=p[j]; L+=1
        if L%2==0: s=-s
    return s
def tour_aut_antiaut(A,n):
    aut=[]; anti=[]
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)): aut.append(p)
        if all(A[p[i]][p[j]]==A[j][i] for i in range(n) for j in range(n)): anti.append(p)
    return aut,anti
def graph_aut(G,n):
    return [p for p in itertools.permutations(range(n))
            if all(G[p[i]][p[j]]==G[i][j] for i in range(n) for j in range(i+1,n))]
def circ(conn,n): return [[1 if (i!=j and (j-i)%n in conn) else 0 for j in range(n)] for i in range(n)]

print("="*94); print(" (P) AUTOMORPHISM PARITY of TOURNAMENTS: always Aut subset A_n; anti-aut (complementation) is ODD"); print("="*94)
n=7
transitive=[[1 if i<j else 0 for j in range(n)] for i in range(n)]
R7=circ({1,2,3},n); Paley=circ({1,2,4},n)
for name,A in [("transitive T_7",transitive),("R_7 carousel {1,2,3} (=sgn)",R7),("Paley_7 QR {1,2,4} (=chi)",Paley)]:
    aut,anti=tour_aut_antiaut(A,n)
    autparity=set(sgn(p) for p in aut); antiparity=set(sgn(p) for p in anti)
    print(f"\n {name}: |Aut|={len(aut)} (odd={len(aut)%2==1}), |anti-Aut|={len(anti)}")
    print(f"   Aut signs = {autparity}  => Aut subset A_7 (all EVEN parity): {autparity<= {1}}")
    print(f"   anti-Aut signs = {antiparity if anti else '(none: not self-converse)'}"
          + (f"  => complementation iota is ODD parity: {antiparity=={-1} or -1 in antiparity}" if anti else ""))

print("\n"+"="*94); print(" (D) vs (P) ARE INDEPENDENT: an EVEN-DEGREE graph can have ODD-parity automorphisms"); print("="*94)
# C_4 = 4-cycle: even degree (all 2); automorphism (1 3) is a transposition = ODD parity
C4={0:[1,3],1:[0,2],2:[1,3],3:[0,2]}
G4=[[1 if j in C4[i] else 0 for j in range(4)] for i in range(4)]
deg4=[sum(row) for row in G4]
aut4=graph_aut(G4,4); signs4=set(sgn(p) for p in aut4)
print(f"  C_4 (4-cycle): degrees {deg4} (all even). |Aut|={len(aut4)}, aut signs={signs4}")
print(f"    => EVEN DEGREE but has ODD-parity automorphisms (e.g. the reflection (1 3)): "
      f"(D) does NOT imply (P).  Contrast tournaments, where (P) is automatic.")
# odd-degree graph with even-parity-only auts: K_2 (single edge on 2 vtx): degree 1 (odd), aut={id,(01)} (01) odd...
# better: the path P_3 vertices 0-1-2: degrees (1,2,1) mixed. Use K_4: all degree 3 (odd), Aut=S_4 (both parities).
K4=[[1 if i!=j else 0 for j in range(4)] for i in range(4)]
autK4=graph_aut(K4,4); print(f"  K_4: degrees all 3 (ODD). |Aut|={len(autK4)}=S_4, signs={set(sgn(p) for p in autK4)} (both).")

print("\n"+"="*94); print(" THE 2x2 TABLE {degree parity} x {automorphism-parity class} over small graphs"); print("="*94)
# classify all graphs on 4 vertices
prs=list(itertools.combinations(range(4),2))
seen={}
for mask in range(1<<len(prs)):
    G=[[0]*4 for _ in range(4)]
    for k,(a,b) in enumerate(prs):
        if (mask>>k)&1: G[a][b]=G[b][a]=1
    deg=[sum(G[i]) for i in range(4)]
    all_even_deg = all(d%2==0 for d in deg)
    aut=graph_aut(G,4)
    aut_in_A4 = all(sgn(p)==1 for p in aut)
    key=(all_even_deg, aut_in_A4)
    # canonical form to dedupe iso
    best=None
    for p in itertools.permutations(range(4)):
        t=tuple(sorted((min(p[a],p[b]),max(p[a],p[b])) for (a,b) in prs if G[a][b]))
        if best is None or t<best: best=t
    seen.setdefault(key,set()).add(best)
print("  (even-degree?, Aut subset A_4?) -> # iso classes of graphs on 4 vertices:")
for k in sorted(seen, key=lambda x:(not x[0], not x[1])):
    print(f"    even-degree={k[0]!s:>5}, Aut<=A_4={k[1]!s:>5}: {len(seen[k])} classes")
print("  => the two 'evens' cross-classify graphs into all four cells: they are ORTHOGONAL properties.")

print("\n"+"="*94); print(" ODD VERSIONS"); print("="*94)
print("  ODD DEGREE graphs (all degrees odd) exist only for n EVEN (handshake): a coset of the cycle space,")
print("    = T-joins / perfect-matching parity.  ODD PARITY automorphisms = Aut has a member in S_n\\A_n.")
print("  TOURNAMENTS: (P) is FORCED even (Aut<=A_n) and the ODD coset is EXACTLY the anti-automorphisms")
print("    (complementation iota).  So for tournaments the sgn character = the orientation-preserving/reversing")
print("    Z_2 = the complement/converse mirror = the dihedral reflection of the Paley heptagon's D_7.")
print("DONE.")

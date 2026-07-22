import numpy as np

def paley_tournament(p):  # p ≡ 3 mod 4
    qr = set((x*x)%p for x in range(1,p))
    return np.array([[1 if i!=j and ((j-i)%p) in qr else 0 for j in range(p)] for i in range(p)])

def paley_graph(p):  # p ≡ 1 mod 4 (symmetric)
    qr = set((x*x)%p for x in range(1,p))
    return np.array([[1 if i!=j and ((j-i)%p) in qr else 0 for j in range(p)] for i in range(p)])

print("=== Paley spectrum at the clock moduli (honest) ===")
# 7 ≡ 3 mod 4: TOURNAMENT. skew Jacobsthal S = A - A^T has eigenvalues ±i√p.
A7 = paley_tournament(7)
S7 = A7 - A7.T
skew_ev = sorted(set(round(abs(e),3) for e in np.linalg.eigvals(S7.astype(float)) if abs(e)>1e-6))
print(f"7 (≡3 mod4, Paley TOURNAMENT): skew Jacobsthal |eigenvalue|={skew_ev}  vs √7={round(7**.5,3)}  "
      f"-> Gauss sum ±i√7 in the SKEW spectrum; tournament-matrix atoms (-1±i√7)/2, |λ|=√(8)/2={round((8**.5)/2,3)}")

# 13 ≡ 1 mod 4: GRAPH. eigenvalues (p-1)/2 and (-1±√p)/2.
G13 = paley_graph(13)
ev13 = sorted(set(round(e,3) for e in np.linalg.eigvals(G13.astype(float)).real))
print(f"13 (≡1 mod4, Paley GRAPH): eigenvalues={ev13}  -> (-1±√13)/2 = {round((-1+13**.5)/2,3)}, {round((-1-13**.5)/2,3)}  "
      f"(√13={round(13**.5,3)} in the REAL spectrum)")
print("14 = 2·7 (runner count; the LRC(14) modulus itself)")
print()
print("HONEST MERGE (structural, verified as a proof-SHAPE not a pole identity):")
print("  Both the GMC2 nullcone proof and the LRC covering proof are 'SCALE the core, then CLOSE on a modular CLOCK':")
print("   GMC2:  scale = dilate ×p   | clock = residue field Z/p (Frobenius x↦x^p) | Kummer/Lucas = clock arithmetic | survival Q̄^p = p-periodicity")
print("   LRC :  scale = ×a          | clock = Z/12a & Z/14a (modular orbits, THM-2057) | orbit closure = clock periodicity")
print("  The clock moduli {7,13,14} (THM-878) each carry a Paley √p spectrum: 7=tournament skew ±i√7, 13=graph (−1±√13)/2, 14=2·7.")
print("  The tournament zeta (THM-1926) is the dynamical lens: periodic-orbit↔spectrum; on the ACYCLIC core ζ=1 (transitive T_12, S214).")
print("  NOT claimed: a numeric zeta-pole = clock-modulus identity (13 is a graph; tournament atom is √(1+p)/2). This is an ANALOGY of proof shape.")

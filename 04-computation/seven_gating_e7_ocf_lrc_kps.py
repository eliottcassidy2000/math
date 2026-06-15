"""
seven_gating_e7_ocf_lrc_kps.py
Rigorous separation of genuine structure vs coincidence for the prime 7 (and 21)
across E7, the OCF forbidden values, and LRC(14).

We do NOT prove new theorems here. We CHECK arithmetic/combinatorial claims so the
verdict catalog is grounded in computation, not assertion.
"""

from math import comb
from itertools import combinations

print("="*70)
print("PART 1: The 7&21 catalog -- arithmetic checks")
print("="*70)

# E7 published invariants (sourced: Wikipedia E7, reflection groups degree tables)
dim_E7 = 133
rank_E7 = 7
roots_E7 = 126
posroots_E7 = 63
exponents_E7 = [1,5,7,9,11,13,17]
degrees_E7 = [2,6,8,10,12,14,18]
coxeter_h = 18
weyl_order_E7 = 2903040
minuscule_56 = 56

print(f"dim E7        = {dim_E7}      7*19 = {7*19}   -> 7 | dim?  {dim_E7%7==0}")
print(f"rank E7       = {rank_E7}")
print(f"roots         = {roots_E7}    = 18*7 = {18*7}  -> 7|126?   {roots_E7%7==0}")
print(f"pos roots     = {posroots_E7}     = 9*7 = {9*7}   -> 7|63?    {posroots_E7%7==0}")
print(f"exponents     = {exponents_E7}  -> is 7 an exponent? {7 in exponents_E7}")
print(f"  middle exponent (h/2 symmetry): exponents sum to {sum(exponents_E7)} ; "
      f"pairs m+m'=h=18? {[ (e, coxeter_h-e) for e in exponents_E7]}")
print(f"degrees       = {degrees_E7}  -> degree 8 = 7+1 present? {8 in degrees_E7}")
print(f"Coxeter h     = {coxeter_h}    (NOT divisible by 7: {coxeter_h%7==0})")
print(f"|Weyl(E7)|    = {weyl_order_E7} ; divisible by 7? {weyl_order_E7%7==0}; "
      f"by 7^2? {weyl_order_E7%49==0}")
# factor weyl order
def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
print(f"  |Weyl(E7)| factorization = {factor(weyl_order_E7)}  (expected 2^10*3^4*5*7)")
print(f"minuscule rep = {minuscule_56} = 2*28 = {2*28} ; = 7*8 = {7*8}  -> BOTH hold: "
      f"56 = 2*28 (antipodal pairs) AND 56 = 7*8 (Fano lines x tripartite)")

print()
print("--- Where does 21 genuinely appear? ---")
print(f"C(7,2) = {comb(7,2)} = dim so(7) (antisymmetric 7x7).  21 | ?")
print(f"  Under G2: 21 = 14 + 7  ->  {14+7==21}  (dim G2=14, vector 7). SOURCED (Baez).")
print(f"C(7,3) = {comb(7,3)} = 35 = dim of 3-forms on R^7 (the G2 associative 3-form lives here).")
print(f"The binomial tower on the Fano 7-set: "
      f"C(7,1),C(7,2),C(7,3) = {comb(7,1)},{comb(7,2)},{comb(7,3)} = 7,21,35.")
print(f"  7+21+35 = {7+21+35} ; note 56-rep weight under SL(2)^7? 56 = ? "
      f"  (56 != 7+21+35={7+21+35}; the 7+21+35 is the EXTERIOR algebra of F_2^7 truncated,"
      f" NOT the 56). So 21 in E7 = dim so(7), a SUBALGEBRA fact, distinct from the 56.")

print()
print("--- 21 in tournaments / strong-component side ---")
print(f"non-strong 6-vertex tournament classes = 21 (sourced, A051337-ish).  "
      f"strong = 35, total = 56.  35+21 = {35+21}.")
print(f"  So on the TOURNAMENT side: 56 (T(6)) = 35 strong + 21 non-strong = C(7,3)+C(7,2).")
print(f"  This MIRRORS the E7 binomial tower numerically. Flagged for verdict.")

print()
print("="*70)
print("PART 2: OCF forbidden 7 -- the Phi_3(2) / Fano / K3 facts (arithmetic)")
print("="*70)
def Phi3(x): return x*x+x+1
print(f"Phi_3(2) = 2^2+2+1 = {Phi3(2)}  = 2^3-1 = {2**3-1} = M_3 (Mersenne) = |PG(2,2)| Fano points")
print(f"I(K3, z) = 1 + 3z  (K3 = triangle conflict graph, 3 vertices pairwise adjacent => "
      f"alpha_0=1, alpha_1=3, alpha_2=0)")
print(f"I(K3, 2) = 1 + 3*2 = {1+3*2}  = the forbidden H value 7.  So 7 = I(K3,2) = Phi_3(2). CHECK.")
print(f"21 = 3 * 7 = {3*7} = 3 * Phi_3(2).  Forbidden because only factorization routes through 7.")
print(f"Note: the '3' in K3 = A2 root system / smallest cyclotomic Phi_3 / the 3-cycle. "
      f"Phi_3(2)=7 ties the SMALLEST cyclotomic to the SMALLEST projective plane.")

print()
print("="*70)
print("PART 3: The 'exceptional cluster gas' -- is the E7 quartic a partition sum gated by 7?")
print("="*70)
# OCF: H = sum_k alpha_k 2^k  (hard-core gas, fugacity 2). degree = clique number of Omega.
# E7 quartic I4: degree-4 SL(2)^7-covariant; Duff-Ferrara: I4 = combination of the 7 Cayley
# hyperdeterminants Det_3 (one per Fano line), each Det_3 is a DEGREE-4 invariant of a 2x2x2 tensor.
# So I4 is a sum/cluster over the 7 Fano lines of degree-4 'local' invariants.
print("OCF partition function:  Z = I(Omega,2) = sum_k alpha_k * 2^k   (hard-core lattice gas).")
print("  - lives at FUGACITY 2; degree in z = independence/clique number of conflict graph Omega.")
print("  - cluster expansion (Mayer): log Z = sum over connected subgraphs (Witt/necklace W_k).")
print()
print("E7 quartic Cartan invariant I4 (Duff-Ferrara 2007, PRD 76 025018):")
print("  - DEGREE 4 in the 56 charges; SL(2)^7-covariant.")
print("  - 56 = 7 lines x 8 (each line = a tripartite 2x2x2 system).")
print("  - I4 = (combination of) the 7 Cayley hyperdeterminants Det(line_i), one per Fano line.")
print("  - each Cayley hyperdeterminant of a 2x2x2 tensor is itself DEGREE 4 (Cayley's = 3-tangle).")
print()
# The KEY DIMENSIONAL COMPARISON: is the OCF degree-2 fugacity 'a partition function' and is the
# E7 quartic also one? The honest test: both are MULTILINEAR invariants summed over a combinatorial
# index set (cliques of Omega // Fano lines). Compare the *structure*, not just the word 'quartic'.
print("STRUCTURAL COMPARISON TABLE")
print(f"{'feature':32} | {'OCF gas Z=I(Omega,2)':30} | {'E7 quartic I4':30}")
print("-"*100)
rows = [
 ("index set",         "cliques/indep sets of Omega",  "7 lines of Fano PG(2,2)"),
 ("local object",      "a k-clique (z^k term)",        "a 2x2x2 tripartite tensor"),
 ("local invariant",   "monomial 2^k",                 "Cayley hyperdet (deg 4)"),
 ("global = sum of",   "alpha_k 2^k over k",           "Det over 7 Fano lines"),
 ("gating prime",      "7=Phi3(2)=I(K3,2) FORBIDDEN",  "7 = #lines = rank = #qubits"),
 ("degree",            "= clique number (varies)",     "fixed 4 (quartic)"),
 ("role of 7",         "forbidden VALUE of Z",         "STRUCTURAL count of summands"),
]
for a,b,c in rows:
    print(f"{a:32} | {b:30} | {c:30}")
print()
print("CRITICAL DISTINCTION:")
print(" - In OCF, 7 is a FORBIDDEN OUTPUT VALUE of the partition function (Z never equals 7).")
print(" - In E7, 7 is the NUMBER OF TERMS/SUMMANDS (Fano lines) building the invariant.")
print(" - These are DIFFERENT ROLES for 7. Same numeral, different semantic slot.")
print(" - BUT both 7's descend from the SAME object: PG(2,2). OCF's 7 = I(K3,2)=|PG(2,2)| as a")
print("   VALUE; E7's 7 = the 7 lines of PG(2,2) as an INDEX SET. The Fano plane is the common root.")

print()
print("="*70)
print("PART 4: LRC(14) = 2*7 -- where is the 7?")
print("="*70)
print("LRC(14): modulus 14 = 2*7. s(t)=sin(pi t/7)/(pi t) vanishes IFF 7|t (the 7-vanishing).")
print("period of sign(s) = 14 = 2*7. 7 = the half-period (the 'forbidden' residue class).")
print("So LRC's 7 = the modulus's odd part = where the singular-series kernel VANISHES.")
print("This is a THIRD distinct role: 7 as a MODULUS / vanishing locus, not a value, not a count.")
print()
print("THREE ROLES OF 7, made explicit:")
print("  OCF:  7 = a forbidden VALUE of H  (7 = I(K3,2) = Phi_3(2))")
print("  E7:   7 = the COUNT of Fano lines / rank / qubits (structural index)")
print("  LRC:  7 = the odd MODULUS / vanishing locus (s(t)=0 iff 7|t)")
print("Common ancestor candidate: 7 = 2^3-1 = Phi_3(2) = |PG(2,2)|. The Fano plane is the only")
print("object that appears LITERALLY in two of three (OCF value via I(K3,2)=|PG(2,2)|; E7 via")
print("octonion/qubit Fano lines). LRC's 7 is just the odd part of 14 -- no Fano structure shown.")

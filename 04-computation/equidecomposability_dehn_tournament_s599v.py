"""EQUINUMEROSITY vs EQUIDECOMPOSABILITY for tournaments (and the Cl2(π/3) scissors-congruence
constant). H(T) = the EQUINUMEROSITY invariant (a count = bijection class). The STRONG-COMPONENT
MULTISET = the EQUIDECOMPOSABILITY invariant (scissors-congruence class: cut into strong pieces,
reassemble). Same H, different piece-multiset = EQUINUMEROUS but NOT equidecomposable = the
'Dehn-distinct' case. Forbidden {7,21} = NO equidecomp class. Cl2(π/3) = ideal-reg-tetrahedron
volume = the log-sin/Mahler scissors invariant at the 3-cycle angle 2π/3. opus-2026-06-03-S599v."""
import math
from itertools import combinations
def Cl2(theta,N=300000): return sum(math.sin(k*theta)/(k*k) for k in range(1,N+1))
def Lob(theta,N=300000): return 0.5*Cl2(2*theta,N)  # Lobachevsky Λ(θ)=½Cl2(2θ)
def main():
    th=math.pi/3
    cl=Cl2(th)
    # log-sin / Mahler form: Cl2(θ) = -∫_0^θ log|2 sin(t/2)| dt
    M=400000; integ=-sum(math.log(abs(2*math.sin((j+0.5)*th/M/2)))*(th/M) for j in range(M))
    print(f"Cl2(π/3): series={cl:.6f};  log-sin integral −∫₀^{{π/3}}log|2sin(t/2)|dt={integ:.6f}  (Mahler/Dehn density)")
    Vtet=3*Lob(th)   # ideal regular tetrahedron volume = 3 Λ(π/3)
    print(f"ideal regular tetrahedron volume = 3·Λ(π/3) = {Vtet:.6f}  (= the maximal hyperbolic simplex; scissors-congruence invariant)")
    print(f"  the 3-cycle (smallest strong tournament) has skew eigenvalues at angle ±2π/3 = the tetrahedron dihedral angle\n")

    # strong H-values by min vertex-size (from s599s)
    strong={3:[3],4:[5],5:[9,11,13,15],6:[15,17,19,23,25,27,29,31,33,37,41,43,45]}
    prims=sorted(set(v for L in strong.values() for v in L))
    # enumerate factorizations of H into strong-values (multisets) = equidecomposability classes
    def facts(H, maxp):
        if H==1: return [()]
        out=[]
        for p in prims:
            if p>maxp or p>H or H%p: continue
            for tail in facts(H//p, p):
                out.append((p,)+tail)
        return out
    print("EQUINUMEROSITY (H) vs EQUIDECOMPOSABILITY (# strong-piece multisets realizing H):")
    print("  H | #equidecomp classes | the strong-piece multisets")
    for H in [1,3,5,7,9,15,21,45,63,189]:
        F=set(facts(H, max(prims)))
        tag = " <-- FORBIDDEN (no equidecomp class)" if len(F)==0 else ("" if len(F)==1 else "  <-- equinumerous but NOT equidecomposable")
        egs=sorted(F)[:5]
        print(f"  {H:3d} | {len(F):2d} | {egs}{tag}")
    print("\n=> H is the COUNT (equinumerosity); the strong-component multiset is the SCISSORS-CONGRUENCE")
    print("   class (equidecomposability); same H w/ multiple multisets = Dehn-distinct; {7,21} have NONE.")
    print("   Cl2(π/3) is the volume invariant of the basic piece (3-cycle/ideal tetrahedron, angle 2π/3).")
if __name__=='__main__': main()

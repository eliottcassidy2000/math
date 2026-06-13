"""
monad-explorer 2026-06-13  —  HYP-2460 / N* probe.

Can a RESONANT product (two triangular-lattice patches glued at a Moser angle w_t)
beat 3N at n in {25,26,27} -- lowering THM-431's ceiling N*<=28?

For each N, every factorization N = a*b, every pair of candidate a- and b-point
triangular-lattice patches, at each Moser rung t in {3,4,7,...} (t Loeschian so the
transverse vectors exist), compute the EXACT Moser-angle unit-distance count
   U = e(A)*b + a*e(B) + Delta_t(A,B).
Report the max and whether it beats 3N. All counts exact (Q(sqrt3,sqrtD)).

NOTE: this is a LOWER-bound search over a curated patch family (not all subsets),
so a non-beat is evidence (consistent with conjecture u(27)=81, N*=28), not proof.
"""
from fractions import Fraction as F
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
from resonant_product_bonus_monad import (
    make_field, ONE, lattice_point, build_product, count_unit_distances,
    edges_in_lattice_graph, eisen_norm, norm_t_displacements, m_alpha,
)

# ---------- candidate triangular-lattice patches of a given size ----------
def hex_region(R):
    r = int((4*R)**0.5)+3
    pts = [(m,n) for m in range(-r,r+1) for n in range(-r,r+1) if eisen_norm(m,n) <= R]
    pts.sort(key=lambda mn: (eisen_norm(*mn), mn))
    return pts

def densest_prefixes(maxsize):
    """greedy: order hex lattice by norm, take prefixes -> dense-ish patches of each size"""
    pool = hex_region(maxsize*2)
    out = {}
    for k in range(1, maxsize+1):
        out.setdefault(k, []).append(pool[:k])
    return out

def parallelograms(maxsize):
    """a x b parallelograms in the triangular lattice (rich in sqrt(t) displacements)"""
    out = {}
    for a in range(1,8):
        for b in range(1,8):
            if a*b>maxsize: continue
            pts = [(i,j) for i in range(a) for j in range(b)]
            out.setdefault(a*b, []).append(pts)
    return out

def triangles(maxsize):
    """triangular blocks T_k = {(i,j): i+j<k} (k rows)"""
    out={}
    for k in range(2,9):
        pts=[(i,j) for i in range(k) for j in range(k-i)]
        if len(pts)<=maxsize:
            out.setdefault(len(pts),[]).append(pts)
    return out

def all_candidates(maxsize):
    cands = {}
    for src in (densest_prefixes, parallelograms, triangles):
        d = src(maxsize)
        for k,lst in d.items():
            cands.setdefault(k, [])
            for p in lst:
                # dedup by frozenset
                if frozenset(p) not in {frozenset(q) for q in cands[k]}:
                    cands[k].append(p)
    return cands

def resonance_bonus(A, B, t):
    Aset=set(A); Bset=set(B); tot=0
    for (am,an) in norm_t_displacements(t):
        mg = sum(1 for (m,n) in A if (m-am,n-an) in Aset)
        mh = sum(1 for (m,n) in B if (m-am,n-an) in Bset)
        tot += mg*mh
    return tot//2

def moser_count_fast(A,B,t):
    eA=edges_in_lattice_graph(A); eB=edges_in_lattice_graph(B)
    generic = eA*len(B)+len(A)*eB
    return generic + resonance_bonus(A,B,t), generic

def factorizations(N):
    return [(a, N//a) for a in range(2, int(N**0.5)+1) if N%a==0]

def loeschian_rungs(tmax):
    return [t for t in range(3, tmax+1) if norm_t_displacements(t)]

if __name__=="__main__":
    MAXSIZE=14
    cands = all_candidates(MAXSIZE)
    rungs = loeschian_rungs(16)   # t with transverse vectors: 3,4,7,9,12,13,...
    print(f"Loeschian rungs (transverse vectors exist): {rungs}")
    print()
    # also pull in the proven small optima as factor candidates via the pure-product
    # K3^[]3 = 81 tie at 27 (THM-433 D): include the Cartesian CUBE explicitly.
    K3=[(0,0),(1,0),(0,1)]
    for N in range(24,31):
        best=( -1, None)
        for (a,b) in factorizations(N):
            for A in cands.get(a,[]):
                for B in cands.get(b,[]):
                    for t in rungs:
                        U,generic = moser_count_fast(A,B,t)
                        if U>best[0]:
                            best=(U,(a,b,t,generic,A,B))
        if best[1] is None:
            print(f"N={N:2d}  3N={3*N:3d}  PRIME -- no two-factor product")
            continue
        U,(a,b,t,generic,A,B)=best
        # pure generic product cap for reference (THM-433): best generic over factorizations
        flag = "  <<< BEATS 3N!" if U>3*N else ("  (ties 3N)" if U==3*N else "")
        print(f"N={N:2d}  3N={3*N:3d}  best resonant-product U={U:3d}  "
              f"[{a}x{b}, rung t={t}, generic={generic}, bonus={U-generic}]{flag}")
        # exact recount of the champion to be safe
        D=4*t-1
        P=build_product(A,B,t,D)
        exactU=count_unit_distances(P,D)
        ndup = len(P)-len({p for p in P})
        assert exactU==U, f"FORMULA MISMATCH at N={N}: formula {U} vs exact {exactU}"
        assert ndup==0, f"coincident points at N={N}: {ndup}"
    print()
    print("Note: pure product cap P(N) (THM-433) ties 3N at 27 (K3^[]3=81) & 30; the")
    print("two-factor resonant search above is over a curated patch family (lower bound).")
    print("All champion counts EXACT-recounted from coordinates; 0 coincident points.")

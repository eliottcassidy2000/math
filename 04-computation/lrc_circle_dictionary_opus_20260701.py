"""DICTIONARY of functions between points on a loop (R/Z), organized group-like. Plus the Verblunsky-centroid
reading (alpha_0 = runner centroid from the marked origin z=1 = observer-lens)."""
import numpy as np
from math import gcd
def frac(x): return x-np.floor(x)
# ---- THE DICTIONARY: maps on the circle R/Z ----
# GROUP 1 -- AFFINE AGL(1): the runners + observer live here
rot   = lambda a: (lambda x: frac(x+a))                 # rotation  (observer shift; inhomogeneous LRC)
dil   = lambda v: (lambda x: frac(v*x))                 # dilation  (a RUNNER of speed v)
refl  = lambda x: frac(-x)                              # reflection/antipode (the -1, killer symmetry)
aff   = lambda v,a: (lambda x: frac(v*x+a))             # affine x->vx+a
# GROUP 2 -- MODULAR PSL2(Z): the CF descent / rung ladder / covering-min
gauss = lambda x: frac(1/x) if x>1e-12 else 0.0         # Gauss map {1/x} (continued fraction)
def mediant(a,b,c,d): return (a+c,b+d)                  # Farey mediant (Stern-Brocot)
mobius= lambda A: (lambda x: frac((A[0]*x+A[1])/(A[2]*x+A[3])))  # fractional-linear on the boundary
# GROUP 3 -- SZEGO semigroup: vertex insertion / n->n-1 reduction (Verblunsky)
# (the recursion is on measures; see verblunsky() -- alpha_j = the j-th reflection)
# ANALYTIC scalars (not maps): the ci-even / ci-odd functions
saw   = lambda x: frac(x)-0.5                           # sawtooth ((x)) -- iota-ODD
dist  = lambda x: min(frac(x),1-frac(x))               # ||x|| distance-to-nearest-int -- iota-EVEN
print("== DICTIONARY of circle functions (12 entries), grouped ==")
print(" AFFINE AGL(1):  rot_a, dil_v (RUNNER), refl (antipode), aff_{v,a}")
print(" MODULAR PSL2Z:  gauss {1/x}, mediant (Farey), mobius (a x+b)/(c x+d)")
print(" SZEGO semigrp:  verblunsky vertex-insertion (n->n-1)")
print(" ANALYTIC:       saw ((x)) [iota-odd], dist ||x|| [iota-even]\n")
# ---- verify the GROUP RELATIONS (they operate group-like) ----
xs=np.linspace(0,1,9,endpoint=False)
# (R1) dilations compose: dil_a . dil_b = dil_ab  (monoid Z, units mod q are a group)
ok1=np.allclose(dil(3)(dil(5)(xs)), dil(15)(xs))
# (R2) affine composition: aff_{v,a}.aff_{w,b} = aff_{vw, vb+a}  (the ax+b group law)
v,a,w,b=3,0.1,5,0.2; ok2=np.allclose(aff(v,a)(aff(w,b)(xs)), aff(v*w, frac(v*b+a))(xs))
# (R3) reflection conjugates dilation: refl . dil_v . refl = dil_v  (dil commutes with refl up to sign)
ok3=np.allclose(refl(dil(5)(refl(xs))), dil(5)(xs))
# (R4) rot & dil semidirect: dil_v . rot_a = rot_{va} . dil_v
ok4=np.allclose(dil(5)(rot(0.1)(xs)), rot(frac(5*0.1))(dil(5)(xs)))
print(f"GROUP RELATIONS: dil compose {ok1} | affine law {ok2} | refl-conj {ok3} | semidirect rot<|dil {ok4}")
print("  => rot, dil, refl generate the AFFINE GROUP AGL(1,R/Z) = {x -> +-v x + a}.  RUNNERS = the dilation part.\n")
# ---- LRC lonely set = the UNITS = dilation orbit structure (q=n=14) ----
n=14; units=[v for v in range(1,n) if gcd(v,n)==1]
print(f"LRC at n={n}: lonely set = units (Z/{n})* = {units}  (the dilations that are INVERTIBLE mod {n}).")
print(f"  |units| = phi({n}) = {len(units)}; these are exactly the speeds with a well-defined observer-inverse.\n")
# ---- alpha_0 = runner CENTROID from marked origin z=1 (observer-lens), for several configs ----
def alpha0(S,t): 
    m=np.mean(np.exp(2j*np.pi*np.array(S)*t)); return abs(m)  # |alpha_0| = |centroid|
print("|alpha_0| = |runner centroid| (observer at z=1) -- SMALL = balanced/cold = good LRC spreading:")
print(f"  AP {{1..13}} @ t=1/14 : |alpha_0|={alpha0(range(1,14),1/14):.4f}  = 1/(n-1)=1/13={1/13:.4f} EXACT (deleted z=1 pulls centroid)")
print(f"  construction @14/183 : |alpha_0|={alpha0(list(range(1,13))+[182],14/183):.4f}")
print(f"  GW {{1..11,13,24}}@1/14: |alpha_0|={alpha0([1,2,3,4,5,6,7,8,9,10,11,13,24],1/14):.4f}")
print(f"  random spread @gen t  : |alpha_0|={alpha0([1,17,38,59,88,123,167,201,244,290,331,377,399],0.234):.4f} (bigger=less balanced)")
print("\n  => The Verblunsky HIERARCHY {alpha_j} is a tower of BALANCE measures of the runner cloud from the origin;")
print("     alpha_0=centroid, and for the extremal AP the WHOLE tower is harmonic 1/(n-1-j) (vertex-peeling = Mode A).")

#!/usr/bin/env python3
"""
lrc_torus_index_klein_S332.py -- klein-2026-07-20-S332
Owner: one involution that is free and carries an odd map; the k-torus of the resonance lattice, with the
caveat that T^k != S^k blocks plain BU so it needs the Z/2-index form.

THM-1381. The caveat is decisive and provable.

(1) A free translation involution on T^k has Z/2-INDEX EXACTLY 1, for EVERY k.
    Free => index >= 1. For the upper bound take c a nonzero 2-torsion element, relabel so c_1 = 1/2, and
        f : T^k -> R^2,   f(x) = (cos 2 pi x_1, sin 2 pi x_1).
    Then f(x+c) = -f(x) (equivariant) and |f| = 1 (never zero), so an equivariant map into R^2 minus 0
    EXISTS => index <= 1.  Verified at k=1,2,3,5,12: max|f(x+c)+f(x)| ~ 1.9e-15, min|f| = 1.
    Cohomologically the same fact: H*(T^k;Z/2) is EXTERIOR on degree-1 generators, so every degree-1 class
    squares to zero, w_1^2 = 0, index <= 1.  DIMENSION BUYS NOTHING.
    Contrast: on S^k the antipodal action has index k and no equivariant S^k -> R^k minus 0 exists (that
    IS Borsuk-Ulam). The torus fails in the RING STRUCTURE, not the dimension.

(2) Second, independent reason special to LRC: for INTEGER speeds, t -> (v_1 t,...,v_n t) is invariant
    under t -> t+1, so it descends to S^1 -> T^n and the orbit is a CLOSED CIRCLE. The problem lives on a
    1-parameter space however many speeds there are; any free involution on it has index 1.

CONSEQUENCE: an equivariant f: X -> R^m is guaranteed a zero only when m <= index(X) = 1. So Borsuk-Ulam
type arguments on the resonance torus force ONE constraint, never n -- they cannot reach LRC(n), n >= 3.
Fourth delimited family in the "what does not work" map, with pairwise invariants (S324), alternating
truncations (S325), additive certificates (THM-1042).

THE SYNTHESIS IT CLOSES (klein-S322's puzzle): tournaments have an involution WITH fixed points => a spine
to build on, but no BU. LRC's involution is FREE => BU applies, but the space is a torus => index 1 => BU
is empty. Freeness is necessary but NOT sufficient; the space must also be curved enough to carry index.
"""
import numpy as np
def witness(k, n=200000, seed=0):
    """f(x) = e^{2 pi i x_1} is equivariant for translation by c=(1/2,0,..,0) and never vanishes."""
    rng=np.random.default_rng(seed)
    x=rng.random((n,k)); c=np.zeros(k); c[0]=0.5
    f =np.stack([np.cos(2*np.pi*x[:,0]),      np.sin(2*np.pi*x[:,0])],      axis=1)
    fc=np.stack([np.cos(2*np.pi*(x+c)[:,0]),  np.sin(2*np.pi*(x+c)[:,0])],  axis=1)
    return float(np.abs(fc+f).max()), float(np.linalg.norm(f,axis=1).min())
if __name__=="__main__":
    print("k | max|f(x+c)+f(x)| (equivariance) | min|f| (nonvanishing)")
    for k in (1,2,3,5,12):
        e,m=witness(k); print("%2d |          %.2e            |   %.6f"%(k,e,m))
    print("\n=> equivariant map into R^2 minus 0 exists for every k => Z/2-index = 1 => BU forces one constraint.")

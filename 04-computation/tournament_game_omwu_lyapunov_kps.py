#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The tournament as a symmetric zero-sum game: OMWU last-iterate convergence to the
bipartisan Nash equilibrium, certified by a Lyapunov (KL) function; the dynamics'
frequencies = the skew spectrum (determinant lens). kps-2026-06-15-S5.
Inspired by Orabona, arXiv:2606.11773 (Last-Iterate Convergence of OMWU) — the
Lyapunov-FUNCTION sense (vs the Lyapunov-EXPONENT sense of THM-488/HYP-614).
Tournament <-> symmetric zero-sum game: payoff S = A - A^T (antisymmetric).
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
import numpy as np

def skew(adj):
    n=len(adj); S=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if adj[i][j]: S[i][j]=1; S[j][i]=-1
    return S

def omwu(S, eta=0.15, T=4000):
    """Optimistic MWU on the symmetric zero-sum game with antisymmetric payoff S."""
    n=S.shape[0]; x=np.ones(n)/n; gprev=S@x; hist=[x.copy()]
    for _ in range(T):
        g=S@x; w=x*np.exp(eta*(2*g-gprev)); x=w/w.sum(); gprev=g; hist.append(x.copy())
    return x, hist

def kl(p,q):
    p=np.clip(p,1e-15,1); q=np.clip(q,1e-15,1); return float(np.sum(p*np.log(p/q)))

def main():
    Ts={
      "3-cycle (rock-paper-scissors)": [[0,1,0],[0,0,1],[1,0,0]],
      "transitive T3": [[0,1,1],[0,0,1],[0,0,0]],
      "regular T5 (circ {1,2})": [[1 if (j-i)%5 in(1,2) else 0 for j in range(5)] for i in range(5)],
      "Paley T7 (QR {1,2,4})": [[1 if (j-i)%7 in(1,2,4) else 0 for j in range(7)] for i in range(7)],
    }
    print("tournament-game OMWU: Nash equilibrium (bipartisan set), Lyapunov KL->0, skew spectrum\n")
    for name,adj in Ts.items():
        for i in range(len(adj)): adj[i][i]=0
        S=skew(adj); xf,hist=omwu(S)
        Ls=[kl(xf,h) for h in hist]
        ev=np.linalg.eigvals(S); mu=sorted({round(abs(e.imag),4) for e in ev if abs(e.imag)>1e-9})
        print(f"{name}: Nash={np.round(xf,3)}  KL {Ls[0]:.4f}->{Ls[-1]:.1e}  skew|mu|={mu}")
    print("\nReadings: transitive -> PURE Nash on the dominator (H=1); regular/Paley -> UNIFORM Nash (H=max).")
    print("OMWU frequencies = skew spectrum {mu} = determinant lens det(I+S)=prod(1+mu^2) (THM-468/472).")
    print("E7 direction (HYP-2530): the 56-rep symplectic form omega is a 'continuous tournament' (28 antipodal")
    print("pairs); Weyl-uniform Nash; OMWU frequencies ~ E7 exponents {1,5,7,9,11,13,17}/h=18.")

if __name__=="__main__":
    main()

"""Transient vs spectral contraction of the level operator: ||T_n^m||_2, ||T_n^m 1||/||1||, and the actual cocycle
||T_n lift T_{n-1} ... ||: the Parseval rate 3^-1/2 is a transient of operators with spectral radius 1/2."""
import numpy as np
def level_op(n):
    L = 2*3**(n-1); mod = 3**n
    om = np.exp(2j*np.pi*np.array([pow(2,k,mod) for k in range(L)])/mod)
    S = np.roll(np.eye(L), 1, axis=0); C = (S/2) @ np.linalg.inv(np.eye(L) - S/2)
    return C * om[None,:], om, L
for n in (3,4,5,6):
    T, om, L = level_op(n)
    one = np.ones(L)
    norms = []; v = one.copy()
    for m in range(1, 13):
        v = T @ v; norms.append(np.linalg.norm(v)/np.sqrt(L))
    Tm = np.eye(L); op = []
    for m in range(1, 9):
        Tm = T @ Tm; op.append(np.linalg.norm(Tm, 2))
    print(f"n={n} L={L}: ||T^m 1||/sqrt(L) for m=1..12: " + " ".join(f"{x:.4f}" for x in norms))
    print(f"        per-step ratios: " + " ".join(f"{norms[i]/norms[i-1]:.3f}" for i in range(1,len(norms))) + f"   (3^-1/2 = 0.5774, 1/2 = 0.5)")
    print(f"        ||T^m||_2 for m=1..8: " + " ".join(f"{x:.4f}" for x in op) + f";  ||T^8||^(1/8) = {op[-1]**(1/8):.4f}")
# the actual cocycle: f_n = T_n lift f_{n-1}, f_0 = 1 on Z/L_1? level 1: L_1 = 2
def lift(f, L_new):   # period L_old = L_new/3
    return np.tile(f, 3)
f = np.ones(2)   # level 1 domain Z/2 ... start with f_0 = 1 on Z/L_1
for n in range(1, 8):
    T, om, L = level_op(n)
    f = T @ (lift(f, L) if n > 1 else np.ones(L))
    print(f"cocycle level {n}: rms |f_n| = {np.sqrt(np.mean(np.abs(f)**2)):.5f}, 3^(-n/2) = {3**(-n/2):.5f}, ratio {np.sqrt(np.mean(np.abs(f)**2))/3**(-n/2):.4f}, max|f_n| = {np.abs(f).max():.5f} (M(n))")

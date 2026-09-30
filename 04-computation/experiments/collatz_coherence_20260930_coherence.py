"""(a) exact energy transfer on the unit cycle: ||T_n g||^2 = sum_xi |Ghat(xi)|^2 E_g(xi), E_g the energy spectrum of omega_n g
    (mean 1/3); the cocycle's per-level covariance with the symbol; (b) the half-turn: f_n(k+L/2) = conj f_n(k) and
    R T_n = conj(T_n) R; (c) the level-2 hexagon operator; (d) every periodic parity word is a cycle of every parity graph G_k;
    (e) motif census of the parity graphs G_sigma (3-node and 5-node induced subdigraphs) for sigma = +, -, random."""
import numpy as np, itertools, math, random
from collections import Counter
# ---------- (a),(b),(c): the cocycle on the full unit cycle up to n = 12
def level(n):
    L = 2*3**(n-1); mod = 3**n
    om = np.exp(2j*np.pi*np.array([pow(2,k,mod) for k in range(L)])/mod)
    return L, om
def apply_T(om, g, A=60):
    L = len(g); out = np.zeros(L, dtype=complex); P = om*g
    for a in range(1, A+1): out += 2.0**(-a)*np.roll(P, a)
    return out
f = np.ones(2); print("(a) energy transfer: level n, ||f_n||^2/||f_{n-1}||^2 (lifted), covariance part = ratio - 1/3, low-band (|xi|<L/20) share of the twisted energy spectrum vs its symbol-weighted share")
prev_norm2 = None
for n in range(1, 13):
    L, om = level(n)
    g = np.tile(f, 3) if n > 1 else np.ones(L)
    tw = om*g; Fh = np.fft.fft(tw); E = np.abs(Fh)**2/L      # sum E = ||tw||^2 = ||g||^2
    xi = np.arange(L); Gh = 1/(2*np.exp(2j*np.pi*xi/L) - 1)          # symbol on Z/L (untruncated kernel)
    f = np.fft.ifft(np.fft.fft(g*om)*np.conj(np.fft.fft(np.array([0]+[2.0**-a for a in range(1,L)])))) if False else apply_T(om, g)
    ratio = np.sum(np.abs(f)**2)/np.sum(np.abs(g)**2)
    pred = np.sum(np.abs(Gh)**2*E)/np.sum(E)
    low = np.abs(((xi+L//2) % L) - L//2) < L/20
    print(f"  n={n:2d} L={L:6d}: ratio {ratio:.6f}  symbol-weighted spectral mean {pred:.6f}  (1/3 = 0.333333; covariance {ratio-1/3:+.5f});  low-band mass {E[low].sum()/E.sum():.4f} (uniform {low.mean():.4f}); half-turn check max|f(k+L/2)-conj f(k)| = {np.abs(np.roll(f, -L//2)-np.conj(f)).max():.1e}; max|f_n| = {np.abs(f).max():.5f}")
# (b) intertwining: R T_n = conj(T_n) R on random vectors
L, om = level(5); v = np.random.randn(L)+1j*np.random.randn(L)
lhs = np.roll(apply_T(om, v), -L//2); rhs = apply_T(np.conj(om), np.roll(v, -L//2))
print("(b) intertwining R T = Tbar R at n=5: max diff", np.abs(lhs-rhs).max(), "; omega(k+L/2) = conj omega(k):", np.abs(np.roll(om,-L//2)-np.conj(om)).max())
# (c) level-2 hexagon operator
L, om = level(2); S = np.roll(np.eye(L), 1, axis=0); C = (S/2)@np.linalg.inv(np.eye(L)-S/2); T2 = C*om[None,:]
ev = np.linalg.eigvals(T2); print("(c) level-2 hexagon: units mod 9 in the order 2^k:", [pow(2,k,9) for k in range(6)], "; eigenvalue moduli", np.round(np.sort(np.abs(ev)),6), "; 63 lambda^6 + lambda^3 - 1 at eigs:", np.abs(63*ev**6+ev**3-1).max())
# (d) periodic words are cycles of every parity graph G_k (sigma=+): x_w = c_w/(2^p - 3^a) is a 2-adic integer; its residues mod 2^k close a walk
def T(x, s): return x//2 if x%2==0 else (3*x+s)//2
words = {"0":[0],"10 (1,2)":[1,0],"1 (-1)":[1],"110 (-5,-7,-10)":[1,1,0],"1111 0 111 000 (-17)":[1,1,1,1,0,1,1,1,0,0,0]}
def periodic_point(w):   # x with T-parity word w periodic (sigma=+), as an exact rational
    from fractions import Fraction
    p=len(w); a=sum(w); c=0
    for j,b in enumerate(w):
        if b: c = 3*c + 2**j
    return Fraction(c, 2**p - 3**a)
print("(d) periodic words as cycles of G_k: for each word the rational periodic point and, for k=1..8, whether the residues mod 2^k of its orbit form a closed walk of G_k (node s -> lifts of T(s) mod 2^(k-1)):")
for name,w in words.items():
    x = periodic_point(w); p=len(w)
    ok=[]
    for k in range(1,9):
        m = 2**k; num, den = x.numerator, x.denominator; inv = pow(den, -1, m); r0 = (num*inv) % m
        # walk: r -> T-step according to the word letter; check parity matches w at each step and closes
        r=r0; good=True
        for j in range(p):
            if r%2 != w[j]: good=False; break
            r = (r//2) % m if w[j]==0 else ((3*r+1)//2) % m
        good = good and (r == r0)
        ok.append(good)
    print(f"   {name:24s} x_w = {x}: closed walks at k=1..8: {ok}")
# (e) motif census of G_sigma at level k
def parity_graph(k, sigma):
    N = 2**k; adj = {s:set() for s in range(N)}
    for s in range(N):
        t = T(s, sigma(s)) % (2**(k-1)) if s%2 else (s//2) % (2**(k-1))
        adj[s] = {t, t + 2**(k-1)}
    return adj
def canon(sub, adj):
    """canonical form of the induced subdigraph on the vertex tuple sub: sorted adjacency matrices over all permutations (small)"""
    n=len(sub); best=None
    for perm in itertools.permutations(range(n)):
        m=tuple(tuple(1 if sub[perm[j]] in adj[sub[perm[i]]] else 0 for j in range(n)) for i in range(n))
        if best is None or m<best: best=m
    return best
def census(adj, size):
    nodes=sorted(adj); c=Counter()
    for sub in itertools.combinations(nodes, size): c[canon(sub, adj)] += 1
    return c
random.seed(1)
for k in (4,5):
    N=2**k
    sig_plus = lambda s: 1; sig_minus = lambda s: -1
    rs = {s: random.choice([1,-1]) for s in range(N)}; sig_rand = lambda s: rs[s]
    cp, cm, cr = parity_graph(k, sig_plus), parity_graph(k, sig_minus), parity_graph(k, sig_rand)
    for size in (3,5):
        if size==5 and k==5: 
            Cp, Cm, Cr = census(cp,5), census(cm,5), census(cr,5)
        else:
            Cp, Cm, Cr = census(cp,size), census(cm,size), census(cr,size)
        same = (Cp == Cm)
        diff = sum(abs(Cp[t]-Cr[t]) for t in set(Cp)|set(Cr))
        print(f"(e) level k={k}, {size}-node induced motifs: types(+) = {len(Cp)}, types(-) = {len(Cm)}, identical(+,-) = {same}; random sigma: types {len(Cr)}, L1 distance to + = {diff} of {sum(Cp.values())} subsets")

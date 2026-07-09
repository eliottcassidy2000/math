"""
mac-mini-2026-07-08-S59 -- VERIFY the exact closed form for the Fourier coefficients of the
uncovered-measure function 𝒲 on T^{k-1} (the shared 𝒲̂-decay constant for THM-664 & opus-S157).

𝒲(phi_1,..,phi_{k-1}) = uncovered measure of the circle by k arcs of length l=1/7 at
{0, phi_1,..,phi_{k-1}} (phase 0 pinned).  CLAIM (derived):
   𝒲̂(n) = (6/7)^z * [prod_{i: n_i != 0} b(n_i)] * Q(N),
where z=#zeros among n_1..n_{k-1}, N = sum_i n_i, and
   b(m) = (1 - e(m*l))/(2*pi*i*m)        (m != 0),   |b(m)| = |sin(pi m l)|/(pi|m|),
   Q(0) = 1-l = 6/7,   Q(N) = (e(-N*l)-1)/(2*pi*i*N)  (N != 0),   |Q(N)|=|sin(pi N l)|/(pi|N|).
KEY: b(m)=0 whenever 7|m (sin(pi m /7)=0); Q(N)=0 whenever 7|N, N!=0.  => 7-commensurate support.

Verify by direct numerical Fourier integral over T^{k-1} for k=3 (2D) and k=4 (3D).
"""
import numpy as np
import cmath
l = 1.0/7.0

def uncovered(phases):
    """exact uncovered measure of circle by arcs [p, p+l) at 'phases' (list, p in [0,1))."""
    ps = sorted(p % 1.0 for p in phases)
    # coverage = union of [p,p+l); uncovered = 1 - covered. Compute via sweep of arc union length.
    intervals = []
    for p in ps:
        a, b = p, p+l
        if b <= 1: intervals.append((a,b))
        else: intervals.append((a,1.0)); intervals.append((0.0,b-1.0))
    intervals.sort()
    covered = 0.0; cur_a, cur_b = None, None
    for a,b in intervals:
        if cur_a is None: cur_a,cur_b = a,b
        elif a <= cur_b: cur_b = max(cur_b,b)
        else: covered += cur_b-cur_a; cur_a,cur_b = a,b
    if cur_a is not None: covered += cur_b-cur_a
    return 1.0 - covered

def b(m):
    if m == 0: return None
    return (1 - cmath.exp(2j*cmath.pi*m*l))/(2j*cmath.pi*m)
def Q(N):
    if N == 0: return 1-l
    return (cmath.exp(-2j*cmath.pi*N*l)-1)/(2j*cmath.pi*N)
def What_formula(n):
    z = sum(1 for x in n if x == 0)
    N = sum(n)
    prod = 1.0+0j
    for x in n:
        if x != 0: prod *= b(x)
    return (1-l)**z * prod * Q(N)

def What_direct(n, k, G):
    """direct 2D/3D Fourier integral of 𝒲 over T^{k-1}, grid G per axis."""
    d = k-1
    axes = [np.arange(G)/G for _ in range(d)]
    total = 0j
    if d == 2:
        for i1 in range(G):
            for i2 in range(G):
                ph = [0.0, axes[0][i1], axes[1][i2]]
                w = uncovered(ph)
                total += w*cmath.exp(-2j*cmath.pi*(n[0]*axes[0][i1]+n[1]*axes[1][i2]))
        return total/G**2
    else:  # d==3
        for i1 in range(G):
            for i2 in range(G):
                for i3 in range(G):
                    ph = [0.0, axes[0][i1], axes[1][i2], axes[2][i3]]
                    w = uncovered(ph)
                    total += w*cmath.exp(-2j*cmath.pi*(n[0]*axes[0][i1]+n[1]*axes[1][i2]+n[2]*axes[2][i3]))
        return total/G**3

print("VERIFY closed form 𝒲̂(n) = (6/7)^z prod b(n_i) Q(N)\n")
print("k=3 (T^2), grid 210:")
for n in [(0,0),(1,0),(0,1),(1,1),(1,-1),(2,1),(3,0),(7,0),(1,6),(2,-2),(1,2)]:
    f = What_formula(n); dct = What_direct(n,3,210)
    ok = abs(f-dct) < 2e-3
    print(f"  n={str(n):8s}: formula={f.real:+.5f}{f.imag:+.5f}i  direct={dct.real:+.5f}{dct.imag:+.5f}i  "
          f"|diff|={abs(f-dct):.2e} {'OK' if ok else 'MISMATCH'}")
print("\nk=4 (T^3), grid 84:")
for n in [(0,0,0),(1,0,0),(1,1,0),(1,1,-2),(1,-1,0),(1,1,1),(2,0,0),(7,0,0),(1,2,-3)]:
    f = What_formula(n); dct = What_direct(n,4,84)
    ok = abs(f-dct) < 4e-3
    print(f"  n={str(n):10s}: formula={f.real:+.5f}{f.imag:+.5f}i  direct={dct.real:+.5f}{dct.imag:+.5f}i  "
          f"|diff|={abs(f-dct):.2e} {'OK' if ok else 'MISMATCH'}")
print(f"\nSANITY: 𝒲̂(0)=(6/7)^k should equal (6/7)^3={float((6/7)**3):.5f}, (6/7)^4={float((6/7)**4):.5f}")
print(f"  formula 𝒲̂(0,0)={What_formula((0,0)).real:.5f}  𝒲̂(0,0,0)={What_formula((0,0,0)).real:.5f}")
print("  (these are the k=3,4 iid means; note (6/7)^z*Q(0)=(6/7)^{k-1}*(6/7)=(6/7)^k)")

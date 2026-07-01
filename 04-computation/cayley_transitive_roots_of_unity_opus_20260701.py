"""The transitive tournament's Cayley spectrum = ROOTS OF UNITY (the LRC runner cloud). Pin the odd/even pattern.
Extend: span/tile-flip spectral effect; the rotational (circulant) tournament -> exact roots too."""
import numpy as np
def cayley_angles(A):
    S=(A-A.T).astype(float); n=A.shape[0]
    U=np.linalg.solve(np.eye(n)+S, np.eye(n)-S)
    return np.sort(np.round(np.angle(np.linalg.eigvals(U)),6))
def transitive(n):
    A=np.zeros((n,n),int)
    for i in range(n):
        for j in range(i+1,n): A[i,j]=1
    return A
print("TRANSITIVE tournament Cayley spectrum vs roots of unity:")
for n in range(3,10):
    ang=cayley_angles(transitive(n))/np.pi   # in units of pi
    nth = sorted(np.round((2*np.arange(n)/n +1)%2 -1,6))     # n-th roots angles/pi: 2k/n mod 2 -> [-1,1)
    nth_roots = sorted(np.round(((2*np.arange(n)/n +1)%2)-1,6))
    # n-th roots of unity angles/pi = 2k/n for k, reduced to (-1,1]
    rootsN = sorted(np.round((((2*np.arange(n))/n +1)%2)-1,6))
    roots2N_odd = sorted(np.round((((2*np.arange(n)+1)/n +1)%2)-1,6))  # primitive-ish 2n-th: (2k+1)/n
    match_n = np.allclose(sorted(np.round(ang,6)), rootsN, atol=1e-5)
    match_2n = np.allclose(sorted(np.round(ang,6)), roots2N_odd, atol=1e-5)
    tag = "= n-th roots of unity" if match_n else ("= {(2k+1)/n}pi (prim 2n-th roots)" if match_2n else "???")
    print(f"  n={n} ({'odd ' if n%2 else 'even'}): angles/pi={list(np.round(ang,3))}   {tag}")
print("\n=> ODD n: transitive Cayley spectrum = the n-th roots of unity {e^(2pi i k/n)} (incl 1 at angle 0).")
print("   EVEN n: = {e^(i pi (2k+1)/n)} = the primitive 2n-th roots (no fixed point 1). The parity split again.")
# rotational/circulant tournament also -> exact roots
def rotational(n):  # i beats i+1,..,i+(n-1)/2 (odd n)
    A=np.zeros((n,n),int)
    for i in range(n):
        for d in range(1,(n-1)//2+1): A[i,(i+d)%n]=1
    return A
for n in [5,7]:
    print(f"  rotational C_{n}: Cayley angles/pi = {list(np.round(cayley_angles(rotational(n))/np.pi,3))}  (circulant => exact roots, the runner cloud)")

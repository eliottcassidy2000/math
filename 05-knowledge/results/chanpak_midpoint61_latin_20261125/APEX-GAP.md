# How far apart must the endpoints of a centered apex flip be?

**PROVED.** Let \(T=2^{p-1}\), \(p\ge3\), and \(E_p=2T^2-T\). In a tournament with an apex and a transitive spine, reverse one apex contact. If its two Hamilton path counts have arithmetic midpoint exactly \(E_p\), let \(D\) be their absolute half-difference. The two flank addresses have separation \(d\) satisfying
\[
                    (2D+1)^2\ \ge\ d^2\ \ge\ 4T+2 .
\]
The half-gap \(D\) is odd. **FINITE-EXACT refinement at exponent 61.** Every such centered apex flip at any order must have
\[
                     D\ge 34{,}376{,}511{,}491 .
\]
An edge crossing exactly at \(E_{61}\) with prescribed endpoints \(E_{61}\pm1\) would have \(D=1\). This gives a quantitative counterpart to the [one-apex endpoint obstruction](SPINE-ARC.md). Existence of an apex flip at the much larger lower bound is open.

## A separation bound for every exponent

Let \(u,w\) be the left address and the reversed complemented right address. Put \(S=u+w,\ d=|w-u|\). The [one-contact formula](../chanpak_apex_flip_local_global_20261122/README.md) gives
\[
 M=2uw+\lceil u/2\rceil+\lceil w/2\rceil+f(u)+f(w)+1,
 \qquad D=|\lceil u/2\rceil-\lceil w/2\rceil|.
\]
Here \(f(0)=0,\ f(2z+1)=f(z),\ f(2z)=f(z)+z+[z\text{ odd}]\).
Both digital reads are even and \(0\le f(z)\le z\). A midpoint at the even integer \(E_p\) between two odd endpoints requires odd \(D\), so \(d>0\). Rounding gives \(d\le2D+1\).

First \(S\ge2T\). For \(S\le2T-3\),
\[
 M\le S^2/2+(S+2)/2+S+1 < 2T^2-T.
\]
For \(S=2T-2\) the positive separation \(d\) is even, hence at least 2. Use \(2uw=(S^2-d^2)/2\) in the same bound to get \(M\le E_p-1\). For \(S=2T-1\), the ceiling sum is \(T\), which is even, and thus \(M\) is odd. Finally, for every remaining sum,
\[
                   d^2\ge S^2+S+2-2E_p\ge4T+2.             \tag{1}
\]
This follows by dropping the two nonnegative digital reads in the midpoint identity.

## An inverse for a fixed address sum

Fix \(S\), let \(L\) be its bit length and \(P=2^L\). Write \(\phi(z)=\sum_i z_i2^{-i-1}\) for the binary radical inverse of the bits of \(z\), and let \(t_i=u_i+w_i\). If \(s_i\) are the bits of \(S\), the ordinary addition carries satisfy
\[
 t_i=s_i+2c_{i+1}-c_i\in\{0,1,2\},\qquad c_0=c_L=0.
\]
Thus each carry sequence fixes \(t_0\) and \(R=\phi(u)+\phi(w)=\sum_i t_i2^{-i-1}\), regardless of which flank owns a particular one. Take \(u\le w\) by symmetry. The reversed-product error gives \(|f(z)-z(1-\phi(z))|\le L+1\) for \(0\le z<P\). Consequently
\[
 \left|d^2-J\right|\le d+4(L+1),\qquad
 J=S^2+3S+t_0+2-2E_p-SR .                                  \tag{2}
\]
Indeed the ceiling sum is \((S+t_0)/2\). Twice the sum of the leading digital terms is \(2S-SR+d(\phi(u)-\phi(w))\). Each \(\phi\) lies in \([0,1)\). [Reversed-product estimate](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md)

Whenever \(d>4(L+1)+1\), the error in (2) is smaller than \(2d-1\). Hence \(J\) is strictly between \((d-1)^2\) and \((d+1)^2\). For each carry sequence it suffices to try \(\lfloor\sqrt{\lfloor J\rfloor}\rfloor\) and its successor. The parity of the midpoint also forces \(S+t_0\equiv2\pmod4\). Checking the original digital recurrence on the deduplicated candidates gives an exact answer on the entire sum slice. The exponentially many assignments of individual one bits are unnecessary.

At \(T=2^{60}\), (1) gives the required large \(d\) on every tested slice. For the 1,024 values \(2T\le S\le2T+1023\), the inverse leaves **695,466** carry sequences after the parity test and **223,038** candidate address pairs across the slices. None has midpoint \(E_{61}\). Thus \(S\ge2T+1024\) and (1) gives \(d\ge68{,}753{,}022{,}980\). Since \(2D+1\ge d\) and \(D\) is odd, the claimed bound follows.

The following self-contained Python 3 block reproduces the counts. A positive control also starts from every address pair with sum in \([90,199]\), separation at least \(4L+8\), and even midpoint. It recovers the supplied separation on all 2,912 such inputs.

~~~python
from math import isqrt

def f(z):
    out = 0
    while z:
        z, rem = divmod(z, 2)
        if rem == 0:
            out += z + (z & 1)
    return out

def candidates(S, E):
    width = S.bit_length()
    P = 1 << width
    states = [(0, 0, 0)]  # carry into current bit, reversed digit sums, low digit
    for i in range(width):
        bit = (S >> i) & 1
        states = [(nxt, 2*rev + t, t if i == 0 else low)
                  for carry, rev, low in states for nxt in (0, 1)
                  if 0 <= (t := bit + 2*nxt - carry) <= 2]
    candidates = set()
    phases = 0
    for carry, rev, low in states:
        if carry or (S + low) % 4 != 2:
            continue
        phases += 1
        J = S*S + 3*S + low + 2 - 2*E - (S*rev + P-1)//P
        if J < 0:
            continue
        k = isqrt(J)
        for d in (k, k+1):
            if d <= S and (S-d) % 2 == 0:
                candidates.add(d)
    return phases, candidates

def scan(last):
    T = 1 << 60
    E = 2*T*T-T
    phases = tests = 0
    for s in range(last+1):
        S = 2*T+s
        count, ds = candidates(S,E)
        phases += count
        tests += len(ds)
        for d in ds:
            u,w=(S-d)//2,(S+d)//2
            mid = 2*u*w + (u+1)//2 + (w+1)//2 + f(u)+f(w)+1
            if mid == E:
                print("hit",s,u,w)
    print("excess through",last,"phases",phases,"exact reads",tests)

planted = 0
for S in range(90, 200):
    for d in range(4*S.bit_length() + 8 + S%2, S+1, 2):
        u, w = (S-d)//2, (S+d)//2
        E = 2*u*w + (u+1)//2 + (w+1)//2 + f(u) + f(w) + 1
        if E % 2 == 0:
            assert d in candidates(S, E)[1]
            planted += 1
assert planted == 2912
scan(1023)
~~~

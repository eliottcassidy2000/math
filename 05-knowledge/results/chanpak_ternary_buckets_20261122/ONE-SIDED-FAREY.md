# A target inverse for the marked Farey receiver

Let \(E=2^{60}(2^{61}-1)\). Form the tagged ordinal sum \(A\to D\to V_n\) with its single shared apex. Here \(A\) is any tagged tournament on one through five vertices, the ordered arc flip is \(1,4\) in the six-vertex block \(D=(16,5,9,19,38,15)\) (zero-based outneighbor masks), and \(V_n\) consists of a transitive spine of length \(n\) with two additional vertices obeying \(\mathrm{spine}\to a\to b\to\mathrm{spine}\). All three blocks permit arbitrary vertex tags. Outside this core append the binary word flanks of the star chart. A **permanent unit flip** adds or subtracts exactly two with the same sign at every pair of outside addresses. The order and tags in the two hosts agree, including both flanks. [Marked chart and join law](../chanpak_decimal_kernel_20261122/E15-FAREY.md)

**PROVED.** For any fixed tagged left prefix and flipped prefix that agree in their first four coordinates, the desired signed difference determines at most four receiver taggings at each \(n\ge 2\), unless every obtainable difference is divisible by \(2^n+1\). These four possibilities are indexed by the two cap tags. The rule finds the receiver before inspecting either outside address.

Use the seven readings \((h,q,c,t,e,r,k)\) of the cited join law. A reversal confined to \(D\) preserves the first four readings of the whole sum exactly when it preserves those of \(D\): the products detecting \(h,q\) have positive factors, and with those fixed the joins scale \(\Delta c,\Delta t\) by positive factors. Permanent unit flips must preserve these coefficients, by taking mixed differences in two addresses for \(q\), and by evaluating the faces at even and odd addresses to separate \(h,c,t,e,r\). For example on either face a putative constant difference requires \(2\Delta s+2\Delta h=0\) at address 2 and \(6\Delta s+4\Delta h=0\) at address 6, forcing both differences to vanish. Then address 1 isolates its parity coefficient. The mixed difference at \((1,1)\) isolates \(q\).

Write \(v,w\) for the original and flipped profiles of \(A\to D\). Joining a fixed \(V_n\) on the right automatically resets the difference in \(e\), while the identical left block forces the difference in \(r\) to be zero. Direct substitution gives
\[
       k_{\rm flipped}-k_{\rm original}
        = G(w_k-v_k)+(w_e-v_e)r_V,\qquad G=2^n+1.
\]
A zero slope permits no difference \(\pm2\) for \(n\ge2\). At nonzero slope, solve this equation for the specific sign. Each untagged cap contributes one unit to \(r_V\). The start counts along the transitive spine, from first to last vertex, are \(2^{n-1},2^{n-2},\dots,1\). Indeed a Hamiltonian path starting at the \(i\)-th spine vertex must use \(a,b\) consecutively. Choose which vertices later than \(i\) lie before \(a\); all other spine vertices appear after \(b\) in increasing order. A path starting at either cap has one possible completion. Accordingly the \(n\)-bit binary expansion of \(r_V-\#\{\text{untagged caps}\}\), when in range, specifies exactly the untagged spine vertices. This proves the inverse.

**FINITE-EXACT.** For every \(n\ge24\) this family contains no permanent unit flip with Hamiltonian counts \(E-1,E+1\) having at least one empty external flank. For \(n\ge52\), two unrestricted external flanks also fail. Smaller spine lengths and the remaining two-flank charts are outside this exclusion.

For \(24\le n\le120\) the replay reads every isomorphism class through order five and all vertex subsets. This gives 344 distinct marked left profiles. Exactly two of the 64 middle tag masks preserve the first four readings and change the profile. Joining yields 688 distinct ordered prefix-profile pairs. Of the 68,436 resulting prefix/receiver choices with the required gap, 36,720 have lower constant at most \(E-1\). The cap inverse produces 3,294 distinct triples (length, start reading, cap mask). The one-sided search returned zero hits after 4,280,161 Python evaluations and 289,715,580 candidates from the accelerated reversed-suffix search. Lower chart constants exceeding \(E-1\) exclude both flanks because all remaining contributions are nonnegative. For \(n\ge121\) even the induced receiver alone has \(2^n+1>E-1\) Hamiltonian paths. Inserting a vertex injects paths by choosing a valid insertion position and recovering the original by deletion.

For positive addresses the lower chart is at least \(k+q+c+t\); the final passing length is 69. For the outstanding charts enumerate the smaller of the two positive addresses through \(\lfloor\sqrt{(E-1-k)/q}\rfloor\), fixing its orientation, and use the interval and residue of the other address at each parity. This gives 10,604 complete positive-pair searches for \(n\ge52\), with zero hits. [Two-word integer inverse](../chanpak_decimal_kernel_20261122/fermat_receiver.py)

On a single face solve \(s x+h f(x)=M\) where \(M=E-1-k-d[x\text{ odd}]\), \((s,d)=(t,r)\) or \((c,e)\). Split \(x=BX+i\) with \(B=2^b\), \(b\ge1\) and \(0\le i<B\) of the chosen parity. The binary recurrence gives
\[
 f(BX+i)=f(i)+(B-1-\operatorname{rev}_b(i))X+f(X)+[i<B/2][X\text{ odd}].
\]
The last two summands together are between zero and \(X\), since \(f(X)\le X\) and \(f(2u+1)=f(u)\). With \(R=M-si-hf(i)\) and \(S=Bs+h(B-1-\operatorname{rev}_b(i))\), the remaining possibilities satisfy
\[
     \left\lceil\frac{R}{S+h}\right\rceil\le X\le
     \left\lfloor\frac{R}{S}\right\rfloor,
     \qquad BsX\equiv R\pmod h.
\]
Testing the exact recurrence on this entire progression exhausts words of every length. The split width changes runtime only. The replay checks the cap transfer at small order against literal tournament path counts and checks the complete inverse on planted and adjacent small readings. [Earlier split inverse](../chanpak_decimal_kernel_20261122/check_farey_axis.py) · [Expanded receiver census](farey_face_tags.py)

For the accelerated replay compile the repository’s [C++ suffix solver](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md) and append `--axis /path/to/farey_axis`. The self-contained, slower Python invocation from the repository root is `python3 -B 05-knowledge/results/chanpak_ternary_buckets_20261122/farey_face_tags.py`.

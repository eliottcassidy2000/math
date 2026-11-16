# Complete the family, keep the dyadic gate: QRH, mod 8 and three moving strands

There are two distinct meanings of “7/8 to 8/8.” In the [September 30 quasi-Riemann paper](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf), 7/8 is a boundary in a zero-free half-plane, not a count of residue classes: replacing it by 8/8 would shrink rather than strengthen the region. A stronger zero-free bound would move the number down, toward 1/2. Below, “8/8” means coverage of the eight classes of a normalized mod-8 readout (e.g. (H−1)/2 taken modulo 8, since H itself is odd).

The [October 5 alternate quasi-Riemann manuscript](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex) passes from twisted mean squares over the Eisenstein field to a completed cubic-theta family and a smaller dual scale. Our analog is strictly a **design transfer**: completing the address box gives an exact Fourier operator norm, and a dyadic test window handles the degenerate conductor levels. Neither a zero-free region nor the theta transformation is assumed here. This extends the complete **odd prime** box theorem in the earlier [repository QRH transfer](../chanpak_qrh_residual_20261115/README.md) to dyadic boxes.

## 1. Exact operator and smoothing at every dyadic level

Let q=2^k, and use a whole residue system for each address. Write
```text
R(x,y)=cxy+F(x)+G(y)pmod q
```
with arbitrary unary functions F,G (including checksum and parity terms from the physical three-center family; its c is 432). Define the Fourier collision sum
```text
S(a)=Σ_{x,ymod q} e_q(aR(x,y)).
```
For g=gcd(ac,q), the bilinear Fourier matrix K_a(x,y)=e_q(acxy) has K_aK_a^*(x,x')=q when x≡x' modulo q/g, and 0 otherwise. The blocks are q times all-one g×g matrices. Thus its exact operator norm is √(qg), and Cauchy on the two unary phase vectors gives
```text
|S(a)| ≤ q^(3/2) √gcd(ac,q).
```
This includes highly imprimitive even characters; treating them as primitive would be wrong at 2.

For h=2^r≤q choose the nonnegative triangular mask W_h(t)=#{i,j in [0,h): i−j≡t mod q}. Its Fourier coefficients are nonnegative squares. If U_j is their sum over **nonzero** frequencies divisible by 2^j, then
```text
U_j= (q h/2^j)-h^2    (0≤j≤k-r),qquad U_j=0    (j≥k-r).
```
Proof: compact the frequency subgroup to the transform on q/2^j points. If h is at most that size, orthogonality gives the total qh/2^j including the zero coefficient h². If h is larger, it is an integer multiple of this size and all nonzero subgroup frequencies vanish. This is the dyadic choice's exact cancellation of conductor levels, not an estimate.

Writing v=v₂(c) for nonzero c, expansion of the weight 2^(v₂(a)/2) by difference of consecutive levels yields
```text
Σ_{a≠0} √gcd(ac,q) · Ŵ_h(a)
 ≤ 2^(v/2)[U_0 + Σ_{j≥1}(2^(j/2)−2^((j−1)/2)) U_j]
 ≤ (1+1/√2) 2^(v/2) qh.
```
The last geometric sum uses U_j≤qh/2^j. Fourier inversion therefore bounds the wrapped weighted number of representations of a target r below by
```text
Σ_{x,y} W_h(R(x,y)−r) ≥ qh² − (1+1/√2)2^(v/2) q^(3/2) h.
```
In particular, **every** r has an attained readout within cyclic distance h−1 whenever h>(1+1/√2)2^{v/2}√q. For c=432 the sufficient constant is 4+2√2. Always using h=q trivially covers the cyclic group regardless of this sufficient condition. For q=8 and c=432 there is no bilinear mixing (432≡0), so this argument does not give exact mod-8 steering of the physical masks; actual digit choices, as in the earlier mod-8 atlas, remain necessary.

**Scale boundary.** A full address box may be taken over *actual* nonnegative representatives, and then its checksum functions really are unary on that box. The resulting radius is a **residue** radius. Completing a short physical box to q representatives without transporting the returned pair back into its original absolute target band invalidates a count assertion. Near addresses of size √N, varying both length-q address ranges can cost Θ(q√N) in the raw product alone, whereas the theorem promises only an O(√q) cyclic residue error. For the physical readout this diameter is not merely an upper risk. The checksum satisfies 0≤p(u)≤u by its recurrence (for odd u let p(u)=p((u−1)/2)). With x and y of order √N, comparing the two endpoints of a growing length-q y-interval at fixed x allows the checksum to cancel at most O(√N+q) of the signed difference, while the positive product difference is 432x(q−1). Consequently a full box has actual diameter Ω(q√N) for growing q=O(√N). A complete-box proof cannot simply place its whole box in the earlier O(N^(1/3)) bracket. No absolute residual contraction is inferred.

## 2. A sharp 6/8 obstruction that has three equal-status strands

For three integer strands x,y,z constrained by x+y+z=0, their symmetric half-energy equals
```text
(x²+y²+z²)/2 = x²+xy+y² = Norm(x−yω),   ω²+ω+1=0.
```
This exact same Eisenstein field enters the analytic QRH source, but only the elementary norm algebra is used. Any interchange of strands, regardless of the amount of time or pairwise hand swaps, leaves this energy unchanged. This does not claim that energy is a braid invariant or that it captures the knot; it says that if we only retain this symmetric energy, interchanging strands alone cannot cure its arithmetic obstruction.

If x,y are not both zero, let t=min(v_2(x),v_2(y)) with v_2(0)=∞. Removing 2^t makes at least one argument odd; the quadratic form is 1 on all three nonzero pairs in F₂². Thus v_2(N)=2t. Conversely every odd residue class modulo 2^k occurs in the norm: take y=1,x=2u, so N=1+2(u+2u²). The map u↦u+2u² is a permutation modulo 2^{k−1}, because its difference factors as (u−v)(1+2(u+v)), with the second factor odd. Scaling both inputs covers all nonzero residues with even valuation, and (0,0) supplies zero. Therefore the complete image modulo 8 is
```text
{0,1,3,4,5,7},
```
missing exactly 2 and 6. This is not a failure of a sampling experiment: it holds over all integer inputs and survives arbitrary hand permutation. The odd readout 1+2N misses 5 mod 8. A genuine second readout offset of +2, where physically realized, does remove the local obstruction: N and N+1 jointly cover every dyadic residue, so the two odd readouts cover all odd classes at each power of two. Merely relabeling the strands cannot act as that offset.

There is also a rigorous **limit to a phase-only 8/8 theorem**. Even where c is a unit modulo 8, legal adversarial unary functions F(x)=cx² and G(y)=cy² give exactly the above six-element image (multiplied by c). Thus an estimate uniform in the unary phases cannot itself prove exact residue surjectivity. One must use either the actual checksum structure or a physically distinct control. The shifted norm example is an arithmetic model, not a claim that these entire norm rays are tournament realizable.

Replay: `python check_qrh_dyadic.py` checks conductor subgroup sums and full norm residues through q=256, using only integer arithmetic.

## 3. Two total-sum sectors give a literal arithmetic 8/8

An alternative normalized signed readout is the symmetric pair interaction
`D(x,y,z)=xy+yz+zx`, for three integers with total S. For any q=2^k, let s=S/3 modulo q, which is well-defined because 3 is odd. Recenter x,y,z by s. Their centered values sum to zero modulo q, so the exact relation is

```text
D(x,y,z) ≡ S²/3 − Norm((x−s)−(y−s)ω)  (mod q).
```

Conversely arbitrary x−s,y−s occur as residues for integer triples of actual sum S, by taking the third coordinate S−x−y. Thus the entire fixed-S image is one translated Eisenstein norm image. At S=0 its mod-8 image is {0,1,3,4,5,7}, just as for the norm; at S=1 it is {0,2,3,4,6,7}. These achieve **8/8 in union**, and the same two sectors hit every residue modulo *every* 2^k: the first supplies every odd residue, and the second, an odd translate of the negated norm, supplies every even residue. This isolates the difference between swapping hand assignments (preserving S) and changing the common anchor (changing S). It is an exact arithmetic steering model, not a claim that every input triple or that change of S is implemented by the tournament gadget.

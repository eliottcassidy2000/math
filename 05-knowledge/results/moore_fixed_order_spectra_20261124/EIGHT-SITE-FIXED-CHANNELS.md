# An eight-site guard: three fixed linear +2 channels and the exponent-31 tail test

## 1. A genuinely different guarded core

In the [six-site module](GEWIRTZ-AND-DUAL-GUARDS.md), keep masks 177 and 181 on vertices 0–5, with their literal arc flip \((0,3)\). Append vertices 6,7 with \(6\to7\). For \(i=0,\ldots,5\) and \(j=6,7\), use \(i\to j\) by default and reverse it exactly when bit \(2i+j-6\) in the integer **189** is set. The common apex beats local sites **2 and 5** (outmask 36). Both eight-site cores now have
\[
 (h,Q,R,J)=(651,8808,2691,6117),\qquad
 H(Y+a)=2313,\ 2315 .
\]
The complement of the apex marking, outmask 219, keeps \(h,Q\), replaces \((R,J)\) by \((6768,2040)\), and reverses the sign of the gap. These are two marked versions of the same underlying core edit. Since h=651, this is not an ordinal source/sink padding of the old six-site h=29 core.

The discovery replay checks all \(2^{12}\) choices for the directions from old sites to the two new sites, fixing \(6\to7\) and the old flipped arc; it checks all 256 apex masks once the same \(h,Q\) tests pass. Of 4096 extensions, 98 preserve \(h,Q\), and 14 marked extensions meet all four guards and absolute gap 2. Twelve have \(h=29\) and the two above have \(h=651\). The scan is complete only **in this extension family**, with the old six-site pair fixed. It does not purport to census eight-site tournaments.

Use the [variable-tail transfer](VARIABLE-TAIL-SELECTOR.md) with \(B=2v,\ j=2u\), carrying \(h,Q,R,J\) as free guard coordinates. The new lower height is
\[
 F_8(v,u)=35232vu+13536v+5382u-1302(q(v)+q(u))+2313 .
\]
The opposite marking transposes \(v,u\). Again the proposed physical mates really differ by one arc, and by two Hamilton paths. In particular with binary payload \(1010\) and tail \(10\), the newly counted pairs are \(250233\to250235\) and \(217617\to217619\) on 15 sites. Independent full subset DP verifies both and all 15 tails of length at most four starting in 1 for the first marking.

## 2. Three *fixed* guarded channels give a universal linear inverse

Take six-site chart B with one guard, six-site chart A with two guards, and the new eight-site chart with tail 10. On a common even payload \(2v\), subtract their constants \(493,743,7695\), respectively, and divide each lower height by two. The resulting readings are
\[
 U=798v-29q,\qquad V=1471v-29q,\qquad N=24384v-651q .
\]
For all integral payloads the exact linear inverse is
\[
 \begin{split}
  v&=33457(V-U)-120(29N-651U),\\
  q&=-202(798v-U)+9(24384v-N).
 \end{split}
\]
Indeed \(V-U=673v,\ 29N-651U=187638v\) and \(33457\cdot673-120\cdot187638=1\). Also \(-202\cdot29+9\cdot651=1\). This gives exact joint inversion modulo **every** \(M\) from physical lower counts modulo \(2M\), and recovers the entire doubled-modulus state from heights modulo \(4M\). Unlike the previous three-channel recipe, no tail depends on the modulus and no nonlinear one-apex measurement is needed. Exact absolute compatibility still demands \(q=q(v)\) and equality of all three readings. The [prefix routing theorem](TWO-COORDINATE-RIGIDITY.md) continues to give every prescribed congruence state.

No pair of fixed *linear* measurements taken from the six- and eight-site modules with a fixed tail can be universal over all moduli. Two six-site measurements both lose the \(q\)-coefficient modulo 29. Any pair involving an eight-site measurement loses one entire row modulo 3: both its digit and slope coefficients are divisible by 3 for any fixed tail in either marking. Thus three are necessary within these two affine families. This is not a lower bound on all possible physical charts. The earlier four-channel result concerned a different fixed list; the new core changes that list itself.

## 3. Exact exponent-31 exclusion in the two variable-tail families

Let \(E_{31}=2^{30}(2^{31}-1)=2305843008139952128\). Neither \(E_{31}-1\) nor \(E_{31}+1\) occurs at either physical endpoint of **either of the two variable-tail guard families** (h=29 or h=651), even allowing both nonnegative parameters and both scalar orderings. This does not exclude either target as a general tournament height or a seam in another family.

The compiled replay is a scaled form of the two-tail inverse, using exact 64-bit arithmetic for this exponent and 128-bit comparisons for monotonicity. For code speed it computes q by concatenating precomputed 16-bit blocks, and tests the necessary height congruence modulo the h-coefficient before exact equality. The reversal-bucket choice directly requires
\[
  h\left(\left\lfloor\frac{\mathrm{lo}}{B}\right\rfloor-1\right)
       >(L+h)(B-1)
\]
at bucket size \(B\); this guarantees monotonicity for every prefix in the candidate interval (the last valid tested \(B\) or a singleton). The exact prefix and feasible-interval receipts are:

| Guard h | tested lower height | smaller parameter range | feasible intervals before residue filter | searched prefix buckets | hits |
|---:|---:|---|---:|---:|---:|
| 29 | \(E_{31}+1\) | 0 through 42,443,372 | 23,629,572 | 338,665,045 | 0 |
| 29 | \(E_{31}-1\) | same | 23,629,572 | 338,665,045 | 0 |
| 29 | \(E_{31}-3\) | same | 23,629,572 | 338,665,045 | 0 |
| 651 | \(E_{31}-1\) | 0 through 8,089,957 | 4,097,914 | 60,815,545 | 0 |

The other two lower-height queries for the latter family are ruled out modulo 3: all its lower heights are \(0\bmod3\). The program independently checks the block digit formula against the direct recurrence on 20,000 deterministic integers. A small-exponent check also agrees with the standalone JavaScript inverse.

Reproduce:

    c++ -O3 -std=c++20 eight-guard.cpp -o eight-guard
    ./eight-guard
    node eight-fixed-channels.js
    c++ -O3 -std=c++17 tail-inverse-large.cpp -o tail-inverse-large
    ./tail-inverse-large 31 1
    ./tail-inverse-large 31 -1
    ./tail-inverse-large 31 -3
    ./tail-inverse-large 31 -1 0 50000000 eight

The largest tests can take several minutes. No extrapolation to exponent 61 is made.

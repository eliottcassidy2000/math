---
id: THM-3729
title: "Rooted Pfaffian response and sign-root deletion average"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every real skew K
  and root vector u, the odd rooted-Pfaffian square sum is exactly
  det(I+K)u^T(I+K)^(-1)u.  For tournaments this identifies THM-1950's open
  strong base with an even/odd parity-energy inequality.  Averaging the odd
  energy over every incident sign root equals one half the sum of all
  vertex-deletion discriminants when n>=2.  Over a full switching class the
  Hamiltonian-path mean is n!/2^(n-1), and it dominates both the constant
  even-energy mean and the odd-energy mean by Hadamard bounds.  The fixed
  all-one-root inequality remains open.
source: root / 2026-08-22
audit: >
  PASS.  A separate Fraction/Gaussian-elimination implementation checked
  every tournament and incident sign root for 2<=n<=4, the n=2 normalization,
  every diagonal cofactor, both deletion/extension averages, and the
  five-vertex hostile.  The main companion exhausts all 33,867 tournaments
  through order six.  A separate switching-gauge implementation checks all
  1,099 labelled switching classes and all 33,866 tournaments of orders two
  through six, including both exact means and both envelopes.  Normal and
  optimized streams byte-match their frozen transcripts.
depends_on:
  - THM-468-tournament-determinant-floor
  - THM-1950-h-ge-disc-reduced-to-strongly-connected
  - THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay
script: 04-computation/tournament_rooted_pfaffian_response_thm3729.py
output: 05-knowledge/results/tournament_rooted_pfaffian_response_thm3729.out
script_sha256: 5dc79e980288341dab35cbec3475d53e095ace1c4b3a4b1cc6d9bad8d546ceda
output_sha256: 78b860d301b990ef74b3c43a8cfe2f6f668877fd9bca4733641573527508c73e
independent_audit_script: 04-computation/tournament_rooted_pfaffian_deletion_average_independent_audit_thm3729.py
independent_audit_output: 05-knowledge/results/tournament_rooted_pfaffian_deletion_average_independent_audit_thm3729.out
independent_audit_script_sha256: 63bac3ee2897220416e0336359d81a9210164e64a979b4855a9b4d727e6bf7df
independent_audit_output_sha256: b1a2a701d11f43ca6f6ad9a28e9470772bbafbe80726a951e7886a15deba3be3
switching_mean_script: 04-computation/tournament_switching_mean_envelope_thm3729.py
switching_mean_output: 05-knowledge/results/tournament_switching_mean_envelope_thm3729.out
switching_mean_script_sha256: ccc630cd5602273bbfcea184d2238b29fb668905a9754826e90339c62de49d4e
switching_mean_output_sha256: 46ed14941c43e0866dad45ab52a1344ba14261f5d8b3d2fcc9f418b53a2c7198
hash_basis: raw LF bytes
---

# THM-3729 -- rooted Pfaffian energy is the exact inverse response

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Rooted identity

Let `K` be a real `n x n` skew matrix, put `M=I+K`, and let
`u in R^n`.  For an odd subset `S` define

```text
p_u(S)=Pf([[K[S],u[S]],[-u[S]^T,0]]).
```

Then

```text
sum_(S odd) p_u(S)^2
 =det(I+K) u^T(I+K)^(-1)u.                            (1)
```

Moreover

```text
s_u:=u^T(I+K)^(-1)u=||(I+K)^(-1)u||^2>=0.             (2)
```

### Proof

A real skew matrix has purely imaginary spectrum, so `M=I+K` is
nonsingular.  Augment `K` by the root:

```text
Khat_u=[[K,u],[-u^T,0]].
```

The Schur complement gives

```text
det(I+Khat_u)=det(M)(1+u^T M^(-1)u).                  (3)
```

The principal-minor expansion of the left side has only even skew minors.
Those not containing the root sum to `det(I+K)`; those containing it are
indexed uniquely by odd `S` and equal `p_u(S)^2`.  Subtracting
`det(I+K)` from (3) proves (1).

For (2), set `x=M^(-1)u`.  Since `u=Mx`,

```text
u^T M^(-1)u=u^T x=x^T(I-K)x=||x||^2,
```

because `x^T Kx=0`.

## 2. Tournament parity energies

For a tournament `T` with skew adjacency matrix `K`, define

```text
E_even(T)=2^(-(n-1)) sum_(S even) Pf(K[S])^2,
E_odd(T;u)=2^(-(n-1)) sum_(S odd) p_u(S)^2.            (4)
```

The principal-minor identity and (1) give

```text
E_even(T)=disc(T),
E_odd(T;u)=disc(T) u^T(I+K)^(-1)u.                    (5)
```

At the fixed all-one root, THM-1950's invariant is therefore

```text
P(T)=max(E_even(T),E_odd(T;1)).                        (6)
```

Its residual strong-tournament base is exactly

```text
H(C)>=max(E_even(C),E_odd(C;1)).                       (7)
```

Equation (7) is still **OPEN**.  If `T` is regular, `K1=0`, so
`E_odd(T;1)=n disc(T)`.

## 3. Sign roots and one-vertex extension

For `u in {+1,-1}^n`, `Khat_u` is the skew matrix of the tournament
`That_u` obtained by adjoining a vertex whose incident signs are `u`.
Taking account of the order-`n+1` normalization in (3),

```text
E_odd(T;u)=2 disc(That_u)-disc(T).                     (8)
```

Every rooted Pfaffian in this sign case is odd: it is a signed sum of an odd
number of perfect matchings.  There are `2^(n-1)` odd subsets.  Together
with THM-468 and (8), this yields

```text
E_odd(T;u) is an integer and E_odd(T;u)>=1.             (9)
```

Switching acts on the rooted object by

```text
(K,u) -> (DKD,Du).                                     (10)
```

The root vector is load-bearing: resetting `Du` to `1` generally changes
the energy.  The companion contains a five-vertex exact control whose
unnormalized even/odd numerators change from `(32,128)` to `(32,32)` when
the covariant root is discarded.

## 4. Exact sign-root deletion average

Assume `n>=2` and average over all `2^n` sign roots.  Since

```text
2^(-n) sum_u u_i u_j=delta_ij,
```

equation (5) gives

```text
2^(-n) sum_u E_odd(T;u)
 =disc(T) tr((I+K)^(-1)).                              (11)
```

Cramer's rule identifies the diagonal cofactor:

```text
((I+K)^(-1))_(v,v)
 =det(I+K[V\{v}])/det(I+K).
```

Because `T-v` has normalization `2^(n-2)`, (11) becomes

```text
2^(-n) sum_u E_odd(T;u)
 =1/2 sum_(v in V(T)) disc(T-v).                       (12)
```

Also

```text
(I+K)^(-1)=(I-K)(I-K^2)^(-1).
```

The second factor is symmetric and commutes with `K`, so the skew term has
trace zero.  Thus the middle member of (11) may equivalently be written
`disc(T)tr((I-K^2)^(-1))`.

Averaging (8) gives the extension form

```text
2^(-n) sum_u disc(That_u)
 =1/2 disc(T)+1/4 sum_v disc(T-v).                     (13)
```

The restriction `n>=2` avoids introducing an unmaintained
`disc(empty)` convention.

## 5. Switching-class mean envelope

Let `S(T)` be the labelled switching class

```text
S(T)={T^D:K(T^D)=DK(T)D}/(D~-D),       |S(T)|=2^(n-1).
```

Every ordering `pi=(v_1,...,v_n)` determines exactly one member of `S(T)`
in which `pi` is a directed Hamiltonian path.  Indeed, its consecutive arcs
impose the path-graph equations

```text
d_(v_i)d_(v_(i+1))K_(v_i,v_(i+1))=1,
```

which determine all signs `d_v` uniquely up to their common negation.
Double counting pairs `(T^D,pi)` therefore gives

```text
2^(-(n-1)) sum_(T' in S(T)) H(T')=n!/2^(n-1).        (14)
```

The covariant rooted action (10) gives
`E_odd(T^D;1)=E_odd(T;D1)`.  The two roots `u` and `-u` have the same
squared energy, so (12) is equivalently the switching-class identity

```text
2^(-(n-1)) sum_(T' in S(T))E_odd(T';1)
 =1/2 sum_v disc(T-v).                               (15)
```

These two exact means satisfy a useful all-scale envelope.  Hadamard's
inequality applied to `I+K` and to every deletion gives

```text
disc(T)<=n^(n/2)/2^(n-1),
1/2 sum_v disc(T-v)
 <=n(n-1)^((n-1)/2)/2^(n-1).                        (16)
```

For every integer `m>=1`, pairing the factors `j` and `m+1-j` proves
`m!>=m^(m/2)`; when `m` is odd the unpaired middle factor is at least
`sqrt(m)`.  Apply this once with `m=n` and once with `m=n-1` in (16):

```text
average_(T' in S(T)) H(T') >= disc(T),
average_(T' in S(T)) H(T')
 >=average_(T' in S(T)) E_odd(T';1).                 (17)
```

This is a mean comparison, not an average of
`max(E_even,E_odd)`, and it is not pointwise in a switched representative.
Thus it identifies a class-level reservoir but does not prove the open
strong-tournament inequality (7).

## 6. Failed raw-square shortcut

One tempting route to (7) first bounds each rooted Pfaffian by the
Hamiltonian-path count of its induced subtournament and then asks for

```text
sum_(S odd) H(T[S])^2 <= 2^(n-1)H(T).                  (18)
```

The second step is false.  In the five-vertex tournament with arcs

```text
1->0, 2->0, 3->0, 0->4, 2->1,
1->3, 4->1, 3->2, 4->2, 4->3,
```

one has

```text
H(T)=15,
sum_(S even)H(T[S])^2=104,
sum_(S odd)H(T[S])^2=272>16H(T)=240,
```

whereas the signed Pfaffian numerators are only `(32,128)`.  The termwise
comparison survives on this witness; the raw-square relaxation (14) is the
first failed implication.  The exact signed/rooted energies, not unsigned
subpath squares, are the strongest survivor.

## 7. Scope and reproduction

The deletion average (12) controls the whole incident-root sign fibre.  It
does not control the fixed `u=1` response in (7), and it proves no
`H>=disc`, LRC, FC, HFC, or JC statement.

```bash
python3 -B 04-computation/tournament_rooted_pfaffian_response_thm3729.py
python3 -B -O 04-computation/tournament_rooted_pfaffian_response_thm3729.py
python3 -B 04-computation/tournament_rooted_pfaffian_deletion_average_independent_audit_thm3729.py
python3 -B -O 04-computation/tournament_rooted_pfaffian_deletion_average_independent_audit_thm3729.py
python3 -B 04-computation/tournament_switching_mean_envelope_thm3729.py
python3 -B -O 04-computation/tournament_switching_mean_envelope_thm3729.py
```

Each normal/optimized pair must agree byte for byte with its frozen
transcript.  **QED.**

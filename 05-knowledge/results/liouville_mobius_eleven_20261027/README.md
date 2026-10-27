# Liouville transport, Möbius guards, and the several roles of eleven

2026-10-27. **PROVED:** elementary identities and finite-prime continuation
theorems below. **FINITE-EXACT:** new arithmetic readouts on the actual 81/27
path-complement pair. **CITED / INHERITED:** the explicitly identified older
results. All-size strong-atom classification, economical closure under unnamed
re-marking, Collatz coverage and LRC(14) remain **OPEN**. No priority claim.

The useful duality is

\[
\mu(n)=\lambda(n)1_{n\text{ squarefree}},\qquad
\mu(ab)=\mu(a)\mu(b)1_{\gcd(a,b)=1}.
\]

Liouville transports a sign through every multiplication. Möbius also enforces
an exclusion guard. This is a literal model of the repo's distinction between
transport and realization, with a completely classifiable continuation carrier.

For a bank of k>=1 distinct primes, and arbitrary multiplication continuations:

| Observed output | Minimal deterministic states | Linear response rank |
|---|---:|---:|
| Liouville alone | 2 signs | 1 |
| Möbius alone | 2^k squarefree supports + 1 failure | 2^k |
| Both | 2^k squarefree supports + 2 failure signs | 2^k+1 |

These ranks concern output responses, without an additional constant mass
coordinate. Deterministic state counts are not numbers of memory bits.
For {2,3,11} the respective state counts are **2, 9, 10**. Keeping both functions
requires preserving the sign even after the squarefree guard has failed.
See [THEORY.md, D1–D4](THEORY.md) for proofs and a distribution-level hostile.

The new tournament test uses the *counted*, uncoloured legal-path operator
`Lf(T)=sum_(P Hamilton in T) f(tau_P T)`, retaining each path's multiplicity.
For the established twelve-vertex hostile pair:

| source | Aut | H | mu(H) | lambda(H) | L(mu composed with H) | L(lambda composed with H) |
|---|---:|---:|---:|---:|---:|---:|
| C3[A,A,A] | 81 | 462915 | 0 | +1 | -2916 | -3321 |
| C3[A,A,A opposite] | 27 | 462915 | 0 | +1 | 297 | -1323 |

Here A=C3->1. Both arithmetic responses distinguish the pair; Möbius reverses
sign. This does not show either response classifies all strong atoms. The
coloured black/blue incidence and joint ports are not replaced by these scalars.
More strongly, both sources occupy the same joint arithmetic failure state:
the exact multiplication carrier is therefore **not closed under complement**.

[ELEVEN.md](ELEVEN.md) recovers the September odd-square theorem, June Euler
ghosts, July three-prime amplifier, October square-zero inventories, level-11
eta work, signed Collatz cycles, and the sharp 4p symmetry threshold. The main
new integrations are:

- A full-rank exclusion response tensor, with bounded marginals and even the
  complete sign distribution insufficient for later joint guards.
- Exact square-dilation transforms between the two summatory functions and
  between their ordinary Euler products.
- The threshold 44=4*11 corresponds to a degree-1 to degree-11 decoration
  change. Modulus11 and Frobenius degree11 play different roles.
- An explicit {2,3,11} circle-projector corner with Möbius value -1/8712;
  its construction is unrestricted interval data, not an LRC source theorem.
- Replacing Möbius by Liouville in inversion passes squarefree-divisor tests
  such as degree11, but fails at degree121 and in the level-11 ghost inverse.

Reproduce from the repo root:

```text
python -B 05-knowledge/results/liouville_mobius_eleven_20261027/verify.py --replay
```

This reruns normal and optimized Python, compiles/runs the complete C++ hostile
enumeration, and audits every saved H-bin with an independent linear sieve.
Independent recursive path counts check both sources and sixteen partners;
the full partner enumeration is one implementation, explicitly so labelled.
[SOURCES.json](SOURCES.json) pins recovered files and their read scopes;
[SHA256SUMS](SHA256SUMS) seals this packet. No executable is committed.

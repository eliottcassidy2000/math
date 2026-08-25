# Source pin: Sun's 2-4-6-8 mixed-binomial conjecture

**Freshness: checked 2026-08-24.** This file separates statement provenance
and fresh counterexample intake from the repo's independent proofs.

## Primary statement and sequence

- [Sun's MathOverflow question 323541](https://mathoverflow.net/questions/323541),
  asked 2019-02-19, fixes the literal domain `w,x,y,z>=2`, notes
  `1/2+1/4+1/6+1/8=25/24`, and records a search through `2*10^12`.
- [OEIS A306477](https://oeis.org/A306477) defines the canonical ordered
  representation count, using one shifted representative for each
  higher-degree zero fibre, and records examples and the same conjecture.
- The MathOverflow thread identifies `20 mod 33` as the unique
  minimum-density class, with probability `16/1089`; the repo independently
  rederives this row.

These sources state and contextualize the conjecture. They do not prove it or
the counterexample.

## Fresh counterexample intake

[Adamczewski's public gist](https://gist.github.com/tadamcz/0c578c8b2b3fb92fe8584bc0725187e3)
reports `896315812331399` and an exact verifier. It is a third-party
2026-08-24 computational source. The repo makes no priority claim from it and
does not use its code as a proof dependency.

## Repo consumers and exclusions

- [THM-4026](../../01-canon/theorems/THM-4026-sun-two-four-six-eight-binomial-counterexample.md)
  is the independent exact counterexample proof.
- [THM-4027](../../01-canon/theorems/THM-4027-sun-two-four-six-eight-universal-modular-solubility.md)
  proves universal local solubility.
- [THM-4028](../../01-canon/theorems/THM-4028-sun-two-four-six-eight-average-order-criticality.md)
  proves the summatory capacity law.
- [HYP-1953](../hypotheses/HYP-1953-additive-basis-normal-form-spectrum.md)
  owns the broader additive-basis frontier.

None proves leastness; priority; finiteness, infinitude, or density of the
exceptional set; a pointwise circle-method law; or an LRC/JC/Rule-30
consequence. A 2026 ResearchGate item linked from OEIS is **RADAR ONLY** and
is not used by any theorem.

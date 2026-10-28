# Audit scope and reproduction

Same-agent mathematical review and alternate-path exact checks.
No independent-agent audit was commissioned. Proofs are in [THEORY.md](THEORY.md).
No conjecture is promoted from these finite checks.

Run from the repository root (Python3; SymPy needed only for prime search):

```powershell
python3 -B 05-knowledge/results/relation_tower_frobenius_ports_20261028/check_fields.py
python3 -B 05-knowledge/results/relation_tower_frobenius_ports_20261028/prime_bank.py
python3 -B 05-knowledge/results/relation_tower_frobenius_ports_20261028/check_signed.py
python3 -B 05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py --max-order 9
python3 -B 05-knowledge/results/relation_tower_frobenius_ports_20261028/audit.py
```

The fourth command took about5 minutes in this run. It uses substantial
memory for1716608 sparse pointed rows and922030 cached deletion parents.
Default maximum order is8; pass9 to reproduce the saved bank. The
executable rejects bounds outside1..9. The all-order mathematical
construction is distinct from this explicit execution cap.

| Receipt | Exact universe / checks | Independent or hostile control |
|---|---|---|
|[fields.json](fields.json)|F4,F16,F64,F256; full alternating-form constraint spaces, all normal-orbit pairs, every three-port phase tuple; all698027 new triple relations;699788 active checks|F64 ambient32768-form exhaustion; F64/F256 full-period collisions; sharp bit-probe witnesses|
|[prime_bank.json](prime_bank.json)|2047 masks over11 coordinate primes,4085 recursive prime certificates|Alternate verifier reconstructs factors and checks Lucas conditions; all698027 integer triple products checked as actual squares|
|[signed.json](signed.json)|8428 signed systems: every s=+/-1 and w in [-3,3] on3-vertex path/triangle and4-vertex path/star; labels0..4|Brute-force distinct integer labelings; odd square-cycle integrality failure; exact eight graceful labelings of the5-path|
|[root_switch.json](root_switch.json)|Every parent attachment orbit and every rerooting through9;6292804 active checks|No child canon/count inputs; exact Burnside check; full independent automorphism/Hamilton audits through8;81/27 hostile|
|[audit.json](audit.json)|900180 alternate checks|Entire native metagraph tc histograms through8, all labelled tournaments/permutations through5, bank type checks,15-cell125-kernel hostile|

The rank statements about Q use its verified identical-row blocks,
not numerical floating row reduction. Nonzero orbit-size entries are
odd, so each disjoint block has rank1 over F2. Q^2=0 at even child
orders does not make Q zero. The node-count formula and characteristic
distinction are inherited theorem consequences, checked at every state.

At order9, the full representative bank and automorphism-order histogram
are computed. Prime-class count and full Hamilton-path histogram are
**not** computed there; `prime_classes:null` and an empty tc histogram
explicitly encode that limit. The independent native-metagraph comparison
stops at8. Actual isomorphism lookup first refines colours and then uses
exact backtracking; equal colours alone are never accepted as an isomorphism.

The prime certificate proof is the classical full-factor Lucas criterion:
if the recorded a has a^(p-1)=1 mod p and
gcd(a^((p-1)/q)-1,p)=1 for every certified prime q dividing p-1, then
modulo each prime divisor l of p its order is divisible by p-1. Thus
p-1 divides l-1, forcing l>=p and p prime. Recursive factors terminate
at2. The search's primality guesses are not the proof. The search found
the bank; no uniform bound on its continuation is inferred.

Proof review checked the following boundaries explicitly: T2 normality
is sufficient, not necessary pair by pair; T3 metadata is fixed and
gluing disconnected pieces needs a relative phase; T4 covers fixed
edge equations and retains global collision tests; T5 preserves the
attachment mask through the same parent isomorphism; T6 relies on
nonuniform equitability and exact internal chart groups. The15-point
geometric bank and15-vertex composite tournament are different objects.

The initial implementation's empty-parent colour refinement needed an
explicit n=0 return; it was fixed before the first census. This was a
program boundary error, not a mathematical count correction. No old
canon claim was changed in this packet. The incoming THM-4840 repairs
at5edadbae58 were read before reuse and are preserved.

Normal and optimized runs of the field, signed-label and alternate-audit
scripts give byte-identical receipts; see [optimized_checks.json](optimized_checks.json).
A separate optimized root-switch replay through8 matches every reported
semantic field and the81/27 witness; execution times are excluded from
that comparison. The order-nine run was not redundantly repeated.
All local Markdown links resolve; [final_validation.json](final_validation.json).
The repository's agent-facing documentation checker also passed.

The concurrent Session86 default checker was replayed after its arrival
at0884018279. It reconstructs ordinary counts through50 and the aggregate
prime cycle index through20 (129 monomials), with independent pair-action
and exponent checks. The [incoming receipt](incoming_a000568.json) agrees
with our constructive class bank through9 and prime counts through8.
It does not independently enumerate those large prime classes; the
prime-index recovery/recomposition is explicitly a structural round trip.

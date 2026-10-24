# Collatz 27: return words, arithmetic cylinders, and the mod-7 boundary

Status: elementary PROVED propositions below; FINITE-EXACT orbit and bounded test data. No novelty claim for parity cylinders, CRT, affine Collatz words, or Paley affine actions. No universal Collatz or square-sum graph theorem follows.

## Inheritance and guards

Read AGENTS.md, START-HERE, RESEARCH-PROTOCOL and ACTIVE-GUARDRAILS in the supplied source checkout. The closest mechanism is the exact affine cylinder, with THM-4512's corrected extra bit for an exact odd-valuation word. The canonical hostile is the odd-step sequence 27 -> 41 -> 31: the child has descended but the original source has not. The corrected near miss is HYP-9282, REFUTED: residue-channel persistence did not establish square-root main terms; octave dependence invalidated its independence model. The least-used useful sidecar here is the endpoint affine carry modulo 7. The finite concept board is clock / dyadic prefix / source debt / affine endpoint / quadratic character.

Source lineage (original relative paths, raw-byte hashes in manifest): procgen_family27_20260926_long_orbit_families.md; entry_20260927_families.md and recursive.md; THM-4512-coefficient-descent-classes-one-member.md; collatz_artin_corrections_20260927_audit.md; collatz_affine_guarded_lifts_20261004.md; THM-4532-tournament-clocks-two-sheet-clocks-all-odd-tournaments-and-the-syracuse-clock.md; HYP-9282-collatz-parity-square-root-residue-biases.md. These are inherited results, not newly established large scans. In particular the old 2^32 family scan is not repeated.

## 1. Three clocks, kept distinct

C(n)=3n+1 for odd n and n/2 for even n. T(n)=(3n+1)/2 for odd n and n/2 for even n. U(n)=oddpart(3n+1) for odd n. Starting from 27, exact direct computations give:

- C: total stopping 111; first value below 27 at time 96; peak 9232 at time 77.
- T: total stopping 70; first value below 27 is 23 at time 59; peak 4616 at time 45.
- U: total stopping 41; first value below 27 is 23 at time 37; peak odd value 3077 at time 32.

The full 41-entry U valuation word, whose sum is 70, and all 42 odd orbit values are retained in normal.json. C and U are computed directly, independently of the T affine recurrence. All clocks stop on their first visit to 1 for these statistics. Fixed-word statements below use the everywhere-defined T map, not an artificially absorbing map at 1.

## 2. Fixed-prefix and inverse-cylinder theorem

For a T parity word w=(e_0,...,e_(k-1)), let s_j be the number of ones before time j. Put b_0=0 and b_(j+1)=3b_j+2^j when e_j=1, b_(j+1)=b_j otherwise. There is exactly one source class r mod 2^k realizing w. On that class,

T^j(n)=(3^s_j n+b_j)/2^j, 0<=j<=k.

Proof: an already fixed j-letter word has one class mod 2^j. Its two lifts mod 2^(j+1) give endpoints differing by odd 3^s_j, so exactly one has each possible next parity. This proves uniqueness and existence inductively. The affine formula follows by substitution. In particular,

T^k(r+2^k t)=T^k(r)+3^s_k t.

The inverse correspondence is y -> (2^k y-b_k)/3^s_k, with the endpoint progression y=T^k(r)+3^s_k t. Positivity remains a separate guard. For the 27 representatives below, 0<27<2^k, hence every positive source in its cylinder has t>=0.

An exact U valuation word with r odd steps and A total divisions requires n congruent to (2^A-b)3^(-r) modulo 2^(A+1). The coarser modulus 2^A allows additional final halvings, as THM-4512 explicitly corrects. For 27's full U word, the exact cylinder is n=27+2^71 t, endpoint 1+2*3^41 t. Its coarse T cylinder is n=27+2^70 t, endpoint 1+3^41 t. Hostile: t=1 in the latter has even endpoint 36472996377170786404 and therefore cannot realize the prescribed exact U word.

## 3. A source-preserving infinite 27-glide family

For every integer t>=0,

n=27+2^59 t has first T descent exactly at step 59,
T^59(n)=23+3^37 t,
first C descent is exactly step 96.

Proof: the exact finite word check verifies 3^s_j>2^j for every 1<=j<59. Because b_j>=0, every lifted positive source remains greater than itself at all these times. At time 59 the multiplier is 3^37<2^59 and the lifted difference is

T^59(n)-n=-4+(3^37-2^59)t<0.

The word has 37 odd steps, so the ordinary clock is 59+37=96. Inserted ordinary odd-step peaks cannot be earlier descents: 3x+1>x. The finite prefix inequalities are an explicit proof certificate, not a sampled extrapolation.

Exact first-descent affine identity:
2^59 T^59(n)=450283905890997363 n+1100931843921811423.

Exact full-word affine identity on n=27 mod2^70:
2^70 T^70(n)=36472996377170786403 n+195820718533800070543.

This second family shares 27's first 70 parity steps, not its total stopping time: for every t>0, T^70(n)=1+3^41 t>1. Neither family is asserted to enter 27's backward tree. Indeed the shortcut predecessors of 27 are only its doubling chain, since a node divisible by 3 has no odd predecessor and doubling preserves divisibility by 3.

## 4. Exact CRT independence, with its scale specified

For any fixed word length k and any fixed odd q, every word occurs with every source residue mod q in exactly one residue class mod q*2^k. Therefore a complete block of q*2^k consecutive integers contains every pair (word, residue) once. Along one word cylinder, any q consecutive lift parameters t hit each residue exactly once. Within [1,X], each such joint class count is X/(q*2^k)+O(1), with absolute error at most 1.

Consequently no mod-7 quadratic-character class can rule out any finite parity prefix. The 27 first-descent cylinder contains infinitely many sources in each of the seven mod-7 classes, with exact equal counts in each complete block of seven lift parameters. This is about integer classes, not equal prime counts and not an independence assertion for growing k or incomplete short intervals. The prior Artin audit specifically warns against using CRT at a modulus much larger than the observation window. No independence of the initial and terminal residues is asserted; they have the affine relation below.

## 5. Positive transport law: Paley orientation, not autonomous Collatz

For q=7, let chi be the quadratic character, with chi(0)=0. A fixed word induces the residue permutation

f_w(x)=3^s 2^(-k) x+b_k 2^(-k) mod7.

Since chi(2)=1 and chi(3)=-1, for every distinct x,y,

chi(f_w(y)-f_w(x))=(-1)^s chi(y-x).

Thus on the Paley tournament x -> y iff chi(y-x)=1, fixed-word transport is an automorphism for even s and an anti-automorphism for odd s. This is a lawful pairwise observable: it remembers the common word, discards integer magnitude, and preserves or reverses only the difference-character predicate. It gives no stopping-time estimate.

For the two concrete 27 words:
- first descent (k,s)=(59,37): f(x)=1-x mod7;
- complete baseline orbit word (70,41): f(x)=-x mod7.

The latter additive carry is divisible by 7. In particular its endpoint character is minus its source character on nonzero residues. This specific character flip is a word identity, not a property of every 70-step orbit. Minimal hostile to an autonomous quotient: 27=34 mod7, but T(27)=41=6 mod7 whereas T(34)=17=3 mod7. The dyadic guard cannot be deleted.

Quadratic reciprocity gives (-7/p)=(p/7) for odd primes p!=7, so the s7 prime split and this finite residue transport can be expressed with the same character. This does not import the s7 supercongruence or its prime-specific conclusions into Collatz. No novelty or external literature attribution is being claimed for this elementary reciprocity identity.

## 6. What prior reset work actually supplies

Inherited entry_20260927_families constructs n=a*8^k-5 with exact initial valuation word (1,2)^k and even-boundary evolution n+5 multiplied by 9/8. Its invariant 2v2(n+5)+3v3(n+5) and prime-to-6 part are phase ledgers, not global descent ranks. Choosing a=2(47*4^h+7)/3^(2k+1), with the displayed integrality guard, makes an actual 2k+1-step odd path to 47. All k>=1 members except 27 then first descend at that step; 27 is the exceptional unpaid source and needs its finite suffix. The recursive note's fixed-coefficient family 4*8^k-5 instead has unpaid growth beyond 2k+4+v2(k) odd steps. It does not inherit the completion theorem merely by sharing its prefix.

The new progression statement in section 3 is a different simple specialization: it retains a fixed source comparison and proves all-height glide stability on one dyadic cylinder. It is not a claim that the inherited infinite completed family, this dyadic family, and the tree Pred*(3077) coincide.

## 7. Boundary with the user's G27

The companion worker identifies G27 as the square-sum graph on vertices 1..27, not a Paley7 graph. No map from its Hamiltonian paths to these Collatz words has been found. Therefore the shared numeral 27 is not a mathematical bridge. The useful common methodology is retaining the growth boundary/endpoint sidecar; the actual mechanisms remain separate.

## Reproduction and audit ledger

python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json

Both modes agree byte for byte. Explicit raising checks, no assert-dependent validation. Exhaustive universe: all 2046 words of lengths 1..10, all seven CRT lifts (14322 prefix checks), all 42 ordered distinct residue pairs for each word (85932 Paley checks). Return-word lifts: t=0..127 plus 2^100+9 for each of k59 and k70. Exact U cylinder: t=0..31. Hostiles: missing final bit, same mod7 source with different successor, inherited child-vs-parent reset debt, and false common total-stopping inference.

Development correction: the first reset hostile accidentally used shortcut T(41) where the inherited example uses odd U(41); this check failed in both modes. It was repaired to the explicit equality T^2(41)=31. No mathematical conclusion relied on the failed clock substitution.

The manifest closes the artifact file set and snapshots the precise inherited inputs. No repository files were changed by this worker; parent publication and independent audit remain separate.

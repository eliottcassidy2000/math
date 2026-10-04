# Modular multiplication, sine/cosine, Fibonacci registers, and guarded routes

**Status:** PROVED elementary identities and theorems, with separately labelled
FINITE-EXACT audits. Design recommendations are OPEN work. This does not prove
Collatz or LRC(14). Date: 2026-10-03. Author: Codex research team.

## 0. Inheritance and working board

The user's anchor is a construction carrying its own route home. The niche is
the information lost by folding modular multiplication tables; the wildcard is
the actual three-colour Zeckendorf diagonal convention. The six live concepts
are modular reflection, signed atoms, guarded route words, tournament blocks,
Fibonacci carries, and first descent.

Closest mechanisms recovered before deriving:

- [Tiling/modular atoms](tiling_modular_atoms_20261003.md): the triangular
  multiplication wedge, A/B centre blocks, and typed p/q maps.
- [Atom route memory](collatz_atom_route_memory_20261003.md): signed siblings,
  fixed-head ternary guard, first-hit tournament code, and safe internal growth.
- [THM-2218, labelled guard-hole Fourier and signed lift energy](../../01-canon/theorems/THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy.md):
  the LRC mechanism is labelled cyclic correlation; magnitudes lose alignment.
  We transfer that mechanism, not its LRC conclusion.
- [THM-1450, odd is sin is skew](../../01-canon/theorems/THM-1450-odd-is-sin-is-skew-the-heptad-spectrum.md):
  skew tournament spectra and reflection parity. The theorem's refuted
  identification of a tournament polynomial with sin(7 theta) is the corrected
  near miss. Different involutions must still be typed separately.
- [Duck/Zeckendorf](duck_zeckendorf_20260925.md),
  [reset colours](reset_20260926_colours.md), and
  [crossroads colours](crossroads_crossing_20260926_colour.md): the diagonal
  carry, a neutral fourth charge, and the confirmed extra-unit stripe with
  unresolved global continuation.
- [THM-4528, golden holonomy](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md)
  is routed through its [audited result note](collatz_golden_holonomy_20261001.md):
  no-adjacent-ones parity words and the polynomial route-tail certificate.
- [Recursive entry family](entry_20260927_recursive.md): the -5 return macro.
  Its retained original-source bound, not an average slope, proves first descent.

The least-used useful sidecars are Fourier phase, the second Fibonacci
register, the first-hit root pointer, and the original-source size budget.
Canonical hostiles below test what happens when each is erased.

## 1. Every positive modulus fits the wedge, but parity has two meanings

Use the explicit convention

    F_n=(x,y,z)=(n,ceil(n/2),2n+1),
    even modulus=2x,   neighbouring odd modulus=z.

Thus x is the *half-size* of the even table. Its own parity alternates the A/B
families. Both moduli use the triangle D_x={(a,b):1<=a<=b<=x}; zero axes are
known separately. Reflect each coordinate to 0..x, preserving its sign, then
transpose if necessary. The signed product reconstructs every residue.

For m>=1 replace a modular product by its faithful character value

    K_m(a,b)=exp(2 pi i ab/m).

This is the unnormalized Fourier matrix. Reflecting either coordinate
conjugates K; transposition leaves K unchanged. Cosine is the reflection-even
part and sine the reflection-odd part. Multiplication by a sends frequency b
to ab: it permutes characters iff gcd(a,m)=1 and otherwise has image size
m/gcd(a,m). In particular doubling is invertible for odd m and collapses pairs
for even m. Moving from m to 2m is a different operation from doubling inside
one fixed residue ring.

The exact addition/multiplication dictionary is

    K_m(a+b,c)=K_m(a,c)K_m(b,c),
    K_m(ab,c)=K_m(a,bc).

Addition becomes multiplication of phases; integer multiplication becomes
frequency substitution. Both identities hold at every modulus. The x/z
choice changes the table geometry and invertibility of2, rather than
assigning addition exclusively to one modulus parity.

Let Rf(a)=f(-a) on real functions on Z/m. Fixed residues satisfy 2a=0.
There are two for even m (0 and m/2), one for odd m (0). Counting fixed points
and two-element orbits proves:

| modulus | reflection-even dimension | reflection-odd dimension | fixed modes |
|---|---:|---:|---|
| 2x | x+1 | x-1 | constant, Nyquist (-1)^a |
| 2x+1 | x+1 | x | constant |

Both table parities contain cosine and sine. What changes is the unpaired
Nyquist mode. The continuous even/cosine and odd/sine convention agrees with
[DLMF 1.8](https://dlmf.nist.gov/1.8); the finite counts above follow directly
from the involution and require no analytic convergence assertion.

For an odd modulus m>=3, its central nonzero-index 2x2 block is

    [[d,-d],[-d,d]] mod m,   d=4^(-1) mod m.

Indeed the two central indices are -1/2 and 1/2 modulo m. This recovers
d=N at m=4N-1 and d=3N+1 at m=4N+1. Writing theta=2 pi d/m, its character
block is

    cos(theta) [[1,1],[1,1]] + i sin(theta) [[1,-1],[-1,1]].

The symmetric and contrast vectors (1,1),(1,-1) have eigenvalues
2 cos(theta), 2i sin(theta). This is an exact local sine/cosine interpretation
of the user's centre observation, including composite odd moduli.

## 2. The sine construction recovers the actual tournament blocks

The symmetric kernel K_m itself supplies no orientation: even its imaginary
part is symmetric in a,b. Instead use a difference observable on vertices Z/m:

    S_ab = sign sin(2 pi (b-a)/m),  a!=b;   S_aa=0.

For m=2x+1, S_ab is +1 exactly when b-a lies in {1,...,x}. There are no
off-diagonal ties. With the chosen positive cyclic generator as orientation
gauge, this is the regular cyclic tournament R_(2x+1). Its preserved target is
the cyclic pair relation; replacing labels by its isomorphism class loses the
chosen origin/generator. The route grammar only needs its order.

The character vector v_r(a)=exp(2 pi ira/m) is an eigenvector with

    lambda_r=2i sum_(d=1)^x sin(2 pi rd/m).

Proof: substitute d=b-a in Sv and pair d with -d. In particular the spectrum
is imaginary, with lambda_0=0. These are precisely the R_(2j+3) blocks in the
existing route code, with x=j+1. Adding two vertices increments j by one;
the ternary guard determines when that graph edit preserves a numeric route.

For even m=2x, the antipodal difference x equals its own negative. Every
translation-invariant skew pair observable therefore vanishes there. A full
tournament requires another ordering choice, breaking that symmetry. The
earlier transitive-clone doubling supplies an order within each clone fibre;
it is not a translation-invariant sine tournament on an even cycle.

There is a second, different intrinsic relation for odd primes p: the quadratic
character of b-a is skew for p=3 mod4 (Paley tournament), symmetric for p=1 mod4
(Paley graph), since chi(-1)=(-1)^((p-1)/2). Thus A/B also controls this
reflection type. This prime relation is not the cyclic half-circle relation;
composite moduli introduce nonunit zero values. It cannot orient all tables.

## 3. A guard mask uses the information discarded by cosine

The complete derivation, compiler and exact audit are in
[Fourier guards](collatz_fourier_guards_20261003.md).
For a compressed prefix (j_1,...,j_p), let

    a_i in {2j_i+1,2j_i+2},  A=sum a_i,
    S=sum_(i=1)^p 3^(p-i) 2^(a_1+...+a_(i-1)).

Its target h is legal exactly when h=2^(-A)S modulo 3^p; each exponent word
supplies one residue. The 2^p residues are distinct: backwards division by 3
forces the parity of the next exponent from the current target residue.
The sibling parameter is a bijection modulo 3^p, so a fixed compressed prefix
has exactly 2^p legal phases in that parameter. Keep positivity and the
first-hit root condition separately: a valid padded 1->1 step is not a
first-arrival certificate.

For head (0,0), terminal 1 and h_j=(4^(j+1)-1)/3, the phases are

    j mod9 in {0,3,4,7}.

They include 75->113->85->1, 151->227->341->1,
19417->14563->21845->1, and 621377->466033->349525->1.
Phase0 is legal for j=9,18,...; j=0 would pad the root and must be removed.
The prior rule fixing both actual head valuations selected only one of these
four classes. This is a fourfold phase enlargement for that specified shape,
not a fourfold enlargement of all certified positive integers.

The minimal Fourier hostile uses only Z/3. Starting at internal j=1, the
allowed increments are {0,2}; starting at j=3 they are {0,1}. These two masks
have identical cosine coefficients and Fourier magnitudes, but +1 is legal
only in the second. Their sine coefficients have opposite signs. The omitted
phase is exactly a direction-sensitive route guard.

There is also a literal sheet mode. For odd target h coprime to 3, signed
predecessors (2^k h-epsilon)/3 require epsilon=s(h)(-1)^k, where s(h)=+1 or -1
is its nonzero residue modulo3. On Z/(2*3^p), this is the Nyquist coordinate;
restricting to k=k_0+2j leaves the odd-modulus guard j modulo3^p. Thus sheet and
guard really can occupy the even and odd parts of one character system.

## 4. Zeckendorf supplies an exact guard reader, with a needed second register

[The saved automaton](zeckendorf_guard_automaton_20261003.md) reads canonical
Zeckendorf digits high to low. Retain X (value) and Y (value with every
Fibonacci weight shifted up once). Appending digit d has the exact update

    (X,Y) -> (Y+d, X+Y+2d).

The preceding digit rejects 11. Modulo a fixed guard modulus this is a finite
automaton processing unbounded words. It evaluates the same arithmetic guard
as the ordinary integer representation; it does not prove that every
integer's Collatz route terminates. For the earlier charge convention,

    b(n)=floor((n+1)/phi^2),  Y=2n-b(n),
    Q(n)=(2Y-3X,2X-Y),       c(n)=(X mod2,Y mod2).

The three nonzero charges cycle under the Fibonacci matrix modulo2. Their
neutral fourth state is necessary. On a fixed-sum diagonal the carry
delta=b(a)+b(b)-b(a+b) in {-1,0,1} is the missing correction to naive XOR;
the integer carry, rather than its parity alone, is needed to reconstruct
the full charge. The three nontrivial Walsh characters of F_2^2 are a useful
Fourier basis for this four-state charge. Their multiplication is XOR in the
dual group, not addition in Z/3; the cyclic permutation of colours is an
automorphism, not the charge's addition law.

For an odd guard modulus M, the pair of registers modulo 2M retains both the
guard residues and this charge by CRT. Specify whether the input integer is
the route source, internal sibling index j, or an exponent; the same reader
works for each, but their residues are not interchangeable.

Hostile: n=2 and n=8 have the same n mod3, charge (0,1), and final digit0.
Their shifted values are 3 and13, unequal modulo3. These three reduced
coordinates cannot replace Y. The saved finite audit checks 70,007 register
states and all 3,174 existing macro positions through this independent path.

The user explicitly confirmed that the intended construction is the
red/black/blue diagonal pattern with extra copies of1. Its historical
35-symbol R/K/B word differs from the Wythoff/lowest-digit candidate
at 35, and the charge convention cannot match even positions2 and4 under
relabeling. The extra-copies-of-1 interpretation remains in its original
scope; a unique global continuation is still unresolved. No new result
silently identifies any of these colour schemes.

The faithful interface is a **marked representation**, not one colour of
the visible integer. Under the inherited ordered alphabet
1_K,1_B,1_R,2,3,5,..., every positive n has two or three nonadjacent-support
representations: Z(n), {K}+Z(n-1), and, when the ordinary red unit is absent
from Z(n-1), {B}+Z(n-1). The marker and its atom position distinguish the
equal-valued units. Normalization can forget that distinction, so a route
record must retain the original marked word when its history matters.

There is now a native reader for these representations. Read the ordinary
part first into (X,Y), retaining r, its red-unit digit. An added black or
blue unit updates the canonical numeric registers by

    (X,Y) -> (X+1,Y+2-r).

Blue requires r=0 by nonadjacency; black permits either r. Keep the marker
as a separate field. This follows from b(X+1)-b(X)=r and Y=2X-b(X).
The resulting pair (value,marker) recovers the complete extended support;
the reader can therefore evaluate modular guards without losing unit type.
The exact audit enumerates all524 representations of n=1..200 and performs
1,572 native-reader checks. It also checks that precisely two representations
in every fibre preserve the full inherited integer charge.

For an actual coloured row retain f(d,c)=1 when diagonal d has colour c.
On the colour coordinate c in {0,1,2}, the three-point Fourier transform
gives a constant channel and two conjugate channels. Equivalently encode
K,B,R as 1,omega,omega^2, where omega=exp(2 pi i/3), *at each position*.
Tensor these with the spatial modular characters to obtain a faithful
Fourier view of the given finite row. Keeping all coefficients retains both
colour and diagonal order; totals or power spectra do not. This C3 colour
transform is distinct from the four-state Walsh transform above.

The colour labels may be permuted as a change of coordinates, but that
permutation is not automatically a symmetry of the atom grammar: K+R is
allowed and K+B is forbidden by nonadjacency. Transport the ordered atom
positions too, or declare the lost structure. Likewise the observed blue
boundary at35 is a legal marked row; the inherited descending-piece grammar
forces red there at row36. Store the row and boundary marker instead of
postulating a fixed infinite diagonal stripe that changes the supplied data.

## 5. Another Fibonacci connection encodes time, not the source's digits

For ordinary Collatz T(n)=3n+1 on odd n and n/2 on even n, consecutive parity
bits cannot both be1. An odd accelerated step with valuation a contributes
the ordinary parity block 1 followed by a zeros. This is the golden-mean
word constraint already recovered in THM-4528's audited result note. The
parity itinerary of n and the Zeckendorf digits of n are different words.

For a *certified* first-hit route of ordinary length L to1, let B(z) be its
finite parity polynomial. Appending the known 100-periodic tail gives

    F_n(z)=B(z)+z^L/(1-z^3),
    P_n(z)=(1-z^3)B(z)+z^L.

This is a finite storage format with its home tail built in. It can be read
in binary, in the golden coordinate, or at Fourier characters. For an mth
root of unity, folding P's coefficients modulo m computes its character
value exactly in the group algebra. A finite set of such evaluations loses
coefficients once their exponents alias; retain P or a compressed route DAG,
the length and root pointer. A spectral fingerprint is not a route proof.

## 6. A constructive gain from a different old cycle

[The -17 return macro](collatz_minus17_return_20261003.md) adapts the retained
source-budget proof of the -5 family. With P=2187,Q=2048 and
W=(1,1,1,2,1,1,4), choose t_m>=1 minimally so

    2^t_m (Q^m-17) > P^m-17,
    beta_m P^m = 17 mod2^t_m.

Every positive b=beta_m mod2^t_m gives n=bQ^m-17 whose exact first descent
occurs after 7m odd steps, to oddpart(bP^m-17). Each m>=2 adds a cylinder
outside the named old 65-cylinder bank and the entire inherited -5 family.
The first is n=4194287 mod8388608, with least representative descending at
step14 to597869. This proves new first descent relative to the specified
bank; it is not by itself a stored route all the way to1 for every parameter.
Compose with a certified lower-source tail when available.

The ledger 7v2(n+17)+11v3(n+17) is preserved by repeated blocks. Its eleven
binary phases and seven odd steps show why every observed triple cannot be
one universal three-state controller. The old -5 sources are3 mod8, these
sources7 mod8: modulo4 merges them, a primitive modulo8 Fourier character
separates them. The size budget still does the descent work.

## 7. Connection audit and next concrete step

| source -> target | map / preserved predicate | loss and needed sidecar | cheapest test |
|---|---|---|---|
| multiplication wedge -> Fourier kernel | residue -> character; full product preserved | cosine forgets reflection sign; keep sine | conjugate guard masks on Z/3 |
| cyclic residues -> tournament block | sign sine of difference; all odd pairs oriented | labels/gauge lost; size suffices only for chosen code | even antipodal tie |
| compressed prefix -> phase mask | all exponent-parity lifts; arithmetic legality | root padding/size not encoded; keep root and parameter | (h,j)=(1,0) |
| Zeckendorf digits -> modular guard | two-register transducer; exact residues | finite colour alone loses shifts; keep Y and digit word | n=2 versus8 |
| route -> parity polynomial | first-hit word plus root tail | finite spectra alias time; keep length, root, polynomial/DAG | coefficient fold and unfold |
| negative cycle -> positive descent family | affine translation plus tail budget | descent does not identify final home; keep lower-source pointer | one-bit budget weakening |

The recommended next implementation is a **typed certificate compiler**:
store a construction rule, exact congruence guard, ordered affine summary,
source range, and certified suffix pointer. Compile the same guard to an
integer residue mask, a phase-preserving Fourier view, and a Zeckendorf digit
automaton. Compare candidate extensions by uncovered source cylinders and
verified first-descent targets. The guard and automaton compilers have now
been implemented; their integration into a growing certificate bank remains
OPEN. Prioritize actual new coverage over another colour-correlation census.

The Fibonacci arithmetic literature provides an algorithmic precedent:
[Ahlbach, Usatine and Pippenger, Efficient Algorithms for Zeckendorf Arithmetic](https://arxiv.org/abs/1207.4497)
proves linear-time addition/subtraction. We do not import a Collatz consequence
from it; the small guard transducer above is independently derived.

## 8. Reproduction and scope

Run the four scripts in `04-computation/experiments/`:

    python modular_sine_tournaments_20261003.py
    python collatz_fourier_guards_20261003.py
    python zeckendorf_guard_automaton_20261003.py
    python collatz_minus17_return_20261003.py

Run from that directory or prepend its path when running from the repo root.
Outputs are saved under the matching names in this results directory.
The modular script uses exact integers: m=1..101, all348,551 table cells,
5,151 multiplication-image counts, 50 odd centre blocks, 2,550 cyclic
differences, 50 even tie controls, 25 odd primes, and 5,000 group-algebra
checks from all certified seeds1..1000. No approximate eigenvalue or
trigonometric calculation is used to establish these checks.

Independent lanes audited the guard phase count, its first-hit exception,
the cyclic sine eigenvalue sign, and the -17 return construction. Detailed
finite universes and hostile controls live with the three specialised notes.

All four saved outputs match fresh replay after LF normalization. SHA-256
hashes below use that same LF-normalized byte convention:

| script/output stem | script SHA-256 | output SHA-256 |
|---|---|---|
| modular_sine_tournaments_20261003 | 7ebef083db6034313f76f93d9f024f55e5dc4f77721dd1a122ff7d50d2ea7e31 | bf7ca77248667031bbd010233fd69f92a6c7f55ed20aa8ee52ccafe8117f152a |
| collatz_fourier_guards_20261003 | 10b267be40d7a43de77860ffde5c315ccbb026029f73793b18c7ff8828715648 | ccb7eb1c21bc27c0dbc4840a959b8e49cf254ae9ee835beb7f3d1f2272197013 |
| zeckendorf_guard_automaton_20261003 | 8cfb6bda1ffe21cc75b0c2ba61db38928039d2e114b03dff41f10c731aff1099 | 98f72ef3d2b5e63fa2920dd9df5bbec59372367fb4667ead23ba6efd67c8cf90 |
| collatz_minus17_return_20261003 | 56661cff3610f7497d8e9f8f65deb4c92d369707b8b785cbe1f0654c6349f33d | ec52a1c20c20f6dc9a46c396a1b4bbcd6aa062f99a0fa3dd0d2e508ae9f9c557 |

Repository documentation audit: the startup surface remains within its byte
budget. The unchanged hypotheses index still has125 maintained lines against
its120-line budget, the pre-existing failure also present at session entry.

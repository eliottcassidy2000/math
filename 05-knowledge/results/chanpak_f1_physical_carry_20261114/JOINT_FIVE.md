# Two-coordinate physical carry, including nonunit reader targets

For every k≥1 and every pair (a,b) modulo 5^k there is an explicitly
chosen strong tournament of one order, \(12k+5\lfloor k/2\rfloor+11\),
with ordinary Hamiltonian count congruent to a and unsigned root
response R_+ congruent to b. The latter is an actual rooted response,
not a second freely assigned integer in a formal transfer.
The construction does not put prescribed **absolute** values in both
coordinates. Its residual tests distinguish that stronger requirement.

A background pass first found an alphabet giving the same statement
for unit b, with 14-vertex digit blocks. Two further changes matter:
use different finite supports in the first and later roles, and use
a reader whose scale is deliberately divisible by five when b is not
a unit. This reaches every pair of residues without losing control
of the ordinary count.

## A positional source with two separate role alphabets

A digit is a literal directed contact mask on the four-vertex
count-five core and two transitive wings, in the notation of
[POSITIONAL_FIVE.md](POSITIONAL_FIVE.md). All digits have
(h,S,T)=(5,0,0). Use the following exposed start menu, listed by
residue in increasing order:

| role | size | f | exposed start labels P |
|---|---:|---:|---|
| first, unit labels only | 11 | 4096 | –, 3806, 3222, 3488, 2734 |
| later, all labels | 12 | 8192 | 6430, 4746, 6932, 5848, 5534 |

For every first label, five physical masks have distinct augmented K
modulo five. For every later label, five masks have distinct converse
end response Q modulo five. The companion script constructs these
45 choices deterministically from the fixed-width masks and verifies
each by literal subset-path dynamic programming. No integer stand-in
for a physical digit is used. The two menus need not coincide:
only K matters when decoding the first place, while only Q matters
in later places.

Start with a fixed four-vertex transitive contact word 1010. Then
P_0=10 and H_0=9. With F_1=4096 and F_j=8192 for j≥2, the ordinal
gluing recurrence for k digits is
\[
 P_j=F_jP_{j-1}+5^{j-1}P'_{j},\qquad
 H_j=5H_{j-1}+5^{j-1}(K_j-5)+P_{j-1}Q_j .
\]
The label \(P'_j\) comes from the indicated role's menu and is the
only choice visible to the first coordinate's recurrence for P.
This is an instance of the [exact common-apex product](../chanpak_movable_reader_returns_20261114/INDEPENDENT_RETURNS.md).

**Separate the two choices.** Prescribe any unit \(w\bmod 5^k\)
for P_k. Subtract \(10\prod_j F_j\) from w, and identify the first
label modulo five using the known multiplication by later
\(F_j\)'s. Subtract the exact weighted label and divide by five.
Continue forward through the roles. These steps give unique actual
P labels, because the F_j are units and the respective menus have
one label per required residue. The first label is a unit, hence
so is each P_j for j≥1.

Now prescribe H_k modulo 5^k. Peel the masks from right to left
inside those already chosen P classes. For j≥2 the invertible
coefficient P_{j−1} makes the five Q residues distinguish the
masks. At j=1, P_0≡0 and the five K residues distinguish them.
Subtract the exact increment, divide by five and continue. For
prescribed absolute (P_k,H_k), check the recomputed P_k and the
final quotient H_0=9, rather than accepting a residue match.

The two role widths are optimal **for these stationary role
requirements** at h=5, S=T=0. The unique-core fact means that
every width-D puncture arises by splitting D padding vertices
around that core and choosing the three contact masks. A complete
census has no full K or Q class for D≤6. At D=7 the P residues
having a complete K class are {1,2,3,4}, and those having a complete
Q class are {1,2,3}; D=8 has all five in both. In particular
it would be wasteful to demand complete K and Q in each mask
class at every position.

## A physical reader for arbitrary valuation

To finish, append an ordinal right reader D whose every vertex
beats the apex. It consists of two certified optional kinds of
tournament blocks and a transitive closing tail:

| block order | outgoing masks | h | f | 5-adic valuation of f |
|---:|---|---:|---:|---:|
| 5 | (24,1,3,6,14) | 11 | 100 | 2 |
| 6 | (10,44,57,48,35,1) | 29 | 320 | 1 |

For 0≤t≤k use \(\lfloor t/2\rfloor\) copies of the first and,
if t is odd, one of the second. Add singleton vertices at the
end to reach fixed reader order \(5\lfloor k/2\rfloor+7\).
There is always at least one singleton. The ordinal products give
\(v_5(f_D)=t\) and \(5\nmid h_D\).

Since all of D beats the apex, its own augmented count is h_D,
its positive out-contact response is zero, and its dual response
is f_D−h_D: the latter is the sum over all nonempty selected
subsets. The main site's signed right endpoint is zero, so the
final full count and rooted response are exactly
\[
             H=h_D H_k+(f_D-h_D)P_k,\qquad
             R_+=5^k h_D+f_D P_k .
\]
The second equality separates the singleton-apex path from the
weighted subset paths whose first old vertex is a contact.
Also \(R_-=-5^k h_D\), because the closing singleton forces every
full path in the puncture to be extendible by the apex. The first
left vertex is beaten by the apex, and the last right vertex
beats it, so the whole graph is strong. Its order is constant
in t for each k.

Given desired b modulo 5^k, let t be its 5-adic valuation, capped
at k. If t<k, solve
\[
            P_k \equiv (b/5^t)(f_D/5^t)^{-1}\pmod {5^{k-t}}.
\]
Any lift to mod 5^k is a unit, so the source decoder supplies
the P labels. If t=k, take any unit P_k; R_+ is already zero
modulo 5^k. Use the actual P_k and solve
\[
                 H_k\equiv h_D^{-1}(a-(f_D-h_D)P_k)\pmod {5^k}
\]
with the backward mask decoder. This covers all \(25^k\) residue
pairs. Simply allowing a nonunit P_k in the original bank would
not work: all P_j for j≥1 would be nonunits, so its last Q digit
could not select an arbitrary ordinary residue. The extra reader
handles the valuation **after** the source has done its selection.

For example, at k=1 the prescribed residues (1,0) give the
23-vertex outputs (28,121,831, 28,650,385) by the exact transfer;
the first is nonzero modulo five while the second is zero.

For absolute positive (H,R_+) at fixed k, infer t from the
capped valuation, solve the displayed R_+ equation by exact
division for P_k, recover the role labels and test the recomputed
integer, then solve and test H_k. This is an exact membership
decoder for the specified bank (with distinguished masks).
No choice of the residue lift alone proves absolute attainment.
The same residue control transfers to an actual positive bridge:
its deletion count is H and its second coordinate is a fixed
dyadic multiple of R_+, a unit rescaling modulo 5^k.

Even the enlarged fixed bank is sparse in actual sizes. Its count
per level is at most \((k+1)25^k\), and the positive cross term
grows at least as a constant times \(8192^k2^{5\lfloor k/2\rfloor}\)
because \(f_D\ge2^{|D|}\) and \(f_D-h_D\ge f_D/2\). Thus this
all-residue result cannot, by itself, deliver full absolute
coverage. Run **python3 replay_joint.py** for the literal digits,
minimum role-width census, exact inverses and residue checks.

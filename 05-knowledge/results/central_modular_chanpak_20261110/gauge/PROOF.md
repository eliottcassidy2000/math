# The affine Collatz gauge and the tournament boundary

**PROVED, elementary and scoped.** This is a small obstruction to a specified
scalar simulation, plus a repair of its intended role. The Chan–Pak
classification and positive-integer Collatz convergence remain **OPEN**.

## G1. The three natural gauges are forced

Let C(x) be ordinary Collatz: 3x+1 for odd positive x, x/2 for even x.
Suppose an affine label L(x)=ax+b, with a>0, is required to obey both

    L(C(x))=(L(x)+1)/2       for every even x,
    L(C(x))=3L(x)+delta      for every odd x.                (1)

Then necessarily b=1 and a=delta+2. Conversely those values satisfy (1).
Indeed the even identity gives b=(b+1)/2; the odd identity then gives
a+1=3+delta. This is an identity on each whole parity class, not a fit to
sample trajectories.

For delta=0,2,4 the unique gauges are therefore

| Odd branch | Gauge | Modulus or label domain | Allowed label reaching a gap |
|---|---|---|---|
| 3L | 2x+1 | odd L>=3 | 13 -> 7 and 41 -> 21 |
| 3L+2 | 4x+1 | L=1 mod4, L>=5 | 41 -> 21 |
| 3L+4 | 6x+1 | L=1 mod6, L>=7 | 13 -> 7 |

The second row is exactly the plus-sign moving-modulus conjugacy in
[central C2](../central/PROOF.md). This explains why a ternary return
appears naturally in the user's tables. The two ternary increments require
different affine gauges; they are not two selectable branches of this one
Collatz map on a common scalar label.

## G2. Exact scalar simulation runs into the known gaps

**THM-1370**, [h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
establishes that 7 and 21 never count Hamiltonian paths in a tournament.
It also gives actual realizations of 13 and 41. Consequently no assignment
of an actual tournament T_x to every positive x can satisfy
H(T_x)=L(x) for any of the three gauges above. For the delta=0 gauge,
L(3)=7 and L(10)=21; for delta=2, L(5)=21; for delta=4, L(1)=7.

Deleting those source labels does not make the proposed simulation closed:
the last column of the table gives an attained label whose even update is
a forbidden label. This proves an obstruction to *these gauges with both
branch laws on their stated domains*. It does not exclude a restricted
domain with additional guards, an address-based encoding, more coordinates,
or a different tournament operation. In particular it is not an argument
against the Chan–Pak conjecture.

The first failed implication would be “an exact arithmetic graph conjugacy
is an exact Hamilton-count realization.” The missing coordinate is the
realization map, including its domain and the physical source preserved by
its operations. Reversing the affine equations cannot supply that map.

## G3. The useful repair: arithmetic on a constructor address

The [fixed-core carry lift C3–C4](../compiler/PROOF.md) supplies a lawful
place for the arithmetic. Given an actual named core and cut, an integer
address B selects a binary flank and therefore an actual tournament with
count Phi(A,B). The odd Collatz update followed by one halving is

    x=2n+1:  (3x+1)/2=3n+2,   (3x-1)/2=3n+1.

Both destination channels are explicit coordinates in the multiplier-three
carry lift. This constructs actual tournaments at the new addresses while
retaining their exact counts. It imposes no false identity H=4x+1 or H=6x+1.

| Connection | Preserved | Lost or still required | Cheapest decisive test |
|---|---|---|---|
| Collatz integers -> affine labels | Both guarded arithmetic laws | Tournament realization | The labels7/21 and boundary13/41 |
| Address -> fixed-core graph | Exact incidences and Hamilton count | Arbitrary old flank names and joint ports | Literal recount after an address change |
| Golden word -> normalized state | Three-coordinate phase at every modulus | Actual normalized-state realization and return readiness | Check M=rU, det U=1, then audit a named return |

The last row uses [Gcontent](../information/PROOF.md): phase can be retained
while multiplying content by5,9,21. The integer inverse of U is a decoder,
not automatically an executable inverse graph operation. To complete an
induction one still needs a positive actual chart switch, with named ports,
whose exact count is3H+2 or3H+4 and whose output lies in a usable next chart.
The carry lift and singular words now give explicit objects on which to
state and test that obligation.

## Reproduction

`python3 -B check.py` checks both branch laws for all x=1..4096 and all
three gauges, the displayed forbidden preimages and boundary arrows, and
the uniqueness equations for integer a=1..20,b=-5..5. The proof above
supplies the all-input quantifiers. The all-order tournament exclusion is
an inherited theorem, not inferred from this finite arithmetic check.

# Literal norm readers with target-selected profiles

The [crossed seam](POLES_TO_PRIMES.md) already gives one positive binary quadratic form for two identical blocks. The phenomenon is not isolated: a directed quotient can read a definite norm or a positive branch of an indefinite norm from the same two integers.

For a nonempty block U put h=H(U), b=f_U(2), u=b−h. Replace marked vertices 0 and 1 in each quotient below by **separate copies** of U. All other vertices stay singletons. Compress consecutive visits to either block into one run. Each quotient has, for a spanning path, only visit shapes (1,1), (1,2), (2,1), (2,2). Their counts are shown in this order. The independent replay enumerates all shapes up to the absolute cutoff e+1, where e is the number of ordinary vertices: the two marked quotient vertices alone are transitive.

| Outgoing masks | Run coefficients | Exact H |
|---|---:|---|
| (2,4,1) | (3,1,1,1) | u²+2h² |
| (6,4,8,3) | (5,2,2,1) | (u+h)²+h² |
| (30,8,58,48,34,3) | (23,11,15,7) | (2h+3u)²−2u² |
| (2,4,9,19,7) | (9,5,3,3) | 4h²+2hu+3u² |
| (54,16,10,19,36,14) | (29,15,15,7) | 6h²+16hu+7u² |

Insert h for one run and b−2h=u−h for two runs to verify each polynomial. For the Pell row literal transitive U of orders 1,2,3 gives H=23,103,431, all prime. The Pell form is positive on every physical input because u≥h>0; it is not a subtraction of signed graph counts. The Gaussian row similarly gives prime H=17 from TT2, and the first gives prime H=11. The crossed diagonal supplies a fourth shape with discriminant −15 in [Prime-safe phase](PRIME_SAFE_PHASE.md).

This changes what it means to choose a tournament for a prime. Instead of first choosing a center and correcting it, solve a norm equation for the requested number and recover the two block readings. In the indefinite example seek x²−2u²=p with x−3u a positive even integer; then h=(x−3u)/2 and b=h+u. The required region x>3u, u≥h bounds u by sqrt(p/7), so this is a finite direct factor-profile test. For p=431 it selects x=23,u=7,h=1,b=8, and TT3 is the physical block. The definite analogues give other filtering regions, not merely other names for the same family. The main open job is to realize a profile on one of these target-selected norm branches, or invent more quotients to read physically easier branches. Arithmetic representability by an unrestricted norm does not suffice.

### A bank with no common quadratic splitting obstruction

A necessary arithmetic screen for a primitive binary quadratic form of discriminant D to represent a prime p not dividing D is that D be a square modulo p. It is not sufficient for a specific integral form or for a physical profile. The bank above passes this **screen** for every prime p>5, in a branch with the correct mod-8 output:

- If p is 1, 3 or 7 modulo 8, use respectively the Gaussian, −8 or Pell row. Their square classes at these primes are −1, −2 and +2.
- If p is 5 modulo 8, the other two rows have discriminants −44 and +88, with square classes −11 and 22. Their quotient is −2 up to a square, a nonsquare in this class, so precisely one passes the square test. (Neither discriminant vanishes since 11 is 3 mod 8.)

For the last two rows, b≡0 mod 8 gives u≡−h and output 5 mod 8. These are not just formal choices: the −44 row gives literal primes 37,709 using TT2,TT4 and the +88 row gives 29,461,7229 using TT1,TT3,TT5. The q=2 reflection congruence forces opposite parity when these two rows are used to produce a prime above 3: even block order for −44, odd for +88. Passing the square test does not prove representation by the chosen form, let alone under the joint block constraints. It removes one reason to insist on a single discriminant in an inverse search.

No identification with the Leech lattice is asserted. The operative analogy is that a fixed finite readout constrains which norm coordinates are available; selecting a congruence class of norms is weaker than finding a vector in the permitted chart. The three exact quotient tables make that distinction usable on tournaments.

Replay: `python3 norm_readers.py` alongside the earlier scripts.

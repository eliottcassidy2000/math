# No complete finite synchronous digit graph for a fixed low-cycle prefix

Session 197. **An all-order obstruction to one precise automaton strategy.** The one-gadget count function for any fixed physical prefix has no synchronous regular binary graph, even after discarding any finite initial interval. This does not rule out a regular or cofinite target range, an automatic selected-witness subrelation, or a target-dependent inverse algorithm.

## 1. Arithmetic and encoding conventions

The reviewed one-seam grammar has a physically reachable rooted prefix (h,e,R), with h>=1, R>=h and 0<=e<=h. Its final gadget is indexed by every integer L>=1, and its actual count is

    F_(h,e,R)(L)=h g(L)+R L−e(L mod2),                 (1)

where

    g(1)=1,
    g(2k)=g(k),
    g(2k−1)=g(k)+k−(k mod2), k>=2.                  (2)

Every such integer gadget is physical. The theorem below uses the arithmetic recurrence, not a relaxation of physical reachability. The inherited exact source is `session197-lowcycle/ARITHMETIC_INVERSE.md`, SHA256 `4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7`.

A synchronous binary graph means a regular language over pairs of bits. Each number is written most-significant-bit first, without leading zeros, then the shorter track is padded on the left by zeros to the longer length. Allowing arbitrary additional common leading-zero columns is equivalent for regularity. Reading both tracks in reverse order is also equivalent, since regular languages are closed under reversal.

**Theorem.** The padded graph of g is not regular. More generally, for every fixed integer h>0, R>0 and e>=0 for which (1) is nonnegative on its domain, the padded graph of F_(h,e,R) is not regular. This includes every physical prefix. The same conclusions hold for an eventual graph restricted to L beyond any fixed finite cutoff.

## 2. A balanced-run slice of the gadget function

Put f(n)=g(n+1), n>=0. Equation(2) gives

    f(0)=1,
    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod2).                           (3)

The last equation also holds at x=0. Hence f(2^a−1)=1. Appending b zero bits to the binary word 1^a adds

    (2^a−1)(1+2+...+2^(b−1))+1,

where the final extra1 comes only from the first, odd prefix. Thus for a,b>=1,

    f((2^a−1)2^b)=(2^a−1)(2^b−1)+2
                   =2^(a+b)−2^a−2^b+3.             (4)

We use only a,b>=3 below, avoiding all short-word boundary ambiguities. Both input and output then have exactly a+b binary digits.

If a=b=k, equation(4) has binary expansion

    1^(k−1) 0^(k−1) 11.                             (5)

If a and b are unequal, set s=min(a,b), t=max(a,b). Its expansion is

    1^(s−1) 0 1^(t−s) 0^(s−2) 11.                  (6)

All displayed internal runs have positive lengths in (6), because t>s and s>=3. These identities follow by summing their positional binary values, so they are exact at arbitrary length.

Consequently, on inputs in the regular language

    1^{at least3} 0^{at least3},

the output belongs to the regular language

    1^{at least2} 0^{at least2} 11

if and only if a=b. No padding columns occur on this slice, since both lengths agree.

## 3. Nonregularity and the fixed-prefix transfer

Suppose the padded graph of f were regular. Intersect it with the two regular track conditions just displayed, and project onto the input track. Closure under intersection and projection would make

    {1^k 0^k : k>=3}

regular. It is not: the prefixes 1^i and 1^j, for distinct i,j>=3, are distinguished by the suffix 0^i, which completes exactly the first one to a word of this language. Infinitely many distinguishable prefixes contradict finite-state recognition. Therefore the graph of f is not regular.

The binary addition relation L=n+1 is synchronous regular. One may recognize it by a finite carry while reading least-significant-bit first, then reverse the language. Regular relations remain regular under product, intersection and existential projection. After projection, take the regular left quotient by the language of common-zero columns, then intersect with canonical pair encodings. This trims any forced common zeros inherited from a longer auxiliary track. Thus a regular graph of g would give a regular graph of f by composition with addition by1, a contradiction.

Now fix h,R,e as in the theorem. The relation

    y+e epsilon=h z+R L,    epsilon=L mod2,           (7)

is synchronous regular: multiplication by a fixed nonnegative integer, addition, and the parity test all use finite carries. This assertion concerns fixed coefficients; it does not introduce an unbounded numerical register. If the graph y=F_(h,e,R)(L) were regular, intersecting it with (7) and projecting away y and epsilon would give precisely z=g(L). Positivity of h makes this recovered z unique. This contradicts the preceding result.

All products and projections may be performed after allowing common leading-zero padding on every track. Projection can leave common leading-zero columns that were forced by a longer auxiliary track; canonical intersection alone could discard the resulting pair. First take the left quotient by the regular language (0,0)*, which deletes any number of common leading-zero columns, and only then intersect with canonical pair encodings. Regular languages are closed under left quotient by a regular language. This yields exactly the canonical recovered graph, so the argument is independent of harmless conventions about common zero columns.

Finally, finite changes preserve regularity. If a graph were regular only beyond a fixed cutoff, adding its finitely many omitted input/output pairs would give a regular full graph. Thus no finite initial exception repairs the obstruction.

## 4. What is and is not excluded

This rules out a finite synchronous automaton recognizing **all exact pairs (L,F(L)) for a fixed prefix**. A finite bank with an explicit prefix tag cannot have a regular complete tagged graph either: fixing one tag would recover a forbidden component graph.

It does not rule out any of the following:

- A regular target range, or a cofinite union of target ranges from finitely many physical prefixes.
- An automaton recognizing a suitable proper subfamily of gadget indices or selecting only some witnesses while still covering desired targets.
- An untagged union of several relations, without a way to recover the individual prefix.
- A finite-dimensional linear representation with unbounded integer registers, a pushdown-like method, nonlinear arithmetic, or a target-dependent search.
- The reviewed congruence inverse or the separate square-root fixed-prefix inverse, whose low-word split grows with the target.

Projection of a nonregular relation can be regular. Thus this result redirects the finite-carry approach; it does not settle the all-order coverage question and does not turn a failure of this strategy into a tournament-spectrum omission.

## 5. Reproduction and relation to the low-word inverse

The verifier is standalone. It checks (3), the arbitrary-length formula and both exact binary forms on a finite control grid, the regular-filter equivalence, and recovery through (7) for several fixed prefixes. Its finite checks corroborate the algebra; the balanced-language proof establishes nonregularity at all lengths.

The separately frozen `session197-lowcycle/SQRT_FIXED_PREFIX_INVERSE.md` derives a compatible low-word identity involving complement-bit reversal and uses a target-dependent split of size O(sqrt(N/R)). That is not a bounded-state synchronous graph and does not conflict with this theorem. No claim about its audit status is needed here.

Run `python3 verify.py` and `python3 -O verify.py`. Neither execution writes files. Their receipts must agree. No graph census, web access, or unproved finite automaton minimization is used.

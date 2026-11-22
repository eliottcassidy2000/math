# A signed root-wrap calculus and a small compiler obstruction

This is not a proof of cofinite Chan–Pak coverage. It isolates where
a physical radix compiler would need more states or larger tapes.
The leading-zero two-coordinate transfer below is the inherited
[rooted migration](../chanpak_bank_fusion_20261117/migration/PROOF.md);
the point here is to analyze signed intermediate stages, including
empty tapes, and a complete short-tape language.

## Signed state and physical operations

Let a tournament \(D\) have a marked root \(a\) that starts exactly
one full path. Write \(h=H(D)\), let \(c\) count ordered two-slot
covers with the first slot a nonempty path starting at \(a\), and
let \(e\) count full paths ending at \(a\). A cover can have an empty
second slot and \(c\ge h\): split every full path just before the
root, putting its suffix in the first slot.

Append an ordinal transitive tail of length \(\ell\), then a fresh
root \(b\) beating the old root and the selected tail tags, and
losing to every other vertex of \(D\). The left-to-right tail tags
are the reverse complement of the \(\ell\)-bit padded address
\(y\), with \(0\le y<2^\ell\). A full path starting at \(b\) must
visit all of \(D\) before the tail, and hence its next vertex
must be the old root. It is unique. If \(D\) was strong, the new
host with nonempty tail is strong exactly when the last tail
vertex beats \(b\), equivalently the leading padded bit of \(y\) is one. With an empty
tail it is strong when \(D\) was nontrivial and strong.
A singleton child gives a nonstrong two-vertex host with empty
tail, but every nonempty tail with leading address bit one is strong.

Let \(P\) denote the binary read \(P(0)=0\),
\(P(2v+1)=P(v)\), \(P(2v)=P(v)+v+[v\ {\rm odd}]\).
For a nonempty tail the exact signed migration is
\[
 \begin{aligned}
 h'&=(1+P(y))h+(y+1)c-[y\ {\rm even}]e,\\
 c'&=(2^\ell-\operatorname{rev}_{\ell}(y))h+2^\ell c,\\
 e'&=[y\ge 2^{\ell-1}]h .
 \end{aligned}
\]
For an empty tail it is instead
\[
                (h',c',e')=(h+c-e,\ h+c,\ h-e).
\]
To derive the first coordinate, attach an apex beating only the
marked root of \(D\): its count is \(h+c-e\) by the singleton-tag
even-path identity (the signed tagged-start mass is \(-e\)).
Apply the exact right-flank reader \(cy+hP(y)+e[y\ {\rm odd}]\).
For the second coordinate, the singleton new-root path contributes
\(h\). Any further first-slot vertex is an admissible tagged start
in the ordinal deletion, whose weighted mass is \(2^\ell c\) plus
\(h\) times the numeric tail tag word. This word has value
\(2^\ell-1-\operatorname{rev}_\ell(y)\). To end at \(b\) after a
nonempty tail, the last tail vertex must beat \(b\), and then there
are \(h\) choices. With empty tail a path ending at \(b\) may use
any old full path except the \(e\) ending at the old root.

These recurrences are physical even when some tapes have a tagged
last vertex, but the following obstruction statement requires
every intermediate host to be strong.

## Complete one-tail compiler has density zero

The only moves with tail length at most one are the empty-tail
move \(Z\) and the untagged singleton-tail move \(A\):
\[
 Z(h,c,e)=(h+c-e,h+c,h-e),\qquad
 A(h,c,e)=(h+2c,h+2c,h).
\]
Starting from the singleton, the first strong move must be \(A\),
giving \((3,3,1)\). All later states obey
\[
 c\ge h,\qquad e\le(h+c)/3,\qquad h/c\ge2/3 .
\]
Indeed \(A\) has equal first two coordinates and
\(3h\le2(h+2c)\). For \(Z\), \(c'-h'=e\ge0\),
\(h'/c'=1-e/(h+c)\ge2/3\), and \(3(h-e)\le2(h+c)-e\)
because \(h-2e\le2c\). Thus \(c\) grows by at least \(5/3\)
under \(Z\), and by at least \(8/3\) under \(A\).
Choose a fixed \(\sigma<1\) sufficiently close to 1 that
\[
 (3/5)^\sigma+(3/8)^\sigma<1 .
\]
Summing the \(-\sigma\) powers of the products of these two
per-letter lower growth factors over all words gives a convergent
geometric sum. Hence the number of words with \(c\le X\) is
\(O(X^\sigma)\), and the same bound holds for absolute height
\(h\le X\), since \(c\le3h/2\).
Therefore this entire physical binary compiler is not cofinite,
even if every word is allowed and collisions are ignored.

## A finite exact decoder for an arbitrary strong rooted target

Each strong nonempty tape has \(y\ge2^{\ell-1}\), and the displayed
formula then gives \(h'\ge (y+1)c\ge2^{\ell-1}h\), since
\((1+P(y))h-e\ge0\). It also gives \(h'\ge3h\): if \(y\) is odd,
use \(h'\ge h+2c\); if \(y\) is even then \(y\ge2\), \(P(y)\ge2\)
and \(h'\ge3h+3c-e\). The last property implies
\(e'\le(h'+c')/3\) on every nonempty strong move.
The empty move preserves this inequality as above, and then
grows \(h\) by at least \(4/3\). Thus a scalar target \(N\)
gives an explicit finite search: traverse rooted tapes after the
singleton, retain only strong gates, truncate tape lengths by
\(1+\lfloor\log_2 N\rfloor\), and prune every state with \(h>N\).
Every nontrivial strong rooted host starting uniquely at its root
decomposes into this grammar by deleting its root and taking the
initial strong component, followed by singleton components.
That standard decomposition is the recursive one used for the
[Fibonacci receivers](../chanpak_e31_physical_20261122/README.md).
Induction on the root deletion proves the search recovers all its
histories. Keeping the parent skeleton gives a literal witness,
not just a numerical state. This is a finite inverse, not a
polynomial-time or coverage claim.

There is also an exact local backward selector when all three
coordinates are given. The empty move gives
\(e=c'-h'\), \(h=e'+e\), \(c=c'-h\). On a nonempty strong
move, \(e'=h\), so \(c\) is obtained by exact division in the
second equation. If \(y\) is even, the first equation uniquely
determines \(e\); otherwise it checks the required height and
allows \(0\le e\le h\) before recursive predecessor tests.
For a scalar target, \(c\le4h\) on every strong state:
on a nonempty tape \(c'\le 2^\ell(h+c)\) and
\(h'\ge 2^{\ell-1}c\); on an empty tape the stronger \(c'\le
3h'/2\) follows from the invariant above. Absolute scalar
membership can therefore also be tested by a finite list of
terminal states \((N,c,e)\) and these exact backward rules.

Neither result addresses arbitrary tournament completions or proves
a cofinite image of the full variable-width rooted grammar.

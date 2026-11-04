# Independent native and reader audits

**Status: FINITE-EXACT PASS**, with the algebraic reader bound below proved
directly. These audits check the two immediate named returns and the narrow
reader inverse. They do not assert repeated closure or count the order21
endpoint Hamilton paths; the separate interface audit does those endpoint
counts.

## Literal native programs

The producer audit_native.py imports no inherited compiler or graph helper. It
reconstructs the actual inherited order17 source from its named suffix and
wrapper, performs every birth from its attachment cut, checks each Hamilton
guard, constructs every off-path complement from unordered pairs, and deletes
only the recorded temporary helper. Every saved intermediate matrix must match.

| Target label | Permanent births | Old arcs changed inside B | Native moves | Temporary births/deletions | Peak order |
|---|---:|---:|---:|---:|---:|
|512423|4|26|52|26|22|
|512425|4|4|8|4|22|

Both endpoint matrices have order21 and retain every original physical name.
Every completed helper gadget changes exactly its requested old arc. All
cross-component arcs, and the restrictions of the named port stars to the
outside of B, agree at those macro endpoints. The intermediate native
complements may change port stars and Hamilton counts. Root39 still has sole
in-port31; prefix50→51 still dominates the entire new component; spare6 survives.
The audit checks the literal component matrices as well as the full endpoints.

The certificate pins a canonical graph payload and the complete native program
payload. It deliberately excludes mutable descriptive metadata. Positive
TT3/native controls, the same-guard involution, and two rejected illegal guards
are included.

## Reader formula and exact inverse

For a binary word, let \(r=\sum_{i:\epsilon_i=0}2^i\), with positions starting
at0. The independent pair-sum expression

\[
H_{\rm word}-1
=\sum_{\substack{i:\epsilon_i=0\\j<i:\epsilon_j=1}}
 \left(2^{i-j-1}+\mathbf1_{i=j+1}\right)
\]

matches literal one-apex graph counts. It also proves the needed bound without
enumeration: the summand over earlier positions at a fixed zero position \(i\)
is the binary value of the earlier bits plus the previous bit, hence is at most
\((2^i-1)+1=2^i\). Therefore \(0\le H_{\rm word}-1\le r\).

For seed fields \(h,e,t,F,p,q\), the reader formula is

\[
h_{\rm out}=h+q(H_{\rm word}-1)+(t-q)r+(q-e)(r\bmod2),
\qquad F_{\rm out}=F+pr.
\]

This follows by the per-bit affine recurrence: relative to \(q\) copies of the
singleton seed, the surplus in \(t\) doubles at each position and contributes
\((t-q)r\); the initial end-count difference contributes only when the first
bit is0. For outer clock \(M=2^s\), write
\(D=N-h-MF,\ K=t-q+Mp,\ d=q-e\). Under \(t\ge q\ge e\ge0\) and \(p>0\),

\[
Kr\le D\le(K+q)r+d,
\qquad
\max\left(1,\left\lceil\frac{D-d}{K+q}\right\rceil\right)
\le r\le\left\lfloor\frac D K\right\rfloor.
\]

A positive \(r\) uniquely specifies a zero-ending word. The empty word is
handled separately. This is the exact arithmetic inverse implemented by the
reader; the bounded audit checks its implementation independently.

The producer audit_reader.py uses fresh forward endpoint DP and subset-cover
sums. It checks93 literal graphs against all eight fields,255 literal one-apex
graphs against the independent pair sum,6141 scalar word controls through
length10, and every target1 through1001 for three seeds. All3003 inverse queries
match brute word menus, including positive and negative targets. The bound on
\(r\) certifies that the ten-letter brute menu contains every possible candidate
at this target cap.

There is a structural boundary beyond the arithmetic. If the seed has at least
two vertices and its retained mark reaches every seed vertex, a nonempty
zero-ending word is strong: all new vertices are reachable through a second
seed vertex, and each new vertex reaches the final zero vertex and then the
mark. A singleton plus the word0 is instead transitive2, with \(h=1,F=2\).
The hostile singleton check is recorded explicitly. Thus a blanket
strong-ending claim needs seed order at least2, or a separate singleton
strongness filter. The actual seeds used for the two targets satisfy this
condition.

Run each producer normally and with python3 -O. Their adjacent normal and
optimized JSON receipts are byte-identical; all checks use explicit exceptions,
not disabled assertions.

## Full depth-three continuation replay

The separate audit_continuation.py extends the fresh native routines to all
14 edges of the recorded binary tree. Each source comes from the previous
edge's actual replayed rows and physical names, indexed by its branch word.
It is never replaced by another graph with the same Hamilton-count label.
The audit checks every birth cut, all944 directed Hamilton guards and native
complements, all472 helper deletions and intermediate hashes, and all final
matrices. The tree has22 permanent births summed over its edges and peak
temporary order25. Every name of each immediate parent survives its children.

Physical ports39,31,6,50,51 keep their identities. The outer root keeps sole
in-port31 and the prefix always begins50→51. Later prefixes can additionally
contain old suffix names, so the audit checks the actual recorded prefix and
its actual strong suffix rather than assuming every endpoint still has TT2.
For example branch42 has prefix50→51→1→2→4. Such changes of membership are
allowed; they are not deletions or relabelings.

The normal and optimized continuation receipts are byte-identical PASS and
pin stable graph and native-program payloads. The audit imports only the
fresh elementary auditor; it neither imports a native compiler nor regenerates
a program. Endpoint Hamilton counts remain the separate interface audit's
responsibility. The result is a verified depth-three tree, not arbitrary-depth
closure.

# The 2-Adic Column Family Structure

**Session:** oracle-2026-05-15
**Depends on:** `the-two-reductions.md`, `everything-is-the-triangle.md`, `project_n7_findings.md` (SC∩SF identity), INV-149 (HYP-217)

---

## The Grid

Every natural number factors uniquely as $n = 2^r \cdot (2k-1)$. This gives a 2D grid:

| | $k=1$ | $k=2$ | $k=3$ | $k=4$ | $\cdots$ |
|---|--------|--------|--------|--------|----------|
| $r=0$ | $1$ | $3$ | $5$ | $7$ | odd numbers |
| $r=1$ | $2$ | $6$ | $10$ | $14$ | first doublings |
| $r=2$ | $4$ | $12$ | $20$ | $28$ | second doublings |
| $r=3$ | $8$ | $24$ | $40$ | $56$ | third doublings |

The **column family** $F_k = \{2^r(2k-1) : r \geq 0\}$ collects all tournament sizes
reachable from the odd base $2k-1$ by repeated doubling. The two grid generators are:

- **Column step** ($n \to n+2$ along top row): moves $(0, k) \to (0, k+1)$
- **Row step** ($n \to 2n$, tournament blowup): moves $(r, k) \to (r+1, k)$

Together they generate all of $\mathbb{N}$ from $n = 1$.

---

## Connection 1: Mode B Is a Column Step

**Mode B** (from `the-two-reductions.md`) takes the staircase $\delta_{n-2}$ to
$\delta_{n-4}$, i.e., a tournament on $n-2$ vertices.

For **odd** $n = 2k-1$: $n - 2 = 2k-3 = 2(k-1)-1$, the top of column $k-1$.

$$\text{Mode B on odd } n = \text{left-column step in the top row.}$$

**Consequence:** All Mode B recursions — Cayley-Dickson descent, meta-graph recursion,
the two-reductions staircase — are the right-to-left traversal of the top row of the
2-adic grid. Every "slow time scale" result in the project is, geometrically, a
column-by-column traversal.

For **even** $n$: Mode B gives $n-2$ which crosses into a different column, breaking
the clean family structure. This is why Mode B theorems are cleanest for odd $n$.

Mode A ($n \to n-1$) crosses columns in a parity-dependent way and has no clean
column-family interpretation — consistent with Mode A being the "fast time scale"
that mixes families.

---

## Connection 2: The Pairs Anomaly — The Seam at $r = 0 \to r = 1$

Define $\mathrm{pairs}(n) = \lfloor n/2 \rfloor$ (maximum matching number of $K_n$).
Under the two generators:

| Operation | $n \to$ | pairs: before $\to$ after |
|-----------|---------|---------------------------|
| Row step ($\times 2$) | $2n$ | $\lfloor n/2 \rfloor \to n$ |
| Column step ($+2$, top row) | $n+2$ | $\lfloor n/2 \rfloor \to \lfloor n/2 \rfloor + 1$ |

Within column $k$:

$$\underbrace{k-1}_{r=0} \;\xrightarrow{??}\; \underbrace{2k-1}_{r=1} \;\xrightarrow{\times 2}\; \underbrace{2(2k-1)}_{r=2} \;\xrightarrow{\times 2}\; \cdots$$

From $r = 1$ onward, pairs **exactly double** each step.
The $r = 0 \to r = 1$ transition gives $k-1 \to 2k-1$ (expected $2(k-1) = 2k-2$):
**always one extra pair**.

**Structural meaning**: the odd tournament $2k-1$ has exactly one unmatched vertex in
any maximum matching. When doubled (node-blowup), this vertex gets a twin, creating
the extra pair. Every other vertex was already paired; the extra 1 is the "birth of
the first twin."

In tiling terms: the staircase for $n = 2k-1$ has an asymmetric apex tile (the tip of
the right triangle), corresponding to the unmatched vertex. After one blowup
($n = 2(2k-1)$), this asymmetry resolves — but marks a structural **seam** at the
$r = 0 \to r = 1$ boundary that distinguishes the odd top-row member from all even
members of the family.

This seam is where SC/SF behavior changes: odd-$n$ behavior (SC requires involution
with fixed point) vs even-$n$ behavior (fixed-point-free involution, $n/2$ orbits).

---

## Connection 3: SC∩SF Count and Column Structure — Coincidence vs. Theorem

**The count identity** $\#(\text{SC} \cap \text{SF})(n) = \#\text{SC}(n-2)$ held for
$n = 4..7$ but **FAILS at $n = 8$**: $\#(\text{SC} \cap \text{SF})(8) = 5 \neq 12 = \#\text{SC}(6)$
(oracle-2026-05-11-S2, exhaustive verification). It is a **small-$n$ coincidence**, not a
general theorem.

In column-family terms: $n = 2k-1$ is at column $k$, $n-2 = 2(k-1)-1$ is at column $k-1$.
The identity would say adjacent top-row columns have related SC∩SF and SC counts — but this
clean ratio breaks down for $k \geq 5$ ($n \geq 9$).

**What the column structure DOES explain:**
- The SF operation is geometrically "column-stripping": it extracts the middle subtournament
  (Mode B), which lies one column to the left. This is a genuine structural connection.
- The **near-regularity** of SC∩SF classes (score deviation $\leq 1$ from regular) DOES
  persist through $n = 8$ — a more robust invariant than the count identity.
- The related ratio $\#\text{SC}(2k)/\#\text{SC}(2k-1) = k/2$ also holds only for $k = 2, 3, 4$
  and fails at $k = 5$ (oracle-2026-05-11-S2).

**Also note**: at $n = 8$ (row 2 of column 2), max $H = 661$ is achieved by SC non-SF,
**breaking the even-$n$ pattern** (which held at $n = 4, 6$ = row 1 of columns 1, 2).
This suggests "max-$H$ by SC∩SF at even $n$" is a **row-1 phenomenon**: holds at
$n = 2, 6, 10, 14, \ldots$ (first doublings of odd bases) but not at row 2 and beyond.

---

## Connection 4: The Ruler Sequence as Canonical Family Ordering

The sequence of 2-adic depths $v_2(1), v_2(2), v_2(3), \ldots = 0, 1, 0, 2, 0, 1, 0, 3, \ldots$
is the **ruler sequence**. It gives a canonical ordering of tournament family depths:
most time on row 0 (odd $n$), with self-similar fractal dips into deeper family rows.

The ruler sequence is already implicit in the Walsh degree formula (OPEN-Q-035, THM-259):
$$\text{Walsh degree} = 2\lfloor(n-1)/2\rfloor = \begin{cases} n-1 & n \text{ odd} \\ n-2 & n \text{ even} \end{cases}$$

This is exactly the row-0 vs row-$r \geq 1$ distinction in the family grid:
the Walsh degree **jumps by 1** exactly at the $r = 0 \to r = 1$ seam (odd to first
even in each column). The Walsh degree anomaly IS the pairs anomaly in disguise —
both are consequences of the structural seam at the first blowup.

S453/THM-371 sharpen the seam by comparing it with the LRC denominator grid.
For odd $m$, reduction $U(2m) \to U(m)$ is a bijection: every unit residue mod
$m$ has exactly one odd lift among $a$ and $a+m$.  Thus the first doubling does
not create a second unit sheet, even though the tournament matching count gains
the missing twin.  The critical local picture is therefore asymmetric:

| structure | odd $m \to 2m$ | later dyadic rows |
|---|---|---|
| tournament pairs | $2\mathrm{pairs}(m)+1$ | honest doubling |
| LRC units | $\phi(2m)=\phi(m)$ | honest doubling |

The same first seam simultaneously creates tournament twin room and LRC
nonunit quotient room.  This is why first-even LRC denominators such as $14$
and $18$ should be treated as parity-seam export problems rather than ordinary
row doublings.

---

## Connection 5: HYP-217 as 2-Adic Orbit Structure

**HYP-217** (INV-149): for circulant tournament $C_n^S$ with $|S| = 2$:
$\beta_2 = 0$ iff $S = \{s, 2s \bmod n\}$ (a "doubling pair").

A doubling pair $\{s, 2s \bmod n\}$ is a **2-adic orbit in $\mathbb{Z}/n\mathbb{Z}$**:
the orbit of $s$ under the map $x \mapsto 2x \bmod n$. For $|S| = 2$, the condition
says $S$ is a complete 2-adic orbit of size 2.

The column family of $n$ is defined by the same map $x \mapsto 2x$ on sizes. The
2-adic orbit in $\mathbb{Z}/n\mathbb{Z}$ is the internal counterpart of the external
column structure. This isomorphism suggests:

**HYP-217 may follow from showing that the column family structure of $n$ induces
a specific chain-complex decomposition whenever the connection set $S$ contains a
full 2-adic orbit**, and that such orbits supply exactly the chains needed to fill
$\ker(\partial_2)$.

For the density threshold generalization (INV-149: $\beta_2 = 0$ for $|S|$ large enough):
the threshold should be related to the number of complete 2-adic orbits needed to
cover $\mathbb{Z}/n\mathbb{Z}$, which depends on $v_2(n)$ (the row in the family grid).

---

## Connection 6: The 2-Adic SIZE Structure vs. the 2-Adic SPECTRAL Structure

The project already has deep results on 2-adic structure in the Walsh amplitudes
(T241, T243, T244):
- $v_2(\text{amplitude})$ follows the ruler sequence across the Walsh spectrum
- Spectral Legendre: $v_2$ spread = $v_2((n-3)!)$
- Total $v_2$ weight = $k^2 + A000788(k-1)$ (new OEIS sequence)

These are about the 2-adic structure of **spectral invariants** for fixed $n$.

The present reflection is about the 2-adic structure of the **tournament size** $n$
itself. These are dual levels:

| Level | Object | 2-adic structure |
|-------|--------|-----------------|
| External | Tournament size $n$ | Column families $n = 2^r(2k-1)$ |
| Internal | Walsh amplitudes for fixed $n$ | Ruler sequence in spectrum |

The connection: the spectral ruler sequence at fixed $n$ and the family ruler sequence
across sizes $n$ are both instances of the same 2-adic self-similarity. The "binary
representation of dimensions controls 2-adic structure" (T241) at the internal level,
while the "2-adic valuation of $n$ controls family structure" at the external level.

---

## Summary Table

| Repo concept | 2-adic column family interpretation |
|---|---|
| Mode B recursion ($n \to n-2$) | Left-column step in top row |
| Mode A recursion ($n \to n-1$) | Column-crossing, parity-dependent |
| Walsh degree jump odd/even | $r=0$ vs $r \geq 1$ seam |
| SC∩SF = SC($n-2$) identity | Adjacent top-row column bijection |
| SF operation | Column-stripping via middle subtournament |
| HYP-217 doubling-closed condition | 2-adic orbit in $\mathbb{Z}/n\mathbb{Z}$ |
| Tournament blowup $T[K_2]$ | Row step: $(r,k) \to (r+1,k)$ |
| Pairs anomaly $\lfloor n/2 \rfloor$ | Unmatched vertex gaining twin at seam |
| Spectral ruler sequence (T241-T244) | Internal 2-adic structure at fixed $n$ |

**Every major Mode B recursion in this project is a column step. Every tournament
blowup is a row step. The fractal self-similarity of the staircase IS the 2-adic
grid, realized geometrically.**

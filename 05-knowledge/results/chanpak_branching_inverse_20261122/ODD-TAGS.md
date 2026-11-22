# Odd-valuation tags and a physical recurrent two-port inverse

**Status: PROVED.** The two output branches are disjoint on actual ordinary Hamilton counts modulo 15, both branches remain coprime to 15, and the *pair* of ordinary readouts has a complete exact inverse. This is not a decoder from the first count alone, and not an all-target Chan–Pak theorem.

## 1. Physical state, borrowed method, and the selector

We use the existing [root migration](../chanpak_bank_fusion_20261117/migration/PROOF.md) with a marked terminal where exactly one full Hamilton path ends. Record $h=H$, $R$ the complementary two-path count with a nonempty designated path ending at the terminal, and $S$ the number of full paths starting at that terminal. Start with the last vertex of a transitive pair, $(h,R,S)=(1,2,0)$.

A migration has width $k$ and integer binary address $u$, padded to width $k$. It uses an ordered new block and a new terminal, with the one return arc from the distinguished old terminal. Let $b_k(u)$ be $k$-bit reversal. The inherited $S=0$ transfer is
$$
 M(k,u)=
 \begin{pmatrix}
 F(u)&u+1\\
 2^k-b_k(u)&2^k
 \end{pmatrix},
 \qquad
 F(0)=1,\quad F(2t+1)=F(t),\quad
 F(2t)=F(t)+t+[t\text{ odd}]\quad(t>0).
$$
Every leading-zero address again gives $S=0$ and preserves the unique full end path. Every leading-one address used as a finish, with the preceding $S=0$, gives ordinary output $F(u)h+(u+1)R$, still with a unique full end path. It makes the graph strong: the new terminal leads through the ordered block that dominates the old graph; the unique old full path ending at its terminal lets all old vertices reach the return arc back to the new terminal.

The inspiration from [Pasten–Vidaux, §§4–6](https://arxiv.org/html/2611.24748v1) is to assign different jobs to two linked arithmetic tests and propagate a local odd valuation through a controlled reader. Their rational-pole predicate and conditional height argument are not invoked. The [linked RH/GRH remark](https://alpo.ge/remark.pdf) suggests accumulating a property with repeated lifts and tracking the information they cost. Here the repeated operation is a literal finite-state tournament migration, with exact integer determinants and modular rank computed below; there is no use of an RH or GRH claim.

## 2. An invertible reset over the integers, rank one at the readout modulus

Set
$$
 A=M(3,2)=\begin{pmatrix}3&3\\6&8\end{pmatrix},
 \qquad B=M(4,4)=\begin{pmatrix}5&5\\14&16\end{pmatrix}.
$$
The ready phase is $(h,R)\equiv(0,1)\pmod {15}$. It is physically seeded by $B(1,2)=(15,46)$. Apply $A$ or $B$ to this seed for the first letter. The two tag phases are
$$
 A(0,1)=(3,8),\qquad B(0,1)=(5,1)\pmod {15}.
$$

The following chronological word gives a common reset $W$. Matrices act on columns, so later operations multiply on the left.

| Step | $(k,u)$ | Matrix | Determinant |
|---:|---:|---|---:|
| 1 | $5,10$ | $\begin{pmatrix}9&11\\22&32\end{pmatrix}$ | 46 |
| 2 | $3,2$ | $A$ | 6 |
| 3 | $5,14$ | $\begin{pmatrix}9&15\\18&32\end{pmatrix}$ | 18 |
| 4 | $5,9$ | $\begin{pmatrix}5&10\\14&32\end{pmatrix}$ | 20 |
| 5 | $4,4$ | $B$ | 10 |

Multiplication gives
$$
 W\equiv\begin{pmatrix}0&0\\1&11\end{pmatrix}\pmod{15}.
$$
Thus $W(3,8)=W(5,1)=(0,1)\pmod{15}$, since $3+11\cdot 8=91$ and $5+11=16$. **The integer map has not lost its inverse:** each listed determinant is nonzero, though its composite has rank one modulo 3 and modulo 5. All five addresses lead with zero, so the state is again eligible to branch. This is a reset of coarse tags, not a rank-one reset of the actual rooted state.

After the first bit, each subsequent bit applies $W$ and then chooses $A$ or $B$. Every nonempty binary word is therefore physically legal. On a ready state $h+R$ is a unit modulo 15, while the two first coordinates are $3(h+R)$ and $5(h+R)$. Accordingly the last letter has the **exact** odd-prime valuation tags
$$
 (v_3(h),v_5(h))=(1,0)\ \text{for }A,\qquad
 (v_3(h),v_5(h))=(0,1)\ \text{for }B.
$$
This is not merely divisibility.

## 3. The ordinary outputs preserve useful tags without forced compositeness

Use two leading-one finish addresses $(5,18)$ and $(6,37)$. Since $F(18)=F(37)=15$, they give the pair of strong ordinary tournament counts
$$
 G=15h+19R,\qquad G'=15h+38R.
$$
The two addresses differ by appending a $1$ bit ($37=2\cdot18+1$); the $F$ coordinate is unchanged and the $R$ coefficient doubles. Both outputs are coprime to 15 on either tag phase. The first read by itself visibly separates the final letters:
$$
 G\equiv 2\pmod {15}\text{ after }A,\qquad
 G\equiv 4\pmod {15}\text{ after }B.
$$
The exact valuation tags survive in the linked read $2G-G'=15h$: they are $(2,1)$ versus $(1,2)$ there. There is no primality assertion, only absence of forced factors 3 and 5.

The **complete inverse on pairs** is explicit. First recover
$$
 h=(2G-G')/15,\qquad R=(G'-G)/19,
$$
rejecting failed positive integer divisions. Select the last letter by $G\bmod15$, equivalently validating the exact tag on $h$. Apply its rational matrix inverse and reject if the result is not an integral positive ready state. If that result is the seed, stop. Otherwise undo $W$, checking integrality and positivity again, and recurse on the preceding tag state. Accept only after reaching the seed; all branch steps strictly increase the first coordinate on positive states, so the exact inverse terminates. The forward replay is a final membership certificate. Surjectivity onto the image and injectivity on physical words follow. This inverse uses both $G$ and $G'$: final-letter separation on the first alone is not a full scalar decoder.

For a scale check, the first-letter tagged states are $A(15,46)=(183,458)$ and $B(15,46)=(305,946)$. They give $(G,G')=(11447,20149)$ and $(22549,40523)$; the mod-15 tags are visible on both first reads. These are actual strong tournaments with a retained physical terminal. The alphabet is binary, but its widths are large relative to one bit, so no density or arbitrary absolute target claim is made.

## 4. Joint CRT steering for every odd modulus

**Theorem.** Let $m$ be any odd positive integer and $(x,y)$ any residue vector modulo $m$. There are infinitely many distinct ordinary Hamilton counts of strong rooted tournaments with $(h,R)\equiv(x,y)\pmod m$ and exactly one full path ending at the marked terminal. In particular one can prescribe exact positive or zero valuations of both coordinates independently at any prescribed finite set of odd primes. This does not prescribe an absolute magnitude, any $2$-adic digits, or an efficient word-length bound.

### A contact-word loop changes one matrix entry modulo $m$

Suppose $m>1$, take $L\ge1$ with $2^L\equiv1\pmod m$, and let $u$ be an odd address padded to width $k$. Append $L$ zero bits, then $L$ one bits. After the zero block, the recurrence gives
$$
 F(2^L u)=F(u)+(2^L-1)u+1.
$$
Appending ones does not change $F$, and restores odd parity. The resulting address and width are
$$
 u'=(2^L u+1)2^L-1,\qquad k'=k+2L.
$$
Thus $u'+1\equiv u+1$ and $2^{k'}\equiv2^k\pmod m$. The old reversed bits retain their exponents, while the appended ones add $2^{k+L}(2^L-1)$ to the reversal, zero modulo $m$. Hence this loop increments $F$ by one modulo $m$ and leaves the other three entries of the transfer unchanged modulo $m$. It preserves the leading bit of the padded word.

Start from address $1$ at width $2$ (padded word $01$). Using one and two loops gives two actual leading-zero addresses whose matrices reduce to
$$
 P_1=\begin{pmatrix}2&2\\2&4\end{pmatrix},\qquad
 P_2=\begin{pmatrix}3&2\\2&4\end{pmatrix}\pmod m.
$$
Their determinants modulo $m$ are 4 and 8, both units. Another physical generator comes by adding one left zero to the address for $P_1$; it has matrix $D P_1$, where $D=\operatorname{diag}(1,2)$ (left padding doubles the full second row, exactly). Every inverse modulo $m$ can be replaced by a positive power of the physical matrix, since the group is finite. In particular $D=(D P_1)P_1^{-1}$ is available on residues.

### Primitive steering and a strong finish

Direct multiplication gives
$$
 J=P_2P_1^{-1}=\begin{pmatrix}2&-1/2\\0&1\end{pmatrix},
 \qquad D J D^{-1}J^{-1}=E_{12}(1/4).
$$
Powers of this commutator give every upper shear over $\mathbb Z/m$, since $1/4$ is a unit and the underlying additive group is cyclic. For $P=E_{12}(-1)P_1$, the vector $P e_1$ is $2e_2$, so conjugating the upper shears by $P$ gives every lower shear. Explicitly, $P=\begin{pmatrix}0&-2\\2&4\end{pmatrix}$ and $e_2^{\mathsf T}P^{-1}=-e_1^{\mathsf T}/2$, so $P E_{12}(t) P^{-1}=E_{21}(-t)$.

The two shear types act transitively on primitive vectors modulo $m$. One can first add a suitable multiple of the lower coordinate to the upper to make the latter a unit (choose separately at each prime, then use CRT). Next eliminate the lower coordinate and normalize the upper by a unit diagonal matrix, expressible in shears:
$$
 \operatorname{diag}(u,u^{-1})=
 (E_{12}(u)E_{21}(-u^{-1})E_{12}(u))
 (E_{12}(-1)E_{21}(1)E_{12}(-1)).
$$
The seed $(1,2)$ is primitive, so positive leading-zero words can steer it to every primitive vector.

For a leading-one strong finish, start from unpadded width $1$ address $1$ (word $1$) and apply the same contact-word loop once. Its matrix on incoming $S=0$ reduces to $\begin{pmatrix}2&2\\1&2\end{pmatrix}$, with determinant 2. Thus it is possible to finish strongly *and* hit any given primitive residue vector by steering to its preimage first.

### Nonprimitive targets and unbounded outputs

For a general target take integer representatives and $d=\gcd(x,y,m)$. First use primitive leading-zero steering to obtain $(h,R)\equiv(1,d-1)\pmod m$. The physical one-bit address $0$ applies the rank-one matrix $\begin{pmatrix}1&1\\2&2\end{pmatrix}$, giving $d(1,2)$ modulo $m$. Now use a primitive steering word, with its strong finish, modulo $m/d$ to take $(1,2)$ to $(x/d,y/d)$, which is primitive there. Apply this physical word to the actual reset state. All transfers are homogeneous integer maps on the incoming $(h,R)$, so the final residue modulo $m$ is $(x,y)$. When $m/d=1$, any strong finish works. This latter use of rank one does intentionally lose the old state, in contrast to the nonsingular tag reset of Section 2.

In either case, before the final finish insert indefinitely many copies of any positive leading-zero matrix word having unit determinant modulo the working odd modulus, raised to its finite group order. This is the identity on residues and strictly increases $h$ and $R$ on positive inputs (the constructed physical $P_1$, or any positive leading-zero migration when the modulus is 1, suffices). The ordinary finish count therefore tends strictly to infinity. For exact local valuation targets, specify residues modulo one power higher than each desired valuation and use CRT.

Replay of the new transfers and finite controls: `python3 check_odd_tags.py`.

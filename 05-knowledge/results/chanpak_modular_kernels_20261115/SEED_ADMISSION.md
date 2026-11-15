# Seed supply: a sharp paid-core admission rule and a persistent signed lift

This note studies the positive offset $r=H(B)-W_B(O)$ and response $M=U_B(O)$ of a nonempty marked pair $(B,O)$. Here $W$ counts two-path covers whose first path ends and second starts in $O$, and $U$ requires just the marked start of the second path. The [inherited closing law (THM-5126)](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md) realizes $r+2^tM$ for every $t\ge1$. A perfect binary response means $M$ is an exact power of two, not just an even integer.

## Exact decomposition: the marked core has to pay

Let $B=A\to C\to S$, with $C$ one strong component. The [inherited strong-component rules (THM-5130(F))](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md) give the following exact admission criterion for $r>0$ and a dyadic response:

1. $O$ is confined to a single SCC $C$; marks in distinct SCCs already give $W\ge H(B)$.
2. All SCCs in $A$ are singletons.
3. $U_C(O)$ and each two-slot optional cover count $P_D(2)$ in the unmarked suffix $S$ are powers of two.

Conversely these conditions and $r_C=H(C)-W_C(O)>0$ suffice. The resulting offset and response are
$$r=r_C\,H(S),\qquad M=U_C(O)\,P_S(2).$$
Indeed $r=r_C H(A)H(S)$, $M=H(A)U_C(O)P_S(2)$, and the odd integer $H(A)\mid M$ must equal one. The two-slot count multiplies along ordinal concatenations, so any odd factors in an unmarked suffix survive. Thus the positive-offset seed problem is exactly a one-marked-strong-core supply problem with a separate semigroup of allowed dyadic two-slot suffix factors. It is not the unrestricted multiplicative semigroup of already realized Hamilton counts.

The literal same-count lifting claim, that every realized $r$ can start from $H(B)=r$ and have dyadic response, is false. If $H(B)=r$, then $W=0$, so $O$ is a singleton because each distinct pair has positive coupling. Every tournament with $H=3$ is one directed triangle surrounded by singleton SCCs. For a directed triangle $C$, each $U_C(y)=3$ and $P_C(2)=12$. Whichever side of $C$ a single mark lies on, its response is divisible by three. More generally, if an ordinal block of $B$ is a directed triangle, its factor of three survives in **every** $U_B(O)$ even for many marks and arbitrary ordinal catalysts: expand $U_{X\to Y}=U_XP_Y(2)+H(X)U_Y$, and use $(H(C),P_C(2),U_C(O))=(3,12,3|O\cap C|)$. Marking can break the triangle's automorphism without breaking its arithmetic obstruction.

Same-count **induced** enlargement with the old selected vertex fixed does not solve this problem. By strict-extension detection via odd-cycle packings, no new vertex can lie in a triangle at unchanged $H$; they are singleton SCCs. Their placement around the selected SCC only multiplies its old response by a power of two.

The first paid offset-three example is small: outgoing rows $(20,1,2,7,14)$, mark mask 17, give $(r,M,H)=(3,32,9)$. The paid cost is necessary: for $H=5$ the strong core is the unique strong four-tournament, whose positive offset-three marks have responses 13 or 14 (before dyadic padding), and $H=7$ is excluded. For $H\ge9$ with $r=3$, at least two marks are required and $M\ge2H\ge18$, so no dyadic response below 32 works. Thus this literal witness is sharp both in base count and response (the response value agrees with the prior type census). Appending $k$ unmarked nearly transitive four-tournaments, each with $(H,P(2))=(5,32)$, gives the explicit family
$$r=3\cdot5^k,\qquad M=32^{k+1}.$$
This is an instance of the exact suffix criterion, not a claim of a new generic dyadic shelf.

There is also a universal budget for paid seeds. For each selected $y$, splitting an actual Hamiltonian path immediately before $y$ shows $U_B(y)\ge H(B)$. If $h=H(B)>r>0$, then $h\ge r+2$ and $|O|\ge2$. A nonempty dyadic closing ending at $N=2^e+r$ therefore obeys
$$2^e=2^t M\ge4(r+2).$$
Every such paid cone lands strictly below $(5/4)2^e$, in the leading binary band **100**. A bare existence theorem for paid seeds cannot by itself handle the other leading-digit bands. This sharpens the familiar positive-offset cone barrier by keeping whether the base has been promoted.

## A closed, physical signed seed lift for every realized count

The obstruction is to an **absolute** dyadic response. There is nonetheless a closed signed transport of arbitrary already realized counts, preserving physical synchronization of two positive trajectories.

Use the four-coordinate state $(r,M,h,E)$ of THM-5130, with $E$ the number of Hamiltonian paths of $B$ ending in $O$. The initial literal pairs are:

| outgoing rows of $B$ | mark mask | $(r,M,h,E)$ |
|---|---:|---:|
| $(18,0,3,7,14)$ | 17 | $(1,24,5,0)$ |
| $(10,0,3,6,15)$ | 16 | $(3,24,3,0)$ |

They agree in $(M,E,h+r)$, but the second offset exceeds the first by two. The common gadget update of THM-5130 has the form
$$
\begin{aligned}
r'&=r_Xr+\kappa_X M,& h'&=H_Xr+P_X(2)M,\\
M'&=u_X(h+r)+\beta_XM-\sigma_X E,& E'&=\epsilon_X(r-h)+\nu_XM .
\end{aligned}
$$
Given any actual nonempty tournament $T$ with $c=H(T)$ and $Q=P_T(2)$, take the gadget $X=T\to TT_3$ and mark the **first two vertices of the terminal triple** only. Then
$$
H_X=c,\quad r_X=-c,\quad \kappa_X=0,\quad
(\sigma_X,\epsilon_X)=(0,0),\quad
P_X(2)=8Q,\quad u_X=6c,\quad \beta_X=30Q,\quad\nu_X=3Q.
$$
To see neutrality, $W_X=2c$. When the third slot is allowed empty, the three-slot $W$ count equals $4Q$: the second slot starts in the marked triple and the copies of $T$ can go in the first and third slots. The terminal triple's marked endpoint count in these optional slots is four, by enumerating the ordered endpoints and position of its last vertex. Subtract $2c$ for an empty third slot; then $\kappa_X=8Q-2W_X-2(4Q-2c)=0$. For $\beta_X$ use twice the optional three-slot start count: the marked first and second positions have $9$ and $6$ completions in $TT_3$ respectively, and the unmarked prefix contributes $Q$. The other identities are immediate from ordinal marks (no Hamiltonian path starts or ends at the marked positions in $X$). These data yield the exact transfer
$$r'=-cr,\quad h'+r'=8QM,\quad
M'=6c(h+r)+30QM,\quad E'=3QM.$$

Thus agreement in $(M,E,h+r)$ is preserved, and each offset is transported by the exact multiplier $-c$ with no common offset contamination. Applying arbitrarily many such physical gadgets with realized counts $c_1,\ldots,c_d$ yields offsets $(-1)^d\prod c_i$ and $3(-1)^d\prod c_i$ with the **same response**. The two literal closing tournaments of any common tail length then differ in their Hamiltonian counts by
$$2(-1)^d\prod_{i=1}^d c_i.$$
This is a repeatable signed lane for every realized multiplier, unlike treating the difference of unrelated host counts as a new seed. It does not turn that signed difference into one absolute count, or claim that the shared response is dyadic. For example the first common response is $36c+720Q$; it is emphatically not a power of two. Each neutral gadget also retains the obstacle $\gcd(c,Q)\mid M'$, in agreement with the triangle barrier.

The literal replay verifies the small paired states, the sharp offset-three case, the transfer against direct subset path counts, and both positive outputs of each example. Chan–Pak and unrestricted paid-core supply remain open.


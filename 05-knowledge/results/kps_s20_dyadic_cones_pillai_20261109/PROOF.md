# THM-5109: dyadic cones, the Pillai chain law and critical fixed-port sheets (proofs)

kind-pasteur-2026-11-09-S20. Base 2 throughout (H = I(Omega, 2) by OCF).

**Priority.**
- Lemma A and Theorem B (sections 1 and 2) are KNOWN. They are THM-4879 (A) and (B), from this machine's own S11
  (2026-11-04, audited), stated there in SEAM's converse convention. Under the converse, gadget position i becomes
  m - i and SEAM's S becomes e. Their R-row coefficient c = 2^(k-z+u) 3^(z-u-1) (z > u) is exactly
  3^(p-q-1) 2^(m-p+q) (q < p) below.
- The proofs are reproduced in the direct convention only because sections 4-5 use them. The pre-audit draft typed
  them as new; that is MISTAKE-750.
- The cone of section 3 is also KNOWN: Session236 offsets section 2, and this machine's THM-5085 (3). The pre-audit
  draft typed it as new; that is MISTAKE-751.
- New in this file:
  - the window bound and the census in section 3;
  - Theorem D;
  - Theorem E.
- Theorem F is an elementary reformulation of known decoders.
- The independent audit (`audit/AUDIT.md`) found no error in a proof. Its corrections, all applied, concerned
  priority, scope and summaries.

## 0. Notation and inherited facts

**Rooted tournaments.** A rooted tournament is a pair (T, v). Its triple is (h, e, R), where:
- h = H(T), with H(empty) = 1;
- e = e_v(T) is the number of Hamiltonian paths of T that end at v;
- R = B_v(T) = sum over S containing v of s_v(T[S]) H(T[V - S]), where s_v(X) is the number of Hamiltonian paths of X
  that start at v. Equivalently, R is the sum of H(T - V(P)) over all paths P of T that start at v, including P = (v).

R is the root response B_t of THM-4865 and B of HYP-9453. In the fixed-port notation of THM-5085, R = E + O and
e = O - E. This is the identification of SEAM.md's (H, S, R) with (h, O - E, E + O) recorded by the S18 audit.
- The identity e = O - E was rechecked on 468 ports of random tournaments of orders 1-9, with 0 failures
  (`check_e_identity.py`).
- The proofs below do not use it. The census and P4 inputs read e from the S19 port statistics as O - E.

**Proposition 0.1 (B >= h).** R >= h.
- *Proof.* Cut every Hamiltonian path Q of T just before v, so Q = Q1 Q2 with Q2 starting at v. Then Q2 is a path from
  v and Q1 is a Hamiltonian path of T - V(Q2). The map Q -> (Q2, Q1) is injective, so
  h <= sum over P from v of H(T - V(P)) = R. QED
- This inequality is stated in HYP-9453; THM-4875 (C) uses R - h >= 0.

**Two-path cover count.**

    Q(X) = sum over S subset V(X) of H(X[S]) H(X[V(X) - S]),    Q(empty) = 1.

- Q(X) is the response of a new source u added above X: R(u => X, u) = Q(X). The paths from u are u followed by a
  path of X. This is Session229's U-formula (`overlap/PROOF.md` section 1, eq. (1)) for a root that beats every block;
  it telescopes to the product of the blocks' Q.
- Q is Session229's proper-graft quantity and the Q of THM-5107 (Session242). Session242 proves 4 | Q for |X| >= 2.
- Q is multiplicative over ordinal sums, Q(X => Y) = Q(X) Q(Y), because H is. This is KNOWN: Session227
  `multiroot/PROOF.md` section 5, and Session229 `dynamic/PROOF.md` eq. (6).

**One-arc join** (SEAM.md; THM-4875 (B)).
- Construction: N = T => G with the single arc v -> u reversed to u -> v, for an entry u of G. The new root is the
  exit z in V(G).
- Transfer: (h, e, R)(N, z) = M(G, u, z) (h, e, R)(T, v), where M = [[h_G, -s_u, E_u], [e_z, -P, E_uz], [R_z, -J_1,
  J_2 + K]].
- First two rows, all path counts in G:
  - h_G = H(G);
  - s_u = number of Hamiltonian paths of G starting at u;
  - E_u = sum of H(G - V(P)) over paths P of G ending at u;
  - e_z = number of Hamiltonian paths of G ending at z;
  - P = number of Hamiltonian paths of G from u to z;
  - E_uz = number of pairs (B1, B2) of disjoint paths covering G with B1 ending at u and B2 nonempty ending at z.
- These two rows are the first two rows of THM-4875 (B), and they follow from the same cut argument as Lemma A below.

## 1. Lemma A (the response row as path sums; KNOWN: THM-4879 (A))

**Lemma A.**

    R(N, z) = B_z(G) h - J_1 e + (J_2 + K) R,

where:
- B_z(G) is the response of z in G;
- J_1 = sum over paths C of G from z avoiding u of s_u(G - V(C));
- J_2 + K = sum over paths C of G from z avoiding u of E_u(G - V(C))
          + sum over paths C of G from z ending at u of Q(G - V(C)).

Here E_u(X) is the E_u of the gadget X with entry u.

**Proof.**
- **Where paths can go.** The only arc from G to T is u -> v. So a path of N that starts at z is one of two kinds:
  - (I) a path C of G from z;
  - (II) a path C of G from z ending at u, followed by v, then a path D of T from v (beginning at v), then a possibly
    empty path C2 of G - V(C).
- **Type (II) concatenations are paths.** Every arc from a vertex of T to a vertex of G other than u is present, and u
  is already used.
- **Type (I) remainders.**
  - If u is not in C, the remainder is the one-arc join of (T, v) with (G - V(C), u). By the first row its H is
    h H(G - V(C)) - e s_u(G - V(C)) + R E_u(G - V(C)).
  - If u is in C, the remainder is the ordinal sum T => (G - V(C)), with H = h H(G - V(C)).
- **Type (II) remainders.** The remainder is (T - V(D)) => (G - V(C) - V(C2)), with H = H(T - V(D)) H(G - V(C) - V(C2)).
  - Summing over D gives R.
  - Summing over C2, the empty path included, gives Q(G - V(C)).
- **Collecting.** Summing type (I) over all C gives h B_z(G), and the remaining terms give the coefficients of e and R
  as displayed. QED

## 2. Theorem B (transitive gadgets; KNOWN: THM-4879 (B))

**Setup.** Let G = TT_m with vertices g_1 < ... < g_m, where g_i -> g_k for i < k. Let u = g_p and z = g_q. Write
(h, e, R) for the base triple.

**Theorem B.**

    h' = h + 2^(p-1) R - [p = 1] e,
    e' = [q = m] h - [p = 1][q = m] e + E_uz R,

where E_uz = 2^(p-1) if q = m > p, E_uz = 2^(q-1) if q < p = m, and E_uz = 0 otherwise.

The new response R' is:

| entry | exit | R' |
|---|---|---|
| p = 1 | q = 1 | 2^(m-1) (h + R) |
| p = 1 | q >= 2 | 2^(m-q) h' |
| p >= 2 | q >= p | 2^(m-q) h' |
| p >= 2 | 2 <= q < p | 2^(m-q) h + 3^(p-q-1) 2^(m-p+q) R |
| p >= 2 | q = 1 | 2^(m-1) h - 2^(m-p) e + 3^(p-2) 2^(m-p+1) R |

**Proof.**

*Gadget statistics.* In TT_m:
- each vertex set carries exactly one path, its increasing order, and every remainder is transitive with H = 1;
- s_u = [p = 1] and E_u = 2^(p-1): a path ending at g_p is any subset of {g_1, ..., g_(p-1)} followed by g_p.

*The value h'.* This is the first row.

*The row e'.*
- e_z = [q = m] and P = [p = 1][q = m].
- E_uz counts pairs (B1, B2). Here B1 is an increasing set ending at g_p, and B2 is the complementary set, which must
  end at g_q. So every vertex after g_q must lie in B1, whose largest vertex is g_p.
  - If q > p, this forces q = m, and B1 is any subset of the first p - 1 vertices plus g_p: 2^(p-1) pairs.
  - If q < p, this forces p = m, and B1 contains every vertex strictly between g_q and g_p plus any subset of the first
    q - 1 vertices: 2^(q-1) pairs.

*The response row* (Lemma A).
- **B_z(G).** The paths from g_q are g_q followed by any subset of later vertices, so B_z(G) = 2^(m-q).
- **Exit at or after the entry (q >= p, p >= 2).**
  - For q > p, no path from z reaches u. Each C avoids u with E_u(G - V(C)) = 2^(p-1). This gives
    J_2 + K = 2^(m-q) 2^(p-1).
  - J_1 = 0 because g_1 survives and precedes u (p >= 2).
  - For q = p, the only path ending at u is (u), with Q(G - u) = 2^(m-1), and J_1 = 0.
  - Both cases give R' = 2^(m-q) h'.
- **Exit before the entry (q < p).**
  - *Paths C from g_q avoiding u.* Each is g_q, then a subset S of the vertices strictly between g_q and g_p, then a
    subset of the vertices after g_p.
    - The remainder keeps q - 1 + (p - q - 1 - |S|) vertices before u, so E_u = 2^(q-1) 2^(p-q-1-|S|).
    - Summing over S and over the 2^(m-p) tails gives 2^(q-1) 3^(p-q-1) 2^(m-p).
  - *Paths C from g_q ending at u.* Each is g_q, then S, then u. The remainder is transitive with m - 2 - |S|
    vertices, so Q = 2^(m-2-|S|). Summing over S gives 2^(m-p+q-1) 3^(p-q-1).
  - Together these give J_2 + K = 3^(p-q-1) 2^(m-p+q).
  - *The coefficient J_1.* u is the source of G - V(C) only if every vertex before u lies in C.
    - If q >= 2, g_1 never does, so J_1 = 0.
    - If q = 1, the paths C covering g_1, ..., g_(p-1), followed by any subset after g_p, give J_1 = 2^(m-p).
- **Entry at the source (p = 1).**
  - For q >= 2: J_1 = 2^(m-q) (u = g_1 is the source of every remainder that contains it), the second sum of J_2 + K
    is empty, and J_2 + K = 2^(m-q). So R' = 2^(m-q)(h - e + R) = 2^(m-q) h'.
  - For q = 1 (z = u): J_1 = 0, and C = (u) contributes Q(G - u) = 2^(m-1). So R' = 2^(m-1)(h + R).

QED

**Checks** (corrected after audit minor 1).
- `check_tt.py`: 24 random bases of order <= 6, every (m, p, q) with m <= 7, 3,360 cases. It compares h' and e' on all
  3,360, and R' on the 2,688 cases with p >= 2 (`out/check_tt.out`).
- `check_tt_p1.py`: 30 bases, m <= 8, 1,080 cases with p = 1. It compares h' and R' (`out/check_tt_p1.out`).
- In all, h' was compared on 4,440 cases, e' on 3,360 and R' on 3,768, with 0 mismatches.
- The audit's own counter checks every row, p = 1 included, on 7,956 cases with 0 mismatches. It also matches
  THM-4879 (B) read verbatim through the converse map (`audit/a01_validate_and_tableB.out`).

## 3. Corollary C (dyadic cones: KNOWN mechanism) and the window bound (new)

**Statement (KNOWN).** Suppose R = 2^j (a *dyadic port*). For every K >= j + 1, let p = K - j + 1 >= 2 and attach TT_m
with any m >= p. This gives H = h + 2^K.

**Priority (audit M1; MISTAKE-751).**
- This is the near-join h + 2^k B of HYP-9453 and THM-4865, specialised to B = 2^j.
- It is already in Session236 `offsets/PROOF.md` section 2 ("Clock shifts A = 2^r"), with the seeds (5, 0, 32) and
  (9, 2, 16) and the rays 2^(m+5) + 5 (69, 133, 261, ...) and 2^(m+4) + 9.
- It is also in this machine's THM-5085 (3) (S19 PROOF sections 2 (b) and 3, "pure dyadic rays"). There, s sinks
  followed by one bridge is literally the one-arc join with TT_(s+1) entered at its sink, that is, the cone with
  p = m.
- Session236 also has the lemma "R = 1 iff T is transitive and rooted at its sink", used in (E2).

**Proposition C.1 (window bound).** Every cone value from a dyadic port lies in (2^K, 3 * 2^(K-1)]. Hence an odd cone
value has binary expansion beginning 10. No cone, from any dyadic port of any tournament, reaches an odd N in
(3 * 2^(K-1), 2^(K+1)): for example 7, 13, 15, 25-31, 49-63.

*Proof.* By Proposition 0.1, h <= R = 2^j, and K >= j + 1 gives 2^j <= 2^(K-1). So 2^K < h + 2^K <= 2^K + 2^(K-1).
QED

Compare THM-5055's even-graph remark that values just below 2^(n-1) are never cones. This was suggested by the audit
(M4).

**Source ports.** If Q(X) = 2^i, then the new source over X is a dyadic port of value H(X), and Q is multiplicative
(both inherited, section 0). The strong tournaments of order 3-8 with Q a power of two are exactly 11 classes
(`out/q_dyadic_strong.out`; reproduced by audit a03). The trivial order-1 tournament also has Q = 2.
- The pairs (5, 32) and (19, 256) are known (Session227 multiroot section 5; Session228/229).
- The values 217, 223 and 599 are new data.

| order | H | Q |
|---|---|---|
| 4 | 5 | 2^5 |
| 6 | 19 | 2^8 (two classes) |
| 8 | 217 | 2^12 (two classes) |
| 8 | 223 | 2^12 (two classes) |
| 8 | 599 | 2^13 (four classes) |

**Consequence.** Let X be an ordinal sum of these tournaments and transitive tournaments. Then H(X) + 2^K is a
Hamiltonian-path count for every K > log_2 Q(X). For instance 5^a 19^b + 2^K for every K >= 5a + 8b + 1; for a = 1,
b = 0 this gives 69, 133, 261, ...

## 4. Theorem D (the chain law)

**Statement.** Take a dyadic port (h, e, 2^j) and a transitive gadget with p >= 2. The exit g_q is again a dyadic port
if and only if q < p and one of the following holds:
- **(i)** q >= 2, p = 2q + j, and h + 3^(q+j-1) is a power of two, 2^t. Then R' = 2^(m-q+t).
- **(ii)** q = 1, and 2^(p-1) h - e + 3^(p-2) 2^(j+1) is a power of two, 2^s. Then R' = 2^(m-p+s).
  - **(ii-a)** When e = 0, this holds if and only if p = j + 2 and h + 3^j is a power of two.

**Proof.**
- **q >= p.** Then R' = 2^(m-q) h' with h' >= 1 + 2^(p-1) odd, so R' is not a power of two.
- **Case (i).** The two summands 2^(m-q) h and 3^(p-q-1) 2^(m-p+q+j) have 2-adic valuations m - q and m - p + q + j.
  - If these differ, the sum is 2^(min) times an odd number at least 3, so it is not a power of two.
  - If they agree, p = 2q + j, and the sum is 2^(m-q)(h + 3^(q+j-1)).
- **Case (ii).** Factor R' = 2^(m-p)(2^(p-1) h - e + 3^(p-2) 2^(j+1)).
- **Case (ii-a).** With e = 0, the valuations of the two summands are p - 1 and j + 1. Unequal valuations give an odd
  cofactor of at least 3. Equal valuations give 2^(j+1)(h + 3^j).

QED

## 5. Theorem E (Pillai rigidity)

### (E1) Pillai near-misses

The solutions of 1 <= 2^x - 3^y <= 2^y in nonnegative integers are exactly (x, y) = (1, 0), (2, 1), (5, 3), (8, 5).
The values 2^x - 3^y are 1, 1, 5, 13.

*Proof.*
- **Ellison's bound.** Ellison (1971; Séminaire de Théorie des Nombres de Bordeaux 1970-71, exposé 12) proved
  |2^x - 3^y| > 2^x e^(-x/10) for all y and all x >= 12 except x = 13, 14, 16, 19, 27. This statement is as quoted in
  M. Waldschmidt, *Perfect powers: Pillai's works and their developments*, arXiv:0908.4031, section 2.
- **Large x.** For x >= 28:
  - 2^x - 3^y > 2^(x (1 - 1/(10 ln 2))) > 2^(0.855 x);
  - 3^y < 2^x gives 2^y < 2^(0.631 x).
  - These contradict 2^x - 3^y <= 2^y.
- **Small x.** `h1_list.py` enumerates x <= 60 exhaustively and finds the four solutions.
- **Sanity check.** It also confirms Ellison's inequality for every 28 <= x <= 60.
- **Audit check.** The audit's own download of Waldschmidt confirms the quote, in section 2.2.2. Numerically the
  inequality fails for x >= 12 exactly at x = 13, 14, 16, 19, 27, up to x = 20,000 (`audit/a05_pillai.out`).

QED

### (E2) Classification of single dyadic cone steps

Every dyadic cone step through a transitive gadget whose new port lies in the gadget is one of the following:
- **(alpha)** 1 -> 9, either from (1, 1, 1) with (p, q) = (4, 2), or from (1, 0, 2) with (p, q) = (3, 1).
- **(beta)** a step of type (ii) with e > 0, p in {3, 5, 7}, and

      e = 2^(p-1) h - d_p 2^(j+1),    (d_3, d_5, d_7) = (1, 5, 13),

  where 0 < e <= h <= 2^j. The output is h' = h + 2^(p+j-1) and R' = 2^(m-p+s), with 2^(s-j-1) = 3^(p-2) + d_p.

*Proof.*
- **Type (i) and type (ii-a).** Write a = q + j - 1 >= j. By B >= h (Proposition 0.1), 2^j >= h, so
  h = 2^t - 3^a <= 2^a. By (E1), (t, a, h) is one of (1, 0, 1), (2, 1, 1), (5, 3, 5), (8, 5, 13).
  - **h = 1.** a = 0 requires q = 1 and j = 0. But R = 1 forces T to be transitive and rooted at its sink, so e = 1
    and type (ii-a) does not apply. The case a = 1 gives exactly the two steps of (alpha).
  - **h = 5 and h = 13.**
    - Since 5 and 13 are prime, a tournament with H in {5, 13} has exactly one nontrivial strong component C, with
      H(C) = H. THM-4810 (A) shows that the strong fibres F_s(5) and F_s(13) are singletons over all orders.
    - With transitive padding, the root responses are:
      - before C: Q(C) 2^i;
      - inside C: B_w(C) 2^i;
      - after C: H(C) 2^i.
    - The census (`out/dyadic_census_le8.out`) contains both singleton classes with a source above them, and gives:
      - H = 5: Q = 2^5 and in-component responses {5, 6, 7, 9};
      - H = 13: response odd parts {7, 9, 13, 15}, with no power of two.
    - So a dyadic port of value 5 has j >= 5, while a step needs j <= a = 3 (q >= 2 would need j <= 2). A
      tournament of value 13 has no dyadic port at all.
- **Type (ii) with e > 0.**
  - Since 0 < e <= h <= 2^j, the power 2^s lies strictly between 3^(p-2) 2^(j+1) and
    2^j (2^(p-1) + 2 * 3^(p-2)).
  - Hence 2^(s-j-1) lies strictly between 3^(p-2) and 3^(p-2) + 2^(p-2).
  - By (E1) with y = p - 2 and the strict upper inequality, (x, y) is one of (2, 1), (5, 3), (8, 5).
  - This gives the three values of p and d_p, and solving the condition for e gives the displayed relation.

QED

### (E3) No two dyadic cone steps chain

No dyadic port admits two consecutive dyadic cone steps through transitive gadgets with both new ports inside the
attached gadgets.

*Proof.*
- **After (alpha).**
  - The output value 9 is not in {1, 5, 13}.
  - The output has e' != 0 only when m = p. Both realisations then give (h', e', j') = (9, 2, 4):
    - from (1, 1, 1) with (m, p, q) = (4, 4, 2), R' = 2^(m-2) + 3 * 2^(m-2) = 2^m = 2^4 and e' = 2^(q-1) R = 2;
    - from (1, 0, 2) with (m, p, q) = (3, 3, 1), R' = 2^(m-1) + 3 * 2^(m-1) = 2^(m+1) = 2^4 and e' = 2^(q-1) R = 2.
  - The (beta) relation e' = 2^(p'-1) 9 - d_(p') 2^(j'+1) then fails for p' = 3, 5, 7: the right-hand sides are 4, -16
    and 160.
  - For m > p the output has e' = 0, so a (beta) step is impossible, and 9 is not in {1, 5, 13}.
- **After (beta).**
  - **A type (i) or (ii-a) second step.** It needs h' in {5, 13}, with h' = h + 2^(p+j-1) and p >= 3. The only
    candidates are (h, j, p) = (1, 0, 3), (5, 1, 3) and (9, 0, 3). Each violates e <= h in the (beta) relation.
  - **A (beta) second step.**
    - It needs e' > 0, so m = p, e' = 2^j and j' = s = x + j + 1.
    - The relation then reads 2^(p'-1) h = 2^j + d_(p') 2^(x+j+2) - 2^(p+p'+j-2).
    - The bracket 1 + d_(p') 2^(x+2) - 2^(p+p'-2) is odd, so oddness of h forces j = p' - 1 and
      h = 1 + d_(p') 2^(x+2) - 2^(p+p'-2).
    - Here 2^x = 3^(p-2) + d_p, so x = 2, 5, 8 for p = 3, 5, 7. The nine pairs (p, p') with p, p' in {3, 5, 7} give:

      | p \ p' | 3 | 5 | 7 |
      |---|---|---|---|
      | 3 | 1 | 17 | -47 |
      | 5 | 65 | 385 | 641 |
      | 7 | 769 | 4097 | 9217 |

    - The constraint h <= 2^j = 2^(p'-1) (that is, h <= 4, 16, 64) leaves only h = 1 at (p, p') = (3, 3). Its first
      step would need e = 2^(p-1) h - d_p 2^(j+1) = 4 - 8 < 0. So no pair survives.

QED

**Abstract triples (audit minor 4).** For abstract triples (h, e, 2^j), constrained only by h <= 2^j, 0 <= e <= h and
e = R (mod 2), the realisability step used against h = 5 and h = 13 in (E2) is unavailable. Three extra first steps
then exist.
- **(5, 0, 8), type (ii-a), p = 5.** The output is h' = 5 + 2^7 = 133, with 2^s = 16 * 5 + 27 * 16 = 2^9, j' = m + 4
  and e' = 8 if m = p, else 0.
- **(13, e, 16), type (i), q = 2, p = 8.** The output is h' = 13 + 2^11 = 2061, with j' = m + 6 and e' = 32 if m = p,
  else 0.
- **(13, 0, 32), type (ii-a), p = 7.** The output is again h' = 2061, with 2^s = 2^14, j' = m + 7 and e' = 32 if m = p.

None of these outputs chains:
- 133 and 2061 are not in {1, 5, 13};
- the (beta) relation would need e' = 2^(p'-1) h' - d_(p') 2^(j'+1) with e' in {8, 32}, and for p' = 3, 5, 7 the right
  side is negative:
  - 532 - 1024, 2128 - 5120 and 8512 - 13312 for h' = 133, j' = 9;
  - 8244 - 32768, 32976 - 163840 and 131904 - 425984 for h' = 2061, j' = 14.

So (E3) holds for abstract triples too.

**Scope (audit M2).** (E2)/(E3) concern the new dyadic port at the exit inside the attached gadget. Outside the gadget,
one-step re-landings at old vertices exist (`audit/a15_scope_witnesses.out`):
- **W1.** TT_2 rooted at its sink, (1, 1, 1), cone through TT_3 at p = 3 to H = 5. The old source of T has
  (e, R) = (0, 32).
- **W2.** (1, 0, 4) cone to 9, with (9, 2, 16) at an old vertex.
- **W3.** (33, 0, 128) cone to 289, with (289, 64, 512) at an old vertex.

With re-rooting allowed (any vertex, or a new source above N), the audit found no second consecutive re-landing. Its
search covered all 181 dyadic bases of order <= 8, TT_m with m <= 4, and outputs of order <= 15
(`audit/a12_reroot_chains.out`). Whether re-rooted dyadic cone chains exist is OPEN.

**Cross-checks.**
- `e_classify.py` sweeps every abstract state (h, e, 2^j) with h <= 2^j, 0 <= e <= h and e = R (mod 2) for small j. It
  confirms the classification and finds no two-step chain.
- Pre-registered P4 searched every realised dyadic base of order <= 8 with m <= 20. It found exactly the steps
  1 -> 9 (two ways) and 33 -> 289. The latter is type (beta) with p = 3: e = 4 = 4 * 33 - 2^7. P4 found no two-step
  chain.

## 6. Theorem F (fixed-port sheets are exactly critical; elementary reformulation of known decoders)

**Priority (audit minor 3).**
- The content is the exact one-core decoder: Session197 `path-polynomials/ONE_CORE_BINARY_INVERSE.md` section 3,
  THM-5096 and THM-5085 (1)-(2).
- The "critical" reading has a quantitative precursor: Session206's zero-root-change mass (ln 2) X ("one free block
  never suffices"), recorded in HYP-9480 (C)/(D).
- What is new is only the tiling formulation.

**Setting.** Fix a port u with h_u = c, an odd number. After a bridge letter B (THM-5085), Delta = E - O = -c. Call such
states the *post-bridge sheet* of c. On the sheet, h and R are odd: h is an H-value, and R = 2R_prev - c.

**Statement.** On the lattice L_c = {(h, R) : h, R odd}, define

    T_0(h, R) = (h + R - c, 2R - c)              (the letter B),
    T_j(h, R) = (h + 2^j R, 2^(j+1) R - c)       (S^j then B, j >= 1).

Then:
- each T_j maps L_c injectively onto {(h', R') in L_c : v_2(R' + c) = j + 1};
- so the images tile L_c minus the line R' = -c, which has no predecessor;
- every point with R' > 0 has exactly one predecessor, under T_j with j = v_2(R' + c) - 1, and that predecessor has
  R = (R' + c)/2^(j+1) > 0.

Forward images of points with R > 0 can have R' < 0, for instance 2R - c when c > 2R. Predecessors of realised points
need not be realised.

**Proof.**
- **j = 0.** R' + c = 2R with R odd. Conversely R = (R' + c)/2, and h = h' - (R' - c)/2. Here (R' - c)/2 is even
  because R' + c = 2 (mod 4) and c is odd, so h is odd.
- **j >= 1.** R' + c = 2^(j+1) R. Conversely R = (R' + c)/2^(j+1), and h = h' - (R' + c)/2, which is odd.

QED

**Remark (heuristic).**
- T_j has determinant 2^(j+1), and the sum over j of 2^(-(j+1)) is 1. In this sense the fixed-port dynamics of a single
  core is exactly critical.
- That is a structural reading of why finite fixed-port banks are expected to miss a positive proportion (HYP-9480 (C),
  HYP-9655, OPEN). Mass enters through root changes, which move to the sheet of the current value.
- No coverage statement is claimed here.

# Exact observer compression and closed six-vertex operations

This addendum continues [NOTE.md](NOTE.md). The operations below are literal strongly connected, fixed-port, induced-source extensions. Their multipliers are not the desired ternary affine branches, so they do not settle the spectrum problem. The scalar wrapper and slot-profile multiplication are inherited from the credited S162 source; the observer formulas below are derived directly.

## 1. A wrapper has a two-coordinate exact observer

Let Y be any nonempty tournament, P=H(Y), R=F_Y(2), and form

    G(Y)=Y→a→b→Y,

with a retained as the observer port. The observer query c beats only a, and the next query d beats only c. Then

    v(G(Y),a)=(P+R, R/2, P+R/2).                         (1)

In particular y=h−x, and the exact profile is recovered by P=h−2x, R=2x.

**Proof.** Only b beats vertices of Y in each of the three query graphs. A Hamilton path therefore visits Y in at most two nonempty runs. If C_2 counts ordered covers of Y by two nonempty paths, then R=2P+C_2. The legal words in the collapsed run symbol Y and the distinguished vertices are as follows.

No query, one run: Yab, abY, bYa.
No query, two runs: YabY.

One query, one run: Yabc, Ycab, abYc, bYca, cabY.
One query, two runs: YabYc, YcabY.

Two queries, one run: Yabdc, Ycabd, Ydcab, abYdc, bYadc, bYcad, bYdca, cabYd, dcabY.
Two queries, two runs: YabYdc, YcabYd, YdcabY.

Each one-run word contributes P and each two-run word C_2. Thus the three queried counts are

    3P+C_2=P+R,
    5P+2C_2=P+2R,
    9P+3C_2=3P+3R.

Taking the observer differences proves (1). R is even because F_Y(2)=2P+C_2 and C_2 is even. ∎

This proves exact observer equality for all gate assignments having the same number of active gates in the construction from NOTE.md. The distinct local apex observers of the three inactive gates do not prevent equality after placing each gate uniformly inside this wrapper.

## 2. Two six-vertex strong operations with commuting observer maps

For the gate B(x,y) from NOTE.md, t=y(1−x), append B(x,y)⊕T_1 to the internal tournament Y, leaving the old wrapper vertices a,b and every old edge unchanged. Call the new wrapped graph G(Y⊕B⊕T_1). This adds exactly six vertices and retains the same physical observer port a. It is strong because the wrapper cycle connects all internal vertices.

Multiplication of slot profiles gives the exact physical transformations

    t=0: (P,R) → (9P,168R),
    t=1: (P,R) → (11P,200R).                              (2)

In observer coordinates, on the plane y=h−x these are

    t=0: (h,x,y) → (9h+318x, 168x, 9h+150x),
    t=1: (h,x,y) → (11h+378x, 200x, 11h+178x).             (3)

The two transformations commute, because (2) is diagonal. All three settings with t=0 give exactly the same new observer for every eligible source. Every output is another source of the same wrapper form, so these operations repeat indefinitely without a state-return assumption or a seed replacement.

## 3. A strong initial source and an exact valuation law

Start with Y=T_2, so P=1, R=4. Then G(Y) has rows (6,4,8,3), port a=2, and observer (5,2,3); it is strong. After s operations, k active and r=s−k inactive, the common exact observer is

    P=9^r11^k,
    R=4·168^r200^k,
    v=(P+R,R/2,P+R/2).                                   (4)

There are 6s+4 vertices and exactly s+1 distinct observer states. The state with k active operations has binom(s,k)3^r physical interface assignments. All are literal induced extensions of the same four-vertex initial source, and the initial port remains fixed.

The first new observers are

    inactive: (681,336,345),
    active:   (811,400,411).

For all s≥0,

    v_3(gcd(h,x,y))=r,
    sum_assignments w^{v_3(gcd(h,x,y))}=(1+3w)^s.          (5)

Indeed if r≥1, h=3^r[3^r11^k+4·56^r200^k], whose bracket is a 3-adic unit; x and y are also divisible by3^r by (4). If r=0, h≡(−1)^s+4(−1)^s=5(−1)^s≠0 mod3. Thus the modular AND detector is compatible with bounded repeated fixed-port source extensions, at the cost of six new vertices per operation. This is a closed source mechanism for its two displayed transformations, not closure under 3h, 3h+2, and 3h+4.

`check.py` verifies the observer formula on the original gate-composition family and checks all four first-step extensions from the four-vertex initial source, including their literal a and aa query graphs. Its conclusions remain active with Python -O.

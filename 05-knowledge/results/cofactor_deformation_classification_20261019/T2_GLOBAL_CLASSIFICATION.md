# Global classification of the t=2 dyadic cofactor deformation

19 October 2026. **PROVED, independently audited, with exact algebra and executable controls.**

This is an arithmetic deformation of the dyadic cofactor map. It is **not** the aliquot map. The proof below concerns this explicitly defined map only.

For a positive even integer n, write n=cp with c a power of two at least 2 and p positive odd. Define

    G(n)=F_2(n)=(c-1)p+4c-1=n-p+4c-1.

All iterates in this note use this same rule with the exact dyadic core of the current integer.

## 1. Main result

Every G-orbit eventually reaches either a fixed point

    c(4c-1),  c=2,4,8,...,

or exactly one of the following six cycles:

    16 -> 78 -> 46 -> 30 -> 22 -> 18 -> 16,
    24 -> 52 -> 54 -> 34 -> 24,
    28 -> 36 -> 42 -> 28,
    32 -> 158 -> 86 -> 50 -> 32,
    40 -> 66 -> 40,
    48 -> 108 -> 96 -> 220 -> 180 -> 150 -> 82 -> 48.

In particular, these are all nontrivial cycles, and no positive even orbit diverges.

The main claim follows from a stronger short-block statement: except for fixed points and the set

    E={2,4,6,8,10,12,16,24,28,32,40,44,48},

there is an iterate G^j(n)<n with 1<=j<=6. E lists exceptional starts for descent, not all cycle members.

## 2. Reduction and notation

Put q=p+1. Then

    G(n)-n=4c-q.

If q>4c there is immediate descent. If q=4c, n is fixed, giving precisely the fixed points stated above. We henceforth assume q<4c, and write

    q=dk,  d=2^{v_2(q)}>=2,  k positive odd.

The independently checked two-step reduction in FINITE_ODD_PART_REDUCTION.md gives G^2(n)<n whenever k>=5. For completeness: growth and k>=5 force d<c; hence the first image has exact core d, and

    d(G^2(n)-n)=(c+d)(4(d-1)-p)+2d+p+1.

Since p>=5d-1 and p+1<4c, its right side is strictly less than

    -(c+d)(d+3)+2d+4c=c(1-d)-d^2-d<0.

Only k=1 and k=3 remain. This is a reduction to two unbounded exponent lanes, rather than an appeal to a finite census.

If d<c, set e=c/d, a power of two at least 2. Direct substitution gives

    n_0=ed(dk-1),
    n_1=d[e(dk+3)-k],
    n_2=e(d-1)(dk+3)+d(4-k)+k-1.                 (2.1)

The core of n_1 is exactly d: e is even and k is odd. Thus n_2 in (2.1) is a genuine second iterate. Later cores are proved individually below.

## 3. The k=3 lane

Growth 3d<4c forces d<=c.

### 3.1. d<c

Equation (2.1) specializes to

    n_0=e(3d^2-d),
    n_2=3e(d^2-1)+d+2.                          (3.1)

If d>=4 and e>=4, the exact core of n_2 is 2, because its first term is divisible by 4 and d+2 is 2 modulo 4. Therefore

    n_3=3e(d^2-1)/2+d/2+8,
    n_3-n_0=-e(3d^2-2d+3)/2+d/2+8<0.

For example, the last expression is at most -6d^2+9d/2+2<0 for d>=4.

If e=2 and d>=8, then

    n_2=6d^2+d-4

has exact core 4. Hence

    n_3=9d^2/2+3d/4+12,
    n_3-n_0=-3d^2/2+11d/4+12<0.

For d=2, (3.1) gives n_0=10e and n_2=9e+4. Thus e>=8 gives two-step descent; e=2 gives the explicit descent from 20 in Section 5; and e=4 gives n_0=40 in E. The remaining pair d=4,e=2 gives n_0=88, also treated in Section 5.

### 3.2. d=c

Here n_0=c(3c-1). For c>=32 the consecutive exact formulas are

    n_1=3c^2,
    n_2=7c^2-4,
    n_3=21c^2/4+12,
    n_4=63c^2/16+24,
    n_5=441c^2/128+52,
    n_6=1323c^2/512+54.

The cores at n_1,n_2,n_3,n_4,n_5 are respectively c^2,4,4,8,4. Their guards follow immediately by reducing after division by those powers of two: c>=4 suffices at n_2, c>=8 at n_3, c>=16 at n_4, and c>=32 at n_5. Therefore

    n_6-n_0=-213c^2/512+c+54<0  (c>=32).

The remaining cores c=2,4 give n_0=10,44 in E; c=8,16 give the explicit six-step descents from 184 and 752 in Section 5. This exhausts k=3.

## 4. The k=1 lane

Growth d<4c permits d<c, d=c, or d=2c, because both d and c are powers of two.

### 4.1. d=2c

Then

    n_0=2c^2-c,
    n_1=c(2c+1),
    n_2=2c^2+3c-2.

The first image has exact core c, and n_2 has exact core 2 for c>=4. Thus for c>=8,

    n_3=c^2+3c/2+6,
    n_3-n_0=-c^2+5c/2+6<0.

The remaining c=2,4 give n_0=6,28 in E.

### 4.2. d=c

For c>=4,

    n_0=c(c-1),
    n_1=c(c+2),
    n_2=c^2+19c/2-2.

The core of n_1 is exactly 2c; for c>=8 the core of n_2 is 2. Hence for c>=16,

    n_3=c^2/2+19c/4+6,
    n_3-n_0=-c^2/2+23c/4+6<0.

The remaining c=2,4 give n_0=2,12 in E; c=8 gives the explicit six-step descent from 56 in Section 5.

### 4.3. d<c and e<d

Set A=(d-1)(d+3)=d^2+2d-3, which is odd. Equations (2.1) give

    n_0=ed(d-1),
    n_2=eA+3d.                                  (4.1)

If e<d, the exact core of n_2 is e: d/e is even. It follows that

    n_3=e(d+1)^2-d^2+d+2-3d/e,
    n_3-n_0=e(3d+1)-d^2+d+2-3d/e.               (4.2)

If e<=d/4, then d>=8 and (4.2) is at most

    -d^2/4+5d/4+2<0.

The only other possibility is e=d/2. In that case

    n_3=d^3/2+3d/2-4.

For d>=16 this has exact core 4, so

    n_4=3d^3/8+9d/8+12,
    n_4-n_0=-d^3/8+d^2/2+9d/8+12<0.

The remaining d=4 gives n_0=24 in E, and d=8 gives the explicit descent from 224 in Section 5.

### 4.4. d<c and e=d

Here c=d^2 and n_0=d^2(d-1). Equation (4.1) gives

    n_2=d^2(d+2).

For d>=4 its exact core is 2d^2. Consequently

    n_3=d^3+10d^2-d/2-2.

For d>=8 the exact core of n_3 is 2, and therefore

    n_4=d^3/2+5d^2-d/4+6,
    n_4-n_0=-d^3/2+6d^2-d/4+6<0  (d>=16).

The remaining d=2,4 give n_0=4,48 in E, while d=8 gives the explicit six-step descent from 448 in Section 5.

### 4.5. d<c and e>d

Put h=e/d=c/d^2, a power of two at least 2. The exact core of n_2 in (4.1) is d, since

    n_2=d[hA+3]

and hA+3 is odd. Thus

    n_0=hd^2(d-1),
    n_3=h(d-1)^2(d+3)+7d-4.                     (4.3)

There are three large-parameter regions.

1. If h=2 and d>=4, n_3 has exact core 2. Consequently

       n_4-n_0=-d^3+3d^2-3d/2+8<0.

   Indeed -d^3+3d^2<=-d^2 for d>=4.

2. If h>=8 and d>=8, n_3 has exact core 4: the first summand in (4.3) is divisible by 8, whereas 7d-4 has core 4. Therefore

       n_4-n_0=-h(d-1)(d-3)^2/4+21d/4+12<0.

   As (d-3)^2>=25 and h>=8, the right side is at most -179d/4+62<0.

3. If h=4 and d>=16, rewriting (4.3) gives

       n_3=4d^3+4d^2-13d+8,

   whose exact core is 8. Therefore

       n_4-n_0=-d^3/2+15d^2/2-91d/8+38<0.

   The first two terms are at most -d^2/2 for d>=16.

The remaining possibilities are controlled as follows.

- d=8,h=4 gives the explicit descent from 1792 in Section 5.
- d=4,h=4,8,16,32 give the explicit descents from 192,384,768,1536 in Section 5.
- For d=4,h>=64, equation (4.3) is n_3=63h+24, of exact core 8. Hence n_4=441h/8+52, of exact core 4, and

      n_5=1323h/32+54,
      n_5-n_0=-213h/32+54<0.

- For d=2,h>=16, (4.3) is n_3=5h+10, of exact core 2. Thus

      n_4=5h/2+12,
      n_4-n_0=12-3h/2<0.

- For d=2,h=2,4,8 the starts are n_0=8,16,32, all in E.

Because h and d are powers of two, the listed regions exhaust e>d.

## 5. All finite descent checks

The following are genuine consecutive G-iterates, with the final entry strictly smaller than the first. No row has more than six steps.

    20 -> 30 -> 22 -> 18,
    56 -> 80 -> 138 -> 76 -> 72 -> 94 -> 54,
    88 -> 108 -> 96 -> 220 -> 180 -> 150 -> 82,
    184 -> 192 -> 444 -> 348 -> 276 -> 222 -> 118,
    192 -> 444 -> 348 -> 276 -> 222 -> 118,
    224 -> 344 -> 332 -> 264 -> 262 -> 138,
    384 -> 892 -> 684 -> 528 -> 558 -> 286,
    448 -> 696 -> 640 -> 1146 -> 580 -> 450 -> 232,
    752 -> 768 -> 1788 -> 1356 -> 1032 -> 934 -> 474,
    768 -> 1788 -> 1356 -> 1032 -> 934 -> 474,
    1536 -> 3580 -> 2700 -> 2040 -> 1816 -> 1620 -> 1230,
    1792 -> 2808 -> 2488 -> 2208 -> 2266 -> 1140.

## 6. Exhaustion of exceptions and global conclusion

The values 16,24,28,32,40,48 in E lie on the six cycles displayed in Section 1. The other seven exceptional starts enter them through

    2 -> 8 -> 38 -> 26 -> 20 -> 30,
    4 -> 18,
    6 -> 10 -> 12 -> 24,
    44 -> 48.

Together these chains include all E, with their endpoints already on a listed cycle.

For any starting value outside the fixed points and E, Sections 2–5 produce a strictly smaller positive even iterate within six steps. Repeat this descent operation until a fixed point or E is reached. An infinite sequence of strictly decreasing positive integers is impossible. This proves eventual entry into the stated fixed-point/cycle list. It also excludes any additional cycle: its least element could neither have a smaller iterate nor belong to an unlisted terminal object.

## 7. Controls and limitations

The accompanying verify_t2_global.py checks symbolic identities, every finite chain and cycle, branch guards and predicted descent steps on an extensive exponent grid, and the short-block claim on an initial interval. It includes guard-boundary controls showing why several tempting smaller thresholds are invalid. These finite checks audit the written proof; they do not establish its unbounded inequalities or exhaustion on their own.

The result closes t=2 only. For general t, the finite odd-coefficient reduction remains valid, but its t residual coefficient lanes retain unbounded exponents. This proof does not establish global periodicity for t>=3, and no aliquot convergence claim follows from the t=2 theorem. Historical priority for this auxiliary-map result has not been established.

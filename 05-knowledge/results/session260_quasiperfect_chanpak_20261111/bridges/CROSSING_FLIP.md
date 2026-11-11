# Every crossing-edge reversal in the two-run cyclic tournament

**Status: PROVED all a,b>=2.** The sole crossing edge whose reversal decreases H is the boundary edge A_a->B_1. Its decrease is 2^(a-1)+2^(b-1)-2. All other A_i->B_j reversals strictly increase H. This extends arithmetic/RESEARCH.md Proposition C from one edge to every crossing rank; no claim about other edge classes is made in this file.

## 1. Exact response

Let T(a,b)=C3(TT_a,TT_b,1), with ordered vertices A_1,...,A_a and B_1,...,B_b, and singleton z. All A vertices beat all B vertices, all B vertices beat z, and z beats all A vertices. Put

    A=2^(a-1), B=2^(b-1), d=a-i, e=j-1,
    X=2^(i-1) 3^d - 2^d,
    Y=2^(b-j) 3^e - 2^e.

Reverse only A_i->B_j. The number L of lost paths and number G of gained paths are

    L=2AB/2^(d+e) - [d=0]B/2^e - [e=0]A/2^d
      +2[d=e=0],                                           (1)

    G=2XY + 2^(d+e+1) - [j=b]2^d - [i=1]2^e.              (2)

Therefore H(new)-H(old)=G-L. Brackets denote indicator functions. All quantities are integers.

## 2. Lost paths

Every old Hamilton path is either z at one end of the increasing A,B word, or

    A_L B_L z A_R B_R,

where the vertices within each displayed block are increasing, A_R is nonempty, and B_L is nonempty. The partitions into L,R are otherwise arbitrary.

For the marked adjacency A_i B_j in the left segment, A_i must be the largest vertex of A_L and B_j the smallest of B_L. There are

    (2^(i-1)-[i=a]) 2^(b-j)

choices: all A above i must go right; when i=a exclude the choice leaving A_R empty. For the marked adjacency in the right segment there are

    2^(i-1)(2^(b-j)-[j=1])

choices, now excluding empty B_L when j=1. Finally both end-z paths use the marked adjacency exactly when i=a and j=1. Their sum is (1).

## 3. Gained paths: three segments, not a scalar insertion

Every gained path uses consecutive B_j,A_i. Contract this pair to a marked symbol w. An outside predecessor of w may be any A vertex, or a B vertex below j; an outside successor may be any B vertex, or an A vertex above i. Neither adjacency z,w nor w,z is allowed.

Deleting both w and z leaves three increasing segments. Write each as an A segment followed by a B segment, allowing empties. There are two separator orders.

### w before z

Write the path as

    A_0 B_0 w A_1 B_1 z A_2 B_2.

The conditions are:

- B_0 consists only of B vertices below j, and may be empty;
- A_1 consists only of A vertices above i, and may be empty;
- B_1 is nonempty;
- either A_2 is nonempty or the whole terminal segment A_2 B_2 is empty.

If A_2 is nonempty, the A assignments number X: the i-1 low A vertices choose slots0 or2, the d high A vertices choose any of slots0,1,2, and the assignments with slot2 empty number 2^d. Independently the B assignments with B_1 nonempty number Y: the e low B vertices choose any of three slots, the b-j high B vertices choose slots1 or2, and the assignments with slot1 empty number 2^e. This contributes XY.

If the whole terminal segment is empty, the A choices number 2^d and the B choices number 2^e-[j=b]. This contributes 2^(d+e)-[j=b]2^d.

### z before w

Write the path as

    A_0 B_0 z A_1 B_1 w A_2 B_2.

Now A_1 is nonempty, B_1 is restricted to low B vertices, A_2 is restricted to high A vertices, and either B_0 is nonempty or the whole initial segment is empty. The nonempty B_0 case again contributes XY. The empty initial-segment case contributes

    (2^d-[i=1])2^e.

Adding the two separator orders proves (2). The partitions reconstruct literal directed paths uniquely.

**Failure boundary of an attractive shortcut.** It is incorrect to count G only by inserting w into Hamilton paths of the residual tournament. A new path may use a low B vertex -> w -> high A vertex; deleting w exposes a forbidden B->A adjacency. The three-segment count keeps this additional path class. A preliminary insertion-only formula failed at a=2,b=3,i=1,j=2 (predicted change9 instead of the actual10); that discarded formula is not used here.

## 4. Global sign classification for crossing edges

At d=e=0, equations (1)-(2) reduce to

    L=2AB-A-B+2,
    G=2(A-1)(B-1)+2,
    G-L=2-A-B<0.

It remains to show G-L>0 when d+e>=1.

### Two or more steps from the boundary

For fixed a, the function

    X(d)=2^(a-1-d)3^d-2^d,       0<=d<=a-1,

is nondecreasing, since for 0<=d<=a-2,

    X(d+1)-X(d)=2^(a-2-d)3^d-2^d >=0.

Thus X>=A-1. Similarly Y>=B-1. Let s=d+e. Using 2^d+2^e<=2^s+1 gives

    G >= 2(A-1)(B-1)+2^(s+1)-2^d-2^e
      >= 2(A-1)(B-1)+2^s-1.

Since A,B>=2,

    2(A-1)(B-1) >= AB/2.

For s>=2, (1) has no positive indicator term and gives L<=2AB/2^s<=AB/2. Hence

    G-L >= 2^s-1 >0.

### Exactly one step from the boundary

For d=1,e=0, direct simplification yields

    G-L=2AB-(5/2)A-4B+8-[a=2].                             (3)

If a=2 then A=2 and (3) equals2. If a>=3 then A>=4 and B>=2, so

    G-L = 2B(A-2)-(5/2)A+8 >= (3/2)A >0.

The case d=0,e=1 is the symmetric statement with A,B exchanged. This proves the strict sign classification for all a,b>=2.

## 5. Consecutive blocks and scope

For b=a+1 the unique decreasing crossing has decrement

    2^(a-1)+2^a-2 = 3*2^(a-1)-2.

Thus the proposed minimum positive decrease over crossing reversals is proved, and there are no other decreasing crossing reversals. Combined with the separate singleton-switch formulas and the positive substitution argument for internal block arcs, this closes the all-edge minimum theorem sought by the parent session. The present file supplies only the crossing-edge component, not an independent derivation of those other components.

For a=b=2, the boundary flip gives11->9. For the consecutive family a>=2, its decrement is at least4, so no crossing flip supplies the missing decrement2. This is a property of this explicit source graph, not an obstruction to realizing the target with a different tournament.

## 6. Small exact controls

`crossing_check.py` compares (1)-(2) against literal subset Hamilton-path DP for every crossing edge of every T(a,b) with 2<=a,b<=5:196 cases. The receipt is `crossing_receipt.json`. The all-order sign theorem is proved above rather than extrapolated from this finite range.

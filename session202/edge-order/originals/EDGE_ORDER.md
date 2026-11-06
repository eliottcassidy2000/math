# A physical two-port family with arbitrarily reversed order bias

Status: proof candidate, not independently audited. This result supplies actual
tournaments and a loss ledger for two-vertex overlap. It proves no spectrum tail.

## 1. Construction and exact sectors

Let p,q>=1. Take vertices x,y,b and two disjoint transitive blocks P,Q of
orders p,q. Orient P -> b -> Q and P -> Q, put x -> y, put
x -> (P union Q) -> y, and put y -> b -> x. These prescriptions give every
edge exactly once. In the private vertex order P,b,Q the private tournament
is transitive.

Let f be the number of Hamilton paths putting x before y, and r the number
putting y before x. Then

    f = 2^p + 2^q + 1,
    r = 2^(p+q),
    H = (2^p+1)(2^q+1).

Indeed x has exactly one inneighbor, b, and y has exactly one outneighbor,
b. A reverse-order Hamilton path must therefore contain the consecutive
subpath y,b,x. Each vertex of P union Q independently goes into the
increasing segment before y or the increasing segment after x. This gives
exactly 2^(p+q) paths.

A forward-order path has one of three disjoint shapes.

* If x is not first, b immediately precedes x and y is last. Before b is
  any subset of P, in increasing order, possibly empty. The other P vertices
  and all Q vertices follow x, in increasing order. There are 2^p choices.
* If y is not last, x is first and b immediately follows y. After b is any
  subset of Q, possibly empty. All P vertices and the other Q vertices lie
  between x and y in increasing order. There are 2^q choices.
* Otherwise x is first, y is last, and the only path is x,P,b,Q,y.

The first two conditions cannot both hold in a forward-order path, since
they would force y,b,x. The last shape exists because p,q>=1. This proves
all three exact formulas.

There is exactly one Hamilton cycle (counted up to cyclic rotation):
x,P,Q,y,b,x. In any cycle x must have predecessor b and y successor b,
forcing the cycle segment y,b,x; the remaining vertices then form their
unique transitive path. In particular these are actual strong C1 graphs.

At p=q=m, f/(f+r) -> 0 and r=(f-1)^2/4. Thus a universal positive lower
bound on the fraction of Hamilton paths following an oriented boundary edge
is false, even within C1. This is an obstruction to a *sector-majority*
shortcut, not an example disproving any specified gluing inequality.

## 2. Arbitrary private modules and the retained product

Replace P,Q by arbitrary nonempty tournaments, keeping the same cross-edge
orientations. Write h_P=H(P), and let c_P count ordered spanning two-path
covers of P allowing either path to be empty. Define h_Q,c_Q similarly.
The same three cases and reverse decomposition give

    f = c_P h_Q + h_P c_Q + h_P h_Q,
    r = c_P c_Q,
    H = (c_P+h_P)(c_Q+h_Q).

For a transitive p-vertex module, h_P=1 and c_P=2^p. The general formula
therefore specializes to Section 1. An arbitrary module need not preserve
the unique-cycle assertion, which was proved only for transitive modules.

This explains the limitation of the family for eventual coverage: the two
boundary sectors become arbitrarily unbalanced, but the scalar count remains
a product of two inherited cyclic-blowup counts. No additive closure follows.

## 3. Finite controls and exact scope

The companion verifier reconstructs the tournaments, counts Hamilton paths,
boundary sectors and Hamilton cycles independently by subset dynamic
programming, and checks the formulas on stated finite inputs. It also checks
the arbitrary-module formula against all labeled tournaments of orders at
most three, paired with each other. These checks validate the implementation;
the all-order statements rest on the case decompositions above.

The scripts are scratch research until an independent audit promotes them.
No numerical evidence here is a theorem of eventual odd-spectrum occupancy.

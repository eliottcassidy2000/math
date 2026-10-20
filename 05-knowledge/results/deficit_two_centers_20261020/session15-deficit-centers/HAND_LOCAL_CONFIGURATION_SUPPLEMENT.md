# A six-case hand replacement for the local configuration census

Status: new optional proof supplement; independent audit requested. It is not needed for the already approved complete D=z^4 proof.

Call a configuration a four-cycle exit if it consists of a directed four-cycle on S and a vertex v outside S beaten by at least three vertices of S. Call its converse a four-cycle entry. The prefix argument applies to an exit; for a self-converse ambient tournament, any entry is carried to an exit by an anti-automorphism.

Claim. Every self-converse strong tournament of order at least six contains a four-cycle exit. Consequently the connectivity-at-most-four reduction can be proved without enumerating all six-vertex tournaments.

## 1. Strong five-vertex cores

Moon's theorem supplies a strong induced five-vertex tournament U. We first show that either U already contains an exit or entry, or U is regular C5 or the nonregular tournament Q5=a->C3->b->a.

A strong four-vertex tournament is unique up to isomorphism. Label it 0->1->2->3->0 with remaining arcs 2->0 and3->1. A directed four-cycle in U supplies such a four-set; call the fifth vertex4.

If vertex4 beats at most one or at least three vertices of the four-set, there is immediately an exit or entry. The remaining possibilities are its six choices of two outneighbors. Four yield the following explicit exits:

    4 beats {0,1}: cycle (1,2,3,4), terminal0, entering vertices2,3,4.
    4 beats {0,2}: cycle (1,4,2,3), terminal0, entering vertices2,3,4.
    4 beats {0,3}: cycle (1,2,4,3), terminal0, entering vertices2,3,4.
    4 beats {1,2}: cycle (0,4,2,3), terminal1, entering vertices0,3,4.

If4 beats {2,3}, all five scores are2, giving the regular C5. If4 beats {1,3}, then a=2,b=1 and the directed triangle on {0,3,4} give exactly Q5=a->C3->b->a. These are the only remaining cases.

## 2. A vertex outside regular C5

Let v be any outside vertex. It either beats at least three C5 vertices or is beaten by at least three, since C5 has odd order five. Every four-vertex induced subtournament of regular C5 is strong. Choose a four-subset containing three vertices on that majority side. It gives an entry or exit. If it is an entry, apply an ambient anti-automorphism to get an exit.

## 3. A vertex outside Q5

Write Q5=a->C3->b->a and let v be outside. For every pair x,y of triangle vertices, {a,b,x,y} is a strong four-set. Unless v has exactly two inneighbors and two outneighbors on each of these three sets, there is already an exit or entry.

Suppose all three counts are balanced. Let t_x indicate x->v and let p be the number of vertices among {a,b} that beat v. For all distinct triangle vertices x,y,

    p+t_x+t_y=2.

Subtracting pairs shows all t_x are equal. There are just two cases.

- C3->v and v->{a,b}. Choose x->y in C3. Then (a,x,y,v) is a directed four-cycle, and x,y,v all beat b. This is an exit with terminal b.
- {a,b}->v and v->C3. Choose x->y in C3. Then (a,v,x,b) is a directed four-cycle, and a,v,x all beat y. This is an exit with terminal y.

Thus an exit or entry always exists; self-converseness converts any entry to an exit. This proves the claim.

## 4. What has and has not changed

At order n>=8 and connectivity at least five, delete the four vertices of the exit cycle. The remainder is strong and has a Hamilton path starting at its terminal; three cycle rotations supply the half-supported clique and force Delta>=4. The all-order cutoff and complete residual census remain exactly as previously approved.

This supplement replaces a local brute-force dependency with four explicit certificates, two recognizable five-vertex exceptions and their elementary extensions. The original independent six-vertex scans remain useful hostile validation. It does not remove the complete finite residual presentation census at the seven remaining orders.

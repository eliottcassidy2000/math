# Independent audit: uniform H4 equality for type (4)(2)

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
Root independently accepts [h4_mixed42](planar_jc48_sep08_h4_mixed42.md): four permutations, each consisting of one 4-cycle and one disjoint transposition, satisfying the declared marked H4 relations must all coincide. The ambient finite set is unbounded in size. Every such equal action is intransitive. Its application to the actual mixed-cusp covering uses the separately proved marked meridian and generation supplier; no general JC(2) conclusion follows.

I read every step of the proof and its standalone implementation. Both normal and optimized executions pass 343,032 gates and reproduce the frozen output byte for byte. An independent support-first census with union-find joint orbits reconstructs the complete 83,160-pair universe and every coupled typed-intersection/joint-orbit table. It does not import the producer implementation.

## Pair reduction and exact inventory

Fix sigma=(1234)(56) on twelve labels. Two permutations of the declared type move at most twelve labels in total. Simultaneous conjugacy can therefore put the first one at sigma and pad the pair's support union with fixed points. This proves the complete pair reduction in any ambient degree. It does not bound the support union of four generators by twelve.

There are binomial(12,4)*6*binomial(8,2)=83,160 distinct partners. The minimum-letter cycle convention removes precisely rotational duplicates and retains all six orientations of each four-set. The unique cycle lengths distinguish the 4-block and 2-block, so no unordered equal-block gauge is needed.

The exact ordinary table has 609 partners, with coupled rows

    (2,0,0,1) -> (3,6): 480
    (2,0,0,2) -> (2,6): 60
    (2,2,2,0) -> (6): 4
    (4,0,0,1) -> (3,4): 60
    (4,0,0,2) -> (2,4): 5.

The left side records the four typed intersections, and the right side the nontrivial joint-orbit sizes. In particular equal 4-blocks in an ordinary pair do not imply equality of its permutations; the five orientations in the final row remain in the proof.

The exact fifth table has 161 partners:

    (2,0,0,2) -> (2,6): 120
    (3,0,0,2) -> (2,5): 24
    (3,1,1,1) -> (6): 16
    (4,0,0,2) -> (2,4): 1.

Here the last row really is the identical partner. Every fifth pair either shares its transposition or has the full cross-mixing row (3,1,1,1). This dichotomy is stronger than an untyped total overlap and is needed later.

The commuting table has 212 partners with precisely the four displayed no-cross-mixing rows of the primary. The independent calculation verifies the joint table itself, not just separate marginal counts of matrices and orbit sizes. Literal cyclic orders are evaluated before recording the invariant; no sufficient matrix condition is promoted to an equivalence.

## Centralizers retain the missing global information

The commuting block assertion has a separate proof. A permutation commuting with a fixes its unique 4-orbit and 2-orbit setwise. On the four-orbit it is a power of that 4-cycle. The square power would introduce two transpositions, impossible for a permutation having only one. Thus its 4-block is either equal to that orbit or disjoint from it. Cross-length mixing is impossible, and the two transposition blocks are equal or disjoint. This remains valid with arbitrarily many additional fixed points.

For an ordinary pair a,b, all nontrivial joint-orbit sizes in the complete inventory are distinct. A permutation d centralizing both must therefore preserve each nontrivial joint orbit individually. It can still permute their common fixed points, and the proof makes no unsupported restriction there.

The centralizer of a transitive permutation group acts semiregularly: an element fixing one point fixes the entire orbit, by commuting with all elements sending that point to the others. Apply this to every power of d. Its effective cyclic order on a preserved orbit divides four, and all its cycles there have that effective order. On three points only identity is possible. On six points order four cannot occur; order two would require three transpositions, whereas d has only one in its complete cycle type. Hence d fixes every joint six-orbit pointwise. The possibility of a nontrivial four-orbit restriction is deliberately preserved.

The short typed obstruction when d fixes all of supp(b) is sound. Ordinary b,c gives at least two points of C4 in B4 and at least one point of C2 in supp(b). Consequently D4 can meet C4 in at most two points, and D2 can meet C2 in at most one. The first inequality forbids the cross-mixing fifth row with 4-block overlap three; the second forbids a shared transposition. All fifth rows are eliminated.

These facts are the proof of the four-generator, unbounded-degree step. The finite centralizer trials are supporting checks, not a substitute for it.

## Complete four-generator case split

If A4=B4, commutation of a,c plus ordinary b,c forces C4=A4: disjointness would violate the required 4-block overlap. Commutation of a,d and the fifth c,d relation similarly forces D4=A4. The fifth table then gives c=d exactly, including their transpositions and orientations. Now b commutes with c and has an ordinary relation with it; cancellation gives b=c. The same argument gives a=b. This covers the non-diagonal ordinary orientations as well as the equal pair.

Otherwise the ordinary joint orbit sizes are only (2,6), (3,6), or (6). In the last two cases, d fixes the whole support of b by the centralizer lemma, contradicting the typed obstruction above. In the (2,6) case, a and b share a transposition T and their four-block union is the fixed six-orbit.

Every remaining implication is paid separately. If C4=A4 then d fixes C4, contradicting its fifth overlap. Hence C4 is disjoint from both A4 and T by commutation with a. The cross-mixing ordinary b,c row (2,2,2,0) would put B2=T inside C4 and is impossible. Every other ordinary row makes C2 meet T. Commutation of a,c then forces C2=T. Finally, C4 has at least two points in B4, which d fixes. Thus C4 and D4 intersect in at most two points, excluding the fifth cross-mixing row. The other fifth alternative gives D2=C2=T.

Only now is the common two-point orbit removed. Its complement is invariant under all four permutations. Their restrictions are single 4-cycles and satisfy the original six relations literally. The proved [single-cycle H4 theorem](planar_jc48_sep08_h4_single_cycle.md) gives equality, and restoring T gives equality of the full permutations. No cycle-extraction map on arbitrary odd-related pairs is invoked.

## Hostiles, positive controls and the geometric consumer

The two displayed six-label hostiles are correct. The ordinary pair has cross-intersection row (2,2,2,0), while its isolated transpositions fail the ordinary relation. The fifth pair has row (3,1,1,1), while its isolated transpositions fail the fifth relation. Six labels are minimal for the declared single permutation type, so these are minimal ambient witnesses to premature extraction. They are not alleged full H4 quadruples.

The equal quadruple is a real positive control and remains intransitive because it already has distinct nontrivial orbits of sizes four and two. The conclusion is equality, not that each generator is identity. The separately proved [actual mixed-cusp supplier](planar_jc48_sep08_mixed_cusp_braid.md) provides the simultaneous marked relations and actual generation when this theorem is used for the curve. A transitive covering action cannot have the equal tuple. No equality of the full complement group with a presentation, full-retention hypothesis or Coxeter quotient is required.

## Independent census and reproduction

I enumerated a six-element moved support first, then its four-block, then all six cyclic orders. This differs from the producer's four-set-first enumeration and independently yields binomial(12,6)*90 partners. Odd words are evaluated by direct nested letter action. Joint orbits are computed by union-find over both generators' edges, independently of the producer's forward-orbit breadth-first calculation. All coupled tables above match exactly.

The following standalone calculation reproduces that independent path:

```python
from itertools import combinations, permutations
from collections import Counter
import json
from hashlib import sha256
D=12
sig=list(range(D))
for cyc in [(0,1,2,3),(4,5)]:
 for j,v in enumerate(cyc):sig[v]=cyc[(j+1)%len(cyc)]
S4={0,1,2,3};S2={4,5}

def relation(t,word):
 for start in range(D):
  left=right=start
  for pos in reversed(range(word)):
   left=(sig if pos%2==0 else t)[left]
   right=(t if pos%2==0 else sig)[right]
  if left!=right:return False
 return True

def orbit_sizes(t):
 # Union-find over all generator edges is independent of the producer BFS.
 par=list(range(D))
 def find(x):
  while par[x]!=x:x=par[x]
  return x
 for i in range(D):
  for j in (sig[i],t[i]):
   a,b=find(i),find(j)
   if a!=b:par[a]=b
 counts=Counter(find(i) for i in range(D))
 return tuple(sorted(v for v in counts.values() if v>1))

hist=[Counter(),Counter(),Counter()];total=0;seen=set()
for support in combinations(range(D),6):
 for block in combinations(support,4):
  two=tuple(x for x in support if x not in block)
  for tail in permutations(block[1:]):
   cycle=(block[0],)+tail;t=list(range(D))
   for j,v in enumerate(cycle):t[v]=cycle[(j+1)%4]
   t[two[0]]=two[1];t[two[1]]=two[0]
   pt=tuple(t)
   if pt in seen:raise RuntimeError('Duplicate')
   seen.add(pt);total+=1
   cell=(len(S4&set(block)),len(S4&set(two)),len(S2&set(block)),len(S2&set(two)))
   flags=[relation(t,3),relation(t,5),all(sig[t[i]]==t[sig[i]] for i in range(D))]
   if any(flags):
    sizes=orbit_sizes(t)
    for h,flag in zip(hist,flags):
     if flag:h[(cell,sizes)]+=1
expected=[
{((2,0,0,1),(3,6)):480,((2,0,0,2),(2,6)):60,((2,2,2,0),(6,)):4,((4,0,0,1),(3,4)):60,((4,0,0,2),(2,4)):5},
{((2,0,0,2),(2,6)):120,((3,0,0,2),(2,5)):24,((3,1,1,1),(6,)):16,((4,0,0,2),(2,4)):1},
{((0,0,0,0),(2,2,4,4)):90,((0,0,0,2),(2,4,4)):90,((4,0,0,0),(2,2,4)):30,((4,0,0,2),(2,4)):2}]
if total!=83160 or any(dict(a)!=b for a,b in zip(hist,expected)):raise RuntimeError((total,hist))
record=[sorted([(list(c),list(o),n) for (c,o),n in h.items()]) for h in hist]
print('Independent support-first joint typed/orbit census: 83160 PASS')
print('Admitted totals:',[sum(h.values()) for h in hist])
print('Joint record sha256:',sha256(json.dumps(record,separators=(',',':')).encode()).hexdigest())
```

It reports 83,160 distinct partners, admitted totals 609,161,212 and joint-record SHA256 `a46705a4d04250dbe036cde9681e1af6994c297b002af53a4b4ab7fb9f14bbda`.

I read the entire producer's source. Its 129,108 literal common-centralizer trials yield 1,892 admitted edges, 720 fixed-six checks and 720 common-transposition propagation cases. It also checks every coupled matrix/orbit row, both forms of the odd relations, and the declared hostile and positive controls. Explicit exceptions preserve the gates under Python optimization.

Reproduce the producer from the worktree root:

    python3 04-computation/planar_jc48_sep08_h4_mixed42.py
    python3 -O 04-computation/planar_jc48_sep08_h4_mixed42.py

Both independently executed modes pass 343,032 gates and equal the frozen 1,263 output bytes.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 9443 | `987f121ef43da7a4ea549c2d4cbd98d7d1c5206a10c0f2504744a2f52b28c9e3` |
| Output and both replays | 1263 | `499481bc6e285848612dff6e35cb00252c45ebb9dbbd435e45b160f180c3e963` |
| Primary before promotion | 15333 | `53f29c452a525087bf1d605f7b596a76b771d6ab629a868390cd04f2088d8689` |

No correction was needed. This audit accepts promotion with the declared cycle type, all six marked relations, arbitrary finite ambient size and exact geometric supplier requirements. Other mixed cycle types and the unrestricted planar Jacobian conjecture remain open.

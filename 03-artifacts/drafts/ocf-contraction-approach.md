# OCF via Contraction: adj(i,j) = H(G_{ij})

## Setup

Tournament T on n vertices with arc i->j. Contract i,j into single vertex w.

**Digraph G = G_{ij}:**
- Vertices: V\{i,j} union {w}, so n-1 vertices
- Arcs between v,u in V\{i,j}: same as T
- v->w iff T[v][i]=1 (v beats i)
- w->v iff T[j][v]=1 (j beats v)

**Digraph G' = G'_{ji}:**
- Same non-{i,j} arcs as T
- v->w iff T[v][j]=1 (v beats j)
- w->v iff T[i][v]=1 (i beats v)

**Key fact:** adj(i,j) = H(G), adj'(j,i) = H(G').

## Arc classification by s_x

For v in V\{i,j}, s_v = 1 - T[v][i] - T[j][v]:

| s_v | T[v][i] | T[j][v] | G arcs v<->w | G' arcs v<->w |
|-----|---------|---------|-------------|--------------|
| -1  | 1       | 1       | v->w, w->v (BOTH) | none |
|  0  | 1       | 0       | v->w only   | v->w only    |
|  0  | 0       | 1       | w->v only   | w->v only    |
| +1  | 0       | 0       | none        | v->w, w->v (BOTH) |

**Observation:** G and G' are IDENTICAL except:
- s=-1 vertices: bidirectional in G, disconnected in G'
- s=+1 vertices: disconnected in G, bidirectional in G'

## Proof framework

1. Paths in G where w's neighbors are all s=0 vertices = paths in G' with
   same property. These cancel in H(G) - H(G').

2. DeltaH = A - B where:
   - A = #{G-paths where w is adjacent to at least one s=-1 vertex}
   - B = #{G'-paths where w is adjacent to at least one s=+1 vertex}

3. A = swap involution's #U_T (unmatched T-paths)
   B = swap involution's #U_T' (unmatched T'-paths)

4. A bidirectional arc v<->w in G corresponds to a 3-cycle {i,j,v} in T
   (since v->i, i->j, j->v). So the s=-1 vertices ARE the third vertices
   of 3-cycles through i->j.

## Connection to cycles

In G, an odd L-cycle through i->j in T becomes an (L-1)-cycle through w.

- 3-cycle {i,j,v} in T -> 2-cycle (bidirectional arc v<->w) in G
- 5-cycle {i,j,a,b,c} in T -> 4-cycle {w,a,b,c} in G
- (2k+1)-cycle in T -> 2k-cycle in G

The bidirectional arcs in G encode EXACTLY the 3-cycles through {i,j}.
The 4-cycles through w in G encode the 5-cycles through {i,j} in T.

## The proof challenge

Show: H(G) - H(G') = DeltaI(Omega(T), 2) for all tournaments T.

Since G and G' only differ in arcs involving w:
- Relate H(G) to sub-tournament structure via bidirectional arcs
- The bidirectional arcs to s=-1 vertices give "flexibility" (choice of direction)
  that should correspond to independence polynomial terms

## n=4 proof via contraction

G is a digraph on 3 vertices {w, a, b}.
If s_a = s_b = 0: G = G', so H(G) = H(G'), delta = 0 = -2(0+0). OK.
If s_a = -1, s_b = 0: G has v_a<->w bidirectional, G' has no arc a<->w.
  H(G) has extra paths using either a->w or w->a.
  H(G') has no path with a adjacent to w.
  Delta = H(G) - H(G') = 2 = -2(-1) OK. (Verified for specific cases.)

This approach may generalize using induction on the number of s != 0 vertices.

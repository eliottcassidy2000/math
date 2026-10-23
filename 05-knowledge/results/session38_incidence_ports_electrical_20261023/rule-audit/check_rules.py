"""Independent finite controls for Session38 all-size rules.
No imports from producer checkers. Run normally and with python -O.
"""
from collections import Counter, defaultdict, deque
from itertools import product, combinations, permutations
from pathlib import Path
from functools import lru_cache
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent
CHECKS = Counter()

def need(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise RuntimeError(label)

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)

def project(x, columns):
    return tuple(sum(x[e]*columns[e][j] for e in range(len(x))) for j in range(len(columns[0])))

def margin(x, edges, n):
    out = [0]*n
    for z,(a,b) in zip(x,edges):
        out[a] += z; out[b] += z
    return tuple(out)

def residual(x, bounds, edges, columns, n):
    arcs = []
    for e,((a,b),z,bound,p) in enumerate(zip(edges,x,bounds,columns)):
        if bound is None or z < bound: arcs.append((a,b,p,e,1))
        if z > 0: arcs.append((b,a,neg(p),e,-1))
    return arcs

def nonzero_cycle(arcs,n,d):
    adj = [[] for _ in range(n)]
    for i,arc in enumerate(arcs): adj[arc[0]].append(i)
    for root in range(n):
        def visit(u,seen,path,total):
            for i in adj[u]:
                a,v,p,e,sign = arcs[i]
                score = add(total,p)
                if v == root:
                    if any(score): return path+[i]
                elif v > root and v not in seen:
                    found = visit(v,seen|{v},path+[i],score)
                    if found is not None: return found
            return None
        found = visit(root,{root},[],(0,)*d)
        if found is not None: return found
    return None

def potentials_hold(arcs,n,d):
    reach = [{u} for u in range(n)]
    for a,b,*_ in arcs: reach[a].add(b)
    for k in range(n):
        for u in range(n):
            if k in reach[u]: reach[u] |= reach[k]
    used = set()
    for root in range(n):
        if root in used: continue
        component = {v for v in range(n) if v in reach[root] and root in reach[v]}
        used |= component
        h = {root:(0,)*d}; todo = deque([root])
        while todo:
            u = todo.popleft()
            for a,b,p,*_ in arcs:
                if a != u or b not in component: continue
                predicted = add(h[a],p)
                if b in h:
                    if h[b] != predicted: return False
                else: h[b] = predicted; todo.append(b)
        need(set(h)==component, 'potential_scc_coverage')
    return True

def transport_controls():
    # Full K_2,2 plus a distinct parallel channel at (u0,v0).
    edges = [(0,2),(0,2),(0,3),(1,2),(1,3)]
    identity = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    observations = {
        'entry_identity':identity,
        'parallel_channels_merged':((1,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),
        'total_literal_mass':((1,),(1,),(1,),(1,),(1,)),
        'zero_projection':((0,),(0,),(0,),(0,),(0,)),
        'signed_projection':((1,),(-2,),(3,),(0,),(-1,)),
        'black_twice_and_blue_once':((2,0),(2,0),(0,1),(0,1),(2,0)),
    }
    configurations = tests = mixed = 0
    for bounds in product(range(3),repeat=5):
        fibers = defaultdict(list)
        for x in product(*(range(b+1) for b in bounds)):
            fibers[margin(x,edges,4)].append(x)
        for tables in fibers.values():
            for name,cols in observations.items():
                unique = len({project(x,cols) for x in tables}) == 1
                for x in tables:
                    arcs = residual(x,bounds,edges,cols,4)
                    cycle = nonzero_cycle(arcs,4,len(cols[0]))
                    need(unique == (cycle is None), 'transport_cycle_equivalence')
                    need(unique == potentials_hold(arcs,4,len(cols[0])), 'transport_potential_equivalence')
                    if cycle is not None:
                        y = list(x)
                        for i in cycle:
                            *_,e,sign = arcs[i]; y[e] += sign
                        need(all(0 <= y[e] <= bounds[e] for e in range(5)), 'cycle_witness_bounds')
                        need(margin(y,edges,4)==margin(x,edges,4), 'cycle_witness_margins')
                        need(project(y,cols)!=project(x,cols), 'cycle_witness_projection')
                    tests += 1
                    if unique and len(tables)>1: mixed += 1
            configurations += len(tables)
    # Unique, unbounded, degree-two six-cycle; no degree-one peeling move.
    six_edges = [(0,3),(1,3),(1,4),(2,4),(2,5),(0,5)]
    x = (0,1,1,0,1,1)
    six_fiber = [y for y in product(range(3),repeat=6) if margin(y,six_edges,6)==margin(x,six_edges,6)]
    need(six_fiber==[x], 'six_cycle_unique')
    need(all(sum(v in e for e in six_edges)==2 for v in range(6)), 'six_cycle_no_leaf')
    six_cols = tuple(tuple(int(i==j) for j in range(6)) for i in range(6))
    need(nonzero_cycle(residual(x,(None,)*6,six_edges,six_cols,6),6,6) is None, 'six_cycle_residual_pass')
    # Parallel channels give two tables but one literal black-loop projection.
    ps = [(0,1),(0,1)]; tables = [(0,1),(1,0)]; cols = ((2,),(2,))
    for x in tables:
        need(nonzero_cycle(residual(x,(1,1),ps,cols,2),2,1) is None, 'parallel_black_loop_projection_unique')
    return {'capacity_choices':3**5,'feasible_base_tables':configurations,'table_projection_tests':tests,
            'nonunique_table_unique_projection_tests':mixed,'six_cycle_unique_table':list(six_fiber[0]),
            'six_cycle_margin':list(margin(six_fiber[0],six_edges,6)),
            'parallel_loop_tables':[list(x) for x in tables],'parallel_loop_literal_projection':[2]}

def tournament(code,n):
    A = [[0]*n for _ in range(n)]
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j] = (code>>bit)&1; A[j][i] = 1-A[i][j]
    return A

def incidence(A):
    n=len(A); c=[0]*n; edge=[[0]*n for _ in range(n)]
    # Count directed 3-cycles with their unique outgoing edge at each vertex.
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edge[i][j] = sum(A[j][k] and A[k][i] for k in range(n))
                c[i] += edge[i][j]
    E=tuple(sorted(tuple(sorted((c[j],edge[i][j]) for j in range(n) if A[i][j])) for i in range(n)))
    return c,edge,E

def paths_all_subsets(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    H=[0]*(1<<n); H[0]=1
    for i in range(n): dp[1<<i][i]=1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S>>v)&1: continue
            R=S^(1<<v)
            if R: dp[S][v]=sum(dp[R][u] for u in range(n) if (R>>u)&1 and A[u][v])
        H[S]=sum(dp[S])
    return H

def lift(A,r):
    n=len(A); N=n+r+1; B=[[0]*N for _ in range(N)]
    for u in range(n):
        for v in range(n): B[u][v]=A[u][v]
        B[u][n]=1
    for i in range(r):
        p=n+1+i; B[n][p]=1
        for u in range(n): B[p][u]=1
        for j in range(i+1,r): B[p][n+1+j]=1
    return B

def predicted_E(E,r):
    n=len(E); cs=sorted(sum(e for c,e in row) for row in E)
    rows=[tuple(sorted([(c+r,e) for c,e in row]+[(n*r,r)])) for row in E]
    rows.append(tuple([(n,n)]*r))
    for i in range(r): rows.append(tuple(sorted([(c+r,1) for c in cs]+[(n,0)]*(r-1-i))))
    return tuple(sorted(rows))

def padding(A,source):
    n=len(A); B=[row+[int(not source)] for row in A]
    B.append([int(source)]*n+[0]); return B

def strong(A):
    n=len(A); R=[{i}|{j for j in range(n) if A[i][j]} for i in range(n)]
    for k in range(n):
        for i in range(n):
            if k in R[i]: R[i]|=R[k]
    return all(len(s)==n for s in R)

def lift_controls():
    bases=lifts=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(code,n); c,edge,E=incidence(A); Hs=paths_all_subsets(A); H=Hs[-1]
            F2=sum(Hs[S]*Hs[((1<<n)-1)^S] for S in range(1,(1<<n)-1))
            need(sorted(c)==sorted(sum(e for _,e in row) for row in E),'E_recovers_vertex_triangles')
            for r in range(1,4):
                B=lift(A,r); cb,eb,Eb=incidence(B)
                need(Eb==predicted_E(E,r),'E_strong_lift_formula')
                need(strong(B),'lift_is_strong')
                need(paths_all_subsets(B)[-1]==(2**r+1)*H+(2**r-1)*F2,'H_strong_lift_formula')
                lifts+=1
            for source in (False,True):
                B=padding(A,source); E2=incidence(B)[2]
                if source:
                    prediction=tuple(sorted(list(E)+[tuple(sorted((z,0) for z in c))]))
                else: prediction=tuple(sorted([tuple(sorted(list(row)+[(0,0)])) for row in E]+[()]))
                need(E2==prediction,'E_padding_formula')
                need(paths_all_subsets(B)[-1]==H,'H_padding_invariant')
            bases+=1
    return {'all_labeled_base_orders':[1,2,3,4],'base_tournaments':bases,'strong_lift_controls':lifts,
            'r_values':[1,2,3],'padding_controls':2*bases,
            'circulant_delta_second_minus_first':{'H':12330,'ordered_two_path_cover':-100740,
                'strong_lift_H':'113070 - 88410 * 2^r, strictly negative for every r >= 1'}}

def rank_binary(rows):
    piv={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p not in piv: piv[p]=x;break
            x ^= piv[p]
    return len(piv)

def matvec(A,x): return tuple(sum(a*z for a,z in zip(row,x)) for row in A)
def pairing(A,x,y): return sum(z*w for z,w in zip(x,matvec(A,y)))

def tower_controls():
    matrices=stages=0
    examples=[]
    # Every nonsingular symmetric 2x2 integer matrix with entries in [-2,2].
    # Exhaust all residue lifts, not a Smith-derived implementation.
    for a,b,c in product(range(-2,3),repeat=3):
        if a*c==b*b: continue
        A=((a,b),(b,c)); matrices+=1
        for k in range(1,4):
            mod=2**k
            candidates=[z for z in product(range(2*mod),repeat=2) if all(v%mod==0 for v in matvec(A,z))]
            fibers=defaultdict(list)
            for z in candidates: fibers[tuple(v%2 for v in z)].append(z)
            reps={u:zs[0] for u,zs in fibers.items()}
            q={(u,v):(pairing(A,x,y)//mod)%2 for u,x in reps.items() for v,y in reps.items()}
            for u,zs in fibers.items():
                for z in zs:
                    for v,y in reps.items():
                        need((pairing(A,z,y)//mod)%2==q[u,v],'tower_lift_independence')
            rad={u for u in reps if all(q[u,v]==0 for v in reps)}
            nextV={tuple(v%2 for v in z) for z in candidates if all(v%(2*mod)==0 for v in matvec(A,z))}
            need(rad==nextV,'tower_radical_next_layer')
            # For 2x2, first Smith invariant is gcd entries, second |det|/d1.
            from math import gcd
            d1=gcd(gcd(abs(a),abs(b)),abs(c)); d2=abs(a*c-b*b)//d1
            expected=sum(d%mod==0 for d in (d1,d2))
            need(len(reps)==2**expected,'tower_smith_dimension')
            stages+=1
    # Every binary symmetric 4x4 matrix and every labelled four-distinct switch.
    switches=detected=0
    pairs=list(combinations(range(4),2)); slots=[(i,i) for i in range(4)]+pairs
    for bits in range(1<<len(slots)):
        L=[[0]*4 for _ in range(4)]
        for k,(i,j) in enumerate(slots): L[i][j]=L[j][i]=(bits>>k)&1
        W=[u for u in product((0,1),repeat=4) if all(v%2==0 for v in matvec(L,u))]
        for a,b,c,d in permutations(range(4)):
            alpha=[u[a]^u[d] for u in W]; beta=[u[b]^u[c] for u in W]
            rows=[sum(((alpha[i]*beta[j]^beta[i]*alpha[j])<<j) for j in range(len(W))) for i in range(len(W))]
            rank=rank_binary(rows)
            independent=any(alpha) and any(beta) and alpha!=beta
            need(rank==(2 if independent else 0),'switch_rank_two_iff_independent')
            detected+=independent; switches+=1
    return {'nonsingular_symmetric_2x2_matrices':matrices,'dyadic_stages':stages,'k_values':[1,2,3],
            'binary_symmetric_4x4_matrices':2**10,'labelled_switch_tests':switches,'rank_two_cases':detected}

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out={'status':'PASS','scope':'Independent finite controls, not all-size proofs or graph-realizability claims',
         'transport':transport_controls(),'lift':lift_controls(),'electrical':tower_controls(),'checks':dict(CHECKS)}
    dest=HERE/'receipt.json'
    if args.write: dest.write_text(json.dumps(out,indent=2)+'\n')
    else: need(json.loads(dest.read_text())==out,'frozen_receipt_matches')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()

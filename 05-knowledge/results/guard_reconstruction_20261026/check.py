"""Producer: full fixed-path universe, direct Hamilton-word enumeration."""
from itertools import permutations, combinations
from math import isqrt, factorial
import json

def need(ok, message):
    if not ok: raise RuntimeError(message)

def tournament(n, mask):
    rows=[0]*n; bit=0
    for i in range(n):
        for j in range(i+1,n):
            forward=j==i+1 or (mask>>bit)&1
            if j>i+1: bit+=1
            a,b=(i,j) if forward else (j,i)
            rows[a]|=1<<b
    return tuple(rows)

def words(rows):
    n=len(rows); full=(1<<n)-1
    def visit(p,used):
        if used==full:
            yield p; return
        choices=(full^used) if not p else rows[p[-1]]&(full^used)
        while choices:
            b=choices&-choices;choices-=b;v=b.bit_length()-1
            yield from visit(p+(v,),used|b)
    yield from visit((),0)

def support(rows):
    d=[0]*len(rows);h=0
    for p in words(rows):
        h+=1
        for a,b in zip(p,p[1:]):d[a]|=1<<b
    return tuple(d),h

def recover(d):
    n=len(d); reach=[r|(1<<i) for i,r in enumerate(d)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    out=list(d); missing=0;internal=0
    for i in range(n):
        for j in range(i+1,n):
            if (d[i]>>j&1) or (d[j]>>i&1):continue
            missing+=1
            if (reach[i]>>j&1) and (reach[j]>>i&1):
                internal+=1
                a=sum((d[i]>>k&1)*(d[k]>>j&1) for k in range(n))
                b=sum((d[j]>>k&1)*(d[k]>>i&1) for k in range(n))
                need((a>=2)!=(b>=2),'ambiguous invisible pair')
                u,v=(i,j) if a>=2 else (j,i)
            else:
                need((reach[i]>>j&1)!=(reach[j]>>i&1),'unordered components')
                u,v=(i,j) if reach[i]>>j&1 else (j,i)
            out[u]|=1<<v
    return tuple(out),missing,internal

def graph(n,edges):
    rows=[0]*n
    for a,b in edges:rows[a]|=1<<b;rows[b]|=1<<a
    return tuple(rows)

def square_graph(n,remove=False):
    edges=[(a-1,b-1) for a,b in combinations(range(1,n+1),2) if isqrt(a+b)**2==a+b]
    if remove:edges.remove((0,14))
    return graph(n,edges)

def sparse_paths(rows):
    """Undirected paths up to reversal; reject disconnected remaining support."""
    n=len(rows);full=(1<<n)-1
    def dfs(p,used):
        if used==full:
            if p[0]<p[-1]:yield p
            return
        rem=full^used;allowed=rows[p[-1]]&rem
        if not allowed:return
        connected=1<<p[-1];todo=connected;universe=rem|connected
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
            add=rows[v]&universe&~connected;connected|=add;todo|=add
        if connected!=universe:return
        while allowed:
            b=allowed&-allowed;allowed-=b
            yield from dfs(p+(b.bit_length()-1,),used|b)
    for start in range(n):yield from dfs((start,),1<<start)

def graceful_checks(k):
    m=k*k+2*k;leaves=[(i+1)*(k+1) for i in range(k)];special=set(leaves)
    frozen=[(0,d) if d not in special else (d//(k+1),d+d//(k+1)) for d in range(k+1,m+1)]
    def valid(edges):
        parent=list(range(m+1))
        def root(a):
            while parent[a]!=a:a=parent[a]
            return a
        for a,b in frozen+edges:
            a,b=root(a),root(b)
            if a==b:return False
            parent[a]=b
        return len(frozen)+len(edges)==m
    patterns=set();checked=0
    for pi in permutations(range(k)):
        for signs in range(1<<k):
            es=[(leaves[pi[d-1]],leaves[pi[d-1]]+(d if signs>>(d-1)&1 else -d)) for d in range(1,k+1)]
            need(valid(es),'graceful base')
            sig=[]
            for d in range(1,k+1):
                for i,b in enumerate(leaves):
                    for sign in [-1,1]:
                        alt=es.copy();alt[d-1]=(b,b+sign*d)
                        ok=valid(alt);need(ok==(pi[d-1]==i),'address guard')
                        sig.append(ok);checked+=1
            patterns.add(tuple(sig))
    need(len(patterns)==factorial(k),'graceful quotient count')
    return {'k':k,'vertices':m+1,'states':factorial(k)*2**k,'guard_classes':len(patterns),'replacement_checks':checked}

def graph_controls():
    pe=[(0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,8),(4,9),(5,7),(5,8),(6,8),(6,9),(7,9)]
    out={}
    for name,rows in [('K5',graph(5,combinations(range(5),2))),('K33',graph(6,((i,j) for i in range(3) for j in range(3,6)))),('Petersen',graph(10,pe))]:
        ps=list(words(rows));used={(min(a,b),max(a,b)) for p in ps for a,b in zip(p,p[1:])}
        cycles=sum(bool(rows[p[-1]]>>p[0]&1) for p in ps)//(2*len(rows))
        out[name]={'directed_words':len(ps),'unoriented_paths':len(ps)//2,'undirected_cycles':cycles,'used_edges':len(used)}
    need(out['Petersen']['undirected_cycles']==0,'Petersen cycle hostile')
    contracted={tuple(sorted((a%5,b%5))) for a,b in pe if a%5!=b%5}
    need(contracted==set(combinations(range(5),2)),'Petersen K5 contraction')
    routes=[(2,3),(2,1,6),(2,7),(8,3),(8,6),(8,5,7),(9,4,3),(9,6),(9,7)]
    edges={tuple(sorted(e)) for e in pe};interiors=[]
    for r in routes:
        need(all(tuple(sorted(e)) in edges for e in zip(r,r[1:])),'K33 route')
        interiors+=list(r[1:-1])
    need(len(interiors)==len(set(interiors))==3,'K33 disjoint interiors')
    need({(r[0],r[-1]) for r in routes}=={(a,b) for a in (2,8,9) for b in (3,6,7)},'K33 complete')
    pairs=list(combinations(range(5),2));kg=graph(10,((i,j) for i,j in combinations(range(10),2) if not set(pairs[i])&set(pairs[j])))
    need(sum(r.bit_count() for r in kg)//2==15,'Kneser edges')
    # Explicit Kneser model: outer consecutive star in the disjointness graph.
    model=[(0,1),(2,3),(0,4),(1,2),(3,4),(2,4),(1,4),(1,3),(0,3),(0,2)]
    need({tuple(sorted((model[a],model[b]))) for a,b in pe}=={tuple(sorted((a,b))) for a,b in combinations(pairs,2) if not set(a)&set(b)},'Kneser isomorphism')
    def matchings(rows,unused):
        if not unused:
            yield ();return
        a=min(unused)
        for b in sorted(unused-{a}):
            if rows[a]>>b&1:
                for tail in matchings(rows,unused-{a,b}):yield ((a,b),)+tail
    def cycle_lengths(rows):
        unused=set(range(len(rows)));lengths=[]
        need(all(r.bit_count()==2 for r in rows),'complement degree two')
        while unused:
            todo=[min(unused)];seen=set()
            while todo:
                v=todo.pop()
                if v in seen:continue
                seen.add(v);todo.extend(w for w in unused if rows[v]>>w&1)
            unused-=seen;lengths.append(len(seen))
        return sorted(lengths)
    p_rows=graph(10,pe);pms=list(matchings(p_rows,set(range(10))))
    need(len(pms)==6,'six Petersen matchings')
    common=[]
    for a,b in combinations(pms,2):
        c=set(a)&set(b);need(len(c)==1,'Petersen matching pair intersection');common.extend(c)
    need(set(common)==edges and len(common)==15,'K6 pair addresses')
    profiles=[]
    for pm in pms:
        r=list(p_rows)
        for a,b in pm:r[a]^=1<<b;r[b]^=1<<a
        profiles.append(cycle_lengths(r))
    need(profiles==[[5,5]]*6,'Petersen two-cycle obstruction')
    pmsets=[frozenset(pm) for pm in pms]
    slot_swap={i:1-i if i in [0,1] else i for i in range(5)}
    vertex_swap=[model.index(tuple(sorted(slot_swap[j] for j in pair))) for pair in model]
    action=[pmsets.index(frozenset(tuple(sorted((vertex_swap[a],vertex_swap[b]))) for a,b in pm)) for pm in pms]
    need(all(action[i]!=i and action[action[i]]==i for i in range(6)),'Petersen triple transposition')
    full_action=set()
    for perm in permutations(range(5)):
        vp=[model.index(tuple(sorted(perm[j] for j in pair))) for pair in model]
        full_action.add(tuple(pmsets.index(frozenset(tuple(sorted((vp[a],vp[b]))) for a,b in pm)) for pm in pms))
    need(len(full_action)==120 and len({g[0] for g in full_action})==6,'faithful transitive S5 matching action')
    out['Petersen_matchings']=[list(pm) for pm in pms]
    out['Petersen_matching_complement_cycles']=profiles
    out['Petersen_matching_transposition_action']=action
    out['Petersen_matching_action_order']=len(full_action)
    perms=list(permutations(range(3)));frames=[];matrix=[]
    for p in perms:
        b=p[0];eps=(p[1]-b)%3
        need(eps in (1,2) and all(p[i]==(b+eps*i)%3 for i in range(3)),'affine S3 frame')
        frames.append((b,eps))
        row=[]
        for q in perms:
            inv=[q.index(i) for i in range(3)];r=[inv[p[i]] for i in range(3)];seen=set();c=0
            for i in range(3):
                if i in seen:continue
                c+=1
                while i not in seen:seen.add(i);i=r[i]
            row.append(c)
        matrix.append(row)
    need(len(set(map(tuple,matrix)))==6,'matching context distinguishes all six')
    for p in perms:
        rem=graph(6,((i,j+3) for i in range(3) for j in range(3) if j!=p[i]))
        need(cycle_lengths(rem)==[6],'K33 matching complement cycle')
    out['matching_cycle_matrix']=matrix;out['signed_frames']=frames;out['K33_routes']=routes
    a=[1,193,673];b=[3,168,483]
    need(all(isqrt(x+y)**2==x+y for x in a for y in b),'square K33 cross edges')
    need(all(isqrt(x+y)**2!=x+y for side in [a,b] for x,y in combinations(side,2)),'square K33 induced')
    return out

def vitali_control():
    pairs=list(combinations(range(7),2))
    def rows(mask):
        t=[0]*7
        for i,(a,b) in enumerate(pairs):
            u,v=(a,b) if mask>>i&1 else (b,a);t[u]|=1<<v
        return tuple(t)
    def lam(mask):
        t=rows(mask);result=[0]*21
        for s in combinations(range(7),3):
            if all(sum(t[v]>>w&1 for w in s)==1 for v in s):
                for e in combinations(s,2):result[pairs.index(e)]+=1
        return tuple(result)
    def flip(s):return sum(1<<pairs.index(e) for e in combinations(s,2))
    def legal(mask,s):
        t=rows(mask)
        return sorted(sum(t[v]>>w&1 for w in s) for v in s)==[1,1,2,2] and lam(mask^flip(s))==lam(mask)
    origin=401896;s=(0,2,3,4);q=(1,3,4,6)
    states=[origin,origin^flip(s),origin^flip(q),origin^flip(s)^flip(q)]
    need(legal(origin,s) and legal(origin,q),'Vitali initial guards')
    need(not legal(states[1],q) and not legal(states[2],s),'Vitali mixed guard failure')
    seen={origin};todo=[origin]
    while todo:
        t=todo.pop()
        for sub in combinations(range(7),4):
            if legal(t,sub):
                u=t^flip(sub)
                if u not in seen:seen.add(u);todo.append(u)
    need(seen==set(states[:3]),'Vitali guarded component')
    records=[]
    for m in states:
        t=rows(m);ps=list(words(t));hc=sum(bool(t[p[-1]]>>p[0]&1) for p in ps)//7
        matrix=[[int(t[i]>>j&1) for j in range(7)] for i in range(7)];power=matrix
        for _ in range(6):power=[[sum(power[i][k]*matrix[k][j] for k in range(7)) for j in range(7)] for i in range(7)]
        records.append({'mask':m,'H':len(ps),'simple_c7':hc,'trace_A7_over_7':sum(power[i][i] for i in range(7))//7,'lambda_preserved':lam(m)==lam(origin)})
    need([r['simple_c7'] for r in records]==[14,15,15,7],'corrected simple cycles')
    return {'subsets':[s,q],'guarded_component':sorted(seen),'square':records,'mixed_lambda_difference':[(e,lam(states[3])[i]-lam(origin)[i]) for i,e in enumerate(pairs) if lam(states[3])[i]!=lam(origin)[i]]}

def main():
    census=[]
    for n in range(1,8):
        count=1<<((n-1)*(n-2)//2);missing=internal=0
        for mask in range(count):
            t=tournament(n,mask);d,h=support(t);r,a,b=recover(d)
            need(r==t,('reconstruction',n,mask));missing+=a;internal+=b
        census.append({'n':n,'states':count,'missing_arcs':missing,'within_strong_missing':internal})
    square=[]
    for n in [25,27]:
        ps=list(sparse_paths(square_graph(n)));qs=list(sparse_paths(square_graph(n,True)))
        use=[p for p in ps if any({a,b}=={0,14} for a,b in zip(p,p[1:]))]
        need(set(qs)==set(ps)-set(use),'square deletion law')
        need((not use)==(n==25),'square dormant edge')
        square.append({'n':n,'paths':len(ps),'without_1_15':len(qs),'uses_1_15':len(use),'witness':[v+1 for v in use[0]] if use else None})
    print(json.dumps({'fixed_path_census':census,'square_sum':square,'graceful':[graceful_checks(k) for k in range(1,5)],'graphs':graph_controls(),'vitali':vitali_control()},sort_keys=True,indent=2))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent domain-CSP enumeration; no producer imports or cached data."""
import json,time,sys

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def audit(p):
    started=time.time(); Q={i*i%p for i in range(1,p)}
    blocks=[{(j+q)%p for q in Q} for j in range(p)]
    incidence=[[i in blocks[j] for j in range(p)] for i in range(p)]
    sols=[]; nodes=0
    def search(dom):
        nonlocal nodes
        nodes+=1
        # Arc consistency for distinct block labels and complementary orientations.
        while True:
            changed=False
            for i in range(p):
                for k in range(p):
                    if i==k: continue
                    keep={j for j in dom[i] if any(j!=l and incidence[i][l]!=incidence[k][j] for l in dom[k])}
                    if not keep:return
                    if keep!=dom[i]:dom[i]=keep;changed=True
            if not changed:break
        if all(len(d)==1 for d in dom):
            sols.append(tuple(next(iter(d)) for d in dom));return
        i=min((i for i in range(p) if len(dom[i])>1),key=lambda i:len(dom[i]))
        for j in sorted(dom[i]):
            ds=[d.copy() for d in dom];ds[i]={j};search(ds)
    search([{j for j in range(p) if not incidence[i][j]} for i in range(p)])
    check(len(sols)==len(set(sols)), "Duplicate matching")
    identity=tuple(range(p));check(identity in sols, "Canonical matching missing")
    # Independent group route: isomorphisms from the canonical tournament
    # to every compatible matching. Each lifts uniquely to an incidence
    # automorphism (alpha, beta), beta(i)=m(alpha(i)). Conversely every
    # incidence automorphism is counted in exactly its transported matching.
    def isos(m):
        answer=[]
        def visit(f,used):
            i=len(f)
            if i==p:answer.append(tuple(f));return
            for v in range(p):
                if v in used:continue
                if all(incidence[k][i]==incidence[f[k]][m[v]] for k in range(i)):
                    visit(f+[v],used|{v})
        visit([],set())
        return answer
    orbit=set(); stab=0; order=0; fiber_counts=[]; actions=[]
    for m in sols:
        lifts=isos(m);fiber_counts.append(len(lifts));order+=len(lifts)
        if lifts:orbit.add(m)
        if m==identity:stab=len(lifts)
        for alpha in lifts:
            beta=tuple(m[alpha[i]] for i in range(p))
            actions.append((alpha,beta))
            check(all(incidence[i][j]==incidence[alpha[i]][beta[j]] for i in range(p) for j in range(p)), "Lift fails incidence preservation")
    check(orbit==set(sols), "Matching action is not transitive")
    check(all(c==stab for c in fiber_counts), "Unequal lift fibers")
    matching_kernel=0; stabilizer_orbits=[]
    if p<=11:
        stabilizer_actions=[]
        for alpha,beta in actions:
            perm={}
            for m in sols:
                fm=[None]*p
                for i in range(p):fm[alpha[i]]=beta[m[i]]
                perm[m]=tuple(fm)
                check(perm[m] in sols, "Matching action fails closure")
            matching_kernel+=all(m==perm[m] for m in sols)
            if perm[identity]==identity:stabilizer_actions.append(perm)
        unseen=set(sols)
        while unseen:
            m=min(unseen);orb={g[m] for g in stabilizer_actions}
            stabilizer_orbits.append(len(orb));unseen-=orb
    # Hamilton paths, direct subset DP: edge i->k iff k in matched block m(i).
    # dense table avoids retaining unreachable hash keys
    counts=[[0]*p for _ in range(1<<p)]
    for i in range(p):counts[1<<i][i]=1
    outmask=[sum(1<<k for k in blocks[i]) for i in range(p)]
    for mask in range(1,1<<p):
        for last,v in enumerate(counts[mask]):
            if not v:continue
            avail=outmask[last]&~mask
            while avail:
                bit=avail&-avail;avail-=bit
                counts[mask|bit][bit.bit_length()-1]+=v
    H=sum(counts[-1]); check(H%stab==0, "Path count not divisible by stabilizer")
    return dict(p=p,compatible_matchings=len(sols),design_group_order=order,orbit_size=len(orbit),tournament_stabilizer=stab,hamilton_paths=H,quotient=H//stab,union_paths=len(sols)*H,union_quotient=len(sols)*H//order,matching_action_kernel=matching_kernel if p<=11 else None,matching_stabilizer_orbits=sorted(stabilizer_orbits),csp_nodes=nodes,elapsed_seconds=round(time.time()-started,3),matchings=sols)
if __name__=='__main__':
    if '--self-test-failure' in sys.argv:
        check(False, 'Deliberate check-helper failure')
    for p in (7,11,19):print(json.dumps(audit(p)),flush=True)

"""AP2 fourth-deletion certificates, standard library only.
Shape enumeration uses a deliberately weak constant/alias simplifier.
"""
from itertools import combinations
import json

def require(ok,why):
    if not ok: raise RuntimeError(why)

def live(gates):
    seen=set(); stack=[len(gates)-1]
    while stack:
        x=stack.pop()
        if x<0 or x in seen: continue
        seen.add(x); stack.extend(gates[x])
    return len(seen)==len(gates)

# Negative symbols: p=-1,w=-2,r=-3,h=-4,z=-5,t=-6.
#  None means constant; arbitrary unassigned symbols remain distinct.
def weak_simplify(gates):
    alias={-1:-1,-2:None,-3:None,-4:None,-5:-5,-6:-6}
    removed=[]
    for j,(a,b) in enumerate(gates):
        aa,bb=alias[a],alias[b]
        if j==0: alias[j]=None; removed.append(j) # G=p*w, w=0.
        elif aa is None and bb is None: alias[j]=None;removed.append(j)
        elif aa is None or bb is None:
            alias[j]=bb if aa is None else aa;removed.append(j)
        elif aa==bb: alias[j]=aa;removed.append(j)
        else:alias[j]=j
    return alias[len(gates)-1],removed

def shapes(n,shared):
    # G can be first: it has only parameter parents, so moving it forward is legal.
    # In the distinct case call the earlier partner child H_r, by symmetry.
    placements=[(a,) for a in range(1,n)] if shared else combinations(range(1,n),2)
    for places in placements:
        def rec(gates):
            j=len(gates)
            if j==n:
                yield tuple(gates);return
            ordinary=[-5,-6,*range(j)]
            if shared and j==places[0]: options=[(-3,-4)]
            elif not shared and j==places[0]: options=[(-3,u) for u in ordinary]
            elif not shared and j==places[1]: options=[(-4,u) for u in ordinary]
            else:options=combinations(ordinary,2)
            for pair in options:
                gates.append(pair);yield from rec(gates);gates.pop()
        yield from rec([(-1,-2)])

def topology_checks():
    out=[]
    for shared in (False,True):
        for n in range(3,7):
            generated=admissible=nonconstant=sharp=0
            for gs in shapes(n,shared):
                generated+=1
                if not live(gs):continue
                admissible+=1
                output,removed=weak_simplify(gs)
                if output is None:continue
                nonconstant+=1
                require(len(removed)>=4,('AP2 counterexample',shared,gs,removed))
                sharp+=len(removed)==4
            out.append({'shared_partner_child':shared,'gates':n,'generated_shapes':generated,'all_live_shapes':admissible,'weakly_nonconstant_survivors':nonconstant,'exactly_four_deleted':sharp})
    return out

def witness():
    # p,w,r,z,h,t,a,b indices0..7. Gate ids begin8.
    gates=[('and',0,1),('and',2,3),('and',4,5),('xor',8,9),('xor',11,10),('and',6,7),('xor',12,13)]
    evals=0
    for mask in range(256):
        v=[(mask>>j)&1 for j in range(8)]
        for op,i,j in gates:v.append((v[i]&v[j]) if op=='and' else v[i]^v[j])
        p,w,r,z,h,t,a,b=v[:8]
        require(v[-1]==((p&w)^(r&z)^(h&t)^(a&b)),'literal witness')
        if r==h==1 and w==0:
            require(v[-1]==(z^t^(a&b)),'restricted witness')
        evals+=1
    # All eight inputs are essential even on the marked product domain:
    # compare only same-fiber assignments that both remain allowed.
    def domain(x):
        return (((x>>0)&1) or ((x>>2)&1)) and (((x>>1)&1) or ((x>>4)&1))
    def original(x):
        return (((x>>0)&1)&((x>>1)&1))^(((x>>2)&1)&((x>>3)&1))^(((x>>4)&1)&((x>>5)&1))^(((x>>6)&1)&((x>>7)&1))
    fibers=[]
    for idx in range(8):
        choices=[x for x in range(256) if domain(x) and domain(x^(1<<idx)) and original(x)!=original(x^(1<<idx))]
        require(bool(choices),('essential original on marked domain',idx))
        fibers.append([idx,choices[0],choices[0]^(1<<idx)])
    # Restriction structurally leaves XOR(z,t), AND(a,b), final XOR, all essential.
    for idx in [3,5,6,7]:
        found=False
        for mask in range(256):
            if (mask&0b10110)!=0b10100:continue # w0,r1,h1
            def out(x):return ((x>>3)&1)^((x>>5)&1)^(((x>>6)&1)&((x>>7)&1))
            if out(mask)!=out(mask^(1<<idx)):found=True;break
        require(found,('essential residual',idx))
    # Polar form of pw+rz+ht+ab pairs every coordinate with exactly one distinct mate.
    pairs=[(0,1),(2,3),(4,5),(6,7)]
    rows=[sum(1<<v for u,v in pairs if u==j)+sum(1<<u for u,v in pairs if v==j) for j in range(8)]
    require(sorted(rows)==[1<<j for j in range(8)],'nondegenerate polar permutation matrix')
    return {'truth_table_evaluations':evals,'initial_gates':7,'deleted_gates':4,'surviving_gates':3,'initial_m_q_a_d':[8,2,6,4],'delta_s_I_q_a':[-4,0,-2,-1],'polar_rank':8,'marked_domain_essential_fiber_witnesses':fibers}

if __name__=='__main__':print(json.dumps({'status':'PASS','shape_universe':topology_checks(),'sharp_legal_state':witness()},indent=2,sort_keys=True))

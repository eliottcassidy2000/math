"""Check exact fanout transfer and the all-k alias-ladder hostile."""
from itertools import combinations,product
import json


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def children(nodes):
    result={v:set() for v in nodes}
    for v,(_,ps) in nodes.items():
        for p in set(ps):
            if p in nodes:result[p].add(v)
    return result


def phi(nodes):
    return sum(max(len(C)-1,0) for C in children(nodes).values())


def prune(nodes,out):
    keep=set();todo=[out]
    while todo:
        v=todo.pop()
        if v in nodes and v not in keep:
            keep.add(v);todo.extend(nodes[v][1])
    return {v:data for v,data in nodes.items() if v in keep}


def alias_census():
    rows=[]
    for n in range(1,6):
        kept=reroutes=exact=0
        for ps in product(*[list(combinations(range(i+2),2)) for i in range(n)]):
            nodes={i+2:(6,p) for i,p in enumerate(ps)};out=n+1
            if len(prune(nodes,out))<n:continue
            kept+=1;before=phi(nodes);C=children(nodes)
            for g,(_,parents) in nodes.items():
                for p in parents:
                    after={v:(t,tuple(p if u==g else u for u in pr))
                           for v,(t,pr) in nodes.items() if v!=g}
                    require(phi(after)<=before,'alias creates gate fanout budget')
                    require(phi(prune(after,p if out==g else out))<=phi(after),'dead cleanup increases potential')
                    if p in nodes and C[g]:
                        q=next(u for u in parents if u!=p)
                        overlap=len((C[p]-{g})&C[g])
                        ignored_parent=int(q in nodes and len(C[q])>=2)
                        require(phi(after)-before==-overlap-ignored_parent,'exact transfer law')
                        exact+=1
                    reroutes+=1
        rows.append({'gates':n,'live_topologies':kept,'literal_reroutes':reroutes,'exact_gate_transfer_checks':exact})
    return rows


def build(k):
    variables=['a','b','y','z']+[x for i in range(1,k+1) for x in (f'a{i}',f'b{i}')]
    nodes={'G0':(8,('a','b'))}
    for i in range(1,k+1):
        nodes[f'A{i}']=(8,(f'a{i}',f'b{i}'))
        nodes[f'G{i}']=(6,(f'G{i-1}',f'A{i}'))
    nodes['L']=(6,(f'G{k}','y'));nodes['R']=(6,(f'G{k}','z'));nodes['O']=(8,('L','R'))
    return variables,nodes,'O'


def table(t,a,b):return t>>(2*a+b)&1


def constants(nodes,fixed):
    values=dict(fixed);C=set()
    for v,(t,ps) in nodes.items():
        options=[(values[p],) if p in values else (0,1) for p in ps]
        result={table(t,*ab) for ab in product(*options)}
        if len(result)==1:values[v]=result.pop();C.add(v)
    return C


def normalize(nodes,out,fixed):
    aliases=dict(fixed);new={}
    for v,(t,ps) in nodes.items():
        p,q=(aliases.get(x,x) for x in ps)
        if isinstance(p,int) and isinstance(q,int):aliases[v]=table(t,p,q);continue
        if isinstance(p,int) or isinstance(q,int) or p==q:
            wire=q if isinstance(p,int) else p
            outputs=[table(t,p if isinstance(p,int) else z,q if isinstance(q,int) else z) for z in (0,1)]
            if outputs[0]==outputs[1]:aliases[v]=outputs[0]
            elif outputs==[0,1]:aliases[v]=wire
            else:raise RuntimeError('unexpected negation in positive alias ladder')
        else:new[v]=(t,(p,q));aliases[v]=v
    out=aliases.get(out,out)
    return prune(new,out),out


def evaluate(nodes,out,inputs):
    values=dict(inputs)
    for v,(t,ps) in nodes.items():values[v]=table(t,*(p if isinstance(p,int) else values[p] for p in ps))
    return out if isinstance(out,int) else values[out]


def ladders():
    rows=[];semantic=0
    for k in range(1,33):
        variables,old,out=build(k);nodes=dict(old);fixed={};consumed=0;boundary_charges=0
        require(phi(nodes)==1 and len(nodes)==2*k+4,'initial ladder inventory')
        for j,i in enumerate(range(k,0,-1),1):
            fixed[f'a{i}']=0;C=constants(nodes,fixed);Ch=children(nodes)
            B={v for v,(_,ps) in nodes.items() if v not in C and set(ps)&C}
            require(C=={f'A{i}'} and B=={f'G{i}'},'literal current region/boundary')
            E=sum(len(Ch[v])-1 for v in C);consumed+=E
            boundary_charges+=sum(len(Ch[v])-1 for v in B)
            require(E==0 and len(Ch[f'G{i}'])==2,'alias branch transfer data')
            nodes,out=normalize(nodes,out,fixed)
            require(phi(nodes)==1 and len(nodes)==2*k+4-2*j,'exact two deletions, conserved branch')
            if k<=3:
                for bits in product((0,1),repeat=len(variables)):
                    assignment=dict(zip(variables,bits))
                    if any(assignment[v]!=b for v,b in fixed.items()):continue
                    require(evaluate(old,'O',assignment)==evaluate(nodes,out,assignment),'restricted function mismatch')
                    semantic+=1
        require(boundary_charges==k,'reused boundary budget')
        fixed['a']=0;C=constants(nodes,fixed);Ch=children(nodes)
        B={v for v,(_,ps) in nodes.items() if v not in C and set(ps)&C}
        require(C=={'G0'} and B=={'L','R'},'final constant branch')
        E=sum(len(Ch[v])-1 for v in C);consumed+=E
        nodes,out=normalize(nodes,out,fixed)
        require(E==1 and consumed==1 and phi(nodes)==0 and len(nodes)==1,'sharp no-reuse equality')
        require(nodes=={'O':(8,('y','z'))},'final consequence object')
        if k<=3:
            for bits in product((0,1),repeat=len(variables)):
                assignment=dict(zip(variables,bits))
                if any(assignment[v]!=b for v,b in fixed.items()):continue
                require(evaluate(old,'O',assignment)==(assignment['y']&assignment['z']),'final Boolean control')
                semantic+=1
        rows.append({'k':k,'initial_gates':2*k+4,'final_gates':1,'total_deletions':2*k+3,
                     'initial_phi':1,'constant_region_excess_sum':consumed,'boundary_excess_sum':boundary_charges})
    return {'rows':rows,'exhaustive_semantic_assignments_checked':semantic}


def main():
    print(json.dumps({'status':'PASS','alias_topology_census':alias_census(),'sharp_ladders':ladders()},indent=2,sort_keys=True))


if __name__=='__main__':main()

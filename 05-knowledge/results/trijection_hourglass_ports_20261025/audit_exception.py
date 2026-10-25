"""Independent literal reconstruction of the uniform boundary certificate."""
import json


def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)


def run(k):
    n=2*k+3
    cycle=tuple(range(k+1))+tuple(range(k+2,n))+(k+1,)
    position={v:i for i,v in enumerate(cycle)}
    diffs={1}|set(range(k+2,2*k+2))
    b={(u,v) for u in range(n) for v in range(n)
       if (position[v]-position[u])%n in diffs}
    p=set(zip(range(n-1),range(1,n)))
    need(p<=b,('boundary P missing',k))
    t=p|{(v,u) for u,v in b-p}
    need({v for v in range(n) if (v,0) in t}==set(range(n-k,n)),('ports',k))
    need(all((n-1-v,n-1-u) in t for u,v in t),('source not blue',k))
    need((2,4) in t,('wrong chord',k))
    a=tuple(range(k+1));u=k+1;v=k+2;w=k+3;s=k+4;z=tuple(range(k+5,n))
    words=[a+(w,u,v,s)+z,a+(w,u)+z+(v,s),a+(u,)+z+(v,w,s)]
    arcs=t
    for i,q in enumerate(words):
        need(sorted(q)==list(range(n)),('word vertices',k,i))
        joins=set(zip(q,q[1:]))
        need(joins<=arcs,('illegal join',k,i))
        r=dict(zip(q,reversed(q)))
        need(any((r[y],r[x]) not in arcs for x,y in arcs),('blue word',k,i))
        arcs=joins|{(y,x) for x,y in arcs-joins}
    need(arcs==b,('wrong boundary target',k))
    return n


if __name__=='__main__':
    orders=[run(k) for k in range(4,51)]
    print(json.dumps({'status':'PASS','k_range':[4,50],
        'boundary_orders':[orders[0],orders[-1]],'literal_certificates':len(orders),
        'directed_black_moves':3*len(orders),'all_sources_reconstructed_from_cyclic_difference_set':True},indent=2))

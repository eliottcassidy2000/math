"""Exact two-anchor reflection fibres and regular-circle unit-graph recurrences."""
from collections import Counter
from itertools import product,combinations
import json
import sympy as s
from return_probe import HERE,need

def eligibility(d):
    k=3*d-1;edges=[]
    for i,j in combinations(range(1,k+1),2):
        if j-i==d:edges.append((i-1,j-1,0))
        elif i+j in (d,5*d):edges.append((i-1,j-1,1))
    return edges
def components(k,edges):
    seen=set();out=[]
    for start in range(k):
        if start in seen:continue
        nodes={start};queue=[start];seen.add(start)
        for a in queue:
            for i,j,sign in edges:
                b=j if a==i else i if a==j else None
                if b is not None and b not in seen:seen.add(b);nodes.add(b);queue.append(b)
        out.append(nodes)
    return out
def main():
    z=s.Symbol('z');E6=1+15*z**2+15*z**4+z**6;receipts=[]
    for d in range(1,13):
        k=3*d-1;edges=eligibility(d);parts=components(k,edges);sizes=sorted(map(len,parts))
        expected=[2]+([3] if d%2==0 else [])+[6]*((d-1)//2)
        need(sizes==sorted(expected),'exact K2/P3/C6 eligibility classification')
        poly=s.expand(z**(3*d+1)*(1+z)**(1+2*int(d%2==0))*E6**((d-1)//2))
        if d<=6:
            # Fix the first reflection sign to remove simultaneous reflection.
            graphs=Counter()
            for tail in range(1<<(k-1)):
                signs=tail<<1;word=0
                for e,(i,j,required) in enumerate(edges):
                    if (((signs>>i)^(signs>>j))&1)==required:word|=1<<e
                graphs[word]+=1
            need(set(graphs.values())=={2**(len(parts)-1)},'uniform metric-phase fibre after global reflection')
            counts=Counter(word.bit_count()+3*d+1 for word in graphs)
            need(sum(count*z**degree for degree,count in counts.items())==poly,'independent complete labelled unit-graph polynomial')
            finite={'phase_orbits':1<<(k-1),'distinct_graphs':len(graphs),'fibre':2**(len(parts)-1)}
        else:finite=None
        # Exact cyclotomic check of the geometric unit-step criterion.
        x=s.Symbol('x');N=6*d;Phi=s.Poly(s.cyclotomic_poly(N,x),x);unit_steps=[]
        for a in range(N):
            if s.rem(s.Poly(x**a+x**((N-a)%N)-1,x),Phi).is_zero:unit_steps.append(a)
        need(unit_steps==[d,5*d],'exact unit chord criterion in the cyclotomic field')
        if d>=3:
            prior=z**(3*(d-2)+1)*(1+z)**(1+2*int(d%2==0))*E6**((d-3)//2)
            need(s.expand(poly-z**6*E6*prior)==0,'all-size two-step polynomial recurrence')
        receipts.append({'d':d,'vertices_with_anchors':3*d+1,'eligibility_component_sizes':sizes,
             'unit_steps':unit_steps,'edge_polynomial':str(poly),'complete_phase_test':finite})
    # A separate four-point hostile: two equal individual anchored charts,
    # but induced unit completion changes a star to a triangle with a leaf.
    A=s.Matrix([0,0]);B=s.Matrix([1,0]);P=s.Matrix([0,1]);examples=[]
    for sign in (1,-1):
        Q=s.Matrix([s.sqrt(3)/2,s.Rational(sign,2)]);pts=[A,B,P,Q]
        unit=[(i,j) for i,j in combinations(range(4),2) if s.simplify((pts[i]-pts[j]).dot(pts[i]-pts[j])-1)==0]
        examples.append({'Q_sign':sign,'unit_edges':unit,'P_Q_squared':str(s.simplify((P-Q).dot(P-Q))),
                         'chromatic_polynomial':str(s.Symbol('q')*(s.Symbol('q')-1)**2*(s.Symbol('q')-2) if sign==1 else s.Symbol('q')*(s.Symbol('q')-1)**3)})
    need(examples[0]['P_Q_squared']=='1' and examples[1]['P_Q_squared']=='3','one actual cross-unit probe separates relative reflection')
    (HERE/'geometric_phases.json').write_text(json.dumps({'status':'FINITE_EXACT_GEOMETRIC_CERTIFICATE','regular_circle_family':receipts,
        'elementary_hostile':examples,'scope':'marked supplied point configurations; graphs are full induced unit-distance completions'},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: two-anchor phase fibres, complete graph counts through d6, exact unit geometry through d12',flush=True)
if __name__=='__main__':main()

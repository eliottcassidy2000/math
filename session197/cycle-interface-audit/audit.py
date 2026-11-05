"""Independent continuation rank, physical profiles, and transfer controls."""
import functools, importlib.util, itertools, json, pathlib
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('independent_rooted', ROOT.parent/'session197-lowcycle-rooted-audit'/'audit.py')
ra = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ra)

def cycles(rows):
    # Fix vertex 0, and recurse over the remaining vertices to close at 0.
    n = len(rows)
    if n < 3: return 0
    full = (1 << n)-1
    @functools.lru_cache(None)
    def walk(mask, last):
        if mask == full: return int(bool(rows[last]&1))
        return sum(walk(mask|(1<<v), v) for v in range(n)
                   if not mask&(1<<v) and rows[last]&(1<<v))
    return walk(1, 0)

def state(rows, r):
    p = ra.profiles(rows)[r]
    return p[:4] + (cycles(rows),)

def literal_HP(rows):
    return sum(all(rows[u]&(1<<v) for u,v in zip(p,p[1:]))
               for p in itertools.permutations(range(len(rows))))

def physical(A, a, G, b):
    # Gadget vertices first, unlike both source and prior audit constructors.
    m = len(G); n = len(A); rows = list(G) + [x<<m for x in A]
    for v in range(m):
        for u in range(n):
            if (u,v)==(a,b): rows[v] |= 1<<(m+u)
            else: rows[m+u] |= 1<<v
    return tuple(rows),b

def determinant(M):
    a = [[Fraction(x) for x in row] for row in M]; det=Fraction(1)
    for i in range(len(a)):
        k=next((k for k in range(i,len(a)) if a[k][i]), None)
        if k is None: return 0
        if k!=i: a[k],a[i]=a[i],a[k];det=-det
        pivot=a[i][i];det*=pivot
        for j in range(i+1,len(a)):
            q=a[j][i]/pivot
            for k in range(i,len(a)):a[j][k]-=q*a[i][k]
    assert det.denominator==1
    return int(det)

def main():
    seed=((2,4,1),0)
    inputs=[seed]+[physical(*seed,G,b) for G,b in [((0,),0),((0,1),0),((0,1,3),0),((8,1,3,6),1)]]
    expected=[(3,1,3,1,1),(5,2,6,1,1),(9,3,9,1,1),(15,3,15,1,1),(41,6,54,3,2)]
    matrix=[]; literal_checks=0
    for (A,a),want in zip(inputs,expected):
        x=state(A,a); assert x==want and len(ra.components(A))==1
        S,s=physical(A,a,(0,),0); E,e=physical(A,a,(0,1),0)
        row=[x[0],state(S,s)[0],state(E,e)[0],x[4],cycles(E)]
        assert row==[literal_HP(A),literal_HP(S),literal_HP(E),ra.literal_cycles(A),ra.literal_cycles(E)]
        literal_checks+=5;matrix.append(row)
    assert matrix==[[3,5,9,1,1],[5,9,17,1,1],[9,15,27,1,1],[15,27,45,1,1],[41,89,149,2,3]]
    assert determinant(matrix)==-24
    # Check factorization into actual states and continuation functionals.
    query_rows=[(1,0,0,0,0),(1,-1,1,0,0),(1,0,2,0,0),(0,0,0,0,1),(0,0,0,1,0)]
    assert [[sum(u*v for u,v in zip(x,q)) for q in query_rows] for x in expected]==matrix
    assert determinant(expected)!=0 and determinant(query_rows)!=0
    collisions=[]
    for A,a in [((32,17,11,19,37,14),3),((32,17,11,35,13,22),2)]:
        p=ra.profiles(A)[a];x=state(A,a)
        assert len(ra.components(A))==1
        assert (p[0],p[1],p[2],p[4],x[4])==(37,5,71,13,4)
        T,r=physical(A,a,(0,),0);y=state(T,r)
        assert y[:3]==(103,32,108) and y[4]==x[3]
        assert literal_HP(A)==37 and ra.literal_cycles(A)==4
        assert literal_HP(T)==103 and ra.literal_cycles(T)==x[3]
        literal_checks+=4
        collisions.append({'input':A,'root':a,'state':x,'deletion_H':p[4],'output':y})
    assert [c['output'][4] for c in collisions]==[9,8]
    # Every labelled rooted gadget through order four, against all five
    # reachable witness inputs. Non-strong and zero-start gadgets included.
    transfer_checks=0;zero_cycle=0;nonunit_d=0
    for n in range(1,5):
        for G in ra.labelled(n):
            pg=ra.profiles(G)
            for b in range(n):
                g,E,Q,t,d,B,L=pg[b]
                for (A,a),x in zip(inputs,expected):
                    h,e,R,s,C=x;T,r=physical(A,a,G,b);actual=state(T,r)
                    want=(g*h-t*e+L*R,E*h-int(n==1)*e,Q*h+B*R,d*s,E*s)
                    assert actual==want,(A,a,G,b,actual,want)
                    transfer_checks+=1;zero_cycle+=int(not actual[4]);nonunit_d+=int(d>1)
    out={'status':'PASS','input_states':expected,'continuation_matrix':matrix,
         'determinant':determinant(matrix),'state_determinant':determinant(expected),
         'query_determinant':determinant(query_rows),'literal_output_checks':literal_checks,
         'physical_transfer_checks':transfer_checks,'zero_cycle_controls':zero_cycle,
         'nonunit_deletion_controls':nonunit_d,'collisions':collisions}
    (ROOT/'AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()

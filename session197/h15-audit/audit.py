import functools,importlib.util,itertools,json,math,pathlib
from fractions import Fraction
ROOT=pathlib.Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
pa=load('prime_audit',ROOT.parent/'session197-prime-atoms-audit'/'audit.py')
ca=load('census_audit',ROOT.parent/'session197-lowcycle-audit'/'check.py')
A=(16,9,3,5,14);B=(24,17,3,6,12);C=(36,1,2,7,15,30);D=(44,1,2,6,15,30)
TRI=(4,1,2);FIVE=(8,1,3,6)

def join(blocks):
    n=sum(map(len,blocks));rows=[];off=0
    for block in blocks:
        later=sum(1<<v for v in range(off+len(block),n))
        rows.extend((r<<off)|later for r in block);off+=len(block)
    return tuple(rows)

def runs(rows,r):
    # Memoized word construction, not the source's permutation/gap enumeration.
    outside=[v for v in range(len(rows)) if v!=r];full=(1<<len(outside))-1
    values=[]
    for K in range(1,len(rows)+1):
        @functools.lru_cache(None)
        def finish(mask,last,k):
            if mask==full and k==K:return 1
            total=0
            for j,v in enumerate(outside):
                if not mask&(1<<j) and rows[last]&(1<<v):total+=finish(mask|(1<<j),v,k)
            if last!=r and k<K and rows[last]&(1<<r):total+=finish(mask,r,k+1)
            return total
        value=finish(0,r,1)
        value+=sum(finish(1<<j,v,0) for j,v in enumerate(outside))
        values.append(value)
    return values

def matching(rows,r):
    edges=[(u,v) for u in range(len(rows)) for v in range(len(rows)) if rows[r]&(1<<u) and rows[v]&(1<<r) and rows[u]&(1<<v)]
    best=0
    for mask in range(1<<len(edges)):
        es=[edges[j] for j in range(len(edges)) if mask&(1<<j)]
        if len({u for u,v in es})==len(es)==len({v for u,v in es}):best=max(best,len(es))
    return best

def clone(rows,r,m):
    labels=[r]*m+[v for v in range(len(rows)) if v!=r];out=[0]*len(labels)
    for i,j in itertools.combinations(range(len(labels)),2):
        if labels[i]==labels[j] or rows[labels[i]]&(1<<labels[j]):out[i]|=1<<j
        else:out[j]|=1<<i
    return tuple(out)

def main():
    expected={ca.canonical(T) for T in (A,B,C,D)};found=set();census=[];min6=10**9;mixed_checks=0
    for n in range(1,9):
        classes=[];mass=0;strong_mass=0
        for line in (ROOT.parent/'session197-lowcycle-audit'/f'classes{n}.txt').read_text().splitlines():
            key,h,c,st=map(int,line.split());rows=tuple((key>>(8*i))&255 for i in range(n))
            if n==6 and st:min6=min(min6,h)
            if st and 3<=n<=6:
                for mask in range(1,(1<<n)-1):
                    child=tuple(row|((1<<n) if not mask&(1<<v) else 0) for v,row in enumerate(rows))+(mask,)
                    hh=pa.data(child)[1][-1]
                    assert hh>=h+2 and hh%2==1 and len(pa.components(child))==1
                    mixed_checks+=1
            if h!=15:continue
            classes.append(rows);weight=math.factorial(n)//ca.autos(rows);mass+=weight
            if st:found.add(ca.canonical(rows));strong_mass+=weight
        if n>=5:
            assert len(classes)==(n-3)*(n-4)
            assert mass==math.factorial(n)*(5*n*n-17*n-32)//15
        else:assert not classes
        census.append({'n':n,'unlabelled_H15':len(classes),'labelled_H15':mass,'labelled_strong_H15':strong_mass})
    assert found==expected and min6==15
    assert [ca.autos(T) for T in (A,B,C,D)]==[3,5,1,1]
    # A separate published check of the strong cutoff, not needed by the
    # source's valid universal-growth proof.
    assert 15**3<5**6
    constructed=[]
    for n in range(5,13):
        graphs=[];mass=Fraction(0)
        for T,aut in zip((A,B,C,D),(3,5,1,1)):
            for p in range(max(0,n-len(T)+1)):
                q=n-len(T)-p;graphs.append(join([(0,)]*p+[T]+[(0,)]*q));mass+=Fraction(1,aut)
        for left,right in ((TRI,FIVE),(FIVE,TRI)):
            for p in range(max(0,n-6)):
                for q in range(n-7-p+1):
                    s=n-7-p-q;graphs.append(join([(0,)]*p+[left]+[(0,)]*q+[right]+[(0,)]*s));mass+=Fraction(1,3)
        assert len(graphs)==(n-3)*(n-4)
        assert len({tuple(tuple(pa.induced(T,c)) for c in pa.components(T)) for T in graphs})==len(graphs)
        assert math.factorial(n)*mass==math.factorial(n)*(5*n*n-17*n-32)//15
        if n<=9:
            for T in graphs:assert pa.data(T)[1][-1]==15
        constructed.append({'n':n,'classes':len(graphs),'labelled':int(math.factorial(n)*mass)})
    functions={(0,0,15)};details=[];actual=0
    for T in (A,B,C,D):
        roots=[]
        for r in range(len(T)):
            a=runs(T,r);nu=matching(T,r)
            assert a[0]==15 and nu<=2 and max(i+1 for i,x in enumerate(a) if x)==nu+1 and not any(a[3:])
            coeff=(a[1]-3*a[2],a[2],15-2*a[1]+3*a[2]);functions.add(coeff)
            values=[]
            for m in range(1,7):
                h=pa.data(clone(T,r,m))[1][-1];assert h==coeff[0]*(1<<m)+coeff[1]*3**m+coeff[2];values.append(h);actual+=1
            roots.append({'root':r,'run_counts':a[:3],'matching':nu,'coefficients':coeff,'values':values})
        details.append({'rows':T,'roots':roots})
    expected_functions={(0,0,15),(9,0,-3),(5,0,5),(6,0,3),(11,0,-7),(5,2,-1),(7,2,-5),(4,6,-11)}
    assert functions==expected_functions
    small={}
    for T in (TRI,FIVE):
        vals=set()
        for r in range(len(T)):
            a=runs(T,r);assert matching(T,r)==1 and not any(a[2:]);vals.add((a[1],a[0]-2*a[1]))
        small[len(T)]=sorted(vals)
    assert small=={3:[(1,1)],4:[(2,1),(3,-1)]}
    for blocks in ((TRI,FIVE),(FIVE,TRI),((0,),A), (B,(0,))):
        T=join(blocks)
        for r in range(len(T)):
            vals=[pa.data(clone(T,r,m))[1][-1] for m in range(1,5)]
            assert any(vals==[a*(1<<m)+b*3**m+c for m in range(1,5)] for a,b,c in functions)
            actual+=4
    doubled=sorted({4*a+9*b+c for a,b,c in functions});assert doubled==[15,25,27,33,37,41,59]
    assert all(8*a+27*b+c>=45 for a,b,c in functions if a or b)
    result={'status':'PASS','independent_census':census,'padding_formula_controls':constructed,
            'strict_mixed_strong_extension_checks':mixed_checks,
            'core_automorphisms':[3,5,1,1],'clone_roots':details,'complete_functions':sorted(functions),
            'small_component_functions':small,'actual_clone_graph_checks':actual,'two_clone_values':doubled}
    (ROOT/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',actual,'actual clone checks')

if __name__=='__main__':main()

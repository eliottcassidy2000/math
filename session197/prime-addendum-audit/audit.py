"""Independent boundary data and sequential nonnegative component transfer."""
import functools,importlib.util,itertools,json,math,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('prime_audit',ROOT.parent/'session197-prime-atoms-audit'/'audit.py')
pa=importlib.util.module_from_spec(spec);spec.loader.exec_module(pa)

def converse(rows):return tuple(sum(1<<v for v in range(len(rows)) if rows[v]&(1<<u)) for u in range(len(rows)))

@functools.lru_cache(None)
def local(rows,mask):
    pp,H,B,start,end=pa.data(rows);qq=pa.data(converse(rows))[0];full=(1<<len(rows))-1
    O=[v for v in range(len(rows)) if mask&(1<<v)]
    return (H[-1],sum(H[S]*H[full^S] for S in range(full+1)),
            sum(pp[x][2] for x in O),sum(qq[x][2] for x in O),
            sum(end[-1][x] for x in O),sum(start[-1][x] for x in O),
            sum(B[y][x] for x in O for y in O if x!=y))

def advance(state,entry):
    d,U,f,a=state;h,q,u,v,A,b,F=entry
    return h*d,q*U+u*d,h*f+F*d+v*U+b*a,A*d

def profile(state):
    d,U,f,a=state;return d+U-f-a,d-a,d+U,d

def assemble(blocks,outs):
    n=sum(map(len,blocks));rows=[];rootout=0;offset=0
    for C,mask in zip(blocks,outs):
        future=sum(1<<v for v in range(offset+len(C),n))
        rows.extend((r<<offset)|future for r in C)
        rootout|=mask<<offset;offset+=len(C)
    return tuple(rows[i]|((1<<n) if not rootout&(1<<i) else 0) for i in range(n))+(rootout,)

def compositions(total,parts):
    if parts==1:yield (total,);return
    for first in range(total+1):
        for rest in compositions(total-first,parts-1):yield (first,)+rest

def classify(bank,d,target):
    k=(1+math.isqrt(1+8*target))//2;ell=max(1,target).bit_length()-1
    found={};examined=matching=0
    for singles in range(k+ell+1):
        # Distribute singleton blocks among gaps, rather than choose core
        # positions and then enumerate global incidence words as the source.
        for gaps in compositions(singles,len(bank)+1):
            blocks=[]
            for i,C in enumerate(bank):blocks.extend([(0,)]*gaps[i]);blocks.append(C)
            blocks.extend([(0,)]*gaps[-1])
            options=[[(mask,mask.bit_count(),local(C,mask)) for mask in range(1<<len(C))] for C in blocks]
            def visit(i,state,used,seen_out,unresolved_zero,outs):
                nonlocal examined,matching
                if i==len(blocks):
                    if not seen_out or unresolved_zero:return
                    examined+=1
                    if state[2]!=target:return
                    matching+=1;U=state[1]
                    if U not in found:
                        T=assemble(blocks,outs);z=profile(state)
                        assert pa.data(T)[0][-1]==z
                        found[U]={'rows':T,'root':len(T)-1,'profile':z}
                    return
                singleton=len(blocks[i])==1
                for mask,pop,entry in options[i]:
                    if used+pop>k:continue
                    if singleton and not mask and not seen_out:continue
                    pending=0 if mask else unresolved_zero+int(singleton)
                    visit(i+1,advance(state,entry),used+pop,seen_out or bool(mask),pending,outs+[mask])
            visit(0,(1,0,0,0),0,False,0,[])
    return examined,matching,found

def main():
    census=ROOT.parent/'session197-lowcycle-audit';types={h:[] for h in (3,5,7,9,11,13)}
    minimum6=10**9;strong6=0
    for n in range(1,9):
        for line in (census/f'classes{n}.txt').read_text().splitlines():
            key,H,C,st=map(int,line.split())
            if not st:continue
            rows=tuple((key>>(8*i))&255 for i in range(n))
            if n==6:
                minimum6=min(minimum6,H)
                # Automorphism order by degree-preserving maps.
                degrees={d:[v for v in range(n) if rows[v].bit_count()==d] for d in sorted({r.bit_count() for r in rows})}
                aut=0
                for images in itertools.product(*(itertools.permutations(c) for c in degrees.values())):
                    p=list(range(n))
                    for cell,image in zip(degrees.values(),images):
                        for u,v in zip(cell,image):p[u]=v
                    aut+=all(bool(rows[u]&(1<<v))==bool(rows[p[u]]&(1<<p[v])) for u in range(n) for v in range(n))
                strong6+=math.factorial(n)//aut
            if H in types:
                assert n<=5;types[H].append(rows)
    assert minimum6==15 and strong6==22320
    assert {h:len(v) for h,v in types.items()}=={3:1,5:1,7:0,9:2,11:1,13:1}
    # Independently, the published strong bound already gives this cutoff.
    assert 13**3<5**5
    roots=0
    for n in range(1,6):
        for rows in pa.source.labelled(n):
            rows=tuple(rows);actual=pa.data(rows)[0]
            for r,z in enumerate(actual):
                vs=[v for v in range(n) if v!=r];S=pa.induced(rows,vs);state=(1,0,0,0)
                for C in pa.components(S):
                    block=pa.induced(S,C);mask=sum(1<<i for i,v in enumerate(C) if rows[r]&(1<<vs[v]))
                    state=advance(state,local(block,mask))
                assert profile(state)==z;roots+=1
    report={'status':'PASS','independent_strong6_labelled_count':strong6,
            'independent_strong6_minimum':minimum6,'strong_core_types':types,
            'sequential_transfer_root_checks_through_n5':roots,'tables':{}}
    expected={(9,1):([20,24,32,44],18092,5),(9,-1):([26,27,34,35,41,45,99],77271,15),
              (11,1):([29,37,47],8493,3),(11,-1):([],8493,0),(13,1):([],8493,0),(13,-1):([],8493,0)}
    extension_checks=0
    for d in (9,11,13):
        banks=[(C,) for C in types[d]]
        if d==9:banks.append((types[3][0],types[3][0]))
        for beta in (1,-1):
            allfound={};examined=matches=0
            for bank in banks:
                a,b,found=classify(bank,d,d-beta);examined+=a;matches+=b;allfound.update(found)
            exp,ea,eb=expected[d,beta]
            assert sorted(allfound)==exp and examined==ea and matches==eb,(d,beta,examined,matches,sorted(allfound))
            report['tables'][f'{d}:{beta}']={'examined':examined,'matching':matches,'bases':sorted(allfound),'witnesses':allfound}
            for U,w in allfound.items():
                T=w['rows'];r=w['root'];h,e,R,dd=w['profile'];f=R+e-h-dd
                for m in (1,2):
                    n=len(T)
                    T=tuple(row|((1<<n) if i!=r else 0) for i,row in enumerate(T))+(1<<r,)
                    assert pa.data(T)[0][r]==((U<<m)+dd-f,dd,dd+(U<<m),dd)
                    extension_checks+=1
            print('independent table',d,beta,examined,matches,sorted(allfound),flush=True)
    allbases={(1,1):[1],(1,-1):[3,5],(3,1):[6],(3,-1):[9,15],(5,1):[11,12,16],(5,-1):[15,39],(7,1):[],(7,-1):[]}
    allbases.update({key:value[0] for key,value in expected.items()})
    allowed=[]
    for p in (3,5,7,11,13,17,19):
        for beta in (-1,1):
            cap=(2*p+2*beta-1)//3
            if any(p in allbases[d,beta] for d in range(1,cap+1,2)):allowed.append((p,beta))
    assert allowed==[(3,-1),(5,-1),(11,1)]
    assert all(23 not in allbases[d,-1] for d in range(1,15,2))
    report['prime_amplitudes_through19']=allowed
    report['independent_actual_K_extension_checks']=extension_checks
    assert extension_checks==28
    hostile=(52,1,2,7,14,30);assert pa.data(hostile)[0][5]==(19,1,73,9)
    report['power_of_two_hostile_profile']=[19,1,73,9]
    (ROOT/'AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',flush=True)

if __name__=='__main__':main()

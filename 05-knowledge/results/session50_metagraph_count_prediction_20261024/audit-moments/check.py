#!/usr/bin/env python3
"""Independent Session50 audit: source bytes, normalization, inference, hostile controls."""
import argparse,csv,hashlib,itertools as it,json,math,subprocess,sys,tempfile
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
from reconstruct import congruence,decode,mul,require,polynomial_bounds
from simple_certificate import certify
ROOT=Path(__file__).resolve().parent


def v2(x):
    if x==0:return math.inf
    x=abs(x);return (x&-x).bit_length()-1


def decode_matrix(code,n):
    rows=[0]*n;shift=n*(n-1)//2
    for a in range(n):
        for b in range(a+1,n):
            shift-=1
            if code>>shift&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
    return tuple(rows)


def transpose(rows):
    n=len(rows);mask=(1<<n)-1
    return tuple(mask^row^(1<<v) for v,row in enumerate(rows))


def path_count(rows):
    # Subset recurrence by the FIRST vertex, unlike inherited endpoint implementation.
    n=len(rows);table={1<<v:[int(v==w) for w in range(n)] for v in range(n)}
    for size in range(2,n+1):
        for vertices in it.combinations(range(n),size):
            mask=sum(1<<v for v in vertices);values=[0]*n
            for v in vertices:
                prior=table[mask^(1<<v)]
                values[v]=sum(prior[w] for w in vertices if rows[v]>>w&1)
            table[mask]=values
    return sum(table[(1<<n)-1])


def signatures(rows):
    degrees=[r.bit_count() for r in rows]
    return [(degrees[v],tuple(sorted(degrees[w] for w in range(len(rows)) if row>>w&1))) for v,row in enumerate(rows)]


def iso_count(rows,target):
    n=len(rows);left,right=signatures(rows),signatures(target)
    if sorted(left)!=sorted(right):return 0
    candidates=[[w for w in range(n) if left[v]==right[w]] for v in range(n)]
    order=sorted(range(n),key=lambda v:len(candidates[v]));image=[-1]*n
    def visit(depth,used):
        if depth==n:return 1
        v=order[depth];total=0
        for w in candidates[v]:
            if used>>w&1:continue
            if any(((rows[v]>>u)&1)!=((target[w]>>image[u])&1) for u in order[:depth]):continue
            image[v]=w;total+=visit(depth+1,used|(1<<w));image[v]=-1
        return total
    return visit(0,0)


def canonical(rows):
    n=len(rows);best=None
    # Canonical search only within degree-ordered labels is sufficient and independent.
    groups=defaultdict(list)
    for v,row in enumerate(rows):groups[row.bit_count()].append(v)
    menus=[list(it.permutations(groups[d])) for d in sorted(groups)]
    for blocks in it.product(*menus):
        order=sum(blocks,());code=0
        for i in range(n):
            for j in range(i+1,n):code=2*code+((rows[order[i]]>>order[j])&1)
        if best is None or code<best:best=code
    return best


def fixed_words(n):
    pairs=[(a,b) for a in range(n) for b in range(a+2,n)]
    for word in range(1<<len(pairs)):
        rows=[(1<<(v+1)) if v+1<n else 0 for v in range(n)]
        for j,(a,b) in enumerate(pairs):
            if word>>j&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        yield tuple(rows)


def small_enumeration():
    results=[]
    for n in range(1,7):
        hist=Counter();representatives={}
        for rows in fixed_words(n):
            code=canonical(rows);hist[code]+=1;representatives[code]=rows
        h_hist=Counter(hist.values());N=sum(hist.values());m=(n-1)*(n-2)//2
        require(N==1<<m,'small fixed word total')
        merged={}
        for code,h in hist.items():
            rows=representatives[code];H=path_count(rows);aut=iso_count(rows,rows)
            bruteH=sum(all(rows[a]>>b&1 for a,b in zip(order,order[1:])) for order in it.permutations(range(n)))
            require(H==bruteH,'independent small Hamiltonian count')
            require(H%2==aut%2==h%2==1 and H==aut*h,'small H/Aut odd fibre')
            other=canonical(transpose(rows));require(other in hist and hist[other]==h,'small converse symmetry')
            key=tuple(sorted({code,other}));merged[key]=sum(hist[c] for c in key)
        for key,mu in merged.items():require((1 if mu%2 else 2)==len(key),'small merged parity decoder')
        k=1
        while 2*k-k.bit_count()<m:k+=1
        moments=[sum(mu*(mu if mu%2 else mu//2)**r for mu in merged.values()) for r in range(k)]
        pred=decode({'n':n,'moments':moments})
        require(pred['prediction']==len(hist),'small finite-moment exact reconstruction')
        results.append({'n':n,'fixed_words':N,'ordinary_classes_verified_after_prediction':len(hist),'merged_nodes':len(merged),'h_histogram':dict(sorted(h_hist.items())),'moments':moments,'k':k,'modulus':pred['modulus'],'prediction':pred['prediction']})
    return results


def polynomial_controls():
    identity_checks=0
    for k in range(1,41):
        Q=[1]
        for j in range(k):Q=mul(Q,[-(2*j+1),1])
        require(Q[0]%2==1,'odd constant')
        e=2*k-k.bit_count()
        for h in range(1,258,2):
            value=sum(c*h**j for j,c in enumerate(Q))
            t=(h-1)//2
            require(value==(1<<k)*math.factorial(k)*(math.comb(t,k) if t>=k else 0),'falling-factorial identity')
            require(value%(1<<e)==0,'universal odd-grid divisibility')
            identity_checks+=1
        require(v2(math.prod(2*k+1-(2*j+1) for j in range(k)))==e,'exact Q valuation')
        # Exact Newton denominator valuation, including j=0.
        value=F(1)
        for j in range(k+1):
            if j:value*=F(-(2*j-1),2*j)
            require(v2(value.denominator)==2*j-j.bit_count(),'Newton denominator valuation')
    exhaustive=0
    for k in range(1,5):
        for c0 in [-1,1]:
            for tail in it.product(range(-2,3),repeat=k):
                coeff=(c0,)+tail
                val=min(v2(sum(c*(2*t+1)**i for i,c in enumerate(coeff))) for t in range(k+1))
                require(val<=2*k-k.bit_count(),'integer polynomial odd-constant optimality control')
                exhaustive+=1
    # Odd-mass multisets, not assumed to arise from tournaments: tests the full generality.
    artificial=0
    for masses in it.product(range(1,12,2),repeat=4):
        for k in range(1,7):
            moments=[sum(h**(r+1) for h in masses) for r in range(k)]
            _,_,M,residue=congruence(moments)
            require(len(masses)%M==residue,'artificial odd-mass congruence')
            artificial+=1
    return {'falling_factorial_checks':identity_checks,'degree_le4_integer_polynomial_sharpness_controls':exhaustive,'artificial_odd_mass_reconstructions':artificial,'Newton_denominator_checks_through_degree':40}


def atlas_audit():
    for pin in json.loads((ROOT/'SOURCE_PINS.json').read_text()):
        b=(ROOT/pin['snapshot']).read_bytes()
        require(len(b)==pin['bytes'] and hashlib.sha256(b).hexdigest()==pin['sha256'],'source byte pin')
    rows=list(csv.DictReader((ROOT/'sources/merged_nodes.tsv.txt').read_text().splitlines(),delimiter='\t'))
    explicit=json.loads((ROOT/'sources/explicit_matrices.json.txt').read_text())
    for slot in explicit['slots']:
        a=decode_matrix(slot['canon'],8)
        require([[int(row>>w&1) for w in range(8)] for row in a]==slot['canonical_matrix'],'frozen explicit matrix decoding')
    masses=[];actual_count=0;self_count=0;maxaut=0
    for index,row in enumerate(rows):
        require(int(row['merged_id'])==index,'sequential source IDs')
        a=decode_matrix(int(row['canon_a']),8);b=decode_matrix(int(row['canon_b']),8)
        H=path_count(a);aut=iso_count(a,a);maxaut=max(maxaut,aut)
        require(H>0 and aut>0 and H%aut==0 and H%2==aut%2==1,'all atlas odd H and Aut')
        h=H//aut
        self_converse=iso_count(a,transpose(a))>0
        require(iso_count(a,transpose(b))==aut,'literal converse pairing')
        require(self_converse==(row['class_a']==row['class_b']),'atlas self-converse flag independently recomputed')
        multiplicity=1 if self_converse else 2;mu=int(row['fibre'])
        require(mu==h*multiplicity,'all atlas H/Aut normalization')
        require((mu if mu%2 else mu//2)==h,'merged mass parity alone gives h')
        masses.append(mu);actual_count+=multiplicity;self_count+=self_converse
    # Decoder payload generated from merged numeric masses ONLY: no class labels/count/hmax.
    moments=[sum(mu*(mu if mu%2 else mu//2)**r for mu in masses) for r in range(12)]
    require(moments[0]==1<<21,'atlas complete fixed-word checksum')
    prediction=decode({'n':8,'moments':moments[:8]})
    short=certify(moments[:8])
    require(F(short['upper'])==F(496641366973194539,15625000000000),'short exact upper evaluation')
    require(F(short['upper'])<short['modulus'] and short['prediction']==prediction['prediction'],'short upper alone isolates positive residue')
    require(prediction['candidate_count']==1,'unique interval congruence representative')
    require(prediction['prediction']==actual_count,'held-back source count compared AFTER prediction')
    # Generic 12-moment theorem uses only [1,N], so bypass optional rational bounds.
    _,e,M,residue=congruence(moments)
    first=1+(residue-1)%M;last=1<<21
    require(M>=last and first==actual_count,'generic 12-moment reconstruction')
    low,up,certificate=polynomial_bounds(moments[:8])
    for h in [F(1),F(3,2),F(2),F(3),F(17),F(661),F(10**9)]:
        pl=sum(F(c)*h**j for j,c in enumerate(certificate['lower_polynomial']))
        pu=sum(F(c)*h**j for j,c in enumerate(certificate['upper_polynomial']))
        require(pl<=1/h<=pu,'global square-identity spot control')
    # Input leakage must not be accepted even if additional values are unused.
    rejected=[]
    for key,value in [('target',actual_count),('hmax',max(mu if mu%2 else mu//2 for mu in masses)),('atlas',rows)]:
        try:decode({'n':8,'moments':moments[:8],key:value})
        except ValueError:rejected.append(key)
        else:raise RuntimeError('unexpected side-channel input admitted')
    # All statistics survive arbitrary atlas order.
    require(moments==[sum(mu*(mu if mu%2 else mu//2)**r for mu in masses[::-1]) for r in range(12)],'row order invariance')
    # Separate-process actual replay with only decoder and 2-key data in temp directory.
    with tempfile.TemporaryDirectory(prefix='s50-blind-decoder-') as temp:
        p=Path(temp);(p/'reconstruct.py').write_bytes((ROOT/'reconstruct.py').read_bytes())
        (p/'input.json').write_text(json.dumps({'n':8,'moments':moments[:8]}))
        out=subprocess.check_output([sys.executable,'-I','reconstruct.py','input.json'],cwd=p,text=True)
        require(json.loads(out)==prediction,'isolated moment-only decoder replay')
        (p/'simple_certificate.py').write_bytes((ROOT/'simple_certificate.py').read_bytes())
        (p/'eight.json').write_text(json.dumps(moments[:8]))
        simpleout=subprocess.check_output([sys.executable,'-I','simple_certificate.py','eight.json'],cwd=p,text=True)
        require(json.loads(simpleout)==short,'isolated short-certificate moment-only replay')
    return {'rows_checked':len(rows),'self_converse_rows':self_count,'every_H_Aut_and_merged_mass_independently_checked':True,'largest_automorphism_order':maxaut,'moment_input':{'n':8,'moments':moments[:8]},'prediction':prediction,'short_certificate':short,'generic_twelve_moment_result':{'exponent':e,'modulus':M,'residue':residue,'prediction':first,'moments':moments},'actual_class_count_held_back_until_prediction':actual_count,'forbidden_side_channel_keys_rejected':rejected,'isolated_moment_only_replay':True,'unmerged_masses_not_replaced_by_merged_masses':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    result={'status':'PASS','source_commit':'f26bc37f4e77f72a7ce4e4f9eabf2b271acc3390','small_enumeration':small_enumeration(),'polynomial_controls':polynomial_controls(),'atlas':atlas_audit(),'scope':'Exact inference from already known finite metagraph mass moments; not an n-only count formula or a demonstrated faster enumeration algorithm.'}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(text)
    else:print(text,end='')
if __name__=='__main__':main()

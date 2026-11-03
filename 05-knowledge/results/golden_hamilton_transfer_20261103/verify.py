"""Independent graph tests and closed-set certificates. Checks survive -O."""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
import argparse, importlib.util, json, sys
import sympy as sp
from transfer import SMALL, LETTERS, STATES, M, R, product_matrix, vector3, apply3, words, targets
from extension import chart, strong
HERE=Path(__file__).resolve().parent
checks=0

def need(ok,why):
    global checks
    checks+=1
    if not ok: raise RuntimeError(why)

def join(rows,block):
    """Actual one-reversed-bridge composition in processing order."""
    m,n=len(rows),len(block)
    if not m:return tuple(block)
    out=[r|(((1<<n)-1)<<m) for r in rows]+[r<<m for r in block]
    out[m-1]&=~(1<<m);out[m]|=1<<(m-1)
    return tuple(out)

def reverse_labels(rows):
    n=len(rows);out=[0]*n
    for i in range(n):
        for j in range(n):
            if rows[i]>>j&1:out[n-1-i]|=1<<(n-1-j)
    return tuple(out)

def actual_complement(rows):
    out=list(rows)
    for i,j in combinations(range(len(rows)),2):
        if j!=i+1:out[i]^=1<<j;out[j]^=1<<i
    return tuple(out)

def source(word,extra=None):
    extra=extra or {};out=[]
    for l in reversed(word):
        r=extra[l] if l in extra else (chart(1,0) if l=='s' else chart(3,0) if l=='c' else chart(4,int(l[1:])))
        n=len(out);m=len(r)
        out=[row|(((1<<m)-1)<<n) for row in out]+[row<<n for row in r]
    return tuple(out)

def literal_target(word,extra=None):return reverse_labels(actual_complement(source(word,extra)))

def stats(rows):
    """Subset DP on actual arcs; C uses a separate marked-start DP."""
    n=len(rows)
    if not n:return (1,0,0)
    lim=1<<n;a=n-1
    end=[[0]*n for _ in range(lim)];start=[[0]*n for _ in range(lim)]
    for v in range(n):end[1<<v][v]=1
    start[1<<a][a]=1
    for mask in range(1,lim):
        for v in range(n):
            if not (end[mask][v] or start[mask][v]):continue
            opt=rows[v]&~mask
            while opt:
                bit=opt&-opt;opt-=bit;u=bit.bit_length()-1
                end[mask|bit][u]+=end[mask][v]
                start[mask|bit][u]+=start[mask][v]
    totals=[sum(x) for x in end];totals[0]=1
    h=totals[-1];b=end[-1][a]
    c=sum(sum(start[mask])*totals[(lim-1)^mask] for mask in range(lim) if mask>>a&1)
    return h-b,b,c

def h(rows):v=stats(rows);return v[0]+v[1]

def add(rows,mask):
    n=len(rows)
    return tuple(r|((1<<n) if not (mask>>i&1) else 0) for i,r in enumerate(rows))+(mask,)

def converse(rows):
    n=len(rows)
    return tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows))

def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            u,v=(i,j) if mask>>k&1 else (j,i)
            rows[u]|=1<<v
        yield tuple(rows)

def check_closed(certificate,cap,bank,absent):
    states={tuple(v) for v in certificate}
    need(len(states)==len(certificate),('duplicate certificate state',cap))
    need((1,0,0) in states,('missing initial state',cap))
    for v in states:
        need(sum(v[:2])<=cap,('cap',cap,v))
        for l,m in bank.items():
            w=apply3(v,m)
            need(sum(w[:2])>=sum(v[:2]),('H decrease',l,v,w))
            need(sum(w[:2])>cap or w in states,('closure',cap,l,v,w))
    for value in absent:need(all(sum(v[:2])!=value for v in states),('absent',value))
    return len(states)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='verification.json');args=ap.parse_args()
    ext=json.loads((HERE/'extension.json').read_text())
    reach=json.loads((HERE/'reachability.json').read_text())
    details=ext['strong5_charts'];bank=dict(SMALL,**{l:tuple(map(tuple,d['matrix'])) for l,d in details.items()})
    extra={l:tuple(d['source_rows']) for l,d in details.items()}
    for l in SMALL:need(product_matrix(M[l],R)==product_matrix(R,SMALL[l]),('quotient',l))
    count=0
    for n in range(11):
        for word in words(n):
            rows=literal_target(word)
            need(rows==targets(word),('literal chart and descriptor disagree',word))
            need(stats(rows)==vector3(word),('three literal coordinates',word,stats(rows),vector3(word)))
            count+=1
    # Three independent actual prefixes reconstruct every block matrix.
    prefixes=[literal_target(('s',)*i) for i in (1,2,3)]
    base=sp.Matrix([stats(r) for r in prefixes]);need(base.det()!=0,'basis')
    for l,mat in bank.items():
        rows=literal_target((l,),extra)
        outputs=sp.Matrix([stats(join(prefix,rows)) for prefix in prefixes])
        need(base.inv()*outputs==sp.Matrix(mat),('independent matrix',l))
        need(strong(source((l,),extra)),('strong source atom',l))
    # General seam theorem, including blocks without a chosen Hamilton spine.
    small=[r for n in range(1,5) for r in all_tournaments(n)]
    seam_checks=0
    for block in small:
        gs=[]
        extended=block
        for j in range(3):
            reverse=stats(reverse_labels(converse(extended)))
            hh=sum(reverse[:2]);start_b=reverse[1];end_cover_b=reverse[2]
            gs.append(sp.Matrix([[hh,hh-start_b,end_cover_b]]))
            extended=join(extended,(0,))
        cc=(gs[2]-gs[0])/2;aa=gs[1]-cc;bb=gs[0]-aa
        general_matrix=sp.Matrix.vstack(aa,bb,cc).T
        for prefix in small:
            v=stats(prefix);result=stats(join(prefix,block))
            need(sum(result[:2])==sum(v[i]*gs[0][i] for i in range(3)), 'single-bridge formula')
            need(sp.Matrix([v])*general_matrix==sp.Matrix([result]),'general three-coordinate closure')
            seam_checks+=1
    # Distinguish every marked Q4 chart by tiny surrounding contexts.
    phases={l:tuple(sum(vector3(w)[:2]) for w in ((l,),('s',l),(l,'s'),('s',l,'s'))) for l in SMALL if l.startswith('q')}
    need(len(set(phases.values()))==5,'five distinct Q4 phases')
    hankel=sp.Matrix([[sum(vector3(('s',)*(i+j))[:2]) for j in range(3)] for i in range(3)])
    need(hankel.det()==-4,'minimal dimension witness')
    # Matrix inequalities prove the invariant cone outside the formal empty unit.
    for l,m in bank.items():
        need(all(x>=0 for row in m for x in row),('nonnegative',l))
        need(all(row[2]>=row[0]+row[1] for row in m),('lower cone',l))
        if l!='s':
            need(all(row[0]+row[1]>=1 for row in m),('H monotonic',l))
            bound=8 if l in SMALL else ext['C_over_H_bound']
            need(all(row[2]<=bound*(row[0]+row[1]) for row in m),('upper cone',l))
        need(all(sum(row)>=1 for row in m),('W nondecrease',l))
        need(sum(m[0])>1 and sum(m[2])>1,('W strict on cone',l))
    certs={}
    for cap,absent in ((37,[7,21,37]),(133,[43,63,133]),(343,[343])):
        certs[str(cap)]=check_closed(reach['closed_certificates'][str(cap)],cap,SMALL,absent)
    certs['extended47']=check_closed(ext['closed_certificate_47'],47,bank,[7,21,47])
    selected={}
    for value,word in ext['full_bank_10001']['selected'].items():
        if word is None:continue
        rows=literal_target(word,extra);v=(1,0,0)
        for l in word:v=apply3(v,bank[l])
        need(stats(rows)==v and sum(v[:2])==int(value),('extension witness',value,word))
        need(strong(rows),('strong target',value))
        selected[value]={'word':word,'rows':rows,'state':v}
    # Enumerate WORDS, not merely collapsed states: finite by strict W growth.
    queue=[((),(1,0,0))];fiftyfive=[]
    for word,v in queue:
        if sum(v[:2])==55:fiftyfive.append((word,v))
        for l,m in SMALL.items():
            w=apply3(v,m)
            if sum(w[:2])<=55:queue.append((word+(l,),w))
    need(len(queue)==99 and len(fiftyfive)==2,'complete H55 presentation census')
    attachments=[]
    for word,v in fiftyfive:
        rows=literal_target(word);counts=[h(add(rows,mask)) for mask in range(1<<len(rows))]
        need(111 not in counts,('golden55 has111 child',word))
        sv=apply3(v,SMALL['s']);ssv=apply3(sv,SMALL['s'])
        attachments.append({'word':word,'state':v,'after_s':sum(sv[:2]),'after_ss':sum(ssv[:2]),
                            'all_attachment_counts':counts})
    prior=json.loads((HERE.parent/'metagraph_primegame_20261102/PADDED_EARS.json').read_text())
    pos=next(x for x in prior['hits'] if x['H']==55 and x['target']==111)
    need(h(tuple(pos['core']))==55 and h(add(pos['core'],pos['mask']))==111,'inherited source-preserving55->111')
    # Literal permutation count: independent of subset DP and transfer.
    permutation_controls=[]
    for word in [('s',)*7,('s','q1','s'),('q2','c')]:
        rows=literal_target(word)
        direct=sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(rows))))
        need(direct==h(rows),('permutations',word))
        permutation_controls.append([word,direct])
    # Retain the actual hostile automorphism pair, outside the claimed grammar.
    prior_dir=HERE.parent/'primegame_port_spectrum_20261103';sys.path.insert(0,str(prior_dir))
    spec=importlib.util.spec_from_file_location('hostile_control',prior_dir/'verify.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    hostile=module.hostile()
    # Aggregate rational GF over marked source WORDS, not distinct targets/classes.
    x=sp.symbols('x');Q=sum((x**len(LETTERS[l])*sp.Matrix(m) for l,m in SMALL.items()),sp.zeros(3))
    gf=sp.cancel((sp.Matrix([[1,0,0]])*(sp.eye(3)-Q).inv()*sp.Matrix([1,1,0]))[0])
    num,den=map(sp.expand,sp.fraction(gf))
    if den.subs(x,0)==-1:num,den=-num,-den
    nc=[int(num.coeff(x,i)) for i in range(sp.degree(num)+1)]
    dc=[int(den.coeff(x,i)) for i in range(sp.degree(den)+1)]
    need(sp.gcd(num,den)==1,'GF reduced')
    totals=[];vecs=[]
    for n in range(26):
        v=[int(n==0),0,0]
        for l,m in SMALL.items():
            cost=len(LETTERS[l])
            if n>=cost:
                w=apply3(vecs[n-cost],m);v=[a+b for a,b in zip(v,w)]
        vecs.append(tuple(v));totals.append(sum(v[:2]))
        rhs=(nc[n] if n<len(nc) else 0)-sum(dc[i]*totals[n-i] for i in range(1,min(n,len(dc)-1)+1))
        need(rhs==totals[n],('GF coefficients',n))
        if n<=10:need(totals[n]==sum(sum(vector3(w)[:2]) for w in words(n)),('word sums',n))
    def sq_edges(n):return [(i,j) for i,j in combinations(range(1,n+1),2) if int((i+j)**.5)**2==i+j]
    need(len(sq_edges(55))==108 and len(sq_edges(56))==111,'square graph55/56')
    need(all(len(sq_edges(n))<=2*(n-1) for n in range(2,56)),'first square surplus')
    receipt={'status':'PASS','checks':checks,'literal_golden_charts':count,'independent_block_matrices':len(bank),
             'general_seam_pairs_all_labelled_orders_1_to_4':seam_checks,
             'five_Q4_contexts':phases,'hankel_det':int(hankel.det()),'closed_certificate_sizes':certs,
             'extended_literal_witnesses':selected,'H55_all_words_count':len(queue),
             'H55_sources':attachments,'positive55to111':pos,'permutation_controls':permutation_controls,
             'hostile81_27':hostile,'aggregate_GF':{'numerator_ascending':nc,'denominator_ascending':dc,'coefficients':totals},
             'square_edge_counts':{'55':108,'56':111}}
    (HERE/args.output).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','checks','literal_golden_charts','independent_block_matrices','closed_certificate_sizes','hostile81_27','aggregate_GF')}))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent structural order-eleven residual audit. No producer imports."""
import itertools as I
import collections as C
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=C.Counter()
def require(condition,label):
    CHECKS[label]+=1
    if not condition: raise ValueError(label)
def arc(t,a,b): return (t[a]>>b)&1
def cyclic(n): return tuple(sum(1<<((v+j)%n) for j in range(1,(n+1)//2)) for v in range(n))
def valid(t):
    n=len(t)
    require(all(0<=r<1<<n and not arc(t,i,i) for i,r in enumerate(t)),'loopless bounded rows')
    require(all(arc(t,i,j)+arc(t,j,i)==1 for i in range(n) for j in range(i)),'tournament pairs')
def refine(t):
    colors=[0]*len(t)
    while True:
        blocks=[[i for i,c in enumerate(colors) if c==a] for a in sorted(set(colors))]
        sig=[(colors[i],)+tuple(sum(arc(t,i,j) for j in b) for b in blocks) for i in range(len(t))]
        names={s:i for i,s in enumerate(sorted(set(sig)))}
        nxt=[names[s] for s in sig]
        if len(names)==len(blocks): return tuple(tuple(i for i,c in enumerate(nxt) if c==a) for a in range(len(names)))
        colors=nxt
def strong(t):
    full=(1<<len(t))-1
    for reverse in (False,True):
        seen=1
        while True:
            nxt=seen
            for v in range(len(t)):
                if seen>>v&1: nxt|=(full^t[v]^(1<<v)) if reverse else t[v]
            if nxt==seen: break
            seen=nxt
        if seen!=full:return False
    return True
def relabel(t,p):
    z=[0]*len(t)
    for i in range(len(t)):
        for j in range(len(t)):
            if arc(t,i,j):z[p[i]]|=1<<p[j]
    return tuple(z)
def automorphisms(t):
    # Only degree is used for initial pruning, not the proposed block presentation.
    degree=[r.bit_count() for r in t]
    targets=[[w for w in range(len(t)) if degree[w]==degree[v]] for v in range(len(t))]
    order=sorted(range(len(t)),key=lambda v:(len(targets[v]),v))
    p=[-1]*len(t); result=[]
    def go(k,used):
        if k==len(t):result.append(tuple(p));return
        v=order[k]
        for w in targets[v]:
            if used>>w&1:continue
            if any(arc(t,v,u)!=arc(t,w,p[u]) for u in order[:k]):continue
            p[v]=w;go(k+1,used|1<<w);p[v]=-1
    go(0,0)
    return result
def defect(t,group):
    cells=refine(t);size={v:len(c) for c in cells for v in c};seen=set();orbits=[]
    for v in range(len(t)):
        if v not in seen:
            o={p[v] for p in group};seen|=o;orbits.append(tuple(sorted(o)))
    require(sum(map(len,orbits))==len(t),'full automorphism orbits partition')
    return sum(size[o[0]]-len(o) for o in orbits),tuple(sorted(map(len,orbits)))
def graph(n,bits):
    t=[0]*n
    for k,(i,j) in enumerate(I.combinations(range(n),2)):
        a,b=(i,j) if bits>>k&1 else (j,i);t[a]|=1<<b
    return tuple(t)
def canonical_small(t):
    # The only unrestricted canonicalization is at order <=8.
    require(len(t)<=8,'canonicalization input at most eight')
    cells=refine(t);best=None
    for pieces in I.product(*(I.permutations(c) for c in cells)):
        order=sum(pieces,());r=tuple(sum(arc(t,v,w)<<j for j,w in enumerate(order)) for v in order)
        if best is None or r<best:best=r
    return best
def substitution(q,blocks):
    offsets=[0]
    for b in blocks:offsets.append(offsets[-1]+len(b))
    t=[0]*offsets[-1]
    for i,b in enumerate(blocks):
        for v,row in enumerate(b): t[offsets[i]+v]=row<<offsets[i]
    for i in range(len(q)):
        for j in range(len(q)):
            if arc(q,i,j):
                for v in range(offsets[i],offsets[i+1]):t[v]|=((1<<len(blocks[j]))-1)<<offsets[j]
    return tuple(t)
def oddparts(n,start=1):
    if not n:yield ()
    for j in range(start,n+1,2):
        for tail in oddparts(n-j,j):yield (j,)+tail

def regular_seven():
    n=7;pairs=list(I.combinations(range(n),2));rows=[0]*n;d=[0]*n;left=[n-1]*n;found=[]
    def go(k):
        if k==len(pairs):
            found.append(tuple(rows));return
        a,b=pairs[k];left[a]-=1;left[b]-=1
        for u,v in ((a,b),(b,a)):
            d[u]+=1;rows[u]|=1<<v
            if all(d[x]<=3<=d[x]+left[x] for x in (a,b)):go(k+1)
            rows[u]^=1<<v;d[u]-=1
        left[a]+=1;left[b]+=1
    go(0);require(len(found)==2640,'exhaustive labelled regular seven count')
    pending=set(found);reps=[]
    while pending:
        t=min(pending);orb={relabel(t,p) for p in I.permutations(range(7))}
        require(orb<=pending,'regular seven complete permutation orbit');pending-=orb
        aut=automorphisms(t);e,s=defect(t,aut)
        require(len(orb)*len(aut)==5040,'regular seven orbit denominator')
        reps.append((t,e,len(aut),s))
    require(sorted((e,a) for _,e,a,_ in reps)==[(0,7),(0,21),(14,3)],'regular seven full group spectrum')
    return reps

def regular_nine():
    data=json.loads((HERE/'inherited/regular-nine.json').read_text())
    parents=[tuple(r['parent']) for r in data['states']]
    index={canonical_small(t):i for i,t in enumerate(parents)}
    require(len(index)==len(parents)==85,'pinned complete almost regular parent slice')
    children=[];reroot=[]
    for t in parents:
        require(C.Counter(r.bit_count() for r in t)=={3:4,4:4},'almost regular parent degrees')
        child=list(t)+[0]
        for i,row in enumerate(t):
            if row.bit_count()==3:child[i]|=1<<8
            else:child[8]|=1<<i
        child=tuple(child);children.append(child);row=C.Counter()
        for v in range(9):
            order=[i for i in range(9) if i!=v]
            deleted=tuple(sum(arc(child,i,j)<<k for k,j in enumerate(order)) for i in order)
            key=canonical_small(deleted);require(key in index,'all nine deletions return to parent slice');row[index[key]]+=1
        reroot.append(row)
    pending=set(range(85));reps=[];hist=C.Counter()
    while pending:
        i=min(pending);block=set(reroot[i]);require(i in block and block<=pending,'regular nine rooted block partition')
        require(all(reroot[j]==reroot[i] for j in block),'regular nine reroot rows identical')
        pending-=block;t=children[i];aut=automorphisms(t);e,s=defect(t,aut)
        require(s==tuple(sorted(reroot[i].values())),'regular nine literal full Aut matches reroot sizes')
        require(e==9*(len(block)-1),'regular nine error formula');reps.append((t,e,len(aut),s));hist[e]+=1
    require(hist=={0:3,18:3,36:2,72:7},'regular nine spectrum rederived')
    require(sum(e*n for e,n in hist.items())==630,'regular nine error630')
    return reps,hist

def large_single_cells(sevens,nines):
    result={};hist=C.Counter();records=[]
    for weights,core in (((9,1,1),nines),((7,3,1),sevens),((7,1,1,1,1),sevens)):
        qsize=len(weights);perms=[]
        for p in I.permutations(range(qsize)):
            if all(weights[i]==weights[p[i]] for i in range(qsize)):perms.append(p)
        pending={graph(qsize,b) for b in range(1<<(qsize*(qsize-1)//2)) if strong(graph(qsize,b))}
        labelled=len(pending);quotients=[];coalesced=0;subtotal=0;errors=C.Counter()
        while pending:
            q=min(pending);orb={relabel(q,p) for p in perms};require(orb<=pending,'weighted quotient complete orbit');pending-=orb
            actual_pattern=None
            for r,e,a,orbits in core:
                blocks=[r]+[cyclic(w) for w in weights[1:]];t=substitution(q,blocks);valid(t)
                require(strong(t),'strong substitution')
                pattern=tuple(sorted(map(len,refine(t))))
                require(actual_pattern is None or actual_pattern==pattern,'internal seven type does not alter weighted cells')
                actual_pattern=pattern;aut=automorphisms(t);actual_e,actual_o=defect(t,aut)
                require(actual_e==e,'uniform substitution preserves error')
                expected_aut=a
                for w in weights[1:]:expected_aut*=w
                if pattern==tuple(sorted(weights)):
                    require(len(aut)==expected_aut,'full substitution Aut product')
                    subtotal+=actual_e;errors[actual_e]+=1;hist[actual_e]+=1
                    records.append({'pattern':list(pattern),'rows':list(t),'aut':len(aut),'error':actual_e,'orbit_sizes':actual_o})
                else:require(weights==(7,1,1,1,1) and pattern==(1,3,7),'only seven-singleton coalescence')
            quotients.append({'rows':q,'presentation_orbit_size':len(orb),'pattern':actual_pattern})
            if actual_pattern!=tuple(sorted(weights)):coalesced+=1
        key='+'.join(map(str,weights));result[key]={'labelled_strong_quotients':labelled,'weighted_quotient_classes':len(quotients),'coalesced_quotients':coalesced,'actual_quotient_classes':len(quotients)-coalesced,'error_histogram':sorted(errors.items()),'error':subtotal,'quotients':quotients}
    require([result[k]['error'] for k in ('9+1+1','7+3+1','7+1+1+1+1')]==[630,28,308],'large-cell residuals')
    require(result['7+1+1+1+1']['weighted_quotient_classes']==24 and result['7+1+1+1+1']['coalesced_quotients']==2,'twenty-four rooted strong five minus two coalescences')
    return result,hist,records

def matrices(d):
    choices=[sum(1<<j for j in s) for s in I.combinations(range(5),d)]
    return {r for r in I.product(choices,repeat=5) if all(sum(x>>j&1 for x in r)==d for j in range(5))}
def shift(m,a,b):
    r=[0]*5
    for i,row in enumerate(m):r[(i+a)%5]=sum(1<<((j+b)%5) for j in range(5) if row>>j&1)
    return tuple(r)
def five_graph(m):
    # singleton10 -> A, B -> singleton10; A -> B matrix M
    t=list(cyclic(5))+[r<<5 for r in cyclic(5)]+[(1<<5)-1]
    for i in range(5):
        t[5+i]|=1<<10
        for j in range(5):
            u,v=(i,5+j) if m[i]>>j&1 else (5+j,i);t[u]|=1<<v
    return tuple(t)
def fives():
    result={};hist=C.Counter();records=[]
    for d in range(6):
        pending=matrices(d);count=len(pending);orbhist=C.Counter();e_hist=C.Counter();classes=0
        while pending:
            m=min(pending);orb={shift(m,a,b) for a in range(5) for b in range(5)}
            require(orb<=pending,'five matrix group orbit disjoint');pending-=orb;t=five_graph(m);valid(t)
            require(strong(t)==(d>0),'five plus singleton strong iff d positive')
            pattern=tuple(sorted(map(len,refine(t))))
            require(pattern==((11,) if d==3 else (1,5,5)),'five d3 coalesces to regular eleven')
            stabilizer=sum(shift(m,a,b)==m for a in range(5) for b in range(5))
            require(stabilizer*len(orb)==25,'five translation orbit denominator');orbhist[stabilizer]+=1;classes+=1
            if d not in (0,3):
                aut=automorphisms(t);e,orbits=defect(t,aut)
                require(len(aut)==stabilizer,'full eleven Aut equals five matrix translations')
                require(e==(40 if stabilizer==1 else 0),'five literal defect')
                e_hist[e]+=1;hist[e]+=1;records.append({'degree':d,'rows':t,'aut':len(aut),'error':e,'orbit_sizes':orbits})
        result[d]={'labelled_matrices':count,'presentation_orbits':classes,'translation_stabilizers':sorted(orbhist.items()),'actual_strong_error_histogram':sorted(e_hist.items())}
    require([result[d]['labelled_matrices'] for d in range(6)]==[1,120,2040,2040,120,1],'five complete matrix universe')
    require(hist=={0:17,40:88},'five true strong class histogram')
    return result,hist,records

def triangle_graph(q,alpha,beta,matchings):
    t=[0]*11
    for i in range(3):
        for x in range(3):t[3*i+x]|=1<<(3*i+(x+1)%3)
    for (i,j),p in zip(I.combinations(range(3),2),matchings):
        for x in range(3):
            for y in range(3):
                # Majority i->j: complement of its matching.
                win=(y!=p[x]) if arc(q,i,j) else (y==p[x])
                a,b=(3*i+x,3*j+y) if win else (3*j+y,3*i+x);t[a]|=1<<b
    t[9]|=1<<10
    for i in range(3):
        for x in range(3):
            for singleton,bits in ((9,alpha),(10,beta)):
                a,b=(3*i+x,singleton) if bits>>i&1 else (singleton,3*i+x);t[a]|=1<<b
    return tuple(t)
def skeleton_key(q,a,b,perms):
    return min((relabel(q,p),sum(((a>>i)&1)<<p[i] for i in range(3)),sum(((b>>i)&1)<<p[i] for i in range(3))) for p in perms)
def triangles():
    allmatch=list(I.permutations(range(3)));perms3=list(I.permutations(range(3)))
    transitive=(0,1,3);cycle=cyclic(3);result={};hist=C.Counter();records=[]
    for name,q in (('transitive',transitive),('cyclic',cycle)):
        qaut=[p for p in perms3 if relabel(q,p)==q];skeletons={};allcases=C.Counter()
        for a in range(8):
            for b in range(8):
                key=skeleton_key(q,a,b,qaut);skeletons.setdefault(key,(a,b))
                for match in I.product(allmatch,repeat=3):
                    t=triangle_graph(q,a,b,match);valid(t);st=strong(t);pattern=tuple(sorted(map(len,refine(t))))
                    require(st==(a!=0 and b!=7),'triangle two-singleton strong criterion')
                    expected=(1,1,9) if name=='cyclic' and a in (0,7) and b in (0,7) else (1,1,3,3,3)
                    require(pattern==expected,'triangle coalescence iff cyclic constant attachments')
                    allcases[(st,pattern)]+=1
        good=0;coalesced=0;e_hist=C.Counter()
        for key,(a,b) in sorted(skeletons.items()):
            if a==0 or b==7:continue
            if name=='cyclic' and a==7 and b==0:coalesced+=1;continue
            good+=1
            require(sum(skeleton_key(q,a,b,[p])==(q,a,b) for p in qaut)==1,'true triangle skeleton has no residual permutation')
            pending={triangle_graph(q,a,b,m) for m in I.product(allmatch,repeat=3)}
            while pending:
                t=min(pending);orb=set()
                for shifts in I.product(range(3),repeat=3):
                    p=tuple(3*i+(x+shifts[i])%3 for i in range(3) for x in range(3))+(9,10);orb.add(relabel(t,p))
                require(orb<=pending,'triangle rotation orbit disjoint');pending-=orb
                aut=automorphisms(t);e,orbits=defect(t,aut)
                require(len(aut)*len(orb)==27,'triangle full Aut orbit denominator')
                require((len(aut),e) in ((1,18),(3,0)),'triangle full group and defect')
                e_hist[e]+=1;hist[e]+=1;records.append({'majority':name,'alpha':a,'beta':b,'rows':t,'aut':len(aut),'error':e,'orbit_sizes':orbits})
        result[name]={'all_attachment_skeletons':len(skeletons),'strong_true_skeletons':good,'strong_coalesced_skeletons':coalesced,'raw_case_histogram':[[str(k),v] for k,v in sorted(allcases.items())],'error_histogram':sorted(e_hist.items())}
    require(result['transitive']['strong_true_skeletons']==49,'forty-nine transitive true skeletons')
    require(result['cyclic']['strong_true_skeletons']==16 and result['cyclic']['strong_coalesced_skeletons']==1,'seventeen cyclic strong skeletons minus one')
    require(hist=={0:780,18:260},'triangle true strong class histogram')
    return result,hist,records

def run():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for item in pins['sources']:
        require(hashlib.sha256((HERE/item['archived_path']).read_bytes()).hexdigest()==item['sha256'],'inherited source hash')
    partitions=list(oddparts(11));require(len(partitions)==12,'complete odd partitions eleven')
    seven=regular_seven();nine,ninehist=regular_nine()
    big,h1,r1=large_single_cells(seven,nine);five,h2,r2=fives();triple,h3,r3=triangles()
    histogram=h1+h2+h3;histogram.pop(0,None)
    require(histogram=={14:24,18:263,36:2,40:88,72:7},'nonregular eleven full positive distribution')
    total=sum(e*n for e,n in histogram.items());require(total==9166,'nonregular strong eleven total9166')
    a=[1,1,1,2,4]
    square=[sum(a[j]*a[k-j] for j in range(k+1)) for k in range(5)]
    require(square==[1,2,3,6,13],'A squared coefficients indexed zero through four')
    strong_errors={7:14,8:0,9:716,10:3992}
    padding=sum(square[11-k]*e for k,e in strong_errors.items())
    require(padding==10314,'SCC padding coefficient10314')
    require((total+padding)%11==10,'upper trace residue ten modulo eleven')
    require(total+padding==19480,'total eleven residual19480')
    return {'status':'PASS','scope':'Exact nonregular strong order eleven only; regular eleven not computed by this audit. Completeness of the85 inherited almost-regular parents is inherited, not proved by reroot closure. No unrestricted order-eleven census.','odd_partitions':partitions,'regular_seven_classes':[{'rows':t,'error':e,'aut':a,'orbit_sizes':s} for t,e,a,s in seven],'regular_nine_histogram':sorted(ninehist.items()),'large_single_cell_sectors':big,'five_cell_sector':five,'triangle_cell_sector':triple,'positive_nonregular_histogram':sorted(histogram.items()),'positive_nonregular_classes':sum(histogram.values()),'nonregular_strong_eleven_error':total,'A_squared_coefficients_0_through_4':square,'scc_padding':padding,'upper_trace_residue_modulo_11':10,'total_eleven_known_residual':total+padding,'regular_eleven_error':'UNKNOWN','full_automorphism_checked_representatives':len(r1)+len(r2)+len(r3),'representatives':r1+r2+r3,'checks':dict(sorted(CHECKS.items()))}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))

"""Complete zero-budget graph classification and two retained-coordinate refinements.

The native companion exhausts all 126-profile words. This wrapper pins the
universe, preserves all three deletion stages, and independently verifies the
18 chamber-only and 31 slot-only complete residuals and the final stop bank.
"""
from pathlib import Path
from math import gcd,comb
from collections import Counter
from itertools import combinations
from hashlib import sha256
import json,subprocess,sys,tempfile

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path.cwd()
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
BASELINE_SHA='89a3e7545cc77467c5f85fbe0bdaab1f71071cff65184bbc1212eb86c9f07177'
gates=0
def need(ok,why):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(why)
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def semantic(v):return sha256(canonical(v)).hexdigest()
def matrix(t,word,table):
    out={}
    for a in set(word):
        for b in set(word):
            e=gcd(a,b);p=tuple(sorted((a//e,b//e)))
            out[a,b]=table.get((t//e,*p),(0,0,0))
    return out
def components(word,M,bit):
    todo=set(range(7));out=[]
    while todo:
        reached={min(todo)};old=set()
        while reached!=old:
            old=set(reached)
            reached|={j for i in old for j in range(7) if j!=i and M[word[i],word[j]][0]&bit}
        todo-=reached;out.append(sorted(reached))
    return sorted(out)
def slot_witness(word,M):
    mult=Counter(word)
    for a in sorted(mult):
        if M[a,a][0]&2:continue
        terms=[[b,mult[b],M[a,b][2]] for b in sorted(mult) if b!=a and M[a,b][2]]
        bound=sum(n*r for b,n,r in terms)
        if mult[a]>bound:return dict(margin=a,multiplicity=mult[a],slot_bound=bound,terms=terms)
    return None
def geometry(p,q):
    L=14*p*q;A=[(max(0,(14*k-1)*q),min(L,(14*k+1)*q)) for k in range(p+1)]
    B=[(max(0,(14*k-1)*p),min(L,(14*k+1)*p)) for k in range(q+1)]
    out=[];i=j=0
    while i<len(A) and j<len(B):
        a,b=max(A[i][0],B[j][0]),min(A[i][1],B[j][1])
        if a<b:out.append((a,b))
        if A[i][1]<B[j][1]:i+=1
        elif A[i][1]>B[j][1]:j+=1
        else:i+=1;j+=1
    return L,[(out[-1][0]-L,out[0][1])]+out[1:-1]
def literal_control(n,p,q,chamber):
    L,I=geometry(p,q);walls=sorted({n*x%L for ab in I for x in ab})
    phases=[2*w+1 for w in walls] if chamber else [2*w for w in walls]
    def count(r):return sum(-((-(2*n*b-r))//(2*L))-(2*n*a-r)//(2*L)-1 for a,b in I)
    r=min(phases,key=count)
    def danger(j,v):
        den=2*L*n;a=(v*(2*L*j+r))%den
        return 14*min(a,den-a)<den
    A={j for j in range(n) if danger(j,p)};B={j for j in range(n) if danger(j,q)}
    need(len(A&B)==count(r),'literal native-grid intersection agrees with interval arithmetic')
    return dict(quotient_clock=n,p=p,q=q,phase2=r,phase_denominator=2*L,intersection=len(A&B),danger_sizes=[len(A),len(B)],marginal_bound=n//7)

def main():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(raw).hexdigest()==PROFILE_SHA,'full profile supplier pin')
    P=json.loads(raw);V=P['levels']['6']['gcds']
    baseline=json.loads((ROOT/'05-knowledge/results/continuing13_20260908_lrc_zero_cuts_certificate.json').read_bytes())['new_scales']
    need(semantic(baseline)==BASELINE_SHA and len(baseline)==7618 and max(baseline)==11935,'complete inherited necessary array')
    clocks=[t for t in baseline if t%7==0]
    need(len(clocks)==1085,'declared complete remaining zero-budget universe')
    need(len(V)==42 and all(d%7 for d in V),'whole inherited alphabet has no seven factor')
    domains=sorted({tuple(d for d in V if t%d==0) for t in clocks})
    need(len(domains)==203,'all full divisor alphabets')
    index={D:i for i,D in enumerate(domains)}
    rows=[(k,c,w) for k,L in P['levels'].items() for c,w in L['profiles']]
    lines=[str(len(rows))]+[' '.join(map(str,[k,c,*w])) for k,c,w in rows]
    lines += [str(len(domains))]+[' '.join(map(str,[len(D),*D])) for D in domains]
    lines += [str(len(clocks))]+[str(t)+' '+str(index[tuple(d for d in V if t%d==0)]) for t in clocks]
    cpp=HERE.with_suffix('.cpp')
    with tempfile.TemporaryDirectory(prefix='continuing14_lrc_') as directory:
        tmp=Path(directory).resolve()
        need(tmp.parent==Path(tempfile.gettempdir()).resolve(),'private temporary output is under the resolved system temporary directory')
        inp=tmp/'input.txt';inp.write_bytes(('\n'.join(lines)+'\n').encode())
        exe=tmp/('native.exe' if sys.platform=='win32' else 'native')
        build=subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe)],capture_output=True)
        need(build.returncode==0,'native compiler: '+build.stderr.decode(errors='replace'))
        run=subprocess.run([str(exe),str(inp),str(tmp)],capture_output=True)
        need(run.returncode==0,'native complete consumer: '+run.stderr.decode(errors='replace'))
        native=run.stdout.decode().replace('\r\n','\n')
        need('RAW_WORDS 9978787 WORD_CLOCKS 630818' in native,'whole unpruned universe and valid evaluations')
        need('CLOSED_WALL_CLOSURES 713 CHAMBER_CLOSURES 721 SLOT_CLOSURES 729' in native,'all three separately retained closure counts')
        domain_rows=[list(map(int,l.split())) for l in (tmp/'domain_counts.txt').read_text().splitlines()]
        wordbanks=[];domain_cert=[]
        for d,count,raw_count in domain_rows:
            data=(tmp/f'words_{d}.json').read_bytes();words=json.loads(data)
            need(len(words)==count and raw_count==comb(len(domains[d])+6,7),'domain whole multiset universe and accepted bank size')
            wordbanks.append(words)
            domain_cert.append(dict(domain=list(domains[d]),word_count=count,unpruned_multisets=raw_count,words_sha256=sha256(data).hexdigest()))
        need(sum(r['word_count'] for r in domain_cert)==415506,'distinct complete accepted profile words')
        zero_rows=[list(map(int,l.split())) for l in (tmp/'quotient_zero_pairs.txt').read_text().splitlines()]
        table={(n,a,b):(bits,closed,opened) for n,a,b,bits,closed,opened in zero_rows}
        need(len(table)==728,'complete nonempty quotient pair bank')
        for n,a,b,bits,c,o in zero_rows:
            need(0<=o<=c and bool(bits&1)==bool(c) and bool(bits&2)==bool(o),'oriented zero-bank counts and graph flags agree')
            need(gcd(a,b)==1 and n%7==0,'native reduced margin typing')
        records=[list(map(int,l.split())) for l in (tmp/'survey_clocks.txt').read_text().splitlines()]
        need([r[0] for r in records]==clocks,'every predeclared clock appears exactly once in order')
        plain=[r[0] for r in records if r[3]==0]
        chamber=[r[0] for r in records if r[3] and r[4]==0]
        slot=[r[0] for r in records if r[4] and r[5]==0]
        need(chamber==[1855,1911,1939,1981,2051,2121,3822,4242],'exact extra marginal-attainment clock bank')
        need(slot==[3738,3906,4074,5712,5964,6216,7308,7812],'exact extra distinct-speed slot bank')
        refinements=0
        for T,D,n,c,o,s,*_ in records:
            for t,d,nn,cc,oo,ss,*_ in records:
                if t>=T or T%t or d!=D:continue
                refinements+=1
                need(c<=cc and o<=oo and s<=ss,'unchanged-alphabet divisible refinement cannot add a residual word')
                small=matrix(t,domains[d],table);large=matrix(T,domains[d],table)
                for a in domains[d]:
                    for b in domains[d]:
                        need(large[a,b][1]<=small[a,b][1] and large[a,b][2]<=small[a,b][2],'complete oriented zero-slot counts decrease under lawful divisible refinement')
        selected=set(chamber+slot+[7056]);details=[]
        for t,d,n,c,o,s,*owners in records:
            need(n==len(wordbanks[d]) and 0<=s<=o<=c<=n,'complete clock-domain linkage and residual nesting')
            if t not in selected:continue
            residual_c=[];residual_o=[];residual_s=[];proof=[]
            for w in wordbanks[d]:
                M=matrix(t,w,table);cc=components(w,M,1);oo=components(w,M,2)
                if len(cc)==1:residual_c.append(w)
                if len(oo)!=1:
                    if len(cc)==1:proof.append(dict(word=w,chamber_components=oo))
                    continue
                residual_o.append(w);bad=slot_witness(w,M)
                if bad:proof.append(dict(word=w,slot_failure=bad))
                else:residual_s.append(w)
            need([len(residual_c),len(residual_o),len(residual_s)]==[c,o,s],'independent Python reconstruction of every refined-clock complete residual')
            details.append(dict(clock=t,closed_residuals=residual_c,chamber_residuals=residual_o,slot_residuals=residual_s,proofs=proof))
        stop=next(r for r in details if r['clock']==7056)
        need(stop['slot_residuals']==[[4,8,8,9,9,18,24],[4,8,9,9,16,18,24],[4,9,9,16,16,18,24]],'exact complete maximal remaining E0-clock residual')
        critical=[list(map(int,l.split())) for l in (tmp/'critical_banks.txt').read_text().splitlines()]
        need([r for r in critical if r[0]==1939]==[[1939,13,321,1,1,1]],'complete wall-only hostile bank')
        need([r for r in critical if r[0]==1953]==[[1953,33,320,3,1,3]],'complete singleton slot bank')
        wall=literal_control(1939,13,321,False)
        good=literal_control(1953,33,320,True)
        need(wall['intersection']==0 and min(wall['danger_sizes'])<277,'wall zero loses marginal attainment')
        need(good['intersection']==0 and good['danger_sizes']==[279,279],'singleton native arm has a saturated zero chamber')
        need(table[(1953,1,3)][2]==1 and gcd(7812,1280)==4 and gcd(7812,132)==12,'one slot is physically attainable but cannot serve three distinct speeds')
        old=literal_control(1890,48,307,True)
        need(old['intersection']==0 and old['danger_sizes']==[270,270],'old7560 depth hostile retains saturated individual zero arms')
        old_counts=[];removed=[]
        for name,stage in [('ordinary_zero_cuts',plain),('marginal_attainment',chamber),('distinct_speed_slots',slot)]:
            removed+=stage;remaining=[t for t in baseline if t not in set(removed)]
            old_counts.append(dict(mechanism=name,removed_clocks=stage,new_scales=remaining,new_scales_sha256=semantic(remaining)))
        new=old_counts[-1]['new_scales'];zero_left=[t for t in new if t%7==0]
        need(len(new)==6889 and max(new)==11934,'exact final general necessary array')
        need(len(zero_left)==356 and max(zero_left)==7056,'exact final zero-budget slice')
        cert=dict(status='PROVED analytic zero-budget consumers; FINITE-EXACT complete remaining E0 classification',
          profile_sha256=PROFILE_SHA,baseline_semantic_sha256=BASELINE_SHA,declared_clocks=clocks,
          native_source_sha256=sha256(cpp.read_bytes()).hexdigest(),native_transcript=native,
          domains=domain_cert,quotient_zero_pair_columns=['n','a','b','graph_bits','oriented_closed_count','oriented_chamber_count'],
          quotient_zero_pairs=zero_rows,clock_columns=['t','domain','word_count','closed_residual_count','chamber_residual_count','slot_residual_count','first_closed_word[7]','first_chamber_word[7]','first_slot_word[7]'],
          clocks=records,stages=old_counts,refined_complete_residuals=details,critical_native_banks=critical,
          controls=[wall,good,old],new_scales=new,new_scales_sha256=semantic(new),remaining_E0=zero_left,
          lawful_divisible_refinement_pairs=refinements,gates=gates,scope='Primitive thirteen distinct speeds; selected-six gcd clock; actual strict graph on complementary seven labels connected; all126 inherited profiles. Weak clearance only. No assertion that a remaining graph or slot assignment has an actual simultaneous ratio/phase realization.')
    data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('Complete E=0 clock classification: PASS')
    print(native.rstrip())
    print('SEPARATE NEW CLOSURES: ordinary zero cuts713; marginal attainment8; distinct-speed slots8')
    print('NECESSARY CLOCKS6889; maximum11934; remaining E0 clocks356; E0 maximum7056')
    print('STOP: full three-word7056 graph/slot residual retained; no joint ratio or phase realization inferred')
    print('NEW_SCALES_SHA256',semantic(new))
    print('CERTIFICATE_SHA256',sha256(data).hexdigest())
    print('Python always-active gates:',gates)
if __name__=='__main__':main()

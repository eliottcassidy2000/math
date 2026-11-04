#!/usr/bin/env python3
"""Independent complete profile, port-involution and physical-bank audit."""
from pathlib import Path
from tempfile import TemporaryDirectory
from itertools import permutations
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'parity_profile_transfer'

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def graph(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return a

def op(a):return [(((1<<len(a))-1)^(1<<i)^r) for i,r in enumerate(a)]
def endpoints(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        if mask&(mask-1)==0:continue
        for i in range(n):
            if mask>>i&1:d[mask][i]=sum(d[mask^(1<<i)][j] for j in range(n) if (mask>>j&1) and (a[j]>>i&1))
    return d

def profile(a):
    n=len(a);E=endpoints(a);H=list(map(sum,E));q=[[0]*(1<<n) for _ in range(n+1)];q[0][0]=1
    for k in range(1,n+1):
        for mask in range(1,1<<n):
            part=mask
            while part:
                q[k][mask]+=H[part]*q[k-1][mask^part];part=(part-1)&mask
    return [0]+[q[k][-1] for k in range(1,n+1)],E

def prime_and_strong(a):
    n=len(a);full=(1<<n)-1
    for mask in range(1,full):
        if mask.bit_count()>1 and all((mask>>v&1) or (a[v]&mask) in (0,mask) for v in range(n)):return False
    reach=[r|(1<<i) for i,r in enumerate(a)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(r==full for r in reach)

def blow(Q,blocks):
    offsets=[];n=0
    for b in blocks:offsets.append(n);n+=len(b)
    a=[0]*n
    for i,B in enumerate(blocks):
        for v,row in enumerate(B):
            a[offsets[i]+v]=row<<offsets[i]
            for j,C in enumerate(blocks):
                if Q[i]>>j&1:a[offsets[i]+v]|=((1<<len(C))-1)<<offsets[j]
    return a

def attach(a,am,cm):
    a=a[:]+[1<<(len(a)+1),0];n=len(a)-2
    for i in range(n):
        if am>>i&1:a[i]|=1<<n
        else:a[n]|=1<<i
        if cm>>i&1:a[i]|=1<<(n+1)
        else:a[n+1]|=1<<i
    return a

def main():
    names=['PROOF.md','verify.py','literal.cpp','normal.json','optimized.json'];pins={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names}
    expected={'PROOF.md':'77994a1898361f12087b49914f991b150b0a0b41b4e9e846318ebe52116cfcdb','verify.py':'99a7fc911900af9640fbb9b569c161b364323cf042d132b3665216117dae3938','literal.cpp':'a9758fc7cea06f046c4045ef74a8fe586d4c702cf3ecd2f6ec33eb955779b9dc','normal.json':'5bb7c5f6c6c2304a5fc01b358b63ae7a6304cd585befddece177fc04a30a6298','optimized.json':'5bb7c5f6c6c2304a5fc01b358b63ae7a6304cd585befddece177fc04a30a6298'}
    need(pins==expected,'reviewed pins');data=json.loads((SRC/'normal.json').read_text());endpoint={}
    for n,code,q in data['profile_controls']:
        a=graph(n,code);qq,E=profile(a);S=endpoints(op(a));need(qq==q,'complete integer cover profiles')
        need(all(sum(row[i] for row in E)%2==S[-1][i]%2 and sum(row[i] for row in S)%2==E[-1][i]%2 for i in range(n)),'all-support/Hamilton endpoint bridge')
        endpoint[n,code]=sum(((E[-1][i]+S[-1][i])&1)<<i for i in range(n))
    for n,code,p,pr in data['quotient_endpoint_records']:
        need(endpoint[n,code]==p,'full endpoint fingerprint table');need(prime_and_strong(graph(n,code))==pr,'prime and strong flags')
    cases=[];expected_values=[];triples=[]
    for control in data['literal_controls']:
        Q=control['quotient'];blocks=control['child_rows'];T=blow(Q,blocks);U=blow(op(Q),blocks)
        E=endpoints(Q)[-1];S=endpoints(op(Q))[-1];p=sum(((x+y)&1)<<i for i,(x,y) in enumerate(zip(E,S)))
        need(p==control['endpoint_parity_mask'],'actual quotient endpoint control')
        for (am,cm),(f,g) in zip(control['block_queries'],control['query_values']):
            if control['linear_law']:
                need(all(len(B)>1 for B in blocks),'nontrivial-child law scope')
                need((f-g)%4==2*(((am^cm)&p).bit_count()%2),'all complete coarse linear controls')
            elif cm==0 and am.bit_count()==1:need((f-g)%4==2*bool(am&p),'singleton-allowed endpoint controls')
        if control['ambient']:T=blow([2,4,1],[T,[0],[0]]);U=blow([2,4,1],[U,[0],[0]])
        cases.append((T,U,control['physical_queries']));expected_values.extend(control['query_values'])
        if not control['ambient'] and len(Q)<=4:
            for a,c in control['block_queries'][:4]:
                for R in (Q,op(Q)):
                    for z in range(len(Q)):triples.append((attach(R,a,c),z))
    bank=data['full_named_ten'];records=bank['basis_records_a_c_FT_FU_coefT_coefU'];queries=[row[:2] for row in records];cases.append((bank['T'],bank['U'],queries));expected_values.extend([row[2:4] for row in records])
    coarse=[(a,c) for a in range(32) for c in range(32)];expand=lambda m:sum(3<<(2*i) for i in range(5) if m>>i&1);cases.append((bank['T'],bank['U'],[(expand(a),expand(c)) for a,c in coarse]));expected_values.extend(bank['all1024_coarse_values'])
    highest=data['multiple_highest_control'];cases.extend([(highest['T'],highest['middle'],[highest['query']]),(highest['middle'],highest['U'],[highest['query']])]);expected_values.extend([x[0] for x in highest['values']])
    need(len(expected_values)*2==10348 and len(triples)==2352,'complete literal and return-word coverage')
    # The singleton boundary must not inherit the full nontrivial-child linear law.
    cases.append(([2,4,1],[4,1,2],[(1,2)]))
    lines=[str(len(cases))]
    for T,U,qs in cases:lines.extend([f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs])
    lines.append(str(len(triples)))
    for a,z in triples:lines.extend([f'{len(a)} {z}',' '.join(map(str,a))])
    with TemporaryDirectory() as td:
        exe=Path(td)/'independent';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'independent.cpp'),'-o',str(exe)],check=True);out=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
    need(len(out)==10350+2352,'independent output size');actual=[out[i:i+2] for i in range(0,10348,2)];need(actual==expected_values,'all complete exact physical counts')
    boundary=out[10348:10350];need((boundary[0]-boundary[1])%4==2,'singleton full-law hostile retained');return_counts=out[10350:];need(all(x%2==0 for x in return_counts),'all literal return swaps even')
    lookup={(r[0],r[1]):r[2:4] for r in records};hist={};recovered=[]
    for row in records:
        a,c=row[:2];co=[0,0];aa=a
        while True:
            cc=c
            while True:
                sign=1 if ((a^aa).bit_count()+(c^cc).bit_count())%2==0 else -1
                for side,v in enumerate(lookup[aa,cc]):co[side]+=sign*v
                if not cc:break
                cc=(cc-1)&c
            if not aa:break
            aa=(aa-1)&a
        need(co==row[4:],'every complete named coefficient');recovered.append([a,c,*co])
        if co[0]!=co[1]:degree=a.bit_count()+c.bit_count();hist[degree]=hist.get(degree,0)+1
    need(hist=={1:16,2:160,3:880,4:1792},'all cubic/quartic difference coordinates')
    for (a,c),values in zip(coarse,bank['all1024_coarse_values']):
        aa,cc=expand(a),expand(c);got=[sum(r[2+s] for r in recovered if r[0]&aa==r[0] and r[1]&cc==r[1]) for s in (0,1)];need(got==values,'complete polynomial at all coarse queries');need((values[0]-values[1])%4==2*(((a^c)&9).bit_count()%2),'full coarse linear congruence')
    for flags in ([],['-O']):need(subprocess.check_output([sys.executable,*flags,str(SRC/'verify.py')])==(SRC/'normal.json').read_bytes(),'producer both modes exact replay')
    need(pins=={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names},'frozen producer')
    return dict(status='PASS',integer_profile_controls=len(endpoint),literal_exact_values=10350,three_run_involution_cases=len(triples),individual_return_words_checked=sum(return_counts),full_named_coordinates=len(records),all_coarse_interpolations=1024,complete_named_nonzero_counts=hist,singleton_boundary_exact_values=boundary,multiple_highest_exact_values=highest['values'],complete_independent_named_coefficients=recovered,source_sha256=pins,scope='All-order endpoint one-query transfer for arbitrary nonempty children and full linear coarse congruence only for nontrivial children; no discarded physical cubic or quartic data.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))

"""Exact append-section tests; standard library; assertions survive python -O."""
from pathlib import Path
import argparse,hashlib,json
T=(1,0,2,1)
def need(ok,message):
    if not ok:raise RuntimeError(message)
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(word,s=T):
    for c in word:s=step(s,c)
    return s
def coords(s):
    D,E,U,d=s
    return D-E-d,U-D+E-d,E,d-E
def signed(s):return 6*s[2]-7*s[0]
def section(s):
    D,E,U,d=s
    return E==2*(U-D) and d==3*(U-D)
def bank(seed,cap):
    need(seed[2]<=cap,'seed retained')
    out={seed:''};queue=[seed]
    for s in queue:
        for c in 'OZJ':
            q=step(s,c)
            need(q[2]>s[2],'physical cover growth')
            if q[2]<=cap and q not in out:
                out[q]=out[s]+c;queue.append(q)
    return out

def predicted(seed,N):
    a,f,e,p=coords(seed);U=seed[2]
    descendants=bank(seed,max(U,2*N-3))
    targets={q[0] for q in descendants if q[0]<=N}
    if f==0 and seed[3]<=N:targets.add(seed[3])
    completions=[]
    for j in range(3):
        q=state('J'*j,seed)
        if section(q):
            n=signed(q);need(n>0,'short physical positivity')
            completions.append((j,n))
            if n<=N:targets.add(n)
    return targets,completions,len(descendants)

def main(target_cap,bank_cap):
    need(target_cap>=101 and bank_cap>=1000,'declared hostile windows')
    basis_checks=0
    for i in range(4):
        s=tuple(int(i==j) for j in range(4));D,E,U,d=s
        q=state('OJJJJ',s)
        need(q==(5*D+12*U,2*D+4*U,6*D+14*U,3*D+6*U),'rank-two transfer')
        need(signed(q)==D,'W target')
        rows=[(2*D+E-2*U,3*D-3*U+d),(D-3*E,D-3*E),(-2*D+2*E+U,-2*D+2*E+U)]
        values=[6*U-7*D,-D+7*E-U,5*D-6*E-2*U]
        for j in range(3):
            q=state('J'*j,s);DD,EE,UU,dd=q
            need((2*DD+EE-2*UU,3*DD-3*UU+dd)==rows[j],'short slice equations')
            need(signed(q)==values[j],'short signed readout')
        q=state('JJJ',s)
        need(signed(q)==3*(U-D)+E,'balanced endpoint transfer')
        basis_checks+=9
    sources=bank(T,bank_cap);face_counts=[0,0,0];min_checks=0
    for s in sources:
        a,f,e,p=coords(s);need(min(a,f,e,p)>=0,'positive source cone')
        short=[]
        for j in range(3):
            q=state('J'*j,s)
            if section(q):
                n=signed(q);need(n>0,'section physical positivity');short.append(n)
                if f>0:need(n==s[0],'positive-f short target equals D')
            if f==0:
                condition=[e==2*p,e==a+p,a==0][j]
                need(section(q)==condition,'exact f-zero face')
                if condition:
                    face_counts[j]+=1
                    need(signed(q)==[p-a,a,p][j],'face target')
        if f==0:
            q=state('ZJJJ',s)
            need(section(q) and signed(q)==s[3],'deletion construction')
            need(q[2]==14*s[2]-8*s[3],'deletion cover identity')
            minimum=min([s[3]]+short)
        else:minimum=s[0]
        need(minimum>0,'minimum positive');min_checks+=1
    words=['','O','Z','ZOZ','OZZ','OZJJ','OJJJJ','ZO','ZJJ','JOJJ']
    comparisons=[]
    for word in words:
        seed=state(word);N=target_cap;B=max(34*N-42,14*seed[2]+6*N)
        descendants=bank(seed,B)
        actual={signed(q) for q in descendants if section(q) and 0<signed(q)<=N}
        expected,completions,tail_count=predicted(seed,N)
        need(actual==expected,('complete target comparison',word,sorted(actual^expected)))
        a,f,e,p=coords(seed)
        minimum=seed[0] if f>0 else min([seed[3]]+[n for _,n in completions])
        need(actual and min(actual)==minimum,('sharp reachable minimum',word))
        if f>0:need(not any(n<seed[0] for n in actual),'no below-D repair')
        comparisons.append(dict(seed_word=word,seed=seed,target_cap=N,complete_cover_cap=B,
            descendants=len(descendants),tail_inventory=tail_count,targets=len(actual),minimum=minimum,
            short_completions=completions,targets_sha256=hashlib.sha256(json.dumps(sorted(actual)).encode()).hexdigest()))
    hostile=[]
    for r in range(31):
        scale=3**r;s=state('Z'+'OJ'*r+'OZ');balanced=state('Z'+'OJ'*r)
        need(s==tuple(scale*x for x in (9,1,11,1)),'infinite hostile source formula')
        need(coords(s)[1]==2*scale and signed(s)==3*scale<s[0],'infinite hostile inequalities')
        witness=state('OJJJJ',balanced)
        need(section(witness) and signed(witness)==signed(s),'same represented target has another witness')
        if r<4:hostile.append(dict(r=r,source=s,signed_target=signed(s),minimum_extension_target=s[0],replacement_witness=witness))
    return dict(status='PASS',basis_checks=basis_checks,complete_bank_cap=bank_cap,complete_bank_profiles=len(sources),
        face_counts=face_counts,minimum_formula_checks=min_checks,complete_target_comparisons=comparisons,
        hostile_family_checked_r=[0,30],hostile_samples=hostile,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dependency_sha256=hashlib.sha256((Path(__file__).parent/'dependencies/SIGNED-SECTION.md').read_bytes()).hexdigest(),
        scope='Exact append-only orbit theorem and complete declared finite windows; all-target occupancy remains open.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--target-cap',type=int,default=151);parser.add_argument('--bank-cap',type=int,default=3000)
    args=parser.parse_args();print(json.dumps(main(args.target_cap,args.bank_cap),indent=2,sort_keys=True))

"""Independent exact checks of the physical norm swap fixed-point theorem."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
PROOF=HERE.parent/'source-norm-capacity'/'PROOF.md'
PIN='2a75ad55f56377cd430d16a05e4ca3be3dceb20c0bd70f0a17b14c892fe90508'
CHECKS=0

def check(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

def move(s,label):
    D,E,U,d=s
    if label=='O':return D,0,2*U,d
    if label=='Z':return D-E+U,d,2*U-d,d
    if label=='J':return D-E+U,D-E,D+U,D
    raise ValueError(label)

def execute(s,word):
    for letter in word:s=move(s,letter)
    return s

def norm(s):
    D,E,U,d=s
    return (D-E+U)**2+2*(U+E)*(U-d)

def main():
    check(hashlib.sha256(PROOF.read_bytes()).hexdigest()==PIN,'frozen main proof')
    cap=5000;seed=(1,0,2,1);single=(1,1,1,1)
    bank={single,seed};queue=[seed]
    for s in queue:
        for letter in 'OZJ':
            t=move(s,letter)
            check(t[2]>s[2],'strict source growth')
            if t[2]<=cap and t not in bank:bank.add(t);queue.append(t)
    actual={s for s in bank if s!=single and s[2]==s[1]+2*s[3]}
    generated={seed};predecessors={}
    for s in bank:
        if s[0]==s[2]:
            t=execute(s,'OJ')
            check(t==(3*s[0],s[0],3*s[0],s[0]),'double balance exact profile')
            check(norm(t)==41*s[0]**2,'fixed point scalar')
            if t[2]<=cap:generated.add(t);predecessors[t]=s[0]
    check(actual==generated,'physical fixed-point image iff in complete cap')
    balanced_heights={t[0]for t in bank if t[0]==t[2]}
    for s in actual:
        D,E,U,d=s
        check((U+E)%2==0 and (U+E)//2==U-d,'literal norm coordinates fixed')
        if s==seed:check(norm(s)==13,'edge exception')
        else:
            check(s==(3*E,E,3*E,E),'classified profile')
            check(E in balanced_heights,'physical balanced antecedent exists')
    for k in range(20):
        B=execute(single,'OJ'*k)
        check(B[0]==B[2],'iterated balanced input')
        s=execute(B,'OJ');check(norm(s)==41*B[0]**2,'unbounded ray controls')
    v=execute(seed,'JOZ')
    check(v==(9,1,11,1)and norm(v)==601 and v[2]!=v[1]+2*v[3],'nonfixed diagonal hostile')
    report={'status':'PASS','checks':CHECKS,'cap':cap,'source_profiles':len(bank),
            'fixed_profiles':len(actual),'fixed_values':sorted({norm(s)for s in actual}),
            'base_proof_sha256':PIN,
            'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [HERE/'PROOF.md',HERE/'verify.py']}}
    (HERE/'VALIDATION.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':main()

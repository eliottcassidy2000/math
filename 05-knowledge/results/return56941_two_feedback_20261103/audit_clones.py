"""Literal independent test of the three named TT3 insertions and run law."""
from pathlib import Path
import sys,json,subprocess
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import need,birth,strong
def hp(T):
    v=subprocess.run([sys.argv[1]],input=' '.join(map(str,[len(T),0,1,*T]))+'\n',text=True,capture_output=True,check=True).stdout.split();return int(v[0])
def clone(T,v,m):
    labels=list(range(len(T)))+[v]*(m-1);n=len(labels);U=[]
    for i,a in enumerate(labels):
        U.append(sum(1<<j for j,b in enumerate(labels) if i!=j and (i<j if a==b else T[a]>>b&1)))
    return tuple(U)
def main():
    rows=[]
    for z in json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text()):
        T=tuple(z['next_rows']);H=z['next_H'];h=z['source_H'];b=z['budget'];n=len(T);controls=[]
        for m in (1,2,3):
            U=clone(T,n-1,m);actual=hp(U);need(actual==(h+b)*(2**m)+h and strong(U),'literal ordered-run clone law');controls.append(actual)
        hits=[json.loads(x) for x in (HERE/f'two_births_{H}.jsonl').read_text().splitlines()];need(len(hits)==3,'three saved insertion orders')
        for x in hits:
            U=birth(birth(T,x['first']),x['second']);need(hp(U)==3*H+x['offset'] and x['offset']==2*b and strong(U),'literal two-birth hit')
            block=[n-1,n,n+1];outside=list(range(n-1))
            need(all((U[a]>>j&1)==(T[n-1]>>j&1) for a in block for j in outside),'same external clone incidence')
            need(sorted(sum(U[a]>>j&1 for j in block) for a in block)==[0,1,2],'transitive three-clone block')
        rows.append(dict(H=H,base_h=h,b=b,clone_counts_m123=controls,two_birth_hits=hits))
    (HERE/'audit_clones.json').write_text(json.dumps(dict(status='FINITE_EXACT_INDEPENDENT_LITERAL_CONTROLS',rows=rows),indent=2)+'\n');print('PASS four clone laws and twelve named births')
if __name__=='__main__':main()

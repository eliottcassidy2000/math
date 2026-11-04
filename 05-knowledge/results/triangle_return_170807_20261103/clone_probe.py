from pathlib import Path
import json,subprocess,sys
HERE=Path(__file__).resolve().parent
def clone(T,v,m):
    labels=[i for i in range(len(T)) for j in range(m if i==v else 1)]
    return tuple(sum(1<<j for j,b in enumerate(labels) if (a==b and i<j) or (a!=b and T[a]>>b&1)) for i,a in enumerate(labels))
def main():
    old=json.loads((HERE/'interface.json').read_text());S=old['source'];T=old['target'];skip={old['source_names'].index(v) for step in old['trace'] for v in step['arc']};records=[]
    for v in range(len(S)):
        if v in skip:continue
        cases=[clone(U,v,m) for m in [1,2,3] for U in [S,T]]
        text=''.join(' '.join(map(str,[len(U),*U]))+'\n' for U in cases)
        hs=list(map(int,subprocess.run([sys.argv[1]],input=text,text=True,capture_output=True,check=True).stdout.split()))
        rec=dict(port=old['source_names'][v],index=v,source_counts=hs[::2],target_counts=hs[1::2],offsets=[b-3*a for a,b in zip(hs[::2],hs[1::2])]);records.append(rec);print(rec,flush=True)
    (HERE/'clone_probe.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()

"""Enumerate inherited ready grammars and send parents to an exact C++ counter."""
from pathlib import Path
import sys,json,subprocess,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ready_return_pentagon_20261103'))
import ready_inverse as inv

def main():
    exe=sys.argv[1];bound=int(sys.argv[2]) if len(sys.argv)>2 else 12
    hits=[];counts=[]
    for n in range(3,bound+1):
      for b in (1,2):
        bank=[inv.unrank(n,j) for j in range(inv.class_counts(n)[n])] if b==1 else inv.two_start(n)
        text=''.join(' '.join(map(str,[n,b,j,*T]))+'\n' for j,T in enumerate(bank))
        result=subprocess.run([exe,'6325'],input=text,text=True,capture_output=True)
        if result.returncode:raise RuntimeError(result.stderr)
        rows=[json.loads(line) for line in result.stdout.splitlines()];hits.extend(rows)
        counts.append(dict(n=n,budget=b,count=len(bank),input_sha256=hashlib.sha256(text.encode()).hexdigest(),hits=len(rows)))
        (HERE/'ready_search.json').write_text(json.dumps(dict(status='FINITE_EXACT_BOUNDED_GRAMMAR',parent_bound=n,counts=counts,hits=hits),indent=2)+'\n')
        print(n,b,len(bank),'new hits',len(rows),'total',len(hits),flush=True)
if __name__=='__main__':main()

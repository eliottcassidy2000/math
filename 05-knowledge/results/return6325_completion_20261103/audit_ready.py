"""Recreate complete finite input streams; count literal children independently."""
from search_ready import HERE,inv
from completion import need
import json,sys,subprocess,hashlib
def main():
    saved=json.loads((HERE/'ready_search.json').read_text());last=json.loads((HERE/'one_ready13.json').read_text())
    checks=[]
    for rec in saved['counts']+[last]:
        n,b=rec['n'],rec['budget']
        bank=[inv.unrank(n,j) for j in range(inv.class_counts(n)[n])] if b==1 else inv.two_start(n)
        text=''.join(' '.join(map(str,[n,b,j,*T]))+'\n' for j,T in enumerate(bank))
        need(len(bank)==rec['count'] and hashlib.sha256(text.encode()).hexdigest()==rec['input_sha256'],'identical complete grammar stream')
        out=subprocess.run([sys.argv[1]],input=text,text=True,capture_output=True)
        need(out.returncode==0,out.stderr)
        hits=[tuple(map(int,s.split())) for s in out.stdout.splitlines()]
        expected=[(h['parent_order'],h['budget'],h['rank']) for h in (last['hits'] if n==13 else saved['hits']) if h['parent_order']==n and h['budget']==b]
        need(hits==expected,'literal-child scan agrees with independent rooted-cover scan')
        checks.append(dict(n=n,budget=b,children=len(bank),hits=hits));print('audit ready',n,b,len(bank),hits,flush=True)
    exceptional=[]
    for a,c in ((3,-1),(6,1),(5,-1),(13,3),(14,3)):
        k=1
        while a*2**k+c<=6325:
            need(a*2**k+c!=6325,'exceptional ready family misses target');k+=1
        exceptional.append(dict(a=a,c=c,first_above=k))
    (HERE/'audit_ready.json').write_text(json.dumps(dict(status='FINITE_EXACT_TWO_ALGORITHMS_SHARED_PROVED_GRAMMAR',checks=checks,exceptional=exceptional),indent=2)+'\n')
if __name__=='__main__':main()

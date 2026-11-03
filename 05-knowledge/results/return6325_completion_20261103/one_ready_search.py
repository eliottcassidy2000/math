"""The complete next parent order in the inherited unique-start grammar."""
from search_ready import HERE,inv
import sys,json,subprocess,hashlib
def main():
    n=13;b=1;bank=[inv.unrank(n,j) for j in range(inv.class_counts(n)[n])]
    text=''.join(' '.join(map(str,[n,b,j,*T]))+'\n' for j,T in enumerate(bank))
    r=subprocess.run([sys.argv[1],'6325'],input=text,text=True,capture_output=True)
    if r.returncode:raise ValueError(r.stderr)
    hits=[json.loads(line) for line in r.stdout.splitlines()]
    (HERE/'one_ready13.json').write_text(json.dumps(dict(n=n,budget=b,count=len(bank),input_sha256=hashlib.sha256(text.encode()).hexdigest(),hits=hits),indent=2)+'\n')
    print('one-ready parent13',len(bank),'hits',len(hits),flush=True)
if __name__=='__main__':main()

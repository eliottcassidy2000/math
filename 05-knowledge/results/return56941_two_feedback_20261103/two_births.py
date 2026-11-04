from pathlib import Path
import json,subprocess,tempfile,sys
HERE=Path(__file__).resolve().parent
def main():
    H=int(sys.argv[1]) if len(sys.argv)>1 else 56935
    z=next(z for z in json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text()) if z['next_H']==H)
    T=z['next_rows'];exe=Path(tempfile.gettempdir())/('two-births-'+str(H)+'.exe')
    subprocess.run(['g++','-std=c++20','-O3',str(HERE/'two_births.cpp'),'-o',str(exe)],check=True)
    text=' '.join(map(str,[len(T),H,*T]))+'\n'
    with (HERE/('two_births_'+str(H)+'.jsonl')).open('w') as out:
        p=subprocess.run([str(exe)],input=text,text=True,stdout=out,check=True)
    hits=[json.loads(x) for x in (HERE/('two_births_'+str(H)+'.jsonl')).read_text().splitlines()]
    print('COMPLETE',H,'hits',len(hits),{b:sum(x['offset']==b for x in hits) for b in (2,4)},flush=True)
if __name__=='__main__':main()

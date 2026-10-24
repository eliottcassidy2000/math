#!/usr/bin/env python3
import argparse,copy,json,pathlib,random,subprocess,tempfile

def load(path,n):
    values=iter(map(int,path.read_text().split())); count=next(values); out=[]
    for _ in range(count):
        rows=[next(values) for _ in range(n)];ng=next(values)
        out.append((rows,[[next(values) for _ in range(n)] for _ in range(ng)]))
    try:next(values);raise RuntimeError('unexpected tokens')
    except StopIteration:return out

def write(path,objs):
    lines=[str(len(objs))]
    for r,g in objs:
        lines.append(' '.join(map(str,r+[len(g)])))
        lines+=[' '.join(map(str,p)) for p in g]
    path.write_text('\n'.join(lines)+'\n')

def relabel(obj,p):
    r,g=obj;n=len(r);rr=[0]*n
    for i in range(n):
        for j in range(n):rr[p[i]]|=((r[i]>>j)&1)<<p[j]
    gg=[]
    for a in g:
        b=[0]*n
        for i in range(n):b[p[i]]=p[a[i]]
        gg.append(b)
    return rr,gg

def main():
    p=argparse.ArgumentParser();p.add_argument('--producer',required=True);p.add_argument('--binary',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    base=pathlib.Path(a.producer).resolve(); binary=str(pathlib.Path(a.binary).resolve()); report=[]
    with tempfile.TemporaryDirectory(prefix='regular-eleven-hostile-') as d:
        root=pathlib.Path(d);orig=load(base/'seven-O2-parents.txt',6);states=(base/'seven-O2-states.tsv').read_text()
        def run(name,objs=orig,text=states,expected=False):
            path=root/'atlas.txt';tab=root/'states.tsv';write(path,objs);tab.write_text(text)
            c=subprocess.run([binary,'5',str(base/'sources/parents-5.txt'),str(path),str(tab),str(root/'result.json')],capture_output=True,text=True)
            if (c.returncode==0)!=expected:raise RuntimeError(f'{name}: incorrect acceptance {c.returncode} {c.stdout} {c.stderr}')
            report.append(dict(name=name,expected_acceptance=expected,returncode=c.returncode,diagnostic=c.stderr.strip() or c.stdout.strip()))
        run('pristine-seven',expected=True)
        rng=random.Random(51011);objs=[]
        for obj in orig:
            p=list(range(6));rng.shuffle(p);objs.append(relabel(obj,p))
        run('independently-relabelled-seven',objs=objs,expected=True)
        x=copy.deepcopy(orig);x[0]=(x[0][0],[list(range(6))]);run('closed-proper-subgroup-rejected',objs=x)
        x=copy.deepcopy(orig);x[1]=relabel(x[0],[1,0,2,3,4,5]);run('nonliteral-isomorphic-duplicate-rejected',objs=x)
        x=copy.deepcopy(orig);x[0][0][0]^=2;run('one-sided-arc-toggle-rejected',objs=x)
        run('parent-deletion-rejected',objs=orig[:-1])
        lines=states.splitlines();f=lines[1].split('\t');a0,b0=f[5].split(',')[0].split(':');f[5]=f[5].replace(a0+':'+b0,a0+':'+str(int(b0)+1),1);x=lines[:];x[1]='\t'.join(f);run('reroot-multiplicity-rejected',text='\n'.join(x)+'\n')
        f=lines[1].split('\t');f[4]=str((int(f[4])+1)%len(orig));x=lines[:];x[1]='\t'.join(f);run('converse-destination-rejected',text='\n'.join(x)+'\n')
        f=lines[1].split('\t');f[1]='999';x=lines[:];x[1]='\t'.join(f);run('component-id-rejected',text='\n'.join(x)+'\n')
        run('state-row-deletion-rejected',text='\n'.join(lines[:-1])+'\n')
    result=dict(status='PASS',controls=report)
    pathlib.Path(a.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()

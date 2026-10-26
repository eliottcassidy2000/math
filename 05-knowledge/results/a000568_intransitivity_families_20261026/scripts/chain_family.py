"""Chain family: TT_n with a chain 1=a0<...<ar=n reversed (S subset of {3..n-2}); checks against the census tables
(DATA dir with t{n}.tsv and the tourn_enum binary in the working directory)."""
import os, itertools, subprocess, csv
DATA=os.environ.get("DATA",".")
def chain_line(n,S):
    chain=[1]+sorted(S)+[n]; rev=set(zip(chain,chain[1:]))
    return ''.join('0' if (i,j) in rev else '1' for i in range(1,n+1) for j in range(i+1,n+1))
for n in range(5,10):
    rows=list(csv.DictReader(open(os.path.join(DATA,f't{n}.tsv')),delimiter='\t'))
    minkeys={r['key'] for r in rows if r['ncomp']=='1' and int(r['c3'])==n-2}
    subsets=[frozenset(c) for k in range(0,n-3) for c in itertools.combinations(range(3,n-1),k)]
    open('chain.txt','w').write('\n'.join(chain_line(n,S) for S in subsets)+'\n')
    subprocess.run(f'./tourn_enum {n} --table chain_t.tsv < chain.txt > /dev/null',shell=True,check=True)
    ch=list(csv.DictReader(open('chain_t.tsv'),delimiter='\t')); keys=[r['key'] for r in ch]
    def gapok(S):
        c=[1]+sorted(S)+[n]; return all(b-a<=2 for a,b in zip(c,c[1:]))
    pred=[gapok(S) for S in subsets]; act=[r['sig'].startswith(f'P{n}:') and r['sig'].count('.')==n for r in ch]
    print(f'n={n}: chains={len(subsets)} strong&c3=n-2: {all(int(r["c3"])==n-2 and r["ncomp"]=="1" for r in ch)}; distinct: {len(set(keys))==len(keys)}; equal to all min-c3 strong classes: {set(keys)==minkeys}; prime iff gaps<=2: {pred==act} (primes {sum(act)})')

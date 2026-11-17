import csv,io,subprocess,argparse
from itertools import permutations,product
BASE='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities'
def canon(rows):
 n=len(rows); cells=[tuple(i for i,v in enumerate(rows)if v.bit_count()==k)for k in range(n)]
 return min(tuple(sum((rows[i]>>j&1)<<k for k,j in enumerate(p))for i in p)
          for q in product(*(permutations(c) for c in cells if c)) for p in [sum(q,())])
def inp(T):return '\n'.join(' '.join(map(str,x))for x in T).encode()
def main():
 cli=argparse.ArgumentParser();cli.add_argument('repo');cli.add_argument('binary');cli.add_argument('--ref',default='HEAD');a=cli.parse_args()
 def read(f):return subprocess.check_output(['git','-C',a.repo,'show',a.ref+':'+BASE+'/'+f])
 layer={(0,)}
 for n in range(1,9):
  if n>=3:
   b=inp(layer)
   if n==7:b=read('n7_representatives.txt')
   if n==8:b=inp([tuple(int(r['code'])>>(8*i)&255 for i in range(8))for r in csv.DictReader(io.StringIO(read('n8_classes.csv').decode()))])
   print(subprocess.check_output([a.binary,str(n)],input=b).decode(),end='',flush=True)
  if n<6:layer={canon(tuple(v+((1-(m>>i&1))<<n)for i,v in enumerate(T))+(m,))for T in layer for m in range(1<<n)}
if __name__=='__main__':main()

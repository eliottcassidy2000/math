#!/usr/bin/env python3
import hashlib,itertools,json,pathlib,subprocess,collections
HERE=pathlib.Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def canon(r):
 d=collections.defaultdict(list)
 for i,z in enumerate(r):d[z.bit_count()].append(i)
 return min(tuple(sum(((r[i]>>j)&1)<<k for k,j in enumerate(p))for i in p)for blocks in itertools.product(*(itertools.permutations(d[k])for k in sorted(d)))for p in [sum(blocks,())])
def call(mode,n=None,bank=None):
 inp=None if bank is None else f'{n} {len(bank)}\n'+''.join(' '.join(map(str,r))+'\n'for r in bank)
 raw=subprocess.check_output([str(HERE/'certify'),mode],input=inp,text=True)
 label='small'if n is None else mode+str(n)
 (HERE/(label+'.out')).write_text(raw)
 return raw

def main():
 bankfile=HERE/'bank8.json';digest=hashlib.sha256(bankfile.read_bytes()).hexdigest()
 need(digest=='55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903','inherited audited bank hash')
 B=json.loads(bankfile.read_text());need(len(B)==6880,'all order8 representatives')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra',str(HERE/'certify.cpp'),'-o',str(HERE/'certify')],check=True)
 small=call('small');base=call('bank',8,B);raw=call('extend',8,B)
 C=sorted({canon(tuple(map(int,l.split()[1:])))for l in raw.splitlines()if l.startswith('LOW ')})
 need(len(C)==545,'complete low order9 class count')
 (HERE/'low9_bank.json').write_text(json.dumps(C,separators=(',',':'))+'\n')
 ten=call('extend',9,C)
 receipt={'status':'PASS','bank8_sha256':digest,'small':small.splitlines(),'summaries':[next(l for l in x.splitlines()if l.startswith('SUMMARY'))for x in(base,raw,ten)],'low_order9_classes':len(C),'targets':[{'order':n,'rows':[l for l in x.splitlines()if l.startswith('TARGET')]}for n,x in [(8,base),(9,raw),(10,ten)]], 'scope':'Strong rooted pairs with (H,F)=(191,1273) or (201,5303) are absent through order10; the all-order reduction is in PROOF.md.'}
 (HERE/'census_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items()if k!='targets'},sort_keys=True))
if __name__=='__main__':main()

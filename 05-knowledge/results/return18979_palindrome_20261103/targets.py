from search_two_entry import HERE,inv,need
import json
def main():
    out=[]
    for z in json.loads((HERE/'one_start_bank.json').read_text())['rows']:
        if z['kind']!='hit':continue
        T=inv.unrank(z['n'],z['rank']);U,r=inv.birth(T,0)
        p=inv.rooted_data(U,r);need(p[0]==z['target'] and p[3]-p[0]==1,'literal one-ready target')
        inv.tables.cache_clear();out.append(dict(H=z['target'],budget=1,rows=U,root=r,inport=0,profile=p,grammar=z))
    data=json.loads((HERE/'two_entry_bank.json').read_text())
    for z in data['hits']:
        c=data['cores'][z['core_id']];word=tuple(z['tail_bits']>>i&1 for i in range(z['tail_length']))
        T,p=inv.funnel(tuple(c['rows']),set(c['ports']),word);U,r=inv.birth(T,p);v=inv.rooted_data(U,r)
        need(v[0]==z['target'] and v[3]-v[0]==2,'literal two-ready target');inv.tables.cache_clear()
        out.append(dict(H=z['target'],budget=2,rows=U,root=r,inport=p,profile=v,grammar=z))
    (HERE/'targets.json').write_text(json.dumps(out,indent=2)+'\n');print('targets',[(x['H'],x['budget'],len(x['rows']),x['profile']) for x in out],flush=True)
if __name__=='__main__':main()

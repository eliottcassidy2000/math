"""Rebuild the pinned object/path/size inventory without scanning contents."""
import json,gzip,subprocess
import history_scan as h
def run():
    raw=h.frozen_git('rev-list','--objects').decode('utf8','replace').splitlines()
    objects=[q.split(' ',1) for q in raw if ' ' in q and q.split(' ',1)[1].lower().endswith(('.md','.txt','.rst'))]
    result=subprocess.check_output(['git','cat-file','--batch-check=%(objectname) %(objecttype) %(objectsize)'],input=('\n'.join(q[0] for q in objects)+'\n').encode(),cwd=h.ROOT).decode().splitlines()
    if len(result)!=len(objects):raise RuntimeError('inventory length')
    out=[]
    for (oid,path),line in zip(objects,result):
        a,kind,size=line.split()
        if a!=oid:raise RuntimeError('object ordering')
        out.append({'oid':oid,'path':path,'type':kind,'bytes':int(size)})
    data=(json.dumps(out,separators=(',',':'))+'\n').encode()
    (h.HERE/'history_universe_manifest.json.gz').write_bytes(gzip.compress(data,mtime=0))
    print(json.dumps({'objects':len(out),'bytes':sum(q['bytes'] for q in out)}))
if __name__=='__main__':run()

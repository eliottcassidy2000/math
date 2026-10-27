from pathlib import Path
s=Path('tc_window_scan.cpp').read_text()
s=s.replace('#define MAXN WORDSIZE\n#include <nauty/nauty.h>\n#include <nauty/nautinv.h>\n','')
a=s.index('static long long autOrder()'); b=s.index('static bool strong()',a)
s=s[:a]+'''// Independent exact orbit-stabilizer counter. Color refinement is only pruning.
static int colors[32], mp[32], used[32];
static bool viable(int a,int b) {
    if(colors[a]!=colors[b] || used[b]) return false;
    for(int z=0;z<n;z++) if(mp[z]>=0 && (((outm[a]>>z)&1)!=((outm[b]>>mp[z])&1))) return false;
    return true;
}
static bool existsAut(int remaining) {
    if(!remaining) return true;
    int pick=-1, best=33; uint32_t choices=0;
    for(int a=0;a<n;a++) if(mp[a]<0) {
        uint32_t cand=0; for(int b=0;b<n;b++) if(viable(a,b)) cand|=1u<<b;
        int c=__builtin_popcount(cand); if(!c) return false;
        if(c<best){best=c;pick=a;choices=cand;if(c==1)break;}
    }
    for(uint32_t cs=choices;cs;cs&=cs-1){int b=__builtin_ctz(cs);mp[pick]=b;used[b]=1;
        if(existsAut(remaining-1)){mp[pick]=-1;used[b]=0;return true;}
        mp[pick]=-1;used[b]=0;
    }
    return false;
}
static long long autOrder() {
    for(int a=0;a<n;a++) colors[a]=__builtin_popcount(outm[a]);
    for(int iter=0;iter<n;iter++){
        std::map<vector<int>,int> palette; int nc[32];
        for(int a=0;a<n;a++){
            vector<int> sig(1+2*n,0); sig[0]=colors[a];
            for(int b=0;b<n;b++) if(a!=b) sig[1+2*colors[b]+((outm[a]>>b)&1)]++;
            if(!palette.count(sig)){int id=palette.size();palette[sig]=id;} nc[a]=palette[sig];
        }
        bool stable=true;for(int a=0;a<n;a++) if(nc[a]!=colors[a]) stable=false;
        for(int a=0;a<n;a++)colors[a]=nc[a];if(stable)break;
    }
    long long order=1;
    for(int base=0;base<n;base++){
        int orbit=1;
        for(int target=base+1;target<n;target++){
            std::fill(mp,mp+32,-1);std::fill(used,used+32,0);
            for(int b=0;b<base;b++){mp[b]=b;used[b]=1;}
            if(!viable(base,target))continue;
            mp[base]=target;used[target]=1;
            if(existsAut(n-base-1))orbit++;
        }
        order*=orbit;
    }
    return order;
}
''' +s[b:]
s=s.replace('uint64_t','u128')
s=s.replace('using std::vector; using std::string;','using std::vector; using std::string;\nusing u128 = unsigned __int128;')
s=s.replace('if (H < minH) minH = H;','if (H > UINT64_MAX) {fprintf(stderr,"COUNT EXCEEDS uint64\\n");return 3;}\n        if (H < minH) minH = H;')
Path('tc_window_scan_portable.cpp').write_text(s)

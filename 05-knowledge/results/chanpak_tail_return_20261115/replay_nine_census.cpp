#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>
using namespace std;
using Row=array<int,10>;
using U=uint64_t;
int n;
U canonical(const Row& a, vector<int> p) {
    while (true) {
        vector<pair<U,int>> sig;
        for(int i=0;i<n;i++) {
            U x=0;
            for(int b:p) x=(x<<4)|__builtin_popcount(unsigned(a[i]&b));
            sig.emplace_back(x,i);
        }
        sort(sig.begin(),sig.end());
        vector<int> q;
        for(int b:p) {
            U old=~U(0);
            for(auto [key,i]:sig) if((b>>i)&1) {
                if(key!=old) q.push_back(0);
                q.back()|=1<<i;old=key;
            }
        }
        if(q==p) break;
        p=move(q);
    }
    if((int)p.size()==n) {
        U x=0;
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
            x=(x<<1)|!!(a[__builtin_ctz(unsigned(p[i]))]&p[j]);
        return x;
    }
    int i=0; while(__builtin_popcount(unsigned(p[i]))==1) i++;
    U best=~U(0);
    for(int b=p[i];b;b&=b-1) {
        int bit=b&-b;
        vector<int> q=p;
        q[i]^=bit;q.insert(q.begin()+i,bit);
        best=min(best,canonical(a,move(q)));
    }
    return best;
}
Row decode(U a,int n) {
    Row r={};
    for(int i=n-2;i>=0;i--) for(int j=n-1;j>i;j--) {
        if(a&1) r[i]|=1<<j;else r[j]|=1<<i;
        a>>=1;
    }
    return r;
}
int ham(const Row& a) {
    int paths[512][9]={};
    for(int m=1;m<(1<<n);m++) for(int b=m;b;b&=b-1) {
        int v=__builtin_ctz(unsigned(b));int s=m^(1<<v);
        if(!s) paths[m][v]=1;
        else for(int t=s&a[v];t;t&=t-1)
            paths[m][v]+=paths[s][__builtin_ctz(unsigned(t))];
    }
    int h=0;
    for(int v=0;v<n;v++)h+=paths[(1<<n)-1][v];
    return h;
}
int main(){
    unordered_set<U> old={0};
    const int expected[]={0,1,1,2,4,12,56,456,6880,191536};
    for(n=1;n<=9;n++) {
        if(n>1) {
            unordered_set<U> next;
            for(U c:old){
                Row r=decode(c,n-1);
                for(int b=0;b<(1<<(n-1));b++) {
                    Row s=r;s[n-1]=b;
                    for(int i=0;i<n-1;i++)if(!((b>>i)&1)) s[i]|=1<<(n-1);
                    next.insert(canonical(s,{(1<<n)-1}));
                }
            }
            old=move(next);
        }
        if((int)old.size()!=expected[n]) return 1;
        int ma=0;
        map<int,int> hits;
        map<int,Row> witness;
        for(U c:old) {
            Row r=decode(c,n);int h=ham(r);
            ma=max(ma,h);
            if((h>=1255&&h<=1271)||(h>=2275&&h<=2291)) {
                hits[h]++;
                if(!witness.count(h)||r<witness[h])witness[h]=r;
            }
        }
        cout<<"order "<<n<<" classes "<<old.size()<<" maxH "<<ma<<endl;
        for(auto [h,k]:hits) {
            cout<<"  "<<h<<" classes "<<k<<" witness";
            for(int i=0;i<n;i++)cout<<" "<<witness[h][i];
            cout<<endl;
        }
    }
    return 0;
}

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

// Independent all-permutation canonicalization. All words for a permutation
// are traversed in Gray order, so each new candidate needs one bit flip.
void words(int n, const string& path) {
    int bits=(n-1)*(n-2)/2, q=n*(n-1)/2, count=1<<bits;
    uint32_t all=(1u<<q)-1;
    vector<uint32_t> lo(count,all),hi(count,0),p(n);
    vector<vector<int>> chord(n,vector<int>(n,-1));int k=0;
    for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)chord[i][j]=k++;
    for(int i=0;i<n;i++)p[i]=i;
    uint64_t permutations=0;
    do {
        uint32_t initial=0;vector<uint32_t> change(bits,0);int pos=q-1;
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,pos--) {
            int a=p[i],b=p[j];bool flipped=a>b;if(flipped)swap(a,b);
            uint32_t mask=1u<<pos;
            if(b==a+1) { if(!flipped)initial|=mask; }
            else { if(flipped)initial|=mask;change[chord[a][b]]=mask; }
        }
        uint32_t value=initial;lo[0]=min(lo[0],value);hi[0]=max(hi[0],value);
        for(uint32_t j=1;j<(uint32_t)count;j++) {
            value^=change[__builtin_ctz(j)];uint32_t w=j^(j>>1);
            lo[w]=min(lo[w],value);hi[w]=max(hi[w],value);
        }
        permutations++;
    } while(next_permutation(p.begin(),p.end()));
    ofstream out(path);
    out << "word\tordinary\tconverse\tmerged\n";
    for(int w=0;w<count;w++)out<<w<<'\t'<<lo[w]<<'\t'<<(all^hi[w])<<'\t'<<min(lo[w],all^hi[w])<<'\n';
    cerr<<"n="<<n<<" permutations="<<permutations<<" words="<<count<<" candidates="<<permutations*count<<'\n';
}

void codes() {
    int n;uint32_t input;
    while(cin>>n>>input) {
        int q=n*(n-1)/2;uint32_t all=(1u<<q)-1,lo=all,hi=0;
        vector<int> p(n);vector<vector<int>> ix(n,vector<int>(n,-1));int k=0;
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)ix[i][j]=k++;
        for(int i=0;i<n;i++)p[i]=i;
        do {
            uint32_t value=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++) {
                int a=p[i],b=p[j];bool flip=a>b;if(flip)swap(a,b);
                value=(value<<1)|(((input>>ix[a][b])&1)^flip);
            }
            lo=min(lo,value);hi=max(hi,value);
        }while(next_permutation(p.begin(),p.end()));
        cout<<lo<<'\t'<<(all^hi)<<'\t'<<min(lo,all^hi)<<'\n';
    }
}
int main(int argc,char** argv) {
    if(argc==2 && string(argv[1])=="codes")codes();
    else if(argc==3)words(stoi(argv[1]),argv[2]);
    else return 1;
}

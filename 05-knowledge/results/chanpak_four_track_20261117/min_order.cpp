#include <iostream>
#include <array>
#include <vector>
#include <utility>
using namespace std;

pair<int,int> profile(const array<unsigned,7>& a,int n) {
    const unsigned N=1u<<n,F=N-1;
    int tab[128][7]{},h[128]; h[0]=1;
    for(int j=0;j<n;j++)tab[1u<<j][j]=1;
    for(unsigned s=1;s<N;s++){
        h[s]=0;
        for(unsigned bits=s;bits;bits&=bits-1){
            int j=__builtin_ctz(bits),v=tab[s][j];h[s]+=v;
            if(!v)continue;
            for(unsigned x=a[j]&~s;x;x&=x-1){unsigned b=x&-x;tab[s|b][__builtin_ctz(b)]+=v;}
        }
    }
    int p=0;
    for(unsigned s=0;s<N/2;s++)p+=2*h[s]*h[F^s];
    return {h[F],p};
}
int main(){
    for(int n=1;n<=7;n++){
        vector<pair<int,int>> e;
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)e.emplace_back(i,j);
        unsigned M=1u<<e.size();
        vector<pair<int,int>> hp(M);
        for(unsigned t=0;t<M;t++){
            array<unsigned,7>a{};
            for(unsigned i=0;i<e.size();i++){
                auto [j,k]=e[i];if(t>>i&1)swap(j,k);a[j]|=1u<<k;
            }
            hp[t]=profile(a,n);
        }
        for(unsigned t=0;t<M;t++)
            for(unsigned k=0;k<e.size();k++)if(!(t>>k&1)){
                auto [h0,p0]=hp[t];auto[h1,p1]=hp[t^(1u<<k)];
                if(p0==p1&&h0!=h1)return 1;
            }
        cout<<"order "<<n<<": "<<M<<" labeled tournaments; no nontrivial fixed-P2 reversal\n";
    }
    const array<unsigned,7>a={26,28,73,112,100,71,3},b={42,28,73,112,101,70,3};
    if(profile(a,7)!=pair(153,1812)||profile(b,7)!=pair(159,1812))return 1;
    unsigned d=0;for(unsigned i=0;i<7;i++)d+=__builtin_popcount(a[i]^b[i]);
    if(d!=4)return 1;
    cout<<"order 7 two-edit example verified\n";
}

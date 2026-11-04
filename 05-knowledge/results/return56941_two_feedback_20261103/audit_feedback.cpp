// Compare the gap compiler with literal tournament subset DP, independently.
#define main search_main
#include "feedback_words.cpp"
#undef main
#include <bit>
#include <random>
#include <stdexcept>
uint64_t literal(const vector<int>&w,int edge){
 int t=w.size(),n=t+2,p=t,q=t+1;vector<unsigned>T(n);
 for(int i=0;i<t;i++)for(int j=i+1;j<t;j++)T[i]|=1u<<j;
 for(int i=0;i<t;i++){if(w[i]&1)T[p]|=1u<<i;else T[i]|=1u<<p;if(w[i]&2)T[q]|=1u<<i;else T[i]|=1u<<q;}
 if(edge==0)T[p]|=1u<<q;else T[q]|=1u<<p;
 vector<uint64_t>d((1u<<n)*n);for(int i=0;i<n;i++)d[(1u<<i)*n+i]=1;
 for(unsigned m=1;m<(1u<<n);m++)for(unsigned z=m;z;z&=z-1){int v=countr_zero(z);auto a=d[m*n+v];for(unsigned x=T[v]&~m;x;x&=x-1)d[(m|(x&-x))*n+countr_zero(x)]+=a;}
 uint64_t h=0;for(int i=0;i<n;i++)h+=d[((1u<<n)-1)*n+i];return h;
}
uint64_t queries=0;
void check(const vector<int>&w,int pos,int c,int m){
 vector<int>v=w;v.insert(v.begin()+pos,m,c);V state{};for(int a=0;a<3;a++)state[36*a]=state[36*a+18]=1;
 for(int x:v)state=step(state,x);
 for(int edge=0;edge<2;edge++){auto h=H(state,edge);auto g=literal(v,edge);if(h[0]!=g||h[1]!=g||h[2]!=g)throw runtime_error("literal discrepancy");queries++;}
}
int main(){init();
 for(int t=0;t<=4;t++)for(unsigned code=0;code<(1u<<(2*t));code++){vector<int>w(t);for(int j=0;j<t;j++)w[j]=(code>>(2*j))&3;for(int pos=0;pos<=t;pos++)for(int c=0;c<4;c++)for(int m=1;m<=3;m++)check(w,pos,c,m);}
 mt19937 rng(56941);for(int j=0;j<120;j++){int t=5+rng()%6;vector<int>w(t);for(auto&x:w)x=rng()%4;check(w,rng()%(t+1),rng()%4,1+rng()%3);}
 // The inherited literal source A=1243, second feedback cut code 5.
 for(int m=1;m<=5;m++){vector<int>w(m,1);for(int i=0;i<11;i++)w.push_back((1-((1243>>i)&1))+2*(i!=1));auto h=literal(w,1);if(h!=3*(1243ull*(1ull<<m)-1061))throw runtime_error("uniform 3H control");queries++;}
 cout<<"{\"status\":\"FINITE_EXACT_LITERAL_MATRIX_CONTROL\",\"queries\":"<<queries<<",\"complete_word_bound\":4,\"random_long_words\":120,\"uniform_3H_bulk_max\":5}\n";
}

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using I=long long;
void require(bool b){if(!b)throw std::runtime_error("independent check failed");}
// Literal path partition: a Hamilton word is an increasing left subsequence,
// the root, and the complementary increasing right subsequence. State records
// the last left and first right incidences (empty/good/bad). No g recurrence.
I gap_count(I L){
 require(L>=1); I a=L-1; int m=0;for(I x=a;x;x>>=1)++m;
 I d[3][3]={{1,0,0},{0,0,0},{0,0,0}};
 for(int i=0;i<m;++i){
  I z[3][3]={}; bool beatsRoot=(a>>i)&1;
  for(int l=0;l<3;++l)for(int r=0;r<3;++r){
   z[beatsRoot?1:2][r]+=d[l][r];
   z[l][r?r:(beatsRoot?2:1)]+=d[l][r];
  }
  for(int l=0;l<3;++l)for(int r=0;r<3;++r)d[l][r]=z[l][r];
 }
 return d[0][0]+d[0][1]+d[1][0]+d[1][1];
}
I scalar(std::array<I,3> p,I L){return p[0]*gap_count(L)+p[2]*L-p[1]*(L&1);}
// Independent branch evaluator uses the literal partition counter at the
// actual L, not the block identity or the source's g lookup table.
std::vector<I> inverse(std::array<I,3> p,I N){
 I K=N/p[2],B=1; while(B*B<K)B*=2;
 std::vector<I> out;
 for(I b=0;b<B && b<K;++b){
  I lo=0,hi=(K-1-b)/B;
  while(lo<hi){I x=(lo+hi)/2;
   if(scalar(p,B*x+b+1)>=N)hi=x;else lo=x+1;
  }
  I L=B*lo+b+1;if(scalar(p,L)==N)out.push_back(L);
 }
 std::sort(out.begin(),out.end());return out;
}
I literal_g(I L){
 I a=L-1;int m=0;for(I x=a;x;x>>=1)++m;I ans=0;
 for(I left=0;left<(1LL<<m);++left){
  int last=-1,first=-1;
  for(int i=0;i<m;++i){if(left>>i&1)last=i;else if(first<0)first=i;}
  if((last<0 || (a>>last&1)) && (first<0 || !(a>>first&1)))++ans;
 }return ans;
}
std::vector<I> gadget(I L){I a=L-1;int m=0;for(I x=a;x;x>>=1)++m;
 std::vector<I>T(m+1);for(int i=0;i<m;++i){int v=i+1;
  if(a>>i&1)T[v]|=1;else T[0]|=1LL<<v;
  for(int j=i+1;j<m;++j)T[v]|=1LL<<(j+1);
 }return T;
}
// Inclusion-exclusion over allowed vertices counts length n-1 walks using
// every vertex. Such a walk is automatically a Hamilton path.
I walk_ie(const std::vector<I>&T){int n=T.size();I ans=0;
 for(int S=1;S<(1<<n);++S){std::vector<I>d(n);for(int v=0;v<n;++v)d[v]=(S>>v)&1;
  for(int k=1;k<n;++k){std::vector<I>z(n);
   for(int v=0;v<n;++v)if(S>>v&1)for(int w=0;w<n;++w)if((S>>w&1)&&(T[w]>>v&1))z[v]+=d[w];d=z;
  }I total=0;for(I x:d)total+=x;ans+=((n-__builtin_popcount((unsigned)S))&1)?-total:total;
 }return ans;
}
std::vector<I> seam(std::vector<I>T,int root,I L){auto G=gadget(L);int n=T.size(),m=G.size();
 for(I&row:T)row|=((1LL<<m)-1)<<n;for(I row:G)T.push_back(row<<n);
 T[root]^=1LL<<n;T[n]|=1LL<<root;return T;
}
int main(){
 I literal=0,blocks=0,monotone=0,inv=0,physical=0;
 for(I L=1;L<=4096;++L){require(gap_count(L)==literal_g(L));++literal;}
 for(int k=1;k<=9;++k){I B=1LL<<k;
  for(I b=0;b<B;++b){I q=0;for(int i=0;i<k;++i)if(!(b>>i&1))q|=1LL<<(k-1-i);
   for(I x=0;x<=128;++x){require(gap_count(B*x+b+1)==gap_count(x+1)+q*x+(q&1)*(x&1)+gap_count(b+1)-1);++blocks;}
  }
 }
 for(I h=1;h<=8;++h)for(I e=0;e<=h;++e)for(I R=h;R<=h+3;++R){std::array<I,3>p={h,e,R};
  std::vector<std::vector<I>>hits(801);
  for(I L=1;L<=800/R;++L){I N=scalar(p,L);if(N<=800)hits[N].push_back(L);}
  for(I N=0;N<=800;++N){require(inverse(p,N)==hits[N]);++inv;}
  for(I K=2;K<=100;++K){I B=1;while(B*B<K)B*=2;
   for(I b=0;b<B&&b<K;++b)for(I L=b+1;L+B<=K;L+=B){require(scalar(p,L)<scalar(p,L+B));++monotone;}
  }
 }
 for(I L=1;L<=128;++L){require(walk_ie(gadget(L))==gap_count(L));++physical;}
 std::array<std::vector<I>,3>prefix={std::vector<I>{2,4,1},std::vector<I>{12,5,8,2},std::vector<I>{18,20,25,3,8}};
 std::array<std::array<I,3>,3>p={std::array<I,3>{3,1,3},std::array<I,3>{5,1,7},std::array<I,3>{13,2,18}};
 int roots[]={0,1,0};I bigL[]={100000007,123456789,987654321};I Ns[]={486878921,1383336817,29960873285};
 for(int i=0;i<3;++i){for(I L=1;L<=32;++L){require(walk_ie(seam(prefix[i],roots[i],L))==scalar(p[i],L));++physical;}
  require(scalar(p[i],bigL[i])==Ns[i]);require(inverse(p[i],Ns[i])==std::vector<I>{bigL[i]});
 }
 std::cout<<"{\"status\":\"PASS\",\"literal_partitions\":"<<literal<<",\"block_identities\":"<<blocks<<",\"monotone_comparisons\":"<<monotone<<",\"complete_inverse_tests\":"<<inv<<",\"walk_ie_graphs\":"<<physical<<",\"large_inverses\":3}\n";
}

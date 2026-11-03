#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>
int main(int argc,char**argv){
 int n=argc>1?std::stoi(argv[1]):7;if(n<3||n>7)return 2;
 int m=n*(n-1)/2;unsigned full=(1u<<m)-1;int idx[7][7];int e=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)idx[i][j]=idx[j][i]=e++;
 std::vector<int> H(full+1),p(n);std::iota(p.begin(),p.end(),0);unsigned words=0;
 do{unsigned on=0,used=0;for(int i=0;i<n-1;i++){unsigned bit=1u<<idx[p[i]][p[i+1]];used|=bit;if(p[i]<p[i+1])on|=bit;}
   unsigned free=full^used,s=free;for(;;){H[on|s]++;if(!s)break;s=(s-1)&free;}words++;
 }while(std::next_permutation(p.begin(),p.end()));
 auto equal=[&](unsigned a,unsigned b){if(H[a]!=H[b])return false;for(int k=0;k<m;k++){unsigned t=1u<<k;if(H[a|t]-H[a&~t]!=H[b|t]-H[b&~t])return false;}return true;};
 std::unordered_map<uint64_t,std::vector<unsigned>> buckets;buckets.reserve(full+1);
 uint64_t classes=0,collisions=0,false_hash_matches=0;unsigned firsta=0,firstb=0;
 for(unsigned a=0;a<=full;a++){
   uint64_t hash=1469598103934665603ull;
   auto mix=[&](int x){hash^=uint32_t(x);hash*=1099511628211ull;};mix(H[a]);
   for(int k=0;k<m;k++){unsigned t=1u<<k;mix(H[a|t]-H[a&~t]);}
   auto &candidates=buckets[hash];bool found=false;
   for(unsigned b:candidates){if(equal(a,b)){found=true;collisions++;if(collisions==1){firsta=b;firstb=a;}break;}else false_hash_matches++;}
   if(!found){candidates.push_back(a);classes++;}
 }
 std::cout<<"{\"n\":"<<n<<",\"orientations\":"<<full+1<<",\"path_cylinders\":"<<words<<",\"first_jet_states\":"<<classes<<",\"duplicate_states\":"<<collisions<<",\"false_hash_matches\":"<<false_hash_matches<<",\"first_collision\":["<<firsta<<","<<firstb<<"],\"exact_collision_comparisons\":true}\n";
}

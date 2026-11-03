#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(){
  for(int n=3;n<=6;n++){
    int index[6][6]={},m=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)index[i][j]=index[j][i]=m++;
    vector<unsigned> cliques[7];
    for(unsigned s=0;s<(1u<<n);s++){
      int k=__builtin_popcount(s);unsigned mask=0;
      if(k<3)continue;
      for(int i=0;i<n;i++)if(s>>i&1u)for(int j=i+1;j<n;j++)if(s>>j&1u)mask|=1u<<index[i][j];
      cliques[k].push_back(mask);
    }
    unsigned all=(1u<<m)-1;vector<unsigned char> tf(1u<<m,1);
    for(unsigned x=0;x<=all;x++)for(unsigned c:cliques[3])if((x&c)==c){tf[x]=0;break;}
    uint64_t redblue=0,c335=0,c336=0;
    for(unsigned r=0;r<=all;r++)if(tf[r]){
      unsigned available=all^r,b=available;
      while(true){
        if(tf[b]){
          redblue++;unsigned g=available^b;bool good5=true,good6=true;
          for(unsigned c:cliques[5])if((g&c)==c){good5=false;break;}
          for(unsigned c:cliques[6])if((g&c)==c){good6=false;break;}
          c335+=good5;c336+=good6;
        }
        if(!b)break;b=(b-1)&available;
      }
    }
    if(n==5 && c336-c335!=1)throw runtime_error("n5 green-clique difference");
    if(n==6 && c336-c335!=1422)throw runtime_error("n6 green-clique difference");
    cout<<"{\"n\":"<<n<<",\"edge_slots\":"<<m<<",\"red_blue_triangle_free\":"<<redblue<<",\"valid_335\":"<<c335<<",\"valid_336\":"<<c336<<"}\n";
  }
}

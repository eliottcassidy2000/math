#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
using I=long long;
int main(int argc,char**argv){
 if(argc!=2)return 2;std::ifstream in(argv[1]);if(!in)return 3;
 constexpr int n=10,N=1024,full=N-1;
 std::array<unsigned,n> row{},pred{};std::array<std::array<I,n>,N> st{},en{};
 std::array<I,N> et{},us{};std::array<I,n> U{};bool got[2]={false,false};
 std::string word;I count=0,hits=0;
 while(std::getline(in,word)){
  if(!word.empty()&&word.back()=='\r')word.pop_back();if(word.size()!=45)return 4;
  row.fill(0);pred.fill(0);int k=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   if(word[k]=='1')row[i]|=1u<<j;else if(word[k]=='0')row[j]|=1u<<i;else return 5;k++;
  }
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(row[j]>>i&1)pred[i]|=1u<<j;
  st={};en={};et.fill(0);U.fill(0);us.fill(0);
  for(int i=0;i<n;i++)st[1<<i][i]=en[1<<i][i]=1;
  for(unsigned m=1;m<N;m++){
   unsigned todo=m;
   while(todo){unsigned b=todo&-todo;todo-=b;int v=__builtin_ctz(b);unsigned rest=m^b;
    if(rest){unsigned a=rest&pred[v];while(a){unsigned bb=a&-a;a-=bb;en[m][v]+=en[rest][__builtin_ctz(bb)];}
     a=rest&row[v];while(a){unsigned bb=a&-a;a-=bb;st[m][v]+=st[rest][__builtin_ctz(bb)];}}
    et[m]+=en[m][v];
   }
  }
  I H=et[full];for(int j=0;j<n;j++)U[j]=st[full][j];
  for(unsigned m=1;m<full;m++)if(et[m]){unsigned a=full^m;while(a){unsigned bb=a&-a;a-=bb;int j=__builtin_ctz(bb);U[j]+=et[m]*st[full^m][j];}}
  count++;
  for(unsigned m=1;m<N;m++){
   unsigned b=m&-m;int v=__builtin_ctz(b);us[m]=us[m^b]+U[v];I u=us[m],odd=u/(u&-u);
   int target=odd==595?0:(odd==11007?1:-1);if(target<0||got[target])continue;
   I tax=0;
   for(unsigned mask=1;mask<full;mask++){
    I a=0,bv=0;unsigned t=mask&m;while(t){unsigned bb=t&-t;t-=bb;a+=en[mask][__builtin_ctz(bb)];}
    if(!a)continue;t=(full^mask)&m;while(t){unsigned bb=t&-t;t-=bb;bv+=st[full^mask][__builtin_ctz(bb)];}tax+=a*bv;
   }
   I D=target==0?-11:-691;if(H-tax!=D)continue;
   got[target]=true;hits++;int j=__builtin_ctzll((unsigned long long)u);
   std::cout<<"{\"kind\":\"hit\",\"catalogue_position\":"<<count<<",\"rows\":[";
   for(int i=0;i<n;i++){if(i)std::cout<<',';std::cout<<row[i];}
   std::cout<<"],\"selected\":"<<m<<",\"H\":"<<H<<",\"defect\":"<<D<<",\"U\":"<<u<<",\"dyadic_shift\":"<<j<<"}\n"<<std::flush;
  }
  if(count%100000==0)std::cerr<<"head target catalogue sources "<<count<<", targets found "<<hits<<'\n';
  if(got[0]&&got[1])break;
 }
 std::cout<<"{\"kind\":\"summary\",\"catalogue_prefix_checked\":"<<count<<",\"targets_found\":"<<hits<<",\"complete_catalogue\":"<<(count==9733056?"true":"false")<<"}\n";
 return 0;
}

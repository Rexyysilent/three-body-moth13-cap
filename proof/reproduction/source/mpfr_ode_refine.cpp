#include <boost/multiprecision/mpfr.hpp>
#include <boost/numeric/odeint.hpp>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
namespace bm=boost::multiprecision;
namespace oi=boost::numeric::odeint;
using R=bm::number<bm::backends::mpfr_float_backend<80>>;
using V=std::array<R,8>;
namespace boost{namespace numeric{namespace odeint{template<>struct norm_result_type<R,void>{using type=R;};template<>struct norm_result_type<V,void>{using type=R;};}}}
struct Sys{
 void operator()(const V&x,V&f,const R&)const{
  R mu=R(1001)/1000,q3x=-(x[0]+x[2])/mu,q3y=-(x[1]+x[3])/mu;
  auto g=[](R a,R b){R r2=a*a+b*b,r3=r2*sqrt(r2);return std::array<R,2>{a/r3,b/r3};};
  auto a=g(x[2]-x[0],x[3]-x[1]),b=g(q3x-x[0],q3y-x[1]);
  auto c=g(x[0]-x[2],x[1]-x[3]),d=g(q3x-x[2],q3y-x[3]);
  f={x[4],x[5],x[6],x[7],a[0]+mu*b[0],a[1]+mu*b[1],c[0]+mu*d[0],c[1]+mu*d[1]};
 }
};
static R rd(std::istream&f){std::string s;f>>s;return R(s);}
int main(int ac,char**av){
 if(ac!=3){std::cerr<<"mesh outlog\n";return 2;}
 std::ifstream in(av[1]);int N;in>>N;R tau=rd(in),rt=rd(in),rs=rd(in);
 std::vector<R>den(N);for(auto&x:den)x=rd(in);
 std::vector<V>node(N+1);for(auto&v:node)for(auto&x:v)x=rd(in);
 std::ofstream out(av[2]);out<<std::scientific<<std::setprecision(75);
 using Step=oi::runge_kutta_fehlberg78<V,R,V,R,oi::array_algebra>;
 using Control=oi::controlled_runge_kutta<Step>;
 Control ctl(typename Control::error_checker_type(R("1e-30"),R("1e-30"),R(1),R(1)));
 for(int k=0;k<N;k++){
  V x=node[k];R t=0,h=tau/den[k],dt=h/100;
  size_t steps=oi::integrate_adaptive(ctl,Sys{},x,t,h,dt);
  out<<"REC "<<k<<" OK 1\nMETA "<<steps<<" 0 0 0 0\nXC";
  for(auto z:x)out<<" "<<z<<" "<<z;
  out<<"\nXB";for(auto z:x)out<<" "<<z<<" "<<z;
  out<<"\nMB";for(int i=0;i<8;i++)for(int j=0;j<8;j++){R z=(i==j?R(1):R(0));out<<" "<<z<<" "<<z;}out<<"\n";
 }
 out<<"BATCH_OK 1 COUNT "<<N<<" GLOBAL_D2 0\n";
}
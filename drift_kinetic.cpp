#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <tuple>

constexpr int NX=51;
constexpr int NY=51;
constexpr int NZ=51;

constexpr double x_0=NX/2;
constexpr double y_0=NY/2;
constexpr double z_0=NZ/2;

constexpr double step_size=1;
constexpr double timestep=0.1;
constexpr double mass=1;
constexpr double charge=1;

#define SQU(x) ((x) * (x))

enum class Direction {x=0,y,z};

struct Field{
  double x,y,z;

  Field():x(0),y(0),z(0){}

  Field(double _x, double _y, double _z):x(_x),y(_y),z(_z){}

  double abs(){
    return sqrt(SQU(x) + SQU(y) + SQU(z));
  }

  auto& operator+=( const Field& rhs ){
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
    return *this;
  }

  auto operator+( const Field& rhs ){
    Field res = *this;
    res+=rhs;
    return res;
  }

  auto& operator-=( const Field& rhs ){
    this->x -= rhs.x;
    this->y -= rhs.y;
    this->z -= rhs.z;
    return *this;
  }

  auto operator-( const Field& rhs ){
    Field res = *this;
    res-=rhs;
    return res;
  }

  auto& operator*=( double rhs ){
    this->x *= rhs;
    this->y *= rhs;
    this->z *= rhs;
    return *this;
  }

  auto operator*( double rhs ){
    Field res = *this;
    res*=rhs;
    return res;
  }
  
  void print(std::string name=""){
    std::cout<<name<<".x = "<<x<<"  ";
    std::cout<<name<<".y = "<<y<<"  ";
    std::cout<<name<<".z = "<<z<<"  ";
    std::cout<<std::endl;
  }
};

struct Particle{
  double x,y,z;
  double vx,vy,vz;

  Particle():x(0),y(0),z(0),vx(0),vy(0),vz(0){}
  Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz):
    x(_x),y(_y),z(_z),vx(_vx),vy(_vy),vz(_vz){}

  void print(std::string name="pt"){
    std::cout<<name<<".x = "<<x<<"  ";
    std::cout<<name<<".y = "<<y<<"  ";
    std::cout<<name<<".z = "<<z<<"  ";
    std::cout<<name<<".vx = "<<vx<<"  ";
    std::cout<<name<<".vy = "<<vy<<"  ";
    std::cout<<name<<".vz = "<<vz<<"  ";
    std::cout<<std::endl;
  }

  double energy(){
    return 0.5 * mass * (SQU(vx) + SQU(vy) + SQU(vz));
  }
};



struct Vecfield{
  std::vector<Field> fieldvec;

  Vecfield(){
    fieldvec.resize(NX*NY*NZ,Field());
  }

  //generate Bfield structure here
  //quadratic, around (x0,y0), 0 in z
  void generate_bfield(){
    //double fac = 1 / (sqrt(SQU(NX) + SQU(NY)));
    for(auto x=0;x<NX;++x){
      for(auto y=0;y<NY;++y){
	for(auto z=0;z<NZ;++z){
	  double xd = static_cast<double>(x)-x_0;
	  double yd = static_cast<double>(y)-y_0;
	  //double zd = static_cast<double>(z)-z_0;
	  //double B = fac * sqrt(SQU(xd) + SQU(yd));
	  double B=0.5;
	  double angle = atan2(yd,xd);
	  double Bx = -1*B*sin(angle);
	  double By = B*cos(angle);
	  val(x,y,z) = Field(Bx,By,0);
	}
      }
    }
  }

  //basic operations
  int lin_index(int x, int y, int z){
    return x + NX*y + z*NX*NY;
  }

  int size() const {
    return fieldvec.size();
  }

  Field& val(int x, int y, int z){
    return fieldvec[lin_index(x,y,z)];
  }
  
  //obtain value of vectorfield via CIC @ (x,y,z)
  Field at(double x, double y, double z){
    if( x>=NX-1 || y>=NY-1 || z>=NZ-1 || x<0 || y<0 || z<0 ){
      std::cout<<"Interpolation index is ooB at (x,y,z) = ("<<x<<","<<y<<","<<z<<")."<<std::endl;
      exit(-1);
    }
    int xi = static_cast<int>(x);
    int yi = static_cast<int>(y);
    int zi = static_cast<int>(z);

    auto& bl = val(xi,yi,zi); //bottom left in front
    auto& br = val(xi+1,yi,zi); //bottom right in front
    auto& tl = val(xi,yi+1,zi); //top left in front
    auto& tr = val(xi+1,yi+1,zi); //top right in front

    auto& blb = val(xi,yi,zi+1); //bottom left in back
    auto& brb = val(xi+1,yi,zi+1); //bottom right in back
    auto& tlb = val(xi,yi+1,zi+1); //top left in back
    auto& trb = val(xi+1,yi+1,zi+1); //top right in back

    double hx = x - xi;
    double hy = y - yi;
    double hz = z - zi;
    double vol = 1 / (step_size*step_size*step_size);

    //weights @ front
    double w_bl = vol*(1-hx)*(1-hy)*hz/vol;
    double w_br = vol*hx*(1-hy)*hz/vol;
    double w_tl = vol*(1-hx)*hy*hz/vol;
    double w_tr = vol*hx*hy*hz/vol;
      
    //weights @ back
    double w_blb = vol*(1-hx)*(1-hy)*(1-hz);
    double w_brb = vol*hx*(1-hy)*(1-hz);
    double w_tlb = vol*(1-hx)*hy*(1-hz);
    double w_trb = vol*hx*hy*(1-hz);

    //double weight_sum = w_bl + w_br + w_tl + w_tr + w_blb + w_brb + w_tlb + w_trb;
    //std::cout<<"weight sum = "<<weight_sum<<std::endl;
    
    return bl*w_bl + br*w_br + tl*w_tl + tr*w_tr + blb*w_blb + brb*w_brb + tlb*w_tlb + trb*w_trb;
  }

  Field at(const Particle& pt) {
    return at(pt.x,pt.y,pt.z);
  }
  
  Field at(const Field& pt) {
    return at(pt.x,pt.y,pt.z);
  }

  auto& operator[]( int index ){
      return fieldvec[index];
  }

  const auto& operator[]( int index ) const {
      return fieldvec[index];
  }

  auto& operator+=( const Vecfield& rhs ){
    if(this->size() != rhs.size()){
      std::cout<<"Cannot add Vecfields of different dimensions"<<std::endl;
      exit(-1);
    }

    for(auto i=0;i<this->size();++i){
      (*this)[i]+=rhs[i];
    }
    return *this;
  }

  auto operator+( const Vecfield& rhs ){
    if(this->size() != rhs.size()) return *this;
    Vecfield result = *this;
    result+=rhs;
    return result;
  }

  auto& operator-=( const Vecfield& rhs ){
    if(this->size() != rhs.size()){
      std::cout<<"Cannot add Vecfields of different dimensions"<<std::endl;
      exit(-1);
    }

    for(auto i=0;i<this->size();++i){
      (*this)[i]-=rhs[i];
    }
    return *this;
  }

  auto operator-( const Vecfield& rhs ){
    if(this->size() != rhs.size()) return *this;
    Vecfield result = *this;
    result-=rhs;
    return result;
  }

  auto& operator*=( double rhs ){
    for(auto i=0;i<this->size();++i){
      (*this)[i]*=rhs;
    }
    return *this;
  }

  auto operator*( double rhs ){
    Vecfield res = *this;
    res *= rhs;
    return res;
  }

  //derivative of coordinate dir_coord in direction dir_deriv
  template <Direction dir_coord, Direction dir_deriv>
  double deriv_vf(int x, int y, int z){
    //oob case
    //std::cout<<"x="<<x<<"  y="<<y<<"  z="<<z<<"  dir_coord="<<static_cast<int>(dir_coord)<<"  dir_deriv="<<static_cast<int>(dir_deriv)<<std::endl;
    if(x<0 || x>=NX || y<0 || y>=NY || z<0 || z>=NZ){
      return 0;
    }

    double div=1/(2*step_size);

    if(dir_deriv==Direction::x){
      if(x>0 && x<(NX-1)){
	const Field& val1 = val(x+1,y,z);
	const Field& val2 = val(x-1,y,z);
	if(dir_coord==Direction::x) return div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return div*(val1.z-val2.z);
      }else if(x==0){
	const Field& val1 = val(x+1,y,z);
	const Field& val2 = val(x,y,z);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else if(x==NX-1){
	const Field& val1 = val(x,y,z);
	const Field& val2 = val(x-1,y,z);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else{
	return 0;
      }
    }
    if(dir_deriv==Direction::y){
      if(y>0 && y<(NY-1)){
	const Field& val1 = val(x,y+1,z);
	const Field& val2 = val(x,y-1,z);
	if(dir_coord==Direction::x) return div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return div*(val1.z-val2.z);
      }else if(y==0){
	const Field& val1 = val(x,y+1,z);
	const Field& val2 = val(x,y,z);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else if(y==NY-1){
	const Field& val1 = val(x,y,z);
	const Field& val2 = val(x,y-1,z);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else{
	return 0;
      }
    }
    if(dir_deriv==Direction::z){
      if(z>0 && z<(NZ-1)){
	const Field& val1 = val(x,y,z+1);
	const Field& val2 = val(x,y,z-1);
	if(dir_coord==Direction::x) return div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return div*(val1.z-val2.z);
      }else if(z==0){
	const Field& val1 = val(x,y,z+1);
	const Field& val2 = val(x,y,z);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else if(z==NZ-1){
	const Field& val1 = val(x,y,z);
	const Field& val2 = val(x,y,z-1);
	if(dir_coord==Direction::x) return 2*div*(val1.x-val2.x);
	if(dir_coord==Direction::y) return 2*div*(val1.y-val2.y);
	if(dir_coord==Direction::z) return 2*div*(val1.z-val2.z);
      }else{
	return 0;
      }
    }
  }

  template <Direction dir>
  double deriv_scalar(const std::vector<double>& vec, int x, int y, int z){
    if(vec.size() != NX*NY*NZ) return 0;

    double div=1/(2*step_size);

    if(dir==Direction::x){
      if(x>0 && x<(NX-1)){
	double val1 = vec[lin_index(x+1,y,z)];
	double val2 = vec[lin_index(x-1,y,z)];
	return div*(val1-val2);
      }else if(x==0){
	double val1 = vec[lin_index(x+1,y,z)];
	double val2 = vec[lin_index(x,y,z)];
	return 2*div*(val1-val2);
      }else if(x==NX-1){
	double val1 = vec[lin_index(x,y,z)];
	double val2 = vec[lin_index(x-1,y,z)];
	return 2*div*(val1-val2);
      }else{
	return 0;
      }
    }
    if(dir==Direction::y){
      if(y>0 && y<(NY-1)){
	double val1 = vec[lin_index(x,y+1,z)];
	double val2 = vec[lin_index(x,y-1,z)];
	return div*(val1-val2);
      }else if(y==0){
	double val1 = vec[lin_index(x,y+1,z)];
	double val2 = vec[lin_index(x,y,z)];
	return 2*div*(val1-val2);
      }else if(y==NY-1){
	double val1 = vec[lin_index(x,y,z)];
	double val2 = vec[lin_index(x,y-1,z)];
	return 2*div*(val1-val2);
      }else{
	return 0;
      }
    }
    if(dir==Direction::z){
      if(z>0 && z<(NZ-1)){
	double val1 = vec[lin_index(x,y,z+1)];
	double val2 = vec[lin_index(x,y,z-1)];
	return div*(val1-val2);
      }else if(z==0){
	double val1 = vec[lin_index(x,y,z+1)];
	double val2 = vec[lin_index(x,y,z)];
	return 2*div*(val1-val2);
      }else if(z==NZ-1){
	double val1 = vec[lin_index(x,y,z)];
	double val2 = vec[lin_index(x,y,z-1)];
	return 2*div*(val1-val2);
      }else{
	return 0;
      }
    }
  }

  Vecfield unit_field(){
    Vecfield new_field=*this;
    std::for_each(new_field.fieldvec.begin(),new_field.fieldvec.end(),[](Field& f){
	double val=f.abs();
	if(val==0) val=1;
	f.x/=val;
	f.y/=val;
	f.z/=val;
      }
    );
    return new_field;
  }

  std::vector<double> abs_val(){
    std::vector<double> field_abs(fieldvec.size(),0);
    for(uint i=0;i<fieldvec.size();++i){
      field_abs[i]=fieldvec[i].abs();
    }
    return field_abs;
  }

  Vecfield grad_field(){
    auto absolute_values = abs_val();
    Vecfield grad_field;

    for(auto x=0; x<NX;++x){
      for(auto y=0; y<NY;++y){
	for(auto z=0; z<NZ;++z){
	  auto& elem = grad_field.val(x,y,z);
	  elem.x = deriv_scalar<Direction::x>(absolute_values,x,y,z);
	  elem.y = deriv_scalar<Direction::y>(absolute_values,x,y,z);
	  elem.z = deriv_scalar<Direction::z>(absolute_values,x,y,z);
	}
      }
    }
    return grad_field;
  }

  Vecfield rot_field(){
    Vecfield rot_field;
    for(auto x=0; x<NX;++x){
      for(auto y=0; y<NY;++y){
	for(auto z=0; z<NZ;++z){
	  auto& elem = rot_field.val(x,y,z);
	  elem.x = deriv_vf<Direction::z,Direction::y>(x,y,z) - deriv_vf<Direction::y,Direction::z>(x,y,z);
	  elem.y = deriv_vf<Direction::x,Direction::z>(x,y,z) - deriv_vf<Direction::z,Direction::x>(x,y,z);
	  elem.z = deriv_vf<Direction::y,Direction::x>(x,y,z) - deriv_vf<Direction::x,Direction::y>(x,y,z);
	}
      }
    }
    return rot_field;
  }

  template <Direction dir>
  void print_xy(std::string fn){
    std::ofstream ofs(fn);
    for(auto y=0;y<NY;++y){
      for(auto x=0;x<NX;++x){
	if(dir==Direction::x) ofs << val(x,y,z_0).x << " ";
	if(dir==Direction::y) ofs << val(x,y,z_0).y << " ";
	if(dir==Direction::z) ofs << val(x,y,z_0).z << " ";
      }
      ofs<<std::endl;
    }
  }

  void print_xy_absval(std::string fn){
    std::vector<double> absolute_values=abs_val();
    std::ofstream ofs(fn);
    for(auto y=0;y<NY;++y){
      for(auto x=0;x<NX;++x){
	ofs << absolute_values[lin_index(x,y,z_0)] << " ";
      }
      ofs<<std::endl;
    }
  }
};


//scalar product
double scal_prod(const Field& f1, const Field& f2){
  return f1.x * f2.x + f1.y * f2.y + f1.z * f2.z;
}

//single cross product
Field cross_prod(const Field& a, const Field& b){
  Field new_field;
  new_field.x = a.y*b.z - a.z*b.y;
  new_field.y = a.z*b.x - a.x*b.z;
  new_field.z = a.x*b.y - a.y*b.x;
  return new_field;
}

//Vecfield cross product
Vecfield cross_prod(const Vecfield& v1, const Vecfield& v2){
  Vecfield newfield;
  for(auto i=0;i<newfield.size();++i){
    newfield[i] = cross_prod(v1[i],v2[i]);
  }
  return newfield;
}



struct Bvars{
  Vecfield B;
  Vecfield grad_B;

  Vecfield b;
  Vecfield rot_b;

  Vecfield b_x_grad_B;

  Bvars(){
    B.generate_bfield();
    grad_B=B.grad_field();
    b=B.unit_field();
    rot_b=b.rot_field();

    b_x_grad_B=cross_prod(b,grad_B);
  }

  Field B_prime(double x, double y, double z, double v_par){
    return B.at(x,y,z) + rot_b.at(x,y,z)*(mass*v_par/charge);
  }

  Field B_prime(const Particle& pt, double v_par){
    return B.at(pt.x,pt.y,pt.z) + rot_b.at(pt.x,pt.y,pt.z)*(mass*v_par/charge);
  }

  Field B_prime(const Field& pt, double v_par){
    return B.at(pt.x,pt.y,pt.z) + rot_b.at(pt.x,pt.y,pt.z)*(mass*v_par/charge);
  }
};

//returns pair of v_par,SQU(v_perp)
double calc_v_par(const Field& b, const Particle& pt){
  return scal_prod(b,Field(pt.vx,pt.vy,pt.vz));
}

//returns pair of v_par,SQU(v_perp)
std::pair<double,double> calc_v_comp(const Field& b, const Particle& pt){
  double v_par = calc_v_par(b,pt);
  double v_perp_squ = SQU(pt.vx) + SQU(pt.vy) + SQU(pt.vz) - SQU(v_par);
  return std::make_pair(v_par,v_perp_squ);
}

double calc_mu(Bvars& bf, const Particle& pt, double v_perp_squ){
  return v_perp_squ*mass/(2*charge*bf.B.at(pt).abs());
}


//calculate d/dt (v_par) = f(v_par)
//TODO add efield later
//TODO improve speed by eliminating unnecessary calculations
double v_par_dot(Bvars& bf, const Particle& pt, double v_par, double mu){
  Field b = bf.b.at(pt);

  Field B_prime = bf.B_prime(pt,v_par);
  double B_prime_par_neg_inv = -1/scal_prod(B_prime,b);
  B_prime *= B_prime_par_neg_inv;

  Field grad_B = bf.grad_B.at(pt);
  grad_B *= mu;

  return scal_prod(B_prime,grad_B);
}


//calculate v_par_new via rk4
double v_par_new_rk4(Bvars& bf, const Particle& pt){
  Field b = bf.b.at(pt);
  auto vels=calc_v_comp(b,pt);
  double& v_par = vels.first;
  double& v_perp_squ = vels.second;
  //mu is constant of motion
  double mu = calc_mu(bf,pt,v_perp_squ);

  double k1 = v_par_dot(bf, pt, v_par, mu);
  double k2 = v_par_dot(bf, pt, v_par+0.5*timestep*k1, mu);
  double k3 = v_par_dot(bf, pt, v_par+0.5*timestep*k2, mu);
  double k4 = v_par_dot(bf, pt, v_par+step_size*k3, mu);

  k1 *= timestep/6;
  k2 *= timestep/3;
  k3 *= timestep/3;
  k4 *= timestep/6;

  return v_par + k1 + k2 + k3 + k4;
}


Field r_dot(Bvars& bf, const Field& start_pos, double v_par, double mu){
  Field b = bf.b.at(start_pos);
  
  Field B_prime = bf.B_prime(start_pos,v_par)*v_par;
  double B_prime_par_inv = 1 / scal_prod(B_prime,b);
  Field b_x_grad_B = bf.b_x_grad_B.at(start_pos) * (mu/charge);

  return (B_prime + b_x_grad_B)*B_prime_par_inv;
}

//calculate new particle position via rk4 and update particle
//TODO make this void once testing is done
Field r_new_rk4(Bvars& bf, Particle& pt){
  Field b = bf.b.at(pt);
  auto vels = calc_v_comp(b,pt);

  double& v_par = vels.first;
  double& v_perp_squ = vels.second;
  //mu is constant of motion
  double mu = calc_mu(bf,pt,v_perp_squ);

  double v_par_new = v_par_new_rk4(bf,pt);
  double v_par_med = 0.5 * (v_par_new + v_par);

  //std::cout<<"v_perp = "<<sqrt(v_perp_squ)<<"  ";
  //std::cout<<"mu = "<<mu<<"  ";
  //std::cout<<std::endl;

  Field start_pos(pt.x,pt.y,pt.z);

  Field k1 = r_dot(bf, start_pos, v_par, mu);
  Field k2 = r_dot(bf, start_pos + k1*(0.5*timestep), v_par_med, mu);
  Field k3 = r_dot(bf, start_pos + k2*(0.5*timestep), v_par_med, mu);
  Field k4 = r_dot(bf, start_pos + k3*timestep, v_par_new, mu);

  k1 *= timestep/6;
  k2 *= timestep/3;
  k3 *= timestep/3;
  k4 *= timestep/6;

  start_pos += k1 + k2 + k3 + k4;

  //FIXME delete this later!
  start_pos.z = pt.z;

  pt.vx = (start_pos.x - pt.x) / timestep;
  pt.vy = (start_pos.y - pt.y) / timestep;
  pt.vz = (start_pos.z - pt.z) / timestep;

  pt.x = start_pos.x;
  pt.y = start_pos.y;
  pt.z = start_pos.z;

  return start_pos;
}

void print_pt_data(const std::vector<Field>& trajectory, const std::vector<double>& energies, std::string fn){
  if(trajectory.size() == energies.size() && !trajectory.empty()){
    std::ofstream ofs(fn);
    ofs<<"#x y z energy"<<std::endl;
    for(uint i=0;i<trajectory.size();++i){
      const Field& f = trajectory[i];
      ofs << f.x << " " << f.y << " " << f.z << " " << energies[i] << std::endl;
    }
  }else{
    std::cout<<"No matching trajectory data to print." << std::endl;
    exit(-1);
  }
}

int main(int argc, char** argv){
  Bvars bf;
#if 0
  //test of CIC interpolation
  Field x1=bf.B.at(x_0,y_0+5,z_0);
  Field x2=bf.B.at(x_0+2.5,y_0+5,z_0+2.5);
  Field x3=bf.B.at(x_0+5,y_0+5,z_0+5);
  std::cout<<"x1.x="<<x1.x<<"  x2.x="<<x2.x<<"  x3.x="<<x3.x<<std::endl;
  std::cout<<"x1.y="<<x1.y<<"  x2.y="<<x2.y<<"  x3.y="<<x3.y<<std::endl;
  std::cout<<"x1.z="<<x1.z<<"  x2.z="<<x2.z<<"  x3.z="<<x3.z<<std::endl;
  std::cout<<"x1.abs="<<x1.abs()<<"  x2.abs="<<x2.abs()<<"  x3.abs="<<x3.abs()<<std::endl;
#endif

#if 1
  //test of bfield components
  bf.B.print_xy<Direction::x>("Bx.dat");
  bf.B.print_xy<Direction::y>("By.dat");
  bf.B.print_xy_absval("B_abs.dat");

  bf.b.print_xy_absval("b_abs.dat");

  bf.grad_B.print_xy<Direction::x>("gradBx.dat");
  bf.grad_B.print_xy<Direction::y>("gradBy.dat");
  bf.grad_B.print_xy_absval("gradB_abs.dat");

  bf.rot_b.print_xy<Direction::x>("rotBx.dat");
  bf.rot_b.print_xy<Direction::y>("rotBy.dat");
  bf.rot_b.print_xy<Direction::z>("rotBz.dat");
  bf.rot_b.print_xy_absval("rotB_abs.dat");
#endif

  bf.rot_b.print_xy<Direction::z>("rot_bz.dat");
  bf.rot_b.print_xy_absval("rot_b_abs.dat");

#if 1
  //test of solver
  Particle pt(x_0+10,y_0,z_0,0,1,0);
  double e0 = pt.energy();
  pt.print();


  //627 steps for full circle at x_0+10 and vy=1
  int period=627;
  int num_periods=10;
  int maxstep = period*num_periods;

  std::vector<Field> trajectories(maxstep);
  std::vector<double> energies(maxstep);
  for(auto i=0;i<maxstep;++i){
    trajectories[i] = r_new_rk4(bf,pt);
    double enew = pt.energy();
    energies[i] = fabs(e0-enew)/e0*100;
    //std::cout<<"step "<<i<<":  ";
    //pt.print();
  }

  print_pt_data(trajectories,energies,"trajectory.dat");
#endif

  return 0;
}

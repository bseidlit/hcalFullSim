//CosmicSpray class
// Author: Daniel Lis
// Brief: Particel generator Class that sources a muon with a vertex and momentum that should mimic real life
// modified by Shuhang on 03/2022: used a new random energy/theta generation algorithm as well as a new vertex generation method. The result agrees with the real data and offers a stright forward muon rate calculation. 

#include "CosmicSpray.h"

#include <g4main/PHG4InEvent.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4Utils.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <iostream>
#include <TSystem.h>
#include "TROOT.h"
#include "TF3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type

#include <vector>  // for vector, vector<>::const_iterator

// Declarations
// fixes y plane of cosmic spray here
double CosmicSpray::_y_fix;
// max and min for x spray plane geoemtry
double CosmicSpray::_x_max;
double CosmicSpray::_x_min;
// max and min for z spray plane geometry
double CosmicSpray::_z_max;
double CosmicSpray::_z_min;

std::string CosmicSpray::_detector_name;
// parameters to check if the path wen thtrough the set_detector_dimensions
// max and min for x in detector geoemtry
double CosmicSpray::_x_max_det;
double CosmicSpray::_x_min_det;
// max and min for z in detector direction
double CosmicSpray::_z_max_det;
double CosmicSpray::_z_min_det;
double CosmicSpray::_y_det;
// physics params
double CosmicSpray::_R_earth;
double CosmicSpray::_d_earth;
// gun energy
double CosmicSpray::_gun_e;

// offset of spray plane
double CosmicSpray::_offx;
double CosmicSpray::_offz;

double CosmicSpray::_Emin;
double CosmicSpray::_Emax;
double CosmicSpray::_thetamin;
double CosmicSpray::_thetamax;
int CosmicSpray::_E_bin;
int CosmicSpray::_theta_bin;
int CosmicSpray::_nbin;
double CosmicSpray::_length;
double CosmicSpray::cdf[1000001];

class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;

// Calculate the gamma factor
double CosmicSpray::GetGamma0(double &energy){
  double mass = 105; // MeV
  double gamma0 = sqrt(1 + TMath::Power(energy/(mass),2));
  return gamma0;
}

// calculate beta
double CosmicSpray::GetBeta(double &gamma){
  double b = sqrt(1-TMath::Power((1/gamma),2));
  return b;
}

double CosmicSpray::EnergyAngularDistributionnew(double en, double theta){
  double energy = TMath::Power(10, en);
  double y0 = 1000; 
  double rho = 0.76;
  //double Alpha = 0.002382;
  double lambda = 120;
  double kappa = 2.645;
  double bp = 0.771;
  double jp = 148.16;
  double Bm = 1.041231831;
  double alpha_ion = 0.0025;
  
  double Ep = (energy + alpha_ion * y0 *(1/cos(theta) - 0.1) ) / rho;

  double Pm =TMath::Power( (0.1 * cos(theta) * (1 - (alpha_ion *(y0 / cos(theta) - 100 ) / (rho * Ep)))), (Bm / ((rho * Ep + 100 * alpha_ion) * cos(theta) ) ) );

  //double p = Alpha * Pm * (TMath::Power(Ep,-kappa)) * lambda * bp * jp / (Ep * cos(theta) * bp * jp);
  double p = 1000 * Pm * (TMath::Power(Ep,-kappa)) * lambda * bp * jp / (Ep * cos(theta) * bp * jp);
  //std::cout<<"energy = "<<energy<<", theta = "<<theta<<", Ep = "<<Ep<<", Pm"<<Pm<<" ,p = "<<p<<std::endl;
  return p;
  
}

double CosmicSpray::EnergyAngularDistribution(double energy, double theta){
  //if (par[0]==21312) std::cout << "Help";
  double en = TMath::Power(10, energy);
  double th = theta;
  double p;
  double R = 600.;
  double c = 29979245800; //cm/s
  double tau = 0.0000022; // seconds
  double h_0 = 1000000; //cm
  double g = GetGamma0(en);
  double b = GetBeta(g);
  double ld = c*tau*g*b;
  double norm = 1/((sqrt(R*R*cos(th)*cos(th) + 2*R + 1) - R*cos(th))*ld);
  double f = exp(-h_0*(sqrt(R*R*cos(th)*cos(th) + 2*R + 1) - R*cos(th))/(ld));
//std::cout<<"g: "<<g<<"b: "<<b<<"ld: "<<ld<<"norm: "<<norm<<"f: "<<f<<std::endl;
  p = 100000*fabs(norm*f);
  return p;
}

CosmicSpray::CosmicSpray(const std::string &name = "COSMICS", const std::string &detector = "FULL", const int &debug = 0)
  : PHG4ParticleGeneratorBase(name)
{
  _x_max = 264.71;
  _x_min = 183.3;
  _z_max = 304.91;
  _z_min = -304.91;
  _x_max_det = 264.71;
  _x_min_det = 183.3;
  _z_max_det = 304.91;
  _z_min_det = -304.91;
  _y_det = 0.;
  _y_fix = _x_max;
  _R_earth = 637100000;
  _d_earth = 1000000;
  _offx = 100.;
  _offz = 100.;
  _E_bin = 1000;
  _theta_bin = 1000;
  _Emin = TMath::Log10(0.5)+3;
  _Emax = 2.5+3;
  _thetamin = 0;
  _thetamax = TMath::Pi()/2;
  _length = 20;
  _nbin = _E_bin * _theta_bin +1;
  int ibin = 0;
  double integral = 0;
  //double cdf[nbin] = {0};
  //doing integral here
  static const int nbin = _nbin;
  for(int i = 0; i < _E_bin; i++){
    for(int j = 0; j < _theta_bin; j++){
      double t = j  * _thetamax/_theta_bin;
      double e = _Emin + i *(_Emax - _Emin)/_E_bin;
      
      integral += CosmicSpray::EnergyAngularDistribution(e,t);
      //std::cout<<integral<<"e:"<<e<<" t:"<<t<<std::endl;
      ibin++;
      cdf[ibin] = integral;
    }
  }
  // std::cout<<"integal is"<<integral<<std::endl;  
  //normalization
  if(integral <= 0) std::cout<<"wrong prob distribution"<<std::endl;
  for(int i = 0; i< nbin;i++){
    cdf[i] /= integral;
    // std::cout<<cdf[i]<<std::endl;
  }

  if (detector == "HCALSECTOR"){
    _detector_name = detector;
    _x_max = 264.71 + _offx;
    _x_min = 183.3 + _offx;
    _z_max = 304.91 + _offz;
    _z_min = -304.91 + _offz;
    _x_max_det = 264.71;
    _x_min_det = 183.3;
    _z_max_det = 304.91;
    _z_min_det = -304.91;
    _y_fix = 200;
  }
  else if(detector == "EMCALSECTOR"){
    _detector_name = detector;
    _x_max = 264.71 + _offx;
    _x_min = 183.3 + _offx;
    _z_max = 304.91 + _offz;
    _z_min = -304.91 + _offz;
    _x_max_det = 264.71;
    _x_min_det = 183.3;
    _z_max_det = 304.91;
    _z_min_det = -304.91;
    _y_fix = 200;
  }
  else {
    _detector_name = "FULLNEW";
  }

  //vertex function

  //_f_vertex = new TF1("_f_vertex","(1/TMath::Power(TMath::Pi(),2))*((TMath::Pi()/2)-fabs(x)+TMath::Pi()/2-sin(fabs(x)+TMath::Pi()/2)*cos(fabs(x)+TMath::Pi()/2))",-1*TMath::Pi(), TMath::Pi());

  // momentum function
  // _f_momentum = new TF2("_f_momentum",CosmicSpray::EnergyAngularDistribution, 0, TMath::Pi()/2,2, 5);

  _debug = debug;
  return;
}
 // add a particle to the generator
void CosmicSpray::add_particle(const std::string &name, const unsigned int num)
{
  _particle_names.clear();
  _particle_names.push_back(std::make_pair(name, num));
  return;
}

// edit the plane the vertex is grabbed from
void CosmicSpray::edit_cosmics_plane(double height, double x_max, double x_min, double z_min, double z_max){
  _y_fix = height;
  _x_max = x_max;
  _x_min = x_min;
  _z_max = z_max;
  _z_min = z_min;
  return;
}

void CosmicSpray::set_detector_dimensions(double x_max, double x_min, double z_min, double z_max){
  _x_max_det = x_max;
  _x_min_det = x_min;
  _z_max_det = z_max;
  _z_min_det = z_min;
  return;
}

void CosmicSpray::set_sprayplane_offset(double x, double z){
  _x_max = _x_max - _offx;
  _x_min =  _x_min - _offx;
  _z_max = _z_max - _offz;
  _z_min = _z_min - _offz;
  _offx = x;
  _offz = z;
  _x_max = _x_max + _offx;
  _x_min = _x_min + _offx;
  _z_max = _z_max + _offz;
  _z_min = _z_min + _offz;
  return;
}

int CosmicSpray::process_event(PHCompositeNode *topNode)
{
  // set_vertex
  if(_debug) std::cout<<"Processing Event"<<std::endl;
  std::string pdgname;
  for (unsigned int i = 0; i < _particle_names.size(); ++i){
    pdgname = _particle_names[i].first;
    if(_debug) std::cout<<"Particle added: "<<pdgname << std::endl;
  }
  gRandom->SetSeed(0);
  int pdgcode = get_pdgcode(pdgname);
  int trackid = 0;
  double gun_t = 0.0;
  double gun_x =0, gun_y =0, gun_z = 0;
  double gun_px = 0, gun_py = 0, gun_pz = 0;
  // bool GoodEvent = true;
  TRandom *randomGen = new TRandom();
  randomGen->SetSeed();
  
  if(true){
    double x1, y1, z1;
    
    x1 = randomGen->Uniform(-1*_x_max, _x_max);
    y1 = randomGen->Uniform(-1*_y_fix, _y_fix);
    z1 = randomGen->Uniform(-1*_z_max, _z_max);
    Double_t r1 = randomGen->Uniform(0, 1);
    //ramdomly generate momentum
    
    int inbin = TMath::BinarySearch(_nbin,cdf,r1);
    int bine = inbin / _theta_bin;
    int binth = inbin - bine*_theta_bin;
    double er = _Emin + (_Emax - _Emin)/_E_bin * bine;
    double tr = _thetamin + (_thetamax - _thetamin)/_theta_bin * binth;
    tr = randomGen->Uniform(tr, tr + (_thetamax - _thetamin)/_theta_bin );
    double pr;
    //std::cout<<"ibin = "<<inbin<<" bin1 = "<<bine<<" bin2 = "<<binth<<"r1 = "<<r1<<std::endl;
    //er = 4;
    _gun_e = TMath::Power(10, er-3);
    // _gun_e = er;
    pr = randomGen->Uniform(0,  2* TMath::Pi());

    gun_px = _gun_e*sin(tr)*sin(pr);
    gun_py = -1*fabs(_gun_e*cos(tr));
    gun_pz = _gun_e*sin(tr)*cos(pr);
    
    TVector3 gun_p(gun_px, gun_py, gun_pz);
    TVector3 dir = gun_p.Unit();

    gun_x = x1 - _length * dir.X();
    gun_y = y1 - _length * dir.Y();
    gun_z = z1 - _length * dir.Z();
    while((gun_y < _y_fix)&&(abs(gun_x)<_x_max)&&(abs(gun_z)<_z_max)){
      gun_x += - _length * dir.X();
      gun_y += - _length * dir.Y();
      gun_z += - _length * dir.Z();
    }
    // std::cout<<x1<<" "<<y1<<" "<<z1<<std::endl;
    //std::cout<<dir.X()<<" "<<dir.Y()<<" "<<dir.Z()<<std::endl;
    //std::cout<<"x = "<<gun_x<<" gun_y = "<<gun_y<<" gun_z = "<<gun_z<<std::endl;
  }
  /***
  if (_detector_name == "FULL"){
    gun_z = randomGen->Uniform(-1*_z_max, _z_max);
    double theta_ran;
    theta_ran = _f_vertex->GetRandom();
    gun_x = _x_max * sin(theta_ran);
    gun_y = _x_max * cos(theta_ran);
    if(_debug) std::cout<<"Vertex: "<<gun_x<<" / " << gun_y <<" / " << gun_z << std::endl;

    double er, tr,t_lim, pr, p_max;
    t_lim = 0;
    if (theta_ran > TMath::Pi()/2){
      t_lim = fabs(theta_ran) - TMath::Pi()/2;
      _f_momentum->SetRange(t_lim, TMath::Pi()/2, 2,5);
    }
    else if (theta_ran < -1*TMath::Pi()/2){
      t_lim = fabs(theta_ran) - TMath::Pi()/2;
      _f_momentum->SetRange(t_lim, TMath::Pi()/2, 2,5);
    }
    else {
      t_lim = TMath::Pi()/2 - fabs(theta_ran);
      _f_momentum->SetRange(0, TMath::Pi()/2, 2,5);
    }
    if (_debug) std::cout << "about to f_momentum" << std::endl;
    //_f_momentum->GetRandom2(tr, er);
    er = 3.5;
    tr = TMath::Pi()/2 -1.3;
    _gun_e = TMath::Power(10, er - 3);
    if (_debug) std::cout << "just f_momentum" << std::endl;

    if (theta_ran > TMath::Pi()/2){
      p_max = TMath::ACos(sin(t_lim)/sin(tr));
    }
    else if (theta_ran < -1*TMath::Pi()/2){
      p_max = TMath::ACos(sin(t_lim)/sin(tr));
    }
    else {
      if (tr < t_lim) p_max = TMath::Pi();
      else p_max = TMath::Pi()/2 + TMath::ASin(sin(t_lim)/sin(tr));
    }
    if (_debug) std::cout << "about to pr" << std::endl;
    pr = (randomGen->Uniform() - 0.5 )*2*p_max;

    if (theta_ran > 0){
      if (pr > 0) pr = pr + TMath::Pi() - p_max;
      else if (pr < 0) pr = pr - TMath::Pi() + p_max;
    }

    gun_px = _gun_e*sin(tr)*sin(pr);
    gun_py = -1*fabs(_gun_e*cos(tr));
    gun_pz = _gun_e*sin(tr)*cos(pr);
    if (_debug) std::cout << " gun_px=" << gun_px << std::endl;
  }

  else {
    _f_momentum->SetRange(0, TMath::Pi()/2, 2, 5);
    while(!GoodEvent){
      gun_z = randomGen->Uniform(_z_min, _z_max);
      gun_x = randomGen->Uniform(_x_min, _x_max);
      gun_y = _y_fix;
      double tr, er;
      _f_momentum->GetRandom2(tr, er);
      _gun_e = TMath::Power(10, er - 3);
      double pr = randomGen->Uniform (0, 2*TMath::Pi());
      gun_py = _gun_e*-1*abs(cos(tr));
      gun_px = _gun_e*abs(sin(tr))*sin(pr);
      gun_pz = _gun_e*abs(sin(tr))*cos(pr);

      double x = gun_x + (_y_fix - _y_det)*tan(tr)*sin(pr);
      double z = gun_z + (_y_fix - _y_det)*tan(tr)*cos(pr);
      if (x < _x_max_det && x > _x_min_det && z < _z_max_det && z > _z_min_det){
        GoodEvent = true;
      }
    }
  }
  ***/
 

  if(_debug) std::cout<<"Momentum: "<<gun_px<<" / "<<gun_py<<" / "<<gun_pz<<std::endl;
  if(_debug)std::cout<<"total mom: "<<_gun_e<<std::endl;
  if(_debug)std::cout<<"Before adding vertex"<<std::endl;
  _InEvent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  int vtxindex = _InEvent->AddVtx(gun_x, gun_y ,gun_z,gun_t);
  if(_debug)std::cout<<"After adding vertex"<<std::endl;

  PHG4Particle *particle = new PHG4Particlev2();
  particle->set_track_id(trackid);
  if(_debug)std::cout<<"track_id: "<<trackid<<std::endl;
  particle->set_vtx_id(vtxindex);
  if(_debug)std::cout<<"vtxindex: "<<vtxindex<<std::endl;
  particle->set_parent_id(0);
  if(_debug)std::cout<<"parent_id "<<std::endl;
  particle->set_name(pdgname);
  if(_debug)std::cout<<"pdgname: "<<pdgname<<std::endl;
  particle->set_pid(pdgcode);
  if(_debug)std::cout<<"pdgcode: "<<pdgcode<<std::endl;
  particle->set_px(gun_px);
  if(_debug)std::cout<<"px"<<std::endl;
  particle->set_py(gun_py);
  if(_debug)std::cout<<"py"<<std::endl;
  particle->set_pz(gun_pz);
  if(_debug)std::cout<<"pz"<<std::endl;
  particle->set_e(_gun_e);
  if(_debug)std::cout<<"ene"<<std::endl;

  _InEvent->AddParticle(vtxindex, particle);
  if(_debug)std::cout<<"left COSMICS"<<std::endl;

  return 0;
}

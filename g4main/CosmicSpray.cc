//CosmicSpray class
// Author: Daniel Lis
// Brief: Particel generator Class that sources a muon with a vertex and momentum that should mimic real life
// modified by Shuhang on 03/2022: now this class serves as a wrapper class that drives "EcoMug" 
#include "CosmicSpray.h"
#include "EcoMug.h"
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



class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;





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

  if (detector == "HCALSECTOR"){
    _detector_name = detector;
    // _x_max = 264.71 + _offx;
    // _x_min = 183.3 + _offx;
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

  gen.SetUseHSphere(); // half-spherical surface generation
  gen.SetHSphereRadius(5.); // half-sphere radius
  gen.SetHSphereCenterPosition({{0., 0., -2.6471}});
  gen.SetMinimumMomentum(0.5);
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


int CosmicSpray::process_event(PHCompositeNode *topNode)
{
  // set_vertex
  if(_debug) std::cout<<"Processing Event"<<std::endl;
  std::string pdgname = "mu-";
  int pdgcode  = 13;
  int trackid = 0;
  double gun_t = 0.0;
  double gun_x =0, gun_y =0, gun_z = 0;
  double gun_px = 0, gun_py = 0, gun_pz = 0;
  // bool GoodEvent = true;
  gen.Generate();
  std::array<double, 3> muon_position = gen.GetGenerationPosition();
  double muon_p = gen.GetGenerationMomentum();
  _gun_e = sqrt(0.105658* 0.105658 + muon_p * muon_p);
  double tr = gen.GetGenerationTheta();
  double pr = gen.GetGenerationPhi();  
  double muon_charge = gen.GetCharge();

  if(muon_charge >0){
    pdgcode = -13;
    pdgname = "mu+";
  }

  gun_px = muon_p*sin(tr)*sin(pr);
  gun_py = -1*fabs(muon_p*cos(tr));
  gun_pz = muon_p*sin(tr)*cos(pr);
  

  gun_x = muon_position[1] * 100;
  gun_y = muon_position[2] * 100;
  gun_z = muon_position[0] * 100;
  

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

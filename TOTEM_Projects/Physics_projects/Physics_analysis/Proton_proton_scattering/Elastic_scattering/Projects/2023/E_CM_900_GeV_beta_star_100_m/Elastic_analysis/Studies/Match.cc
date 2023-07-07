#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include "math.h"
#include "math.h"

using namespace std;

#include "stdlib.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TFitResult.h"

#include "/afs/cern.ch/work/f/fnemes/main_workspace/2023/Elastic_analysis/E_CM_900_GeV_beta_star_100_m/CMSSW_10_6_17/src/ElasticAnalysis/ElasticAnalyzer/interface/ElasticAnalyzer.h"

#include <iomanip>

const double Ly_near = 19.1 ;
const double Ly_far  = 17.4 ; 

const double Lx_near = 8.8228247 ;
const double Lx_far  = 5.8873862 ;

const double vy_near = 2.47 ;
const double vy_far  = 3.13 ;

const double vx_near = -3.3797 ;
const double vx_far  = -3.0494 ;

const string RP_near_name = "XRPV.D6R5.B1" ;
const string RP_far_name  = "XRPV.B6R5.B1" ;

const string RP_near_hor_name = "XRPH.D6R5.B1" ;
const string RP_far_hor_name  = "XRPH.B6R5.B1" ;

const double n_sigma = 3.0 ;

Double_t my_exponential_distribution(Double_t *x, Double_t *par)
{
      return par[0] * exp(par[1] * (x[0])) ;
}

Double_t theta_star_rad_from_t_GeV2(double t_GeV2, double beam_momentum_GeV)
{
        return asin(sqrt(-t_GeV2) / beam_momentum_GeV) ;
}

const double PI_TIMES_2 = (2.0 * TMath::Pi()) ;
const double beam_momentum_GeV = 450 ;

const double Constant = 2.363e7 ;
const double Slope = -14.0 ;
//const double Slope = -1.69642e+01 ;
// const double Slope = -1.69642 ;

// const string aperture_filename="/afs/cern.ch/eng/lhc/optics/V6.503/aperture/aperture.b1.madx" ;
const string aperture_filename="data/lhc/apertures_after_IP5_B1.txt" ;
const string optics_file_name = "/afs/cern.ch/work/f/fnemes/main_workspace_git/projects/LHCRealOpticsTools/OpticsProject/Project_files/2018/E_CM_900_GeV/Valentina/twiss_beam1_900_beta11_thetay0.tfs" ;

class TAperture
{ 
  public:

  string name ;
  TH2D *histogram ;
  TH2D *histogram_vertex ;
  TH2D *histogram_x_y_mm ;
  TH1D *histogram_t ;
  TH1D *histogram_phi ;

  double rectangle_half_width ;
  double rectangle_half_height ;
  double ellipse_semi_axis_hor ;
  double ellipse_semi_axis_ver ;
  double edge ;
  double hedge ;
  
  double Lx, Ly, vx, vy ;
  
  TAperture(string name, double Lx, double Ly, double vx, double vy, double rectangle_half_width, double rectangle_half_height, double ellipse_semi_axis_hor, double ellipse_semi_axis_ver, double edge, double hedge):
  name(name), Lx(Lx), Ly(Ly), vx(vx), vy(vy), rectangle_half_width(rectangle_half_width), rectangle_half_height(rectangle_half_height), ellipse_semi_axis_hor(ellipse_semi_axis_hor), ellipse_semi_axis_ver(ellipse_semi_axis_ver), edge(edge), hedge(hedge)
  {
    histogram = new TH2D(name.c_str(), name.c_str(), 100, -1000e-6, 1000e-6, 100, -1000e-6, 1000e-6) ;
    histogram_vertex = new TH2D((name + "_vertex").c_str(), (name + "_vertex").c_str(), 100, -1000e-6, 1000e-6, 100, -1000e-6, 1000e-6) ;
    histogram_t = new TH1D((name + "_t").c_str(), (name + "_t").c_str(), 1000, 0, 10) ;
    histogram_phi = new TH1D((name + "_phi").c_str(), (name + "_phi").c_str(), 100, 0, 10) ;

    histogram_x_y_mm = new TH2D("histogram_x_y_mm", "histogram_x_y_mm", 100, -1e-2, 1e-2,  100, -1e-2, 1e-2) ;
  }
  
  int getrpDecId(double) ;
} ;

int TAperture::getrpDecId(double y)
{
  int value = 0 ;

  if(name.compare("XRPH.D6R5.B1") == 0) value = 103 ;
  else if(name.compare("XRPH.B6R5.B1") == 0) value = 123 ;
  else if(name.compare("XRPV.D6R5.B1") == 0)
  {
    value = (y > 0) ? 104 : 105 ;
  }
  else if(name.compare("XRPV.B6R5.B1") == 0)
  {
    value = (y > 0) ? 124 : 125 ;
  }
  
  // if(value != 0) cout << "valuetest " << value << endl ;

  return value ;
}

TF1 *t_GeV2_distribution = NULL ;
TRandom3 myrand ;

class TProton
{
  public:

  double theta_x_star ;
  double theta_y_star ;
  double x_star ;
  double y_star ;

  double theta_x_star_pert ;
  double theta_y_star_pert ;

  double theta_x_star_reco ;
  double theta_y_star_reco ;
  double x_star_reco ;
  double y_star_reco ;

  TProton() {} ;
  void Reconstruct() ;
} ;


void TProton::Reconstruct()
{
  double x_near = (Lx_near * theta_x_star_pert) + (vx_near * x_star) ;
  double x_far  = (Lx_far *  theta_x_star_pert) + (vx_far *  x_star) ;
  double y_near = (Ly_near * theta_y_star_pert) + (vy_near * y_star) ;
  double y_far  = (Ly_far *  theta_y_star_pert) + (vy_far *  y_star) ;

  double determinant_x = ((Lx_near * vx_far) - (Lx_far * vx_near)) ;
  double determinant_y = ((Ly_near * vy_far) - (Ly_far * vy_near)) ;

  theta_x_star_reco = ((x_near * vx_far) - (x_far * vx_near)) / determinant_x ;
  x_star_reco =      -((x_near * Lx_far) - (x_far * Lx_near)) / determinant_x ;

  theta_y_star_reco = ((y_near * vy_far) - (y_far * vy_near)) / determinant_y ;
  y_star_reco =      -((y_near * Ly_far) - (y_far * Ly_near)) / determinant_y ;

  double theta_star_reco = sqrt((theta_x_star_reco*theta_x_star_reco) + (theta_y_star_reco*theta_y_star_reco)) ;
  double minus_t_GeV2_reco = - beam_momentum_GeV * beam_momentum_GeV * theta_star_reco * theta_star_reco ;
}

class TProtonPair
{
  public:

  double minus_t_GeV2 ;
  Double_t phi_IP5_rad ;    

  double theta_star_rad ;

  TProton p_b1, p_b2 ;

  double beam_divergence_x_b1 ;
  double beam_divergence_y_b1 ;

  double beam_divergence_x_b2 ;
  double beam_divergence_y_b2 ;
  
  bool reconstructed ;

  double theta_x_star_reco ;
  double theta_y_star_reco ;
  double x_star_reco ;
  double y_star_reco ;

  TProtonPair() ;
  void GenerateParticles() ;
  bool TestApertures(vector<TAperture *> &) ;
  bool TestAperturesOneProton(TProton, vector<TAperture *> &, bool) ;
  
  virtual void TestDetectorPair(map<unsigned int, RP_struct_type>::iterator, map<unsigned int, RP_struct_type>::iterator, unsigned int, unsigned int) ;  
} ;

const double vertex_size_m = 200e-6 ;

const double beam_divergence_x_rad = 50e-6 ;
const double beam_divergence_y_rad = 27e-6 ;

TProtonPair::TProtonPair()
{
  GenerateParticles() ;
}

void TProtonPair::GenerateParticles()
{
  reconstructed = false ;

  minus_t_GeV2 = t_GeV2_distribution->GetRandom() ;
  phi_IP5_rad = gRandom->Uniform(0, PI_TIMES_2) ;    

  theta_star_rad = theta_star_rad_from_t_GeV2(-minus_t_GeV2, beam_momentum_GeV) ;

  p_b1.theta_x_star =  (theta_star_rad * cos(phi_IP5_rad)) ;
  p_b1.theta_y_star =  (theta_star_rad * sin(phi_IP5_rad)) ;

  p_b2.theta_x_star = -(theta_star_rad * cos(phi_IP5_rad)) ;
  p_b2.theta_y_star = -(theta_star_rad * sin(phi_IP5_rad)) ;

  double factor = 1.0 ;
  double factor_vtx = 1.0 ;

  p_b1.x_star = factor_vtx * vertex_size_m * myrand.Gaus() ;
  p_b1.y_star = factor_vtx * vertex_size_m * myrand.Gaus() ;

  p_b2.x_star = p_b1.x_star ;
  p_b2.y_star = p_b1.y_star ;
  
  beam_divergence_x_b1 = factor * beam_divergence_x_rad * myrand.Gaus() ;
  beam_divergence_y_b1 = factor * beam_divergence_y_rad * myrand.Gaus() ;

  beam_divergence_x_b2 = factor * beam_divergence_x_rad * myrand.Gaus() ;
  beam_divergence_y_b2 = factor * beam_divergence_y_rad * myrand.Gaus() ;

  p_b1.theta_x_star_pert = p_b1.theta_x_star + beam_divergence_x_b1 ;
  p_b1.theta_y_star_pert = p_b1.theta_y_star + beam_divergence_y_b1 ;

  p_b2.theta_x_star_pert = p_b2.theta_x_star + beam_divergence_x_b2 ;
  p_b2.theta_y_star_pert = p_b2.theta_y_star + beam_divergence_y_b2 ;
  
  theta_x_star_reco = 0 ;
  theta_y_star_reco = 0 ;
  x_star_reco = 0 ;
  y_star_reco = 0 ;
  
}

bool TProtonPair::TestApertures(vector<TAperture *> &vector_apertures) 
{
  bool test_p_b1 = TestAperturesOneProton(p_b1, vector_apertures, true) ;
  bool test_p_b2 = TestAperturesOneProton(p_b2, vector_apertures, false) ;
  
  if(test_p_b1 && test_p_b2)
  {
    p_b1.Reconstruct() ;
    p_b2.Reconstruct() ;

    theta_x_star_reco = (p_b1.theta_x_star_reco - p_b2.theta_x_star_reco) / 2.0 ;
    theta_y_star_reco = (p_b1.theta_y_star_reco - p_b2.theta_y_star_reco) / 2.0 ;

    x_star_reco = (p_b1.x_star_reco + p_b2.x_star_reco) / 2.0 ;
    y_star_reco = (p_b1.y_star_reco + p_b2.y_star_reco) / 2.0 ;
    
    reconstructed = true ;
  }

}

const double dx_threshold_between_vertical_and_horizontal_mm = 0.2 ;

map<string, TH2D *> map_of_THorizontal_and_vertical_xy_histogram ;

void TProtonPair::TestDetectorPair(map<unsigned int, RP_struct_type>::iterator it1, map<unsigned int, RP_struct_type>::iterator it2, unsigned int detector_1, unsigned int detector_2)
{
  if((it1->first == detector_1) && (it2->first == detector_2))
  {
    stringstream ss_1 ;
    stringstream ss_2 ;

    ss_1 << it1->first ;
    ss_2 << it2->first ;

    string key_for_coords = ss_1.str() + "_" + ss_2.str() ;
    string key_for_coords_1 = ss_1.str() + "_" + ss_2.str() + "_" + ss_1.str() ;
    string key_for_coords_2 = ss_1.str() + "_" + ss_2.str() + "_" + ss_2.str()  ;
    string key_for_coords_3 = ss_1.str() + "_" + ss_2.str() + "_corr_x_dx"  ;
    string key_for_coords_4 = ss_1.str() + "_" + ss_2.str() + "_corr_y_dy"  ;

    string name_x = "dx_" + ss_1.str() + "_" + ss_2.str() ;
    string name_y = "dy_" + ss_1.str() + "_" + ss_2.str() ;
    
    // cout << name_x << " " << name_y << endl ;
    
    const bool idealized_case = true ; // comparing the exact same coordinate with some offsets

    if((fabs(it2->second.x - it1->second.x) < dx_threshold_between_vertical_and_horizontal_mm) && (fabs(it2->second.y - it1->second.y) < dx_threshold_between_vertical_and_horizontal_mm))
    {
    
      double pert_alpha_rad = 0.0e-2 ;
      double pert_a_mm = 2.0 ;
      double pert_b_mm = 0.0 ;

      double pert_hor_x = 0.0 ;
      double pert_hor_y = 0.0 ;

      if(idealized_case)
      {
        pert_hor_x = (cos(pert_alpha_rad) * it2->second.x) - (sin(pert_alpha_rad) * it2->second.y) + pert_a_mm ;
        pert_hor_y = (sin(pert_alpha_rad) * it2->second.x) + (cos(pert_alpha_rad) * it2->second.y) + pert_b_mm ;
      }
      else
      {
        pert_hor_x = (cos(pert_alpha_rad) * it1->second.x) - (sin(pert_alpha_rad) * it1->second.y) + pert_a_mm ;
        pert_hor_y = (sin(pert_alpha_rad) * it1->second.x) + (cos(pert_alpha_rad) * it1->second.y) + pert_b_mm ;
      }

      map_of_THorizontal_and_vertical_xy_pairs_to_match[key_for_coords].push_back(new THorizontal_and_vertical_xy_pairs_to_match(pert_hor_x, pert_hor_y, it2->second.x, it2->second.y)) ;

      // cout << "hello " << key_for_coords_3 << endl ; 
      map_of_THorizontal_and_vertical_xy_histogram[key_for_coords_1]->Fill(pert_hor_x, pert_hor_y) ;
      map_of_THorizontal_and_vertical_xy_histogram[key_for_coords_2]->Fill(it2->second.x, it2->second.y) ;

      map_of_THorizontal_and_vertical_xy_histogram[key_for_coords_3]->Fill(it1->second.x, it2->second.x - it1->second.x) ;
      map_of_THorizontal_and_vertical_xy_histogram[key_for_coords_4]->Fill(it1->second.y, it2->second.y - it1->second.y) ;
      
      // if(it2->second.x < -4) cout << "to_be_saved " << key_for_coords << " " <<  it1->second.x << " " <<  it1->second.y << " " <<  it2->second.x << " " <<  it2->second.y << " " <<  endl ;
      // cout << " dx: " << (it2->second.x - it1->second.x) << " dy: " << (it2->second.y - it1->second.y) << endl ;

      string name_x2 = "xy_" + ss_1.str() + "_if_" + ss_1.str() + "_" + ss_2.str() ;
      string name_y2 = "xy_" + ss_2.str() + "_if_" + ss_1.str() + "_" + ss_2.str() ;

      // cout << name_x2 << " " << name_y2 << endl ;

    }
  }
}

bool TProtonPair::TestAperturesOneProton(TProton proton, vector<TAperture *> &vector_apertures, bool save_result) 
{
  bool test = true ;
  bool test_horizontals = true ;

  map<unsigned int, RP_struct_type> map_RPs ;
 
  for(int j = 0 ; j < vector_apertures.size() ; ++j)
  {
    double x_pos_meter = ((vector_apertures[j]->Lx * proton.theta_x_star_pert) + (vector_apertures[j]->vx * proton.x_star)) ;
    double y_pos_meter = ((vector_apertures[j]->Ly * proton.theta_y_star_pert) + (vector_apertures[j]->vy * proton.y_star)) ;

    if((fabs(x_pos_meter) < vector_apertures[j]->rectangle_half_width) &&
       (fabs(y_pos_meter) < vector_apertures[j]->rectangle_half_height)  )
    {
      // ((x/A)^2) + ((y/B)^2) == 1

      double x_scaled =  (x_pos_meter / vector_apertures[j]->ellipse_semi_axis_hor) ;

      double y_scaled = sqrt(1 - (x_scaled*x_scaled)) * vector_apertures[j]->ellipse_semi_axis_ver ;

      if(fabs(y_pos_meter) < y_scaled)
      {
        if((vector_apertures[j]->edge == 0) || (fabs(y_pos_meter) > vector_apertures[j]->edge))
        if(!test_horizontals || (test_horizontals && ((vector_apertures[j]->hedge == 0) || (x_pos_meter > vector_apertures[j]->hedge))))
        {
          // if(vector_apertures[j]->hedge != 0) cout << "hedge " << vector_apertures[j]->name << endl ;

          if(save_result)
          {
            vector_apertures[j]->histogram->Fill(proton.theta_x_star_pert, proton.theta_y_star_pert) ;
            vector_apertures[j]->histogram_vertex->Fill(proton.x_star, proton.y_star) ;
            vector_apertures[j]->histogram_t->Fill(minus_t_GeV2) ;
            vector_apertures[j]->histogram_phi->Fill(phi_IP5_rad) ;
            vector_apertures[j]->histogram_x_y_mm->Fill(x_pos_meter, y_pos_meter) ;

            int rpDecId = vector_apertures[j]->getrpDecId(y_pos_meter) ;

            if(rpDecId != 0)
            {
              RP_struct_type my_RP_struct ;
            
              // cout << "weirdtest: " <<  vector_apertures[j]->getrpDecId(y_pos_meter) << endl ;
            
              my_RP_struct.validity = kTRUE ;
              my_RP_struct.rpDecId = rpDecId ;
              my_RP_struct.x = x_pos_meter * 1e3 ;
              my_RP_struct.y = y_pos_meter * 1e3 ;
              my_RP_struct.thx = 0 ;
              my_RP_struct.thy = 0 ;
    
              map_RPs[my_RP_struct.rpDecId] = my_RP_struct ;
            }
            
          }
        }
        else
        {
          test = false ;
        }

      }
      else
      {
        test = false ;
      }
    }
    else
    {
      test = false ;
    }
  }

  bool test_all_apertures_before_overlap_test = false  ;

  if(!test_all_apertures_before_overlap_test || test)
  {  
    for(map<unsigned int, RP_struct_type>::iterator it = map_RPs.begin() ; it != map_RPs.end() ; ++it)
    for(map<unsigned int, RP_struct_type>::iterator it2 = it ; it2 != map_RPs.end() ; ++it2)
    {
      if(it == it2) continue ;
  
      TestDetectorPair(it, it2, 3, 4) ;
      TestDetectorPair(it, it2, 3, 5) ;

      TestDetectorPair(it, it2, 23, 24) ;
      TestDetectorPair(it, it2, 23, 25) ;

      TestDetectorPair(it, it2, 103, 104) ;
      TestDetectorPair(it, it2, 103, 105) ;

      TestDetectorPair(it, it2, 123, 124) ;
      TestDetectorPair(it, it2, 123, 125) ;
    }
  }
  
  
  return test ;
}

void MinimizeHorizontalVerticalPair(string key, vector<THorizontal_and_vertical_xy_pairs_to_match *> &vector_of_corrd)
{
  points = &vector_of_corrd ;

  MinuitFit(key) ;
}

void Minimize()
{
  for(map<string, vector<THorizontal_and_vertical_xy_pairs_to_match *>>:: iterator it = map_of_THorizontal_and_vertical_xy_pairs_to_match.begin() ; it != map_of_THorizontal_and_vertical_xy_pairs_to_match.end() ; ++it)
  {
    cout << endl << endl << "################" << it->first << " " << it->second.size() << endl ;

    actual_detector_combination = it->first ;
    string name_of_hist = "chi2_contribution_" + actual_detector_combination ;
    map_of_hists[actual_detector_combination] = new TH1F(name_of_hist.c_str(), name_of_hist.c_str(), 1000, 0, 100) ;

    MinimizeHorizontalVerticalPair(it->first, it->second) ;
    
    test_hist_ver->SaveAs("test_hist_ver.root") ;
    test_hist_hor->SaveAs("test_hist_hor.root") ;

  }
  
}

void test_aperture(vector<TAperture *> &vector_apertures)
{


  t_GeV2_distribution->SetParameters(Constant, Slope) ;
  t_GeV2_distribution->SetNpx(100000) ;
  
  TH1D *histogram_t = new TH1D("histogram_t", "histogram_t", 1000, 0, 10) ;
  TH1D *histogram_phi = new TH1D("histogram_phi", "histogram_phi", 100, 0, 10) ;

  TH2D *histogram_theta_x_star_rad = new TH2D("histogram_theta_x_star_rad", "histogram_theta_x_star_rad", 100, -1400e-6, 1400e-6,  100, -1400e-6, 1400e-6) ;
  TH2D *histogram_theta_x_y_star_rad = new TH2D("histogram_theta_x_y_star_rad", "histogram_theta_x_y_star_rad", 100, -1400e-6, 1400e-6,  100, -1400e-6, 1400e-6) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_104_103"] = new TH2D("103_104_103", "103_104_103", 100, -20, 20,  100, -20, 20) ;
  map_of_THorizontal_and_vertical_xy_histogram["103_104_104"] = new TH2D("103_104_104", "103_104_104", 100, -20, 20,  100, -20, 20) ;

  map_of_THorizontal_and_vertical_xy_histogram["123_124_123"] = new TH2D("123_124_123", "123_124_123", 100, -20, 20,  100, -20, 20) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_124_124"] = new TH2D("123_124_124", "123_124_124", 100, -20, 20,  100, -20, 20) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_105_103"] = new TH2D("103_105_103", "103_105_103", 100, -20, 20,  100, -20, 20) ;
  map_of_THorizontal_and_vertical_xy_histogram["103_105_105"] = new TH2D("103_105_105", "103_105_105", 100, -20, 20,  100, -20, 20) ;

  map_of_THorizontal_and_vertical_xy_histogram["123_125_123"] = new TH2D("123_125_123", "123_125_123", 100, -20, 20,  100, -20, 20) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_125_125"] = new TH2D("123_125_124", "123_125_124", 100, -20, 20,  100, -20, 20) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_104_corr_x_dx"] = new TH2D("103_104_corr_x_dx", "103_104_corr_x_dx", 100, -20, 20,  100, -0.2, 0.2) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_124_corr_x_dx"] = new TH2D("123_124_corr_x_dx", "123_124_corr_x_dx", 100, -20, 20,  100, -0.2, 0.2) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_105_corr_x_dx"] = new TH2D("103_105_corr_x_dx", "103_105_corr_x_dx", 100, -20, 20,  100, -0.2, 0.2) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_125_corr_x_dx"] = new TH2D("123_125_corr_x_dx", "123_125_corr_x_dx", 100, -20, 20,  100, -0.2, 0.2) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_104_corr_y_dy"] = new TH2D("103_104_corr_y_dy", "103_104_corr_y_dy", 100, -20, 20,  100, -0.2, 0.2) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_124_corr_y_dy"] = new TH2D("123_124_corr_y_dy", "123_124_corr_y_dy", 100, -20, 20,  100, -0.2, 0.2) ;

  map_of_THorizontal_and_vertical_xy_histogram["103_105_corr_y_dy"] = new TH2D("103_105_corr_y_dy", "103_105_corr_y_dy", 100, -20, 20,  100, -0.2, 0.2) ;
  map_of_THorizontal_and_vertical_xy_histogram["123_125_corr_y_dy"] = new TH2D("123_125_corr_y_dy", "123_125_corr_y_dy", 100, -20, 20,  100, -0.2, 0.2) ;

	// for(int i = 0 ; i < 1e6 ; ++i)
	
	int factor = (8600 / 3) ;
	
	for(int i = 0 ; i < (16e4 / factor) ; ++i)
	{
    TProtonPair pp ;
   
    histogram_t->Fill(pp.minus_t_GeV2) ;
    histogram_phi->Fill(pp.phi_IP5_rad) ;
  
    pp.TestApertures(vector_apertures) ;
    
    if(pp.reconstructed)
    {
      histogram_theta_x_star_rad->Fill(pp.p_b1.theta_x_star, pp.theta_x_star_reco) ;
      histogram_theta_x_y_star_rad->Fill(pp.p_b1.theta_x_star_reco, pp.theta_y_star_reco) ;
    }
  }
  
  cout << "######### Minimize start" << endl ;
  Minimize() ;
  cout << "######### Minimize end" << endl ;

	gStyle->SetPadGridX(kTRUE);
	gStyle->SetPadGridY(kTRUE);

  TCanvas c ;
  
  for(int j = 0 ; j < vector_apertures.size() ; ++j)
  {
    vector_apertures[j]->histogram->GetXaxis()->SetTitle("#theta_{x}* (rad)") ;
    vector_apertures[j]->histogram->GetYaxis()->SetTitle("#theta_{y}* (rad)") ;

    vector_apertures[j]->histogram->Draw("colz") ;
    c.SaveAs(("plots/apeture_test/aperture_test_" + vector_apertures[j]->name + ".root").c_str()) ;

    vector_apertures[j]->histogram_vertex->Draw("colz") ;
    c.SaveAs(("plots/apeture_test/aperture_test_" + vector_apertures[j]->name + "_vertex.root").c_str()) ;

    vector_apertures[j]->histogram_t->Draw("") ;
    c.SaveAs(("plots/apeture_test/aperture_test_" + vector_apertures[j]->name + "_t.root").c_str()) ;

    vector_apertures[j]->histogram_phi->Draw("") ;
    c.SaveAs(("plots/apeture_test/aperture_test_" + vector_apertures[j]->name + "_phi.root").c_str()) ;

    vector_apertures[j]->histogram_x_y_mm->Draw("colz") ;
    c.SaveAs(("plots/apeture_test/aperture_test_" + vector_apertures[j]->name + "_x_y.root").c_str()) ;
    
    cout << "statistics " << vector_apertures[j]->name << " " << vector_apertures[j]->histogram->GetEntries() << endl ;
  }
  
  for(map<string, TH2D *>::iterator it = map_of_THorizontal_and_vertical_xy_histogram.begin() ; it != map_of_THorizontal_and_vertical_xy_histogram.end() ; ++it)
  {
    it->second->SaveAs(("plots/apeture_test/xy_overlap_test_" + it->first + ".root").c_str()) ;
  }
  
  histogram_t->SaveAs("histogram_t.root") ;
  histogram_phi->SaveAs("histogram_phi.root") ;
  histogram_theta_x_star_rad->SaveAs("histogram_theta_x_star_rad.root") ;
  histogram_theta_x_y_star_rad->SaveAs("histogram_theta_x_y_star_rad.root") ;
  
  c.cd() ;
  
  histogram_t->Draw("") ;

  for(int j = 0 ; j < vector_apertures.size() ; ++j)
  {
    vector_apertures[j]->histogram_t->Draw("same") ;
  }
  
  c.SaveAs("histogram_t_combined.root") ;
}

void read_apertures()
{

  ifstream apertures(aperture_filename.c_str()) ;
  vector<TAperture *> vector_apertures ;
  
  string name ;
  double rectangle_half_width = 0 ;
  double rectangle_half_height = 0 ;
  double ellipse_semi_axis_hor = 0 ;
  double ellipse_semi_axis_ver = 0 ;
  
  const int RE12_position = 13 ;
  const int RE34_position = 23 ;
  const int BSY_position =  29 ;

  const double LARGE_VALUE = 1e99 ;

  double theta_x_star_limit_due_rectangle_half_width_min = LARGE_VALUE ;
  double theta_x_star_limit_due_ellipse_semi_axis_hor_min = LARGE_VALUE ;

  double theta_y_star_limit_due_rectangle_half_height_min = LARGE_VALUE ;
  double theta_y_star_limit_due_ellipse_semi_axis_ver_min = LARGE_VALUE ;

  string theta_x_star_limit_due_rectangle_half_width_min_name = "" ;
  string theta_x_star_limit_due_ellipse_semi_axis_hor_min_name = "" ;

  string theta_y_star_limit_due_rectangle_half_height_min_name = "" ;
  string theta_y_star_limit_due_ellipse_semi_axis_ver_min_name = "" ;
  
  while(apertures >> name >> rectangle_half_width >> rectangle_half_height >>  ellipse_semi_axis_hor >> ellipse_semi_axis_ver)
  {
    // cout << name << endl ;

    double edge = 0 ;
    double hedge = 0 ;
    
    ifstream optics_file(optics_file_name.c_str()) ;
    
    const string magnet_name = "\"" + name + "\"" ;
    string word ;
    
    bool found = false ;
    
    while(optics_file >> word)
    {
      if(word.compare(magnet_name) == 0)
      {
        found = true ;
        
        double Lx = 0 ;
        double Ly = 0 ;

        double vx = 0 ;
        double vy = 0 ;

        double bsy = 0 ;

        double value = 0 ;
        
        for(int i = 0 ; i < 29 ; ++i)
        {
          optics_file >> value ;
          // cout << value << " " ;

          if(i == (RE12_position - 1)) Lx = value ;
          if(i == (RE34_position - 1)) Ly = value ;

          if(i == (RE12_position - 2)) vx = value ;
          if(i == (RE34_position - 2)) vy = value ;
          
          if(i == (BSY_position - 1))
          {
            bsy = value ;
          
            if((name.compare(RP_near_name) == 0) || (name.compare(RP_far_name) == 0))
            {
              edge = n_sigma * bsy ;
              cout << "megvanedge " << name << " " << edge << " " << endl ;
            }

            if((name.compare(RP_near_hor_name) == 0) || (name.compare(RP_far_hor_name) == 0))
            {
              hedge = n_sigma * bsy ;
              cout << "megvanedge " << name << " " << edge << " " << endl ;
            }
          }
        }
        
        // cout << name << "     " << Lx  << " " << Ly << endl ;
        
        TAperture *aperture = new TAperture(name, Lx, Ly, vx, vy, rectangle_half_width, rectangle_half_height,  ellipse_semi_axis_hor, ellipse_semi_axis_ver, edge, hedge) ;
        vector_apertures.push_back(aperture) ;
        
        double theta_x_star_limit_due_rectangle_half_width = (rectangle_half_width / Lx) ;
        double theta_x_star_limit_due_ellipse_semi_axis_hor = (ellipse_semi_axis_hor / Lx) ;

        double theta_y_star_limit_due_rectangle_half_height = (rectangle_half_height / Ly) ;
        double theta_y_star_limit_due_ellipse_semi_axis_ver = (ellipse_semi_axis_ver / Ly) ;

        cout << name << " theta_x_star_limit_due_rectangle_half_width_min_name:  " << theta_x_star_limit_due_rectangle_half_width_min << endl ;
        cout << name << " theta_x_star_limit_due_ellipse_semi_axis_hor_min_name: " << theta_x_star_limit_due_ellipse_semi_axis_hor_min << endl ;

        cout << name << " theta_y_star_limit_due_rectangle_half_height_min_name: " << theta_y_star_limit_due_rectangle_half_height_min << endl ;
        cout << name << " theta_y_star_limit_due_ellipse_semi_axis_ver_min_name: " << theta_y_star_limit_due_ellipse_semi_axis_ver_min << endl << endl ;
        
        if(theta_x_star_limit_due_rectangle_half_width_min > theta_x_star_limit_due_rectangle_half_width)   { theta_x_star_limit_due_rectangle_half_width_min_name = name ; theta_x_star_limit_due_rectangle_half_width_min = theta_x_star_limit_due_rectangle_half_width ; } ;
        if(theta_x_star_limit_due_ellipse_semi_axis_hor_min > theta_x_star_limit_due_ellipse_semi_axis_hor) { theta_x_star_limit_due_ellipse_semi_axis_hor_min_name = name ;theta_x_star_limit_due_ellipse_semi_axis_hor_min = theta_x_star_limit_due_ellipse_semi_axis_hor ; } ;

        if(theta_y_star_limit_due_rectangle_half_height_min > theta_y_star_limit_due_rectangle_half_height) { theta_y_star_limit_due_rectangle_half_height_min_name = name ; theta_y_star_limit_due_rectangle_half_height_min = theta_y_star_limit_due_rectangle_half_height ; } ;
        if(theta_y_star_limit_due_ellipse_semi_axis_ver_min > theta_y_star_limit_due_ellipse_semi_axis_ver) { theta_y_star_limit_due_ellipse_semi_axis_ver_min_name = name ; theta_y_star_limit_due_ellipse_semi_axis_ver_min = theta_y_star_limit_due_ellipse_semi_axis_ver ; } ;
        
        // test_aperture(name, Lx, Ly, rectangle_half_width, rectangle_half_height,  ellipse_semi_axis_hor, ellipse_semi_axis_ver) ;
      }
    }
    
    optics_file.close() ;
    
    if(!found)
    {
      cout << name << " not found" << endl ;
    }
  }
  
  cout << endl ;
  cout << "#################" << endl ;  
  cout << endl ;

  cout << "theta_x_star_limit_due_rectangle_half_width_min_name:  " <<  theta_x_star_limit_due_rectangle_half_width_min_name  << theta_x_star_limit_due_rectangle_half_width_min << endl ;
  cout << "theta_x_star_limit_due_ellipse_semi_axis_hor_min_name: " <<  theta_x_star_limit_due_ellipse_semi_axis_hor_min_name << theta_x_star_limit_due_ellipse_semi_axis_hor_min << endl ;
  
  cout << "theta_y_star_limit_due_rectangle_half_height_min_name: " <<  theta_y_star_limit_due_rectangle_half_height_min_name << theta_y_star_limit_due_rectangle_half_height_min << endl ;
  cout << "theta_y_star_limit_due_ellipse_semi_axis_ver_min_name: " <<  theta_y_star_limit_due_ellipse_semi_axis_ver_min_name << theta_y_star_limit_due_ellipse_semi_axis_ver_min << endl ;

  test_aperture(vector_apertures) ;
}

void plot_contour()
{
  string contour_45b_56t_L = "[[-200e-6, 550e-6], [-500e-6, 500e-6], [-560e-6, 220e-6], [-400e-6, 170e-6], [+200e-6, 170e-6], [+440e-6, 260e-6], [+370e-6, 500e-6], [+340e-6, 550e-6]]" ;
  string contour_45t_56b_L = "[[-200e-6, 170e-6], [+400e-6, 170e-6], [+520e-6, 200e-6], [+400e-6, 480e-6], [-430e-6, 480e-6], [-550e-6, 400e-6], [-600e-6, 200e-6]]" ;
  
  ifstream mycontours("data/contours_of_Jan/contours.txt") ;
  string word ;
  string contour_name, an_equal_sign ;
  double x_coord, y_coord ;
  string x_coord_str, y_coord_str ;
  
  mycontours >> contour_name ;
  mycontours >> an_equal_sign ;

  TGraph *graph = new TGraph() ;
  
  string output_string = "" ;
  
  while(mycontours >> x_coord >> y_coord)
  {
    graph->AddPoint(x_coord, y_coord) ;
  }

  mycontours.close() ;
  mycontours.open("data/contours_of_Jan/contours.txt") ;

  mycontours >> contour_name ;
  mycontours >> an_equal_sign ;

  while(mycontours >> x_coord_str >> y_coord_str)
  {
    output_string += "[" + x_coord_str + ", " + y_coord_str + "], " ;
  }
  
  mycontours.close() ;
  
  cout << "[" << output_string << "]" << endl ;
  
  TCanvas c ;
  
  graph->Draw("alp") ;
  graph->SetMarkerStyle(20) ;
  graph->SaveAs("plots/contour.root") ;
  c.SaveAs("plots/contour.pdf") ;

}

int main()
{
  // gErrorIgnoreLevel = 6001 ;

  const int main_scenario_plot_fiducial_cuts = 1 ;
  const int main_scenario_plot_apertures = 2 ;
  const int main_scenario_optics_test = 3 ;

  int main_scenario = main_scenario_plot_apertures ;

  t_GeV2_distribution = new TF1("t_GeV2_distribution", my_exponential_distribution, 0.0, 7.0, 2) ;

  // gStyle->SetOptStat(0) ;
  // gStyle->SetOptFit(1111) ;

  TF1 *t_GeV2_distribution = new TF1("t_GeV2_distribution", my_exponential_distribution, 0.0, 7.0, 2) ;
  t_GeV2_distribution->SetParameters(Constant, Slope) ;
  t_GeV2_distribution->SetNpx(100000) ;

  //   TH1D *hist_minus_t_GeV2 = new TH1D("hist_minus_t_GeV2", "hist_minus_t_GeV2", 100, 0, 2) ;
  //   TH1D *hist_minus_t_GeV2_reco = new TH1D("hist_minus_t_GeV2_reco", "hist_minus_t_GeV2_reco", 100, 0, 2) ;
  //   TH1D *hist_minus_t_GeV2_reco_diff_rel_percent = new TH1D("hist_minus_t_GeV2_reco_diff_rel_percent", "hist_minus_t_GeV2_reco_diff_rel_percent", 100, 1e-2, -1e-2) ;

  TH1D *hist_minus_t_GeV2 = new TH1D("hist_minus_t_GeV2", "hist_minus_t_GeV2", 100, 0, 2) ;
  TH1D *hist_minus_t_GeV2_reco = new TH1D("hist_minus_t_GeV2_reco", "hist_minus_t_GeV2_reco", 100, 0, 2) ;

  TH1D *hist_minus_t_GeV2_reco_diff_rel_percent = new TH1D("hist_minus_t_GeV2_reco_diff_rel_percent", "hist_minus_t_GeV2_reco_diff_rel_percent", 100, 1e-2, -1e-2) ;

  TH2D *hist = new TH2D("hist", "hist", 100, -5e+2, 5e+2, 100, -1000, 1000) ;

  TH2D *hist_compare_theta_y_star = new TH2D("hist_compare_theta_y_star", "hist_compare_theta_y_star", 100, -16.0e-4, 16.0e-4,  100, -16.0e-4, 16.0e-4) ;
  TH2D *hist_compare_y_star = new TH2D("hist_compare_y_star", "hist_compare_y_star", 100, -16.0e-4, 16.0e-4,  100, -16.0e-4, 16.0e-4) ;

  TH2D *hist_compare_theta_x_star = new TH2D("hist_compare_theta_x_star", "hist_compare_theta_x_star", 100, -16.0e-4, 16.0e-4,  100, -16.0e-4, 16.0e-4) ;
  TH2D *hist_compare_x_star = new TH2D("hist_compare_x_star", "hist_compare_x_star", 100, -16.0e-4, 16.0e-4,  100, -16.0e-4, 16.0e-4) ;

  TH2D *hist_theta_x_y_star_reco = new TH2D("hist_theta_x_y_star_reco", "hist_theta_x_y_star_reco", 100, -1000e-6, 1000e-6, 100, -1000e-6, 1000e-6) ;

  // Best effect
  const double delta_n = -2.0 * 0.0 ;
  const double delta_f = -3.0 * 0.0 ;

  // 

  const double Lx_near_reco = Lx_near + delta_n ;
  const double Lx_far_reco = Lx_far   + delta_f ;

  const double Ly_near_reco = Ly_near + delta_n ;
  const double Ly_far_reco = Ly_far   + delta_f ;
  

  if(main_scenario == main_scenario_plot_fiducial_cuts)
  {
    plot_contour() ;
  }
  else if(main_scenario == main_scenario_plot_apertures)
  {

    read_apertures() ;

  }
  else if(main_scenario == main_scenario_optics_test)
  {

    for(int i = 0 ; i < 1e6 ; ++i)
    {

      const double minus_t_GeV2 = t_GeV2_distribution->GetRandom() ;
      const Double_t phi_IP5_rad = gRandom->Uniform(0, PI_TIMES_2) ;    

      double theta_star_rad = theta_star_rad_from_t_GeV2(-minus_t_GeV2, beam_momentum_GeV) ;

      double theta_x_star = (theta_star_rad * cos(phi_IP5_rad)) ;
      double theta_y_star = (theta_star_rad * sin(phi_IP5_rad)) ;

	    double y_star = 1.0 * 200e-6 * myrand.Gaus() ;
	    double x_star = 1.0 * 200e-6 * myrand.Gaus() ;

	    double beam_divergence_x = 0.0 * 20e-6 * myrand.Gaus() ;
	    double beam_divergence_y = 0.0 * 20e-6 * myrand.Gaus() ;

	    double theta_x_star_pert = theta_x_star + beam_divergence_x ;
	    double theta_y_star_pert = theta_y_star + beam_divergence_y ;

	    double x_near = (Lx_near * theta_x_star_pert) + (vx_near * x_star) ;
	    double x_far =  (Lx_far *  theta_x_star_pert) + (vx_far *  x_star) ;
	    double y_near = (Ly_near * theta_y_star_pert) + (vy_near * y_star) ;
	    double y_far =  (Ly_far *  theta_y_star_pert) + (vy_far *  y_star) ;

	    // cout << y_far << endl ;
    //		if(fabs(x_far)< 3e-3)
    //		if(fabs(x_near)< 3e-3)

	    // if((y_far < -1.e-3) && (y_far > -10.0e-3) )
      //	if(fabs(x_near) < 3e-3)
      //	if(fabs(x_far) < 3e-3)
      //	if((y_far < -4.0e-3) && (y_far > -10.0e-3) && (y_near < (-4.0e-3*1.4*(Ly_far/Ly_near))) && (y_near > (-10.0e-3*(Ly_far/Ly_near))))
	    {

        double determinant_x = ((Lx_near_reco * vx_far) - (Lx_far_reco * vx_near)) ;
        double determinant_y = ((Ly_near_reco * vy_far) - (Ly_far_reco * vy_near)) ;

		    double theta_y_star_reco = ((y_near * vy_far) - (y_far * vy_near)) / determinant_y ;
		    // double theta_y_star_reco = ((y_near / Ly_near_reco) + (y_far / Ly_far_reco)) / 2.0 ;
		    double y_star_reco =      -((y_near * Ly_far_reco) - (y_far * Ly_near_reco)) / determinant_y ;

		    double theta_x_star_reco = ((x_near * vx_far) - (x_far * vx_near)) / determinant_x ;
		    double x_star_reco =      -((x_near * Lx_far_reco) - (x_far * Lx_near_reco)) / determinant_x ;

        double theta_star_reco = sqrt((theta_x_star_reco*theta_x_star_reco) + (theta_y_star_reco*theta_y_star_reco)) ;
        double minus_t_GeV2_reco = - beam_momentum_GeV * beam_momentum_GeV * theta_star_reco * theta_star_reco ;

        hist_minus_t_GeV2->Fill(minus_t_GeV2) ;
        hist_minus_t_GeV2_reco->Fill(-minus_t_GeV2_reco) ;

        hist_minus_t_GeV2_reco_diff_rel_percent->Fill(((minus_t_GeV2 + minus_t_GeV2_reco) / minus_t_GeV2) * 100.0) ;

		    hist_compare_theta_x_star->Fill(theta_x_star_pert, theta_x_star_reco) ;
		    hist_compare_x_star->Fill(x_star, x_star_reco) ;

		    hist_compare_theta_y_star->Fill(theta_y_star_pert, theta_y_star_reco) ;
		    hist_compare_y_star->Fill(y_star, y_star_reco) ;

        hist_theta_x_y_star_reco->Fill(theta_x_star_reco, theta_y_star_reco) ;

		    hist->Fill(theta_y_star_reco * 1e6, y_star_reco * 1e6) ;
	    }

    }

    gStyle->SetPalette(1) ;

    TCanvas c ;

    gStyle->SetTitleFillColor(kWhite) ;
    gStyle->SetStatColor(kWhite) ;
    gStyle->SetFillColor(kWhite) ;

    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);

    hist->SetFillColor(kWhite) ;
    hist->SetTitle("") ;

    hist->Draw("colz") ;

    hist->GetXaxis()->SetTitle("#theta_{x,reco}* [#murad]") ;
    hist->GetYaxis()->SetTitle("x_{reco}* [#mum]") ;
    hist->GetYaxis()->SetTitleOffset(1.2) ;

    TLatex *latex = new TLatex() ;
    latex->SetTextFont(132) ;
    latex->SetNDC();
    latex->DrawText(0.24, 0.92, "simulation preliminary") ;

    latex->SetTextFont(22) ;
    latex->DrawText(0.12, 0.92, "TOTEM") ;
    latex->Draw("same") ;


    c.SetFillColor(kWhite) ;

    //	gStyle->SetOptStat(0);

    c.SetGridx() ;
    c.SetGridy() ;	

    hist->GetXaxis()->SetAxisColor(17);
    hist->GetYaxis()->SetAxisColor(17);	

    c.SetLogz() ;
    c.SaveAs("plots/hist.pdf") ;

    hist_compare_theta_x_star->Draw("colz") ;
    c.SaveAs("plots/hist_compare_theta_x_star.pdf") ;

    hist_compare_x_star->Draw("colz") ;
    c.SaveAs("plots/hist_compare_x_star.pdf") ;

    hist_compare_theta_y_star->Draw("colz") ;
    c.SaveAs("plots/hist_compare_theta_y_star.pdf") ;

    hist_compare_y_star->Draw("colz") ;
    c.SaveAs("plots/hist_compare_y_star.pdf") ;

    hist_theta_x_y_star_reco->Draw("colz") ;
    c.SaveAs("plots/hist_theta_x_y_star_reco.pdf") ;
    c.SaveAs("plots/hist_theta_x_y_star_reco.root") ;

    // cout >> hist_minus_t_GeV2->Integral() ;
    hist_minus_t_GeV2->Draw("colz") ;
    hist_minus_t_GeV2->GetYaxis()->SetRangeUser(0.1, 1e7) ;
    c.SetLogy() ;
    c.SaveAs("plots/hist_minus_t_GeV2.pdf") ;

    hist_minus_t_GeV2_reco->Draw("colz") ;
    hist_minus_t_GeV2_reco->GetYaxis()->SetRangeUser(0.1, 1e7) ;
    c.SetLogy() ;
    c.SaveAs("plots/hist_minus_t_GeV2_reco.pdf") ;

    hist_minus_t_GeV2_reco_diff_rel_percent->Draw("colz") ;
    c.SetLogy(0) ;
    c.SaveAs("plots/hist_minus_t_GeV2_reco_diff_rel_percent.pdf") ;
  
  }
  else
  {
    cout << "Unknown scenario: " << main_scenario << endl ;
  }
}	

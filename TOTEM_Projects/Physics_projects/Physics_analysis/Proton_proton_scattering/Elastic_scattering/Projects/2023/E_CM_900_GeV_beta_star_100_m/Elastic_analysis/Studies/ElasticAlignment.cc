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

TH2D *flipHistogram(TH2D *histogram)
{
	int xbins = histogram->GetXaxis()->GetNbins() ;
	int ybins = histogram->GetYaxis()->GetNbins() ;

	int x_axis_min = histogram->GetXaxis()->GetXmin() ;
	int x_axis_max = histogram->GetXaxis()->GetXmax() ;

	int y_axis_min = histogram->GetYaxis()->GetXmin() ;
	int y_axis_max = histogram->GetYaxis()->GetXmax() ;

	TH2D *histogram_rotated = new TH2D("test", "test", ybins, y_axis_min, y_axis_max, xbins, x_axis_min, x_axis_max) ;
	
	for(int i = 0 ; i < xbins ; ++i)
	for(int j = 0 ; j < ybins ; ++j)
	{
		double bincontent = histogram->GetBinContent(i, j) ;
		double bin_unc = histogram->GetBinError(i, j) ;
	
		histogram_rotated->SetBinContent(j, i, bincontent) ;
		histogram_rotated->SetBinError(j, i, bin_unc) ;
	}
	
	// histogram_rotated->SaveAs("histogram_rotated.root") ;
	
	return histogram_rotated ;
}


TFile *file_LBRT = NULL ;
TFile *file_LTRB = NULL ;

void combine_and_fit(string name, double low, double high, bool swap = false)
{
	TCanvas *c = new TCanvas() ;

	TH2D *LBRT = ((TH2D *)(file_LBRT->Get(name.c_str()))) ;
	TH2D *LTRB = ((TH2D *)(file_LTRB->Get(name.c_str()))) ;
	
	c->SetTitle(LBRT->GetName()) ;
	cout << "myname " << LBRT->GetName() << endl ;
	c->SetTitle("") ;
	
	// LBRT->SaveAs("delete_it.root") ;

	TH1D *px1 = NULL ;
	TH1D *px2 = NULL ;
	
	if(swap)
	{
		px1 = LBRT->ProjectionX("px1", 600, 660) ;
		px2 = LTRB->ProjectionX("px2", 366, 425) ;	
	}
	else
	{
		px1 = LTRB->ProjectionX("px1", 600, 660) ;
		px2 = LBRT->ProjectionX("px2", 366, 425) ;	
	}
	
	double px1_mean = px1->GetMean() ;
	double px2_mean = px2->GetMean() ;

	TFitResultPtr ptr1 = px1->Fit("gaus","S", "", px1_mean + low, px1_mean + high) ;
	TFitResultPtr ptr2 = px2->Fit("gaus","S", "", px2_mean + low, px2_mean + high) ;
	
	double mean_1 = ptr1->Parameter(1) ;
	double mean_2 = ptr2->Parameter(1) ;

	double avg = (mean_1 + mean_2) / 2.0 ;
	double var = fabs(avg - mean_1) ;
	
	cout << "fit_resultx mean: " << mean_1 << endl ;
	cout << "fit_resultx mean: " << mean_2 << endl ;
	cout << "fit_resultx avg: " << avg << endl ;
	cout << "fit_resultx var: " << var << endl ;
	cout << "fit_resultx" << endl ;

	px1->Draw() ;
	px1->GetXaxis()->SetRangeUser(-10.0, 10.0) ;
	px2->Draw("same") ;

	c->SaveAs(("plots/ElasticAlignment/proj_" + name + ".root").c_str()) ;
	c->SaveAs(("plots/ElasticAlignment/proj_" + name + ".pdf").c_str()) ;
	
	delete px1 ;
	delete px2 ;
	
	delete c ;

}

void combine_and_fit_y(string name, double low, double high, int bins, int fill)
{
	TCanvas *c = new TCanvas() ;
	
	string filename = "data/acceptance_cuts_for_vertical_alignment/acceptance_cuts_for_vertical_alignment.txt" ;
	
	ifstream acceptance_cuts(filename.c_str()) ;
	
	double LBRT_cut_y_mm_low = 0 ;
	double LBRT_cut_y_mm_high = 0 ;

	double LTRB_cut_y_mm_low = 0 ;
	double LTRB_cut_y_mm_high = 0 ;
	
	// LBRT
	
	string arm = "" ;
	string position = "" ;
	double cut_position_rad = 0 ;
	string diagonal = "" ;
	double Ly_near = 0.0 ;
	double Ly_far = 0.0 ;
	double far_cut_mm = 0.0 ;
	double near_cut_mm = 0.0 ;
	
	string selected_arm = "" ;
	string selected_rp = "" ;
	
	if(name.compare("P0079_PlotsCollection_x_mm_y_mm_near_left_aligned") == 0)
	{
		selected_arm = "left" ;
		selected_rp =  "near" ;
	} 
	else if(name.compare("P0080_PlotsCollection_x_mm_y_mm_far_left_aligned") == 0)
	{
		selected_arm = "left" ;
		selected_rp =  "far" ;
	}
	else if(name.compare("P0081_PlotsCollection_x_mm_y_mm_near_right_aligned") == 0)
	{
		selected_arm = "right" ;
		selected_rp =  "near" ;
	}
	else if(name.compare("P0082_PlotsCollection_x_mm_y_mm_far_right_aligned") == 0)
	{
		selected_arm = "right" ;
		selected_rp =  "far" ;
	}
	else
	{
		exit(1) ;
	}

	
	while(acceptance_cuts >> arm >> position >> cut_position_rad >> diagonal >> Ly_near >> Ly_far >> far_cut_mm >> near_cut_mm)
	{
		// cout << " " <<   arm << " " <<  position << " " <<  cut_position_rad << " " <<  diagonal << " " <<  Ly_near << " " <<  Ly_far << " " <<  far_cut_mm << " " <<  near_cut_mm << " " <<  endl ;

		if(diagonal.compare("LBRT") == 0)
		{
			double sign_and_unit = 1e3 ;

			if(selected_arm.compare("left")) sign_and_unit *= -1.0 ;

			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("near") == 0) && (position.compare("low") == 0)) LBRT_cut_y_mm_low = sign_and_unit * near_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("near") == 0) && (position.compare("high") == 0)) LBRT_cut_y_mm_high = sign_and_unit * near_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("far") == 0) && (position.compare("low") == 0)) LBRT_cut_y_mm_low = sign_and_unit * far_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("far") == 0) && (position.compare("high") == 0)) LBRT_cut_y_mm_high = sign_and_unit * far_cut_mm ;
		}

		if(diagonal.compare("LTRB") == 0)
		{
			double sign_and_unit = 1e3 ;

			if(selected_arm.compare("right")) sign_and_unit *= -1.0 ;

			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("near") == 0) && (position.compare("low") == 0)) LTRB_cut_y_mm_low = sign_and_unit * near_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("near") == 0) && (position.compare("high") == 0)) LTRB_cut_y_mm_high = sign_and_unit * near_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("far") == 0) && (position.compare("low") == 0)) LTRB_cut_y_mm_low = sign_and_unit * far_cut_mm ;
			if((arm.compare(selected_arm) == 0) && (selected_rp.compare("far") == 0) && (position.compare("high") == 0)) LTRB_cut_y_mm_high = sign_and_unit * far_cut_mm ;
		}
	}

	cout <<  "LBRT_cut_y_mm_low " << LBRT_cut_y_mm_low << endl ;
	cout <<  "LBRT_cut_y_mm_high " << LBRT_cut_y_mm_high << endl ;

	cout << "LTRB_cut_y_mm_low " << LTRB_cut_y_mm_low << endl ;
	cout << "LTRB_cut_y_mm_high " << LTRB_cut_y_mm_high << endl ;
	cout << endl ;
	

	// return ;

	TH2D *LBRT = ((TH2D *)(file_LBRT->Get(name.c_str()))) ;
	TH2D *LTRB = ((TH2D *)(file_LTRB->Get(name.c_str()))) ;
	
	TH2D *combined = (TH2D *)LBRT->Clone(name.c_str()) ;
	combined->Add(LTRB) ;
	
	combined->SaveAs(("plots/ElasticAlignment/combined_" + name + ".root").c_str()) ;

	cout << "vertmyname " << LBRT->GetName() << endl ;

	TH1D *py1 = NULL ;
	TH1D *py2 = NULL ;
	TH1D *py3 = NULL ;
	
	int bin1 = combined->GetXaxis()->FindBin(-2.0) ;
	int bin2 = combined->GetXaxis()->FindBin( 2.0) ;
	cout << "bin1 " << bin1 << endl ;
	cout << "bin2 " << bin2 << endl ;
	
	// py1 = combined->ProjectionY("py1", bin1, bin2) ; 
	// py1 = combined->ProjectionY("py1", bin1, bin2) ; 
	py1 = combined->ProjectionY("py1") ; 
	py2 = LBRT->ProjectionY("py2") ; 
	py3 = LTRB->ProjectionY("py3") ; 

	py1->SetLineColor(kRed) ;

	//py2->Draw() ;
	//py3->Draw("same") ;
	py1->Draw("") ;
	
	py1->SaveAs(("plots/ElasticAlignment/projy_" + name + "_hist_before.root").c_str()) ;
	
	// py1->GetXaxis()->SetRangeUser(-30.0, 30.0) ;
	
	int nbins = py1->GetXaxis()->GetNbins() ;
	
	for(int i = 0 ; i < nbins ; ++i)
	{
		double biny = py1->GetXaxis()->GetBinCenter(i) ;
		// cout << biny << endl ;

		if((LBRT_cut_y_mm_low < 0) && ((LBRT_cut_y_mm_high > biny) || ((biny > LBRT_cut_y_mm_low) && (biny < LTRB_cut_y_mm_low)) || (LTRB_cut_y_mm_high < biny)))
		{
			// cout << "here I am2" << endl ;
			// if(fill == 7301)
			{
				py1->SetBinContent(i, 0) ;
				py1->SetBinError(i, 0) ;
			}
		}
		
		if((LBRT_cut_y_mm_low > 0) && ((LTRB_cut_y_mm_high > biny) || ((biny > LTRB_cut_y_mm_low) && (biny < LBRT_cut_y_mm_low)) || (LBRT_cut_y_mm_high < biny)))
		{
			cout << "here I am" << endl ;

			// if(fill == 7301)
			{
				py1->SetBinContent(i, 0) ;
				py1->SetBinError(i, 0) ;
			}
		}




	}
	
	TFitResultPtr ptr = py1->Fit("gaus", "S", "", low, high) ;
	py1->GetFunction("gaus")->SetLineColor(kBlue) ;
	double mean = ptr->Parameter(1) ;
	cout << "fit_resulty mean: " << name << " " << mean << endl ;
	
	py1->SetLineColor(kBlue) ;

	py1->SaveAs(("plots/ElasticAlignment/projy_" + name + "_hist.root").c_str()) ;
	c->SaveAs(("plots/ElasticAlignment/projy_" + name + ".root").c_str()) ;
	c->SaveAs(("plots/ElasticAlignment/projy_" + name + ".pdf").c_str()) ;

	delete py1 ;
	delete py2 ;
	delete py3 ;

}

main()
{
	gStyle->SetLineScalePS(.3) ;
	gStyle->SetOptFit(1111);

	// file_LBRT = TFile::Open("/afs/cern.ch/work/f/fnemes/tmp/pp/E_CM_900_GeV_beta_star_11_m/Analysis_output_files/7301/Diagonals/DIAGONAL_LEFT_BOTTOM_RIGHT_TOP/All_root_files_to_define_cuts/Generic.root") ;
	// file_LTRB = TFile::Open("/afs/cern.ch/work/f/fnemes/tmp/pp/E_CM_900_GeV_beta_star_11_m/Analysis_output_files/7301/Diagonals/DIAGONAL_LEFT_TOP_RIGHT_BOTTOM/All_root_files_to_define_cuts/Generic.root") ;
	
	int fill = 7302 ;

	if(fill == 7301)
	{
		file_LBRT = TFile::Open("/eos/cms/store/group/phys_diffraction/fnemes/E_CM_900_GeV/Beta_star_11_m/Analysis_output_files/7301/Diagonals/DIAGONAL_LEFT_BOTTOM_RIGHT_TOP/All_root_files_to_define_cuts_aligned_first_optics_test/Generic.root") ;
		file_LTRB = TFile::Open("/eos/cms/store/group/phys_diffraction/fnemes/E_CM_900_GeV/Beta_star_11_m/Analysis_output_files/7301/Diagonals/DIAGONAL_LEFT_TOP_RIGHT_BOTTOM/All_root_files_to_define_cuts_aligned_first_optics_test/Generic.root") ;

		combine_and_fit_y("P0079_PlotsCollection_x_mm_y_mm_near_left_aligned", -10.0, 10.0, 90, fill) ;
		combine_and_fit_y("P0080_PlotsCollection_x_mm_y_mm_far_left_aligned", -10.0, 10.0,  90, fill) ;
		combine_and_fit_y("P0081_PlotsCollection_x_mm_y_mm_near_right_aligned", -9.0, 9.0,  70, fill) ;
		combine_and_fit_y("P0082_PlotsCollection_x_mm_y_mm_far_right_aligned",  -9.0, 9.0,  70, fill) ;
	}
	else if(fill == 7302)
	{
		file_LBRT = TFile::Open("/eos/cms/store/group/phys_diffraction/fnemes/E_CM_900_GeV/Beta_star_11_m/Analysis_output_files/7302/Diagonals/DIAGONAL_LEFT_BOTTOM_RIGHT_TOP/All_root_files_to_define_cuts_aligned_first_optics_test/Generic.root") ;
		file_LTRB = TFile::Open("/eos/cms/store/group/phys_diffraction/fnemes/E_CM_900_GeV/Beta_star_11_m/Analysis_output_files/7302/Diagonals/DIAGONAL_LEFT_TOP_RIGHT_BOTTOM/All_root_files_to_define_cuts_aligned_first_optics_test/Generic.root") ;

		combine_and_fit_y("P0079_PlotsCollection_x_mm_y_mm_near_left_aligned", -10.0, 10.0, 90, fill) ;
		combine_and_fit_y("P0080_PlotsCollection_x_mm_y_mm_far_left_aligned", -10.0, 10.0,  90, fill) ;
		combine_and_fit_y("P0081_PlotsCollection_x_mm_y_mm_near_right_aligned", -9.0, 9.0,  70, fill) ;
		combine_and_fit_y("P0082_PlotsCollection_x_mm_y_mm_far_right_aligned",  -9.0, 9.0,  70, fill) ;
	}


	// combine_and_fit("P0079_PlotsCollection_x_mm_y_mm_near_left_aligned", -2.0, 2.0) ;
	// combine_and_fit("P0080_PlotsCollection_x_mm_y_mm_far_left_aligned", -1.5, 1.5) ;
	// combine_and_fit("P0081_PlotsCollection_x_mm_y_mm_near_right_aligned", -2.0, 2.0, true) ;
	// combine_and_fit("P0082_PlotsCollection_x_mm_y_mm_far_right_aligned", -1.5, 1.5, true) ;
}

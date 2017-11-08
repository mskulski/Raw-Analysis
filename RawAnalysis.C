#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TCut.h"
#include <TRandom3.h>
#include <TF1.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include <sstream> 
#include <fstream>
#include <stdlib.h>
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
using namespace std;

void RawAnalysis(Int_t runno,Int_t stable_reg)
{
	vector<Double_t> event;
	vector<vector<Double_t>> scalers(32);
	vector<vector<Double_t>> scalers_add(32);
	vector<Double_t> concentration;
	vector<Double_t> timeplot;
	
	vector<Double_t> currents0;
	vector<Double_t> currents1;
	stringstream runs;
	Int_t counts = 0;
	
	runs << "run" << runno << ".root";
	TString file = runs.str();
	TFile *file1 = TFile::Open(file);
	
	TCanvas* can = new TCanvas("Canvas","",900,900);
	evtTree->Draw(">>evtlist","95 > sqrt((Pos_X-2458)*(Pos_X-2458))","goff");
	if(evtlist)
	{
		Int_t ntot = evtlist->GetN();
	}
	Int_t time = (scalerTree->GetEntries()-1)*2;
	
	for(int m=0;m<(time/2)+1;m++)
	{
		scalerTree->GetEvent(m);
		timeplot.push_back(m*2);
		for(int k=0;k<32;k++)
		{
			Double_t r = scalerTree->GetLeaf("scalers")->GetValue(k);			
			scalers[k].push_back(r);
		}	
	}
	
	Double_t scalersum;
	Double_t scalerstotal;
	
	for(int k=0;k<32;k++)
		{
			scalersum=0;
			scalerstotal=0;
			for(int n=0;n<scalers[k].size();n++)
			{
				scalerstotal += scalers[k][n];
			}
			for(int m=0;m<scalers[k].size();m++)
			{
				scalersum += scalers[k][m];
				scalers_add[k].push_back(scalersum/scalerstotal);
			}
		}
	
	Int_t totals =0;
	for(int l=0;l<(time/2)+1;l++)
	{
		counts += scalers[8][l];
		Double_t events = 0;
		Double_t events_corrected = 0;
		if(evtlist)
		{
			for(int s=0;s<ntot;s++)
			{
				Int_t eno = evtlist->GetEntry(s);
				if(eno > counts && eno <= counts + scalers[8][l+1])
				{
					events++;
					totals++;
				}
			}
			if(scalers[8][l] != 0)
			{
				events_corrected = (scalers[9][l]/scalers[8][l])*(events);
			}
			else events_corrected = events;
			
			event.push_back(events_corrected);
		}
	}
	
	stringstream reg0;
	stringstream reg1;
	reg0 << "reg0_current_run" << runno << ".txt";
	reg1 << "reg1_current_run" << runno << ".txt";
	TString reg0file = reg0.str();
	TString reg1file = reg1.str();
	ifstream reg0currentfile;
	ifstream reg1currentfile;
	
	Int_t dT=2;
	Int_t dt=0;
	Int_t Nt=0;
	Double_t current=0;
	Double_t currentsum=0;
	TString null;
	
	reg0currentfile.open(reg0file);
	reg1currentfile.open(reg1file);
	while(1)
	{
		reg0currentfile >> dt >> null >> current >> endl;
		if(dt<=dT)
		{
			Nt++;
			currentsum+=current;
		}
		else if(dt>dT)
		{
			currents0.pushback(currentsum/Nt);
			Nt=1;
			currentsum=current;
			dT+=2;
		}
		if(!reg0currentfile.good())
		{break;}
	}
	
	dT=2;
	dt=0;
	Nt=0;
	current=0;
	currentsum=0;
	
	while(1)
	{
		reg1currentfile >> dt >> null >> current >> endl;
		if(dt<=dT)
		{
			Nt++;
			currentsum+=current;
		}
		else if(dt>dT)
		{
			currents1.pushback(currentsum/Nt);
			Nt=1;
			currentsum=current;
			dT+=2;
		}
		if(!reg1currentfile.good())
		{break;}
	}
	reg0currentfile.close();
	reg1currentfile.close();
	
	Double_t concentration_sum=0;
	for(int m=0;m<(time/2)+1;m++)
	{
		timeplot.push_back(m*2);
		if(stable_reg==0)
		{
		concentration_sum+=(1.602e-4*event[m]/currents0[m]);
		concentration.push_back(concentration_sum/(m+1));
		}
		else if(stable_reg==1)
		{
		concentration_sum+=(1.602e-4*event[m]/currents1[m]);
		concentration.push_back(concentration_sum/(m+1));
		}
	}
	Double_t concentration_final = concentration[time/2];
	
	std::ofstream myfile;
	myfile.open("RawConcentrations.txt", std::ofstream::app);
	myfile << runno << " " << concentration_final << " " << totals << endl;
	myfile.close();
	
	TCanvas* concentration_canvas = new TCanvas("Concentration","Concentration",3600,2700);
	TGraph* concentration_graph = new TGraph((time/2)+1,&timeplot[0],&concentration[0]);
	concentration_graph->SetTitle("Concentration vs. Time, run "+runno+", 14C/12C = "+concentration_final);
	concentration_graph->GetXaxis()->SetTitle("Time (seconds)");
	concentration_graph->GetXaxis()->CenterTitle();
	concentration_graph->GetXaxis()->SetTitleOffset(1.4);
	concentration_graph->GetXaxis()->SetLabelSize(0.025);
	concentration_graph->GetYaxis()->SetTitle("Concentration (10^-15)");
	concentration_graph->GetYaxis()->CenterTitle();
	concentration_graph->GetYaxis()->SetTitleOffset(1.4);
	concentration_graph->GetYaxis()->SetLabelSize(0.025);
	concentration_graph->Draw("AL");
	
	TImage *img = TImage::Create();
	img->FromPad(concentration_canvas);
	stringstream imgname;
	imgname << "run"<< runno << ".png";
	TString imgf = imgname.str();
	img->WriteImage(imgf);
}

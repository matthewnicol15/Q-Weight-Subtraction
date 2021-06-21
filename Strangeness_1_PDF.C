{
TFile f1("/mnt/f/PhD/Analysis_Output/RGA/Skim4/Inbending/Strangeness_1/PID/Strangeness_1_RGA_Skim4_e_Kp_FD_ppi_no_FD_Inbending_190521_01.root");

TH1F *hmass_kp = (TH1F*)f1.Get("hmass_kp");
TH1F *hmiss = (TH1F*)f1.Get("hmiss_1");

// Functions to fit Q weight plots
auto* func4=new TF1("func4","gaus(0)+gaus(3)+pol2(6)",0.36,0.65);
auto* func5=new TF1("func5","gaus(0)",0.36,0.65);
auto* func6=new TF1("func6","gaus(0)",0.36,0.65);
auto* func8=new TF1("func8","[0]*exp(-pow(x-[1],2)/(2*pow([2],2)))",0.36,0.65);
auto* func7=new TF1("func7","pol2(0)",0.36,0.65);


// Setting parameters before fitting for 2nd order polynomial and gaus
func4->SetParameter(0,hmass_kp->GetMaximum());
func4->FixParameter(1,0.49368);
func4->SetParameter(2,0.05);
func4->SetParameter(3,hmass_kp->GetMaximum() / 5.0);
func4->FixParameter(4,0.49368);
func4->SetParameter(5,0.02);
func4->SetParameter(6,hmass_kp->GetMaximum() / 2.0);
func4->SetParameter(7,-1);
func4->SetParameter(8,0);

// Making sure the amplitude is positive to avoid negative Q weight values
func4->SetParLimits(0,0,100000);
func4->SetParLimits(3,0,100000);
//
func4->SetLineColor(kRed);
hmass_kp->Fit("func4","R");



func5->FixParameter(0,func4->GetParameter(0));
func5->FixParameter(1,func4->GetParameter(1));
func5->FixParameter(2,func4->GetParameter(2));
func6->FixParameter(0,func4->GetParameter(3));
func6->FixParameter(1,func4->GetParameter(4));
func6->FixParameter(2,func4->GetParameter(5));
func7->FixParameter(0,func4->GetParameter(6));
func7->FixParameter(1,func4->GetParameter(7));
func7->FixParameter(2,func4->GetParameter(8));

cout<<func4->GetParameter(2) / func4->GetParameter(5)<<endl;

func5->SetLineColor(kBlue);
func6->SetLineColor(kGreen);
func7->SetLineColor(kOrange);

auto c1=new TCanvas("c1","",800,800);
hmass_kp->Draw();
func5->Draw("same");
func6->Draw("same");
func7->Draw("same");




}

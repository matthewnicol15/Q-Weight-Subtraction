{
TFile f1("/mnt/f/PhD/Analysis_Output/RGA/Skim4/Inbending/Strangeness_1/PID/Strangeness_1_RGA_Skim4_e_Kp_FD_ppi_no_FD_Inbending_190521_01.root");

TH1F *hmass_kp = (TH1F*)f1.Get("hmass_kp");
TH1F *hmiss = (TH1F*)f1.Get("hmiss_1");

// auto* func1=new TF1("func1","expo(0)",0.8,1.0);
auto* func1=new TF1("func1","pol2(0)+gaus(3)+gaus(6)",0.38,0.7);
// auto* func1=new TF1("func1","expo(0)+gaus(2)",0.4,0.7);
// auto* func1=new TF1("func1","pol1(0)",0.6,1.0);
// auto* func1=new TF1("func1","pol1(0)+gaus(2)",0.6,1.4);
auto* func2=new TF1("func2","pol2(0)",0.38,0.7);
// // auto* func2=new TF1("func2","pol4(0)",0.8,1.4);
auto* func3=new TF1("func3","gaus(0)",0.38,0.7);
auto* func4=new TF1("func4","gaus(0)",0.38,0.7);
//

func1->SetParameter(3,40000);
func1->SetParameter(4,0.5);
func1->SetParameter(5,0.02);
func1->SetParameter(6,5000);
func1->SetParameter(7,0.5);
func1->SetParameter(8,0.04);
//
func1->SetLineColor(kRed);
hmass_kp->Fit("func1","RB");
auto c1=new TCanvas("c1","",800,800);
hmass_kp->Draw();
// func1->Draw("same");

func2->FixParameter(0,func1->GetParameter(0));
func2->FixParameter(1,func1->GetParameter(1));
func2->FixParameter(2,func1->GetParameter(2));
// // func2->FixParameter(3,func1->GetParameter(3));
// // func2->FixParameter(4,func1->GetParameter(4));
func3->FixParameter(0,func1->GetParameter(3));
func3->FixParameter(1,func1->GetParameter(4));
func3->FixParameter(2,func1->GetParameter(5));
func4->FixParameter(0,func1->GetParameter(6));
func4->FixParameter(1,func1->GetParameter(7));
func4->FixParameter(2,func1->GetParameter(8));
//
func2->SetLineColor(kBlue);
func3->SetLineColor(kGreen);
func4->SetLineColor(kBlack);
func2->Draw("same");
func3->Draw("same");
func4->Draw("same");

// auto *c2=new TCanvas("c2",800,800);


}

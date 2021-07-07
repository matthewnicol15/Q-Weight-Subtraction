// Analysis split into bins of photon energy
// Select energy range below

Int_t Energy_Bin = 0 ; // Set to the energy bin you want to look at
// Photon Energy Bins [GeV]
// [0]  0-3,
// [1]  3-4
// [2]  4-4.5
// [3]  4.5-5
// [4]  5-5.5
// [5]  5.5-6
// [6]  6-6.5
// [7]  6.5-7
// [8]  7-7.5
// [9]  7.5-8
// [10]  8-8.5
// [11]  8.5-9
// [12]  9-9.5
// [13]  9.5-10.5

Int_t Nearest_Neighbour_Size; // Current size of vector of nearest neighbours

// Set sigma ratios here depending on photon energy bin
Double_t sigma_ratio_1 = 3.2, sigma_ratio_2 = 3.0, sigma_ratio_3 = 2.8, sigma_ratio_4 = 2.6;

#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>
#include <TCanvas.h>
#include <chrono>
#include <sstream>
#include <string>



// Macro name
void Strangeness_1_RGA_Q_Weight(){
  cout<<Energy_Bin<<endl;

  auto start = std::chrono::high_resolution_clock::now();

  // Read file with information on vectors
  // gROOT->ProcessLine(".L ~/work/Macros/Q_Weight/Loader.C+");
  gROOT->ProcessLine(".L /mnt/f/PhD/Macros/Loader.C+");

  std::ostringstream Input_File_Name, Output_File_Name; // String for root file names
  std::ostringstream Input_TTree_Name, Output_TTree_Name; // String for TTree names

  // Input_File_Name<<"/shared/storage/physhad/JLab/mn688/Trees/Dibaryon/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Split_test_160621_"<<Energy_Bin<<".root";
  Input_File_Name<<"/mnt/f/PhD/Trees/Dibaryon/RGA/Strangeness_1/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Split_test_160621_"<<Energy_Bin<<".root";
  Input_TTree_Name<<"Energy_"<<Energy_Bin;

  Output_File_Name<<"/mnt/f/PhD/Trees/Dibaryon/RGA/Strangeness_1/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Friend_060721_"<<Energy_Bin<<".root";
  // Output_File_Name<<"/shared/storage/physhad/JLab/mn688/Trees/Dibaryon/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Friend_060721_"<<Energy_Bin<<".root";
  Output_TTree_Name<<"Energy_"<<Energy_Bin;

  // Read input root file and assign it to 'f'
  TFile *fin = new TFile(Input_File_Name.str().c_str());
  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)fin->Get(Input_TTree_Name.str().c_str());

  // Recording the calculated mass of kaons for Q weights
  Double_t mass_kp, mass_kp_other;
  Int_t fitStatus1, fitStatus2;

  // Creating TTree friend
  // File to store TTree friend
  TFile *fout = new TFile(Output_File_Name.str().c_str(),"recreate");
  // TFile *ff = new TFile("/shared/storage/physhad/JLab/mn688/Trees/Dibaryon/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Friend_500_010621_04y.root","recreate");
  // TTree friend is copy of original TTree
  TTree t2("t2","it's a tree!");


  // // Functions to fit Q weight plots
  // auto* func1=new TF1("func1","gaus(0)+pol1(3)",0.36,0.65);
  // auto* func2=new TF1("func2","gaus(0)",0.36,0.65);
  // auto* func3=new TF1("func3","pol1(0)",0.36,0.65);

  // Functions to fit Q weight plots
  auto* func4=new TF1("func4","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow(3.2*[2],2))) + pol1(5)",0.36,0.65);
  auto* func5=new TF1("func5","gaus(0)",0.36,0.65);
  auto* func6=new TF1("func6","gaus(0)",0.36,0.65);
  auto* func7=new TF1("func7","pol1(0)",0.36,0.65);

  // Functions to fit Q weight plots
  auto* func8=new TF1("func8","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow(3.0*[2],2))) + pol1(5)",0.36,0.65);
  auto* func9=new TF1("func9","gaus(0)",0.36,0.65);
  auto* func10=new TF1("func10","gaus(0)",0.36,0.65);
  auto* func11=new TF1("func11","pol1(0)",0.36,0.65);

  // Functions to fit Q weight plots
  auto* func12=new TF1("func12","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow(2.8*[2],2))) + pol1(5)",0.36,0.65);
  auto* func13=new TF1("func13","gaus(0)",0.36,0.65);
  auto* func14=new TF1("func14","gaus(0)",0.36,0.65);
  auto* func15=new TF1("func15","pol1(0)",0.36,0.65);

  // Functions to fit Q weight plots
  auto* func16=new TF1("func16","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow(2.6*[2],2))) + pol1(5)",0.36,0.65);
  auto* func17=new TF1("func17","gaus(0)",0.36,0.65);
  auto* func18=new TF1("func18","gaus(0)",0.36,0.65);
  auto* func19=new TF1("func19","pol1(0)",0.36,0.65);

  // Integrals to calculate Q weights
  Double_t sig_1=-100, back_1=10, sig_2=-100, back_2=10, sig_3=-100, back_3=10, sig_4=-100, back_4=10;
  // Q weights
  Double_t Q_Weight_1a=-10, Q_Weight_1b=-10, Q_Weight_2a=-10, Q_Weight_3a=-10, Q_Weight_4a=-10;


  // Creating components to read from TTree
  // Set any vectors to 0 before reading from the TTree
  // Event information
  TLorentzVector *readbeam=NULL;  // Information on the beam
  TLorentzVector *readtarget=NULL; // Information on the target
  Double_t start_time; // Event start time
  Int_t readrunno, readrunno_2;
  Int_t readeventno, readeventno_2;
  Int_t readtriggerno;
  // Number of given particle or charge track in each event
  Int_t readchargetracks; // Number of positive or negative charge tracks
  Int_t readothertracks; // Number of particles excluding p, e^-, pions or kaons
  Int_t region; // which region the particles go in (FT, FD, CD)


  // Particle information
  vector<TLorentzVector> *v_p4=0;   // 4-vectors for the detected particles
  vector<TLorentzVector> *v_vertex=0;   // Vertex information for particles
  vector<double> *v_path=0;   // Measured path of particles
  vector<double> *v_time=0;   // Measured time of flight of particles
  vector<double> *v_beta=0;   // Beta measured from FTOF
  vector<double> *v_energy=0;   // Energy of particle
  vector<double> *v_charge=0;   // Charge of particle from Drift Chambers
  vector<double> *v_PID=0;   // PID from FTOF
  vector<double> *v_chi2PID=0;   // Chi^2 of the PID
  vector<Int_t> *v_region=0; // region particle goes in

  // Setting the branch addresses to read from original TTree
  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("vertex",&v_vertex);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);
  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("region",&v_region);
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);

  // Adding the Q Weights to the TTree friend
  t2.Branch("fitStatus1",&fitStatus1);
  t2.Branch("fitStatus2",&fitStatus2);
  t2.Branch("mass_kp",&mass_kp);
  t2.Branch("Q_Weight_1a",&Q_Weight_1a);
  t2.Branch("Q_Weight_1b",&Q_Weight_1b);
  t2.Branch("Q_Weight_2a",&Q_Weight_2a);
  t2.Branch("Q_Weight_3a",&Q_Weight_3a);
  t2.Branch("Q_Weight_4a",&Q_Weight_4a);
  t2.Branch("eventno",&readeventno_2);
  t2.Branch("runno",&readrunno_2);
  // t2.Branch("Q_Weight_2a",&Q_Weight_2a);
  // t2.Branch("Q_Weight_2b",&Q_Weight_2b);


  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  // Create histograms here
  auto* hmiss_mass_all=new TH1F("miss_all","MM^2(e' K^{+} p #pi^{-});MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",800,-1,3);
  auto* hmiss_momentum_all=new TH1F("hmiss_momentum_all","P(B + T - e' - K^{+} - p - #pi^{-});P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",800,-1,3);
  auto* hinv_lambda=new TH1F("hinv_lambda","Invariant mass of p #pi^{-};M(p #pi^{-}) [GeV];Counts",500,1,2);
  auto* hinv_missing_lambda=new TH2D("hinv_missing_lambda","Missing #Lambda against invariant #Lambda;M(p #pi^{-}) [GeV];MM(e' K^{+}) [GeV]",200,0.5,2.0,200,0.5,2.06);
  auto* hangular_distribution=new TH2D("hangular_distribution","Theta vs Phi for PID 0 assuming it's #pi^{-};#phi [deg];#theta [deg]",200,-200,200,200,-200,200);
  auto* hangular_distribution_momentum_reconstructed=new TH2D("hangular_distribution_momentum_reconstructed","Theta vs P for reconstructed #pi^{-};P [GeV];#theta [deg]",200,0,11,200,0,200);
  auto* hangular_distribution_momentum_detected=new TH2D("hangular_distribution_momentum_detected","Theta vs P for #pi^{-};P [GeV];#theta [deg]",200,0,11,200,0,200);
  auto* hmiss_1=new TH1F("hmiss_1","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_2=new TH1F("hmiss_1_2","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_sig=new TH1F("hmiss_1_sig","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_back=new TH1F("hmiss_1_back","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_phi=new TH1F("hmiss_1_phi","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_2=new TH1F("hmiss_2","MM^{2}(e' K^{+} p);MM^{2}(e' K^{+} p) [GeV^{2}];Counts",800,-1,3);
  auto* hmiss_3=new TH1F("hmiss_3","MM(e' p);MM(e' p) [GeV];Counts",800,-1,3);
  auto* hmass_kp=new TH1F("hmass_kp","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hkaon_pion_sig=new TH2F("hkaon_pion_sig","kaon vs pion;MM(e'K+) [GeV];MM(e' #pi^{+}) [GeV]",200,0,3,200,0,3);
  auto* hkaon_pion_lowback=new TH2F("hkaon_pion_lowback","kaon vs pion;MM(e'K+) [GeV];MM(e' #pi^{+}) [GeV]",200,0,3,200,0,3);
  auto* hkaon_pion_highback=new TH2F("hkaon_pion_highback","kaon vs pion;MM(e'K+) [GeV];MM(e' #pi^{+}) [GeV]",200,0,3,200,0,3);
  auto* h_delta_beta_kp=new TH2F("h_delta_beta_kp","",200,0,11,200,-1,1);
  auto* h_delta_beta_pr=new TH2F("h_delta_beta_pr","",200,0,11,200,-1,1);
  auto* h_delta_beta_pim=new TH2F("h_delta_beta_pim","",200,0,11,200,-1,1);
  auto* h_photon_energy=new TH1F("h_photon_energy","",400,0,11);
  auto* h_cos_theta_kaonp=new TH1F("h_cos_theta_kaonp","",400,-2,2);
  auto* h_mass_kp_Qweights_1a=new TH2D("h_mass_kp_Qweights_1a","Q weights missing mass",400,0.3,0.7,1000,0,1);
  auto* h_mass_kp_Qweights_1b=new TH2D("h_mass_kp_Qweights_1b","Q weights missing mass",400,0.3,0.7,1000,0,1);
  auto* h_mass_kp_Qweights_2a=new TH2D("h_mass_kp_Qweights_2a","Q weights missing mass",400,0.3,0.7,1000,0,1);
  auto* h_mass_kp_Qweights_2b=new TH2D("h_mass_kp_Qweights_2b","Q weights missing mass",400,0.3,0.7,1000,0,1);
  auto* hregion=new TH1F("hregion","Regions;Region;Counts",3,1,4);


  // Create vectors of TLorentzVectors to store information of
  // all particles of a given type (important when you have more than 1
  // of any particle in the same event)
  vector<TLorentzVector> v_el, v_el_other;  // e^-
  vector<TLorentzVector> v_pim, v_pim_other; // pi^-
  vector<TLorentzVector> v_pr, v_pr_other; // protons
  vector<TLorentzVector> v_kp, v_kp_other; // K^+
  vector<TLorentzVector> v_othertracks; // Any other particles are assigned to this

  // Q weight values
  TLorentzVector kaon_boost_com, kaon_boost_com_other; // kaon boosted in the COM reference frame
  Double_t Cos_Theta_Kp_COM, Cos_Theta_Kp_COM_other; // cos theta of kaons in COM frame
  Double_t photon_energy, photon_energy_other; // photon energy
  Double_t distance; // distance between missing mass points for Q weights
  Double_t cos_theta_range = 2, Momentum_range = 7; // Ranges of values for cos theta and photon energy
  vector<Double_t> v_NN_d, v_NN_mKp; // vectors containing information on the nearest neighbours

  // TLorentzVectors for individual particles
  TLorentzVector el, el_other;
  TLorentzVector pim, pim_other;
  TLorentzVector pr, pr_other;
  TLorentzVector kp, kp_other;

  // Vertex information for different particle types
  vector<TLorentzVector> v_vertex_el; // e^-
  vector<TLorentzVector> v_vertex_pim; // pi^-
  vector<TLorentzVector> v_vertex_pr; // protons
  vector<TLorentzVector> v_vertex_kp; // K^+
  TLorentzVector vertex_el;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;

  // These are used to define the missing masses later
  TLorentzVector missall;
  TLorentzVector miss1, miss1_other;
  TLorentzVector misspion;
  TLorentzVector miss2;
  TLorentzVector miss3;
  TLorentzVector lambda, Lambda_other;
  TLorentzVector photon, photon_other;
  TLorentzVector COM, COM_other;
  TVector3 COM_3, COM_3_other;


  // After information is read from the TTree, particles are identified using
  // PID and assigned to that particle type
  // e^-
  vector<Double_t> v_beta_tof_el;  // Beta measured
  Double_t beta_tof_el;
  vector<Double_t> v_P_el;  // Momentum measured
  Double_t P_el;
  vector<Double_t>v_path_el; // Path measured
  Double_t path_el;
  vector<Double_t>v_TOF_el;  // TOF measured
  Double_t TOF_el;
  vector<Double_t> v_beta_calc_el;  // Beta calculated
  Double_t beta_calc_el;
  vector<Double_t> v_delta_beta_el;  // Beta calculated - beta measured
  Double_t delta_beta_el;
  vector<Double_t> v_vertex_time_el;  // Vertex time calculated
  Double_t vertex_time_el;
  vector<Double_t> v_region_el; // region hit
  Double_t region_el;

  // pi^-
  vector<Double_t> v_beta_tof_pim;  // Beta measured
  Double_t beta_tof_pim;
  vector<Double_t> v_P_pim;  // Momentum measured
  Double_t P_pim;
  vector<Double_t>v_path_pim; // Path measured
  Double_t path_pim;
  vector<Double_t>v_TOF_pim;  // TOF measured
  Double_t TOF_pim;
  vector<Double_t> v_beta_calc_pim;  // Beta calculated
  Double_t beta_calc_pim;
  vector<Double_t> v_delta_beta_pim;  // Beta calculated - beta measured
  Double_t delta_beta_pim;
  vector<Double_t> v_vertex_time_pim;  // Vertex time calculated
  Double_t vertex_time_pim;
  vector<Double_t> v_region_pim; // region hit
  Double_t region_pim;

  // protons
  vector<Double_t> v_beta_tof_pr;  // Beta measured
  Double_t beta_tof_pr;
  vector<Double_t> v_P_pr;  // Momentum measured
  Double_t P_pr;
  vector<Double_t>v_path_pr; // Path measured
  Double_t path_pr;
  vector<Double_t>v_TOF_pr;  // TOF measured
  Double_t TOF_pr;
  vector<Double_t> v_beta_calc_pr;  // Beta calculated
  Double_t beta_calc_pr;
  vector<Double_t> v_delta_beta_pr;  // Beta calculated - beta measured
  Double_t delta_beta_pr;
  vector<Double_t> v_vertex_time_pr;  // Vertex time calculated
  Double_t vertex_time_pr;
  vector<Double_t> v_region_pr; // region hit
  Double_t region_pr;

  // K^+
  vector<Double_t> v_beta_tof_kp;  // Beta measured
  Double_t beta_tof_kp, beta_tof_kp_other;
  vector<Double_t> v_P_kp;  // Momentum measured
  Double_t P_kp;
  vector<Double_t>v_path_kp; // Path measured
  Double_t path_kp;
  vector<Double_t>v_TOF_kp;  // TOF measured
  Double_t TOF_kp;
  vector<Double_t> v_beta_calc_kp;  // Beta calculated
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;  // Beta calculated - beta measured
  Double_t delta_beta_kp;
  vector<Double_t> v_vertex_time_kp;  // Vertex time calculated
  Double_t vertex_time_kp;
  vector<Double_t> v_region_kp; // region hit
  Double_t region_kp;
  vector<Double_t> v_mass_kp; // region hit


  Double_t c=30;  // Speed of light used for calculating vertex time

  // Reads the total number of entries in the TTree
  Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  // Long64_t nentries = 1000000;
  cout<<t1->GetEntries()<<endl;
  // This is used to print out the percentage of events completed so far
  Int_t Percentage = nentries/100;
  // This loops over all the entries in the TTree
  for(Long64_t i = 0; i < nentries; i++){
    if((i/10)%10==0 && i%10==0)cout<<i<<endl;
    // Creating the histogram for kaon mass of nearest neighbours
    auto* h_miss1_NN=new TH1F("h_miss1_NN","NN missing mass",100,0.3,0.7);


    // Get this entry from the TTree
    t1->GetEntry(i);

    // Clear nearest neighbour vectors
    v_NN_d.clear();
    v_NN_mKp.clear();

    // Getting the run and event number for the second TTree
    readeventno_2 = readeventno;
    readrunno_2 = readrunno;

    // Setting Q weights each event
    Q_Weight_1a = -10;
    Q_Weight_1b = -10;
    // Q_Weight_2a = -10;
    // Q_Weight_2b = -10;

    // This prints out the percentage of events completed so far
    // if (i % Percentage == 0){
    //   fprintf (stderr, "%lld\r", i/Percentage);
    //   fflush (stderr);
    // }

    // All the vectors must be cleared at the start of each event entry
    // e^-
    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();
    v_region_el.clear();

    // pi^-
    v_pim.clear();
    v_beta_tof_pim.clear();
    v_P_pim.clear();
    v_path_pim.clear();
    v_TOF_pim.clear();
    v_beta_calc_pim.clear();
    v_delta_beta_pim.clear();
    v_vertex_pim.clear();
    v_vertex_time_pim.clear();
    v_vertex_pim.clear();
    v_region_pim.clear();

    // protons
    v_pr.clear();
    v_beta_tof_pr.clear();
    v_P_pr.clear();
    v_path_pr.clear();
    v_TOF_pr.clear();
    v_beta_calc_pr.clear();
    v_delta_beta_pr.clear();
    v_vertex_pr.clear();
    v_vertex_time_pr.clear();
    v_vertex_pr.clear();
    v_region_pr.clear();

    // K^+
    v_kp.clear();
    v_beta_tof_kp.clear();
    v_P_kp.clear();
    v_path_kp.clear();
    v_TOF_kp.clear();
    v_beta_calc_kp.clear();
    v_delta_beta_kp.clear();
    v_vertex_kp.clear();
    v_vertex_time_kp.clear();
    v_vertex_kp.clear();
    v_region_kp.clear();
    v_mass_kp.clear();


    // This reads the number of particles in the current entry
    Int_t Nparticles = v_p4->size();

    // This loops over all the particles in the current entry
    for(Int_t j=0; j<Nparticles; j++){

      // Filling histogram to show the different regions hit
      hregion->Fill(v_region->at(j));
      // Can put a selection on which regions particles are hitting

      // Checking PID and assigning particles
      // e^-
      if(v_PID->at(j)==11){
        // Setting the 4-vector and assigning mass from PDG
        el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(11)->Mass());
        TOF_el = v_time->at(j); // Measured time
        path_el = v_path->at(j); // Measured path
        beta_tof_el = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_el = el.Rho()/(sqrt((pow(el.Rho(),2))+(pow(el.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_el = beta_calc_el-beta_tof_el;
        vertex_time_el = TOF_el - path_el / (beta_tof_el*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);
        region_el = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_el.push_back(el);
        v_beta_tof_el.push_back(beta_tof_el);
        v_P_el.push_back(P_el);
        v_path_el.push_back(path_el);
        v_TOF_el.push_back(TOF_el);
        v_beta_calc_el.push_back(beta_calc_el);
        v_delta_beta_el.push_back(delta_beta_el);
        v_vertex_time_el.push_back(vertex_time_el);
        v_vertex_el.push_back(vertex_el);
        v_region_el.push_back(region_el);
      }

      // pi^-
      else if(v_PID->at(j)==-211){
        // Setting the 4-vector and assigning mass from PDG
        pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-211)->Mass());
        TOF_pim = v_time->at(j); // Measured time
        path_pim = v_path->at(j); // Measured path
        beta_tof_pim = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pim = pim.Rho()/(sqrt((pow(pim.Rho(),2))+(pow(pim.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pim = beta_calc_pim-beta_tof_pim;
        vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);
        region_pim = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pim.push_back(pim);
        v_beta_tof_pim.push_back(beta_tof_pim);
        v_P_pim.push_back(P_pim);
        v_path_pim.push_back(path_pim);
        v_TOF_pim.push_back(TOF_pim);
        v_beta_calc_pim.push_back(beta_calc_pim);
        v_delta_beta_pim.push_back(delta_beta_pim);
        v_vertex_time_pim.push_back(vertex_time_pim);
        v_vertex_pim.push_back(vertex_pim);
        v_region_pim.push_back(region_pim);
      }

      // protons
      else if(v_PID->at(j)==2212){
        // Setting the 4-vector and assigning mass from PDG
        pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2212)->Mass());
        TOF_pr = v_time->at(j); // Measured time
        path_pr = v_path->at(j); // Measured path
        beta_tof_pr = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pr = pr.Rho()/(sqrt((pow(pr.Rho(),2))+(pow(pr.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pr = beta_calc_pr-beta_tof_pr;
        vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);
        region_pr = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pr.push_back(pr);
        v_beta_tof_pr.push_back(beta_tof_pr);
        v_P_pr.push_back(P_pr);
        v_path_pr.push_back(path_pr);
        v_TOF_pr.push_back(TOF_pr);
        v_beta_calc_pr.push_back(beta_calc_pr);
        v_delta_beta_pr.push_back(delta_beta_pr);
        v_vertex_time_pr.push_back(vertex_time_pr);
        v_vertex_pr.push_back(vertex_pr);
        v_region_pr.push_back(region_pr);
      }

      // K^+
      else if(v_PID->at(j)==321){
        // Setting the 4-vector and assigning mass from PDG
        kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());

        TOF_kp = v_time->at(j); // Measured time
        path_kp = v_path->at(j); // Measured path
        beta_tof_kp = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_kp = kp.Rho()/(sqrt((pow(kp.Rho(),2))+(pow(kp.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_kp = beta_calc_kp-beta_tof_kp;
        vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
        region_kp = v_region->at(j);
        mass_kp = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_kp.push_back(kp);
        v_beta_tof_kp.push_back(beta_tof_kp);
        v_P_kp.push_back(P_kp);
        v_path_kp.push_back(path_kp);
        v_TOF_kp.push_back(TOF_kp);
        v_beta_calc_kp.push_back(beta_calc_kp);
        v_delta_beta_kp.push_back(delta_beta_kp);
        v_vertex_time_kp.push_back(vertex_time_kp);
        v_vertex_kp.push_back(vertex_kp);
        v_region_kp.push_back(region_kp);
        v_mass_kp.push_back(mass_kp);
      }
    }

    // Here you can apply conditions on the events you want to analyse
    if(v_kp.size()==1 && v_el.size()==1 && v_pr.size()==1 && v_pim.size()==1){
      // cout<<"mass "<<v_mass_kp.at(0)<<endl;
      // If outside interested mass range then save entry but don't bother fitting
      if(v_mass_kp.at(0) > 0.65){
        Q_Weight_1a = 0;
        Q_Weight_1b = 0;
        Q_Weight_2a = 0;
        Q_Weight_3a = 0;
        Q_Weight_4a = 0;

        fitStatus1 = -1;
        fitStatus2 = -1;
      }

      else{

        // Select which region you want the particles to go in
        // if(v_region_kp.at(0) == 1 && v_region_pr.at(0) == 1 && v_region_pim.at(0) == 1){

        // Missing mass of e' K^{+}, looking for lambda ground state
        miss1 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);
        // Missing mass of e' K^{+} p, looking for pi^{-}
        miss2 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0);
        // Missing mass of e', looking for phi meson background
        miss3 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pr.at(0);
        // Invariant mass of p pi^{-}, looking for lambda ground state
        lambda = v_pr.at(0) + v_pim.at(0);
        // Missing mass of all detected particles, should have peak at 0
        missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0) - v_pim.at(0);
        // Determining the photon TLorentzVector
        photon = (TLorentzVector)*readbeam - v_el.at(0);
        // COM used for boosting
        COM = v_el.at(0) + v_kp.at(0) + v_pr.at(0) + v_pim.at(0);

        // Filling missing mass histograms
        hmass_kp->Fill(mass_kp);
        hmiss_1->Fill(miss1.M());
        hmiss_2->Fill(miss2.M2());
        hmiss_3->Fill(miss3.M());

        // Checking the delta beta for each particle
        h_delta_beta_kp->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
        h_delta_beta_pr->Fill(v_pr.at(0).Rho(),v_delta_beta_pr.at(0));
        h_delta_beta_pim->Fill(v_pim.at(0).Rho(),v_delta_beta_pim.at(0));


        // Boosting the kaon in COM reference frame
        COM_3 = COM.BoostVector();
        kaon_boost_com = v_kp.at(0);
        kaon_boost_com.Boost(-1.0*COM_3);

        // Calculating cos theta for boosted kaon
        Cos_Theta_Kp_COM = cos(kaon_boost_com.Theta());

        // Getting the photon energy
        photon_energy = photon.E();


        // Plotting the cos(theta) distribution of K+
        h_cos_theta_kaonp->Fill(Cos_Theta_Kp_COM);

        // Plotting the photon energy
        h_photon_energy->Fill(photon_energy);

        // Looking at the angular distribution of detected pi^{-}
        hangular_distribution_momentum_detected->Fill(v_pim.at(0).Rho(), v_pim.at(0).Theta()*TMath::RadToDeg());
        // Invariant against missing mass of lambda
        hinv_missing_lambda->Fill(lambda.M(),miss1.M());
        // Filling invariant and missing mass histograms
        hinv_lambda->Fill(lambda.M());
        hmiss_mass_all->Fill(missall.M2());
        hmiss_momentum_all->Fill(missall.Rho());

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Determining nearest neighbours for Q weight calculations

        // Only looking at nearest neighbours for events with good missing mass values
        // if(miss1.M() > 0.7 && miss1.M() < 1.6){

        Int_t size=0;
        // Loop over all other events
        for(Long64_t m = 0; m < nentries; m++){
          if(m==i)continue;

          Nearest_Neighbour_Size = v_NN_d.size();


          // Gets information on event m
          t1->GetEntry(m);

          v_el_other.clear();
          v_kp_other.clear();
          v_pr_other.clear();
          v_pim_other.clear();

          Int_t Nparticles_other = v_p4->size();

          // This loops over all the particles in the current entry
          for(Int_t p=0; p<Nparticles_other; p++){

            if(v_PID->at(p) == 11){
              el_other.SetXYZM(v_p4->at(p).Px(), v_p4->at(p).Py(), v_p4->at(p).Pz(), db->GetParticle(11)->Mass());
              v_el_other.push_back(el_other);
            }

            else if(v_PID->at(p) == 321){
              kp_other.SetXYZM(v_p4->at(p).Px(), v_p4->at(p).Py(), v_p4->at(p).Pz(), db->GetParticle(321)->Mass());
              v_kp_other.push_back(kp_other);
              beta_tof_kp_other = v_beta->at(p);
            }

            else if(v_PID->at(p) == 2212){
              pr_other.SetXYZM(v_p4->at(p).Px(), v_p4->at(p).Py(), v_p4->at(p).Pz(), db->GetParticle(2212)->Mass());
              v_pr_other.push_back(pr_other);
            }

            else if(v_PID->at(p) == -211){
              pim_other.SetXYZM(v_p4->at(p).Px(), v_p4->at(p).Py(), v_p4->at(p).Pz(), db->GetParticle(-211)->Mass());
              v_pim_other.push_back(pim_other);
            }
          }

          // Checking for events with 1 K+, 1 e-, 1 proton and 1 pi-
          if(v_el_other.size() == 1 && v_kp_other.size() == 1 && v_pr_other.size() == 1 && v_pim_other.size() == 1){

            // Invariant mass of proton pi^-
            Lambda_other = v_pr_other.at(0) + v_pim_other.at(0);

            // Checking if event has invariant mass close to that of ground state lambda
            // if(Lambda_other.M() < 1.14){

            // Missing mass of e' K^{+}, looking for lambda ground state
            miss1_other = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el_other.at(0) - v_kp_other.at(0);
            // Determining the four vector for the photon
            photon_other = (TLorentzVector)*readbeam - v_el_other.at(0);
            // Determining the COM frame of reference
            COM_other = v_el_other.at(0) + v_kp_other.at(0) + v_pim_other.at(0) + v_pr_other.at(0);

            mass_kp_other = sqrt((pow(v_kp_other.at(0).Rho(),2) / (pow(beta_tof_kp_other,2))) - pow(v_kp_other.at(0).Rho(),2));

            // cout<<v_kp_other.size()<<" "<<v_beta->size()<<endl;
            if(mass_kp_other > 0.65)continue;

            // Making the 3 vector for COM frame of reference
            COM_3_other = COM_other.BoostVector();
            // Getting the kaon information before boosting
            kaon_boost_com_other = v_kp_other.at(0);
            // Boosting the other kaon into COM frame of reference
            kaon_boost_com_other.Boost(-1.0*COM_3_other);

            // Get cos theta and gamma energy here
            Cos_Theta_Kp_COM_other = cos(kaon_boost_com_other.Theta());
            photon_energy_other = photon_other.E();

            // Calculate difference between cos thetas
            distance = (pow((Cos_Theta_Kp_COM - Cos_Theta_Kp_COM_other) / cos_theta_range,2)) /*+
            (pow((v_kp.at(0).Rho() - v_kp_other.at(0).Rho()) / Momentum_range,2))*/;


            // If there are no entries yet then just push back the values
            if(Nearest_Neighbour_Size<1){
              v_NN_d.push_back(distance);
              v_NN_mKp.push_back(mass_kp_other);
            }
            // If there are entries
            else{

              // Loop over current stored nearest neighbours, starting at largest distance
              for(Int_t k = Nearest_Neighbour_Size - 1; k >= 0; k --){

                // Skip events with distance greater than the max in vector if there are already 500
                if(Nearest_Neighbour_Size == 500 && distance > v_NN_d.at(499)) break;



                // Creating iterator for beginning of vector for smallest distance
                vector<double>::iterator itPos2 = v_NN_d.begin();
                vector<double>::iterator itPos2MM = v_NN_mKp.begin();

                if(distance < v_NN_d.at(0)){

                  // Inserting values into vector of nearest neighbours
                  v_NN_d.insert(itPos2, distance);
                  v_NN_mKp.insert(itPos2MM, mass_kp_other);

                  // Removing the largest value if there are over 500 in the vector
                  if(Nearest_Neighbour_Size==501){
                    v_NN_d.pop_back();
                    v_NN_mKp.pop_back();
                  }

                  break;
                }



                // Checking if distance is larger than the current value in vector
                if(distance > v_NN_d.at(k)){

                  // Creating iterator for the position to insert distance into vector
                  vector<double>::iterator itPos = v_NN_d.begin() + k + 1;
                  vector<double>::iterator itPosMM = v_NN_mKp.begin() + k + 1;

                  // Pushing back values to the next position
                  v_NN_d.insert(itPos, distance);
                  v_NN_mKp.insert(itPosMM, mass_kp_other);

                  // Removing the largest value if there are over 500 in the vector
                  if(Nearest_Neighbour_Size==501){
                    v_NN_d.pop_back();
                    v_NN_mKp.pop_back();
                  }

                  break;
                }
              }
            }
          }
          // }
          // }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Plotting and fitting the nearest neighbours

        // If there are values for the nearest neighbours
        if(Nearest_Neighbour_Size > 0){
          Int_t k = Nearest_Neighbour_Size;
          // Loop over all the nearest neighbour values
          for(Int_t q = 0; q < k; q++){
            // Fill the histogram with the missing mass values for the nearest neighbours
            h_miss1_NN->Fill(v_NN_mKp.at(q));
            // cout<<v_NN_mKp.at(q)<<endl;
          }

          // Setting parameters before fitting for 3rd order polynomial and gaus
          // func1->SetParameter(0,h_miss1_NN->GetMaximum());
          // func1->FixParameter(1,0.49368);
          // func1->SetParameter(2,0.025);
          // func1->SetParameter(3,h_miss1_NN->GetMaximum() / 2.0);
          // func1->SetParameter(4,-1);
          // func1->SetParameter(5,0);
          // func1->SetParameter(6,0);
          //
          // // Making sure the amplitude is positive to avoid negative Q weight values
          // func1->SetParLimits(0,h_miss1_NN->GetMaximum() / 3.0,h_miss1_NN->GetMaximum() * 1.2);
          // func1->SetParLimits(2,0.01, 0.09);
          // func1->SetParLimits(5,0, 100);

          // Setting parameters before fitting for 2nd order polynomial and gaus
          func4->SetParameter(0,h_miss1_NN->GetMaximum());
          func4->FixParameter(1,0.49368);
          func4->SetParameter(2,0.01);
          func4->SetParameter(3,h_miss1_NN->GetMaximum() / 5.0);
          func4->FixParameter(4,0.49368);
          func4->SetParameter(5,h_miss1_NN->GetMaximum() / 2.0);
          func4->SetParameter(6,-1);

          // Making sure the amplitude is positive to avoid negative Q weight values
          func4->SetParLimits(0,0.5,60);
          // func4->SetParLimits(2,0.001,0.012);
          func4->SetParLimits(3,0.5,60);

          // Setting parameters before fitting for 2nd order polynomial and gaus
          func8->SetParameter(0,h_miss1_NN->GetMaximum());
          func8->FixParameter(1,0.49368);
          func8->SetParameter(2,0.01);
          func8->SetParameter(3,h_miss1_NN->GetMaximum() / 5.0);
          func8->FixParameter(4,0.49368);
          func8->SetParameter(5,h_miss1_NN->GetMaximum() / 2.0);
          func8->SetParameter(6,-1);

          // Making sure the amplitude is positive to avoid negative Q weight values
          func8->SetParLimits(0,0.5,60);
          // func8->SetParLimits(2,0.001,0.012);
          func8->SetParLimits(3,0.5,60);

          // Setting parameters before fitting for 2nd order polynomial and gaus
          func12->SetParameter(0,h_miss1_NN->GetMaximum());
          func12->FixParameter(1,0.49368);
          func12->SetParameter(2,0.01);
          func12->SetParameter(3,h_miss1_NN->GetMaximum() / 5.0);
          func12->FixParameter(4,0.49368);
          func12->SetParameter(5,h_miss1_NN->GetMaximum() / 2.0);
          func12->SetParameter(6,-1);

          // Making sure the amplitude is positive to avoid negative Q weight values
          func12->SetParLimits(0,0.5,60);
          // func12->SetParLimits(2,0.001,0.012);
          func12->SetParLimits(3,0.5,60);

          // Setting parameters before fitting for 2nd order polynomial and gaus
          func16->SetParameter(0,h_miss1_NN->GetMaximum());
          func16->FixParameter(1,0.49368);
          func16->SetParameter(2,0.01);
          func16->SetParameter(3,h_miss1_NN->GetMaximum() / 5.0);
          func16->FixParameter(4,0.49368);
          func16->SetParameter(5,h_miss1_NN->GetMaximum() / 2.0);
          func16->SetParameter(6,-1);

          // Making sure the amplitude is positive to avoid negative Q weight values
          func16->SetParLimits(0,0.5,60);
          // func16->SetParLimits(2,0.001,0.012);
          func16->SetParLimits(3,0.5,60);

          h_miss1_NN->Sumw2();
          fitStatus1 =  h_miss1_NN->Fit("func4","RQ");
          fitStatus1 =  h_miss1_NN->Fit("func8","RQ");
          fitStatus1 =  h_miss1_NN->Fit("func12","RQ");
          fitStatus1 =  h_miss1_NN->Fit("func16","RQ");

          func5->FixParameter(0,func4->GetParameter(0));
          func5->FixParameter(1,func4->GetParameter(1));
          func5->FixParameter(2,func4->GetParameter(2));
          func6->FixParameter(0,func4->GetParameter(3));
          func6->FixParameter(1,func4->GetParameter(4));
          func6->FixParameter(2,func4->GetParameter(2) * 3.2);
          func7->FixParameter(0,func4->GetParameter(5));
          func7->FixParameter(1,func4->GetParameter(6));

          func9->FixParameter(0,func8->GetParameter(0));
          func9->FixParameter(1,func8->GetParameter(1));
          func9->FixParameter(2,func8->GetParameter(2));
          func10->FixParameter(0,func8->GetParameter(3));
          func10->FixParameter(1,func8->GetParameter(4));
          func10->FixParameter(2,func8->GetParameter(2) * 3.0);
          func11->FixParameter(0,func8->GetParameter(5));
          func11->FixParameter(1,func8->GetParameter(6));

          func13->FixParameter(0,func12->GetParameter(0));
          func13->FixParameter(1,func12->GetParameter(1));
          func13->FixParameter(2,func12->GetParameter(2));
          func14->FixParameter(0,func12->GetParameter(3));
          func14->FixParameter(1,func12->GetParameter(4));
          func14->FixParameter(2,func12->GetParameter(2) * 2.8);
          func15->FixParameter(0,func12->GetParameter(5));
          func15->FixParameter(1,func12->GetParameter(6));

          func17->FixParameter(0,func16->GetParameter(0));
          func17->FixParameter(1,func16->GetParameter(1));
          func17->FixParameter(2,func16->GetParameter(2));
          func18->FixParameter(0,func16->GetParameter(3));
          func18->FixParameter(1,func16->GetParameter(4));
          func18->FixParameter(2,func16->GetParameter(2) * 2.6);
          func19->FixParameter(0,func16->GetParameter(5));
          func19->FixParameter(1,func16->GetParameter(6));

          // Calculating the signal and background for 1st fit
          sig_1 = func5->Eval(mass_kp) + func6->Eval(mass_kp);
          back_1 = func7->Eval(mass_kp);
          // Calculating the signal and background for 1st fit
          sig_2 = func9->Eval(mass_kp) + func10->Eval(mass_kp);
          back_2 = func11->Eval(mass_kp);
          // Calculating the signal and background for 1st fit
          sig_3 = func13->Eval(mass_kp) + func14->Eval(mass_kp);
          back_3 = func15->Eval(mass_kp);
          // Calculating the signal and background for 1st fit
          sig_4 = func17->Eval(mass_kp) + func18->Eval(mass_kp);
          back_4 = func19->Eval(mass_kp);

          // Calculating Q weight for 1st fit using 2 methods
          Q_Weight_1a = sig_1 / (sig_1 + back_1);
          Q_Weight_1b = 1 - (back_1 / h_miss1_NN->GetBinContent(h_miss1_NN->FindBin(mass_kp)));
          // Calculating Q weight for 1st fit using 2 methods
          Q_Weight_2a = sig_2 / (sig_2 + back_2);
          // Q_Weight_2b = 1 - (back_2 / h_miss1_NN->GetBinContent(h_miss1_NN->FindBin(mass_kp)));
          // Calculating Q weight for 1st fit using 2 methods
          Q_Weight_3a = sig_3 / (sig_3 + back_3);
          // Q_Weight_3b = 1 - (back_3 / h_miss1_NN->GetBinContent(h_miss1_NN->FindBin(mass_kp)));
          // Calculating Q weight for 1st fit using 2 methods
          Q_Weight_4a = sig_4 / (sig_4 + back_4);
          // Q_Weight_4b = 1 - (back_4 / h_miss1_NN->GetBinContent(h_miss1_NN->FindBin(mass_kp)));

        }

        // Setting Q weight to 0 for events without "good" missing masses
        else{
          Q_Weight_1a = 0;
          Q_Weight_1b = 0;
          Q_Weight_2a = 0;
          // Q_Weight_2b = 0;
          Q_Weight_3a = 0;
          Q_Weight_4a = 0;
        }

      }
      h_mass_kp_Qweights_1a->Fill(mass_kp,Q_Weight_1a);
      h_mass_kp_Qweights_1b->Fill(mass_kp,Q_Weight_1b);
      // h_mass_kp_Qweights_2a->Fill(mass_kp,Q_Weight_2a);
      // h_mass_kp_Qweights_2b->Fill(mass_kp,Q_Weight_2b);
      // if(Q_Weight_1a < 0) cout<<"Bad! "<<mass_kp<<" "<<Q_Weight_1a<<" "<<sig_1<<" "<<back_1<<endl;
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      if(lambda.M() < 1.14){
        hmiss_1_2->Fill(miss1.M());


        if(miss3.M()> 0.9 && miss3.M() < 1.1) hmiss_1_phi->Fill(miss1.M());
        // Looking at the sidebands from the mass of kaons
        if(mass_kp > 0.454404 && mass_kp < 0.535267){
          hmiss_1_sig->Fill(miss1.M());
          hkaon_pion_sig->Fill(miss1.M(),misspion.M());
        }

        else if(mass_kp > 0.575698 && mass_kp < 0.616129){
          hmiss_1_back->Fill(miss1.M());
          hkaon_pion_highback->Fill(miss1.M(),misspion.M());

        }
        else if(mass_kp > 0.373542 && mass_kp < 0.413973){
          hmiss_1_back->Fill(miss1.M());
          hkaon_pion_lowback->Fill(miss1.M(),misspion.M());

        }

      }

      // } // Selecting events with kaons and protons hitting FD
    } // Selecting events with 1 e, 1 K^{+} and 1 p

    // auto *c1=new TCanvas("c1","",800,800);
    // h_miss1_NN->Draw();
    // func1->SetLineColor(kYellow);
    // func2->SetLineColor(kBlack);
    // func3->SetLineColor(kOrange);
    // func1->Draw("same");
    // func2->Draw("same");
    // func3->Draw("same");
    //
    // auto *c2=new TCanvas("c2","",800,800);
    // h_miss1_NN->Draw();
    // func5->SetLineColor(kBlue);
    // func6->SetLineColor(kBlack);
    // func7->SetLineColor(kGreen);
    // // func4->Draw("same");
    // func5->Draw("same");
    // func6->Draw("same");
    // func7->Draw("same");
    //
    // auto *c3=new TCanvas("c3","",800,800);
    // h_miss1_NN->Draw();
    // func9->SetLineColor(kBlue);
    // func10->SetLineColor(kBlack);
    // func11->SetLineColor(kGreen);
    // // func4->Draw("same");
    // func9->Draw("same");
    // func10->Draw("same");
    // func11->Draw("same");
    //
    // auto *c4=new TCanvas("c4","",800,800);
    // h_miss1_NN->Draw();
    // func13->SetLineColor(kBlue);
    // func14->SetLineColor(kBlack);
    // func15->SetLineColor(kGreen);
    // // func4->Draw("same");
    // func13->Draw("same");
    // func14->Draw("same");
    // func15->Draw("same");
    //
    // auto *c5=new TCanvas("c5","",800,800);
    // h_miss1_NN->Draw();
    // func17->SetLineColor(kBlue);
    // func18->SetLineColor(kBlack);
    // func19->SetLineColor(kGreen);
    // // func4->Draw("same");
    // func17->Draw("same");
    // func18->Draw("same");
    // func19->Draw("same");

    // auto *c3=new TCanvas("c3","",800,800);
    // h_distance->Draw();

    h_miss1_NN->Delete();
    t2.Fill();
  } // Event loop



  // Save root file with TTree friend
  fout->Write();
  // Close file with TTree friend
  fout->Close();

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count();
}

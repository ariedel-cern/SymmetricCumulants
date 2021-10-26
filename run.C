/**
 * File              : run.C
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 26.10.2021
 * Last Modified Date: 26.10.2021
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

R__ADD_INCLUDE_PATH($GRID_UTILITY_SCRIPTS)

Int_t run() {

  TFile *file = new TFile("SymmetricCumulants.root", "RECREATE");

  // pdf for multiplicity flucutations
  Int_t NumberOfParticles0 = 900;
  Int_t NumberOfParticles1 = 1100;
  TF1 *PDFmultiplicity =
      new TF1("Multiplicity", "1", NumberOfParticles0, NumberOfParticles1);
  PDFmultiplicity->Write();

  // pdf for phi distribution
  // set flow harmonics
  const Int_t N = 6;
  std::vector<Double_t> harmonics = {};
  for (int i = 0; i < N; ++i) {
    harmonics.push_back(0.05 + i * 0.01);
  }
  // write the formula
  // f(phi) = 1+2*Sum_i^N v_i cos(i*(phi-psi_i))
  TString formula = "1+2*(";

  for (int i = 0; i < N; i++) {
    formula += Form("[%d]*TMath::Cos(%d*(x-[%d]))", 2 * i, i + 1, 2 * i + 1);
    if (i != N - 1) {
      formula += TString("+");
    }
  }

  formula += ")";

  // std::cout << formula << std::endl;

  Double_t phi0 = 0;
  Double_t phi1 = TMath::TwoPi();
  Int_t binsPhi = 360;
  TF1 *PDFphi = new TF1("phi", formula, phi0, phi1);

  // set flow harmonics
  // symmetry planes will be set run by run
  for (std::size_t i = 0; i < harmonics.size(); i++) {
    PDFphi->SetParameter(2 * i, harmonics.at(i));
  }
  PDFphi->Write();

  // symmetric cumulants we want to compute
  std::vector<std::vector<Int_t>> SC = {{2, 3}, {2, 4}, {2, 3, 4}};

  // custom seed
  UInt_t Seed = 2020201;

  // configure task for uniform acceptance
  AliAnalysisTaskAR *task = new AliAnalysisTaskAR("SymmetricCumulant");
  task->SetMCOnTheFly(kTRUE);
  task->SetCustomSeed(Seed);
  task->SetSymmetricCumulants(SC);
  task->SetTrackControlHistogramBinning(kPHI, binsPhi, phi0, phi1);
  task->SetEventControlHistogramBinning(kMULQ,
                                        NumberOfParticles1 - NumberOfParticles0,
                                        NumberOfParticles0, NumberOfParticles1);

  // set pdfs
  task->SetMCMultiplicityPdf(PDFmultiplicity);

  // add all tasks to the analysis manager in a loop
  std::vector<AliAnalysisTaskAR *> tasks = {task};

  // connect to some output to silence error messages
  AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  TString OutputFile("AnalysisResults.root:AnalysisResults");
  // AliAnalysisDataContainer *cinput = nullptr;
  AliAnalysisDataContainer *coutput = nullptr;

  // loop over all tasks
  for (auto T : tasks) {
    mgr->AddTask(T);
    cout << "Added to manager: " << T->GetName() << endl;
    // cinput = mgr->GetCommonInputContainer();
    coutput =
        mgr->CreateContainer(T->GetName(), TList::Class(),
                             AliAnalysisManager::kOutputContainer, OutputFile);
    // mgr->ConnectInput(T, 0, cinput);
    mgr->ConnectOutput(T, 1, coutput);
  }

  mgr->InitAnalysis();

  // run analysis
  task->UserCreateOutputObjects();

  // loop over all events
  Int_t NumberOfEvents = 10000;
  Double_t Psi = 0.;
  for (int i = 0; i < NumberOfEvents; ++i) {

    std::cout << "Processing Event " << i << std::endl;

    // randomize psi
    Psi = gRandom->Uniform(TMath::TwoPi());
    for (std::size_t i = 0; i < harmonics.size(); i++) {
      PDFphi->SetParameter(2 * i + 1, Psi);
    }
    task->SetMCKinematicPdf(kPHI, PDFphi);

    task->UserExec(nullptr);
  }
  task->Terminate(nullptr);

  task->GetEventControlHistogram(kSIM, kMULQ, kAFTER)->Write();
  task->GetTrackControlHistogram(kSIM, kPHI, kBEFORE)->Write();

  // std::vector<Double_t> output;
  // task->GetSymmetricCumulantValues(&output);
  // for (auto v : output) {
  //   std::cout << v << std::endl;
  // }
  TList *cor = task->GetFinalResultProfilesList();
  cor->Write();

  TList *sc = task->GetFinalResultSymmetricCumulantsList();
  sc->Write();

  file->Close();
  return 0;
}

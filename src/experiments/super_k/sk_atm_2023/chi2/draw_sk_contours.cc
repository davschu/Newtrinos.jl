{
    TCanvas* c = new TCanvas("c", "SK Results 2023", 800, 400);
    c->Divide(2,1);

    TFile* f = TFile::Open("sk_q13-constrained_2023.root", "read");
    /*
        KEY: TGraph   SK2023_q13-constrained_dcp_chi2_no;1    Graph
        KEY: TGraph   SK2023_q13-constrained_dcp_chi2_io;1    Graph
        KEY: TGraph   SK2023_q13-constrained_m23_s23_68percent_no;1   Graph
        KEY: TGraph   SK2023_q13-constrained_m23_s23_68percent_io;1   Graph
        KEY: TGraph   SK2023_q13-constrained_m23_s23_90percent_no;1   Graph
        KEY: TGraph   SK2023_q13-constrained_m23_s23_90percent_io;1   Graph
    */

    // Left figure: delta_{CP} Delta chi^2
    c->cd(1);
    TGraph* gCP_no = (TGraph*)f->Get("SK2023_q13-constrained_dcp_chi2_no");
    TGraph* gCP_io = (TGraph*)f->Get("SK2023_q13-constrained_dcp_chi2_io");
    gCP_no->SetLineColor(434);
    gCP_io->SetLineColor(802);
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gCP_no);
    mg->Add(gCP_io);
    mg->Draw("al");
    mg->SetTitle("SK Results 2023;#delta_{CP};#Delta#chi^{2}");

    TLegend* l = new TLegend(0.1, 0.8, 0.4, 0.9);
    l->AddEntry(gCP_no, "Normal", "L");
    l->AddEntry(gCP_io, "Inverted", "L");
    l->Draw();

    // Right figure: atmospheric mixing parameter contours 
    c->cd(2);
    TGraph* gAtm_68_no = (TGraph*)f->Get("SK2023_q13-constrained_m23_s23_68percent_no");
    TGraph* gAtm_68_io = (TGraph*)f->Get("SK2023_q13-constrained_m23_s23_68percent_io");
    TGraph* gAtm_90_no = (TGraph*)f->Get("SK2023_q13-constrained_m23_s23_90percent_no");
    TGraph* gAtm_90_io = (TGraph*)f->Get("SK2023_q13-constrained_m23_s23_90percent_io");
    gAtm_68_no->SetLineColor(434);
    gAtm_68_io->SetLineColor(802);
    gAtm_68_no->SetLineStyle(2);
    gAtm_68_io->SetLineStyle(2);
    gAtm_90_no->SetLineColor(434);
    gAtm_90_io->SetLineColor(802);
    TMultiGraph* mg2 = new TMultiGraph();
    mg2->Add(gAtm_68_no);
    mg2->Add(gAtm_68_io);
    mg2->Add(gAtm_90_no);
    mg2->Add(gAtm_90_io);
    mg2->Draw("al");
    mg2->SetTitle("SK Results 2023;sin^{2}#theta_{23};|#Delta m^{2}_{32,31}| [eV^{2}]");

    TLegend* l2 = new TLegend(0.1, 0.8, 0.4, 0.9);
    l2->AddEntry(gAtm_90_no, "Normal", "L");
    l2->AddEntry(gAtm_90_io, "Inverted", "L");
    TLine* l68 = new TLine();
    l68->SetLineStyle(2);
    TLine* l90 = new TLine();
    l68->SetLineColor(kBlack);
    l90->SetLineColor(kBlack);
    l2->AddEntry(l68, "68%", "L");
    l2->AddEntry(l90, "90%", "L");
    l2->Draw();
}

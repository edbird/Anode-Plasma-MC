#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "HistogramWrapper.hpp"

class HistogramWrapper2
{

    public:

    HistogramWrapper2(TFile* tfile)
        : histogram{TH1F("name2", "name2", 10, 0.0, 1.0)}
        , tfile{tfile}
    {
        histogram.Fill(0.1);
    }

    ~HistogramWrapper2()
    {
        TCanvas *c = new TCanvas("cname2", "cname2", 800, 600);
        histogram.Draw();
        c->SaveAs("canvas2.png");
        delete c;

        tfile->cd();
        histogram.Write();
        tfile->Close();
    }

    TH1F histogram;
    TFile *tfile;

};

int main()
{
    TFile *tfile4 = new TFile("file4.root", "recreate");
    HistogramGroupProperties gp("name-base", "");
    HistogramProperties p2("name4", 10, 0.0, 1.0);
    HistogramGroupFloat g(tfile4, gp, p2);
    std::string address{g.Add("name4_2", 10, 0.0, 1.0)};
    g.Ref(address).Fill(0.1);
    g.Write();
    g.Canvas();

    TFile *tfile3 = new TFile("file3.root", "recreate");
    HistogramProperties p("name3", 10, 0.0, 1.0);
    HistogramWrapperFloat w(tfile3, p);
    w.Ref().Fill(0.1);
    w.Write();
    w.Canvas();

    TFile *tfile2 = new TFile("file2.root", "recreate");
    HistogramWrapper2 w2(tfile2);


    TFile *tfile = new TFile("file.root", "recreate");

    TH1F histogram{TH1F("name", "name", 10, 0.0, 1.0)};
    histogram.Fill(0.1);
    histogram.Write();

    TCanvas *c = new TCanvas("cname", "cname", 800, 600);
    histogram.Draw();
    c->SaveAs("canvas.png");
    delete c;


    tfile->Close();


}

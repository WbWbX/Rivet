#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MatrixN.hh"


namespace Rivet {


  /// Implements the cuts described in https://arxiv.org/pdf/1702.06996.pdf
  class CMS_WbWb_internal : public Analysis {
  public:

    /// Minimal constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_WbWb_internal);


    /// @name Analysis methods
    //@{


    /// Set up projections and book histograms
    void init() {

      // Cuts
      fs_cut       = (Cuts::abseta < 5) and (Cuts::pT > 0.0*MeV);
      lepton_cut   = (Cuts::abseta<2.4 && Cuts::pT > 25.*GeV);
      allj_cut     = (Cuts::pt>25*GeV && Cuts::abseta<4.7);
      bjet_cut     = (Cuts::pt>25*GeV && Cuts::abseta<2.5);

      // Complete final state
      FinalState fs(fs_cut);

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);


      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);

      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, lepton_cut, true);
      declare(dressed_leptons, "DressedLeptons");

      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      declare(FastJets(fsForJets, FastJets::ANTIKT, 0.4,
                       JetAlg::Muons::ALL, JetAlg::Invisibles::NONE), "Jets");

      //Projection for MET
      declare(MissingMomentum(fsForJets), "MissingET");
      
      //booking of histograms
      book(_cut_flow, "cut_flow", 10,0.5,10.5);
      book(histos["total"],    "total",         1, 0.0, 2.0);
      book(histos["wpt"],    "wpt",         40, 0.0, 280.0);
      book(histos["weta"],   "weta",        40,-5,5);
	  book(histos["njets"],  "njets",       6, 0.5, 6.5);
	  book(histos["nbjets"],  "nbjets",       4, 0.5, 4.5);
	  book(histos["nljets"],  "nljets",       6, -0.5, 5.5);
	  book(histos["bm1eta"], "bm1eta",        40,-5,5);
	  book(histos["bp1eta"], "bp1eta",        40,-5,5);
	  book(histos["b1pt"],   "b1pt",        40, 25.0, 225.0);
	  book(histos["b2pt"],   "b2pt",        40, 25.0, 225.0);
	  book(histos["lj1pt"], "lj1pt",        40, 25.0, 225.0);
	  book(histos["lj2pt"], "lj2pt",        40, 25.0, 225.0);
	  book(histos["wj1pt"],  "wj1pt",       40, 0.0, 250.0);
	  book(histos["topmass_pos"],   "topmass_pos",   40,100,300);
	  book(histos["topmass_neg"],   "topmass_neg",   40,100,300);
	  book(histos["topmass_recpos"],   "topmass_recpos",   40,100,300);
	  book(histos["topmass_recneg"],   "topmass_recneg",   40,100,300);
	  book(histos["charge_bpos"],   "charge_bpos",   40,-1,1);
	  book(histos["charge_bneg"],   "charge_bneg",   40,-1,1);
    }

    
    //run the selection and fill the histograms
    void analyze(const Event& event) {
		
	  int cut_num = 1;
      _cut_flow->fill(cut_num++);  
       histos["total"]->fill(1.0);	 
	   
	  //lepton
      const vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      if(leptons.size()!=1) vetoEvent;
	  _cut_flow->fill(cut_num++);
      const DressedLepton *lepton( &(leptons[0]) );
	  
	  // Truth level bquarks (select all):
	  Particles bquarks;
	  for (ConstGenParticlePtr p : HepMCUtils::particles(event.genEvent())) {
		  if (abs(p->pdg_id()) == PID::BQUARK) bquarks += Particle(p);
	  }
	  
	  // jets
      const Jets all_jets = applyProjection<FastJets>(event, "Jets").jetsByPt(allj_cut);
      std::vector<const Jet *> bjets;
      std::vector<const Jet *> nonbjets;
      std::vector<float> q_bjets;
      std::vector<float> qreco_bjets;
      for (const Jet& jet : all_jets) {

        // check for jet-lepton overlap -> do not consider for selection
        if (deltaR(jet, *lepton) < 0.4) continue;
        
		// truth charge
		float q = 0, maxdR = 0.4;	
		for (const Particle& b : bquarks) { // parton level association
			if (deltaR(jet, b) <maxdR) {q = b.charge()/b.abscharge(); maxdR = deltaR(jet, b);}
		}
		// reconstructed charge
		float k_power = 1; // weight as pT^k
		float qrec=getJetCharge(&jet, k_power);
        
        if(jet.bTagged(bjet_cut)) { bjets.push_back(&jet);    q_bjets.push_back(q);  qreco_bjets.push_back(qrec); }
        else                      { nonbjets.push_back(&jet);}
		
      }
		
      //b-jet cut
      if(bjets.size()==0) vetoEvent;
      _cut_flow->fill(cut_num++);
	  
       // skip events with bad association
       bool passchargeid( fabs(q_bjets[0]) );
       if( !passchargeid ) vetoEvent;
       _cut_flow->fill(cut_num++);
	  
       // Select events with bjet_pt>50GeV (otherwise, we don't see mt peak with POWHEG samples)
       if( bjets[0]->pt()<50*GeV ) vetoEvent;	  
       _cut_flow->fill(cut_num++);

       // Select events with at least 1 light jet
       //if(nonbjets.size()==0) vetoEvent;
       //_cut_flow->fill(cut_num++);
	  
      //bool passbj2( bjets.size()==1 || bjets[1]->pt()<50*GeV );
      //bool passj2( nonbjets.size()<2 || nonbjets[1]->pt()<50*GeV);

      //missing momentum
      const Vector3& met_p3 = apply<MissingMomentum>(event, "MissingET").vectorMissingPt();

      //reconstruct the W boson
      const FourMomentum lepton_p4 = lepton->mom();
      double pz = findZcomponent(lepton_p4, met_p3);
      FourMomentum neutrino_p4(sqrt(sqr(met_p3.x()) + sqr(met_p3.y()) + sqr(pz)), met_p3.x(), met_p3.y(), pz);
      FourMomentum w_p4(lepton_p4+neutrino_p4);
      bool passw(fabs(w_p4.eta())<4.0 && w_p4.pt()<120*GeV);
      bool passwj1(true);
      float wj1pt(-9999.);
      if(nonbjets.size()>0) {
        FourMomentum wj1(w_p4+nonbjets[0]->mom());
        wj1pt=wj1.pt();
        passwj1 &= (wj1pt<140*GeV);
      }

      //top quark is used to define off-shellness
      std::vector<float> tmass;
      for(size_t i=0; i<bjets.size();i++) {
        FourMomentum top(w_p4+bjets[i]->mom());
        tmass.push_back(top.mass());
      }

      //fill histograms  
	  if(q_bjets[0]>0) {
		histos["topmass_pos"]->fill(tmass[0]/GeV);
		histos["bp1eta"]->fill(bjets[0]->eta());
		histos["charge_bpos"]->fill(qreco_bjets[0]);
	  }
	  else {
		histos["topmass_neg"]->fill(tmass[0]/GeV);
		histos["bm1eta"]->fill(bjets[0]->eta());
		histos["charge_bneg"]->fill(qreco_bjets[0]);
	  }	  
	  if(qreco_bjets[0]>0) {
		histos["topmass_recpos"]->fill(tmass[0]/GeV);
	  }
	  else {
		histos["topmass_recneg"]->fill(tmass[0]/GeV);
	  }

	  histos["njets"]->fill(all_jets.size());
	  histos["nbjets"]->fill(bjets.size());
	  histos["nljets"]->fill(nonbjets.size());
	  histos["wpt"]->fill(w_p4.pt()/GeV);
	  histos["weta"]->fill(w_p4.eta());
	  histos["b1pt"]->fill(bjets[0]->pt()/GeV);
	  if(bjets.size()>1) histos["b2pt"]->fill(bjets[1]->pt()/GeV);
	  if(nonbjets.size()>0) {
		histos["lj1pt"]->fill(nonbjets[0]->pt()/GeV);
		histos["wj1pt"]->fill(wj1pt/GeV);
	  }
	  if(nonbjets.size()>1) histos["lj2pt"]->fill(nonbjets[1]->pt()/GeV);

    } // end analyze

    //solve the quadratic equation to find the pz of the neutrino based on the W mass constraint and the lepton+MET system
    //this was just copied from https://rivet.hepforge.org/analyses/MC_TTBAR.html
    double findZcomponent(const FourMomentum& lepton, const Vector3& met) const {
      // estimate z-component of momentum given lepton 4-vector and MET 3-vector
      double pz_estimate;
      double m_W = 80.399*GeV;
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.x() + lepton.py() * met.y());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.perp() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pz_estimate = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pz_estimate = quad[0];
        else                          pz_estimate = quad[1];
      }
      return pz_estimate;
    }

    //jet charge based on the DELPHI paper
    float getJetCharge(const Jet *jet,float k=0.5) {

      float q(0.), total(0.);
      for (const Particle& p : jet->particles(Cuts::pT > 1.*GeV)) {
        if (p.charge() == 0) continue;
        q += p.charge() * pow( p.pt()/GeV, k );
		total += pow( p.pt()/GeV, k );
      }

      return (total ? q / total : total);
    }


    //scale to cross section
    void finalize() {
      std::cout << "Cross section = " << crossSection()/picobarn << std::endl;
      std::cout << " sumOfWeights() = " << sumOfWeights() << std::endl;
      //double sf (1.0);
      //if(sumOfWeights()!=0) sf = (crossSection()/picobarn) / sumOfWeights();
      //for (auto hist : histos) { scale(hist.second, sf); }
	  //asymm(histos["topmass_neg"],histos["topmass_pos"], _topmass_asymm);
	  scale(_cut_flow,1/crossSection());

    }

    //@}


  private:

    // @name Histogram data members
    //@{
    Cut fs_cut,lepton_cut,allj_cut,bjet_cut;

    map<string, Histo1DPtr> histos;
	Histo1DPtr _cut_flow;
	//Scatter2DPtr _topmass_asymm;

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_WbWb_internal);

}

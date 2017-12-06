#include "Muon.h"
#include "Electron.h"
#include "Jet.h"
#include "EventID.h"
#include "GeneratorInfo.h"

#include "tqgamma_fNtuple.hh"
#include <string>
#include <vector>
#include <iostream>

#include "Muon.cc"
#include "Electron.cc"
#include "Jet.cc"
#include "GeneratorInfo.cc"
#include "EventID.cc"
#include "Lepton.cc"
#include "Candidate.cc"
#include "CandidateWithID.cc"

using namespace std;

/*
  1. Pt фотона
  2. b-tag дискриминатор? o_O
  3. Pt b-jet
  4. Pt muon (tcy)
  5. cos(p_t, p_y)
  6. dR(b-jet, y)
  7. dR(muon, y)
  8. charge of lepton (tuy)
  9. jet multiplicity

*/

float get_dR(float eta_a, float phi_a, float eta_b, float phi_b){
  return sqrt( pow(eta_a - eta_b, 2) + pow(phi_a - phi_b, 2) );
}

int analyse_single_file(string inp_name, string out_name, Long64_t max_events){
  // msg("analyse_single_file() ... ", inp_name, out_name, max_events);

  // SETUP OUTPUT TTREE
  TFile * file_out = new TFile( out_name.c_str(), "RECREATE");
  TTree * tree_out = new TTree("ttree", "ttree");
  float G_PT, BJ_BTAG, BJ_PT, MU_PT, COS_T_G, dR_BJ_G, dR_Mu_G, MU_CHARGE, J_MULT;
  vector<pec::Muon> muons; 
  auto muonsPtr = &muons;
  vector<pec::Electron> electrons;
  auto elePtr = &electrons;
  vector<pec::Jet> jets;
  auto jetPtr = &jets;
  vector<pec::Candidate> storeMETs;
  auto metPtr = &storeMETs;
  pec::EventID id;
  auto eventIDPtr = &id;
  pec::GeneratorInfo genInfo;
  auto genInfoPtr = &genInfo;
  tree_out->Branch("pecMuons/muons",         &muonsPtr);
  tree_out->Branch("pecElectrons/electrons",     &elePtr);
  tree_out->Branch("pecJetMET/jets",          &jetPtr);
  tree_out->Branch("METs",          &metPtr);
  tree_out->Branch("eventId",          &eventIDPtr);
  tree_out->Branch("generator",          &genInfoPtr);

  tree_out->Branch("var@G_PT",      &G_PT);
  tree_out->Branch("var@BJ_BTAG",   &BJ_BTAG);
  tree_out->Branch("var@BJ_PT",     &BJ_PT);
  tree_out->Branch("var@MU_PT",     &MU_PT);
  tree_out->Branch("var@COS_T_G",   &COS_T_G);
  tree_out->Branch("var@dR_BJ_G",   &dR_BJ_G);
  tree_out->Branch("var@dR_Mu_G",   &dR_Mu_G);
  tree_out->Branch("var@MU_CHARGE", &MU_CHARGE);
  tree_out->Branch("var@J_MULT",    &J_MULT);

  // INPUT TTREE
  TFile * file = TFile::Open( inp_name.c_str() );
  tqgamma_Event event;
  event.Init(file);
  Long64_t nevents = 0;

  pec::Muon muon;
  pec::Electron ele;
  pec::Jet jet;

  if (max_events == -1)
    max_events = event.fChain_Event->GetEntries();

  for(Long64_t nevent = 0; nevent < max_events; nevent++) {
    if (nevent % 1000 == 0)
        cout << nevent << endl;
    
    event.GetEntry( nevent );


////////////////////////////////////////////////////////////////////////
// M U O N S
////////////////////////////////////////////////////////////////////////

    muons.clear();
    for (unsigned i = 0; i < event.MuonTight_size; ++i)
    {
      //pat::Muon const &mu = srcMuons->at(i);
      muon.Reset();
      
      // Set four-momentum. Mass is ignored
      muon.SetPt(event.MuonTight_PT[i]);
      muon.SetEta(event.MuonTight_Eta[i]);
      muon.SetPhi(event.MuonTight_Phi[i]);

      muon.SetCharge(event.MuonTight_Charge[i]);
      
      // Relative isolation with delta-beta correction [1]
      //[1] https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=22#Muon_Isolation
      //auto const &isoR04 = mu.pfIsolationR04();
      //muon.SetRelIso((isoR04.sumChargedHadronPt +
      // max(isoR04.sumNeutralHadronEt + isoR04.sumPhotonEt - 0.5 * isoR04.sumPUPt, 0.)) / mu.pt());
      muon.SetRelIso(event.MuonLoose_IsolationVar[i]);  
        
      // Moun identification bits [1]. Note this does not imply selection on isolation or
      //kinematics
      //[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2?rev=22#Muon_Identification
      muon.SetBit(0, 1);
      muon.SetBit(1, 1);
      muon.SetBit(2, 1);
      
      
      // Evaluate user-defined selectors if any
      unsigned const nUsedBits = 3;
      
      //for (unsigned i = 0; i < muSelectors.size(); ++i)
      //  muon.SetBit(nUsedBits + i, muSelectors[i](mu));
      
      // The muon is set up. Add it to the vector
      muons.emplace_back(muon);
    }


////////////////////////////////////////////////////////////////////////
// E L E C T R O N S
////////////////////////////////////////////////////////////////////////

    electrons.clear();
    for (unsigned i = 0; i < event.ElectronLoose_size; ++i)
    {
        ele.Reset();
        
        
        // Set four-momentum. Mass is ignored
        ele.SetPt(event.ElectronLoose_PT[i]);
        ele.SetEta(event.ElectronLoose_Eta[i]);
        ele.SetPhi(event.ElectronLoose_Phi[i]);
        
        ele.SetCharge(event.ElectronLoose_Charge[i]);
        
        
        // Isolation is calculated by a dedicated method
        //ele.SetRelIso(CalculateRhoIsolation(el, *rho));
        ele.SetRelIso(event.ElectronLoose_IsolationVar[i]);
        
        
        // Set pseudorapidity of the associated supercluster
        //ele.SetEtaSC(el.superCluster()->eta());
        
        bool isTightEle = false;
        for (unsigned j = 0; j < event.ElectronTight_size; ++j) {
            double dr = get_dR(event.ElectronLoose_Eta[i], event.ElectronLoose_Phi[i], event.ElectronTight_Eta[j], event.ElectronTight_Phi[j]);
            //cout << "DR: " << dr << endl;
            if (dr < 0.4) {
                isTightEle = true;
                break;
            }
        }
        ele.SetBooleanID(0, 1);
        //cout << isTightEle << endl;
        ele.SetBooleanID(1, isTightEle);
        
        /*
        // Copy non-triggering MVA ID stored as a userFloat
        // ele.SetContinuousID(0,
        //   el.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"));
        unsigned nUsedContIDs = 2;
        
        
        // Copy embedded ID decisions
        unsigned const nEmbeddedBoolIDs = embeddedBoolIDLabels.size();
        unsigned const nEmbeddedContIDs = embeddedContIDLabels.size();
        
        for (unsigned i = 0; i < nEmbeddedBoolIDs; ++i)
            ele.SetBooleanID(i, (el.electronID(embeddedBoolIDLabels.at(i)) > 0.5f));
            //^ Since pat::Electron::electronID returns a float, need to be accurate with the
            //conversion to a boolean value
        
        for (unsigned i = 0; i < nEmbeddedContIDs; ++i)
            ele.SetContinuousID(nUsedContIDs + i,
              el.electronID(embeddedContIDLabels.at(i)));
        
        nUsedContIDs += nEmbeddedContIDs;
        
        
        // Copy additional ID decisions from the maps
        Ptr<pat::Electron> const elPtr(srcElectrons, i);
        
        for (unsigned i = 0; i < boolIDMaps.size(); ++i)
            ele.SetBooleanID(nEmbeddedBoolIDs + i, (*boolIDMaps.at(i))[elPtr]);
        
        for (unsigned i = 0; i < contIDMaps.size(); ++i)
            ele.SetContinuousID(nUsedContIDs + i, (*contIDMaps.at(i))[elPtr]);
        
        
        // Evaluate loose selection on impact parameters [1]. It is implemented as in [2-3].
        //[1] https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Offline_selection_criteria
        //[2] https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDxyCut.cc#L58-L68
        //[3] https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDzCut.cc#L58-L68
        bool passIPCuts;
        double const d0 = std::abs(el.gsfTrack()->dxy(firstPV.position()));
        double const dz = std::abs(el.gsfTrack()->dz(firstPV.position()));
        
        if (std::abs(el.superCluster()->eta()) < 1.479)
            passIPCuts = (d0 < 0.05 and dz < 0.10);  // units are cm
        else
            passIPCuts = (d0 < 0.10 and dz < 0.20);
        
        ele.SetBit(0, passIPCuts);
        
        
        // Evaluate user-defined selectors if any
        unsigned const nUsedBits = 1;
        
        for (unsigned i = 0; i < eleSelectors.size(); ++i)
            ele.SetBit(nUsedBits + i, eleSelectors[i](el));
        */
        
        // The electron is set up. Add it to the vector.
        electrons.emplace_back(ele);
    }


////////////////////////////////////////////////////////////////////////
// J E T S
////////////////////////////////////////////////////////////////////////


    jets.clear();
    
    for (unsigned int i = 0; i < event.JetPUPPI_size; ++i)
    {
        jet.Reset();
        
        //reco::Candidate::LorentzVector const &rawP4 = j.correctedP4("Uncorrected");
        
        jet.SetPt(event.JetPUPPI_PT[i]);
        jet.SetEta(event.JetPUPPI_Eta[i]);
        jet.SetPhi(event.JetPUPPI_Phi[i]);
        jet.SetM(event.JetPUPPI_Mass[i]);
        
        /*
        if (not rawJetMomentaOnly)
        {
            if (runOnData)
            {
                storeJet.SetCorrFactor(1. / j.jecFactor("Uncorrected"));
                //^ Here jecFactor("Uncorrected") returns the factor to get raw momentum starting
                //from the corrected one. Since in fact the raw momentum is stored, the factor is
                //inverted
            }
            else
            {
                double const jerFactorNominal = j.userFloat("jerFactorNominal");
                double const jerFactorUp = j.userFloat("jerFactorUp");
                double const jerFactorDown = j.userFloat("jerFactorDown");
                
                storeJet.SetCorrFactor(1. / j.jecFactor("Uncorrected") * jerFactorNominal);
                //^ See the comment for real data concerning the inverted JEC factor
                storeJet.SetJECUncertainty(j.userFloat("jecUncertainty"));
                
                // For JER the variation is not necessarily symmetric. Save the largest
                //variation. Information about the sign of the variation is preserved, and the
                //stored uncertainty is negative if "up" variation of JER decreases the smearing
                //factor
                if (std::abs(jerFactorUp - jerFactorNominal) >
                  std::abs(jerFactorDown - jerFactorNominal))
                    storeJet.SetJERUncertainty(jerFactorUp / jerFactorNominal - 1.);
                else
                    storeJet.SetJERUncertainty(1. - jerFactorDown / jerFactorNominal);
            }
        }
        */
        
        //storeJet.SetArea(j.jetArea());
        //jet.SetCharge(event.JetPUPPI_Charge[i]);
        
        
        // Save b-tagging discriminators
        jet.SetBTag(pec::Jet::BTagAlgo::CSV , (event.JetPUPPI_DeepCSV[i] & 7) > 0);
        jet.SetBTag(pec::Jet::BTagAlgo::CMVA, (event.JetPUPPI_MVAv2[i] & 7) > 0);
        //jet.SetBTagDNN(event.JetPUPPI_DeepCSV[i]);
        /*
        j.bDiscriminator("pfDeepCSVJetTags:probbb"),
          j.bDiscriminator("pfDeepCSVJetTags:probb"),
          j.bDiscriminator("pfDeepCSVJetTags:probcc"),
          j.bDiscriminator("pfDeepCSVJetTags:probc"),
          j.bDiscriminator("pfDeepCSVJetTags:probudsg"));
        */
        
        // Save pileup ID
        //[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID?rev=29#Information_for_13_TeV_data_anal
        //storeJet.SetPileUpID(j.userFloat("pileupJetId:fullDiscriminant"));
        
        /*
        if (not runOnData)
        {
            storeJet.SetFlavour(j.hadronFlavour(), j.partonFlavour(),
              (j.genParton() ? j.genParton()->pdgId() : 0));
            storeJet.SetBit(0, bool(j.userInt("hasGenMatch")));
        }
        
        
        // Loose PF jet ID [1]. Accessors to energy fractions take into account JEC, so there is no
        //need to undo the corrections.
        //[1] https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016?rev=1
        bool passPFID = false;
        double const absEta = std::abs(rawP4.Eta());
        
        if (absEta <= 2.7)
        {
            bool const commonCriteria = (j.neutralHadronEnergyFraction() < 0.99 and
              j.neutralEmEnergyFraction() < 0.99 and
              (j.chargedMultiplicity() + j.neutralMultiplicity() > 1));
            
            if (absEta <= 2.4)
                passPFID = (commonCriteria and j.chargedHadronEnergyFraction() > 0. and
                  j.chargedMultiplicity() > 0 and j.chargedEmEnergyFraction() < 0.99);
            else
                passPFID = commonCriteria;
        }
        else if (absEta <= 3.)
            passPFID = (j.neutralMultiplicity() > 2 and
              j.neutralHadronEnergyFraction() < 0.98 and j.neutralEmEnergyFraction() > 0.01);
        else
            passPFID = (j.neutralMultiplicity() > 10 and j.neutralEmEnergyFraction() < 0.9);
        
        storeJet.SetBit(1, passPFID);
        
        
        // User-defined selectors if any. The first two bits have already been used.
        for (unsigned i = 0; i < jetSelectors.size(); ++i)
            storeJet.SetBit(i + 2, jetSelectors[i](j));
        
        
        // The jet is set up. Add it to the vector
        storeJets.emplace_back(storeJet);
        
        
        // Update the partial T1 MET correction
        auto const deltaT1JetP4 = -(j.p4() - j.correctedP4("L1FastJet"));
        metT1Corr += TVector2(deltaT1JetP4.Px(), deltaT1JetP4.Py());
        */
    }


    
    storeMETs.clear();
    pec::Candidate storeMET;
    //^ Will reuse this object to fill the vector of METs
    
    // Nominal MET (type-I corrected)
    storeMET.Reset();
    storeMET.SetPt(event.PuppiMissingET_MET[0]);
    storeMET.SetPhi(event.PuppiMissingET_Phi[0]);
    storeMETs.emplace_back(storeMET);


////////////////////////////////////////////////////////////////////////

    id.SetRunNumber(event.Run);
    id.SetLumiSectionNumber(event.Lumi);
    id.SetEventNumber(event.Event);


////////////////////////////////////////////////////////////////////////
// O T H E R
//////////////////////////////////////////////////////////////////////// 


    // msg_progress( float(nevent) / max_events );

    //   float G_PT, BJ_BTAG, BJ_PT, MU_PT, COS_T_G, dR_BJ_G, dR_Mu_G, MU_CHARGE, J_MULT;
    //   1. Pt фотона

    // 2. b-tag дискриминатор? o_O
    if(event.JetPUPPI_size){
      BJ_BTAG = event.JetPUPPI_MVAv2[0];
    } else BJ_BTAG = -1;

    // 3. Pt b-jet
    if(event.JetPUPPI_size){
      BJ_PT = event.JetPUPPI_PT[0];
    } else BJ_PT = -1;

    //  4. Pt muon (tcy)
    if(event.MuonTight_size){
      MU_PT = event.MuonTight_PT[0];
    } else MU_PT = -1;

    //  5. cos(p_t, p_y)
    COS_T_G = -1;

    //   6. dR(b-jet, y)

    //   8. charge of lepton (tuy)
    if(event.MuonTight_size){
      MU_CHARGE = event.MuonTight_Charge[0];
    } else MU_CHARGE = -10;

    //   9. jet multiplicity
    J_MULT = event.JetPUPPI_size;

    tree_out->Fill();
  }

  file_out->cd();
  tree_out->Write();
  // msg("analyse_single_file() ... ", inp_name, out_name, max_events, "done!");

  return 0;
}


int s1_create_var_trees(){

  // create part
  analyse_single_file("MiniEvents.root", "vars.root", -1);

  return 0;
}


int main() {
    s1_create_var_trees();
    return 0;
}




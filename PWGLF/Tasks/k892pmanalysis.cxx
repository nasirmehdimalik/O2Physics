// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file k892analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// modifeid by  Nasir Mehdi Malik  for reconstruction of K(892)*pm
#include <CCDB/BasicCCDBManager.h>
#include <TLorentzVector.h>

#include <cmath>
#include <array>
#include <cstdlib>
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct k892pmanalysis {
  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  // Configurables

  Configurable<bool> isData{"isData", 1, "set 0 for mc"};

  // event
  Configurable<double> cfgCutVertex{"cfgCutvertex", 10.0f, "vertex z"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for daughter tracks"};
  Configurable<int> cfgMinCrossedRows{"cfgMinCrossedRows", 70, "min crossed rows"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "pt cut on daughters track"};

  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for track"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for track"};

  // vo
  Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.8f, "rapidty"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA between  V0 Daughters pi plus and pi minus"};
  Configurable<float> dcanegtopv{"dcanegtopv", .3, "DCA Neg pi minus To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .3, "DCA Pos pi plus To PV"};
  Configurable<float> dcav0topv{"dcav0topv", .3, "DCA v0 To PV"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "V0 Cosine of Pointing Angle"};
  Configurable<float> v0radius{"v0radius", 0.5, "V0 radius"};
  Configurable<float> lifetimecut{"lifetimecut", 20., "lifetimecut"};

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Mass QA (quick check)

    AxisSpec multAxis = {3000, 0, 3000, "mult"};
    AxisSpec ptAxiss = {100, 0.0f, 10.0f, "#it{p}_{T}  (GeV/#it{c})"};
    AxisSpec invmassAxis = {900, 0.6, 1.5, "{M}_{#pion #K^{0}} (GeV/#it{c}^2)"};
    // 3d histogram
    if (isData) {

      histos.add("h3k892pminvmass", "Invariant mass of K(892)pm", kTHnSparseD, {multAxis, ptAxiss, invmassAxis});
      histos.add("h3k892pminvmassME", "Invariant mass of K(892)pm mixed event", kTHnSparseD, {multAxis, ptAxiss, invmassAxis});

      histos.add("collision/hCentrality", "Centrality distribution", kTH1F, {{3000, 0.0, 3000.0}});
      histos.add("collision/hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("collision/hNcontributor", "Number of primary vertex contributor", kTH1F, {{1000, 0.0f, 1000.0f}});
    }
    // bachlor track
    histos.add("Bachtrack/hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("Bachtrack/hDCAxy", "DCAxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("Bachtrack/hDCAz", "DCAz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("Bachtrack/hNsigmapionTPC_withoutNsigmacut", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("Bachtrack/hNsigmapionTOF_withoutNSigmacut", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("Bachtrack/hNsigmapionTPC", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("Bachtrack/hNsigmapionTOF", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("Bachtrack/hNsigmapion_TPC_TOF_MAP", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add<TH1>("Bachtrack/tpcNClsCR", "tpcnclscrossedrows; #tpcnclscrossedrows; counts", HistType::kTH1F, {{400, 50.0, 150.0}});

    // vo track

    histos.add("v0track/hRapidity", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("v0track/hV0Radius", "v0 radius", kTH1F, {{500, 0.0f, 100.0f, "cm"}});
    histos.add("v0track/hV0CosPA", "V0 Cosine of Pointing Angle", kTH1F, {{200, 0.9f, 1.2f}});
    histos.add("v0track/hDCAPosToPV", "dca pos to pv", kTH1F, {{100, -1.0f, 1.0f, "cm"}});
    histos.add("v0track/hDCANegToPV", "dca neg to pv", kTH1F, {{100, -1.0f, 1.0f, "cm"}});
    histos.add("v0track/hDCAV0Dau", "dca v0 dau", kTH1F, {{100, 0.0f, 4.0f, ""}});
    histos.add("v0track/hDCAV0toPV", "dca v0 t0 pV", kTH1F, {{250, 0.0f, 0.5f, "cm"}});

    if (!isData) {
      histos.add("h3k892pminvmass_gen", "gen Invariant mass of K(892)pm", kTHnSparseD, {ptAxiss, invmassAxis});
      histos.add("h3k892pminvmass_rec", "rec Invariant mass of K(892)pm", kTHnSparseD, {ptAxiss, invmassAxis});
    }
  }

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massK0 = RecoDecay::getMassPDG(kK0);
  double mass{0.};
  double recMass{0.};
  double rapidity;
  double pT{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  // static constexpr float defaultLifetimeCuts[1] = {20.};

  template <typename T>
  bool selectpion(const T& candidate, const char* pname = "sE")
  {
    histos.fill(HIST("Bachtrack/hNsigmapionTPC_withoutNsigmacut"), candidate.pt(), candidate.tpcNSigmaPi());
    histos.fill(HIST("Bachtrack/hNsigmapionTOF_withoutNSigmacut"), candidate.pt(), candidate.tofNSigmaPi());
    if (candidate.tpcNClsCrossedRows() > cfgMinCrossedRows &&
        TMath::Abs(candidate.dcaXY()) < cfgCutDCAxy &&
        TMath::Abs(candidate.dcaZ()) < cfgCutDCAz &&
        TMath::Abs(candidate.eta()) < cfgCutEta &&
        TMath::Abs(candidate.pt()) < cfgCutPT) {

      if ((candidate.pt() < 0.3) && (std::abs(candidate.tpcNSigmaPi()) > 6.0))
        return false;

      if ((candidate.pt() >= 0.3) && (candidate.pt() < 0.4) && (std::abs(candidate.tpcNSigmaPi()) > 4.0))
        return false;

      if ((candidate.pt() >= 0.4) && (std::abs(candidate.tpcNSigmaPi()) > 3.0))
        return false;
      if (strcmp(pname, "sE") == 0) {
        histos.fill(HIST("Bachtrack/hNsigmapionTPC"), candidate.pt(), candidate.tpcNSigmaPi());
        histos.fill(HIST("Bachtrack/hNsigmapionTOF"), candidate.pt(), candidate.tofNSigmaPi());
        histos.fill(HIST("Bachtrack/hNsigmapion_TPC_TOF_MAP"), candidate.tpcNSigmaPi(), candidate.tofNSigmaPi());
        histos.fill(HIST("Bachtrack/hDCAxy"), candidate.dcaXY());
        histos.fill(HIST("Bachtrack/hDCAz"), candidate.dcaZ());
        histos.fill(HIST("Bachtrack/hEta"), candidate.eta());
        histos.fill(HIST("Bachtrack/tpcNClsCR"), candidate.tpcNClsCrossedRows());
      }

      return true;
    }
    return false;
  }
  template <typename c0, typename T0>
  bool selectV0(const c0& cc, const T0& candidate, const char* pname = "sE")
  {
    if (candidate.v0radius() > v0radius &&
        candidate.v0cosPA(cc.posX(), cc.posY(), cc.posZ()) > v0cospa &&
        candidate.dcav0topv(cc.posX(), cc.posY(), cc.posZ()) < dcav0topv &&
        candidate.distovertotmom(cc.posX(), cc.posY(), cc.posZ()) * massK0 < lifetimecut &&
        TMath::Abs(candidate.yK0Short()) < cfgCutRapidity &&
        TMath::Abs(candidate.pt()) > cfgCutPT &&
        TMath::Abs(candidate.dcapostopv()) > dcapostopv &&
        TMath::Abs(candidate.dcanegtopv()) > dcanegtopv &&
        TMath::Abs(candidate.dcaV0daughters()) < dcav0dau) {
      if (strcmp(pname, "sE") == 0) {
        histos.fill(HIST("v0track/hV0CosPA"), candidate.v0cosPA(cc.posX(), cc.posY(), cc.posZ()));
        histos.fill(HIST("v0track/hDCAV0toPV"), candidate.dcav0topv(cc.posX(), cc.posY(), cc.posZ()));
        histos.fill(HIST("v0track/hV0Radius"), candidate.v0radius());
        histos.fill(HIST("v0track/hDCAPosToPV"), candidate.dcapostopv());
        histos.fill(HIST("v0track/hDCANegToPV"), candidate.dcanegtopv());
        histos.fill(HIST("v0track/hRapidity"), candidate.yK0Short());
        histos.fill(HIST("v0track/hDCAV0Dau"), candidate.dcaV0daughters());
      }
      return true;
    }
    return false;
  }

  template <typename TracksType, typename TracksV0>
  void fillHistograms(float multiplicity, const TracksType& trk1, const TracksV0& trk2, const char* pname)
  {

    pvec0 = array{trk1.px(), trk1.py(), trk1.pz()};
    pvec1 = array{trk2.px(), trk2.py(), trk2.pz()};
    auto arrMom = array{pvec0, pvec1};

    mass = RecoDecay::m(arrMom, array{massPi, massK0});
    pT = RecoDecay::pt(array{trk1.px() + trk2.px(), trk1.py() + trk2.py()});
    rapidity = RecoDecay::y(array{trk1.px() + trk2.px(), trk1.py() + trk2.py(), trk1.pz() + trk2.pz()}, mass);

    if (std::abs(rapidity) < 0.5) {

      if (strcmp(pname, "sE") == 0)
        histos.fill(HIST("h3k892pminvmass"), multiplicity, pT, mass);

      if (strcmp(pname, "mE") == 0)
        histos.fill(HIST("h3k892pminvmassME"), multiplicity, pT, mass);
    }
  }

  using V0TrackCandidate = aod::V0Datas;
  using TrackCandidate = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

  void processData(EventCandidates::iterator const& collision, TrackCandidate const& tracks, V0TrackCandidate const& V0tracks)
  {

    //    if(!collision.sel8()) return;
    if (TMath::Abs(collision.posZ()) > cfgCutVertex)
      return;
    float multiplicity = collision.multTPC();

    histos.fill(HIST("collision/hCentrality"), multiplicity);
    histos.fill(HIST("collision/hNcontributor"), collision.numContrib());
    histos.fill(HIST("collision/hVtxZ"), collision.posZ());

    for (auto track1 : tracks) {
      if (!selectpion(track1))
        continue;

      for (auto track2 : V0tracks) {
        if (!selectV0(collision, track2))
          continue;

        fillHistograms(multiplicity, track1, track2, "sE");
      }
    }
  }
  PROCESS_SWITCH(k892pmanalysis, processData, "Process Event for data", true);

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 1500}, "multiplicity axis for histograms"};

  ConfigurableAxis axisVertex{"axisVertex", {VARIABLE_WIDTH, -10, 10}, "vertex axis for bin"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default
  Pair<EventCandidates, TrackCandidate, V0TrackCandidate, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1};

  void processMixedEvent(EventCandidates const& collisions, TrackCandidate const& tracks, V0TrackCandidate const& V0tracks)
  {

    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      // if(!c1.sel8()) return;
      // if(!c2.sel8()) return;

      float multiplicity = c1.multTPC();
      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!selectpion(t1, "mE"))
          continue;
        if (!selectV0(c2, t2, "mE"))
          continue;

        fillHistograms(multiplicity, t1, t2, "mE");
      }
    }
  }

  PROCESS_SWITCH(k892pmanalysis, processMixedEvent, "Process EventMixing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892pmanalysis>(cfgc, TaskName{"lf-k892pmanalysis"})};
}

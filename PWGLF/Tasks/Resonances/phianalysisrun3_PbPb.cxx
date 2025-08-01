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
// Preliminary QA analysis task for resonances
// (1) For Run3
// (2) Event and track selection need to be optimized
// (3) particle = 0 --> phi
// (4) particle = 1 --> kstar
// (5) particle = 2 --> lambdastar
// (6) 4 process function (a) Data same event (b) Data mixed event (c) MC generated (d) MC reconstructed

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;

struct phianalysisrun3_PbPb {
  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  RCTFlagsChecker rctChecker;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmacutTPC{"nsigmacutTPC", 2.0f, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmacutTOF{"nsigmacutTOF", 2.0f, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> fillOccupancy{"fillOccupancy", true, "fill Occupancy"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> additionalEvSel1{"additionalEvSel1", true, "Additional evsel1"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", true, "Additional evsel4"};
  Configurable<bool> additionalEvSel5{"additionalEvSel5", true, "Additional evsel5"};
  Configurable<bool> additionalEvSel6{"additionalEvSel6", true, "Additional evsel6"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<int> nBkgRotations{"nBkgRotations", 3, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0f * TMath::Pi() / 6.0f, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0f * TMath::Pi() / 6.0f, "Maximum of rotation"};
  Configurable<bool> pdgcheck{"pdgcheck", true, "pdgcheck"};
  Configurable<bool> reco{"reco", true, "reco"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0, 3.5, 5.67, 7.45, 8.85, 10.0, 11.21, 12.26, 13.28, 14.23, 15.27}, "Binning of the impact parameter axis"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 3000, "Occupancy cut"};
  Configurable<int> centestimator{"centestimator", 0, "Select multiplicity estimator: 0 - FT0C, 1 - FT0A, 2 - FT0M, 3 - FV0A, 4 - PVTracks"};

  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy distribution", kTH1F, {{500, 0, 50000}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
    if (!isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRot", "Invariant mass of Phi meson Rotation", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassSame", "Invariant mass of Phi meson same", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h2PhiRapidity", "phi meson Rapidity", kTH2F, {{200, 0.0f, 20.0f}, {200, -4, 4}});
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("EL1", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("EL2", "MC Event statistics", kTH1F, {centAxis});
      histos.add("ES1", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("ES3", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("ES2", "MC Event statistics", kTH1F, {centAxis});
      histos.add("ES4", "MC Event statistics", kTH1F, {centAxis});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiGen1", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("Centrec", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("Centgen", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("hVtxZgen", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hVtxZrec", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h3PhiRec3", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3Phi1Rec3", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3PhiGen3", "Phi meson Gen", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixedMC", "Invariant mass of Phi meson Mixed", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassSameMC", "Invariant mass of Phi meson same", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRotMC", "Invariant mass of Phi meson Rotation", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h2PhiGen1", "Phi meson gen", kTH2F, {ptAxis, impactParAxis});
      histos.add("h1PhiRec1", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1Phimassgen", "Phi meson gen", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimasssame", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassmix", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrot", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phi1massrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phipt", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hOccupancy1", "Occupancy distribution", kTH1F, {{500, 0, 50000}});
      histos.add("h1PhifinalRec", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1Phifinalgenmass", "Phi meson gen mass", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h3PhifinalRec", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h1PhifinalGen", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h2PhifinalGen", "Phi meson Gen", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("hMC1", "MC Event statistics", kTH1F, {{15, 0.0f, 15.0f}});
      histos.add("Centrec1", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("Centgen1", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("h1PhiRecsplit1", "Phi meson Rec split", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterRec", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterGenCen", "Impact parameter of generated MC events", kTH2F, {impactParAxis, centAxis});
      histos.add("hImpactParameterRecCen", "Impact parameter of generated MC events", kTH2F, {impactParAxis, centAxis});
      histos.add("TOF_Nsigma_MC", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      histos.add("TPC_Nsigma_MC", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      histos.add("TOF_Nsigma1_MC", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      histos.add("TPC_Nsigma1_MC", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      histos.add("trkDCAxy", "DCAxy distribution of positive kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      histos.add("trkDCAz", "DCAxy distribution of negative kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
      if (doprocessEvtLossSigLossMC) {
        histos.add("QAevent/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
        histos.add("QAevent/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
        histos.add("QAevent/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{120, 0.0f, 120.0f}, impactParAxis});
        histos.add("QAevent/phigenBeforeEvtSel", "phi before event selections", kTH2F, {ptAxis, impactParAxis});
        histos.add("QAevent/phigenAfterEvtSel", "phi after event selections", kTH2F, {ptAxis, impactParAxis});
      }
    }

    histos.add("hEta", "eta of kaon track candidates", HistType::kTH2F, {{200, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("hPhi", "phi of kaon track candidates", HistType::kTH2F, {{65, 0, 6.5}, {200, 0.0f, 20.0f}});

    // DCA QA
    // DCA histograms: separate for positive and negative kaons, range [-1.0, 1.0]
    histos.add("QAbefore/trkDCAxy_pos", "DCAxy distribution of positive kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAbefore/trkDCAxy_neg", "DCAxy distribution of negative kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAbefore/trkDCAz_pos", "DCAz distribution of positive kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAbefore/trkDCAz_neg", "DCAz distribution of negative kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});

    histos.add("QAbefore/trkDCAxypt_pos", "DCAxy distribution of positive kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAbefore/trkDCAxypt_neg", "DCAxy distribution of negative kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAbefore/trkDCAzpt_pos", "DCAz distribution of positive kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAbefore/trkDCAzpt_neg", "DCAz distribution of negative kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});

    histos.add("QAafter/trkDCAxy_pos", "DCAxy distribution of positive kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAafter/trkDCAxy_neg", "DCAxy distribution of negative kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAafter/trkDCAz_pos", "DCAz distribution of positive kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});
    histos.add("QAafter/trkDCAz_neg", "DCAz distribution of negative kaon track candidates", HistType::kTH1F, {{150, -1.0f, 1.0f}});

    histos.add("QAafter/trkDCAxypt_pos", "DCAxy distribution of positive kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAafter/trkDCAxypt_neg", "DCAxy distribution of negative kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAafter/trkDCAzpt_pos", "DCAz distribution of positive kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    histos.add("QAafter/trkDCAzpt_neg", "DCAz distribution of negative kaon track candidates", HistType::kTH2F, {{150, -1.0f, 1.0f}, {200, 0.0f, 20.0f}});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Mapka_all_pos", "TOF + TPC Combined PID for positive Kaon;#sigma_{TOF}^{K^{+}};#sigma_{TPC}^{K^{+}}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAbefore/TOF_TPC_Mapka_all_neg", "TOF + TPC Combined PID for negative Kaon;#sigma_{TOF}^{K^{-}};#sigma_{TPC}^{K^{-}}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});

    histos.add("QAbefore/TOF_Nsigma_all_pos", "TOF NSigma for positive Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{K^{+}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAbefore/TOF_Nsigma_all_neg", "TOF NSigma for negative Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{K^{-}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});

    histos.add("QAbefore/TPC_Nsigma_all_pos", "TPC NSigma for positive Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{K^{+}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAbefore/TPC_Nsigma_all_neg", "TPC NSigma for negative Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{K^{-}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});

    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Mapka_all_pos", "TOF + TPC Combined PID for positive Kaon;#sigma_{TOF}^{K^{+}};#sigma_{TPC}^{K^{+}}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAafter/TOF_TPC_Mapka_all_neg", "TOF + TPC Combined PID for negative Kaon;#sigma_{TOF}^{K^{-}};#sigma_{TPC}^{K^{-}}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});

    histos.add("QAafter/TOF_Nsigma_all_pos", "TOF NSigma for positive Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{K^{+}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAafter/TOF_Nsigma_all_neg", "TOF NSigma for negative Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{K^{-}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});

    histos.add("QAafter/TPC_Nsigma_all_pos", "TPC NSigma for positive Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{K^{+}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAafter/TPC_Nsigma_all_neg", "TPC NSigma for negative Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{K^{-}}", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
  }

  double massKa = o2::constants::physics::MassKPlus;
  double rapidity;
  double genMass, recMass, resolution;
  ROOT::Math::PxPyPzMVector phiMother, daughter1, daughter2;
  double mass{0.};
  double massrotation{0.};
  double pT{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  array<float, 3> pvec1rotation;
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!isNoTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmacutTPC) {
      return true;
    }
    if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmacutTPC) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmacutTPC) {
      return true;
    }
    if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmacutTPC && std::abs(candidate.tofNSigmaKa()) < nsigmacutTOF) {
      return true;
    }
    return false;
  }
  template <typename CollType>
  bool myEventSelections(const CollType& collision)
  {
    if (std::abs(collision.posZ()) > cfgCutVertex)
      return false;

    if (!collision.sel8())
      return false;

    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;

    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;

    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;

    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    if (additionalEvSel5 && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return false;
    if (additionalEvSel6 && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      return false;
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy))
      return false;

    return true;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }
  template <typename T1, typename T2>
  void fillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, float massd1, float massd2)
  {
    pvec0 = std::array<float, 3>{candidate1.px(), candidate1.py(), candidate1.pz()};
    pvec1 = std::array<float, 3>{candidate2.px(), candidate2.py(), candidate2.pz()};
    auto arrMom = std::array<std::array<float, 3>, 2>{pvec0, pvec1};

    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    mass = RecoDecay::m(arrMom, std::array<float, 2>{massd1, massd2});

    pT = RecoDecay::pt(std::array<float, 2>{
      candidate1.px() + candidate2.px(),
      candidate1.py() + candidate2.py()});

    rapidity = RecoDecay::y(std::array<float, 3>{
                              candidate1.px() + candidate2.px(),
                              candidate1.py() + candidate2.py(),
                              candidate1.pz() + candidate2.pz()},
                            mass);

    constexpr float kRapidityCut = 0.5;
    constexpr int kOppositeCharge = 0;

    // default filling
    if (std::abs(rapidity) < kRapidityCut && track1Sign * track2Sign < kOppositeCharge) {

      if (unlike) {
        histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
        histos.fill(HIST("h2PhiRapidity"), pT, rapidity);
      }
      if (mix) {
        histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
      }
    }
  }
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcacutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  // using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::McCollisionLabels>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFV0As>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                    aod::McTrackLabels>>;
  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "multiplicity  for bin"};

  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  SliceCache cache;

  using BinningTypeVertexContributor1 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeVertexContributor2 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeVertexContributor3 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor4 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;
  ROOT::Math::PxPyPzMVector phiMesonMother, kaonPlus, kaonMinus;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 2.5);
    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 3.5);
    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 4.5);
    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 5.5);
    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 6.5);
    if (additionalEvSel5 && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 7.5);
    if (additionalEvSel6 && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 8.5);
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    histos.fill(HIST("hEvtSelInfo"), 9.5);
    float multiplicity{-1};
    const int kCentFT0C = 0;
    const int kCentFT0A = 1;
    const int kCentFT0M = 2;
    const int kCentFV0A = 3;

    if (centestimator == kCentFT0C) {
      multiplicity = collision.centFT0C();
    } else if (centestimator == kCentFT0A) {
      multiplicity = collision.centFT0A();
    } else if (centestimator == kCentFT0M) {
      multiplicity = collision.centFT0M();
    } else if (centestimator == kCentFV0A) {
      multiplicity = collision.centFV0A();
    }

    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      int track1Sign = track1.sign(); // or track1.charge(), assuming it returns ±1

      histos.fill(HIST("hEta"), track1.eta(), track1.pt());
      histos.fill(HIST("hPhi"), track1.phi(), track1.pt());
      if (track1Sign > 0) { // Positive kaon
        histos.fill(HIST("QAbefore/TPC_Nsigma_all_pos"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        histos.fill(HIST("QAbefore/TOF_Nsigma_all_pos"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        histos.fill(HIST("QAbefore/trkDCAxy_pos"), track1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz_pos"), track1.dcaZ());
        histos.fill(HIST("QAbefore/trkDCAxypt_pos"), track1.dcaXY(), track1.pt());
        histos.fill(HIST("QAbefore/trkDCAzpt_pos"), track1.dcaZ(), track1.pt());
        histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all_pos"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
      } else if (track1Sign < 0) { // Negative kaon
        histos.fill(HIST("QAbefore/TPC_Nsigma_all_neg"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        histos.fill(HIST("QAbefore/TOF_Nsigma_all_neg"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        histos.fill(HIST("QAbefore/trkDCAxy_neg"), track1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz_neg"), track1.dcaZ());
        histos.fill(HIST("QAbefore/trkDCAxypt_neg"), track1.dcaXY(), track1.pt());
        histos.fill(HIST("QAbefore/trkDCAzpt_neg"), track1.dcaZ(), track1.pt());
        histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all_neg"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
      }

      auto track1ID = track1.globalIndex();
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        bool unlike = true;
        bool mix = false;
        if (!ispTdepPID && selectionPID(track1) && selectionPID(track2)) {
          int track1Sign = track1.sign(); // Assuming `charge()` gives +1 or -1

          if (track1Sign > 0) { // Positive kaon
            histos.fill(HIST("QAafter/TPC_Nsigma_all_pos"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/TOF_Nsigma_all_pos"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/trkDCAxy_pos"), track1.dcaXY());
            histos.fill(HIST("QAafter/trkDCAz_pos"), track1.dcaZ());
            histos.fill(HIST("QAafter/trkDCAxypt_pos"), track1.dcaXY(), track1.pt());
            histos.fill(HIST("QAafter/trkDCAzpt_pos"), track1.dcaZ(), track1.pt());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_all_pos"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          } else if (track1Sign < 0) { // Negative kaon
            histos.fill(HIST("QAafter/TPC_Nsigma_all_neg"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/TOF_Nsigma_all_neg"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/trkDCAxy_neg"), track1.dcaXY());
            histos.fill(HIST("QAafter/trkDCAz_neg"), track1.dcaZ());
            histos.fill(HIST("QAafter/trkDCAxypt_neg"), track1.dcaXY(), track1.pt());
            histos.fill(HIST("QAafter/trkDCAzpt_neg"), track1.dcaZ(), track1.pt());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_all_neg"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          }

          fillinvMass(track1, track2, multiplicity, unlike, mix, massKa, massKa);
        }

        if (ispTdepPID && selectionPIDpTdependent(track1) && selectionPIDpTdependent(track2)) {
          int track1Sign = track1.sign(); // Same assumption as above

          if (track1Sign > 0) { // Positive kaon
            histos.fill(HIST("QAafter/TPC_Nsigma_all_pos"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/TOF_Nsigma_all_pos"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/trkDCAxy_pos"), track1.dcaXY());
            histos.fill(HIST("QAafter/trkDCAz_pos"), track1.dcaZ());
            histos.fill(HIST("QAafter/trkDCAxypt_pos"), track1.dcaXY(), track1.pt());
            histos.fill(HIST("QAafter/trkDCAzpt_pos"), track1.dcaZ(), track1.pt());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_all_pos"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          } else if (track1Sign < 0) { // Negative kaon
            histos.fill(HIST("QAafter/TPC_Nsigma_all_neg"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/TOF_Nsigma_all_neg"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("QAafter/trkDCAxy_neg"), track1.dcaXY());
            histos.fill(HIST("QAafter/trkDCAz_neg"), track1.dcaZ());
            histos.fill(HIST("QAafter/trkDCAxypt_neg"), track1.dcaXY(), track1.pt());
            histos.fill(HIST("QAafter/trkDCAzpt_neg"), track1.dcaZ(), track1.pt());
            histos.fill(HIST("QAafter/TOF_TPC_Mapka_all_neg"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          }

          fillinvMass(track1, track2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processSameEvent, "Process Same event", false);
  void processMixedEvent1(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor1 binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor1> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctChecker(c1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctChecker(c2)) {
        continue;
      }
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel1 && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if (additionalEvSel4 && (!c1.selection_bit(aod::evsel::kIsGoodITSLayersAll) || !c2.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
        continue;
      }
      if (additionalEvSel5 && (!c1.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) || !c2.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel6 && (!c1.selection_bit(aod::evsel::kNoCollInRofStandard) || !c2.selection_bit(aod::evsel::kNoCollInRofStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();

      if (fillOccupancy && (occupancy1 > cfgCutOccupancy || occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      float multiplicity;
      multiplicity = c1.centFT0C();
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent1, "Process Mixed event", false);
  void processMixedEvent2(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor2 binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor2> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctChecker(c1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctChecker(c2)) {
        continue;
      }
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel1 && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if (additionalEvSel4 && (!c1.selection_bit(aod::evsel::kIsGoodITSLayersAll) || !c2.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
        continue;
      }
      if (additionalEvSel5 && (!c1.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) || !c2.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel6 && (!c1.selection_bit(aod::evsel::kNoCollInRofStandard) || !c2.selection_bit(aod::evsel::kNoCollInRofStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();

      if (fillOccupancy && (occupancy1 > cfgCutOccupancy || occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      float multiplicity;
      multiplicity = c1.centFT0A();
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent2, "Process Mixed event", false);
  void processMixedEvent3(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor3 binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor3> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctChecker(c1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctChecker(c2)) {
        continue;
      }
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel1 && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if (additionalEvSel4 && (!c1.selection_bit(aod::evsel::kIsGoodITSLayersAll) || !c2.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
        continue;
      }
      if (additionalEvSel5 && (!c1.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) || !c2.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel6 && (!c1.selection_bit(aod::evsel::kNoCollInRofStandard) || !c2.selection_bit(aod::evsel::kNoCollInRofStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();

      if (fillOccupancy && (occupancy1 > cfgCutOccupancy || occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      float multiplicity;
      multiplicity = c1.centFT0M();
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent3, "Process Mixed event", false);
  void processMixedEvent4(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor4 binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor4> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      if (rctCut.requireRCTFlagChecker && !rctChecker(c1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctChecker(c2)) {
        continue;
      }
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel1 && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if (additionalEvSel4 && (!c1.selection_bit(aod::evsel::kIsGoodITSLayersAll) || !c2.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
        continue;
      }
      if (additionalEvSel5 && (!c1.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) || !c2.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel6 && (!c1.selection_bit(aod::evsel::kNoCollInRofStandard) || !c2.selection_bit(aod::evsel::kNoCollInRofStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();

      if (fillOccupancy && (occupancy1 > cfgCutOccupancy || occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      float multiplicity;
      multiplicity = c1.centFV0A();
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          fillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent4, "Process Mixed event", false);
  void processRotEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    float multiplicity{-1};
    if (cfgMultFT0)
      multiplicity = collision.centFT0C();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    for (const auto& track1 : tracks) {

      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("QAbefore/TPC_Nsigma_all"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
      histos.fill(HIST("QAbefore/TOF_Nsigma_all"), track1.tofNSigmaKa(), multiplicity, track1.pt());
      histos.fill(HIST("QAbefore/trkDCAxy"), track1.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAz"), track1.dcaZ());
      histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());

      auto track1ID = track1.globalIndex();
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (track1.sign() * track2.sign() < 0) {
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }
        phiMesonMother = kaonPlus + kaonMinus;
        if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassSame"), multiplicity, phiMesonMother.pt(), phiMesonMother.M());
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (track1.sign() * track2.sign() < 0) {
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              kaonPlus = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            phiMesonMother = kaonPlus + kaonMinus;
            if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h3PhiInvMassRot"), multiplicity, phiMesonMother.pt(), phiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRotEvent, "Process Rot event", false);
  void processMC(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("hMC"), 2);
      return;
    }
    for (const auto& RecCollision : RecCollisions) {
      histos.fill(HIST("hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }
      if (timFrameEvsel && (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        histos.fill(HIST("hMC"), 5);
        continue;
      }
      if (additionalEvSel2 && (!RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup) || !RecCollision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      int occupancy = RecCollision.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
        continue;
      }
      if (std::abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      auto oldindex = -999;
      auto rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (const auto& track1 : rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track1)) {
          continue;
        }
        if (ispTdepPID && !selectionPIDpTdependent(track1)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (const auto& track2 : rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track2)) {
            continue;
          }
          if (ispTdepPID && !selectionPIDpTdependent(track2)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!selectionPair(track1, track2)) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }
          const auto mctrack1 = track1.mcParticle();
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (!(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kKPlus)) {
            continue;
          }
          for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (std::abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (pdgcheck && std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kPhi) {
                continue;
              }
              if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
                continue;
              }
              if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() * track2.sign() < 0) {
                kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              phiMesonMother = kaonPlus + kaonMinus;

              if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
                continue;
              }
              histos.fill(HIST("h1PhiRec1"), phiMesonMother.pt());
              histos.fill(HIST("h2PhiRec2"), phiMesonMother.pt(), centrality);
              histos.fill(HIST("h1Phimassrec"), phiMesonMother.M());
              histos.fill(HIST("h3PhiRec3"), phiMesonMother.pt(), centrality, phiMesonMother.M());
              histos.fill(HIST("Centrec"), centrality);
            }
          }
        }
      }
      // loop over generated particle
      for (const auto& mcParticle : GenParticles) {
        if (std::abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (pdgcheck && mcParticle.pdgCode() != o2::constants::physics::kPhi) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        const size_t kExpectedDaughterCount = 2;

        if (kDaughters.size() != kExpectedDaughterCount) {

          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            if (!genacceptancecut) {
              daughtp = true;
            }
            kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && std::abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            kaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          phiMesonMother = kaonPlus + kaonMinus;
          histos.fill(HIST("h1PhiGen"), phiMesonMother.pt());
          histos.fill(HIST("h2PhiGen2"), phiMesonMother.pt(), centrality);
          histos.fill(HIST("Centgen"), centrality);
          histos.fill(HIST("h1Phimassgen"), phiMesonMother.M());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMC, "Process Reconstructed", false);
  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {

    histos.fill(HIST("hMC"), 0.5);
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 1.5);
    }
    float imp = mcCollision.impactParameter();
    histos.fill(HIST("hImpactParameterGen"), imp);
    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    auto multiplicity = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      int occupancy = collision.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
        continue;
      }
      histos.fill(HIST("hOccupancy1"), occupancy);
      multiplicity = collision.centFT0C();
      histos.fill(HIST("Centgen"), multiplicity);
      histos.fill(HIST("hVtxZgen"), collision.mcCollision().posZ());
      histos.fill(HIST("hImpactParameterGenCen"), imp, multiplicity);

      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
      histos.fill(HIST("hMC"), 2.5);
    }
    selectedEvents.resize(nevts);

    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    histos.fill(HIST("EL1"), imp);
    histos.fill(HIST("EL2"), multiplicity);
    if (reco && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("ES1"), imp);
    histos.fill(HIST("ES2"), multiplicity);
    for (const auto& mcParticle : mcParticles) {
      const double kMaxAcceptedRapidity = 0.5;

      if (std::abs(mcParticle.y()) >= kMaxAcceptedRapidity) {

        continue;
      }
      if (pdgcheck && mcParticle.pdgCode() != o2::constants::physics::kPhi) {
        continue;
      }
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      const size_t kExpectedNumberOfDaughters = 2;

      if (kDaughters.size() != kExpectedNumberOfDaughters) {

        continue;
      }
      auto daughtp = false;
      auto daughtm = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
          daughtp = true;
          kaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
          daughtm = true;
          kaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        }
      }
      if (daughtp && daughtm) {
        phiMesonMother = kaonPlus + kaonMinus;
        histos.fill(HIST("h1PhiGen"), phiMesonMother.pt());
        histos.fill(HIST("h2PhiGen2"), phiMesonMother.pt(), multiplicity);
        histos.fill(HIST("h2PhiGen1"), phiMesonMother.pt(), imp);
        histos.fill(HIST("h1Phimassgen"), phiMesonMother.M());
        histos.fill(HIST("h3PhiGen3"), phiMesonMother.pt(), multiplicity, phiMesonMother.M());
      }
    }
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processGen, "Process Generated", false);
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex || !collision.sel8()) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    auto multiplicity = collision.centFT0C();
    histos.fill(HIST("Centrec"), multiplicity);
    histos.fill(HIST("hVtxZrec"), collision.posZ());
    float imp = collision.mcCollision().impactParameter();
    histos.fill(HIST("hImpactParameterRec"), imp);
    histos.fill(HIST("hImpactParameterRecCen"), imp, multiplicity);
    histos.fill(HIST("ES3"), imp);
    histos.fill(HIST("ES4"), multiplicity);
    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }
      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        if (!(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kKPlus)) {
          continue;
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);

        phiMother = daughter1 + daughter2;
        histos.fill(HIST("h1Phi1massrec"), phiMother.M());
        histos.fill(HIST("h3Phi1Rec3"), phiMother.pt(), multiplicity, phiMother.M());
        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            const double kMaxRapidityCut = 0.5;

            if (std::abs(mothertrack1.y()) >= kMaxRapidityCut) {
              continue;
            }

            if (pdgcheck && std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kPhi) {
              continue;
            }
            if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
              continue;
            }
            if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
              continue;
            }

            histos.fill(HIST("TPC_Nsigma_MC"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("TOF_Nsigma_MC"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            if (track1.sign() * track2.sign() < 0) {
              kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            phiMesonMother = kaonPlus + kaonMinus;

            if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h1PhiRec1"), phiMesonMother.pt());
            histos.fill(HIST("h2PhiRec2"), phiMesonMother.pt(), multiplicity);
            histos.fill(HIST("h1Phimassrec"), phiMesonMother.M());
            histos.fill(HIST("h3PhiRec3"), phiMesonMother.pt(), multiplicity, phiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRec, "Process Reconstructed", false);
  void processSameEventMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!collision.sel8()) {
      return;
    }
    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    if (additionalEvSel5 && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    if (additionalEvSel6 && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    float multiplicity{-1};
    multiplicity = collision.centFT0C();
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      auto track1ID = track1.globalIndex();
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (track1.sign() * track2.sign() < 0) {
          kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }
        phiMesonMother = kaonPlus + kaonMinus;
        if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassSameMC"), multiplicity, phiMesonMother.pt(), phiMesonMother.M());
        histos.fill(HIST("h1Phimasssame"), phiMesonMother.M());
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (track1.sign() * track2.sign() < 0) {
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              kaonPlus = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              kaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            phiMesonMother = kaonPlus + kaonMinus;
            if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h3PhiInvMassRotMC"), multiplicity, phiMesonMother.pt(), phiMesonMother.M());
            histos.fill(HIST("h1Phimassrot"), phiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processSameEventMC, "Process Same event", false);
  void processMixedEventMC(EventCandidatesMC const& recCollisions, TrackCandidatesMC const& RecTracks, aod::McParticles const&)
  {

    auto tracksTuple = std::make_tuple(RecTracks);
    BinningTypeVertexContributor1 binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeVertexContributor1> pairs{binningOnPositions, cfgNoMixedEvents, -1, recCollisions, tracksTuple, &cache};

    for (const auto& [c1, tracks1, c2, tracks2] : pairs) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel1 && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      if (additionalEvSel4 && (!c1.selection_bit(aod::evsel::kIsGoodITSLayersAll) || !c2.selection_bit(aod::evsel::kIsGoodITSLayersAll))) {
        continue;
      }
      if (additionalEvSel5 && (!c1.selection_bit(aod::evsel::kNoCollInTimeRangeStandard) || !c2.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel6 && (!c1.selection_bit(aod::evsel::kNoCollInRofStandard) || !c2.selection_bit(aod::evsel::kNoCollInRofStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy1 > cfgCutOccupancy)) {
        continue;
      }
      if (fillOccupancy && (occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      auto multiplicity = c1.centFT0C();
      for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        histos.fill(HIST("hMC"), 6.5);
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && (!selectionPID(t1) || !selectionPID(t2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(t1) || !selectionPIDpTdependent(t2))) {
          continue;
        }
        if (t1.sign() * t2.sign() < 0) {
          kaonPlus = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          kaonMinus = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa);
        }
        phiMesonMother = kaonPlus + kaonMinus;
        if (std::abs(phiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassMixedMC"), multiplicity, phiMesonMother.pt(), phiMesonMother.M());
        histos.fill(HIST("h1Phimassmix"), phiMesonMother.M());
      }
    }
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEventMC, "Process Mixed event MC", true);
  void processGen1(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    histos.fill(HIST("hMC1"), 0.5);
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC1"), 1.5);
    }
    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    auto multiplicity = -1.0;
    histos.fill(HIST("Centgen1"), multiplicity);
    histos.fill(HIST("hMC1"), 2.5);
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      histos.fill(HIST("hMC1"), 3.5);
      if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 4.5);
      if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 5.5);
      if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 6.5);
      if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 7.5);
      if (additionalEvSel5 && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 8.5);
      if (additionalEvSel6 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 9.5);
      int occupancy = collision.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
        continue;
      }
      histos.fill(HIST("hMC1"), 10.5);
      multiplicity = collision.centFT0C();
      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    histos.fill(HIST("hMC1"), 11.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("hMC1"), 12.5);
    for (const auto& mcParticle : mcParticles) {
      const double kMaxRapidityCut = 0.5;

      if (std::abs(mcParticle.y()) >= kMaxRapidityCut) {
        continue;
      }

      if (mcParticle.pdgCode() != o2::constants::physics::kPhi) {
        continue;
      }
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      const size_t kExpectedNumberOfDaughters = 2;

      if (kDaughters.size() != kExpectedNumberOfDaughters) {

        continue;
      }
      auto daughtp = false;
      auto daughtm = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
          daughtp = true;
        } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
          daughtm = true;
        }
      }
      if (daughtp && daughtm) {
        histos.fill(HIST("h1PhifinalGen"), mcParticle.pt());
        histos.fill(HIST("h2PhifinalGen"), mcParticle.pt(), multiplicity);
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processGen1, "Process Generated", false);
  void processRec1(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex || !collision.sel8()) {
      return;
    }
    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }

    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }

    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    if (additionalEvSel5 && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    if (additionalEvSel6 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    auto multiplicity = collision.centFT0C();
    histos.fill(HIST("Centrec1"), multiplicity);
    histos.fill(HIST("hMC1"), 13.5);
    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }
      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        if (!(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kKPlus)) {
          continue;
        }
        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            const double kMaxRapidityCut = 0.5;

            if (std::abs(mothertrack1.y()) >= kMaxRapidityCut) {
              continue;
            }

            if (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kPhi) {
              continue;
            }
            if (!(selectionPID(track1) && selectionPID(track2))) {
              continue;
            }
            histos.fill(HIST("TPC_Nsigma1_MC"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("TOF_Nsigma1_MC"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("trkDCAxy"), track1.dcaXY());
            histos.fill(HIST("trkDCAz"), track1.dcaZ());
            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1PhiRecsplit1"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            std::array<float, 3> pvec0 = {track1.px(), track1.py(), track1.pz()};
            std::array<float, 3> pvec1 = {track2.px(), track2.py(), track2.pz()};
            std::array<std::array<float, 3>, 2> arrMomrec = {pvec0, pvec1};

            auto motherP = mothertrack1.p();
            auto motherE = mothertrack1.e();
            genMass = std::sqrt(motherE * motherE - motherP * motherP);
            recMass = RecoDecay::m(arrMomrec, array{massKa, massKa});

            histos.fill(HIST("h1PhifinalRec"), mothertrack1.pt());
            histos.fill(HIST("h3PhifinalRec"), mothertrack1.pt(), multiplicity, recMass);
            histos.fill(HIST("h1Phifinalgenmass"), genMass);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRec1, "Process Reconstructed", false);
  void processEvtLossSigLossMC(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {

    // Event loss estimation
    auto impactPar = mcCollision.impactParameter();
    histos.fill(HIST("QAevent/hImpactParameterGen"), impactPar);

    bool isSel = false;
    auto centrality = -999.;
    for (const auto& RecCollision : recCollisions) {
      if (!myEventSelections(RecCollision))
        continue;
      centrality = RecCollision.centFT0C();
      isSel = true;
    }

    if (isSel) {
      histos.fill(HIST("QAevent/hImpactParameterRec"), impactPar);
      histos.fill(HIST("QAevent/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      const double kMaxRapidity = 0.5;

      if (std::abs(mcPart.y()) >= kMaxRapidity || std::abs(mcPart.pdgCode()) != o2::constants::physics::kPhi)

        continue;

      // signal loss estimation
      histos.fill(HIST("QAevent/phigenBeforeEvtSel"), mcPart.pt(), impactPar);
      if (isSel) {
        // signal loss estimation
        histos.fill(HIST("QAevent/phigenAfterEvtSel"), mcPart.pt(), impactPar);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processEvtLossSigLossMC, "Process Signal Loss, Event Loss", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisrun3_PbPb>(cfgc, TaskName{"phianalysisrun3_PbPb"})};
}

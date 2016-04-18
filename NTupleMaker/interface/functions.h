#ifndef NTupleMakerFunctions_h
#define NTupleMakerFunctions_h

#include "TMath.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

const double MuMass = 0.105658367;
const double tauMass = 1.776;
const double electronMass = 0;
const double muonMass = 0.10565837;
const double pionMass = 0.1396;
using namespace std;


vector <int > runvect_;
vector <int > lumivect_;
vector <int > eventvect_;
struct Event {
    std::string name;
    int run;
    int lumi;
    int eventrn;
};

std::vector<std::string>
split(const std::string &s, char delim = ':')
{

    std::vector<std::string> elems;

    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
            elems.push_back(item);
    }
    return elems;
}
  


bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


int binNumber(float x, int nbins, float * bins) {

  int binN = -1;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double PtoPhi(double Px, double Py) {
  return TMath::ATan2(Py,Px);
}

double PtoPt(double Px, double Py) {
  return TMath::Sqrt(Px*Px+Py*Py);
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;

}

double mT (TLorentzVector v1, TLorentzVector Metv){

  double dPhimT=dPhiFrom2P( v1.Px(), v1.Py(), Metv.Px(),  Metv.Py() );

  double MT = TMath::Sqrt(2*v1.Pt()*Metv.Pt()*(1-TMath::Cos(dPhimT)));
  return MT;
}



double DeltaPhi(TLorentzVector METV, TLorentzVector LeptonV){
  TLorentzVector Ws = METV + LeptonV;

  //Delta phi between W and Lep
  //standard root defintion (+ fabs)takes care of getting values between 0 and pi
  double DelPhiWLep = fabs(Ws.DeltaPhi(LeptonV));
  //alternative definiton with the same result, if you want to cross check
  Double_t DelPhiWLepAlt = (Ws.Phi() - LeptonV.Phi());
  if (DelPhiWLepAlt > TMath::Pi()) DelPhiWLepAlt -= 2*TMath::Pi();
  if (DelPhiWLepAlt <= -TMath::Pi()) DelPhiWLepAlt += 2*TMath::Pi();
  DelPhiWLepAlt = fabs(DelPhiWLepAlt);
		
  return DelPhiWLep;

}



bool electronMvaIdTight(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.73) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.57) passed = true;
  }
  else {
    if (mva>0.05) passed = true;
  }

  return passed;

}

bool electronMvaIdLoose(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.35) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.20) passed = true;
  }
  else {
    if (mva>-0.52) passed = true;
  }

  return passed;

}

bool electronMvaIdWP80(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.253;
    else 
      passed = mva > 0.965;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > 0.081;
    else
      passed = mva > 0.917;
  }
  else {
    if (pt<10)
      passed = mva > -0.081;
    else
      passed = mva > 0.683;
  }

  return passed;

}

bool electronMvaIdWP90(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.483;
    else 
      passed = mva > 0.933;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > -0.267;
    else
      passed = mva > 0.825;
  }
  else {
    if (pt<10)
      passed = mva > -0.323;
    else
      passed = mva > 0.337;
  }

  return passed;

}


bool electronVetoTight(double SuperClusterEta, double eta, double phi, double full5, double hOverE, double d0, double dZ, double ooE, double pfISO, double nMissing, bool convVeto) {


  bool passed = false;


  if (fabs(SuperClusterEta) <= 1.479 ){
	
    if (  fabs(eta)<0.013625 && fabs(phi)< 0.230374 && fabs(full5) < 0.011586 && hOverE < 0.181130 && fabs(d0) < 0.094095 && fabs(dZ) < 0.713070 && ooE < 0.295751 && pfISO < 0.158721 && nMissing <= 2 && convVeto) 
      passed = true;
  }

  else if (      1.479  < fabs(SuperClusterEta)  && fabs(SuperClusterEta) < 2.5 ){
 
    if ( fabs(eta)<0.011932 && fabs(phi)< 0.255450 && fabs(full5) < 0.031849 && hOverE < 0.223870 && fabs(d0) < 0.342293 && fabs(dZ) < 0.953461 && ooE < 0.155501 && pfISO < 0.177032 && nMissing <= 3 && convVeto) 

      passed = true;
  }
  else  
    passed = false;

  return passed;
}


void ComputeMetFromHadRecoil(float Hparal,
			     float Hperp,
			     float genVPx, 
			     float genVPy,
			     float visVPx,
			     float visVPy,
			     float & metX,
			     float & metY) {
  
  float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
  float unitX = genVPx/genVPt;
  float unitY = genVPy/genVPt;

  float unitPhi = TMath::ATan2(unitY,unitX);
  float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
  float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

  float det = unitX*unitPerpY - unitY*unitPerpX;
  float Hx = (Hparal*unitPerpY - Hperp*unitY)/det;
  float Hy = (Hperp*unitX - Hparal*unitPerpX)/det;

  metX = -Hx - visVPx;
  metY = -Hy - visVPy;

}

void ComputeHadRecoilFromMet(float metX,
			     float metY,
			     float genVPx, 
			     float genVPy,
			     float visVPx,
			     float visVPy,
			     float & Hparal,
			     float & Hperp) {

  float genVPt = TMath::Sqrt(genVPx*genVPx+genVPy*genVPy);
  float unitX = genVPx/genVPt;
  float unitY = genVPy/genVPt;

  float unitPhi = TMath::ATan2(unitY,unitX);
  float unitPerpX = TMath::Cos(unitPhi+0.5*TMath::Pi());
  float unitPerpY = TMath::Sin(unitPhi+0.5*TMath::Pi());

  float Hx = -metX - visVPx;
  float Hy = -metY - visVPy;

  Hparal = Hx*unitX + Hy*unitY;
  Hperp = Hx*unitPerpX + Hy*unitPerpY;

}


struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject, myobjectX;


namespace genTools{
  TLorentzVector genZ(const AC1B& analysisTree){
    TLorentzVector genZ; genZ.SetXYZM(0,0,0,91.2);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      if (analysisTree.genparticles_pdgid[igen]==23||analysisTree.genparticles_pdgid[igen]==22) {
	if (analysisTree.genparticles_fromHardProcess[igen])
	  genZ.SetXYZT(analysisTree.genparticles_px[igen],
		       analysisTree.genparticles_py[igen],
		       analysisTree.genparticles_pz[igen],
		       analysisTree.genparticles_e[igen]);
      }
    }
    return genZ;
  }

  TLorentzVector genL(const AC1B& analysisTree){
    TLorentzVector genL; genL.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
     
    bool isMuon = 0;
    bool isElectron = 0;
    bool isChargedLepton = 0;
    bool isNeutrino = 0;
    bool fromHardProcessFinalState = 0;
    bool isDirectHardProcessTauDecayProduct = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
      isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
      isChargedLepton = isMuon || isElectron;
      isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	fabs(analysisTree.genparticles_pdgid[igen])==14||
	fabs(analysisTree.genparticles_pdgid[igen])==16;
      fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen] && analysisTree.genparticles_status[igen]==1;
      isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

      /*if((fromHardProcessFinalState && isChargedLepton) || (isDirectHardProcessTauDecayProduct && (!isNeutrino)))
	genL += genPart;
    }

    return genL;*/

      if(fromHardProcessFinalState && isChargedLepton)
	genL += genPart;
    }
    
    bool fromHardProcess = 0;
    bool isLastCopy = 0;
    for (unsigned int itau=0; itau<analysisTree.gentau_count; ++itau) {
      genPart.SetXYZT(analysisTree.gentau_visible_px[itau],
		      analysisTree.gentau_visible_py[itau],
		      analysisTree.gentau_visible_pz[itau],
		      analysisTree.gentau_visible_e[itau]);
      
      fromHardProcess = analysisTree.gentau_fromHardProcess[itau];
      isLastCopy = analysisTree.gentau_isLastCopy[itau];

      if (fromHardProcess && isLastCopy)
	genL += genPart;
    }
      
    return genL;
      
  }

  TLorentzVector genNu(const AC1B& analysisTree){
    TLorentzVector genNu; genNu.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    bool isNeutrino = 0;
    bool isPrompt = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	fabs(analysisTree.genparticles_pdgid[igen])==14||
	fabs(analysisTree.genparticles_pdgid[igen])==16;
      isPrompt = analysisTree.genparticles_isPrompt[igen]||
	analysisTree.genparticles_isPromptTauDecayProduct[igen];
      
      if (analysisTree.genparticles_status[igen]==1&&isPrompt) {
	if (isNeutrino) 
	  genNu += genPart;
      }
    }

    return genNu;
  }

  TLorentzVector genV(const AC1B& analysisTree){
    TLorentzVector genV; genV.SetXYZM(0,0,0,0);
    TLorentzVector genPart; genPart.SetXYZM(0,0,0,0);
    
    bool isMuon = 0;
    bool isElectron = 0;
    bool isChargedLepton = 0;
    bool isNeutrino = 0;
    bool fromHardProcessFinalState = 0;
    bool isDirectHardProcessTauDecayProduct = 0;
    
    for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
      genPart.SetXYZT(analysisTree.genparticles_px[igen],
		      analysisTree.genparticles_py[igen],
		      analysisTree.genparticles_pz[igen],
		      analysisTree.genparticles_e[igen]);
      
      isMuon = fabs(analysisTree.genparticles_pdgid[igen])==13;
      isElectron = fabs(analysisTree.genparticles_pdgid[igen])==11;
      isChargedLepton = isMuon || isElectron;
      isNeutrino = fabs(analysisTree.genparticles_pdgid[igen])==12||
	fabs(analysisTree.genparticles_pdgid[igen])==14||
	fabs(analysisTree.genparticles_pdgid[igen])==16;
      fromHardProcessFinalState = analysisTree.genparticles_fromHardProcess[igen] && analysisTree.genparticles_status[igen]==1;
      isDirectHardProcessTauDecayProduct = analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen];

      /*if((fromHardProcessFinalState && (isChargedLepton || isNeutrino)) || isDirectHardProcessTauDecayProduct)
	genV += genPart;
    }

    if (genV.Pt()<0.1)
      genV.SetXYZM(0.1,0.1,0.,0.);
    
      return genV;*/


      if(fromHardProcessFinalState && (isChargedLepton || isNeutrino))
	genV += genPart;
    }
    
    bool fromHardProcess = 0;
    bool isFirstCopy = 0;
    for (unsigned int itau=0; itau<analysisTree.gentau_count; ++itau) {
      genPart.SetXYZT(analysisTree.gentau_px[itau],
		      analysisTree.gentau_py[itau],
		      analysisTree.gentau_pz[itau],
		      analysisTree.gentau_e[itau]);
      
      fromHardProcess = analysisTree.gentau_fromHardProcess[itau];
      isFirstCopy = analysisTree.gentau_isFirstCopy[itau];

      if (fromHardProcess && isFirstCopy)
	genV += genPart;
    }
   
    if (genV.Pt()<0.1)
      genV.SetXYZM(0.1,0.1,0.,0.);
    
    return genV;  
  }

  int nJetsHad(const AC1B& analysisTree){
    int njetshad = 0;
    bool isChargedLepton = 0;
    bool fromHardProcess = 0;
    
    float genEta = 0.;
    float genPhi = 0.;
    
    for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

      if (analysisTree.pfjet_pt[jet]<=30.) continue;
      float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
      if (absJetEta >= 4.7) continue;
      
      // jetId
      float energy = analysisTree.pfjet_e[jet];
      energy *= analysisTree.pfjet_energycorr[jet];
      float chf = analysisTree.pfjet_chargedhadronicenergy[jet]/energy;
      float nhf = analysisTree.pfjet_neutralhadronicenergy[jet]/energy;
      float phf = analysisTree.pfjet_neutralemenergy[jet]/energy;
      float elf = analysisTree.pfjet_chargedemenergy[jet]/energy;
      float muf = analysisTree.pfjet_muonenergy[jet]/energy;
      float chm = analysisTree.pfjet_chargedmulti[jet];
      float nm = analysisTree.pfjet_neutralmulti[jet];
      float npr = analysisTree.pfjet_chargedmulti[jet] + analysisTree.pfjet_neutralmulti[jet];
      //bool isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      bool isPFJetId = false;
      if (absJetEta<=3.0)
	isPFJetId = (nhf < 0.99 && phf < 0.99 && npr > 1) && (absJetEta>2.4 || (chf>0 && chm > 0 && elf < 0.99));
      else
	isPFJetId = phf < 0.9 && nm > 10;
      //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99 && muf < 0.8) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      //isPFJetId = (npr>1 && phf<0.99 && nhf<0.99) && (absJetEta>3.0 || (elf<0.99 && chf>0 && chm>0));
      
      if (!isPFJetId) continue;

      int overlap = 0;
      
      for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {	
	isChargedLepton = ( fabs(analysisTree.genparticles_pdgid[igen])==11 ||
			    fabs(analysisTree.genparticles_pdgid[igen])==13);
	fromHardProcess = analysisTree.genparticles_fromHardProcess[igen];

	if (!isChargedLepton) continue;
	if (!fromHardProcess) continue;

	genEta = PtoEta( analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen]); 
	genPhi = PtoPhi( analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen]);
	
	float dR = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			  genEta, genPhi);
	if (dR<=0.5) continue;
	
	overlap = 1;
	break;
      }

      if (overlap == 1) continue;

      for (unsigned int itau=0; itau<analysisTree.gentau_count; ++itau) {
	fromHardProcess = analysisTree.gentau_fromHardProcess[itau];
	
      	if (!fromHardProcess) continue;

	genEta = PtoEta( analysisTree.gentau_visible_px[itau], analysisTree.gentau_visible_py[itau], analysisTree.gentau_visible_pz[itau]);
	genPhi = PtoPhi( analysisTree.gentau_visible_px[itau], analysisTree.gentau_visible_py[itau]);
			 	
	float dR = deltaR(analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet],
			  genEta, genPhi);
	if (dR<=0.5) continue;
	
	overlap = 1;
	break;
      }

      if (overlap == 1) continue;
      
      njetshad++;
    }

    return njetshad;
  } 
  
  enum RecoilCorrectionsMethod{QuantileRemap=1, MeanResolution};

  int RecoilCorrections( RecoilCorrector& corr, int method,
			 float met, float metphi,
			 float vx, float vy,
			 float lx, float ly,
			 int njets,
			 float& metcorr, float& metphicorr ){
    float metx = met*TMath::Cos(metphi);
    float mety = met*TMath::Sin(metphi);
    float metcorrx = metx;
    float metcorry = mety;
    if (method == 1)
      corr.Correct(metx, mety, vx, vy, lx, ly, njets, metcorrx, metcorry);
    else if(method == 2)
      corr.CorrectByMeanResolution(metx, mety, vx, vy, lx, ly, njets, metcorrx, metcorry);

    metcorr = TMath::Sqrt(metcorrx * metcorrx + metcorry * metcorry);
    metphicorr = TMath::ATan2(metcorry, metcorrx);
            
    return method;
  }

  float mt( float lpt, float lphi, float met, float metphi){
    return sqrt(2*lpt*met*(1.-TMath::Cos(lphi-metphi)));
  }

  float pzetamiss( float zx, float zy, float met, float metphi){
    return zx*met*TMath::Cos(metphi)+zy*met*TMath::Sin(metphi);
  }
}
#endif

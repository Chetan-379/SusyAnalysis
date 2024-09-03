#define NtupleVarsTProxy_cxx

#include "NtupleVarsTProxy.h"
#include <TH2.h>
#include <TStyle.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myLV;


//#pragma link C++ class NtupleVarsTProxy+;
double NtupleVarsTProxy::DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double NtupleVarsTProxy::DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}
double NtupleVarsTProxy::TransMass(double phi1, double phi2, double pt1, double pt2){
  double l_=-1;
  if(phi1 == -9999.0 && phi2 == -9999.0 && pt1==-9999.0 && pt2==-9999.0 ) return l_;
  
  double dphi = DeltaPhi(phi1, phi2);
  double Cos= 1-cos(dphi);
  return std::sqrt(2*pt1*pt2*Cos);
}
double NtupleVarsTProxy::MinDr(TLorentzVector v1,vector<TLorentzVector> v2)
{
  double dr = 60;
  for(int j=0;j<v2.size();j++)
    { if(dr>=v1.DeltaR(v2[j]))
	{ dr = v1.DeltaR(v2[j]);}
    }
  return dr;
}

double NtupleVarsTProxy::MinDr_myLV(myLV v1,vector<myLV> v2)
{
  double dr = 60;
  for(int j=0;j<v2.size();j++)
    { if(dr>=DeltaR(v1.Eta(),v1.Phi(),v2[j].Eta(),v2[j].Phi()))
	{ dr = DeltaR(v1.Eta(),v1.Phi(),v2[j].Eta(),v2[j].Phi());}
    }
  return dr;
}

double NtupleVarsTProxy::MinDr2(vector<TLorentzVector> v1,TLorentzVector v2)
{
  double dr = 60;
  for(int j=0;j<v1.size();j++)
    { if(dr>=v1[j].DeltaR(v2))
	{ dr = v1[j].DeltaR(v2);}
    }
  return dr;
}

double NtupleVarsTProxy::MinDr2_myLV(vector<myLV> v1,myLV v2)
{
  double dr = 60;
  for(int j=0;j<v1.size();j++)
    {if(dr>=DeltaR(v1[j].Eta(),v1[j].Phi(),v2.Eta(),v2.Phi())) dr=DeltaR(v1[j].Eta(),v1[j].Phi(),v2.Eta(),v2.Phi());
    }
  return dr;
}

double NtupleVarsTProxy::getCrossSection(std::string process_name)
{
  std::map<std::string, float>::iterator it = cross_sectionValues.find(process_name);
  //std::map<std::string, float>::iterator it = cross_sectionValues.begin();
  //cout << it << endl;
  if(it !=cross_sectionValues.end()){
    return it->second;
  }
  else {
    return 0;
  }
}

double NtupleVarsTProxy::getEventWeight(TString s_Process, double xsec)
{
  //double xsec = getCrossSection(process_name);
  //double xsec = constXsec;
  //if(xsec>0)
  if(s_Process.Contains("2018.WGJets_MonoPhoton_PtG-40to130UL") ||
     s_Process.Contains("2018.WGJets_MonoPhoton_PtG-130UL")||
     s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-40to130UL") ||
     s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-130UL") ||
     s_Process.Contains("2017.WGJets_MonoPhoton_PtG-40to130UL")||
     s_Process.Contains("2017.WGJets_MonoPhoton_PtG-130UL")||
     s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-130UL")||
     s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-40to130UL")) {
    //std::cout << "getEventWeight " << xsec*59.83*1000.0 << std::endl;
    return xsec * 59.83 * 1000.0;}
  else
    return Weight * 59.83 * 1000.0;
}

void NtupleVarsTProxy::sortTLorVec(vector<TLorentzVector> *vec){
  TLorentzVector temp;
  for(int i=1;i<vec->size();i++){
    for(int j=i;j<vec->size();j++){
      if( (*vec)[i-1].Pt() < (*vec)[j].Pt() ){
	temp = (*vec)[i-1];
	(*vec)[i-1] = (*vec)[j];
	(*vec)[j] = temp;
      }
    }
  }
}

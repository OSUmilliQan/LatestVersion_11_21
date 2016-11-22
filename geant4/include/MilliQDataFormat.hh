#ifndef MilliQDataFormat_h
#define MilliQDataFormat_h 1

#include "g4root.hh"
#include "G4Types.hh"
#include "MilliQEventAction.hh"

#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

template<typename T> std::vector<T> ptree_array(const boost::property_tree::ptree pt, const std::string sEntry, const G4double unit) {

  const std::string sConfig = pt.get<std::string>(sEntry);

  std::vector<T> v;
  std::stringstream ss(sConfig);
  std::string item;

  while(std::getline(ss, item, ',')) v.push_back(boost::lexical_cast<T>(item) * unit);

  return v;
}

class MilliQDataFormat {
 public:

  MilliQDataFormat();
  virtual ~MilliQDataFormat();

  virtual void CreateDataMembers() {
    CreateAnalysisManager();
    CreateHistograms();
    CreateNtuples();
  }

  enum NtupleIds { kAll, kCAEN, kDEDX, kNumNtupleIds };

  inline void SetVerboseLevel(G4int v) { fVerboseLevel = v; }
  inline G4int GetVerboseLevel()       { return fVerboseLevel; }

 private:

  virtual void CreateAnalysisManager();
  virtual void CreateHistograms();
  virtual void CreateNtuples();

  G4AnalysisManager* analysisManager;

  G4int fVerboseLevel;
};

#endif

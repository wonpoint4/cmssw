#ifndef DQM_RPCMonitorClient_DQMDaqInfo_H
#define DQM_RPCMonitorClient_DQMDaqInfo_H

// system include files
#include <iostream>
#include <fstream>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "CondFormats/RunInfo/interface/RunSummary.h"
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"

class RPCDaqInfo : public DQMEDHarvester {
public:
  explicit RPCDaqInfo(const edm::ParameterSet &);
  ~RPCDaqInfo() override;

protected:
  void beginJob() override;
  void dqmEndLuminosityBlock(DQMStore::IBooker &,
                             DQMStore::IGetter &,
                             edm::LuminosityBlock const &,
                             edm::EventSetup const &) override;       //performed in the endLumi
  void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override;  //performed in the endJob

private:
  void myBooker(DQMStore::IBooker &);

  edm::ESGetToken<RunInfo, RunInfoRcd> runInfoToken_;

  bool init_;

  MonitorElement *DaqFraction_;
  MonitorElement *DaqMap_;
  constexpr static int kNWheels = 5;
  MonitorElement *daqWheelFractions[kNWheels];
  constexpr static int kNDisks = 10;
  MonitorElement *daqDiskFractions[kNDisks];

  std::pair<int, int> FEDRange_;

  int numberOfDisks_, NumberOfFeds_;
};

#endif

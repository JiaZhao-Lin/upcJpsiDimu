// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIzy35dIzaochendIcmsJpsiDimudIjpsiAnaCodedIanaDatadIanaEvt_C_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/zy35/zaochen/cmsJpsiDimu/jpsiAnaCode/anaData/./anaEvt.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_anaEvt_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./anaEvt.C",
0
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/include",
"/cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/include/Geant4",
"/cvmfs/sft.cern.ch/lcg/releases/jsonmcpp/3.9.1-72770/x86_64-centos7-gcc10-opt/include",
"/cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/src/cpp",
"/cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/include",
"/cvmfs/sft.cern.ch/lcg/releases/Python/3.8.6-3199b/x86_64-centos7-gcc10-opt/include/python3.8",
"/cvmfs/sft.cern.ch/lcg/releases/R/3.6.3-0846b/x86_64-centos7-gcc10-opt/lib64/R/include",
"/cvmfs/sft.cern.ch/lcg/releases/R/3.6.3-0846b/x86_64-centos7-gcc10-opt/lib64/R/library/RInside/include",
"/cvmfs/sft.cern.ch/lcg/releases/R/3.6.3-0846b/x86_64-centos7-gcc10-opt/lib64/R/library/Rcpp/include",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/etc/",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/etc//cling",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/include/",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/include",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.24.00-a725e/x86_64-centos7-gcc10-opt/include/",
"/home/zy35/zaochen/cmsJpsiDimu/jpsiAnaCode/anaData/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "anaEvt_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "anaEvt_C_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./anaEvt.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"anaEvt", payloadCode, "@",
"bookHistos", payloadCode, "@",
"goodMuPair", payloadCode, "@",
"grabNeutronNum", payloadCode, "@",
"hCen_Final", payloadCode, "@",
"hCenvsTrig", payloadCode, "@",
"hDeltaEtavsDetaPhi4Cosmic", payloadCode, "@",
"hDeltaPtvsDeltaEtavsRap", payloadCode, "@",
"hHFMinusvsHFPlusvsCen", payloadCode, "@",
"hHFMinusvsHFPlusvsCen_Sel", payloadCode, "@",
"hMuPtvsEtavsRap", payloadCode, "@",
"hMvsAsyPhivsRap", payloadCode, "@",
"hMvsAsyPhivsRap_NeuDir", payloadCode, "@",
"hMvsPtvsRap", payloadCode, "@",
"hMvsPtvsRap_NeuDir", payloadCode, "@",
"hMvsPtvsRap_WS", payloadCode, "@",
"hNegMuEtavsPosMuEtavsRap", payloadCode, "@",
"hNegMuPtvsPosMuPtvsRap", payloadCode, "@",
"hNeuNumMinusvsNeuNumPlus", payloadCode, "@",
"hNtrkHP_2SoftMuons", payloadCode, "@",
"hNtrkHPvsNtrkofflinevsCen", payloadCode, "@",
"hNtrkHPvsNtrkofflinevsCen_Sel", payloadCode, "@",
"hPreScalevsTrig", payloadCode, "@",
"hRawCen", payloadCode, "@",
"hTrigMuPtvsEtavsRap", payloadCode, "@",
"hVzvsVyvsVx", payloadCode, "@",
"hVzvsVyvsVx_Sel", payloadCode, "@",
"hZDCMinusvsZDCPlus", payloadCode, "@",
"hZDCMinusvsZDCPlus_LS", payloadCode, "@",
"hZDCMinusvsZDCPlus_Only2MuTrk", payloadCode, "@",
"hZDCMinusvsZDCPlus_Only2MuTrk_LS", payloadCode, "@",
"hZDCMinusvsZDCPlus_Sel", payloadCode, "@",
"hZDCMinusvsZDCPlus_Sel_LS", payloadCode, "@",
"hZDCvsNeuNum", payloadCode, "@",
"hnEvts", payloadCode, "@",
"shiftDeltaPhi", payloadCode, "@",
"shiftToPi", payloadCode, "@",
"writeHistos", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("anaEvt_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_anaEvt_C_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_anaEvt_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_anaEvt_C_ACLiC_dict() {
  TriggerDictionaryInitialization_anaEvt_C_ACLiC_dict_Impl();
}

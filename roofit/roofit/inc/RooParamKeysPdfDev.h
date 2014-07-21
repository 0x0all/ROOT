/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooParamKeysPdfDev.h,v 1.10 2007/05/11 09:13:07 verkerke Exp $
 * Authors:                                                                  *
 *   GR, Gerhard Raven,   UC San Diego,        raven@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_PARAM_KEYS_DEV
#define ROO_PARAM_KEYS_DEV

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooSetProxy.h"
#include "RooObjCacheManager.h"

class RooRealVar;
class RooDataHist;
class RooHistPdf;
class RooChangeTracker;

class TH1;

class RooParamKeysPdfDev : public RooAbsPdf {
public:
  enum Mirror { NoMirror, MirrorLeft, MirrorRight, MirrorBoth,
		MirrorAsymLeft, MirrorAsymLeftRight,
		MirrorAsymRight, MirrorLeftAsymRight,
		MirrorAsymBoth };
  RooParamKeysPdfDev() ;
  RooParamKeysPdfDev(
    const char *name, const char *title,
    RooAbsReal& x, RooAbsReal& deltax, 
    RooDataSet& data, Mirror mirror= NoMirror, Double_t rho=1, Int_t nPoints=1000
  );
  RooParamKeysPdfDev(
    const char *name, const char *title,
    RooAbsReal& x, RooAbsReal& deltax, double centralValue, RooAbsReal& multiplicativeShift,
    RooDataSet& data, Mirror mirror= NoMirror, Double_t rho=1, Int_t nPoints=1000
  );
  RooParamKeysPdfDev(
    const char *name, const char *title,
    const RooArgSet& xSet, const RooArgList& deltaxList, 
    const std::vector<double>& centralValues, 
    const RooArgList& multiplicativeShiftList, 
    RooDataSet& data,
    const std::vector<int>& mirror, const std::vector<Double_t>& rho,
    Int_t nPoints=1000
  );
  
  RooParamKeysPdfDev(const RooParamKeysPdfDev& other, const char* name=0);
  virtual TObject* clone(const char* newname) const {return new RooParamKeysPdfDev(*this,newname); }
  virtual ~RooParamKeysPdfDev();

  void LoadDataSet( RooDataSet& data);

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  RooAbsReal* createIntegral(const RooArgSet& iset, const RooArgSet* nset=0, const RooNumIntConfig* cfg=0, const char* rangeName=0) const ;  

protected:

  Int_t _nDim;
  
  //RooRealProxy _x ;
  RooArgSet   _histObsList ; // List of observables defining dimensions of histogram
  RooSetProxy _xSet;         // List of observables

  TIterator*  _histObsIter ; //! 
  TIterator*  _xSetIter;     //!


  RooDataHist*   _dataHist;     // pointer to underlying histogram

  //RooRealProxy _deltax ;
  RooListProxy _deltaxList;

  double* _centralValues; //[_nDim]
  //RooRealProxy _multiplicativeShift;
  RooListProxy _multiplicativeShiftList;

  Double_t evaluate() const;

private:

  
  Double_t evaluateFull(Double_t* x) const;

  Int_t      _nEvents; 
  Double_t ** _dataPts;   //[_nEvents][_nDim]
  Double_t *  _dataWgts;  //[_nEvents]
  Double_t ** _weights;   //[_nEvents][_nDim]
  Double_t   _sumWgt;
  
  Int_t _nPoints;
  Int_t _nPointsAll;
  
  Double_t* _lookupTable; // [_nPointsAll]
  
  Double_t fgaus(Double_t* x, Double_t* sigma) const;

  Bool_t *_mirrorLeft;  // [_nDim]
  Bool_t *_mirrorRight; // [_nDim]
  Bool_t *_asymLeft; // [_nDim]
  Bool_t *_asymRight; // [_nDim]

  // cached info on variable

  //Char_t _varName[128];
  std::string* _varNames; //! do not persistify
  
  //Double_t _lo, _hi, _binWidth;
  Double_t* _loArr; // [_nDim]
  Double_t* _hiArr; // [_nDim]
  Double_t* _binWidths; // [_nDim]

  Double_t* _rho; // [_nDim]
  
  TMatrixD* _rotMat;
  
  ClassDef(RooParamKeysPdfDev,4) // One-dimensional non-parametric kernel estimation p.d.f.
};

#endif

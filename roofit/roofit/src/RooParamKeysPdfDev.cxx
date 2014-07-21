/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooParamKeysPdfDev.cxx 44644 2012-06-11 11:47:21Z moneta $
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
#include "RooFit.h"

#include <math.h>
#include "Riostream.h"
#include "TMath.h"
#include "TMatrixDSymEigen.h"

#include "RooParamKeysPdfDev.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooLinearVar.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooChangeTracker.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include "RooHistPdf.h"
#include "RooDataHist.h"

using namespace std;

ClassImp(RooParamKeysPdfDev)


//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooParamKeysPdfDev implements a one-dimensional kernel estimation p.d.f which model the distribution
// of an arbitrary input dataset as a superposition of Gaussian kernels, one for each data point,
// each contributing 1/N to the total integral of the p.d.f.
// <p>
// If the 'adaptive mode' is enabled, the width of the Gaussian is adaptively calculated from the
// local density of events, i.e. narrow for regions with high event density to preserve details and
// wide for regions with log event density to promote smoothness. The details of the general algorithm
// are described in the following paper: 
// <p>
// Cranmer KS, Kernel Estimation in High-Energy Physics.  
//             Computer Physics Communications 136:198-207,2001 - e-Print Archive: hep ex/0011057
// <p>
// END_HTML
//


//_____________________________________________________________________________
RooParamKeysPdfDev::RooParamKeysPdfDev() : _dataHist(0), _centralValues(0),
  _nEvents(0), _dataPts(0), _dataWgts(0),_weights(0), _sumWgt(0.), _nPoints(1000), _lookupTable(0),
  _mirrorLeft(0), _mirrorRight(0), _asymLeft(0), _asymRight(0), 
  _loArr(0), _hiArr(0), _binWidths(0), _rho(0), _rotMat(0)
{ 

  // coverity[UNINIT_CTOR]
  //cout<<"empty c'tor"<<endl;
  _xSetIter = _xSet.createIterator();
  _histObsIter = _histObsList.createIterator();

}


//_____________________________________________________________________________
RooParamKeysPdfDev::RooParamKeysPdfDev(const char *name, const char *title,
				       RooAbsReal& x, RooAbsReal& deltax, RooDataSet& data,
				       Mirror mirror, Double_t rho, Int_t nPoints) :
  RooAbsPdf(name,title),
  _nDim(1),
  _xSet("xSet","Dependent List", this),
  _dataHist(0),
  _deltaxList("deltaxList","Dependent List",this),
  _centralValues(new double[1]),
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),
  _nPoints(nPoints),
  _nPointsAll(TMath::Power(_nPoints+1,_nDim)),
  _mirrorLeft(new Bool_t[_nDim]),
  _mirrorRight(new Bool_t[_nDim]),
  _asymLeft(new Bool_t[_nDim]),
  _asymRight(new Bool_t[_nDim]),
  _varNames(new std::string[_nDim]),
  _loArr(new Double_t[_nDim]),
  _hiArr(new Double_t[_nDim]),
  _binWidths(new Double_t[_nDim]),
  _rho(new Double_t[_nDim]),
  _rotMat(new TMatrixD(_nDim,_nDim))
{
  //cout<<"1-D c'tor"<<endl;

  _centralValues[0]=0.0;
  _rho[0]=rho;

  _xSet.add(x);
  _xSetIter = _xSet.createIterator();

  _histObsList.addClone(_xSet);
  _histObsIter = _histObsList.createIterator();

  ////cout<<"RooParamKeysPdfDev c'tor with xSet="<<_xSet<<endl;
  //_xSet.Print("V");

  _deltaxList.add(deltax);

  // cache stuff about x
  _varNames[0]=x.GetName();
  
  RooRealVar* real= (RooRealVar*)(&x);

  Double_t lo=real->getMin();
  Double_t hi=real->getMax();

  _loArr[0]=lo;
  _hiArr[0]=hi;
  _binWidths[0]=(hi-lo)/(_nPoints-1);

  _mirrorLeft[0] =mirror==MirrorLeft || mirror==MirrorBoth || mirror==MirrorLeftAsymRight;
  _mirrorRight[0]=mirror==MirrorRight || mirror==MirrorBoth || mirror==MirrorAsymLeftRight;
  _asymLeft[0]   =mirror==MirrorAsymLeft || mirror==MirrorAsymLeftRight || mirror==MirrorAsymBoth;
  _asymRight[0]  =mirror==MirrorAsymRight || mirror==MirrorLeftAsymRight || mirror==MirrorAsymBoth;

  _lookupTable=new Double_t[_nPointsAll];

  // form the lookup table
  LoadDataSet(data);

}


//_____________________________________________________________________________
RooParamKeysPdfDev::RooParamKeysPdfDev(const char *name, const char *title,
				       RooAbsReal& x, RooAbsReal& deltax, 
				       double centralValue, RooAbsReal& multiplicativeShift, 
				       RooDataSet& data,
				       Mirror mirror, Double_t rho, Int_t nPoints) :
  RooAbsPdf(name,title),
  _nDim(1),
  _xSet("xSet","Dependent List", this),
  _dataHist(0),
  _deltaxList("deltaxList","Dependent List", this),
  _centralValues(new double[_nDim]),
  _multiplicativeShiftList("multiplicativeShiftList","Dependent",this),  
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),
  _nPoints(nPoints),
  _nPointsAll(TMath::Power(_nPoints+1,_nDim)),
  _mirrorLeft(new Bool_t[_nDim]),
  _mirrorRight(new Bool_t[_nDim]),
  _asymLeft(new Bool_t[_nDim]),
  _asymRight(new Bool_t[_nDim]),
  _varNames(new std::string[_nDim]),
  _loArr(new Double_t[_nDim]),
  _hiArr(new Double_t[_nDim]),
  _binWidths(new Double_t[_nDim]),
  _rho(new Double_t[_nDim]),
  _rotMat(new TMatrixD(_nDim,_nDim))
{
  //cout<<"1-D c'tor"<<endl;

  _centralValues[0]=centralValue;
  _rho[0]=rho;

  _xSet.add(x);
  _xSetIter = _xSet.createIterator();

  _histObsList.addClone(_xSet);
  _histObsIter = _histObsList.createIterator();

  //cout<<"RooParamKeysPdfDev c'tor with xSet="<<_xSet<<endl;
  //_xSet.Print("V");

  _deltaxList.add(deltax);
  _multiplicativeShiftList.add(multiplicativeShift);

  // cache stuff about x
  _varNames[0]=x.GetName();
  
  RooRealVar* real= dynamic_cast<RooRealVar*>(&x);

  Double_t lo=real->getMin();
  Double_t hi=real->getMax();

  _loArr[0]=lo;
  _hiArr[0]=hi;
  _binWidths[0]=(hi-lo)/(_nPoints-1);

  _mirrorLeft[0] =mirror==MirrorLeft || mirror==MirrorBoth || mirror==MirrorLeftAsymRight;
  _mirrorRight[0]=mirror==MirrorRight || mirror==MirrorBoth || mirror==MirrorAsymLeftRight;
  _asymLeft[0]   =mirror==MirrorAsymLeft || mirror==MirrorAsymLeftRight || mirror==MirrorAsymBoth;
  _asymRight[0]  =mirror==MirrorAsymRight || mirror==MirrorLeftAsymRight || mirror==MirrorAsymBoth;

  _lookupTable=new Double_t[_nPointsAll];
  
  // form the lookup table
  LoadDataSet(data);
  
}

//_____________________________________________________________________________
RooParamKeysPdfDev::RooParamKeysPdfDev(const char *name, const char *title,
				       const RooArgSet& xSet, const RooArgList& deltaxList, 
				       const std::vector<double>& centralValues, const RooArgList& multiplicativeShiftList, 
				       RooDataSet& data,
				       const std::vector<int>& mirror, const std::vector<Double_t>& rho,
				       Int_t nPoints) :
  RooAbsPdf(name,title),
  _nDim(rho.size()),
  _xSet("xSet","Dependent List", this),
  _dataHist(0),
  _deltaxList("deltaxList","Dependent List", this),
  _centralValues(new double[_nDim]),
  _multiplicativeShiftList("multiplicativeShiftList","Dependent",this),  
  _nEvents(0),
  _dataPts(0),
  _dataWgts(0),
  _weights(0),
  _nPoints(nPoints),
  _nPointsAll(TMath::Power(_nPoints+1,_nDim)),
  _mirrorLeft(new Bool_t[_nDim]),
  _mirrorRight(new Bool_t[_nDim]),
  _asymLeft(new Bool_t[_nDim]),
  _asymRight(new Bool_t[_nDim]),
  _varNames(new std::string[_nDim]),
  _loArr(new Double_t[_nDim]),
  _hiArr(new Double_t[_nDim]),
  _binWidths(new Double_t[_nDim]),
  _rho(new Double_t[_nDim]),
  _rotMat(new TMatrixD(_nDim,_nDim))
{
  //cout<<"2-D c'tor"<<endl;

  _xSet.add(xSet);
  _xSetIter = _xSet.createIterator();

  _histObsList.addClone(_xSet);
  _histObsIter = _histObsList.createIterator();

  RooAbsArg* x;
  for (Int_t idim=0; (x = (RooAbsArg*)_xSetIter->Next()); ++idim) {

    //_xSet.add(*x);    
    _varNames[idim]=std::string(x->GetName());
    RooRealVar* var=dynamic_cast<RooRealVar*>(x);

    _loArr[idim] = var->getMin();
    _hiArr[idim] = var->getMax();
    _binWidths[idim] = (_hiArr[idim]-_loArr[idim])/(_nPoints-1);

    RooAbsArg* delta = (RooAbsArg*)deltaxList.at(idim);
    _deltaxList.add(*delta);

    RooAbsArg* shift = (RooAbsArg*)multiplicativeShiftList.at(idim);
    _multiplicativeShiftList.add(*shift);

    Mirror m=(Mirror)mirror[idim];
    _mirrorLeft[idim]  = m==MirrorLeft      || m==MirrorBoth          || m==MirrorLeftAsymRight;
    _mirrorRight[idim] = m==MirrorRight     || m==MirrorBoth          || m==MirrorAsymLeftRight;
    _asymLeft[idim]    = m==MirrorAsymLeft  || m==MirrorAsymLeftRight || m==MirrorAsymBoth;
    _asymRight[idim]   = m==MirrorAsymRight || m==MirrorLeftAsymRight || m==MirrorAsymBoth;
    
    _centralValues[idim]=centralValues[idim];

    _rho[idim]=rho[idim];
  }
  
  // form the lookup table
  _lookupTable=new Double_t[_nPointsAll];
      
  LoadDataSet(data);
  
}



//_____________________________________________________________________________
RooParamKeysPdfDev::RooParamKeysPdfDev(const RooParamKeysPdfDev& other, const char* name):
  RooAbsPdf(other,name), 
  _nDim(other._nDim),
  _xSet("xSet",this,other._xSet), 
  _dataHist(other._dataHist),
  _deltaxList("deltaxList",this,other._deltaxList),
  _centralValues(new double[_nDim]),
  _multiplicativeShiftList("multiplicativeShiftList",this,other._multiplicativeShiftList),
  _nEvents(other._nEvents),
  _dataPts(0), _dataWgts(0), _weights(0), _sumWgt(0.),
  _nPoints(other._nPoints),
  _nPointsAll( other._nPointsAll ),
  _mirrorLeft(new Bool_t[_nDim]),
  _mirrorRight(new Bool_t[_nDim]),
  _asymLeft(new Bool_t[_nDim]),
  _asymRight(new Bool_t[_nDim]),
  _varNames(new std::string[_nDim]),
  _loArr(new Double_t[_nDim]),
  _hiArr(new Double_t[_nDim]),
  _binWidths(new Double_t[_nDim]),
  _rho(new Double_t[_nDim]),
  _rotMat(new TMatrixD(*other._rotMat))
  //_dhist(other._dhist), 
  //_hist(0)
{

  //cout<<"copy c'tor"<<endl;

  //cout<<"RooParamKeysPdfDev copy c'tor with xSet="<<_xSet<<endl;
  //_xSet.Print("V");

  _xSetIter = _xSet.createIterator();

  _histObsList.addClone(other._histObsList); 
  _histObsIter = _histObsList.createIterator();

  RooAbsArg* arg;
  for (int i=0; (arg = (RooAbsArg*)_xSetIter->Next()); i++) {

    _centralValues[i]=other._centralValues[i];
    _loArr[i]        =other._loArr[i];
    _hiArr[i]        =other._hiArr[i];
    _binWidths[i]    =other._binWidths[i];
    _mirrorLeft[i]   =other._mirrorLeft[i]; 
    _mirrorRight[i]  =other._mirrorRight[i];
    _asymLeft[i]     =other._asymLeft[i]; 
    _asymRight[i]    =other._asymRight[i];
    _rho[i]          =other._rho[i];
    _varNames[i]     = arg->GetName();
  }

  // copy over data and weights... not necessary, commented out for speed ...
  // _dataPts = new Double_t*[_nEvents];
  // _dataWgts = new Double_t[_nEvents];
  // _weights = new Double_t*[_nEvents];  
  // for (Int_t i= 0; i<_nEvents; i++) {
  //   _dataPts[i]=new Double_t[_nDim];
  //   _dataWgts[i]=other._dataWgts[i];
  //   _weights[i]=new Double_t[_nDim];
  //   for (Int_t j=0;j<_nDim;j++) {
  //     _dataPts[i][j]= other._dataPts[i][j];
  //     _weights[i][j]= other._weights[i][j];
  //   }
  // }
  

  _lookupTable=new Double_t[_nPointsAll];
  // copy over the lookup table
  for (Int_t i= 0; i<_nPointsAll; i++) {
    _lookupTable[i]= other._lookupTable[i];  
  }

}


//_____________________________________________________________________________
RooParamKeysPdfDev::~RooParamKeysPdfDev() {

  delete [] _dataPts;
  delete [] _dataWgts;
  delete [] _weights;
  delete [] _varNames;
  delete [] _loArr;
  delete [] _hiArr;
  delete [] _binWidths;

  delete [] _lookupTable;

  delete [] _mirrorLeft;
  delete [] _mirrorRight;
  delete [] _asymLeft;
  delete [] _asymRight;

  delete [] _rho;

  delete _rotMat;

  delete _histObsIter ;
  delete _xSetIter;

  //delete _dataHist;
}


void
//_____________________________________________________________________________
RooParamKeysPdfDev::LoadDataSet( RooDataSet& data) {

  //
  // Calculates _nEvents, _weights, _dataPts and _dataWgts for use by evaluateFull method
  // 

  delete [] _dataPts;
  delete [] _dataWgts;
  delete [] _weights;


  // get nEvents, extend to either side if any dimension has mirroring
  _nEvents = (Int_t)data.numEntries();

  vector<Long_t> mirrors; // the configuration is encoded in a Long_t, with each digit being a 0, 1, or 2 (NoMirror, MirrorLeft, MirrorRight)

  if (_nDim==1) {

    if (_mirrorLeft[0])  mirrors.push_back((Long_t)MirrorLeft);
    if (_mirrorRight[0]) mirrors.push_back((Long_t)MirrorRight);
  }
  else if (_nDim==2) {
    
    // first add mirror points for 2nd dimension with 1st dimension not mirrored
    if (_mirrorLeft[1])  mirrors.push_back( (Long_t) 10*MirrorLeft );  // (0,-)
    if (_mirrorRight[1]) mirrors.push_back( (Long_t) 10*MirrorRight ); // (0,+)
    
    // now add mirror left for 1st dimension, and all possibilities for 2nd dimension
    if (_mirrorLeft[0]) {
      mirrors.push_back((Long_t) MirrorLeft);                                         // (-,0)
      if (_mirrorLeft[1])  mirrors.push_back( (Long_t) MirrorLeft + 10*MirrorLeft );  // (-,-)
      if (_mirrorRight[1]) mirrors.push_back( (Long_t) MirrorLeft + 10*MirrorRight ); // (-,+)    
    }
    
    // now add mirror right for 1st dimension, and all possibilities for 2nd dimension
    if (_mirrorRight[0]) {
      mirrors.push_back((Long_t) MirrorRight);                                         // (+,0)
      if (_mirrorLeft[1])  mirrors.push_back( (Long_t) MirrorRight + 10*MirrorLeft );  // (+,-)
      if (_mirrorRight[1]) mirrors.push_back( (Long_t) MirrorRight + 10*MirrorRight ); // (+,+)    
    }      
  
  }
  else if (mirrors.size()>0) {
    std::cerr<<"not ready to do more than 2-dimensions yet with mirroring"<<std::endl;
    exit(3);
  }

  _nEvents += (Int_t) mirrors.size() * data.numEntries();

  // make new arrays for data and weights to fill
  _dataPts  = new Double_t*[_nEvents];
  _dataWgts = new Double_t[_nEvents];
  _weights  = new Double_t*[_nEvents];
  _sumWgt=0.;

  for (int ievt=0; ievt<_nEvents; ievt++) {
    _dataPts[ievt]  = new Double_t[_nDim];
    _weights[ievt]  = new Double_t[_nDim];
  }

  Double_t* x0=new Double_t[_nDim];
  Double_t* x1=new Double_t[_nDim];
  Double_t* x2=new Double_t[_nDim];

  const RooArgSet* values=data.get();
  vector<RooRealVar*> dVars(_nDim);
  for (Int_t idim=0;idim<_nDim;idim++) { 
    dVars[idim]=(RooRealVar*)values->find(_varNames[idim].c_str());
    x0[idim]=x1[idim]=x2[idim]=0.;
  }

  // matrix for ??
  TMatrixD mat(_nDim,_nDim);
  mat.Zero();
    
  Int_t i, idata=0;
  for (i=0; i<data.numEntries(); i++) {

    data.get(i);

    // event weights
    double weight=data.weight();
    _sumWgt += weight;

    _dataWgts[idata] = weight;

    // get dataPts and dataWgts, fill intermediate values for determining sigma
    for (Int_t idim=0;idim<_nDim;idim++) {
      
      double val_i=dVars[idim]->getVal();

      for (int jdim=0;jdim<_nDim;jdim++) {

	double val_j=dVars[jdim]->getVal();

	mat(idim,jdim) += val_i * val_j * weight;
      }
      
      _dataPts[idata][idim] = val_i;

      x0[idim] += weight ; 
      x1[idim] += weight * val_i; 
      x2[idim] += weight * val_i * val_i;
    }

    idata++;
  
    for (int imirror=0; imirror<(Int_t)mirrors.size(); imirror++) {
      
      Long_t val=mirrors[imirror];
      for (Int_t idim=0; idim<_nDim; idim++) {
	Double_t xval=dVars[idim]->getVal();
	Int_t dimMirror=val%10;
	if (dimMirror==0) _dataPts[idata][idim] = xval;
	else if (dimMirror==MirrorLeft) xval = 2.*_loArr[idim] - xval;
	else if (dimMirror==MirrorRight) xval = 2.*_hiArr[idim] - xval;
	val -= dimMirror;
	val /= 10;
	
	_dataPts[idata][idim] = xval;
	
	x0[idim] += weight;
	x1[idim] += weight * dVars[idim]->getVal(); 
	x2[idim] += weight * dVars[idim]->getVal() * dVars[idim]->getVal();

	for (Int_t jdim=0; jdim<_nDim; jdim++) {
	  mat(idim,jdim) += dVars[idim]->getVal() * dVars[jdim]->getVal() * weight;
	}
      }
      
      _dataWgts[idata] = data.weight();

      _sumWgt+=weight;

      idata++;      
    }
  }

  //
  // calculate _weights
  //

  Double_t *meanv=new Double_t[_nDim];
  Double_t *sigmav=new Double_t[_nDim];
  for (Int_t idim=0; idim<_nDim; idim++) {
    meanv[idim]  = x1[idim]/x0[idim];
    sigmav[idim] = sqrt(x2[idim]/x0[idim] - meanv[idim]*meanv[idim]);
  }

  TMatrixDSym covMatRho(_nDim); // covariance matrix times rho parameters
  for (Int_t j=0; j<_nDim; j++) {
    for (Int_t k=0; k<_nDim; k++) { 
      //covMat(j,k)    = mat(j,k)/x0[j] - meanv[j]*meanv[k];
      covMatRho(j,k) = (mat(j,k)/x0[j] - meanv[j]*meanv[k]) * _rho[j] * _rho[k];
    }
  }


  // find decorrelation matrix and eigenvalues (R)
  TMatrixDSymEigen evCalculatorRho(covMatRho);

  _rotMat->Zero();
  TMatrixD rotMat= evCalculatorRho.GetEigenVectors();
  *_rotMat = rotMat.T(); // transpose

  const TVectorD& sigmaRotVec = evCalculatorRho.GetEigenValues();
  
  Double_t sigmav1[_nDim];
  for (Int_t idim=0; idim<_nDim; idim++) {
    sigmav1[idim]=sqrt(sigmaRotVec[idim]);
    _rho[idim]=1.;
  }
  
  Double_t sigmaAvgR=1.;
  for (Int_t j=0; j<_nDim; j++) { sigmaAvgR *= sigmav1[j]; }
  sigmaAvgR = TMath::Power(sigmaAvgR, 1./_nDim) ;

  // calculate variables
  Double_t power=1./Double_t(_nDim+4);
  Double_t h[_nDim];
  Double_t hmin[_nDim];
  Double_t norm[_nDim];
  Double_t weights[_nDim];
  for (Int_t idim=0; idim<_nDim; idim++) {
    h[idim]   =
      TMath::Power( Double_t(4)/Double_t(_nDim+2)/Double_t(_nEvents), power) * _rho[idim];
    hmin[idim]    = h[idim]*sigmav1[idim]*sqrt(2.)/10.0;
    norm[idim]    = h[idim]*sqrt(sigmav1[idim])/(2.0*sqrt(3.0));
    weights[idim] = h[idim]*sigmav1[idim];
  }


  for(Int_t j=0;j<_nEvents;++j) {

    Double_t gaus=fgaus(_dataPts[j],weights);
    Double_t denominator = (_nDim==1) ? sqrt(gaus) : TMath::Power(gaus,1./Double_t(2*_nDim));
    for (Int_t idim=0; idim<_nDim; idim++) {
      _weights[j][idim]=norm[idim]/denominator;
      if (_weights[j][idim]<hmin[idim]) {
	_weights[j][idim]=hmin[idim];
      }
    }
  }


  delete x0;
  delete x1;
  delete x2;

 
  Double_t xVar[_nDim];
  /*
  Int_t indices[_nDim];
  for (Int_t idim=0;idim<_nDim;idim++) {
    indices[idim]=0;
    xVar[idim]=_loArr[idim];
  }
    
  for (int ipt=0;ipt<_nPointsAll;++ipt) {
    
    //for (int jdim=_nDim-1; jdim>=0; jdim--) {
    for (int jdim=0; jdim<_nDim; jdim++) {
      bool doBreak=true;
      xVar[jdim]=_loArr[jdim]+(Double_t(indices[jdim])*_binWidths[jdim]);
      indices[jdim]=indices[jdim]+1;
      if (indices[jdim]==_nPoints+1) { 
	indices[jdim]=0;
	doBreak=false; 
      }
      if (doBreak) {
	xVar[jdim-1]=_loArr[jdim-1]+(Double_t(indices[jdim-1])*_binWidths[jdim-1]);	
	break;
      }
    }
    
    _lookupTable[ipt]=evaluateFull( xVar );
    if (ipt%10000==0) //cout<<"Processing "<<ipt<<"/"<<_nPointsAll<<endl;
  }  
  */

  Char_t histname[128]; snprintf(histname,128,"%s_hist",GetName());
  TH1* hist(0);

  if (_nDim==1) {
    hist = new TH1F(histname,histname,_nPoints,_loArr[0],_hiArr[0]);
    //Double_t xVar[1];
    for (int ipt=0; ipt<_nPointsAll; ipt++) {
      xVar[0]=_loArr[0]+(Double_t(ipt)*_binWidths[0]);
      double val=evaluateFull(xVar)/_binWidths[0];
      ((TH1F*)hist)->SetBinContent(ipt+1,val);
      if (ipt%1000==0) cout<<"Processing "<<ipt<<"/"<<_nPointsAll<<endl;
    }
  }
  else if (_nDim==2) {
    
    hist = new TH2F(histname,histname,_nPoints,_loArr[0],_hiArr[0],_nPoints,_loArr[1],_hiArr[1]);

    int ipt0=0;
    for (int ipt=0; ipt<_nPoints+1; ipt++) {
      xVar[0]=_loArr[0]+(Double_t(ipt)*_binWidths[0]);
      for (int jpt=0; jpt<_nPoints+1; jpt++) {
  	xVar[1]=_loArr[1]+(Double_t(jpt)*_binWidths[1]);
  	((TH2F*)hist)->SetBinContent(ipt+1,jpt+1,evaluateFull(xVar)/_binWidths[0]/_binWidths[1]);
	if (ipt0%10000==0) cout<<"Processing "<<ipt0<<"/"<<_nPointsAll<<endl;
	ipt0++;
      }
    }
  }
  
  Char_t dhistname[128]; snprintf(dhistname,128,"%s_dhist",GetName());
  _dataHist = new RooDataHist(dhistname,dhistname,_xSet,hist,1.);

}



//_____________________________________________________________________________
Double_t RooParamKeysPdfDev::evaluate() const {
  
  /*
  double val = 0.0;

  bool forcePositive = false; // this is just to suppress warning for values outside of range.

  Double_t xVar[_nDim];
  Double_t dx[_nDim];
  Int_t indices[_nDim];
  Int_t igridpt0=0;
  Int_t dim=1;

  _xSetIter->Reset();
  RooRealVar* x;
  for (int idim=0; (x = (RooRealVar*)_xSetIter->Next()); ++idim) {

    indices[idim]=0;

    RooAbsReal* delta=(RooAbsReal*)_deltaxList.at(idim);
    RooAbsReal* multiplicativeShift=(RooAbsReal*)_multiplicativeShiftList.at(idim);

    double deltax = delta->getVal()-_centralValues[idim];
    if( multiplicativeShift ) 
      deltax += _centralValues[idim]*(multiplicativeShift->getVal()-1.0);
    
    xVar[idim]=x->getVal()-deltax;

    Int_t ifloor = (Int_t)floor((x->getVal()-deltax-_loArr[idim])/_binWidths[idim]);
    if (ifloor<0)          { ifloor=0;          forcePositive = true; }
    if (ifloor>_nPoints-1) { ifloor=_nPoints-1; forcePositive = true; }

    dx[idim] = (Double_t(x->getVal()-deltax)-(_loArr[idim]+ifloor*_binWidths[idim]))/_binWidths[idim];

    igridpt0 += ifloor*dim;

    dim *= _nPoints+1;// dim=TMath::Power(_nPoints+1,idim)
    
  } 

  // set up initial (2 ^ _nDim) grid from which linear interpolation will be done
  Int_t ngridpoints=1<<_nDim; // 2^_nDim
  Double_t* gridvalues=new Double_t[ngridpoints];
  for (Int_t igridpt=0;igridpt<ngridpoints;igridpt++) {
      
    Int_t ipt=igridpt0;
    Int_t dim1=1;
    for (Int_t idim=0;idim<_nDim;idim++) {
      ipt+=((igridpt>>idim)%2)*dim1;     
      dim1*=_nPoints+1; // dim=TMath::Power(_nPoints+1,idim)
    } 
    gridvalues[igridpt]=_lookupTable[ipt];
  }


  // now iterate, interpolate and reduce grid by a power of 2 each time until only 1 value remaining
  for (Int_t iter=_nDim-1;iter>=0;iter--) {

    int new_ngridpoints=1<<iter; // 2^iter

    Double_t *new_gridvalues=new Double_t[new_ngridpoints];
    int igridpt=0;    
    for (Int_t j=0;j<iter+1;j++) {  
      Double_t upper=gridvalues[igridpt+new_ngridpoints];
      Double_t lower=gridvalues[igridpt];
      new_gridvalues[igridpt] = (upper-lower)*dx[iter]+lower;
      igridpt++;
    }
    delete [] gridvalues;
    gridvalues=new_gridvalues;
  }

  // last dimension
  val = gridvalues[0];
  delete [] gridvalues;
    
  if( forcePositive && val<0.0 ) val = 0.0; 
  */
  // Calculate current value of object

  /*
  // Transfer values from   
  if (_xSet.getSize()>0) {
    _histObsIter->Reset() ;
    _xSetIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
      parg = (RooAbsArg*)_xSetIter->Next() ;
      if (harg != parg) {
	parg->syncCache() ;
	harg->copyCache(parg,kTRUE) ;
	if (!harg->inRange(0)) {
	  return 0 ;
	}
      }
    }
  }
  */

  //cout<<"RooParamKeysPdfDev::evaluate()"<<endl;

  _xSetIter->Reset();
  _histObsIter->Reset();

  const RooRealVar* parg;
  RooRealVar* harg;

  int idim=0;
  while ((parg=(RooRealVar*)_xSetIter->Next())) {

    harg=(RooRealVar*)_histObsIter->Next();

    RooAbsReal* delta=(RooAbsReal*)_deltaxList.at(idim);
    RooAbsReal* multiplicativeShift=(RooAbsReal*)_multiplicativeShiftList.at(idim);
    
    double deltax = delta->getVal()-_centralValues[idim];
    if( multiplicativeShift ) 
      deltax += _centralValues[idim]*(multiplicativeShift->getVal()-1.0);

    harg->setVal(parg->getVal()-deltax);

    idim++;
    //cout<<"harg: "<<harg->getVal();
  }

  //cout<<"calling dh->weight"<<endl;

  Double_t ret =  _dataHist->weight(_histObsList,3,kFALSE,kFALSE) ;  

  //cout<<", ret="<<ret<<endl;

  if (ret<0) {
    ret=0 ;
  }  
  return ret;

}


//_____________________________________________________________________________
Double_t RooParamKeysPdfDev::evaluateFull( Double_t* xVar ) const {
  
  Double_t y=0;

  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));

  bool asymLeft=false;
  bool asymRight=false;
  for (Int_t idim=0; idim<_nDim; idim++) {
    if (_asymLeft[idim]) asymLeft=true;
    if (_asymRight[idim]) asymRight=true;
  }

  for (Int_t ievt=0;ievt<_nEvents;++ievt) {

    TVectorD dx(_nDim);
    for (int idim=0; idim<_nDim; idim++) 
      dx[idim] = xVar[idim]-_dataPts[ievt][idim];

    dx *= *_rotMat;

    Double_t rfull=1.;
    for (int idim=0; idim<_nDim; idim++) {
      Double_t chi=dx[idim]/_weights[ievt][idim];
      rfull *= exp(-0.5*chi*chi)/(_weights[ievt][idim]*sqrt2pi);      
    }
    y += _dataWgts[ievt]*rfull;

    // if mirroring the distribution across either edge of
    // the range ("Boundary Kernels"), pick up the additional
    // contributions

    // asymmetric right
    if (asymRight) {
      rfull=1.;
      for (Int_t idim=0; idim<_nDim; idim++) {
	if (_asymLeft[idim]) {	
	  Double_t chi=(xVar[idim]-(2*_loArr[idim]-_dataPts[ievt][idim]))/_weights[ievt][idim];
	  rfull *= exp(-0.5*chi*chi)/(_weights[ievt][idim]*sqrt2pi);
	}
	y-=_dataWgts[ievt]*rfull;
      }
    }
    
    // asymmetric right
    if (asymLeft) {
      rfull=1.;
      for (Int_t idim=0; idim<_nDim; idim++) {
	if (_asymRight[idim]) {
	  Double_t chi=(xVar[idim]-(2*_hiArr[idim]-_dataPts[ievt][idim]))/_weights[ievt][idim];
	  rfull *= exp(-0.5*chi*chi)/(_weights[ievt][idim]*sqrt2pi);
	}
	y-=_dataWgts[ievt]*rfull;
      }
    }
    
  }
  
  double val=y/_sumWgt;
  return val;

}
  

//_____________________________________________________________________________
Double_t RooParamKeysPdfDev::fgaus(Double_t* x,Double_t* sigmav) const {
  
  Double_t c=Double_t(1)/Double_t(2);

  Double_t y=0;
  for (Int_t ievt=0;ievt<_nEvents;++ievt) {

    TVectorD dx(_nDim);
    for (int idim=0; idim<_nDim; idim++) {      
      dx[idim]=x[idim]-_dataPts[ievt][idim];
    }
 
    // rotate to decorrelated frame
    dx *= *_rotMat;

    Double_t rfull=1.;
    for (Int_t idim=0; idim<_nDim; idim++) {    
      Double_t r = dx[idim]/sigmav[idim];
      rfull *= exp(-c*r*r);
    }
    y+=rfull;
  }

  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));
  for (Int_t idim=0; idim<_nDim; idim++) y/=(sqrt2pi*sigmav[idim]);

  return y/_nEvents;
}

namespace {
    bool fullRange(const RooAbsArg& x ,const char* range)  {
      if (range == 0 || strlen(range) == 0 ) return true;
      const RooAbsRealLValue *_x = dynamic_cast<const RooAbsRealLValue*>(&x);
      if (!_x) return false;
      return ( _x->getMin(range) == _x->getMin() && _x->getMax(range) == _x->getMax() ) ; 
    }
}

//________________________________________________________________________
RooAbsReal* RooParamKeysPdfDev::createIntegral(const RooArgSet& iset, const RooArgSet* nset, const RooNumIntConfig* cfg, const char* rangeName) const
{
  return RooAbsPdf::createIntegral(iset,nset,cfg,rangeName);

  //cout<<"in createIntegral"<<endl;

  if (_nDim==1) { cout<<"returning RooAbsPdf::createIntegral"<<endl; return RooAbsPdf::createIntegral(iset,nset,cfg,rangeName); }

  if (iset.getSize()==2) { 

    Double_t sum=_dataHist->sum(kTRUE,kFALSE) ;

    string varname=Form("%s_allInt",GetName());

    return new RooRealVar(varname.c_str(),varname.c_str(),sum);

  }  

  // create 1-D hist
  //cout<<"iset="<<iset<<endl;
  //cout<<"nset="<<nset<<endl;

  RooArgSet newXSet;
  _xSetIter->Reset();
  _histObsIter->Reset();
  RooRealVar* parg;
  RooRealVar* harg;
  Int_t idim=0;
  while ((parg=(RooRealVar*)_xSetIter->Next())) {
    //cout<<"parg="<<parg<<endl;
    harg=(RooRealVar*)_histObsIter->Next();
    if (parg==iset.first()) {
      //cout<<"parg=iset.first"<<endl;
      break;
    }
  }
 
  //cout<<"parg: "<<parg<<endl; parg->Print();
  //cout<<"harg: "<<harg<<endl; harg->Print();

  string dhname=Form("%s_%s",GetName(),parg->GetName());
  
  Double_t lo=parg->getMin();
  Double_t hi=parg->getMax();
  Double_t binwidth=(hi-lo)/(_nPoints-1);
  
  TH1F* hist=new TH1F((dhname+"_hist").c_str(),(dhname+"_hist").c_str(),_nPoints,parg->getMin(),parg->getMax());
  for (Int_t ipt=0;ipt<_nPoints;ipt++) {
    double xval=parg->getMin()+binwidth*(Double_t)ipt;    
    parg->setVal(xval);
    double yval=_dataHist->sum(iset,_histObsList,kTRUE,kTRUE);
    hist->SetBinContent(ipt+1,yval);

  }
  //cout<<"creating newDataHist"<<endl;
  RooDataHist* newDataHist = new RooDataHist((dhname+"_dhist").c_str(),(dhname+"_dhist").c_str(),RooArgList(*harg),hist,1.0);

  //cout<<"creating histPdf"<<endl;
  RooHistPdf* histPdf = new RooHistPdf((dhname+"_dhist").c_str(),(dhname+"_dhist").c_str(),RooArgSet(*harg),*newDataHist,3);

  histPdf->Print();

  return histPdf;
}

//_____________________________________________________________________________
Int_t RooParamKeysPdfDev::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{

  //cout<<"RooParamKeysPdfDev::getAnalyticalIntegral, allVars="<<allVars<<", analVars="<<analVars<<endl;

  return 0;

   // Determine integration scenario. If no interpolation is used,
  // RooHistPdf can perform all integrals over its dependents
  // analytically via partial or complete summation of the input
  // histogram. If interpolation is used on the integral over
  // all histogram observables is supported


  // First make list of pdf observables to histogram observables
  // and select only those for which the integral is over the full range
  RooArgList hobsl(_histObsList),pobsl(_xSet) ;
  RooArgSet allVarsHist ;
  TIterator* iter = allVars.createIterator() ;
  RooAbsArg* pdfobs ;
  while((pdfobs=(RooAbsArg*)iter->Next())) {
    Int_t idx = pobsl.index(pdfobs) ;
    if (idx>=0) {
      RooAbsArg* hobs = hobsl.at(idx) ;
      if (hobs && fullRange( *hobs, rangeName ) ) {
	allVarsHist.add(*hobs) ;
      }
    }
  }
  delete iter ;
  //cout<<"allVarsHist: "<<allVarsHist<<endl;

  // Simplest scenario, integrate over all dependents
  RooAbsCollection *allVarsCommon = allVarsHist.selectCommon(_histObsList) ;  
  Bool_t intAllObs = (allVarsCommon->getSize()==_histObsList.getSize()) ;
  delete allVarsCommon ;
  if (intAllObs) {
    analVars.add(allVars) ;
    cout<<"RooParamKeysPdfDev::getAnalyticalIntegral, returning 1000, analVars="<<analVars<<endl;
    return 1000 ;
  }

  // Disable partial analytical integrals if interpolation is used
//   if (_intOrder>0) {
//     return 0 ;
//   }

  // Find subset of _histObsList that integration is requested over
  //cout<<"histObsList="<<_histObsList<<endl;

  RooArgSet* allVarsSel = (RooArgSet*) allVarsHist.selectCommon(_histObsList) ;
  if (allVarsSel->getSize()==0) {
    delete allVarsSel ;
    cout<<"RooParamKeysPdfDev::getAnalyticalIntegral returning 0"<<endl;
    return 0 ;
  }

  // Partial integration scenarios.
  // Build unique code from bit mask of integrated variables in depList
  Int_t code(0),n(0) ;
  iter = _histObsList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (allVarsHist.find(arg->GetName())) {
      code |= (1<<n) ;
      analVars.add(*pobsl.at(n)) ;
    }
    n++ ;
  }
  delete iter ;
  cout<<"RooParamKeysPdfDev::getAalyticalIntegral returning code="<<code<<endl; 
  return code ;
 
}

//_____________________________________________________________________________
Double_t RooParamKeysPdfDev::analyticalIntegral(Int_t code, const char* ) const 
{
   // Return integral identified by 'code'. The actual integration
  // is deferred to RooDataHist::sum() which implements partial
  // or complete summation over the histograms contents

  //cout<<"RooParamKeysPdfDev::analyticalIntegral"<<endl;

  _xSetIter->Reset();
  _histObsIter->Reset();
  
  const RooRealVar* parg;
  RooRealVar* harg;
  
  int idim=0;
  while ((parg=(RooRealVar*)_xSetIter->Next())) {
    
    harg=(RooRealVar*)_histObsIter->Next();
    
    RooAbsReal* delta=(RooAbsReal*)_deltaxList.at(idim);
    RooAbsReal* multiplicativeShift=(RooAbsReal*)_multiplicativeShiftList.at(idim);
    
    double deltax = delta->getVal()-_centralValues[idim];
    if( multiplicativeShift ) 
      deltax += _centralValues[idim]*(multiplicativeShift->getVal()-1.0);
    
    harg->setVal(parg->getVal()-deltax);
    
    idim++;
  }

  // Simplest scenario, integration over all dependents
  if (code==1000) {    
    return _dataHist->sum(kTRUE,kFALSE) ;
  }

  // Partial integration scenario, retrieve set of variables, calculate partial sum
  RooArgSet intSet ;
  TIterator* iter = _histObsList.createIterator() ;
  RooAbsArg* arg ;
  Int_t n(0) ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (code & (1<<n)) {
      intSet.add(*arg) ;
    }
    n++ ;
  }
  delete iter ;
  //cout<<"intSet (before): "<<endl;
  //intSet.Print("V");

  /*
  // WVE must sync hist slice list values to pdf slice list
  // Transfer values from   
  if (_xSet.getSize()>0) {
    _histObsIter->Reset() ;
    _xSetIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
      parg = (RooAbsArg*)_xSetIter->Next() ;
      if (harg != parg) {
	parg->syncCache() ;
	harg->copyCache(parg,kTRUE) ;
      }
    }
  }  
  */


  //cout << "intSet = " << intSet << endl ;
  //cout << "slice position = " << endl ;
  //_histObsList.Print("v") ;
  
  Double_t ret =  _dataHist->sum(intSet,_histObsList,kTRUE,kTRUE) ;

//    cout << "intSet = " << intSet << endl ;
//    cout << "slice position = " << endl ;
//    _histObsList.Print("v") ;
//    cout << "RooHistPdf::ai(" << GetName() << ") code = " << code << " ret = " << ret << endl ;

  return ret ;
}






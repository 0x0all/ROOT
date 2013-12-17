// @(#)root/meta:$Id$
// Author: Fons Rademakers   08/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBaseClass.h"
#include "TClass.h"
#include "TInterpreter.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Each class (see TClass) has a linked list of its base class(es).    //
//  This class describes one single base class.                         //
//  The base class info is obtained via the CINT api.                   //
//     see class TCling.                                                 //
//                                                                      //
//  The base class information is used a.o. in to find all inherited    //
//  methods.                                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


ClassImp(TBaseClass)

//______________________________________________________________________________
TBaseClass::TBaseClass(BaseClassInfo_t *info, TClass *cl) : TDictionary()
{
   // Default TBaseClass ctor. TBaseClasses are constructed in TClass
   // via a call to TCling::CreateListOfBaseClasses().

   fInfo     = info;
   fClass    = cl;
   fClassPtr = 0;
   if (fInfo) SetName(gCling->BaseClassInfo_FullName(fInfo));
}

//______________________________________________________________________________
TBaseClass::~TBaseClass()
{
   // TBaseClass dtor deletes adopted CINT BaseClassInfo object.

   gCling->BaseClassInfo_Delete(fInfo);
}

//______________________________________________________________________________
void TBaseClass::Browse(TBrowser *b)
{
   // Called by the browser, to browse a baseclass.

   TClass *c = GetClassPointer();
   if (c) c->Browse(b);
}

//______________________________________________________________________________
TClass *TBaseClass::GetClassPointer(Bool_t load)
{
   // Get pointer to the base class TClass.

   if (!fClassPtr) {
      if (fInfo) fClassPtr = TClass::GetClass(gCling->BaseClassInfo_ClassInfo(fInfo),load);
      else fClassPtr = TClass::GetClass(fName, load);
   }
   return fClassPtr;
}

//______________________________________________________________________________
Int_t TBaseClass::GetDelta() const
{
   // Get offset from "this" to part of base class.

   return (Int_t)gCling->BaseClassInfo_Offset(fInfo);
}

//______________________________________________________________________________
const char *TBaseClass::GetTitle() const
{
   // Get base class description (comment).

   TClass *c = ((TBaseClass *)this)->GetClassPointer();
   return c ? c->GetTitle() : "";
}

//______________________________________________________________________________
ROOT::ESTLType TBaseClass::IsSTLContainer()
{
   // Return which type (if any) of STL container the data member is.

   if (!fInfo) return ROOT::kNotSTL;
   const char *type = gCling->BaseClassInfo_TmpltName(fInfo);
   if (!type) return ROOT::kNotSTL;

   if (!strcmp(type, "vector"))   return ROOT::kSTLvector;
   if (!strcmp(type, "list"))     return ROOT::kSTLlist;
   if (!strcmp(type, "deque"))    return ROOT::kSTLdeque;
   if (!strcmp(type, "map"))      return ROOT::kSTLmap;
   if (!strcmp(type, "multimap")) return ROOT::kSTLmultimap;
   if (!strcmp(type, "set"))      return ROOT::kSTLset;
   if (!strcmp(type, "multiset")) return ROOT::kSTLmultiset;
   return ROOT::kNotSTL;
}

//______________________________________________________________________________
Long_t TBaseClass::Property() const
{
   // Get property description word. For meaning of bits see EProperty.

   return gCling->BaseClassInfo_Property(fInfo);
}

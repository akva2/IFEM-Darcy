// $Id$
//==============================================================================
//!
//! \file MixedDarcy.C
//!
//! \date oct 20 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mixed Darcy flow problems.
//!
//==============================================================================

#include "MixedDarcy.h"

#include "AnaSol.h"
#include "BlockElmMats.h"
#include "ElmNorm.h"
#include "FiniteElement.h"
#include "Function.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "TimeDomain.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <cmath>
#include <ext/alloc_traits.h>
#include <iostream>
#include <memory>
#include <vector>


MixedDarcy::MixedDarcy (unsigned short int n, int torder) :
  Darcy(n, torder)
{
  pp = 1;
  cc = 2;
  cp = 4;
  npv = 2;
  sourceC = nullptr;
}


LocalIntegral* MixedDarcy::getLocalIntegral (size_t nen, size_t, bool neumann) const
{
  BlockElmMats* result = new BlockElmMats(2, 1);

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->resize(5, 3);
  result->redim(pp, nen, 1, 1);
  result->redim(cc, nen, 1, 1);
  result->redimOffDiag(cp, 0);
  result->redimNewtonMat();

  return result;
}


bool MixedDarcy::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                          const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!this->Darcy::evalInt(elmInt, fe, time, X))
    return false;

  if (!dispersivity || !porosity)
    return false;

  const double D = (*this->dispersivity)(X);
  const double phi = (*this->porosity)(X);

  WeakOps::Laplacian(elMat.A[cc], fe, D);

  double cn = fe.N.dot(elmInt.vec[0],1,2);

  Vec3 perm = this->getPermeability(X);
  WeakOps::Laplacian(elMat.A[cp], fe, perm[0]*cn);

  if (sourceC)
    WeakOps::Source(elMat.b[cc], fe, (*sourceC)(X));

  if (bdf.getActualOrder() > 0 && elmInt.vec.size() > 1) {
    double c = 0.0;
    for (int t = 1; t <= bdf.getOrder(); t++) {
      double val = fe.N.dot(elmInt.vec[t],1,2);
      c -= val * phi * bdf[t] / time.dt;
    }
    WeakOps::Source(elMat.b[cc], fe, c);
    WeakOps::Mass(elMat.A[cc], fe, phi*bdf[0] / time.dt);
  }

  return true;
}


bool MixedDarcy::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                          const Vec3& X, const Vec3& normal) const
{
  if (!this->Darcy::evalBou(elmInt,fe,X,normal))
    return false;

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!dispersivity)
    return false;

  const double D = (*this->dispersivity)(X);

  Vec3 perm = this->getPermeability(X);

  for (size_t i = 1; i <= fe.N.size(); ++i)
    for (size_t j = 1; j <= fe.N.size(); ++j)
      for (int k = 1; k <= nsd; ++k)
        elMat.A[cp](i,j) += (fe.N(j)*perm[k-1]*fe.dNdX(j,k)*normal[k-1] - D*fe.dNdX(j,k))*fe.N(i)*fe.detJxW;

  return true;
}


bool MixedDarcy::evalSol2 (Vector& s, const Vectors& eV,
                           const FiniteElement& fe, const Vec3& X) const
{
  if (!this->Darcy::evalSol2(s,eV,fe,X))
    return false;

  // Evaluate the concentration gradient
  RealArray temp;
  if (eV.empty() || !fe.dNdX.multiply(eV.front(),temp,1.0,0.0,true,2,1,1))
  {
    std::cerr <<" *** MixedDarcy::evalSol: Invalid solution vector.\n"
              <<"     size(eV) = "<< (eV.empty() ? 0 : eV.front().size())
              <<" size(dNdX) = "<< fe.dNdX.rows() <<","<< fe.dNdX.cols()
              << std::endl;
    return false;
  }

  for (const double& v : temp)
    s.push_back(v);

  return true;
}


std::string MixedDarcy::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return "p&&c";

  if (i >= 2)
    return "";

  static const char* s[2] = {"p", "c"};
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string MixedDarcy::getField2Name (size_t i, const char* prefix) const
{
  if (nsd == 2 && i > 1)
    ++i;

  if (i >= 6) return "";

  static const char* s[6] = {"p,x","p,y","p,z",
                             "c,x","c,y","c,z"};

  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* MixedDarcy::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new MixedDarcyNorm(*const_cast<MixedDarcy*>(this),
                              asol->getScalarSecSol(0),
                              asol->getScalarSecSol(1));

  else
    return new MixedDarcyNorm(*const_cast<MixedDarcy*>(this));
}


MixedDarcyNorm::MixedDarcyNorm (MixedDarcy& p, VecFunc* a, VecFunc* c)
  : DarcyNorm(p,a), anac(c)
{
}


bool MixedDarcyNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3& X) const
{
  if (!this->DarcyNorm::evalInt(elmInt,fe,X))
    return false;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const Darcy& problem = static_cast<const Darcy&>(myProblem);
  const double D = problem.getDispersivity(X);

  // Evaluate the concentration field gradient
  Vector dCh(fe.dNdX.cols());
 fe.dNdX.multiply(pnorm.vec.front(),dCh,1.0,0.0,true,2,1,1);

  pnorm[H1_Ch] += D*dCh.dot(dCh)*fe.detJxW;

  Vector dC(fe.dNdX.cols());
  if (anac) {
    dC.fill((*anac)(X).ptr(),fe.dNdX.cols());
    pnorm[H1_C] += D*dC.dot(dC)*fe.detJxW;
    Vector error =  dC-dCh;
    pnorm[H1_E_Ch] += D*error.dot(error)*fe.detJxW;
  }

  size_t ip = this->getNoFields(1);
  size_t f = 2;
  for (const Vector& psol : pnorm.psol) {
    if (!psol.empty())
    {
      // Evaluate projected concentration field
      Vector dCr(fe.dNdX.cols());

      for (size_t j = 0; j < fe.dNdX.cols(); j++)
        dCr[j] = psol.dot(fe.N,j+fe.dNdX.cols(),nrcmp);

      // Integrate the energy norm a(c^r,c^r)
      pnorm[ip+H1_Cr] += D*dCr.dot(dCr)*fe.detJxW;
      // Integrate the estimated error in energy norm a(c^r-c^h,c^r-c^h)
      Vector error = dCr - dCh;
      pnorm[ip+H1_Cr_Ch] += D*error.dot(error)*fe.detJxW;

      if (anac)
      {
        // Integrate the error in the projected solution a(c-c^r,c-c^r)
        error = dC - dCr;
        pnorm[ip+H1_E_Cr] += D*error.dot(error)*fe.detJxW;
      }
    }
    ip += this->getNoFields(f++);
  }

  return true;
}

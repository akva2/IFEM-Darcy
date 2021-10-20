// $Id$
//==============================================================================
//!
//! \file MixedDarcy.h
//!
//! \date Oct 20 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for mixed-formulation Darcy flow problems.
//!
//==============================================================================

#ifndef _MIXED_DARCY_H_
#define _MIXED_DARCY_H_

#include "Darcy.h"


/*!
  \brief Class representing the integrand of the mixed Darcy problem.
*/

class MixedDarcy : public Darcy
{
public:
  //! \brief The constructor initializes all pointers to zero.
  explicit MixedDarcy(unsigned short int n, int torder = 0);
  //! \brief Empty destructor.
  virtual ~MixedDarcy() {}

  //! \brief Defines the concentration source function.
  void setCSource(RealFunc* s) override { sourceC = s; }

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral contribution object for given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time,
               const Vec3& X) const override;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  using IntegrandBase::evalSol2;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? 2*nsd : 2; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix) const override;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

protected:
  RealFunc* sourceC;       //!< Concentration source function
};


/*!
  \brief Class representing the integrand of mixed Darcy energy norms.
*/

class MixedDarcyNorm : public DarcyNorm
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical pressure gradient (optional)
  //! \param[in] c The analytical concentration gradient (optional)
  explicit MixedDarcyNorm(MixedDarcy& p,
                          VecFunc* a = nullptr,
                          VecFunc* c = nullptr);

  //! \brief Empty destructor.
  virtual ~MixedDarcyNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

private:
  VecFunc* anac; //!< Analytical concentration gradient
};

#endif

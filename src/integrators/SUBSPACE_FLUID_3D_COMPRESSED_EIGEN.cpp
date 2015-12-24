/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2013 Theodore Kim
 */
// SUBSPACE_FLUID_3D_COMPRESSED_EIGEN.cpp: implementation of the SUBSPACE_FLUID_3D_COMPRESSED_EIGEN class.
//
//////////////////////////////////////////////////////////////////////

#include "SUBSPACE_FLUID_3D_COMPRESSED_EIGEN.h"
#include "BIG_MATRIX.h"
#include "MATRIX_COMPRESSION_DATA.h"
#include "COMPRESSION_REWRITE.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::SUBSPACE_FLUID_3D_COMPRESSED_EIGEN(int xRes, int yRes, int zRes, const string& reducedPath, unsigned int* boundaries, bool usingIOP, bool loadNothing) :
  FLUID_3D_MIC(xRes, yRes, zRes, 4, boundaries), _reducedPath(reducedPath), _discardThreshold(-1), _usingIOP(usingIOP)
{
  _snapshotPath = string("");

  // set up the boundary conditions
  if (boundaries != NULL)
  {
    _domainBcFront  = boundaries[0];
    _domainBcBack   = boundaries[1];
    _domainBcLeft   = boundaries[2];
    _domainBcRight  = boundaries[3];
    _domainBcTop    = boundaries[4];
    _domainBcBottom = boundaries[5];
  }
  else
  {
    _domainBcFront  = 1; 
    _domainBcBack   = 1;
    _domainBcLeft   = 1;
    _domainBcRight  = 1;
    _domainBcTop    = 0;
    _domainBcBottom = 0;
  }

  if (!loadNothing)
  {
    if (_usingIOP) { 
      cout << "Calling initOutOfCoreIOP(). " << endl;
      initOutOfCoreIOP();
    }

    else { 
      cout << "Calling vanilla initOutOfCore. " << endl;
      initOutOfCore();
    }

    readAdvectionCubature();
  }
  else
  {
    // init the peeled dimensions
    _xPeeled = _xRes - 2;
    _yPeeled = _yRes - 2;
    _zPeeled = _zRes - 2;
    _slabPeeled = _xPeeled * _yPeeled;
  }

}


//////////////////////////////////////////////////////////////////////
// initialize the peeled version where there is a separate basis
// for each stage
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::initOutOfCore()
{
  // init the peeled dimensions
  _xPeeled = _xRes - 2;
  _yPeeled = _yRes - 2;
  _zPeeled = _zRes - 2;
  _slabPeeled = _xPeeled * _yPeeled;

  bool pcaBuilt = fileExists(_reducedPath + string("U.final.componentX")) &&
                  fileExists(_reducedPath + string("U.final.componentY")) &&
                  fileExists(_reducedPath + string("U.final.componentZ")) &&
                  fileExists(_reducedPath + string("U.preproject.matrix")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentX")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentY")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentZ")) &&
                  fileExists(_reducedPath + string("U.pressure.matrix"));

  // check if pre-built matrices exist
  bool filesBuilt = fileExists(_reducedPath + string("projected.A.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptof.matrix")) &&
                    fileExists(_reducedPath + string("projected.vtod.matrix")) &&
                    fileExists(_reducedPath + string("damping.peeled.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptov.matrix")) &&
                    fileExists(_reducedPath + string("inverseProduct.matrix"));

  if (!filesBuilt) {
    cout << "Error: pre-built matrices do not yet exist!" << endl;
    cout << "You should run ./buildProducts first." << endl;
    exit(1);
  }
  else
  {
    string filename;
    
    filename = _reducedPath + string("projected.ptof.matrix");
    EIGEN::read(filename, _preprojectToFinal);
    cout << "Read projected.ptof.matrix into _preprojectToFinal. " << endl;

    filename = _reducedPath + string("projected.vtod.matrix");
    EIGEN::read(filename, _reducedVelocityToDivergence);

    filename = _reducedPath + string("projected.ptov.matrix");
    EIGEN::read(filename, _reducedPressureToVelocity);

    filename = _reducedPath + string("projected.A.matrix");
    EIGEN::read(filename, _reducedA);

    filename = _reducedPath + string("damping.peeled.matrix");
    EIGEN::read(filename, _dampingMatrixReduced);
    
    filename = _reducedPath + string("inverseProduct.matrix");
    bool success = EIGEN::read(filename, _inverseProduct);

    if (!success)
    {
      // Needs everything prior to be built already
      MatrixXd inverse = _reducedA.inverse();
      _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
      filename = _reducedPath + string("inverseProduct.matrix");
      EIGEN::write(filename, _inverseProduct);
      TIMER::printTimings();
    }
  }
}

//////////////////////////////////////////////////////////////////////
// initialize the peeled version where there is a separate basis
// for each stage, using IOP
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::initOutOfCoreIOP()
{
  // init the peeled dimensions
  _xPeeled = _xRes - 2;
  _yPeeled = _yRes - 2;
  _zPeeled = _zRes - 2;
  _slabPeeled = _xPeeled * _yPeeled;

  bool pcaBuilt = fileExists(_reducedPath + string("U.final.componentX")) &&
                  fileExists(_reducedPath + string("U.final.componentY")) &&
                  fileExists(_reducedPath + string("U.final.componentZ")) &&
                  fileExists(_reducedPath + string("U.preproject.matrix")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentX")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentY")) &&
                  fileExists(_reducedPath + string("U.preadvect.componentZ")) &&
                  fileExists(_reducedPath + string("U.pressure.matrix"));

  // check if pre-built matrices exist
  bool filesBuilt = fileExists(_reducedPath + string("projected.A.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptof.matrix")) &&
                    fileExists(_reducedPath + string("projected.vtod.matrix")) &&
                    fileExists(_reducedPath + string("damping.peeled.matrix")) &&
                    fileExists(_reducedPath + string("projected.ptov.matrix")) &&
                    fileExists(_reducedPath + string("inverseProduct.matrix"));

  if (!filesBuilt) {
    cout << "Error: pre-built matrices do not yet exist!" << endl;
    cout << "You should run ./buildProducts first." << endl;
    exit(1);
  }
  else
  {
    string filename;
    
    filename = _reducedPath + string("projected.ptof.matrix");
    // EIGEN::read(filename, _preprojectToPreadvect);
    EIGEN::read(filename, _preprojectToFinal);

    filename = _reducedPath + string("projected.vtod.matrix");
    EIGEN::read(filename, _reducedVelocityToDivergence);

    filename = _reducedPath + string("projected.ptov.matrix");
    EIGEN::read(filename, _reducedPressureToVelocity);

    filename = _reducedPath + string("projected.A.matrix");
    EIGEN::read(filename, _reducedA);

    filename = _reducedPath + string("damping.peeled.matrix");
    EIGEN::read(filename, _dampingMatrixReduced);
    
    filename = _reducedPath + string("inverseProduct.matrix");
    bool success = EIGEN::read(filename, _inverseProduct);

    if (!success)
    {
      // Needs everything prior to be built already
      MatrixXd inverse = _reducedA.inverse();
      _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
      filename = _reducedPath + string("inverseProduct.matrix");
      EIGEN::write(filename, _inverseProduct);
      TIMER::printTimings();
    }
  }
}

SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::~SUBSPACE_FLUID_3D_COMPRESSED_EIGEN()
{
}

//////////////////////////////////////////////////////////////////////
// The reduced solver, with peeled boundaries, 
// with cubature enabled
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::stepReorderedCubatureStam()
  {
    TIMER functionTimer(__FUNCTION__);
    Real goalTime = 0.1;
    Real currentTime = 0;

    // compute the CFL condition
    _dt = goalTime;

    // wipe forces
    _force.clear();

    // wipe boundaries
    _velocity.setZeroBorder();

    // compute the forces
    addBuoyancy(_heat.data());
    _velocity.axpy(_dt, _force);

    _force.clear();
    addVorticity();
    _velocity.axpy(_dt, _force);

    // VectorXd flat;

    advectHeatAndDensityStam();
    // cout << "did advect heat and density" << endl;
    
    VECTOR::printVertical = false;

    // VECTOR3_FIELD_3D velocityPeeled = _velocity.peelBoundary();
    // VECTOR velocityFlatPeeled = velocityPeeled.flattened();
    // velocityFlatPeeled.write("velocityFlatPeeledUnstable.field");

    TIMER projectionTimer("Velocity projection");

    PeeledCompressedProjectTransformNoSVD(_velocity, &_U_preadvect_data, &_qDot);
    // cout << "did peeled compressed project" << endl;
    // VECTOR qDot = EIGEN::convert(_qDot);
    // cout << "qDot preadvect: " << qDot << "; " << endl;

    projectionTimer.stop();

    reducedAdvectStagedStamFast();
    // cout << "did reduced advect staged" << endl;

    // qDot = EIGEN::convert(_qDot);
    // cout << "qDot postadvect: " << qDot << "; " << endl;
    
    TIMER diffusionProjectionTimer("Velocity projection");
    diffusionProjectionTimer.stop();

    reducedPeeledDiffusion();
    // cout << "did reduced peeled diffusion" << endl;

    // qDot = EIGEN::convert(_qDot);
    // cout << "qDot postdiffuse: " << qDot << "; " << endl;

    reducedStagedProject();

    // cout << "did reduced staged project" << endl;

    // qDot = EIGEN::convert(_qDot);
    // cout << "qDot post reduced project: " << qDot << "; " << endl;
     
    PeeledCompressedUnprojectTransform(&_U_final_data, _qDot, &_velocity);
    // VECTOR unprojectFlat = _velocity.flattened();
    // unprojectFlat.write("unprojectFlatUnstable.vector");

    // qDot = EIGEN::convert(_qDot);
    // cout << "qDot post unproject: " << qDot << "; " << endl;

    // do the full space unprojection
    TIMER unprojectionTimer("Velocity unprojection");

    currentTime += _dt;

    cout << " Simulation step " << _totalSteps << " done. " << endl;

    _totalTime += goalTime;
    _totalSteps++;

    // diff the current sim results against ground truth
    diffGroundTruth();
  }

//////////////////////////////////////////////////////////////////////
// The reduced solver, with IOP, with peeled boundaries, 
// with cubature enabled
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::stepWithObstacle()
{
  TIMER functionTimer(__FUNCTION__);
  Real goalTime = 0.1;
  Real currentTime = 0;

  // compute the CFL condition
  _dt = goalTime;

  // wipe forces
  _force.clear();

  // wipe boundaries
  _velocity.setZeroBorder();

  // compute the forces
  addBuoyancy(_heat.data());
  _velocity.axpy(_dt, _force);

  _force.clear();
  addVorticity();
  _velocity.axpy(_dt, _force);

  TIMER projectionTimer("Velocity projection");

  // project into the subspace
  PeeledCompressedProjectTransformNoSVD(_velocity, &_U_preadvect_data, &_qDot);
  cout << "finished projection! " << endl;
  projectionTimer.stop();

  // then advect

  // full-space advect heat and density
  advectHeatAndDensityStam();
  cout << "finished advect heat and density!" << endl;

  // reduced advect velocity
  //reducedAdvectStagedStamFast();
  //cout << "finished reduced advection!" << endl;

  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  reducedAdvectCompressionFriendly();
  cout << "finished compression-friendly advection!" << endl;
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  
  //cout << " post advection qDot: " << _qDot << endl;
  //exit(0);

  // then diffuse 
  TIMER diffusionProjectionTimer("Reduced diffusion");
  reducedPeeledDiffusion();
  diffusionProjectionTimer.stop();

  // do IOP
  reducedSetZeroSphere();

  // then pressure project
  reducedStagedProject();


  // come back to full space
  TIMER unprojectionTimer("Velocity unprojection");
  PeeledCompressedUnprojectTransform(&_U_final_data, _qDot, &_velocity);
  unprojectionTimer.stop();

  currentTime += _dt;

  cout << " Simulation step " << _totalSteps << " done. " << endl;

  _totalTime += goalTime;
  _totalSteps++;

  // diff the current sim results against ground truth
  diffGroundTruth();
}

//////////////////////////////////////////////////////////////////////
// do a full-rank advection of heat and density
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::advectHeatAndDensityStam()
{
  TIMER functionTimer(__FUNCTION__);

  if(_domainBcLeft == 0) 
    _velocity.copyBorderX();
  else 
    _velocity.setZeroX();

  if(_domainBcTop == 0) 
    _velocity.copyBorderY();
  else 
    _velocity.setZeroY();

  if(_domainBcFront == 0) 
    _velocity.copyBorderZ();
  else 
    _velocity.setZeroZ();

  _density.swapPointers(_densityOld);
  _heat.swapPointers(_heatOld);

  const Real dt0 = _dt / _dx;
  VECTOR3_FIELD_3D::advect(dt0, _velocity, _densityOld, _density);
  VECTOR3_FIELD_3D::advect(dt0, _velocity, _heatOld, _heat);
  
  _density.setZeroBorder();
  _heat.setZeroBorder();
}

//////////////////////////////////////////////////////////////////////
// do a reduced zeroing out of the sphere interior for IOP
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedSetZeroSphere()
{
  TIMER functionTimer(__FUNCTION__);

  cout << "_qDot.size in reducedZeroSphere: " << _qDot.size() << endl;
  if (_reducedIOP.cols() == 0) {
    string filename;
    filename = _reducedPath + string("U.iop.subspace.matrix");
    if (fileExists(filename)) {
      EIGEN::read(filename, _reducedIOP);
    }
  }

  cout << "_reducedIOP.cols: " << _reducedIOP.cols() << endl;
  _qDot = _reducedIOP * _qDot;

  
}
//////////////////////////////////////////////////////////////////////
// perform reduced order diffusion with separate boundary slabs
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedPeeledDiffusion() 
{
  TIMER functionTimer(__FUNCTION__);

  // apply the reduced damping to the middle
  _qDot = _dampingMatrixReduced * _qDot;
}

//////////////////////////////////////////////////////////////////////
// get the projection error of a vector with respect to a basis
//////////////////////////////////////////////////////////////////////
Real SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::projectionError(const MatrixXd& basis, const VectorXd& v)
{
  VectorXd after = basis * (basis.transpose() * v);

  return (v - after).norm();
}

//////////////////////////////////////////////////////////////////////
// compute pressure-to-velocity matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::computePressureToVelocity()
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  int peeledDims = xPeeled * yPeeled * zPeeled;
  int totalCells = 3 * peeledDims;
  _pressureToVelocity.resize(totalCells, peeledDims);

  Real invDx = 1.0f / _dx;
  Real halfInvDx = 0.5f / _dx;
  int index = _slabSize + _xRes + 1;
  for (int z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (int y = 1; y < _yRes - 1; y++, index += 2)
      for (int x = 1; x < _xRes - 1; x++, index++)
      {
        int peeledIndex = (x - 1) + (y - 1) * xPeeled + (z - 1) * slabPeeled;
        int peeledIndex3 = 3 * peeledIndex;

        if (x == 1)
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex + 1) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex)     -= -invDx;
        }
        else if (x == _xRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex)     -=  invDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex - 1) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex + 1) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 0, peeledIndex - 1) -= -halfInvDx;
        }

        if (y == 1)
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex + xPeeled) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex)           -= -invDx;
        }
        else if (y == _yRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex)           -=  invDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex - xPeeled) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex + xPeeled) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 1, peeledIndex - xPeeled) -= -halfInvDx;
        }

        if (z == 1)
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex + slabPeeled) -=  invDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex)              -= -invDx;
        }
        else if (z == _zRes - 2)
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex)              -=  invDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex - slabPeeled) -= -invDx;
        }
        else
        {
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex + slabPeeled) -=  halfInvDx;
          _pressureToVelocity(peeledIndex3 + 2, peeledIndex - slabPeeled) -= -halfInvDx;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// do a staged reduced order pressure projection
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedStagedProject()
{
  TIMER functionTimer(__FUNCTION__);
  cout << "_preprojectToFinal size: " << '(' << _preprojectToFinal.rows() 
    << ", " << _preprojectToFinal.cols() << ")\n";

  cout << "_inverseProduct size: " << '(' << _inverseProduct.rows()
    << ", " << _inverseProduct.cols() << ")\n";

  _qDot = _preprojectToFinal * _qDot + _inverseProduct * _qDot;
}
//////////////////////////////////////////////////////////////////////
// do a staged reduced order pressure projection for IOP
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedStagedProjectIOP()
{
  TIMER functionTimer(__FUNCTION__);
  cout << "_preprojectToPreadvect.cols: " << _preprojectToPreadvect.cols() << endl;
  cout << "_qDot.size: " << _qDot.size() << endl;
  cout << "_inverseProduct rows, cols: " << "(" << _inverseProduct.rows() << ", " << _inverseProduct.cols() << ")" << endl;
  _qDot = _preprojectToPreadvect * _qDot + _inverseProduct * _qDot;
}

//////////////////////////////////////////////////////////////////////
// diff the current sim results against ground truth
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::diffGroundTruth()
{
  TIMER functionTimer(__FUNCTION__);
  if (_fullRankPath.size() == 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " No snapshot path has been set to compare to! " << endl;
    exit(0);
  }

  char buffer[256];
  string path = _fullRankPath;
  sprintf(buffer, "%sfluid.%04i.fluid3d", path.c_str(), _totalSteps - 1);
  string filename(buffer);
  FLUID_3D_MIC ground;
  ground.readGz(filename);

  cout << "=====================================================" << endl;
  cout << " Ground truth difference for frame " << _totalSteps - 1 << endl;
  cout << "=====================================================" << endl;
  VectorXd diff = _velocity.peelBoundary().flattenedEigen() - ground.velocity().peelBoundary().flattenedEigen();
  _velocityErrorAbs.push_back(diff.norm());
  _velocityErrorRelative.push_back(diff.norm() / _velocity.peelBoundary().flattened().norm2());
  cout << " velocity abs error:      " << _velocityErrorAbs.back() << endl;
  cout << " velocity relative error: " << _velocityErrorRelative.back() << endl;

  diff = _density.peelBoundary().flattenedEigen() - ground.density().peelBoundary().flattenedEigen();
  _densityErrorAbs.push_back(diff.norm());
  _densityErrorRelative.push_back(diff.norm() / _density.peelBoundary().flattened().norm2());
  cout << " density abs error:      " << _densityErrorAbs.back() << endl;
  cout << " density relative error: " << _densityErrorRelative.back() << endl;
}

//////////////////////////////////////////////////////////////////////
// get the sub-basis associated with a cell
//////////////////////////////////////////////////////////////////////
MatrixXd SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::cellBasisPeeled(const MatrixXd& U, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // decompose into x,y,z
  const int decompose = index;
  const int z = decompose / _slabPeeled;
  const int y = (decompose % _slabPeeled) / _xPeeled;
  const int x = (decompose % _slabPeeled) % _xPeeled;

  assert(x >= 0);
  assert(x < _xRes - 2);
  assert(y >= 0);
  assert(y < _yRes - 2);
  assert(z >= 0);
  assert(z < _zRes - 2);
  assert(3 * index < U.rows());
  
  return EIGEN::getRows(3 * index, 3, U); 
}

//////////////////////////////////////////////////////////////////////
// get the sub-basis associated with a cell, compressed version
//////////////////////////////////////////////////////////////////////
MatrixXd SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::cellBasisCompressedPeeled(MATRIX_COMPRESSION_DATA& U_data, const int index)
{
  TIMER functionTimer(__FUNCTION__);
  // decompose into x,y,z
  const int decompose = index;
  const int z = decompose / _slabPeeled;
  const int y = (decompose % _slabPeeled) / _xPeeled;
  const int x = (decompose % _slabPeeled) % _xPeeled;

  COMPRESSION_DATA* dataX = U_data.get_compression_dataX();
  const VEC3I& dims = dataX->get_dims();
  const int xRes = dims[0];
  const int yRes = dims[1];
  const int zRes = dims[2];
  const int numRows = 3 * xRes * yRes * zRes;
  const int numCols = dataX->get_numCols();

  assert(x >= 0);
  assert(x < _xRes - 2);
  assert(y >= 0);
  assert(y < _yRes - 2);
  assert(z >= 0);
  assert(z < _zRes - 2);
  assert(3 * index < numRows);
  
  MatrixXd result(3, numCols); 
  GetSubmatrixNoSVD(3 * index, &U_data, &result); 
  return result;
  // return EIGEN::getRows(3 * index, 3, U); 
}


//////////////////////////////////////////////////////////////////////
// advect a single cell
//////////////////////////////////////////////////////////////////////

VectorXd SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::advectCellStamPeeled(MATRIX_COMPRESSION_DATA& U_data, const MatrixXd& cellU,
    Real dt, const VectorXd& qDot, int index, MatrixXd* submatrix)
{
  TIMER functionTimer(__FUNCTION__);
  // peeled coordinates were passed in -- need to promote to full grid
  const int decompose = index;
  const int z = decompose / _slabPeeled + 1;
  const int y = (decompose % _slabPeeled) / _xPeeled + 1;
  const int x = (decompose % _slabPeeled) % _xPeeled + 1;

  // get the velocity to backtrace
  VectorXd v = cellU * qDot;

  // backtrace
  const VEC3F velocity(v[0], v[1], v[2]);
  Real xTrace = x - dt * velocity[0];
  Real yTrace = y - dt * velocity[1];
  Real zTrace = z - dt * velocity[2];

  // clamp backtrace to grid boundaries
  
  // keeping this comment block here for reference
  xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
  xTrace = (xTrace > _xRes - 2.5) ? _xRes - 2.5 : xTrace;
  yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
  yTrace = (yTrace > _yRes - 2.5) ? _yRes - 2.5 : yTrace;
  zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
  zTrace = (zTrace > _zRes - 2.5) ? _zRes - 2.5 : zTrace;

  // locate neighbors to interpolate --
  // since we're in peeled coordinates, the lookup needs to be modified slightly
  
  // keeping this comment block here for reference
  const int x0 = (int)xTrace - 1;
  const int x1 = x0 + 1;
  const int y0 = (int)yTrace - 1;
  const int y1 = y0 + 1;
  const int z0 = (int)zTrace - 1;
  const int z1 = z0 + 1;

  // get interpolation weights
  const Real s1 = (xTrace - 1) - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = (yTrace - 1) - y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = (zTrace - 1) - z0;
  const Real u0 = 1.0f - u1;

  const int z0Scaled = z0 * _slabPeeled;
  const int z1Scaled = z1 * _slabPeeled;
  const int y0Scaled = y0 * _xPeeled;
  const int y1Scaled = y1 * _xPeeled;

  const int i000 = 3 * (x0 + y0Scaled + z0Scaled);
  const int i010 = 3 * (x0 + y1Scaled + z0Scaled);
  const int i100 = 3 * (x1 + y0Scaled + z0Scaled);
  const int i110 = 3 * (x1 + y1Scaled + z0Scaled);
  const int i001 = 3 * (x0 + y0Scaled + z1Scaled);
  const int i011 = 3 * (x0 + y1Scaled + z1Scaled);
  const int i101 = 3 * (x1 + y0Scaled + z1Scaled);
  const int i111 = 3 * (x1 + y1Scaled + z1Scaled);

  // NOTE: it spends most of its time (+50%) here
  
  // COMPRESSION_DATA* dataX = U_data.get_compression_dataX();
  // const int totalCols = dataX->get_numCols();

  TIMER multiplyTimer0("multiplyTimer0");

  GetSubmatrixNoSVD(i000, &U_data, submatrix);
  const VectorXd v000 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i010, &U_data, submatrix);
  const VectorXd v010 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i100, &U_data, submatrix);
  const VectorXd v100 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i110, &U_data, submatrix);
  const VectorXd v110 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i001, &U_data, submatrix);
  const VectorXd v001 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i011, &U_data, submatrix);
  const VectorXd v011 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i101, &U_data, submatrix);
  const VectorXd v101 = (*submatrix) * qDot;

  GetSubmatrixNoSVD(i111, &U_data, submatrix);
  const VectorXd v111 = (*submatrix) * qDot;

  multiplyTimer0.stop();

  const Real w000 = u0 * s0 * t0;
  const Real w010 = u0 * s0 * t1;
  const Real w100 = u0 * s1 * t0;
  const Real w110 = u0 * s1 * t1;
  const Real w001 = u1 * s0 * t0;
  const Real w011 = u1 * s0 * t1;
  const Real w101 = u1 * s1 * t0;
  const Real w111 = u1 * s1 * t1;
  
  // interpolate
  // (indices could be computed once)
  //
  // NOTE: it's deceptive to think this cuts down on
  // multiplies, since they will all occur on a VECTOR,
  // not just a scalar
  //
  //return u0 * (s0 * (t0 * v000 + t1 * v010) +
  //             s1 * (t0 * v100 + t1 * v110)) +
  //       u1 * (s0 * (t0 * v001 + t1 * v011) +
  //             s1 * (t0 * v101 + t1 * v111));
  return w000 * v000 + w010 * v010 + w100 * v100 + w110 * v110 +
         w001 * v001 + w011 * v011 + w101 * v101 + w111 * v111;
  }





//////////////////////////////////////////////////////////////////////
// read in a cubature scheme
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::readAdvectionCubature()
{
  string filename = _reducedPath + string("cubature");
  cout << " Trying to read in cubature file  " << filename.c_str() << " ... "; flush(cout);

  // read in the cubature file
  FILE* file;
  file = fopen(filename.c_str(), "rb");

  if (file == NULL)
  {
    cout << " No cubature file " << filename.c_str() << " found! " << endl;
    cout << " Maybe you're still training? " << endl;
    return;
  }

  // read dimensions
  int size;
  fread((void*)&size, sizeof(int), 1, file);

  // read out the key tet indices
  int cellIndex;
  for (int x = 0; x < size; x++)
  {
    fread((void*)&cellIndex, sizeof(int), 1, file);
    _keyAdvectionCells.push_back(cellIndex);
  }

  // read out the key tet weights
  double weight;
  for (int x = 0; x < size; x++)
  {
    fread((void*)&weight, sizeof(double), 1, file);
    _keyAdvectionWeights.push_back(weight);
  }

  VECTOR weightVector(_keyAdvectionWeights);
  VECTOR::printVertical = false;
  cout << " Read in " << _keyAdvectionCells.size() << " cubature points " << endl;
  fclose(file);

  cout << " done." << endl;

  // look for precached advection matrices
  filename = _reducedPath + string("advection.cubature.cache");
  file = fopen(filename.c_str(), "rb");

  int totalPoints = _keyAdvectionCells.size();

  // if there's no cache, create one
  if (file == NULL)
  {
    //fclose(file);
    cout << " No advection cubature cache found. Precaching all the advection cubature matrices. " << endl;

    // make the before cache
    _advectionCubatureBefore.clear();
    // TODO: insert the decoder for U.preadvect

  
    COMPRESSION_DATA compression_data0;
    COMPRESSION_DATA compression_data1;
    COMPRESSION_DATA compression_data2;

    // const char* filename = "U.preadvect.SVD.data";
    // ReadSVDData(filename, &compression_data0);

    string preadvectFile = _reducedPath + string("U.preadvect.component0");
    int* allData0 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data0);
    preadvectFile = _reducedPath + string("U.preadvect.component1");
    int* allData1 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data1);
    preadvectFile = _reducedPath + string("U.preadvect.component2");
    int* allData2 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data2);


    _U_preadvect_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
        &compression_data0, &compression_data1, &compression_data2); 

    _U_preadvect_data.dct_setup(-1);
    _U_preadvect_data.init_cache();

    // EIGEN::read(preadvectFile, _preadvectU);
    //
    for (int x = 0; x < totalPoints; x++)
    {
      const int index = _keyAdvectionCells[x];
      _advectionCubatureBefore.push_back(cellBasisCompressedPeeled(_U_preadvect_data, index));
    }
    // _preadvectU.resize(0,0);

    // make the after cache
    _advectionCubatureAfter.clear();
    string prediffuseFile = _reducedPath + string("U.prediffuse.matrix");
    EIGEN::read(prediffuseFile, _prediffuseU);
    for (int x = 0; x < totalPoints; x++)
    {
      const int index = _keyAdvectionCells[x];
      const Real weight = _keyAdvectionWeights[x];
      
      const MatrixXd cellU = cellBasisPeeled(_prediffuseU, index);
      _advectionCubatureAfter.push_back(weight * cellU.transpose());
    }
    _prediffuseU.resize(0,0);

    // write the results out to a file
    filename = _reducedPath + string("advection.cubature.cache");
    file = fopen(filename.c_str(), "wb");

    if (file == NULL)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Couldn't open advection cubature cache file " << filename.c_str() << "!!! " << endl;
      exit(0);
    }

    for (int x = 0; x < totalPoints; x++)
      EIGEN::write(file, _advectionCubatureBefore[x]);
    for (int x = 0; x < totalPoints; x++)
      EIGEN::write(file, _advectionCubatureAfter[x]);
    fclose(file);
  }
  // if there's already a cache, read it in
  else
  {
    cout << " Found advection cubature cache! Reading ... " << flush;
    _advectionCubatureBefore.resize(_keyAdvectionCells.size());
    _advectionCubatureAfter.resize(_keyAdvectionCells.size());

    for (int x = 0; x < totalPoints; x++)
      EIGEN::read(file, _advectionCubatureBefore[x]);
    for (int x = 0; x < totalPoints; x++)
      EIGEN::read(file, _advectionCubatureAfter[x]);
    fclose(file);

    cout << " done. " << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// do Stam-style advection using cubautre
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedAdvectStagedStamFast()
{
  TIMER functionTimer(__FUNCTION__);
  VectorXd final(_advectionCubatureAfter[0].rows());
  final.setZero();
  const Real dt0 = _dt / _dx;

  vector<VectorXd> finals(_keyAdvectionCells.size());
  int totalPoints = _keyAdvectionCells.size();
  
  assert(_advectionCubatureAfter.size() > 0);
  assert(_advectionCubatureBefore.size() > 0);

  MATRIX_COMPRESSION_DATA& data = (*this).U_final_data();
  COMPRESSION_DATA* dataX = data.get_compression_dataX();
  const int numCols = dataX->get_numCols();
  
  MatrixXd submatrix(3, numCols);

  // loop through the cubature points
#pragma omp parallel
#pragma omp for  schedule(static)
  for (int x = 0; x < totalPoints; x++)
  {
    const int index = _keyAdvectionCells[x];
    // finals[x] = _advectionCubatureAfter[x] * advectCellStamPeeled(_preadvectU, _advectionCubatureBefore[x], dt0, _qDot, index);
    finals[x] = _advectionCubatureAfter[x] * advectCellStamPeeled(_U_preadvect_data, _advectionCubatureBefore[x], dt0, _qDot, index, &submatrix);
  }

  for (int x = 0; x < totalPoints; x++)
    final += finals[x];
  
  _qDot = final;
}

/*
//////////////////////////////////////////////////////////////////////
// build the staged bases and the projected matrices
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::buildOutOfCoreMatrices()
{
  TIMER functionTimer(__FUNCTION__);

  cout << "=====================================================" << endl;
  cout << " Building huge projected matrices " << endl;
  cout << "=====================================================" << endl;

  // in preparation for writing results to files
  string filename;

  TIMER pressureTimer("Pressure projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  int totalCells = (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseA(totalCells, totalCells);
  buildFlatA(sparseA, _obstacles);
  cout << " Projecting the pressure matrix ... " << flush;
  _reducedA = sparseA.projectVerySparse(_pressureU, _pressureU);
  sparseA.clear();
  cout << " done." << endl;

  filename = _reducedPath + string("projected.A.matrix");
  EIGEN::write(filename, _reducedA);
  _Asparse.clear();
  purge();
  pressureTimer.stop();
  TIMER::printTimings();

  TIMER vtodTimer("Velocity to div projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.preproject.matrix");
  EIGEN::read(filename, _preprojectU);

  // build reduced velocity to divergence
  // Need: _pressureU and _preprojectU
  computeVelocityToDivergence();
  cout << " Projecting velocity to divergence ... " << flush;
  _reducedVelocityToDivergence = _velocityToDivergence.project(_pressureU, _preprojectU);
  cout << " done." << endl;

  filename = _reducedPath + string("projected.vtod.matrix");
  EIGEN::write(filename, _reducedVelocityToDivergence);
  _velocityToDivergence.clear();
   
  // stomp matrices to make room in memory
  _pressureU.resize(0,0);
  purge();
  vtodTimer.stop();
  TIMER::printTimings();

  TIMER dampingTimer("Damping matrix projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);

  // NEW
  cout << " Projecting damping matrix ... " << flush;
  int totalCellsD = 3 * (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  SPARSE_MATRIX_ARRAY sparseD(totalCellsD, totalCellsD);
  buildPeeledDampingMatrixFlat(sparseD);
  _dampingMatrixReduced = sparseD.projectVerySparse(_preprojectU, _prediffuseU);
  filename = _reducedPath + string("damping.peeled.matrix");
  EIGEN::write(filename, _dampingMatrixReduced);
  sparseD.clear();
  cout << "done. " << endl;

  // stomp matrices to make room in memory
  _prediffuseU.resize(0,0);
  purge();
  dampingTimer.stop();
  TIMER::printTimings();

  TIMER preprojectTimer("Preprojection projection");
  // read in the needed matrices

  
  // filename = _reducedPath + string("U.final.matrix");
  // EIGEN::readBig(filename, _U);
  // EIGEN::read(filename, _U);
  
  
  int* allDataX = NULL;
  int* allDataY = NULL;
  int* allDataZ = NULL;
  DECOMPRESSION_DATA decompression_dataX;
  DECOMPRESSION_DATA decompression_dataY;
  DECOMPRESSION_DATA decompression_dataZ;
  filename = _reducedPath + string("U.final.componentX");
  ReadBinaryFileToMemory(filename.c_str(), allDataX, decompression_dataX);
  filename = _reducedPath + string("U.final.componentY");
  ReadBinaryFileToMemory(filename.c_str(), allDataY, decompression_dataY);
  filename = _reducedPath + string("U.final.componentZ");
  ReadBinaryFileToMemory(filename.c_str(), allDataZ, decompression_dataZ);
  MATRIX_COMPRESSION_DATA U_data(allDataX, allDataY, allDataZ,
      decompression_dataX, decompression_dataY, decompression_dataZ);



  // need: preproject and U
  EIGEN::transposeProduct(_U, _preprojectU, _preprojectToFinal);
  filename = _reducedPath + string("projected.ptof.matrix");
  EIGEN::write(filename, _preprojectToFinal);

  // stomp matrices to make room in memory
  _preprojectU.resize(0,0);
  purge();
  preprojectTimer.stop();
  TIMER::printTimings();

  TIMER ptovTimer("Pressure to velocity projection");
  // read in the needed matrices
  filename = _reducedPath + string("U.pressure.matrix");
  EIGEN::read(filename, _pressureU);

  // build reduced pressure to velocity
  // Need: pressureU and U
  computePressureToVelocity();
  _reducedPressureToVelocity = _pressureToVelocity.project(_U, _pressureU);
  filename = _reducedPath + string("projected.ptov.matrix");
  EIGEN::write(filename, _reducedPressureToVelocity);

  // stomp pressure, just in case
  _pressureU.resize(0,0);
  purge();
  ptovTimer.stop();
  TIMER::printTimings();

  // Needs everything prior to be built already
  MatrixXd inverse = _reducedA.inverse();
  _inverseProduct = _reducedPressureToVelocity * inverse * _reducedVelocityToDivergence;
  filename = _reducedPath + string("inverseProduct.matrix");
  EIGEN::write(filename, _inverseProduct);
  TIMER::printTimings();

  // clear out all memory, just to be sure
  stompAllBases();
  purge();
  
  cout << " Done building matrices " << endl;
  TIMER::printTimings();
}
*/

//////////////////////////////////////////////////////////////////////
// check of a file exists
//////////////////////////////////////////////////////////////////////
bool SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::fileExists(const string& filename)
{
  FILE* file;
  file = fopen(filename.c_str(), "rb");
  
  if (file == NULL)
    return false;

  fclose(file);
  return true;
}

//////////////////////////////////////////////////////////////////////
// stomp all loaded bases
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::stompAllBases()
{
  _pressureU.resize(0,0);
  _preprojectU.resize(0,0);
  _prediffuseU.resize(0,0);
  // _preadvectU.resize(0,0);
  // _U.resize(0,0);
}

//////////////////////////////////////////////////////////////////////
// init the cache in the matrix compression data members
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::initCompressionData()
{
  _U_final_data.dct_setup(-1);
  _U_final_data.init_cache();
  _U_preadvect_data.dct_setup(-1);
  _U_preadvect_data.init_cache();
}

//////////////////////////////////////////////////////////////////////
// stomp the other matrices and load the ones needed for cubature 
// training
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::loadCubatureTrainingBases()
{
  string filename;

  // TODO: insert the decoder for U.preadvect 
   
  COMPRESSION_DATA compression_data0;
  COMPRESSION_DATA compression_data1;
  COMPRESSION_DATA compression_data2;

  // filename = string("U.preadvect.SVD.data");
  // ReadSVDData(filename.c_str(), &compression_data0);

  string preadvectFile = _reducedPath + string("U.preadvect.component0");
  int* allData0 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data0);
  preadvectFile = _reducedPath + string("U.preadvect.component1");
  int* allData1 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data1);
  preadvectFile = _reducedPath + string("U.preadvect.component2");
  int* allData2 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data2);

  _U_preadvect_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
      &compression_data0, &compression_data1, &compression_data2); 

  _U_preadvect_data.dct_setup(-1);
  _U_preadvect_data.init_cache();

  // EIGEN::read(filename, _preadvectU);

  filename = _reducedPath + string("U.prediffuse.matrix");
  EIGEN::read(filename, _prediffuseU);
}

//////////////////////////////////////////////////////////////////////
// load the bases needed for cubature runtime
//
// the path variable is there in case we want to load off the SSD
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::loadReducedRuntimeBases(string path)
{
  TIMER functionTimer(__FUNCTION__);
  if (path.length() == 0)
    path = _reducedPath;

  string filename;
  
  // TODO: insert the decoder for U.preadvect
  COMPRESSION_DATA compression_data0;
  COMPRESSION_DATA compression_data1;
  COMPRESSION_DATA compression_data2;

  // filename = string("U.preadvect.SVD.data");
  // ReadSVDData(filename, &compression_data0);

  string preadvectFile = _reducedPath + string("U.preadvect.component0");
  int* allData0 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data0);
  preadvectFile = _reducedPath + string("U.preadvect.component1");
  int* allData1 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data1);
  preadvectFile = _reducedPath + string("U.preadvect.component2");
  int* allData2 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data2);

  _U_preadvect_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
      &compression_data0, &compression_data1, &compression_data2); 

  _U_preadvect_data.dct_setup(-1);
  _U_preadvect_data.init_cache();


  COMPRESSION_DATA* preadvect_dataX = _U_preadvect_data.get_compression_dataX();
  const VEC3I& dims = preadvect_dataX->get_dims();
  int numRows = 3 * dims[0] * dims[1] * dims[2]; 

  // EIGEN::read(filename, _preadvectU);
  // if (_preadvectU.rows() > 1000000)

  if (numRows > 1000000) 
    purge();
  
  COMPRESSION_DATA final_compression_data0;
  COMPRESSION_DATA final_compression_data1;
  COMPRESSION_DATA final_compression_data2;

  // filename = string("U.final.SVD.data");
  // ReadSVDData(filename, &compression_data0);

  string finalFile = _reducedPath + string("U.final.component0");
  allData0 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data0);
  finalFile = _reducedPath + string("U.final.component1");
  allData1 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data1);
  finalFile = _reducedPath + string("U.final.component2");
  allData2 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data2);


  _U_final_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
      &final_compression_data0, &final_compression_data1, &final_compression_data2); 

  _U_final_data.dct_setup(-1);
  _U_final_data.init_cache();

  // EIGEN::read(filename, _U);
  
  COMPRESSION_DATA* final_dataX = _U_final_data.get_compression_dataX();
  const VEC3I& dimsFinal = final_dataX->get_dims();
  numRows = 3 * dimsFinal[0] * dimsFinal[1] * dimsFinal[2]; 

  // EIGEN::read(filename, _preadvectU);
  // if (_U.rows() > 1000000)
  if (numRows > 1000000) 
    purge();

  TIMER::printTimings();
}

//////////////////////////////////////////////////////////////////////
// load the IOP bases needed for cubature runtime
//
// the path variable is there in case we want to load off the SSD
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::loadReducedIOP(string path, bool debug)
{

  // we need to load preproject, preadvect, and final (??)
  // nope, just preadvect and final!!
  
 TIMER functionTimer(__FUNCTION__);
  if (path.length() == 0)
    path = _reducedPath;

  string filename;
  
  COMPRESSION_DATA compression_data0;
  COMPRESSION_DATA compression_data1;
  COMPRESSION_DATA compression_data2;

  // filename = string("U.preadvect.SVD.data");
  // ReadSVDData(filename, &compression_data0);

  string preadvectFile = _reducedPath + string("U.preadvect.component0");
  int* allData0 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data0);
  preadvectFile = _reducedPath + string("U.preadvect.component1");
  int* allData1 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data1);
  preadvectFile = _reducedPath + string("U.preadvect.component2");
  int* allData2 = ReadBinaryFileToMemory(preadvectFile.c_str(), &compression_data2);

  compression_data0.set_debug(debug);
  compression_data1.set_debug(debug);
  compression_data2.set_debug(debug);

  _U_preadvect_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
      &compression_data0, &compression_data1, &compression_data2); 

  _U_preadvect_data.dct_setup(-1);
  _U_preadvect_data.init_cache();


  COMPRESSION_DATA* preadvect_dataX = _U_preadvect_data.get_compression_dataX();
  const VEC3I& dims = preadvect_dataX->get_dims();
  int numRows = 3 * dims[0] * dims[1] * dims[2]; 

  // EIGEN::read(filename, _preadvectU);
  // if (_preadvectU.rows() > 1000000)

  if (numRows > 1000000) 
    purge();

  filename = path + string("U.final.matrix");
  
  COMPRESSION_DATA final_compression_data0;
  COMPRESSION_DATA final_compression_data1;
  COMPRESSION_DATA final_compression_data2;

  // filename = string("U.final.SVD.data");
  // ReadSVDData(filename, &compression_data0);

  string finalFile = _reducedPath + string("U.final.component0");
  allData0 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data0);
  finalFile = _reducedPath + string("U.final.component1");
  allData1 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data1);
  finalFile = _reducedPath + string("U.final.component2");
  allData2 = ReadBinaryFileToMemory(finalFile.c_str(), &final_compression_data2);


  final_compression_data0.set_debug(debug);
  final_compression_data1.set_debug(debug);
  final_compression_data2.set_debug(debug);

  _U_final_data = MATRIX_COMPRESSION_DATA(allData0, allData1, allData2,
      &final_compression_data0, &final_compression_data1, &final_compression_data2); 

  _U_final_data.dct_setup(-1);
  _U_final_data.init_cache();

  // EIGEN::read(filename, _U);
  
  COMPRESSION_DATA* final_dataX = _U_final_data.get_compression_dataX();
  const VEC3I& dimsFinal = final_dataX->get_dims();
  numRows = 3 * dimsFinal[0] * dimsFinal[1] * dimsFinal[2]; 

  // EIGEN::read(filename, _preadvectU);
  // if (_U.rows() > 1000000)
  if (numRows > 1000000) 
    purge();

  
  /*
  filename = path + string("U.preproject.matrix");
  UallDataX = NULL;
  UallDataY = NULL;
  UallDataZ = NULL;
  // rewrite into Udecompression_data
  filename = _reducedPath + string("U.preproject.componentX");
  ReadBinaryFileToMemory(filename.c_str(), UallDataX, Udecompression_dataX);
  filename = _reducedPath + string("U.preproject.componentY");
  ReadBinaryFileToMemory(filename.c_str(), UallDataY, Udecompression_dataY);
  filename = _reducedPath + string("U.preproject.componentZ");
  ReadBinaryFileToMemory(filename.c_str(), UallDataZ, Udecompression_dataZ);
  MATRIX_COMPRESSION_DATA U_preproject_data(UallDataX, UallDataY, UallDataZ,
      Udecompression_dataX, Udecompression_dataY, Udecompression_dataZ);
  _U_preproject_data = U_preproject_data;
  _U_preproject_data.init_cache();
  */

  TIMER::printTimings();
}
//////////////////////////////////////////////////////////////////////
// compute the velocity-to-divergence matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::computeVelocityToDivergence()
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  int peeledDims = xPeeled * yPeeled * zPeeled;
  int totalCells = 3 * peeledDims;

  Real coeff = -_dx * 0.5;
  _velocityToDivergence.resize(peeledDims, totalCells);
  _velocityToDivergence.clearAndStompSparsity();

  for (int z = 0; z < zPeeled; z++)
    for (int y = 0; y < yPeeled; y++)
      for (int x = 0; x < xPeeled; x++)
      {
        int index = x + y * xPeeled + z * slabPeeled;
        int right = index + 1;
        int left  = index - 1;
        int up   = index + xPeeled;
        int down = index - xPeeled;
        int far  = index + slabPeeled;
        int near = index - slabPeeled;

        // this needs to be redone if Neumann is being used
        //
        // the zero will get constant folded away -- it's just there
        // as a reminder.
        if (x != xPeeled - 1)
          _velocityToDivergence(index, 3 * right + 0) += coeff;
        else if (_domainBcLeft != 0)
          _velocityToDivergence(index, 3 * index + 0) += -coeff;
        else if (_domainBcLeft == 0)
          _velocityToDivergence(index, 3 * left + 0) += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (x != 0)
          _velocityToDivergence(index, 3 * left + 0)  += -coeff;
        else if (_domainBcLeft != 0)
          _velocityToDivergence(index, 3 * index + 0) += coeff;
        else if (_domainBcLeft == 0)
          _velocityToDivergence(index, 3 * right + 0) += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (y != yPeeled - 1)
          _velocityToDivergence(index, 3 * up + 1)   += coeff;
        else if (_domainBcTop == 0)
          _velocityToDivergence(index, 3 * down + 1) += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (y != 0)
          _velocityToDivergence(index, 3 * down + 1) += -coeff;
        else if (_domainBcTop == 0)
          _velocityToDivergence(index, 3 * up + 1)   += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }
        
        if (z != zPeeled - 1)
          _velocityToDivergence(index, 3 * far + 2)  += coeff;
        else if (_domainBcFront != 0)
          _velocityToDivergence(index, 3 * index + 2) += -coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }

        if (z != 0)
          _velocityToDivergence(index, 3 * near + 2) += -coeff;
        else if (_domainBcFront != 0)
          _velocityToDivergence(index, 3 * index + 2)  += coeff;
        else
        {
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
          cout << " FELL THROUGH! " << endl;
          exit(0);
        }
      }
}

//////////////////////////////////////////////////////////////////////
// build the peeled damping matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::buildPeeledDampingMatrixFlat(SPARSE_MATRIX_ARRAY& peeledDampingMatrix)
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Building flat peeled damping matrix ... ";flush(cout);
  const Real w = 0.9;
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  const Real wMinus = 1.0 - w;
  const Real sixth = (1./ 6.)* w;

  // initially everybody is identity
  for (int x = 0; x < peeledDampingMatrix.rows(); x++)
    peeledDampingMatrix(x,x) = 1.0;

  // populate the rest
  for (int z = 0; z < zPeeled; z++) 
    for (int y = 0; y < yPeeled; y++)
      for (int x = 0; x < xPeeled; x++) 
      {
        int index = 3 * (x + y * xPeeled + z * slabPeeled);
        int left  = 3 * ((x + 1) + y * xPeeled + z * slabPeeled);
        int right = 3 * ((x - 1) + y * xPeeled + z * slabPeeled);
        int up   = 3 * (x + (y + 1) * xPeeled + z * slabPeeled);
        int down = 3 * (x + (y - 1) * xPeeled + z * slabPeeled);
        int near = 3 * (x + y * xPeeled + (z + 1) * slabPeeled);
        int far  = 3 * (x + y * xPeeled + (z - 1) * slabPeeled);

        for (int i = 0; i < 3; i++)
        {
          peeledDampingMatrix(index + i, index + i) = wMinus;

          if (x != xPeeled - 1)
            peeledDampingMatrix(index + i, left + i) = sixth;

          if (x != 0)
            peeledDampingMatrix(index + i, right + i) = sixth;

          if (y != yPeeled - 1)
            peeledDampingMatrix(index + i, up + i)    = sixth;

          if (y != 0)
            peeledDampingMatrix(index + i, down + i)  = sixth;

          if (z != zPeeled - 1)
            peeledDampingMatrix(index + i, near + i)  = sixth;

          if (z != 0)
            peeledDampingMatrix(index + i, far + i)   = sixth;
        }
      }

  cout << " done. " << endl;
}

//////////////////////////////////////////////////////////////////////
// build a sparse version of the Poisson matrix
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::buildFlatA(SPARSE_MATRIX_ARRAY& sparseA, unsigned char* skip)
{
  int xPeeled = _xRes - 2;
  int yPeeled = _yRes - 2;
  //int zPeeled = _zRes - 2;
  int slabPeeled = xPeeled * yPeeled;
  //int peeledDims = xPeeled * yPeeled * zPeeled;

  int x, y, z, index;
  Real Aoff = 1.0;
  index = _slabSize + _xRes + 1;
  for (z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
      {
        // if the cell is a variable
        if (!skip[index])
        {
          int peeledIndex = (x - 1) + (y - 1) * xPeeled + (z - 1) * slabPeeled;

          // set the matrix to the Poisson stencil in order
          if (!skip[index + _xs]) sparseA(peeledIndex, peeledIndex) += Aoff;
          if (!skip[index - _xs]) sparseA(peeledIndex, peeledIndex) += Aoff;
          if (!skip[index + _ys]) sparseA(peeledIndex, peeledIndex) += Aoff;
          if (!skip[index - _ys]) sparseA(peeledIndex, peeledIndex) += Aoff;
          if (!skip[index + _zs]) sparseA(peeledIndex, peeledIndex) += Aoff;
          if (!skip[index - _zs]) sparseA(peeledIndex, peeledIndex) += Aoff;

          if (!skip[index - 1] && x != 1)
            sparseA(peeledIndex, peeledIndex - 1) = -Aoff;
          if (!skip[index + 1] && x != _xRes - 2)
            sparseA(peeledIndex, peeledIndex + 1) = -Aoff;
          if (!skip[index - _xRes] && y != 1)
            sparseA(peeledIndex, peeledIndex - xPeeled) = -Aoff;
          if (!skip[index + _xRes] && y != _yRes - 2)
            sparseA(peeledIndex, peeledIndex + xPeeled) = -Aoff;
          if (!skip[index - _slabSize] && z != 1)
            sparseA(peeledIndex, peeledIndex - slabPeeled) = -Aoff;
          if (!skip[index + _slabSize] && z != _zRes - 2)
            sparseA(peeledIndex, peeledIndex + slabPeeled) = -Aoff;
        }
      }
}

//////////////////////////////////////////////////////////////////////
// do Stam-style advection using cubautre
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::reducedAdvectCompressionFriendly()
{
  TIMER functionTimer(__FUNCTION__);
  VectorXd final(_advectionCubatureAfter[0].rows());
  final.setZero();
  const Real dt0 = _dt / _dx;

  vector<VectorXd> finals(_keyAdvectionCells.size());
  int totalPoints = _keyAdvectionCells.size();

  assert(_advectionCubatureAfter.size() > 0);
  assert(_advectionCubatureBefore.size() > 0);

  MATRIX_COMPRESSION_DATA& data = (*this).U_final_data();
  COMPRESSION_DATA* dataX = data.get_compression_dataX();
  const int numCols = dataX->get_numCols();
  const VEC3I& dims = dataX->get_dims();

  // the blocks that are needed
  //
  // the first int is the block number
  // the second is all the data you need once you have the block decoded
  multimap<int, CUBATURE_DATA> requestedBlocks;

  // accumulate the blocks you need
  for (int x = 0; x < totalPoints; x++)
  {
    const int index = _keyAdvectionCells[x];
    accumAdvectRequests(_advectionCubatureBefore[x], dt0, _qDot, x, index, dims, requestedBlocks);
  }

  // prepare an array for the results of the advection
  vector<VectorXd> advected(totalPoints);
  for (unsigned int x = 0; x < totalPoints; x++)
  {
    advected[x].resize(3);
    advected[x].setZero();
  }

  // decompress each block in turn
  for (auto iter = requestedBlocks.begin(); iter != requestedBlocks.end(); iter++)
  {
    // unpack the block data
    //int currentBlock = (*iter).first;
    CUBATURE_DATA cubatureData = (*iter).second;
    int cubaturePoint = cubatureData.index;
    int ixxx= cubatureData.ixxx;
    Real wxxx= cubatureData.wxxx;

    // get the block index for the block we're looking at
    int blockNumber, blockIndex;
    ComputeBlockNumber(ixxx, dims, &blockNumber, &blockIndex);

    // caching should work ideally now. Since the multimap sorts by the block number,
    // you will only get a miss when you're done using the block for the whole timestep
    MatrixXd submatrix(3, numCols);
    //GetSubmatrixNoSVD(ixxx, &_U_preadvect_data, &submatrix);
    GetSubmatrixNoSVDSparse(ixxx, &_U_preadvect_data, &submatrix);
    
    // compute the Uq product
    const VectorXd vxxx = submatrix * _qDot;
    
    // add it to the current VEC3 for the cubature point
    advected[cubaturePoint] += wxxx * vxxx;
  }

  for (int x = 0; x < totalPoints; x++)
    finals[x] = _advectionCubatureAfter[x] * advected[x];

  for (int x = 0; x < totalPoints; x++)
    final += finals[x];
  
  _qDot = final;
}

//////////////////////////////////////////////////////////////////////
// accumulate decoder requests for a single cell
//////////////////////////////////////////////////////////////////////
void SUBSPACE_FLUID_3D_COMPRESSED_EIGEN::accumAdvectRequests(const MatrixXd& cellU, 
      const Real dt, const VectorXd& qDot, const int cubatureIndex, const int index, 
      const VEC3I& dims, multimap<int, CUBATURE_DATA>& requestedBlocks)
{
  //TIMER functionTimer(__FUNCTION__);
  // peeled coordinates were passed in -- need to promote to full grid
  const int decompose = index;
  const int z = decompose / _slabPeeled + 1;
  const int y = (decompose % _slabPeeled) / _xPeeled + 1;
  const int x = (decompose % _slabPeeled) % _xPeeled + 1;

  // get the velocity to backtrace
  VectorXd v = cellU * qDot;

  // backtrace
  const VEC3F velocity(v[0], v[1], v[2]);
  Real xTrace = x - dt * velocity[0];
  Real yTrace = y - dt * velocity[1];
  Real zTrace = z - dt * velocity[2];

  // clamp backtrace to grid boundaries
  
  // keeping this comment block here for reference
  xTrace = (xTrace < 1.5) ? 1.5 : xTrace;
  xTrace = (xTrace > _xRes - 2.5) ? _xRes - 2.5 : xTrace;
  yTrace = (yTrace < 1.5) ? 1.5 : yTrace;
  yTrace = (yTrace > _yRes - 2.5) ? _yRes - 2.5 : yTrace;
  zTrace = (zTrace < 1.5) ? 1.5 : zTrace;
  zTrace = (zTrace > _zRes - 2.5) ? _zRes - 2.5 : zTrace;

  // locate neighbors to interpolate --
  // since we're in peeled coordinates, the lookup needs to be modified slightly
  
  // keeping this comment block here for reference
  const int x0 = (int)xTrace - 1;
  const int x1 = x0 + 1;
  const int y0 = (int)yTrace - 1;
  const int y1 = y0 + 1;
  const int z0 = (int)zTrace - 1;
  const int z1 = z0 + 1;

  // get interpolation weights
  const Real s1 = (xTrace - 1) - x0;
  const Real s0 = 1.0f - s1;
  const Real t1 = (yTrace - 1) - y0;
  const Real t0 = 1.0f - t1;
  const Real u1 = (zTrace - 1) - z0;
  const Real u0 = 1.0f - u1;

  const int z0Scaled = z0 * _slabPeeled;
  const int z1Scaled = z1 * _slabPeeled;
  const int y0Scaled = y0 * _xPeeled;
  const int y1Scaled = y1 * _xPeeled;

  const int i000 = 3 * (x0 + y0Scaled + z0Scaled);
  const int i010 = 3 * (x0 + y1Scaled + z0Scaled);
  const int i100 = 3 * (x1 + y0Scaled + z0Scaled);
  const int i110 = 3 * (x1 + y1Scaled + z0Scaled);
  const int i001 = 3 * (x0 + y0Scaled + z1Scaled);
  const int i011 = 3 * (x0 + y1Scaled + z1Scaled);
  const int i101 = 3 * (x1 + y0Scaled + z1Scaled);
  const int i111 = 3 * (x1 + y1Scaled + z1Scaled);

  //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
  //cout << " i000: " << i000 << endl;

  const Real w000 = u0 * s0 * t0;
  const Real w010 = u0 * s0 * t1;
  const Real w100 = u0 * s1 * t0;
  const Real w110 = u0 * s1 * t1;
  const Real w001 = u1 * s0 * t0;
  const Real w011 = u1 * s0 * t1;
  const Real w101 = u1 * s1 * t0;
  const Real w111 = u1 * s1 * t1;

  CUBATURE_DATA cubatureData;
  pair<int, CUBATURE_DATA> requestPair;
  int blockNumber, blockIndex;

  ComputeBlockNumber(i000, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i000;
  cubatureData.wxxx = w000;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i010, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i010;
  cubatureData.wxxx = w010;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i100, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i100;
  cubatureData.wxxx = w100;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i110, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i110;
  cubatureData.wxxx = w110;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i001, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i001;
  cubatureData.wxxx = w001;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i011, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i011;
  cubatureData.wxxx = w011;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i101, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i101;
  cubatureData.wxxx = w101;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);

  ComputeBlockNumber(i111, dims, &blockNumber, &blockIndex);
  cubatureData.index = cubatureIndex;
  cubatureData.ixxx = i111;
  cubatureData.wxxx = w111;
  requestPair = pair<int, CUBATURE_DATA>(blockNumber, cubatureData);
  requestedBlocks.insert(requestPair);
}

/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cylinder.h"
#include "cylinder.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "poiseuille.h"
#include "poiseuille.hh"
#ifdef HDF5
#include "io/xdmfDataOutput.h"
#include "io/xdmfDataOutput.hh"
#include "io/hdfWrapper.h"
#include "io/hdfWrapper.hh"
#endif


using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR MRTD2Q9Descriptor
typedef MRTdynamics<T, DESCRIPTOR> BackgroundDynamics;

// Class to define a domain from a native 2D boolean mask
template <typename T>
class MaskShapeDomain2D : public DomainFunctional2D {
public:
    // Constructor: takes a native 2D mask (vector of vector of bool)
    MaskShapeDomain2D(const std::vector<std::vector<bool>>& mask_)
        : mask(mask_)
    { }

    // operator(): returns true if the mask is active at a given lattice index
    virtual bool operator()(plint iX, plint iY) const override
    {
        if (iY >= 0 && iY < static_cast<plint>(mask.size()) &&
            iX >= 0 && iX < static_cast<plint>(mask[iY].size()))
        {
            return mask[iY][iX];
        }
        return false; // out of bounds -> inactive
    }

    // clone() required by Palabos
    virtual MaskShapeDomain2D<T>* clone() const override
    {
        return new MaskShapeDomain2D<T>(*this);
    }

private:
    std::vector<std::vector<bool>> mask;
};

/// A functional, used to initialize a pressure boundary to constant density
template <typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_) : density(density_) { }
    T operator()(plint, plint) const
    {
        return density;
    }

private:
    T density;
};

/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
void defineGeometry(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters,
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition, Array<plint, 2> forceIds)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    // Define the different boundary boxes
    Box2D inlet(0, 0, 1, ny-2);
    Box2D outlet(nx-1, nx-1, 1, ny-2);
    Box2D bottomWall(0, nx-1, 0, 0);
    Box2D topWall(0, nx-1, ny-1, ny-1); // Careful here to not have this wall step on the outlet

    // Implement boundary conditions
    // So far: top and bottom walls are treated as inlets
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, bottomWall);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, topWall);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, outlet, boundary::outflow);

    // Define the value of the imposed velocity on all inlet BC
    T u = parameters.getLatticeU();
    setBoundaryVelocity(lattice, lattice.getBoundingBox(), Array<T, 2>(u, (T)0.));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1., Array<T, 2>(u, (T)0.));
    
    // Instead of plain BounceBack, use the dynamics MomentumExchangeBounceBack,
    //   to compute the momentum exchange, and thus, the drag and the lift on
    //   the obstacle, locally during collision.

    /////////// Code with two embedded vectors ///////

    // Define the mask 
    MultiScalarField2D<bool> boolMaskPLB(nx, ny);
    plb_ifstream ifile("geometry.dat");
    ifile >> boolMaskPLB;

    // Convert to native 2D mask
    std::vector<std::vector<bool>> boolMask(ny, std::vector<bool>(nx));
    for (plint y = 0; y < ny; ++y)
        for (plint x = 0; x < nx; ++x)
            boolMask[y][x] = boolMaskPLB.get(x, y);

    //Define the dynamics and initialize momentum exchange with the class defined above
    defineDynamics(
       lattice, lattice.getBoundingBox(), new MaskShapeDomain2D<T>(boolMask),
       new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds));
    initializeMomentumExchange(
       lattice, lattice.getBoundingBox(), new MaskShapeDomain2D<T>(boolMask));
    ////////////////////////////////////////////////
    
    lattice.initialize();
}

void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(
        createFileName("u", iter, 6), *computeVelocityNorm(lattice), imSize, imSize);
}

void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);

    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx / dt);
    vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
    vtkOut.writeData<float>(*computePressure(lattice), "pressure", dx / dt);
    vtkOut.writeData<float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);

}

void writeHDF5(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters, plint iter)
{
    
#ifdef HDF5
    ParallelXdmfDataWriter2D xdmfOut("NACA0012");
#endif
    
#ifdef HDF5
    xdmfOut.writeDataField<T>(*computeVelocity(lattice), "velocity");
#endif  


}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    // Defines the flow parameters
    IncomprFlowParam<T> parameters(
        (T)1e-2,  // uMax
        (T)1000.,  // Re
        150,       // N
        10.,       // lx
        4.        // ly
    );
    const T logT = (T)0.01;


#ifndef PLB_REGRESSION
    const T imSave = (T)0.5;
    const T vtkSave = (T)1.;
    const T maxT = (T)20.1;
#else
    const T maxT = (T)0.5;
#endif

    writeLogFile(parameters, "Flow parameters");
    
    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        parameters.getNx(), parameters.getNy(),
        new BackgroundDynamics(parameters.getOmega())); // This currently implements MRT dynamics
    lattice.initialize();


    // The drag and lift acting on the obstacle are computed with help of the
    // internal statistics object of the lattice. For this purpose, they
    // need to be registered first, as it is done in the following lines.
    Array<plint, 2> forceIds;
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
        createInterpBoundaryCondition2D<T, DESCRIPTOR>();
    // boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    defineGeometry(lattice, parameters, *boundaryCondition, forceIds);


    //plb_ofstream ofileVelocity("velocityProfiles.dat");
    plb_ofstream ofileDrag("dragProfiles.dat");
    plb_ofstream ofileLift("liftProfiles.dat");

    // Main loop over time iterations.
    for (plint iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {
#ifndef PLB_REGRESSION
        if (iT % parameters.nStep(imSave) == 0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
        }
#endif

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "step " << iT << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT * parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "; av energy=" << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho=" << getStoredAverageDensity<T>(lattice)
                  << "; drag=" << lattice.getInternalStatistics().getSum(forceIds[0]) 
                  << "; lift=" << lattice.getInternalStatistics().getSum(forceIds[1]) << endl;
            
                  
            // Computation of lift and drag 
            ofileDrag << setprecision(10) << lattice.getInternalStatistics().getSum(forceIds[0]) << endl;
            ofileLift << setprecision(10) << lattice.getInternalStatistics().getSum(forceIds[1]) << endl;
            //pcout << "Saving HDF5 file ..." << endl;
            //writeHDF5(lattice, parameters, iT);
            

        }
        if (iT % parameters.nStep(vtkSave) == 0 && iT > 0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }
    }

    delete boundaryCondition;
}


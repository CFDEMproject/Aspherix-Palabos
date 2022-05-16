/*
  simulations of the experiments performed in

  Ten Cate, A., et al. "Particle imaging velocimetry experiments and
  lattice-Boltzmann simulations on a single sphere settling under gravity."
  Physics of Fluids (1994-present) 14.11 (2002): 4012-4025.

  parameters for the experiments in SI units : (rho_f, mu_f, u_inf)
  E1: 970 ; 0.373 ; 0.038
  E2: 965 ; 0.212 ; 0.06
  E3: 962 ; 0.113 ; 0.091
  E4: 960 ; 0.058 ; 0.128

 */

#define LBDEM_USE_WEIGHING

#include "palabos3D.h"
#include "palabos3D.hh"

#include "ibCompositeDynamics3D.h"
#include "ibDataWritingFunctionals3D.h"
#include "ibProcessors3D.h"
#include "physunits.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "aspherixSocketWrapper.h"

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::D3Q19Descriptor
#define DYNAMICS IBcompositeDynamics<T,DESCRIPTOR>(new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()))

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters,
              PhysUnits3D<T> const& units, plint iter)
{

  MultiScalarField3D<T> tmp(lattice.getNx(),lattice.getNy(),lattice.getNz());

  T p_fact = units.getPhysForce(1)/pow(units.getPhysLength(1),2)/3.;

  std::string fname(createFileName("vtk", iter, 6));

  VtkImageOutput3D<T> vtkOut(fname, units.getPhysLength(1));
  vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", units.getPhysVel(1));

  MultiScalarField3D<T> p(*computeDensity(lattice));
  subtractInPlace(p,1.);
  vtkOut.writeData<float>(p,"pressure",p_fact );

  pcout << "wrote " << fname << std::endl;
}

void writeGif(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
  const plint imSize = 600;

  Box3D slice(0,lattice.getNx(),
              lattice.getNy()/2,lattice.getNy()/2,
              0,lattice.getNz());

  MultiScalarField3D<T> vel(*computeVelocityNorm(lattice,slice));

  ImageWriter<T> imageWriter("leeloo");
  imageWriter.writeGif(createFileName("uNorm", iter, 6),
                       vel,
                       0., 0.02,
                       imSize, imSize );

}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);


    plint N(0);
    T uMax(0.),rho_f(0.),mu_f(0.),v_inf(0.), maxT(0.);
    std::string outDir;

    try {
        global::argv(1).read(N);
        global::argv(2).read(uMax);
        global::argv(3).read(rho_f);
        global::argv(4).read(mu_f);
        global::argv(5).read(v_inf);
        global::argv(6).read(maxT);
        global::argv(7).read(outDir);
    } catch(PlbIOException& exception) {
        pcout << "Error the parameters are wrong. The structure must be :\n";
        pcout << "1 : grid points along particle diameter\n";
        pcout << "2 : uMax\n";
        pcout << "3 : rho_f\n";
        pcout << "4 : mu_f\n";
        pcout << "5 : expected v_inf\n";
        pcout << "6 : maximal run time\n";
        pcout << "7 : outDir\n";
        exit(1);
    }


    std::string lbOutDir(outDir), demOutDir(outDir);
    lbOutDir.append("tmp/"); demOutDir.append("post/");
    global::directories().setOutputDir(lbOutDir);


    const T nu_f = mu_f/rho_f;
    const T lx = 0.1, ly = 0.1, lz = 0.16;

    T r_ = 0.015/2.;

    PhysUnits3D<T> units(2.*r_,v_inf,nu_f,lx,ly,lz,N,uMax,rho_f);
    units.setLbOffset(0.,0.,0.);


    IncomprFlowParam<T> parameters(units.getLbParam());

    const T vtkT = 0.002;
    const T logT = 0.02;

    const plint maxSteps = units.getLbSteps(maxT);
    const plint vtkSteps = max<plint>(units.getLbSteps(vtkT),1);
    const plint logSteps = max<plint>(units.getLbSteps(logT),1);

    writeLogFile(parameters, "3D sedimenting sphere");

    pcout << "-----------------------------------\n";
    pcout << "grid size: "
          << parameters.getNx() << " "
          << parameters.getNy() << " "
          << parameters.getNz() << std::endl;
    pcout << "-----------------------------------" << std::endl;

    pcout << "-----------------------------------\n";
    pcout << "props: rho_f\n"
          << rho_f << std::endl;
    pcout << "-----------------------------------" << std::endl;

    writeLogFile(parameters, "3D diagonal cavity");

    T omega = parameters.getOmega();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(omega) );

//=====
//BCs
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        //= createInterpBoundaryCondition3D<T,DESCRIPTOR>();
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), boundary::dirichlet);

    T u = std::sqrt((T)2)/(T)2 * units.getLbVel(0.); //* units.getLbVel(v_inf);
    initializeAtEquilibrium(lattice, everythingButTopLid, (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T) 1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

//=====



    defineDynamics(lattice,lattice.getBoundingBox(),new DYNAMICS);

    lattice.periodicity().toggleAll(false);


    T dt_phys = units.getPhysTime(1);
    pcout << "omega: " << parameters.getOmega() << "\n"
          << "dt_phys: " << dt_phys << "\n"
          << "Re : " << parameters.getRe() << "\n"
          << "vtkSteps: " << vtkSteps << "\n"
          << "grid size: "
          << parameters.getNx() << " "
          << parameters.getNy() << " "
          << parameters.getNz() << std::endl;

    lattice.initialize();


    clock_t start = clock();

    bool const useHyperthreading = true;
    AspherixSocketWrapper asx(global::mpi().getRank(),useHyperthreading,nu_f,rho_f);
    asx.initComm();

    double DEMts = asx.getDEMts();
    int nCGs = asx.getNumCG();
    int cg = asx.getCG()[0];

    // ==================
    // initialize
    // build vectors of push and pull names/types
    //setPushPullProperties();

    // define push (from DEM to CFD) properties
    asx.addPushProperty("radius",AspherixSocketWrapper::PropertyType::SCALAR_ATOM);
    asx.addPushProperty("x",AspherixSocketWrapper::PropertyType::VECTOR_ATOM);
    asx.addPushProperty("v",AspherixSocketWrapper::PropertyType::VECTOR_ATOM);

    asx.addPullProperty("dragforce",AspherixSocketWrapper::PropertyType::VECTOR_ATOM);

    asx.createProperties();
    asx.setParticleShapeType("sphere");

    // recv ok on comm setup fom DEM
    asx.confirmComm();

    int const maxStepsTmp = 3;
    // Loop over main time iteration.
    for (plint iT=0; iT<=maxStepsTmp; ++iT) {

        ParticleData<T>::ParticleDataArrayVector x_lb, v_lb;
        ParticleData<T>::ParticleDataScalarVector r_lb;

        bool const isLastExchange = iT == maxStepsTmp;
        asx.beginExchange(isLastExchange);

        // handle BB
        double limits[6] = {-100,100,-100,100,-100,100};
        asx.exchangeDomain(limits);

        // this relies on the fact that there is exactly one block on each lattice
        plint iBlock = lattice.getLocalInfo().getBlocks()[0];
        std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
        Box3D localBB = blockmap[iBlock];

        double r=0;
        double x[3] = {0.,0.,0.}, v[3] = {0.,0.,0.}, omega[3] = {0.,0.,0.};

        int nP_(0);
        asx.receiveData(nP_);

        while(asx.getNextParticleData(r,x,v))
        {
            r = units.getLbLength(r);
            r_lb.push_back(r);
            for(int i=0;i<3;i++)
            {
                x[i] = units.getLbPosition(x[i]);
                x_lb.push_back(Array<T,3>(x[0],x[1],x[2]));
                v[i] = units.getLbVel(v[i]);
                v_lb.push_back(Array<T,3>(v[0],v[1],v[2]));
            }

            pcout << "received "
                  << "r " << r << " | "
                  << "x " << x[0] << " " << x[1] << " " << x[2] << " | "
                  << "v " << v[0] << " " << v[1] << " " << v[2] << " | "
                  << std::endl;
            // set other properties (currently not communicated)
            plint const id = 1;  // TODO - we would need particle ID to have more than 1 particle in sim

            SetSingleSphere3D<T,DESCRIPTOR> *sss = 0;
            sss = new SetSingleSphere3D<T,DESCRIPTOR>(x,v,omega,x,r,id,false);
            Box3D sss_box = sss->getBoundingBox();

            // only go over part that lies on local processor
            // to avoid unnecessary communication overhead
            Box3D sss_box_intersect(0,0,0,0,0,0);

            bool boxes_intersect = intersect(sss_box,localBB,sss_box_intersect);

            if(boxes_intersect)
                applyProcessingFunctional(sss,sss_box_intersect,lattice);
            else
                delete sss;
        }

        // this one returns modif::staticVariables and forces an update of those along processor
        // boundaries
        applyProcessingFunctional(new AttributeFunctional<T,DESCRIPTOR>(),lattice.getBoundingBox(),lattice);

        //============================

        if( iT%vtkSteps == 0){
            writeVTK(lattice,parameters,units,iT);
            // writeGif(lattice,iT);
        }

        lattice.collideAndStream();

        plint const n_force = nP_*3;
        std::vector<T> force(n_force),torque(n_force);

        if(nP_ > 0)
        {
            SumForceTorque3D<T,DESCRIPTOR> *sft = new SumForceTorque3D<T,DESCRIPTOR>(x_lb,
                                                                                     &force.front(),&torque.front()
                );

            // this relies on the fact that there is exactly one block on each processor
            plint iBlock = lattice.getLocalInfo().getBlocks()[0];
            std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
            Box3D localBB = blockmap[iBlock];
            applyProcessingFunctional(sft,localBB, lattice);
        }

        //=========
        // gather local data send_data
        //gatherData(send_data,nP_); // TODO HERE
        for(int i=0;i<nP_;i++)
        {
            double f[3] = {
                units.getPhysForce(force[3*i+0]),
                units.getPhysForce(force[3*i+1]),
                units.getPhysForce(force[3*i+2])
            };

            asx.addNextSendParticle(f);
        }
        //=========

        asx.sendData();

        if(iT%logSteps == 0){
            clock_t end = clock();
            T time = ((T)difftime(end,start))/((T)CLOCKS_PER_SEC);
            T mlups = ((T) (lattice.getNx()*lattice.getNy()*lattice.getNz()*((T)logSteps)))/time/1e6;
            pcout << "time: " << time << " " ;
            pcout << "calculating at " << mlups << " MLU/s" << std::endl;
            start = clock();

        }


    }

}

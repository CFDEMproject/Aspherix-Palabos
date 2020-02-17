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

//#include "plb_ib.h"

#ifndef PLB_IB_H
#define PLB_IB_H

#include "ibCompositeDynamics3D.h"
#include "ibDataWritingFunctionals3D.h"
#include "ibProcessors3D.h"
#include "physunits.h"

#endif // PLB_IB_H


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include "aspherix_cosim_socket.h"

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


    const T g = 9.81;
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

    AspherixCoSimSocket sock_(false,global::mpi().getRank());

    pcout << "receiving DEMts ..." << std::endl;
    double DEMts;
    sock_.read_socket(&DEMts, sizeof(double));
    pcout << "receiving DEMts done." << std::endl;

    pcout << "sending couple_nevery ..." << std::endl;
    int couple_nevery(1);
    sock_.write_socket(&couple_nevery, sizeof(int));
    pcout << "sending couple_nevery done." << std::endl;

    pcout << "receiving nr nCGs ..." << std::endl;
    int nCGs;
    sock_.read_socket(&nCGs, sizeof(int));
    pcout << "receiving nr nCGs done. nCGs=" << nCGs << std::endl;

    pcout << "receiving cg ..." << std::endl;
    int cg;
    sock_.read_socket(&cg, sizeof(int));
    pcout << "receiving cg done. cg=" << cg << std::endl;

    // ==================
    // initialize
    // build vectors of push and pull names/types
    //setPushPullProperties();

    // define push (from DEM to CFD) properties
    sock_.pushBack_pushNameList(std::string("radius"));
    sock_.pushBack_pushTypeList(std::string("scalar-atom"));
    sock_.pushBack_pushNameList(std::string("x"));
    sock_.pushBack_pushTypeList(std::string("vector-atom"));
    sock_.pushBack_pushNameList(std::string("v"));
    sock_.pushBack_pushTypeList(std::string("vector-atom"));

    // define pull (from CFD to DEM) properties
    //sock_.pushBack_pullNameList(std::string("Ksl"));
    //sock_.pushBack_pullTypeList(std::string("scalar-atom"));
    //sock_.pushBack_pullNameList(std::string("uf"));
    //sock_.pushBack_pullTypeList(std::string("vector-atom"));
    sock_.pushBack_pullNameList(std::string("dragforce"));
    sock_.pushBack_pullTypeList(std::string("vector-atom"));

    // build byte pattern of the push (from DEM to CFD) & pull (from CFD to DEM) properties
    sock_.buildBytePattern();

    // send push (from DEM to CFD) & pull (from CFD to DEM) property names and types
    pcout << "send push/pull ..." << std::endl;
    sock_.sendPushPullProperties();
    pcout << "send push/pull done." << std::endl;

    // send particleShapeType to DEM
    string shapeTypeName("sphere");
    size_t send_dataSize = shapeTypeName.size()+1;
    char *send_data = new char[send_dataSize];
    strcpy(send_data,const_cast<char*>(shapeTypeName.c_str()));
    sock_.sendData(send_dataSize,send_data);
    // ==================


    // Loop over main time iteration.
    for (plint iT=0; iT<maxSteps; ++iT) {
    //for (plint iT=0; iT<10; ++iT) {

        //===
        // exchange status
        // TODO here is the point where we can tell DEM to stop by sending SocketCodes::request_quit
        pcout << "send SocketCodes..." << std::endl;
        sock_.exchangeStatus(::SocketCodes::start_exchange,::SocketCodes::start_exchange);
        pcout << "send SocketCodes done." << std::endl;

        // send flag to proceed
        pcout << "sending info on command models ..." << std::endl;
        ::SocketCodes hh=::SocketCodes::start_exchange;
        sock_.write_socket(&hh, sizeof(SocketCodes));
        pcout << "sending info on command models done." << std::endl;

        // handle BB
        pcout << "send BB..." << std::endl;
        bool useBB(true);
        double limits[6];
        //setBounds(limits,useBB);
        limits[0]=-100;
        limits[1]=100;
        limits[2]=-100;
        limits[3]=100;
        limits[4]=-100;
        limits[5]=100;
        sock_.exchangeDomain(useBB,limits);
        pcout << "send BB done." << std::endl;

        // init rcv_data pointer
        size_t rcv_dataSize(0);
        char *rcv_data=nullptr;

        // receive rcv_data from DEM to CFD
        pcout << "receive data..." << std::endl;
        sock_.rcvData(rcv_dataSize,rcv_data);
        pcout << "receive data done." << std::endl;

        // set nr of particles
        int nP_ = rcv_dataSize/sock_.get_rcvBytesPerParticle();

        //============================
        // distribute received rcv_data to local structure
        //distributeData(rcv_data,nP_);  // TODO HERE

        // this relies on the fact that there is exactly one block on each lattice
        plint iBlock = lattice.getLocalInfo().getBlocks()[0];
        std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
        Box3D localBB = blockmap[iBlock];

        // loop all particles
        int index_from(0);
        int len(0);
        int h=sock_.get_rcvBytesPerParticle();
        int h2=sock_.get_pushBytesPerPropList().size();
        for(int i=0;i<nP_;i++)
        {
            T r; r=0;
            T x[3],v[3],omega[3];
            for(int k=0;k<3;k++) x[k]=v[k]=omega[k]=0;
            plint id;

            // loop all properties
            for (int j = 0; j < h2; j++)
            {
                index_from = i*h + sock_.get_pushCumOffsetPerProperty()[j];
                len = sock_.get_pushBytesPerPropList()[j];

                // cut out part of string
                char* str = new char[len+1];
                memmove(str, rcv_data + index_from, len);
                str[len] = '\0';

                if(sock_.get_pushTypeList()[j]=="scalar-atom")
                {
                    double* b=(double*)(str);
                    pcout << "receiving " << sock_.get_pushNameList()[j] << " = " << b[0] << endl;

                    int fieldID(-1);
                    if(sock_.get_pushNameList()[j]=="radius")
                    {
                        r = units.getLbLength(b[0]);
                    }
                    else pcout<<"ERROR: unknown fieldID!!! pushNameList[j]=\n" << sock_.get_pushNameList()[j] << " of type" << sock_.get_pushTypeList()[j] << endl;
                }
                else if(sock_.get_pushTypeList()[j]=="vector-atom")
                {
                    double* b=(double*)(str);
                    pcout << "receiving " << sock_.get_pushNameList()[j] << " = " << b[0] << "," << b[1] << "," << b[2] << "." << endl;

                    int fieldID(-1);
                    if(sock_.get_pushNameList()[j]=="v")
                    {
                        for(int k=0;k<3;k++)
                        {
                            v[k] = units.getLbVel(b[k]);
                        }
                    }
                    else if(sock_.get_pushNameList()[j]=="x")
                    {
                        for(int k=0;k<3;k++)
                        {
                            x[k] = units.getLbPosition(b[k]);
                        }
                    }
                    else pcout<<"unknown fieldID!!! pushNameList[j]=\n" << sock_.get_pushNameList()[j] << " of type" << sock_.get_pushTypeList()[j] << endl;
                }
                delete [] str;
            }

            // set other properties (currently not communicated)
            id = 1;  // TODO - we would need particle ID to have more than 1 particle in sim
            for(plint i=0;i<3;i++)  omega[i] = units.getLbFreq(0.);

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

        size_t send_dataSize = nP_*sock_.get_sndBytesPerParticle();
        char *send_data = new char[send_dataSize];

        //getForcesFromLattice(lattice,wrapper,sock_,units);  //we want to have a verion of this w/o wrapper
        static std::vector<T> force,torque;
        {
            static typename ParticleData<T>::ParticleDataArrayVector x_lb;

            plint const nPart = nP_;//wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
            plint const n_force = nPart*3;


            if(nPart > 0)
            {
                //===
                // TODO - can we move this to where pos is available?
                if(nPart > x_lb.size()){
                  for(plint iPart=0;iPart<x_lb.size();iPart++){
                    x_lb[iPart][0] = units.getLbPosition(0.05);
                    x_lb[iPart][1] = units.getLbPosition(0.05);
                    x_lb[iPart][2] = units.getLbPosition(0.08);
                  }
                  for(plint iPart = x_lb.size();iPart < nPart; iPart++)
                    x_lb.push_back( Array<T,3>( units.getLbPosition(0.05),
                                                units.getLbPosition(0.05),
                                                units.getLbPosition(0.08) ) );
                } else{
                  for(plint iPart=0;iPart<nPart;iPart++){
                    x_lb[iPart][0] = units.getLbPosition(0.05);
                    x_lb[iPart][1] = units.getLbPosition(0.05);
                    x_lb[iPart][2] = units.getLbPosition(0.08);
                  }
                }
                //===

                if(n_force > force.size()){
                  for(plint i=0;i<force.size();i++){
                    force[i] = 0;
                    torque[i] = 0;
                  }
                  for(plint i=force.size();i<n_force;i++){
                    force.push_back(0.);
                    torque.push_back(0.);
                  }
                } else {
                  for(plint i=0;i<n_force;i++){
                    force[i] = 0;
                    torque[i] = 0;
                  }
                }
                SumForceTorque3D<T,DESCRIPTOR> *sft = new SumForceTorque3D<T,DESCRIPTOR>(x_lb,
                                                                                         &force.front(),&torque.front()
                                                                                         );

                // this relies on the fact that there is exactly one block on each processor
                plint iBlock = lattice.getLocalInfo().getBlocks()[0];
                std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
                Box3D localBB = blockmap[iBlock];
                applyProcessingFunctional(sft,localBB, lattice);
            }

        }

        //=========
        // gather local data send_data
        //gatherData(send_data,nP_); // TODO HERE
        h=sock_.get_sndBytesPerParticle();
        h2=sock_.get_pullNameList().size();
        for(int i=0;i<nP_;i++)
        {
            for (int j = 0; j < h2; j++)
            {
                char *ptr = nullptr;

                int index_from = i*h + sock_.get_pullCumOffsetPerProperty()[j];
                int len = sock_.get_pullBytesPerPropList()[j];

                if(sock_.get_pullTypeList()[j]=="scalar-atom")
                {
                    double* v=new double[1];
                    int fieldID(-1);
                    if(sock_.get_pullNameList()[j]=="Ksl")
                    {
                        pcout<<"ERROR: NOT IMPLEMENTED" << endl;
                    }
                    else pcout<<"ERROR: unknown fieldID!!! pullNameList[j]=\n" << sock_.get_pullNameList()[j] << " of type" << sock_.get_pullTypeList()[j] << endl;

                    ptr = (char*)v;
                    memcpy(&send_data[index_from], ptr, len);
                    delete [] v;
                }
                else if(sock_.get_pullTypeList()[j]=="vector-atom")
                {
                    double* v=new double[3];
                    int fieldID(-1);
                    if(sock_.get_pullNameList()[j]=="dragforce")
                    {
                        for(int k=0;k<3;k++) v[k]=units.getPhysForce(force[3*i+k]); // TODO: we do not yet receive a force here!?
                        pcout << "sending " << sock_.get_pullNameList()[j] << " = " << v[0] << "," << v[1] << "," << v[2] << "." << endl;
                    }
                    else if(sock_.get_pullNameList()[j]=="uf")
                    {
                        pcout<<"ERROR: NOT IMPLEMENTED" << endl;
                    }
                    else pcout<<"ERROR: unknown fieldID!!! pullNameList[j]=\n" << sock_.get_pullNameList()[j] << " of type" << sock_.get_pullTypeList()[j] << endl;

                    ptr = (char*)v;
                    memcpy(&send_data[index_from], ptr, len);
                    delete [] v;
                }
            }
        }
        //=========

        // send snd_data to DEM from CFD
        sock_.sendData(send_dataSize,send_data);
        pcout << "sock_.sendData(send_dataSize,send_data) done." << std::endl;

        delete[] send_data;

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

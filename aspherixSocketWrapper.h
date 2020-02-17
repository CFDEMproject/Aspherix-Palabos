#ifndef ASPHERIX_SOCKET_WRAPPER_H
#define ASPHERIX_SOCKET_WRAPPER_H

#include "aspherix_cosim_socket.h"


class AspherixSocketWrapper {
public:
    AspherixSocketWrapper(int const nRank)
        : sock_(false,nRank),
          demTS(0.),
          nCGs(0),
          cg(nullptr)
    {}

    void initComm()
    {
        sock_.read_socket(&DEMts, sizeof(double));
        int const couple_nevery = 1; // hardcoded for the moment
        sock_.write_socket(&couple_nevery, sizeof(int));
        sock_.read_socket(&nCGs, sizeof(int));

        if(cg) delete[] cg;
        cg = new int[nCGs];

        sock_.read_socket(&cg, sizeof(int)*nCGs);
    }


    double getDEMts() const { return demTS; }
    int getNumCG() const { return nCGs; }
    int const * const getCG() const { return cg; }

private:
    AspherixCoSimSocket sock_;

    double demTS;
    int nCGs;
    int  *cg;

};


#endif

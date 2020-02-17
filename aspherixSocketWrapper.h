#ifndef ASPHERIX_SOCKET_WRAPPER_H
#define ASPHERIX_SOCKET_WRAPPER_H

#include <map>
#include <string>
#include <vector>

#include "aspherix_cosim_socket.h"

class AspherixSocketWrapper {
public:
    enum class PropertyType { SCALAR_ATOM, VECTOR_ATOM };

    AspherixSocketWrapper(int const nRank)
        : sock_(false,nRank),
          demTS(0.),
          nCGs(0),
          cg(nullptr)
    {}

    void initComm()
    {
        sock_.read_socket(&demTS, sizeof(double));
        int couple_nevery = 1; // hardcoded for the moment
        sock_.write_socket(&couple_nevery, sizeof(int));
        sock_.read_socket(&nCGs, sizeof(int));

        if(cg) delete[] cg;
        cg = new int[nCGs];

        sock_.read_socket(&cg, sizeof(int)*nCGs);
    }

    void addPushProperty(std::string const &name, PropertyType const type)
    {
        sock_.pushBack_pushNameList(name);
        sock_.pushBack_pushTypeList(propTypeToStr(type));
    }

    void addPullProperty(std::string const &name, PropertyType const type)
    {
        sock_.pushBack_pullNameList(name);
        sock_.pushBack_pullTypeList(propTypeToStr(type));
    }

    void createProperties()
    {
        sock_.buildBytePattern();
        sock_.sendPushPullProperties();
    }

    void setParticleShapeType(std::string const &pst)
    {
        size_t sendSize = pst.size()+1;
        char *c = new char[sendSize];
        strcpy(c,pst.c_str());
        sock_.sendData(sendSize,c);
        delete[] c;
    }

    void beginExchange()
    {
        sock_.exchangeStatus(::SocketCodes::start_exchange,::SocketCodes::start_exchange);
        ::SocketCodes hh=::SocketCodes::start_exchange;
        sock_.write_socket(&hh, sizeof(SocketCodes));
    }

    void exchangeDomain(double *limits)
    {
        bool const useBB = true;
        sock_.exchangeDomain(useBB,limits);
    }

    void exchangeData()

    double getDEMts() const { return demTS; }
    int getNumCG() const { return nCGs; }
    int const * const getCG() const { return cg; }

    AspherixCoSimSocket sock_;

private:
    double demTS;
    int nCGs;
    int  *cg;

    std::string propTypeToStr(PropertyType const type) const
    {
        switch(type)
        {
        case PropertyType::SCALAR_ATOM:
            return "scalar-atom";
        case PropertyType::VECTOR_ATOM:
            return "vector-atom";
        }
    }
};


#endif

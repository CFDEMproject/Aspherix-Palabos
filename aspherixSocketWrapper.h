#ifndef ASPHERIX_SOCKET_WRAPPER_H
#define ASPHERIX_SOCKET_WRAPPER_H

#include <map>
#include <string>
#include <vector>

#include "aspherix_cosim_socket.h"

class AspherixSocketWrapper {
public:
    enum class PropertyType { SCALAR_ATOM, VECTOR_ATOM };

    AspherixSocketWrapper(int const nRank, bool const _hyperthreading)
        : sock_(false,nRank),
          hyperthreading(_hyperthreading),
          demTS(0.),
          nCGs(0),
          cg(nullptr),
          send_dataSize(0),
          rcv_dataSize(0),
          send_data(nullptr),
          rcv_data(nullptr),
          nPart(0),
          iPart_send(0),
          iPart_rcv(0),
          firstStep(true)
    {}

    void initComm()
    {
        std::string commProtocolVersion_("0c84OrF9\0");
        size_t nSend = commProtocolVersion_.size()+1;
        char *versionSend = new char[nSend];
        strcpy(versionSend, const_cast<char*>(commProtocolVersion_.c_str()));

        sock_.sendData(nSend,versionSend);
        delete[] versionSend;

        bool socketOK = false;
        sock_.read_socket(&socketOK,sizeof(bool));

        if(!socketOK)
        {
            std::cerr << "socket library versions don't match" << std::endl;
            exit(1);
        }

        sock_.read_socket(&demTS, sizeof(double));
        int couple_nevery = 1; // hardcoded for the moment
        sock_.write_socket(&couple_nevery, sizeof(int));

        sock_.write_socket(&hyperthreading,sizeof(bool));

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
        if(!hyperthreading && !firstStep)
        {
            // start DEM
            sock_.exchangeStatus(::SocketCodes::ping,::SocketCodes::ping);
            // wait for DEM to finish
            sock_.exchangeStatus(::SocketCodes::ping,::SocketCodes::ping);
        }
        firstStep = false;

        sock_.exchangeStatus(::SocketCodes::start_exchange,::SocketCodes::start_exchange);

        ::SocketCodes hh=::SocketCodes::start_exchange;
        sock_.write_socket(&hh, sizeof(SocketCodes));
    }

    void exchangeDomain(double *limits)
    {
        bool const useBB = true;
        sock_.exchangeDomain(useBB,limits);
    }

    void receiveData(int &nP)
    {
        rcv_dataSize = 0;
        rcv_data = nullptr;

        sock_.rcvData(rcv_dataSize,rcv_data);

        nPart = rcv_dataSize/sock_.get_rcvBytesPerParticle();
        nP = nPart;
        iPart_rcv = 0;
    }

    bool getNextParticleData(double &r, double x[3], double v[3])
    {
        if(iPart_rcv >= nPart)
            return false;

        int h=sock_.get_rcvBytesPerParticle();
        int h2=sock_.get_pushBytesPerPropList().size();
        // loop all properties
        for (int j = 0; j < h2; j++)
        {
            int const index_from = iPart_rcv*h + sock_.get_pushCumOffsetPerProperty()[j];
            int const len = sock_.get_pushBytesPerPropList()[j];

            // cut out part of string
            char* str = new char[len+1];
            memmove(str, rcv_data + index_from, len);
            str[len] = '\0';

            if(sock_.get_pushTypeList()[j]=="scalar-atom")
            {
                double* b=(double*)(str);

                if(sock_.get_pushNameList()[j]=="radius")
                {
                    r = b[0];
                }
            }
            else if(sock_.get_pushTypeList()[j]=="vector-atom")
            {
                double* b=(double*)(str);

                if(sock_.get_pushNameList()[j]=="v")
                {
                    for(int k=0;k<3;k++)
                    {
                        v[k] = b[k];
                    }
                }
                else if(sock_.get_pushNameList()[j]=="x")
                {
                    for(int k=0;k<3;k++)
                    {
                        x[k] = b[k];
                    }
                }
            }
            delete [] str;

        }
        iPart_rcv++;
        return true;
    }

    void addNextSendParticle(double *f)
    {
        if(iPart_send == 0)
        {
            send_dataSize = nPart*sock_.get_sndBytesPerParticle();
            send_data = new char[send_dataSize];
        }

        // we assume for the moment that only the drag force is being transfered
        // can be easily extended to account for torque as well
        int const h = sock_.get_sndBytesPerParticle();
        int const h2 = sock_.get_pullNameList().size();

        for(int j=0;j<h2;j++)
        {
            int const index_from = iPart_send*h + sock_.get_pullCumOffsetPerProperty()[j];
            int const len = sock_.get_pullBytesPerPropList()[j];
            if(sock_.get_pullTypeList()[j]=="vector-atom" && sock_.get_pullNameList()[j]=="dragforce")
            {
                memcpy(&send_data[index_from],(char*)f,len);
            }
        }

        iPart_send++;
    }

    void sendData()
    {
        sock_.sendData(send_dataSize,send_data);
        delete[] send_data;
        send_data = nullptr;
        iPart_send = 0;
    }

    double getDEMts() const { return demTS; }
    int getNumCG() const { return nCGs; }
    int const * const getCG() const { return cg; }

    AspherixCoSimSocket sock_;

private:
    bool hyperthreading;
    double demTS;
    int nCGs;
    int  *cg;

    // init rcv_data pointer
    size_t send_dataSize, rcv_dataSize;
    char *send_data, *rcv_data;

    int nPart, iPart_send, iPart_rcv;
    bool firstStep;

    std::string propTypeToStr(PropertyType const type) const
    {
        switch(type)
        {
        case PropertyType::SCALAR_ATOM:
            return "scalar-atom";
        case PropertyType::VECTOR_ATOM:
            return "vector-atom";
        }
        return "";
    }
};


#endif

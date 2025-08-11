#ifndef ASPHERIX_SOCKET_WRAPPER_H
#define ASPHERIX_SOCKET_WRAPPER_H

#include <array>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "aspherix_cosim_socket.h"
#include "aspherix_cosim_field.h"
#include "aspherix_cosim_settings.h"

constexpr auto kRecv = CoSimSocket::SyncDirection::kServerToClient;
constexpr auto kSend = CoSimSocket::SyncDirection::kClientToServer;

class AspherixSocketWrapper {
public:
    enum class PropertyType { SCALAR_ATOM, VECTOR_ATOM };

    AspherixSocketWrapper(int const nRank, bool const _hyperthreading,
                          double const _cfd_viscosity, double const _cfd_density)
        : sock_(std::make_shared<CoSimSocket::AspherixCoSimSocket>(CoSimSocket::Mode::kClient, nRank)),
          hyperthreading(_hyperthreading),
          demTS(0.),
          cfd_viscosity(_cfd_viscosity),
          cfd_density(_cfd_density),
          iPart_send(0),
          iPart_rcv(0),
          firstStep(true)
    {}

    void processFields()
    {
        pull_size_one_point_ = 0;
        push_size_one_point_ = 0;
        for (auto &field : fields_)
        {
            // const auto field_data_length = field.data_length();
            // if (needs_reallocation_)
            //     >allocateFieldData(field_index, 0.0, field_data_length, num_points_);
            if (field.isClientToServer())
            {
                field.setOffset(push_size_one_point_);
                push_size_one_point_ += field.dataTypeSize();
            }
            else if (field.isServerToClient())
            {
                field.setOffset(pull_size_one_point_);
                pull_size_one_point_ += field.dataTypeSize();
            }
        }
    }

    void initializeCommunication()
    {
        using namespace CoSimSocket;

        const std::string commProtocolVersion_ = "HaikuCup";
        auto settings = CoSimSettings(sock_);
        settings.addSetting("client_protocol_id", commProtocolVersion_, kSend);
        settings.sync();

        bool socketOK = false;
        settings.clearSettings();
        settings.addSetting("matching_protocol_ids", socketOK, kRecv);
        settings.sync();
        socketOK = settings.getSetting<bool>("matching_protocol_ids");

        if(!socketOK)
        {
            std::cerr << "socket library versions don't match" << std::endl;
            exit(1);
        }

        // Receive data from Aspherix
        using array3 = std::array<double, 3>;
        array3 gravity_DEM;
        int n_particle_templates = 0;

        settings.clearSettings();
        settings.addSetting("DEMts", demTS, kRecv);
        settings.addSetting("gravity", gravity_DEM, kRecv);
        settings.addSetting("num_particle_templates", n_particle_templates, kRecv);
        settings.sync();

        demTS = settings.getSetting<double>("DEMts");
        gravity_DEM = settings.getSetting<array3>("gravity");
        n_particle_templates = settings.getSetting<int>("num_particle_templates");

        // Send data to Aspherix (hardcoded for the moment)
        int couple_nevery = 1;
        double start_time = 0.;
        bool const_density = true;
        bool trigger_DEM_output = false;          // write independently
        double write_interval = 0.;               // write interval -- unused but for compliance
        std::vector<std::string> boundary_names;  // unused

        settings.clearSettings();
        settings.addSetting("couple_nevery", couple_nevery, kSend);
        settings.addSetting("cfd_start_time", start_time, kSend);
        settings.addSetting("is_const_density", const_density, kSend);
        settings.addSetting("cfd_density", cfd_density, kSend);
        settings.addSetting("cfd_viscosity", cfd_viscosity, kSend);
        settings.addSetting("hyperthreading", hyperthreading, kSend);
        settings.addSetting("trigger_output", trigger_DEM_output, kSend);
        settings.addSetting("write_interval_CFD", write_interval, kSend);
        settings.addSetting("boundary_names", boundary_names, kSend);
        settings.sync();

        // check status against DEM
        if (hyperthreading)
            sock_->exchangeStatus(SocketCodes::kStartExchange, SocketCodes::kStartExchange);
        else
            sock_->exchangeStatus(SocketCodes::kPing, SocketCodes::kPing);

        [[maybe_unused]] const int num_coupled_boundaries = sock_->readValue<std::size_t>();
        cg = sock_->readData<double>();

        const std::string solverName("cfdemSolverPiso");
        settings.clearSettings();
        settings.addSetting("particle_shape_type", particle_shape_, kSend);
        settings.addSetting("fields", fields_, kSend);
        settings.addSetting("solver_name", solverName, kSend);
        settings.sync();

        auto solver_code = sock_->exchangeStatus(SocketCodes::kStartExchange);
        if (solver_code != SocketCodes::kStartExchange)
        {
            if (solver_code == SocketCodes::kInvalid)
                std::cerr << "The solver you are using ('" << solverName << "') is not "
                    << "compatible with your Aspherix version. "
                    << "Check the error message from your DEM code for more details."
                    << std::endl;
            else
                std::cerr << "Socket status codes dont match! This is fatal." << std::endl;
            exit(1);
        }

        // Check for errors on DEM side
        SocketCodes code = sock_->exchangeStatus();
        if (code != SocketCodes::kStartExchange)
        {
            if (code == SocketCodes::kInvalid)
                // there was an error on the DEM side, report
                std::cerr << "Not all requested properties are available on the DEM side. "
                    << "Check the error message from your DEM code for more details."
                    << std::endl;
            else
                // status codes dont match --> something else went wrong
                std::cerr << "Socket status codes dont match! This is fatal." << std::endl;
            exit(1);
        }

        sock_->exchangeStatus(SocketCodes::kStartExchange, SocketCodes::kPing);
    }

    void closeComm()
    {
        using CoSimSocket::SocketCodes;
        // if(hyperthreading)
        //     sock_->exchangeStatus(SocketCodes::kCloseConnection,SocketCodes::kStartExchange);
        // else
        //     sock_->exchangeStatus(SocketCodes::kCloseConnection,SocketCodes::ping);

        std::cerr << "************ CLOSING!!!!! ***************" << std::endl;

        sock_->exchangeStatus(SocketCodes::kRequestQuit,SocketCodes::kStartExchange);
    }

    void addField(std::string const &name, PropertyType const type, const CoSimSocket::SyncDirection direction)
    {
        using namespace CoSimSocket;
        auto [data_type, data_length, object] = propTypeToTuple(type);
        fields_.emplace_back(CoSimField(name, "particles", data_type, data_length, object, direction));
    }

    void setParticleShapeType(std::string const &pst)
    {
        particle_shape_ = pst;
    }

    void beginExchange(bool const isLastExchange)
    {
        using namespace CoSimSocket;

        if(!hyperthreading && !firstStep)
        {
            // start DEM
            sock_->exchangeStatus(SocketCodes::kPing,SocketCodes::kPing);
            // wait for DEM to finish
            sock_->exchangeStatus(SocketCodes::kPing,SocketCodes::kPing);
        }
        firstStep = false;

        // communicate simulation end
        SocketCodes code = isLastExchange ? SocketCodes::kRequestQuit : SocketCodes::kStartExchange;
        sock_->exchangeStatus(code, SocketCodes::kStartExchange);

        sock_->exchangeStatus(SocketCodes::kStartExchange, SocketCodes::kPing);
    }

    void exchangeDomain(std::array<double, 6> limits)
    {
        using namespace CoSimSocket;

        sock_->exchangeStatus(SocketCodes::kStartExchange, SocketCodes::kStartExchange);
        sock_->exchangeStatus(SocketCodes::kBoundingBoxUpdate, SocketCodes::kPing);
        sock_->exchangeValue<CoSimSocket::SyncDirection::kSend>(limits);
    }

    void receiveData()
    {
        using namespace CoSimSocket;

        const std::size_t n_containers = sock_->readValue<std::size_t>();
        auto server_code = sock_->exchangeStatus(SocketCodes::kStartExchange);
        if (server_code == SocketCodes::kPing)
        {
            num_points_ = sock_->readValue<std::size_t>();
            const auto size_one_point = sock_->readValue<std::size_t>();
            sock_->readData<char>(recv_buffer_);

            assert(size_one_point == pull_size_one_point_);

            send_buffer_.resize(num_points_ * push_size_one_point_);
            recv_buffer_.resize(num_points_ * pull_size_one_point_);
        }
        else
        {
            std::cerr << "Unexpected status code! This is fatal." << std::endl;
            exit(1);
        }
        iPart_rcv = 0;
    }

    bool getNextParticleData(double &r, double x[3], double v[3])
    {
        if(iPart_rcv >= num_points_)
            return false;

        for (const auto& field : fields_)
        {
            if (field.isClientToServer())
                continue;

            const int index_from = iPart_rcv * pull_size_one_point_ + field.offset();
            const void* ptr_to_value = recv_buffer_.data() + index_from;
            const double *value = reinterpret_cast<const double*>(ptr_to_value);

            if (field.name() == "radius")
                r = value[0];
            else
            {
                double *ptr_to_data = nullptr;
                if (field.name() == "x")
                    ptr_to_data = &x[0];
                else if (field.name() == "v")
                    ptr_to_data = &v[0];

                for (std::size_t k = 0; k < field.data_length(); ++k)
                {
                    ptr_to_data[k] = value[k];
                }
            }
        }
        iPart_rcv++;
        return true;
    }

    void addNextSendParticle(double *f)
    {
        // we assume for the moment that only the drag force is being transfered
        // can be easily extended to account for torque as well
        for (const auto& field : fields_)
        {
            if (field.isServerToClient())
                continue;

            const auto data_index_from = iPart_send * push_size_one_point_ + field.offset();
            if (field.name() == "dragforce")
                memcpy(&send_buffer_[data_index_from], f, field.dataTypeSize());
        }
        iPart_send++;
    }

    void sendData()
    {
        num_points_ = sock_->readValue<std::size_t>();
        const auto size_one_point = sock_->readValue<std::size_t>();
        sock_->writeData(send_buffer_);
        iPart_send = 0;
    }

    double getDEMts() const { return demTS; }
    int getNumCG() const { return static_cast<int>(cg.size()); }
    const std::vector<double> getCG() const { return cg; }
    std::size_t getNumParticles() const { return num_points_; }

private:
    std::shared_ptr<CoSimSocket::AspherixCoSimSocket> sock_;

    bool hyperthreading;
    double demTS, cfd_viscosity, cfd_density;
    std::vector<double> cg;
    std::string particle_shape_;

    std::vector<CoSimSocket::CoSimField> fields_;
    std::size_t pull_size_one_point_;
    std::size_t push_size_one_point_;

    std::size_t num_points_;

    std::vector<char> recv_buffer_;
    std::vector<char> send_buffer_;

    std::size_t iPart_send, iPart_rcv;
    bool firstStep;

    std::tuple<CoSimSocket::DataType, std::size_t, CoSimSocket::DataObject> propTypeToTuple(PropertyType const type) const
    {
        using namespace CoSimSocket;
        DataType data_type = DataType::kDouble;
        std::size_t data_length = 0;
        DataObject object = DataObject::kParticle;

        switch(type)
        {
        case PropertyType::SCALAR_ATOM:
            data_length = 1;
            break;
        case PropertyType::VECTOR_ATOM:
            data_length = 3;
            break;
        }
        return {data_type, data_length, object};
    }
};


#endif

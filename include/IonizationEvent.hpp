#ifndef IONIZATIONEVENT_HPP
#define IONIZATIONEVENT_HPP


#include "vector3.hpp"

// TODO: check that this class does the following
// stores charge as electron integer type if set is called with
// an intger argument
// stores charge as double coulomb type if set is called with
// a double argument
// stores charge as integer whenever possible unless overflow
// this isn't relevant unless we are adding etc
class ElectronicCharge
{

    enum class ElectronicChargeMode
    {
        DOUBLE,
        INTEGER
    };


    public:

    explicit
    ElectronicCharge(const double charge)
        : _mode_{ElectronicChargeMode::DOUBLE}
        , _coulomb_{charge}
    {
    }

    explicit
    ElectronicCharge(const uint64_t charge)
        : _mode_{ElectronicChargeMode::INTEGER}
        , _electron_{charge}
    {
    }

    ElectronicCharge& operator+=(const ElectronicCharge& other)
    {
        if(_mode_ == ElectronicChargeMode::INTEGER && other._mode_ == ElectronicChargeMode::INTEGER)
        {
            // compute the amount of free space
            uint64_t free_space{_uint64_t_max_ - _electron_};
            if(other._electron_ <= free_space)
            {
                _electron_ += other._electron_;
            }
            else
            {
                _coulomb_ = electron_to_coulomb(_electron_) + other.get_charge_coulomb();
                _mode_ = ElectronicChargeMode::DOUBLE;
            }
        }
        else if(_mode_ == ElectronicChargeMode::DOUBLE && other._mode_ == ElectronicChargeMode::DOUBLE)
        {
            _coulomb_ += other._coulomb_;
        }
        else if((_mode_ == ElectronicChargeMode::INTEGER && other._mode_ == ElectronicChargeMode::DOUBLE) || (_mode_ == ElectronicChargeMode::DOUBLE && other._mode_ == ElectronicChargeMode::INTEGER))
        {
            _coulomb_ = get_charge_coulomb() + other.get_charge_coulomb();
            _mode_ == ElectronicChargeMode::DOUBLE;
        }
        else
        {
            throw __func__;
        }
        return *this;
    }

    double GetChargeCoulomb() const
    {
        return get_charge_coulomb();
    }

    uint64_t GetChargeElectron() const
    {
        return get_charge_electron();
    }

    private:
    inline
    double get_charge_coulomb() const
    {
        if(_mode_ == ElectronicChargeMode::INTEGER)
        {
            return electron_to_coulomb(_electron_);
        }
        else if(_mode_ == ElectronicChargeMode::DOUBLE)
        {
            return _coulomb_;
        }
        else
        {
            throw __func__;
        }
    }

    inline
    uint64_t get_charge_electron() const
    {
        if(_mode_ == ElectronicChargeMode::INTEGER)
        {
            return _electron_;
        }
        else if(_mode_ == ElectronicChargeMode::DOUBLE)
        {
            if(check_overflow(_coulomb_))
            {
                throw __func__;
            }
            return coulomb_to_electron(_coulomb_);
        }
        else
        {   
            throw __func__;
        }

    }

    private:
    bool check_overflow(const double coulomb) const
    {
        return coulomb > _max_electron_count_uint64_as_double_;
    }

    
    // set coulomb from coulomb
    void set_charge_coulomb(const double charge)
    {
        _coulomb_ = charge;
        _mode_ = ElectronicChargeMode::DOUBLE;
    }
    
    // set coulomb from electron
    void set_charge_coulomb(const uint64_t charge)
    {
        _coulomb_ = electron_to_coulomb(charge); 
        _mode_ = ElectronicChargeMode::DOUBLE;
    }

    // set electron from electron
    void set_charge_electron(const uint64_t charge)
    {
        _electron_ = charge;
        _mode_ = ElectronicChargeMode::INTEGER;
    }

    // set electron from coulomb
    void set_charge_electron(const double charge)
    {
        _electron_ = coulomb_to_electron(charge);
        _mode_ = ElectronicChargeMode::INTEGER;
    }
    

    inline
    uint64_t coulomb_to_electron(const double charge) const
    {
        return std::round(charge / _electron_charge_fundamental_);
    }

    inline
    double electron_to_coulomb(const uint64_t charge) const
    {
        return _electron_charge_fundamental_ * (double)charge;
    }
    

    public:
    void SetStorageMode(const ElectronicChargeMode mode)
    {
        if(mode == _mode_)
        {
            return;
        }
        else if(mode == ElectronicChargeMode::INTEGER)
        {
            //SetChargeElectron(_coulomb_);
            if(check_overflow(_coulomb_))
            {
                set_charge_electron(_coulomb_);
            }
        }
        else if(mode == ElectronicChargeMode::DOUBLE)
        {
            //SetChargeCoulomb(_electron_);]
            set_charge_coulomb(_electron_);
        }
        else
        {
            throw;
        }

        _mode_ = mode;
        return;
    }

    void SetChargeCoulomb(const double charge)
    {
        if(_mode_ == ElectronicChargeMode::DOUBLE)
        {
            _coulomb_ = charge;
        }
        else if(_mode_ == ElectronicChargeMode::INTEGER)
        {
            if(check_overflow(charge))
            {
                throw __func__;
            }
            _electron_ = coulomb_to_electron(charge);
        }
        else
        {
            throw __func__;
        }
    }

    void SetChargeElectron(const uint64_t charge)
    {
        if(_mode_ == ElectronicChargeMode::DOUBLE)
        {
            _coulomb_ = electron_to_coulomb(charge);
        }
        else if(_mode_ == ElectronicChargeMode::INTEGER)
        {
            _electron_ = charge;
        }
        else
        {
            throw __func__;
        }
    }

    private:

    // how to measure
    ElectronicChargeMode _mode_;

    // in unit of coulomb
    double _coulomb_;

    // in unit of e
    uint64_t _electron_;

    // physics constant e, fundamental unit of electronic charge in coulomb
    static constexpr const double _electron_charge_fundamental_{1.6021766208e-19};
    // the maximum charge that can be stored in both double format
    // and uint64_t format without overflow in uint64_t
    static constexpr const uint64_t _uint64_t_max_{std::numeric_limits<uint64_t>::max()}; // 9223372036854775807
    static constexpr const double _max_electron_count_uint64_as_double_{((double)_uint64_t_max_ - 0.5) * _electron_charge_fundamental_};

    public:
    static constexpr const double ELECTRON_CHARGE{-_electron_charge_fundamental_};
    static constexpr const double ELECTRON_VOLT{_electron_charge_fundamental_};
};


class IonizationEvent
{

    public:

    IonizationEvent(const vector3<double>& position)
        : _position_{position}
        , _charge_((uint64_t)0)
    {
    }

    IonizationEvent(const vector3<double>& position, const uint64_t num_electron)
        : _position_{position}
        , _charge_(num_electron)
    {
    }

    const vector3<double>& GetPosition() const
    {
        return _position_;
    }

    uint64_t GetCharge() const
    {
        return _charge_.GetChargeElectron();
    }


    private:

    vector3<double> _position_;
    ElectronicCharge _charge_;

};


class Ion
{

    public:

    Ion(const double mass, const double charge)
        : _mass_{mass}
        , _charge_{charge}
    {
    }

    Ion(const double mass, const double charge, const vector3<double> position, const vector3<double> momentum)
        : _mass_{mass}
        , _charge_{charge}
        , _position_{position}
        , _momentum_{momentum}
    {
    }

    ////////////////////////////////////////////////////////////////////////////
    // public functions
    ////////////////////////////////////////////////////////////////////////////

    virtual
    double Charge() const
    {
        return _charge_;
    }

    virtual
    double Mass() const
    {
        return _mass_;
    }

    vector3<double> Momentum() const
    {
        return _momentum_;
    }

    vector3<double> Position() const
    {
        return _position_;
    }

    void SetMomentum(const vector3<double>& momentum)
    {
        _momentum_ = momentum;
    }

    void SetPosition(const vector3<double>& position)
    {
        _position_ = position;
    }

    double KE() 
    {
        calc_ke();
        return _kinetic_energy_;
    }

    double Velocity() const
    {
        return velocity();
    }

    /*
    double VelocitySquare() const
    {
        return velocity2();
    }
    */

    // step function
    // step particle forwards by distance length in direction of
    // momentum and update momentum magnitude
    // this is used to step the particle through when simulating collisions
    // using exponentially distributed random variables
    /*
    void Step(const double length, const double potential_change)
    {
        std::cout << "_position_=" << _position_ << std::endl;
        std::cout << "length=" << length << std::endl;
        std::cout << "_momentum_=" << _momentum_ << std::endl;
        std::cout << "Unit=" << _momentum_.Unit() << std::endl;
        _position_ += length * _momentum_.Unit();
        std::cout << "_position_=" << _position_ << std::endl;
        double potential_energy_change = _charge_ * potential_change;
        std::cout << "potential_energy_change=" << potential_energy_change << std::endl;
        double kinetic_energy{KE() - potential_energy_change}; // subtract because -PE = +KE
        std::cout << "_momentum_=" << _momentum_ << std::endl;
        _momentum_ = std::sqrt(2.0 * kinetic_energy * _mass_) * _momentum_.Unit();
        std::cout << "_momentum_=" << _momentum_ << std::endl;
    }
    */

    void SetKE(const double kinetic_energy)
    {
        _momentum_ = std::sqrt(2.0 * kinetic_energy * _mass_) * _momentum_.Unit();
    }

    private:

    ////////////////////////////////////////////////////////////////////////////
    // private functions
    ////////////////////////////////////////////////////////////////////////////

    void calc_ke()
    {
        _kinetic_energy_ = 0.5 * _mass_ * velocity2();
    }

    double velocity2() const
    {
        return ((1.0 / _mass_) * _momentum_).Mod2();
    }

    double velocity() const
    {
        return ((1.0 / _mass_) * _momentum_).Mod();
    }

    ////////////////////////////////////////////////////////////////////////////
    // member variables
    ////////////////////////////////////////////////////////////////////////////

    // charge of ion
    double _charge_;
    
    // kinetic variables
    // position, momentum and kinetic energy
    double _mass_;
    vector3<double> _position_;
    vector3<double> _momentum_;
    double _kinetic_energy_;

};


class Electron : public Ion
{

    public:

    Electron(const double mass, const vector3<double> position, const vector3<double> momentum)
        : Ion(mass, ElectronicCharge::ELECTRON_CHARGE, position, momentum)
    {
    }

    ~Electron()
    {
    }

    // Note: charge variable is unused
    /*
    double Charge() const
    {
        return IonizationEvent::ELECTRON_CHARGE;
    }
    */


    private:




};

#endif // IONIZATIONEVENT_HPP

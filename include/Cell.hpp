#ifndef CELL_HPP
#define CELL_HPP


#include "vector3.hpp"
#include "Generator.hpp"
#include "Wire.hpp"
#include "EndCap.hpp"
#include "IonizationEvent.hpp"

#include "function_debug.hpp"

#include "HistogramWrapper.hpp"

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>


const double CELL_ENDCAP_LENGTH{0.1}; // 100 mm
const double WIRE_RADIUS{4.0e-5}; // 40 um

class Cell
{

    public:

    // volume is the world volume used for generating points for events
    Cell(const double length, const double radius, const double anode_voltage, const vector3<double>& position, const vector3<double>& volume)
        //: _position_(0.0, 0.0, 0.0)
        //, _direction_(0.0, 0.0, 1.0) // anode wire points along z
        //: _length_{length}
        //, _radius_{radius}
        : _cube_(vector3<double>(2.0 * radius, 2.0 * radius, length))
        , _wire_anode_(vector3<double>(radius, radius, 0.0), vector3<double>(0.0, 0.0, 1.0), length, WIRE_RADIUS, anode_voltage)
        , _tfile_{new TFile("Cell_output.root", "recreate")}
        , _histogram_group_(_tfile_, HistogramGroupProperties("", "./canvas"), HistogramProperties(100, 0.0, 1.0))
        , properties("h_ionization_count_wrapper", 10, 0.0, 1.0e1)
        , h_ionization_count_wrapper(_tfile_, properties)
        , h_total_ionization_count_wrapper(_tfile_, HistogramProperties("h_total_ionization_count", 20, 0.0, 20.0))
        , h_exponential_wrapper(_tfile_, HistogramProperties("h_exponential_wrapper", 100, 0.0, 1.0e-1))
        , h_ionization_position(_tfile_, HistogramProperties("h_ionization_position", 100, 0.0, 2.0 * radius))
    {

        DEBUG_MESSAGE(function_debug_arguments(__PRETTY_FUNCTION__, "volume", volume.String()));

        // set position of cube (volume of cell)
        _cube_.SetPosition(position);

        init_wire(anode_voltage);
        vector3<double> endcap_position(0.0, 0.0, 0.0);
        endcap_position += _cube_.Position();
        init_endcap(endcap_position, vector3<double>(0.0, 0.0, 1.0), CELL_ENDCAP_LENGTH, radius);

        // initialize the electric field
        init_electric_field();

        //f = new TFile("Cell_output.root");
        //t = new TTree("celloutput");
        //h_event_pos_x = new TH1F("h_event_pos_x", "h_event_pos_x", 100, -5.0 * radius, 5.0 * radius);
        //h_event_pos_y = new TH1F("h_event_pos_y", "h_event_pos_y", 100, -5.0 * radius, 5.0 * radius);
        //h_event_pos_z = new TH1F("h_event_pos_z", "h_event_pos_z", 100, -5.0 * length, 5.0 * length);

        _histogram_group_.Add("h_event_pos_x", 100, 0.0, volume.GetX());
        _histogram_group_.Add("h_event_pos_y", 100, 0.0, volume.GetY());
        _histogram_group_.Add("h_event_pos_z", 100, 0.0, volume.GetZ());

        //_tfile_->SetDirectory(0);

        h_perpendicular_distance = new TH1F("h_perpendicular_distance", "h_perpendicular_distance", 100, 0.0, 0.05);

        h_cell_voltage = new TH1F("h_cell_voltage", "h_cell_voltage", 100000, 0.0, 0.020);

        for(Int_t i{1}; i <= h_cell_voltage->GetNbinsX(); ++ i)
        {
            Double_t radial_position{h_cell_voltage->GetBinCenter(i)};
            Double_t voltage{V_radial(radial_position)};
            h_cell_voltage->SetBinContent(i, voltage);
        }

        h_event_initial_voltage = new TH1F("h_event_initial_voltage", "h_event_initial_voltage", 100, 0.0, 2000.0);

        //h_ionization_count = new TH1F("h_ionization_count", "h_ionization_count", 100, 0.0, 2.0e4);
    }

    ~Cell()
    {
        h_ionization_count_wrapper.Write();
        h_ionization_count_wrapper.Canvas();
        
        h_total_ionization_count_wrapper.Write();
        h_total_ionization_count_wrapper.Canvas();


        h_exponential_wrapper.Write();
        h_exponential_wrapper.Canvas();

        h_ionization_position.Write();
        h_ionization_position.Canvas();
        
        h_event_initial_voltage->Write();
        TCanvas* c_event_initial_voltage = new TCanvas("c_event_initial_voltage", "c_event_initial_voltage", 800, 600);
        c_event_initial_voltage->SetLogy();
        h_event_initial_voltage->Draw();
        c_event_initial_voltage->SaveAs("c_event_initial_voltage.png");
        delete c_event_initial_voltage; 

        h_cell_voltage->Write();
        TCanvas* c_cell_voltage = new TCanvas("c_cell_voltage", "c_cell_voltage", 800, 600);
        h_cell_voltage->Draw();
        c_cell_voltage->SaveAs("c_cell_voltage.png");
        delete c_cell_voltage; 

        h_perpendicular_distance->Write();
        TCanvas* c = new TCanvas("c_perpendicular_distance", "c_perpendicular_distance", 800, 600);
        h_perpendicular_distance->Draw();
        c->SaveAs("c_perpendicular_distance.png");
        delete c;

        //Canvas(h_event_pos_x, "h_event_pos_x");
        //Canvas(h_event_pos_y, "h_event_pos_y");
        //Canvas(h_event_pos_z, "h_event_pos_z");

        _histogram_group_.Write();
        _histogram_group_.Canvas();
        
        //h_event_pos_x->Write();
        //h_event_pos_y->Write();
        //h_event_pos_z->Write();


        //t->Write();
        _tfile_->Close();

        //delete _tfile_;
    }

    
    /*
    void Canvas(const TH1* histogram, const std::string& name) const
    {
        //h->SetStats(0); // move to init function

        TCanvas *c = new TCanvas(name.c_str(), name.c_str(), 800, 600);
        histogram->Draw();
        c->SaveAs(name + std::string(".png"));
        delete c;
    }
    */
    
    void init_electric_field()
    {

    }

    // convert radial distance from anode wire to voltage (electric potential)
    // radial position is relative to the anode wire position
    double V_radial(const double radial_position)
    {
        
        //DEBUG_MESSAGE(function_debug_arguments(__PRETTY_FUNCTION__, "radial_position", radial_position));

        // wire voltage
        const double V0{_wire_anode_.GetVoltage()};

        // wire radius
        const double wire_radius{_wire_anode_.GetCylinder().Radius()};
        // cell radius
        // TODO: assumed radial field configuration THIS IS WRONG
        const double radius{_cube_.Size().GetX() / 2.0}; 
        
        // is position is inside wire, return voltage of wire
        // TODO: no check before this function call to check if position is
        // valid (inside wire)
        if(radial_position <= wire_radius)
        {
            return V0;
        }
        else if(radial_position >= radius)
        {
            return 0.0;
        }
        
        // compute potential if outside wire
        const double ln_R{std::log(radius)};
        const double ln_r0{std::log(wire_radius)};
        const double ln_r{std::log(radial_position)};

        //DEBUG_MESSAGE(function_debug_locals(__PRETTY_FUNCTION__, "wire_radius", wire_radius, "radius", radius, "ln_R", ln_R, "ln_r0", ln_r0, "ln_r", ln_r));

        return V0 * ((ln_R - ln_r) / (ln_R - ln_r0));

    }

    // for electron, solve for radial position, given a required energy gain
    // and an initial radial position
    // energy drop specified in electron volt
    bool solve_radial(double& output_radial_position, const double radial_position, const double energy_drop)
    {
        // electron charge
        const double e{ElectronicCharge::ELECTRON_CHARGE};

        // anode wire voltage constant
        const double V0{_wire_anode_.GetVoltage()};

        // delta_E / (e * V0)
        //const double k{energy_drop / (e * V0)}; // joule
        const double k{energy_drop / (1.0 * V0)}; // eV

        // wire radius
        const double r0{_wire_anode_.GetCylinder().Radius()};

        // maximum radius
        const double R{_cube_.Size().GetX() / 2.0};

        // radius ratio
        const double r_frac{r0 / R};

        // multiplier
        const double mult{std::pow(r_frac, k)};

        // TODO: mult is a constant - make static, member of Cell

        // solution
        const double r_final{radial_position * mult};

        //DEBUG_MESSAGE(function_debug_locals(__PRETTY_FUNCTION__, "e", e, "V0", V0, "k", k, "r0", r0, "R", R, "r_frac", r_frac, "mult", mult, "radial_position", radial_position, "r_final", r_final));

        // check if solution exists / is valid
        if(r_final < r0)
        {
            return false;
        }
        
        // solution exists
        output_radial_position = r_final;

        return true;


        // TODO: inefficient if voltage already calculated

        //const double voltage_initial{V_radial(radial_position)};
        //const double energy_initial{e * voltage_initial};
        //const double energy_final{};
    }
    
    bool solve_radial_and_energy(double& output_radial_position, double& output_energy, const double radial_position, const double energy, const double energy_drop)
    {
        // solve radial position
        if(solve_radial(output_radial_position, radial_position, energy_drop) == false)
        {
            // solution does not exist
            return false;
        }
        
        // solution exists
        // compute energy gain
        // output = input + gain due to potential - ionization
        // all units are eV
        const double energy_gain{V_radial(output_radial_position) - V_radial(radial_position)};
        //std::cout << "the energy gain is " << energy_gain << " eV" << std::endl;
        //std::cin.get();
        // TODO: This always returns 24.something eV, as this is the ionization energy
        // an over-constrained problem?
        // what about when initial energy non-zero?

        // TODO
        // not yet using gas mean free path information!

        output_energy = energy + energy_gain - energy_drop;

        return true;
    }

    // position is relative to the global origin
    double electric_potential(Generator& generator, vector3<double> position)
    {
        // TODO: should use get_placement_relative(), a function to convert
        // to global position using a heirarchy of objects, eg; a cell contains
        // a wire, and the position of the wire is relative to the cell
        // so the position of the wire relative to the global origin is given
        // by the sum
        vector3<double> cell_position{_cube_.Position()};
        vector3<double> wire_position{_wire_anode_.GetCylinder().Position()};
        vector3<double> wire_direction{_wire_anode_.GetCylinder().Direction()};
        vector3<double> delta_position{position - (cell_position + wire_position)};
        // TODO: derive this algorithm and optimize
        double distance_along{dot(delta_position, wire_direction)};
        vector3<double> perpendicular{delta_position - distance_along * wire_direction};
        double perpendicular_distance{perpendicular.Length()};

        /*
        std::cout << "wire_position=" << wire_position + cell_position << std::endl;
        std::cout << "wire_direction=" << wire_direction << std::endl;
        std::cout << "delta_position=" << delta_position << std::endl;
        std::cout << "distance_along=" << distance_along << std::endl;
        std::cout << "perpendicular=" << perpendicular << std::endl;
        std::cout << "perpendicular_distance=" << perpendicular_distance << std::endl;
        */

        h_perpendicular_distance->Fill(perpendicular_distance);
        
        double voltage{V_radial(perpendicular_distance)};
        h_event_initial_voltage->Fill(voltage);
        // solve to find radial position of next ionization event
        //double next_radial{0.0};
        //double next_energy;

        // https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable1.htm
        const double HELIUM_IONIZATION_ENERGY{24.587387}; // eV 
        const double ionization_energy_ev{85.7}; // eV
        const double ionization_energy{85.7 * ElectronicCharge::ELECTRON_VOLT}; // fake ionization energy eV //{HELIUM_IONIZATION_ENERGY};
        //int ionization_count{0};

        // He mean free path for electrons
        //const double HELIUM_MEAN_FREE_PATH{1.0e-6}; // 1 um
        const double HELIUM_MEAN_FREE_PATH{1.0e-2}; // change to mean of ~ 100 steps

        // structure containing electrons and properties
        const double electron_mass{9.10938356e-31};
        std::vector<Electron> electron;
        // TODO: position is not a true position, using x value as radial
        // distance
        // kinetic energy = 0.0
        // momentum is zero
        electron.push_back(Electron(electron_mass, vector3<double>(perpendicular_distance, 0.0, 0.0), vector3<double>(0.0, 0.0, 0.0)));

        // process all electrons
        for(std::size_t ix{0}; ix < electron.size(); ++ ix)
        {
            //std::cout << "electron " << ix + 1 << " of " << electron.size() << std::endl;

            // count number of new ionization events
            std::size_t ionization_count{0};
            
            const double initial_position{electron[ix].Position().GetX()};
            //std::cout << "initial potential energy is " << _wire_anode_.GetVoltage() - V_radial(initial_position) << " which is enough for " << (int)((_wire_anode_.GetVoltage() - V_radial(initial_position)) / ionization_energy_ev) <<  " ionization events" << std::endl;


            // loop until electron reaches anode wire
            for(;;)
            {

                // draw exponentially distributed random number
                double distance{generator.GetRandomExponential(1.0 / HELIUM_MEAN_FREE_PATH)};
                h_exponential_wrapper.Ref().Fill(distance);

                // initial and final radial positions
                double initial_position{electron[ix].Position().GetX()};
                double next_position{initial_position - distance};

                if(next_position < _wire_anode_.GetCylinder().Radius())
                {
                    break;
                }

                // calculate potential_drop
                double initial_potential{V_radial(initial_position)};
                double next_potential{V_radial(next_position)};
                double potential_change{next_potential - initial_potential};
                //DEBUG_MESSAGE(function_debug_locals(__PRETTY_FUNCTION__, "initial_p", initial_potential, "next_p", next_potential));
                //std::cout << "potential_change=" << potential_change << std::endl;

                // move electron to this radial position
                //electron[ix].Step(distance, potential_change);
                electron[ix].SetPosition(vector3<double>(next_position, 0.0, 0.0));
                double potential_energy_change{electron[ix].Charge() * potential_change};
                double intermediate_kinetic_energy{electron[ix].KE() - potential_energy_change};
                
                // only required for debug statement
                double intermediate_momentum_magnitude{std::sqrt(2.0 * intermediate_kinetic_energy * electron[ix].Mass())};
                electron[ix].SetMomentum(intermediate_momentum_magnitude * vector3<double>(-1.0, 0.0, 0.0));
                //std::cout << "distance=" << distance << " delta_potential=" << potential_change << " eV ke=" << electron[ix].KE() / ElectronicCharge::ELECTRON_VOLT << " eV" << std::endl;

                // check if kinetic energy is great enough to ionize
                // convert to SI unit
                if(intermediate_kinetic_energy >= ionization_energy)
                {
                    h_ionization_position.Ref().Fill(initial_position);

                    // TODO: correct collision formula

                    // set kinetic energy to half the remaining energy after ionization
                    double next_kinetic_energy{0.5 * (intermediate_kinetic_energy - ionization_energy)};
                    double next_momentum_magnitude{std::sqrt(2.0 * next_kinetic_energy * electron[ix].Mass())};

                    electron[ix].SetMomentum(next_momentum_magnitude * vector3<double>(-1.0, 0.0, 0.0));

                    // next position and momentum for new electron created in ionization process
                    vector3<double> next_momentum{electron[ix].Momentum()};
                    vector3<double> next_position{electron[ix].Position()};

                    // add new electron to vector
                    electron.push_back(Electron(electron_mass, next_position, next_momentum));

                    ++ ionization_count;

                    //std::cout << "ionization event ke=" << electron[ix].KE() / ElectronicCharge::ELECTRON_VOLT << " eV" << std::endl;
                }


                //std::cin.get();
            }

            h_ionization_count_wrapper.Ref().Fill(ionization_count);

            //std::cout << "ionization_count=" << ionization_count << std::endl;
            //std::cin.get();



            /*
            // set the initial radial value of this electron
            //double initial_radial{*it};
            //double initial_radial{electron_radial_position[ix]};
            //double initial_energy{electron_energy[ix]};
            double initial_radial{electron[ix].Position().Mod()};
            double initial_energy{electron[ix].KE()};
            //std::cout << "initial_radial=" << initial_radial << std::endl;
            //std::cin.get();

            double ionization_energy_minus_kinetic_energy{ionization_energy - initial_energy};
            if(ionization_energy_minus_kinetic_energy < 0.0) ionization_energy_minus_kinetic_energy = 0.0;

            // solve for radial position
            while(solve_radial_and_energy(next_radial, next_energy, initial_radial, initial_energy, ionization_energy_minus_kinetic_energy))
            //while(solve_radial(next_radial, initial_radial, ionization_energy_minus_kinetic_energy))
            {
                // if there was a solution, then an ionization event occured
                // TODO: store the kinetic energy / velocity as well as position
                //++ ionization_count;
                //std::cout << "next ionization occurs at position: " << next_radial << std::endl;

                // compute remaining energy
                //next_energy{initial_energy + energy_gain};
                //energy_gain = 0.0;

                // add new electron to vector at this position
                //electron_radial_position.push_back(next_radial);
                //electron_energy.push_back(next_energy);
                vector3<double> next_momentum(-1.0, 0.0, 0.0);
                next_momentum *= std::sqrt(2.0 * next_energy * electron_mass);
                electron.push_back(Electron(electron_mass, vector3<double>(next_radial, 0.0, 0.0), next_momentum));

                // reset the while loop, try to find another solution
                initial_radial = next_radial;
                next_radial = 0.0;
                initial_energy = next_energy;
                next_energy = 0.0;

            }
            // when this loop ends, some new electrons may have been added to
            // the vector

            //std::cout << "NEXT ELECTRON" << std::endl;
        
            //std::cout << "electron: ix=" << ix << " ionization_count=" << ionization_count << std::endl;
            */
        }

        h_total_ionization_count_wrapper.Ref().Fill(electron.size() - 1);

        return voltage;

    }

    vector3<double> electric_field(vector3<double> position)
    {
        // subtract relative position of cell volume
        position -= _cube_.Position();

        // anode wire position
        // this is relative to the cube position
        vector3<double> anode_wire_position{_wire_anode_.GetCylinder().Position()};

        vector3<double> field_e;

        double electric_field_magnitude{0.0}; // TODO

        field_e = anode_wire_position - position;
        field_e *= electric_field_magnitude;

    }

    
    void SetPosition(const vector3<double> position)
    {
        std::cout << "TODO" << std::endl;
        //_position_ = position;
    }

    IonizationEvent GenerateIonizationEvent(Generator& generator, const vector3<double>& volume) const
    {

        std::size_t counter{0};
        const std::size_t COUNTER_MAX{100000};
        

        vector3<double> event_position;
        for(;; ++ counter)
        {
            // generate event within volume
            event_position = generator.GetRandomPosition(volume);

            /*_histogram_group_.Ref("h_event_pos_x").Get()->Fill(event_position.GetX());
            _histogram_group_.Ref("h_event_pos_y").Get()->Fill(event_position.GetY());
            _histogram_group_.Ref("h_event_pos_z").Get()->Fill(event_position.GetZ());
            */
            
            // check if event is within cell volume
            // TODO: don't need this because no "track"
            // when generating muon tracks, this was required
            // now it is not

            if(_cube_.PointIntersectionTest(event_position))
            {
                //_histogram_group_.Ref("h_event_pos_x").Ref().Fill(event_position.GetX());
                //_histogram_group_.Ref("h_event_pos_y").Ref().Fill(event_position.GetY());
                //_histogram_group_.Ref("h_event_pos_z").Ref().Fill(event_position.GetZ());
                _histogram_group_.Ref("h_event_pos_x").Fill(event_position.GetX());
                _histogram_group_.Ref("h_event_pos_y").Fill(event_position.GetY());
                _histogram_group_.Ref("h_event_pos_z").Fill(event_position.GetZ());

                break;
            }

            if(counter > COUNTER_MAX) break;

        }
        return IonizationEvent(event_position, (uint64_t)1);
    }


    private:

    
    void init_wire(const double anode_voltage)
    {
        // initialized anode wire
        //_wire_anode_.Init(_position_, _direction_, _length_, _radius_, anode_voltage);
        //_wire_anode_.SetVoltage(anode_voltage);
        
        // initialize ground wire vector
        _wire_ground_.clear();

        vector3<double> position_origin{_wire_anode_.GetCylinder().Position()};
        vector3<double> direction{_wire_anode_.GetCylinder().Direction()};
        //double length{_cube_.Size().GetZ()};
        double length{_wire_anode_.GetCylinder().Length()};
        double radius{0.5 * _cube_.Size().GetX()};

        for(int y{-1}; y <= 1; ++ y)
        {
            for(int x{-1}; x <= 1; ++ x)
            {
                if(x == 0 && y == 0) continue;

                //vector3<double> position{_position_};
                const vector3<double> x_delta(radius, 0.0, 0.0);
                const vector3<double> y_delta(0.0, radius, 0.0);
                vector3<double> position{position_origin};
                position += x * x_delta;
                position += y * y_delta;
                _wire_ground_.push_back(Wire(position, direction, length, WIRE_RADIUS, 0.0));
            }
        }
    }

    void init_endcap(vector3<double>& position, const vector3<double> direction, const double length, const double radius)
    {
        _endcap_0_.Init(position, direction, length, radius, 0.0);
        position += vector3<double>(0.0, 0.0, length);
        _endcap_1_.Init(position, -direction, length, radius, 0.0); // TODO radius should be slightly smaller
    }
    


    //vector3<double> _position_; // vector to anode wire start position
    //vector3<double> _direction_; // unit vector in direction along anode wire
    //double _length_;
    //double _radius_;
    Geometry::Cube _cube_; // cell volume
    const double _ACTIVE_RADIUS_FRACTION_{0.9}; // the active radius is 
    // TODO

    Wire _wire_anode_;
    /*
    Wire _wire_ground_0_;
    Wire _wire_ground_1_;
    Wire _wire_ground_2_;
    Wire _wire_ground_3_;
    Wire _wire_ground_4_;
    Wire _wire_ground_5_;
    Wire _wire_ground_6_;
    Wire _wire_ground_7_;
    */
    std::vector<Wire> _wire_ground_;


    EndCap _endcap_0_;
    EndCap _endcap_1_;


    // ROOT data output
    TFile *_tfile_;
    //TTree *t;
    //TH1F *h_event_pos_x;
    //TH1F *h_event_pos_y;
    //TH1F *h_event_pos_z;
   
    mutable HistogramGroupFloat _histogram_group_;

    mutable TH1F *h_perpendicular_distance;
    
    mutable TH1F *h_cell_voltage;
    
    mutable TH1F *h_event_initial_voltage;

    //mutable TH1F *h_ionization_count;
    HistogramProperties properties;
    HistogramWrapperFloat h_ionization_count_wrapper;
    HistogramWrapperFloat h_total_ionization_count_wrapper;
    HistogramWrapperFloat h_exponential_wrapper;
    HistogramWrapperFloat h_ionization_position;

};






#endif // CELL_HPP

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
const double WIRE_RADIUS{1.0e-6}; // 1 um

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
        , _histogram_group_(_tfile_, HistogramGroupProperties("", "."), HistogramProperties(100, 0.0, 1.0))
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

        h_perpendicular_distance = new TH1F("h_perpendicular_distance", "h_perpendicular_distance", 100, -0.1, 0.5);

    }

    ~Cell()
    {
        h_perpendicular_distance->Write();

        TCanvas* c = new TCanvas("c_perpendicular_distance", "c_perpendicular_distance", 800, 600);
        h_perpendicular_distance->Draw();
        c->SaveAs("c_perpendicular_distance.png");
        delete c;

        //Canvas(h_event_pos_x, "h_event_pos_x");
        //Canvas(h_event_pos_y, "h_event_pos_y");
        //Canvas(h_event_pos_z, "h_event_pos_z");

        _histogram_group_.Canvas();

        //h_event_pos_x->Write();
        //h_event_pos_y->Write();
        //h_event_pos_z->Write();

        _histogram_group_.Write();

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
    double V_radiual(const double radial_position)
    {
        const double radius{_cube_.Size().GetX() / 2.0};
        const double wire_radius{_wire_anode_.GetCylinder().Radius()};
        const double ln_R{std::log(radius)};
        const double ln_r0{std::log(wire_radius)};
        const double ln_r{std::log(radial_position)};
        const double V0{_wire_anode_.GetVoltage()};
        return V0 * ((ln_R - ln_r) / (ln_R - ln_r0));
    }

    // for electron, solve for radial position, given a required energy gain
    // and an initial radial position
    bool solve_radial(double& output_radial_position, const double radial_position, const double energy_drop)
    {
        const double e{ElectronicCharge::ELECTRON_CHARGE};
        // TODO check solution exists
    }

    // position is relative to the global origin
    double electric_potential(vector3<double> position)
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

        std::cout << "wire_position=" << wire_position + cell_position << std::endl;
        std::cout << "wire_direction=" << wire_direction << std::endl;
        std::cout << "delta_position=" << delta_position << std::endl;
        std::cout << "distance_along=" << distance_along << std::endl;
        std::cout << "perpendicular=" << perpendicular << std::endl;
        std::cout << "perpendicular_distance=" << perpendicular_distance << std::endl;

        h_perpendicular_distance->Fill(perpendicular_distance);

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

};






#endif // CELL_HPP

#ifndef TRACKERGAS_HPP
#define TRACKERGAS_HPP



class TrackerGas
{


    private:

        
    double _ionization_energy_; // first ionization energy of gas atom / molecule
    double _gamma_mean_free_path_; // assume all energies of photons have same mean free path
    double _gas_fraction_ionizable_;
    double _gas_fraction_quenching_;

};


#endif

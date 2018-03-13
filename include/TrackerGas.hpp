#ifndef TRACKERGAS_HPP
#define TRACKERGAS_HPP



class TrackerGas
{


    public:

    // electron mean free path:
    // https://agenda.infn.it/getFile.py/access?contribId=27&resId=0&materialId=slides&confId=4542

    TrackerGas()
        : _electron_mean_free_path_{2.7e-6} // 2.7 um in argon
    {
    }


    private:

        
    double _ionization_energy_; // first ionization energy of gas atom / molecule
    double _gamma_mean_free_path_; // assume all energies of photons have same mean free path
    double _electron_mean_free_path_; // assume all energies of electrons have same mean free path
    double _gas_fraction_ionizable_;
    double _gas_fraction_quenching_;

};


#endif // TRACKERGAS_HPP

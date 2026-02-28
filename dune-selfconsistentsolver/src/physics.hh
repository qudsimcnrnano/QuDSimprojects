
#ifndef QUDSIM_PHYSICS_HH
#define QUDSIM_PHYSICS_HH

#include <cmath>
#include <string>

namespace QuDSim {

namespace Constants {
    constexpr double q0    = 1.0;
    constexpr double m0    = 1.0;
    constexpr double hbar  = 1.0;
    constexpr double eps0  = 1.0/(4.0*M_PI);
    constexpr double Hartree_to_eV  = 27.211;
    constexpr double eV_to_Hartree  = 1.0/27.211;
    constexpr double Bohr_to_nm     = 0.05291;
    constexpr double nm_to_Bohr     = 1.0/0.05291;
    constexpr double Bohr_to_cm     = 5.291e-9;
    constexpr double Hartree_to_K   = 3.157e5;
    constexpr double kB_eV          = 8.61735e-5;
    constexpr double pi = M_PI;
}

inline double eV_to_Ha(double eV) { return eV * Constants::eV_to_Hartree; }
inline double Ha_to_eV(double Ha) { return Ha * Constants::Hartree_to_eV; }
inline double nm_to_au(double nm) { return nm * Constants::nm_to_Bohr; }
inline double au_to_nm(double au) { return au * Constants::Bohr_to_nm; }
inline double cm3_to_au3(double n_cm3) {
    double au = Constants::Bohr_to_cm;
    return n_cm3 * au * au * au;
}

struct Material {
    std::string name;
    double eps_r;        // relative permittivity
    double Eg;           // band gap (eV)
    double chi;          // electron affinity (eV)
    double me_eff;       // electron effective mass (m0)
    double mh_eff;       // hole effective mass (m0)
    double VBO_electron; // conduction band offset vs Si (eV)
    double VBO_hole;     // valence band offset vs Si (eV)
};

namespace Materials {
    const Material Silicon = {"Silicon", 11.7, 1.12, 4.05, 0.26, 0.49, 0.0, 0.0};
    const Material SiO2    = {"SiO2",    3.9, 9.0, 0.95, 0.5, 0.5, 3.1, 4.5};
    const Material HfO2    = {"HfO2",   25.0, 5.68, 2.05, 0.2, 0.15, 1.5, 3.4};
    const Material Al2O3   = {"Al2O3",   9.0, 6.4, 1.35, 0.4, 0.4, 2.8, 2.6};

    inline Material getMaterialByRegion(int regionID) {
        switch(regionID) {
            case 200: case 201: case 202: return Silicon;
            case 400: return SiO2;
            case 501: return Al2O3;
            case 502: return HfO2;
            default:  return Silicon;
        }
    }
}

struct SimulationParams {
    double Na = 1e18;
    double Nd = 0.0;
    double ni = 1.5e10;
    double T = 300.0;
    double Vg = 0.0;
    double Vfb = -0.88;
    int    max_iterations = 100;
    double tolerance      = 1.0e-6;
    double mixing_alpha   = 0.3;
    int num_eigenvalues = 10;
    bool solve_electrons = true;
    bool solve_holes     = true;

    double beta, Ef, Np0, Pf;

    void initialize() {
        using namespace Constants;
        double Na_au = cm3_to_au3(Na);
        double ni_au = cm3_to_au3(ni);
        beta = Hartree_to_eV / (kB_eV * T);
        Np0 = 0.5 * (Na_au + std::sqrt(Na_au * Na_au + 4.0 * ni_au * ni_au));
        Ef = (1.0 / beta) * std::log(Np0 / ni_au);
        Pf = (1.0 / beta) * std::log(Na_au / ni_au);
    }

    double getVss() const {
        double Eg_Ha = eV_to_Ha(Materials::Silicon.Eg);
        return Eg_Ha / 2.0 + Pf + eV_to_Ha(Vg);
    }
};

} // namespace QuDSim
#endif

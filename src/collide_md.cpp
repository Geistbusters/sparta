#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "collide_md.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "fix_vibmode.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

//#include "modify.h"

#include "units.h"
#include "math_utils.h"
#include "collision_data.h"
#include "md_integrator.h"
#include "rovib_utils.h"
#include "jv_utils.h"
#include "init_utils.h"
#include "output_manager.h"
#include "qct_collisions.h"
#include "sampling_utils.h"
#include "input_parser.h"

#include <filesystem>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono> // For time-based seeding



using namespace SPARTA_NS;
using namespace MathConst;



enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define MAXLINE 1024
    //CollisionProcessor::CollisionProcessor() : verbosity_(2) {
    //    integrator_.setVerbosity(verbosity_);
    	// Enable surface hopping for nonadiabatic trajectories
   //     integrator_.enableSurfaceHopping(false);
    //    integrator_.setInitialSurface(0);  // Start on ground state (surface 0)
   // }
#ifdef SPARTA_MD_BUILD

CollideMD::CollideMD(SPARTA *sparta, int narg, char **arg) : Collide(sparta, narg, arg){
	//SPARTA STUFF
 	DEBUG_PRINT("DEBUG: Constructor called with narg=%d\n" << narg);
	DEBUG_PRINT( "DEBUG: mixture pointer = " << mixture );
	DEBUG_PRINT("DEBUG: sizeof(Params) = %zu\n" << sizeof(Params));
	DEBUG_PRINT("DEBUG: nparams = %d\n" << nparams);
	DEBUG_PRINT("DEBUG: broadcast size = %zu bytes\n" <<  nparams*nparams*sizeof(Params));	
	DEBUG_PRINT( "DEBUG: sparta pointer = " << sparta );

	/*if (mixture) {
	DEBUG_PRINT( "DEBUG: mixture->nspecies = " << mixture->nspecies );
	} else {
	DEBUG_PRINT( "ERROR: mixture is NULL!" );
	}*/  

	//. SPARTA STUFF

	if (narg < 3) error->all(FLERR,"Illegal collide command");

	// optional args

	relaxflag = CONSTANT;

	int iarg = 3;
	while (iarg < narg) {
		if (strcmp(arg[iarg],"relax") == 0) {
		if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
		if (strcmp(arg[iarg+1],"constant") == 0) relaxflag = CONSTANT;
		else if (strcmp(arg[iarg+1],"variable") == 0) relaxflag = VARIABLE;
		else error->all(FLERR,"Illegal collide command");
		iarg += 2;
		} else error->all(FLERR,"Illegal collide command");
	}

	// proc 0 reads file to extract params for current species
	// broadcasts params to all procs

	nparams = particle->nspecies;
	if (nparams == 0)
		error->all(FLERR,"Cannot use collide command with no species defined");

	

	memory->create(params,nparams,nparams,"collide:params");
	memory->create(prefactor,nparams,nparams,"collide:prefactor");


	// if (comm->me == 0) read_param_file(arg[2]);
	// MPI_Bcast(params[0],nparams*nparams*sizeof(Params),MPI_BYTE,0,world);
	double b_max_ang = 6.0;      // Your desired b_max in Angstroms
	double b_max_m = b_max_ang * 1.0e-10;  // Convert to meters
	double d_collision = 2.0 * b_max_m;    // Collision diameter = 2 * b_max

	for (int i = 0; i < nparams; i++) {
		for (int j = 0; j < nparams; j++) {
			params[i][j].diam = d_collision;       
			params[i][j].omega = 0.5;              // Hard sphere
			params[i][j].tref = 300.0;             
			params[i][j].alpha = 1.0;              
			
			// HS cross section set to: σ = π * (5Å)²
	prefactor[i][j] = M_PI * b_max_m * b_max_m;
		}
	}


	//memory->create(params,nparams,nparams,"collide:params");
	//if (comm->me == 0) read_param_file(arg[2]);
	//MPI_Bcast(params[0],nparams*nparams*sizeof(Params),MPI_BYTE,0,world);

	// allocate per-species prefactor array

	//memory->create(prefactor,nparams,nparams,"collide:prefactor");
	

	// pes_path_ = "./";
	// verbosity_ = 2;
	// // Configure MD integrator
	// integrator_.setVerbosity(verbosity_);
	// integrator_.setPESPath(pes_path_);
	

	// if (integration_type_ == 10) {
	// 			integrator_.enableSurfaceHopping(true);
	// }

	// if (!integrator_.initializePES("O3_quintet")) {
	// 	error->all(FLERR, "Failed to initialize PES for MD collisions");
	// }
		

	//. CRDS STUFF
	  std::string input_file = "crds_input.txt"; 
    auto CRDSparams = parse_crds_input(input_file);
    config_ = populate_config(CRDSparams, 0);

    std::cout << " CONFIG(integration type): " << config_.integration_type << std::endl;
    std::cout << " CONFIG(initial surface): " << config_.initial_surface << std::endl;
    std::cout << " CONFIG(ADIABATIC): " << config_.adiabatic << std::endl;
    // if (!config_.adiabatic){
    //   integrator_.num_states=6;
    // }

    // Tell the OutputManager to use the config's verbosity
    OutputManager::set_verbosity(config_.verbosity); // (Or whatever the method is named in output_manager.h)
    
    // Tell the Integrator to use the config's verbosity
    integrator_.setVerbosity(config_.verbosity);

    // 2. Initialize the PES once!
    //pes_path_ = "./";
    //integrator_.setVerbosity(2);
    //integrator_.setPESPath(pes_path_);
    
    // Use the system name from your config (or hardcode "O3_quintet" if needed)
    if (!integrator_.initializePES("O3_quintet")) {
        error->all(FLERR, "Failed to initialize PES for MD collisions");
    }

	if (config_.integration_type >= 10) {
        integrator_.enableSurfaceHopping(true);
    }
	
	printf("DEBUG: CollideMD constructor end\n");


}


	void CollideMD::init(){
	//printf("DEBUG: typeid(update).name() = %s\n", typeid(update).name());  
	if (nparams != particle->nspecies)
    error->all(FLERR,"MD parameters do not match current species");

  Collide::init();
}

#else
CollideMD::CollideMD() {
 // if (narg < 3) error->all(FLERR,"Illegal collide command");  

  pes_path_ = "./";
  verbosity_ = 2;
  // Configure MD integrator
  integrator_.setVerbosity(verbosity_);
  integrator_.setPESPath(pes_path_);


  if (integration_type_ == 10) {
            integrator_.enableSurfaceHopping(true);
   }

  if (!integrator_.initializePES()) {
	  std::cerr <<  "Failed to initialize PES for MD collisions" <<std::endl;
  }
}
#endif

CollideMD::~CollideMD() {
}

//*********************************************************************************//
//									           //
//                                                                                 //
/////// PREPARE FUNCTION                                                           //
//                                                                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////






bool CollideMD::prepare(CRDSConfig& config, CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2, int integration_type) {
            
        
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//																				  //
	//                                                                                //
	//PASS DATA IN FROM SPARTA PARTICLES	                                          //
	//                                                                                //
	//                                                                                //
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	DEBUG_PRINT( "\n=== Testing prepare() ===" );
	collision.clear();
  if(!config.adiabatic){
    collision.integration.n_states = integrator_.getNumStates();
  }else{
     collision.integration.n_states = 1;
  }
   
    bool success = prepareFromParticles(config, collision, p1, p2);
    
    if (!success) return false;

    // Set up the integrator based on the prepared collision data
    //integrator.initializePES("O3_quintet"); // or dynamically set
    //integrator.setInitialSurface(collision.initial_state.electronic_state);
    
    if (integration_type >= 10) {
        integrator_.enableSurfaceHopping(true);
        
    }
    
    return true;
}


bool CollideMD::prepareFromParticles(CRDSConfig& config, CollisionData& collision, 
                                     Particle::OnePart* p1, Particle::OnePart* p2) {
    
    // 1. Get species names dynamically from SPARTA
    std::string p1_name = particle->species[p1->ispecies].id;
    std::string p2_name = particle->species[p2->ispecies].id;

    // 2. Calculate initial relative velocity (p1 - p2)
    std::vector<double> rel_vel(3, 0.0);
    for (int i = 0; i < 3; i++) {
        rel_vel[i] = (p1->v[i] - p2->v[i]) * Units::MPS_TO_AUVEL;
    }

    // 3. Figure out which one is the diatom so we grab the correct J, V, and Elec states
    bool is_p1_diatom = (p1_name.length() >= 2);
    Particle::OnePart* sparta_atom = is_p1_diatom ? p2 : p1;
    Particle::OnePart* sparta_diatom = is_p1_diatom ? p1 : p2;

    // 4. Extract Electronic State mapping
    double eelec_atom = sparta_atom->eelec * Units::J_TO_HARTREE;
    double eelec_diatom = sparta_diatom->eelec * Units::J_TO_HARTREE;
    int ies;
    MDIntegrator::getEState(ies, eelec_atom, eelec_diatom);

    // 5. Extract Rotational (J) and Vibrational (V) states
    // (If both were somehow atoms, xj and xv would default to 0.0 from SPARTA)
    double xj = sparta_diatom->xj;
    double xv = sparta_diatom->xv;

    // 6. Pass everything over the bridge to the MD library!
    // This utilizes your new QCT::prepareFromPreSampled function
    return QCT::prepareFromPreSampled(config, collision, p1_name, p2_name, 
                                      rel_vel, ies, xj, xv, comm->me);
}


// bool CollideMD::prepareFromParticles(CRDSConfig& config, CollisionData& collision, 
//                                      Particle::OnePart* p1, Particle::OnePart* p2) {
    
//     collision.clear();
//     //DEBUG_PRINT( "\n=== Preparing SPARTA particle collision ===" );
//     std::cout<< "\n=== Preparing SPARTA particle collision ===" <<std::endl;

//     try {
//         // 1. Determine which SPARTA particle is the atom and which is the diatom
//         char* p1_name = particle->species[p1->ispecies].id;
//         char* p2_name = particle->species[p2->ispecies].id;
        
//         bool is_p1_diatom = (strlen(p1_name) >= 2);
//         bool is_p2_diatom = (strlen(p2_name) >= 2);

//         Particle::OnePart* sparta_atom = nullptr;
//         Particle::OnePart* sparta_diatom = nullptr;
//         int a1 = -1, a2 = -1, a3 = -1;

//         if (is_p1_diatom && !is_p2_diatom) {
//             sparta_diatom = p1;
//             sparta_atom = p2;
//             a1 = 0; a2 = 1; // Atoms 0 and 1 belong to the diatom
//             a3 = 2;         // Atom 2 is the free atom
//         } else if (!is_p1_diatom && is_p2_diatom) {
//             sparta_atom = p1;
//             sparta_diatom = p2;
//             a3 = 0;         // Atom 0 is the free atom
//             a1 = 1; a2 = 2; // Atoms 1 and 2 belong to the diatom
//         } else {
//             DEBUG_PRINT("Only atom-diatom collisions are supported currently.");
//             return false;
//         }

//         // 2. Size arrays and map indices
//         collision.initial_state.atoms.resize(3);
//         collision.initial_state.diatoms.resize(1);
        
//         collision.initial_state.diatoms[0].atom1_index = a1;
//         collision.initial_state.diatoms[0].atom2_index = a2;

//         // 3. Set Masses
//         //double m_atom = Units::MASSES[sparta_atom->ispecies] * Units::AMU_TO_AUMASS;
//         //double m_diatom = Units::MASSES[sparta_diatom->ispecies] * Units::AMU_TO_AUMASS;
//         // Grab the molecular weight (in amu) directly from SPARTA's species array!
//         double m_atom = particle->species[sparta_atom->ispecies].molwt * Units::AMU_TO_AUMASS;
//         double m_diatom = particle->species[sparta_diatom->ispecies].molwt * Units::AMU_TO_AUMASS;

//         collision.initial_state.atoms[a1].mass = m_diatom / 2.0; 
//         collision.initial_state.atoms[a2].mass = m_diatom / 2.0; 
//         collision.initial_state.atoms[a3].mass = m_atom;
        
//         collision.mu_coll = (m_atom * m_diatom) / (m_atom + m_diatom);
//         collision.mu_dia = (m_diatom / 2.0) * (m_diatom / 2.0) / m_diatom;


//         std::string atom_name = particle->species[sparta_atom->ispecies].id;
//         std::string diatom_name = particle->species[sparta_diatom->ispecies].id;

//         // =======================================================
//         // NEW: 2. Populate the Config Bookkeeping
//         // =======================================================
//         config.species.clear();
//         config.species.push_back(p1_name);
//         config.species.push_back(p2_name);
//         config.diatom_name = diatom_name;
//         // MD parser to set the chemical system (e.g. "O3")
//         config.chemicalSystem = getChemicalSystemName(p1_name, p2_name); 

//         // =======================================================
//         // NEW: 3. Dynamically Extract Atomic Numbers
//         // =======================================================
//         // MD parser to extract the base element (e.g., "O2" -> "O" -> 8)
//         std::map<std::string, int> atom_elements = parseParticle(atom_name);
//         std::map<std::string, int> diatom_elements = parseParticle(diatom_name);
//         std::map<std::string, int> periodic_table = getPeriodicOrder();

//         // Safely grab the first element's atomic number
//         int z_atom = periodic_table[atom_elements.begin()->first];
//         int z_diatom = periodic_table[diatom_elements.begin()->first];

//         // Assign them to the states
//         collision.initial_state.atoms[a1].atomic_number = z_diatom;
//         collision.initial_state.atoms[a2].atomic_number = z_diatom;
//         collision.initial_state.atoms[a3].atomic_number = z_atom;
        
//         collision.initial_state.diatoms[0].atomic_numbers[0] = z_diatom;
//         collision.initial_state.diatoms[0].atomic_numbers[1] = z_diatom;



//         // 4. Extract Kinematics (Translational Energy)
//         std::vector<double> rel_vel(3, 0.0);
//         for (int i = 0; i < 3; i++) {
//             rel_vel[i] = (sparta_atom->v[i] - sparta_diatom->v[i]) * Units::MPS_TO_AUVEL;
//             collision.i_sparta_coll_axis[i] = rel_vel[i]; // Save raw axis for SPARTA update later
//         }
        
//         double vrel_sq = dotProduct(rel_vel, rel_vel);
//         collision.initial_state.etr = 0.5 * collision.mu_coll * vrel_sq;
        
//         // Normalize the saved SPARTA axis
//         double vrel_mag = sqrt(vrel_sq);
//         for (int i = 0; i < 3; i++) {
//             collision.i_sparta_coll_axis[i] /= vrel_mag;
//         }

//         // 5. Extract Internal States
//         collision.initial_state.diatoms[0].erot = sparta_diatom->erot * Units::J_TO_HARTREE;
//         collision.initial_state.diatoms[0].evib = sparta_diatom->evib * Units::J_TO_HARTREE; // Add ZPE here if your model requires it
//         collision.initial_state.diatoms[0].xj = sparta_diatom->xj;
//         collision.initial_state.diatoms[0].xv = sparta_diatom->xv;
        
//         // Propagate down to main molecular state view
//         collision.initial_state.erot = collision.initial_state.diatoms[0].erot;
//         collision.initial_state.evib = collision.initial_state.diatoms[0].evib;

//         // 6. Extract Electronic State
//         double eelec_atom = sparta_atom->eelec * Units::J_TO_HARTREE;
//         double eelec_diatom = sparta_diatom->eelec * Units::J_TO_HARTREE;
//         int ies;
//         MDIntegrator::getEState(ies, eelec_atom, eelec_diatom);
//         collision.initial_state.electronic_state = ies;
//         collision.integration.initial_electronic_state = ies;

//         // 7. Setup Integration Controls (Using config object)
//         collision.integration.max_time_fs = 50000.0;
//         collision.initial_separation = config.initial_separation * Units::ANG_TO_BOHR;
//         collision.integration.bond_check = config.bond_check_distance * Units::ANG_TO_BOHR;
//         collision.integration.integration_type = config.integration_type;

//         if (collision.integration.integration_type > 9) {
//             collision.integration.surface_hopping_enabled = true;
//             collision.integration.timestep_fs = std::isnan(config.timestep_fs) ? 0.005 : config.timestep_fs;
//         } else {
//             collision.integration.surface_hopping_enabled = false;
//             collision.integration.timestep_fs = std::isnan(config.timestep_fs) ? 0.05 : config.timestep_fs;
//         }

//         // 8. Generate Random Numbers for Geometry Setup
//         unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + config.rankNumber * 10;
//         std::mt19937 ic_rng(seed);
//         std::uniform_real_distribution<double> dist(0.0, 1.0);
        
//         collision.rand1 = dist(ic_rng); collision.rand2 = dist(ic_rng);
//         collision.rand3 = dist(ic_rng); collision.rand4 = dist(ic_rng);
//         collision.rand5 = dist(ic_rng); collision.rand6 = dist(ic_rng);
//         collision.rand7 = dist(ic_rng); collision.rand8 = dist(ic_rng);
//         collision.rand9 = dist(ic_rng); collision.rand10 = dist(ic_rng);

//         // 9. Sample Impact Parameter (Relies on config cross-sections)
//         bool sample_success = sampleImpactParameter(config, collision, 0);
//         if (!sample_success) {
//             DEBUG_PRINT("Failed to sample impact parameter!");
//             return false;
//         }

//         // 10. Perform Geometrical Placement!
//         bool setupSuccess = QCT::setupADCollision(config, collision, ic_rng);
        
//         return setupSuccess;

//     } catch (const std::exception& e) {
//         collision.error_msg = std::string("prepareFromParticles() exception: ") + e.what();
//         DEBUG_PRINT("Error: " << collision.error_msg);
//         return false;
//     }
// }


// bool CollideMD::prepareFromParticles(CRDSConfig& config, CollisionData& collision, 
//                                      Particle::OnePart* p1, Particle::OnePart* p2) {
    
//     collision.clear();
//     DEBUG_PRINT( "\n=== Preparing SPARTA particle collision ===" );

//     try {
//         // 1. Determine which SPARTA particle is the atom and which is the diatom
//         char* p1_name = particle->species[p1->ispecies].id;
//         char* p2_name = particle->species[p2->ispecies].id;
        
//         bool is_p1_diatom = (strlen(p1_name) >= 2);
//         bool is_p2_diatom = (strlen(p2_name) >= 2);

//         Particle::OnePart* sparta_atom = nullptr;
//         Particle::OnePart* sparta_diatom = nullptr;
//         int a1 = -1, a2 = -1, a3 = -1;

//         if (is_p1_diatom && !is_p2_diatom) {
//             sparta_diatom = p1;
//             sparta_atom = p2;
//             a1 = 0; a2 = 1; // Atoms 0 and 1 belong to the diatom
//             a3 = 2;         // Atom 2 is the free atom
//         } else if (!is_p1_diatom && is_p2_diatom) {
//             sparta_atom = p1;
//             sparta_diatom = p2;
//             a3 = 0;         // Atom 0 is the free atom
//             a1 = 1; a2 = 2; // Atoms 1 and 2 belong to the diatom
//         } else {
//             DEBUG_PRINT("Only atom-diatom collisions are supported currently.");
//             return false;
//         }

//         // 2. Size arrays and map indices
//         collision.initial_state.atoms.resize(3);
//         collision.initial_state.diatoms.resize(1);
        
//         collision.initial_state.diatoms[0].atom1_index = a1;
//         collision.initial_state.diatoms[0].atom2_index = a2;

//         // 3. Set Masses
//         double m_atom = Units::MASSES[sparta_atom->ispecies] * Units::AMU_TO_AUMASS;
//         double m_diatom = Units::MASSES[sparta_diatom->ispecies] * Units::AMU_TO_AUMASS;

//         collision.initial_state.atoms[a1].mass = m_diatom / 2.0; 
//         collision.initial_state.atoms[a2].mass = m_diatom / 2.0; 
//         collision.initial_state.atoms[a3].mass = m_atom;
        
//         collision.mu_coll = (m_atom * m_diatom) / (m_atom + m_diatom);
//         collision.mu_dia = (m_diatom / 2.0) * (m_diatom / 2.0) / m_diatom;

//         // 4. Extract Kinematics (Translational Energy)
//         std::vector<double> rel_vel(3, 0.0);
//         for (int i = 0; i < 3; i++) {
//             rel_vel[i] = (sparta_atom->v[i] - sparta_diatom->v[i]) * Units::MPS_TO_AUVEL;
//             collision.i_sparta_coll_axis[i] = rel_vel[i]; // Save raw axis for SPARTA update later
//         }
        
//         double vrel_sq = dotProduct(rel_vel, rel_vel);
//         collision.initial_state.etr = 0.5 * collision.mu_coll * vrel_sq;
        
//         // Normalize the saved SPARTA axis
//         double vrel_mag = sqrt(vrel_sq);
//         for (int i = 0; i < 3; i++) {
//             collision.i_sparta_coll_axis[i] /= vrel_mag;
//         }

//         // 5. Extract Internal States
//         collision.initial_state.diatoms[0].erot = sparta_diatom->erot * Units::J_TO_HARTREE;
//         collision.initial_state.diatoms[0].evib = sparta_diatom->evib * Units::J_TO_HARTREE; // Add ZPE here if your model requires it
//         collision.initial_state.diatoms[0].xj = sparta_diatom->xj;
//         collision.initial_state.diatoms[0].xv = sparta_diatom->xv;
        
//         // Propagate down to main molecular state view
//         collision.initial_state.erot = collision.initial_state.diatoms[0].erot;
//         collision.initial_state.evib = collision.initial_state.diatoms[0].evib;

//         // 6. Extract Electronic State
//         double eelec_atom = sparta_atom->eelec * Units::J_TO_HARTREE;
//         double eelec_diatom = sparta_diatom->eelec * Units::J_TO_HARTREE;
//         int ies;
//         MDIntegrator::getEState(ies, eelec_atom, eelec_diatom);
//         collision.initial_state.electronic_state = ies;
//         collision.integration.initial_electronic_state = ies;

//         // 7. Setup Integration Controls (Using config object)
//         collision.integration.max_time_fs = 50000.0;
//         collision.initial_separation = config.initial_separation * Units::ANG_TO_BOHR;
//         collision.integration.bond_check = config.bond_check_distance * Units::ANG_TO_BOHR;
//         collision.integration.integration_type = config.integration_type;

//         if (collision.integration.integration_type > 9) {
//             collision.integration.surface_hopping_enabled = true;
//             collision.integration.timestep_fs = std::isnan(config.timestep_fs) ? 0.005 : config.timestep_fs;
//         } else {
//             collision.integration.surface_hopping_enabled = false;
//             collision.integration.timestep_fs = std::isnan(config.timestep_fs) ? 0.05 : config.timestep_fs;
//         }

//         // 8. Generate Random Numbers for Geometry Setup
//         unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + config.rankNumber * 10;
//         std::mt19937 ic_rng(seed);
//         std::uniform_real_distribution<double> dist(0.0, 1.0);
        
//         collision.rand1 = dist(ic_rng); collision.rand2 = dist(ic_rng);
//         collision.rand3 = dist(ic_rng); collision.rand4 = dist(ic_rng);
//         collision.rand5 = dist(ic_rng); collision.rand6 = dist(ic_rng);
//         collision.rand7 = dist(ic_rng); collision.rand8 = dist(ic_rng);
//         collision.rand9 = dist(ic_rng); collision.rand10 = dist(ic_rng);

//         // 9. Sample Impact Parameter (Relies on config cross-sections)
//         bool sample_success = sampleImpactParameter(config, collision, 0);
//         if (!sample_success) {
//             DEBUG_PRINT("Failed to sample impact parameter!");
//             return false;
//         }

//         // 10. Perform Geometrical Placement!
//         bool setupSuccess = setupADCollision(config, collision, ic_rng);
        
//         return setupSuccess;

//     } catch (const std::exception& e) {
//         collision.error_msg = std::string("prepareFromParticles() exception: ") + e.what();
//         DEBUG_PRINT("Error: " << collision.error_msg);
//         return false;
//     }
// }


// bool CollideMD::prepareFromParticles(CRDSConfig& config, CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2) {

// 	#ifdef PAR_MD_BUILD
//     	std::ostream& out = OutputManager::get_output_stream();
//     #else
//     	std::ostream& out = OutputManager::get_output_stream();
//     #endif
// 	out << "\n=== Preparing collision from SPARTA particles() ===" << std::endl;
    
// 	collision.clear();
        
// 	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// 	//										                                          //
// 	//                                                                                //
// 	//PASS DATA IN FROM PARTICLES	                                                  //
// 	//                                                                                //
// 	//                                                                                //
// 	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    

// 	DEBUG_PRINT( "\n=== Testing prepare() ===" );
// 	// Determine collision type and number of atoms
// 	OutputManager::get_output_stream(OutputManager::StreamType::RANK_SPECIFIC) << "Particle 1 species: " << p1->ispecies << std::endl;
// 	OutputManager::get_output_stream(OutputManager::StreamType::RANK_SPECIFIC) << "Particle 2 species: " << p2->ispecies << std::endl;
	
// 	//int n_atoms = 0;
// 	//int n_diatoms = 0;
	   
// 	// if (p1->ispecies < 1000) {
// 	// 	// It's an atom
// 	// 	n_atoms += 1;
// 	// } else {
// 	// 	// It's a diatom
// 	// 	n_atoms += 2;
// 	// 	n_diatoms += 1;
// 	// }

// 	// // Count atoms in particle 2
// 	// if (p2->ispecies < 1000) {
// 	// 	// It's an atom
// 	// 	n_atoms += 1;
// 	// } else {
// 	// 	// It's a diatom
// 	// 	n_atoms += 2;
// 	// 	n_diatoms += 1;
// 	// }
		
// 	// Look up the species names from SPARTA
// 	char* p1_name = particle->species[p1->ispecies].id;
// 	char* p2_name = particle->species[p2->ispecies].id;

// 	bool is_p1_diatom = (strlen(p1_name) >= 2);
// 	bool is_p2_diatom = (strlen(p2_name) >= 2);

// 	int n_atoms = 0;
// 	int n_diatoms = 0;
// 	int a1 = -1, a2 = -1, a3 = -1, a4 = -1;

// 	// Use the booleans to set sizes and map the atom indices
// 	if (is_p1_diatom && !is_p2_diatom) {
// 		// p1 = Diatom, p2 = Atom
// 		n_atoms = 3; 
// 		n_diatoms = 1;
// 		a1 = 0; a2 = 1; // Atoms 0 and 1 belong to p1 (the diatom)
// 		a3 = 2;         // Atom 2 is p2 (the free atom)
// 	} 
// 	else if (!is_p1_diatom && is_p2_diatom) {
// 		// p1 = Atom, p2 = Diatom
// 		n_atoms = 3; 
// 		n_diatoms = 1;
// 		a3 = 0;         // Atom 0 is p1 (the free atom)
// 		a1 = 1; a2 = 2; // Atoms 1 and 2 belong to p2 (the diatom)
// 	} 
// 	else if (!is_p1_diatom && !is_p2_diatom) {
// 		// p1 = Atom, p2 = Atom
// 		n_atoms = 2; 
// 		n_diatoms = 0;
// 		a1 = 0; a2 = 1; 
// 	} 
// 	else if (is_p1_diatom && is_p2_diatom) {
// 		// p1 = Diatom, p2 = Diatom
// 		n_atoms = 4; 
// 		n_diatoms = 2;
// 		a1 = 0; a2 = 1; // Diatom 1
// 		a3 = 1; a4 = 2; // Diatom 2
// 	}

// 	collision.initial_state.atoms.resize(n_atoms);
// 	collision.initial_state.diatoms.resize(n_diatoms);

// 	// If there's a diatom, tell it which atoms it owns
// 	if (n_diatoms == 1) {
// 		collision.initial_state.diatoms[0].atom1_index = a1;
// 		collision.initial_state.diatoms[0].atom2_index = a2;
// 	}

// 		//Size arrays and set masses accordingly 
// 		collision.initial_state.atoms.resize(n_atoms);
		
// 	if (n_diatoms!=1){
// 		OutputManager::get_output_stream() << "Only atom-diatom collisions right now " << std::endl;
// 		return false;	    
// 	}else{
// 		collision.initial_state.diatoms.resize(n_diatoms);
// 	}

// 	double m1 = get_species_mass(p1->ispecies)*Units::AMU_TO_AUMASS;
// 	double m2 = get_species_mass(p2->ispecies)*Units::AMU_TO_AUMASS;
       
// 	for (int i =0; i<3;i++){
// 		collision.initial_state.atoms[i].mass = m1;
// 	}
// 	collision.initial_state.diatoms[0].mass = m2;

// 	//Set up indices for placing particles in arrays based on initial arrangement, always 0 for PESs with three identical atoms
// 	//collision.initial_arrangement=0;
// 	//int a1,a2,a3;
// 	//MDIntegrator::getAtomInds(a2,a2,a3,collision.initial_arrangement);

// 	//int diatom_ind;

// 	OutputManager::get_output_stream() << "Particle 1 mass: " << collision.initial_state.atoms[a3].mass << std::endl;
//     OutputManager::get_output_stream() << "Particle 2 mass: " << collision.initial_state.diatoms[0].mass << std::endl;
//     OutputManager::get_output_stream() << "Total atoms (including those in diatoms): " << n_atoms << std::endl;
//     OutputManager::get_output_stream() << "Total diatoms: " << n_diatoms << std::endl;
            
// 	// Determine relative momenta of sparta particles
// 	std::vector<double> rel_vel(3,0.0);
// 	std::vector<double> rel_mom(3,0.0);
// 	std::vector<double> sparta_pair_eelec(2,0.0);
// 	double sparta_erot, sparta_xj; 
// 	double sparta_evib, sparta_xv; 


// 	// Use local sparta particle objects to pass in data
// 	//const CollidingParticle* sparta_atom,* sparta_diatom;
// 	Particle::OnePart* sparta_atom,* sparta_diatom;

// 	if(p1->ispecies == 1 && p2->ispecies != 1){
// 		sparta_atom = p2;
// 		sparta_diatom = p1;
// 	}else if (p1->ispecies != 1 && p2->ispecies == 1){
// 		sparta_atom = p1;
// 		sparta_diatom = p2;
// 	}
		
// 		for (int i=0;i<3;i++){
// 				DEBUG_PRINT( "P1V P2V " << p1->v[i] << " " <<  p2->v[i] );
// 			rel_vel[i]= sparta_atom->v[i] - sparta_diatom->v[i];

// 				DEBUG_PRINT( "REL VEL " << i << " " <<  rel_vel[i] );
// 			rel_vel[i] *= Units::MPS_TO_AUVEL;
// 				DEBUG_PRINT( "REL VEL " << i << " " <<  rel_vel[i] );


// 			sparta_pair_eelec[0] = sparta_atom->eelec;

// 			sparta_erot  = sparta_diatom->erot;
// 			sparta_evib  = sparta_diatom->evib;
// 			sparta_xj  = sparta_diatom->xj;
// 			sparta_xv  = sparta_diatom->xv;
// 			sparta_pair_eelec[1]  = sparta_diatom->eelec;			


// 		}

// 		DEBUG_PRINT( "0.5									 Surface Hopping  : " << collision.integration.surface_hopping_enabled );


// 	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// 	//										  //
// 	//                                                                                //
// 	//SETUP COLLISIONDATA STRUCTURE WITH DATA FROM PARTICLES	                  //
// 	//                                                                                //
// 	//                                                                                //
// 	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


// 	try {
// 			//Collision reduced mass
// 			collision.mu_coll = m1*m2/(m1+m2);
// 			collision.mu_dia = collision.initial_state.atoms[a1].mass*collision.initial_state.atoms[a2].mass/(collision.initial_state.atoms[a1].mass+collision.initial_state.atoms[a2].mass);
			
// 			// Store original collision axis in SPARTA Cartesian coordinates
// 			for(int j=0;j<3;j++){
// 				rel_mom[j] = rel_vel[j]*collision.mu_coll;
// 					DEBUG_PRINT( "REL MOM " << j << " " <<  rel_mom[j] );
// 				collision.i_sparta_coll_axis[j] = rel_vel[j]/sqrt(dotProduct(rel_vel,rel_vel));
// 			}

// 				DEBUG_PRINT( "REL MOM " <<  sqrt(dotProduct(rel_mom,rel_mom)) );
// 				DEBUG_PRINT( "SPARTA COLLISION AXIS " << collision.i_sparta_coll_axis[0] << " " << collision.i_sparta_coll_axis[1]<< " " << collision.i_sparta_coll_axis[2] );
// 				DEBUG_PRINT( "EROT " <<     sparta_erot );
// 				DEBUG_PRINT( "EVIB " <<     sparta_evib );
// 				DEBUG_PRINT( "EELEC " <<    sparta_pair_eelec[0] << " " <<  sparta_pair_eelec[1]  );
			
			
			
// 			DEBUG_PRINT( "COLL MU " << collision.mu_coll );
			

// 			collision.initial_state.etr = dotProduct(rel_mom,rel_mom)/(2*collision.mu_coll);
				
// 			DEBUG_PRINT( "COLL EREL_TR " <<  collision.initial_state.etr );
			
			
// 			// Determine initial electronic state (Surface Number) 
// 			int ies;
// 			MDIntegrator::getEState(ies, collision.initial_state.atoms[a3].eelec , collision.initial_state.diatoms[0].eelec);
// 			collision.initial_state.electronic_state = ies;
		
// 				collision.initial_state.electronic_state = ies ; 
// 				collision.final_state.electronic_state = collision.initial_state.electronic_state ; 
			
// 			DEBUG_PRINT( "Initial electronic state is : " << ies );
// 			DEBUG_PRINT( "1									 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
			
			
// 			collision.initial_state.erot = sparta_erot*Units::J_TO_HARTREE;

// 			std::string diatom_name = get_species_name(sparta_diatom->ispecies);
// 			OutputManager::get_output_stream() <<  "DIATOM NAME : " << diatom_name  << std::endl;
// 			double zpe = getZPE(diatom_name, ies);
			
// 			collision.initial_state.evib = (sparta_evib*Units::J_TO_HARTREE)+zpe;

// 			collision.initial_state.diatoms[0].erot = sparta_erot*Units::J_TO_HARTREE;
// 			collision.initial_state.diatoms[0].evib = sparta_evib*Units::J_TO_HARTREE;

// 			collision.initial_state.diatoms[0].xj = sparta_xj;
// 			collision.initial_state.diatoms[0].xv = sparta_xv;

// 			collision.initial_state.atoms[a3].eelec = sparta_pair_eelec[0]*Units::J_TO_HARTREE;
// 			collision.initial_state.diatoms[0].eelec = sparta_pair_eelec[1]*Units::J_TO_HARTREE;
// 			collision.initial_state.eelec =  collision.initial_state.atoms[a3].eelec + collision.initial_state.diatoms[0].eelec;	
			
			
// 			// collision.integration.max_time_fs = 50000.0;
// 			// collision.impactParameter = 0.0; // 0 Angstrom
// 			// collision.initial_separation = config.initial_separation*Units::ANG_TO_BOHR; // 
// 			// collision.bond_check = config.bond_check_distance*Units::ANG_TO_BOHR; //     
// 			// collision.integration.integration_type = config.integration_type; 
// 			// collision.n_states=1;
			
				
// 			//     // Set default parameters for collision data object 
// 			// if(!config.adiabatic){
// 			// 	collision.integration.surface_hopping_enabled = true;
// 			// 	//integrator.enableSurfaceHopping(true);
// 			// 	collision.integration.timestep_fs = 0.005;
// 			//     DEBUG_PRINT( "2 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
// 			// }else{
// 			// 	collision.integration.timestep_fs = 0.01;
// 			// }

// 			// Set parameters for collision data object 
// 			collision.integration.max_time_fs = 50000.0;
// 			collision.intergration.initial_separation = config.initial_separation*Units::ANG_TO_BOHR; // ~15 Angstrom
// 			collision.integration.bond_check = config.bond_check_distance*Units::ANG_TO_BOHR; // ~10 Angstrom
				
// 			collision.integration.integration_type = config.integration_type;

// 			// Nonadiabatic dynamics (integration_type > 9) requires smaller timestep
// 			if(collision.integration.integration_type > 9){
// 				collision.integration.surface_hopping_enabled = true;
// 				// Use user-specified timestep if provided, otherwise default to 0.005 fs
// 				if (std::isnan(config.timestep_fs)) {
// 					collision.integration.timestep_fs = 0.005;
// 				} else {
// 					collision.integration.timestep_fs = config.timestep_fs;
// 				}
// 				DEBUG_PRINT( "2 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
// 				OutputManager::get_output_stream(OutputManager::StreamType::OVERALL, 2) << "Timestep: " << collision.integration.timestep_fs << " fs" << std::endl;
// 			}else{
// 				collision.integration.surface_hopping_enabled = false;
// 				// Use user-specified timestep if provided, otherwise default to 0.05 fs
// 				if (std::isnan(config.timestep_fs)) {
// 					collision.integration.timestep_fs = 0.05;
// 				} else {
// 					collision.integration.timestep_fs = config.timestep_fs;
// 				}
// 			}



// 			DEBUG_PRINT( "3 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
				
				
			
// 			// Set up atoms based on collision type
// 			bool p1_is_atom = (p1->ispecies < 1000);
// 			bool p2_is_atom = (p2->ispecies < 1000);

// 			if (p1_is_atom && p2_is_atom) {
// 				// Atom + Atom collision
// 				//setupAtomAtomCollision(collision, p1, p2);
// 				OutputManager::get_output_stream() << "ATOM ATOM NOT INITIAILZED YET" << std::endl;
// 			} else if ((p1_is_atom && !p2_is_atom) || (!p1_is_atom && p2_is_atom)) {
// 				// Atom + Diatom collision
// 				setupAtomDiatomCollision(config, collision);
// 			} else if (!p1_is_atom && !p2_is_atom) {
// 				// Diatom + Diatom collision
// 				//setupDiatomDiatomCollision(collision, p1, p2);
// 				OutputManager::get_output_stream() << "DIATOM DIATOM NOT INITIAILZED YET" << std::endl;
// 			}



// 				// if (p1->ispecies == 0 && p2->ispecies == 0) {
// 				//     // O + O collision
// 				//     setupOOCollision(collision, p1, p2);
// 				// } else if ((p1->ispecies == 0 && p2->ispecies == 1) || 
// 				//            (p1->ispecies == 1 && p2->ispecies == 0)) {
// 				//     // O2 + O collision
// 				//     setupO2OCollision(collision, p1, p2);
// 				// } else if (p1->ispecies == 1 && p2->ispecies == 1) {
// 				//     // O2 + O2 collision
// 				//     setupO2O2Collision(collision, p1, p2);
// 				// }
				

// 			DEBUG_PRINT( "INITIAL THETA/PHI " << collision.imd_ca_theta << " " << collision.imd_ca_phi );


// 			DEBUG_PRINT( "4 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
// 				OutputManager::get_output_stream(OutputManager::StreamType::RANK_SPECIFIC) << "Collision prepared successfully" << std::endl;
				
// 	#if !defined(SPARTA_MD_BUILD) && !defined(PAR_MD_BUILD) 
// 				int goOn;
// 				OutputManager::get_output_stream() << "\n -1 to exit here...";
// 				std::cin >> goOn;
// 			if(goOn==-1){
// 				exit(0);
// 			}
// 	#endif	
// 			collision.success = true;
// 				return true;
					
// 		} catch (const std::exception& e) {
// 					collision.error_msg = std::string("prepare() exception: ") + e.what();
// 					OutputManager::get_output_stream() << "Error: " << collision.error_msg << std::endl;
// 					return false;
// 		}




// 	}



// 	// 		DEBUG_PRINT( "\n=== Testing prepare() ===" );
// 	// 	    // Determine collision type and number of atoms
// 	//             std::cout << "Particle 1 species: " << p1->ispecies << std::endl;
// 	//             std::cout << "Particle 2 species: " << p2->ispecies << std::endl;
		
// 	//             int n_atoms = 0;
// 	//             int n_diatoms = 0;
// 	// 	    collision.p1Type = p1->ispecies;
// 	// 	    collision.p2Type = p2->ispecies;
				
// 	// 	    if (p1->ispecies == 0){ n_atoms += 1;      // O atom
// 	// 	    }else if (p1->ispecies == 1){ n_atoms += 2; // O2 molecule
// 	// 	    	n_diatoms += 1;
// 	// 	    }


// 	//             if (p2->ispecies == 0){ n_atoms += 1;      // O atom  
// 	// 	    }else if (p2->ispecies == 1) {n_atoms += 2; // O2 molecule
// 	// 	    	n_diatoms += 1; 
// 	// 	    }

// 	// 	//Size arrays and set masses accordingly 
// 	// 	    collision.initial_state.atoms.resize(n_atoms);
				
// 	// 	    if (n_diatoms!=1){
// 	// 		    std::cout << "Only atom-diatom collisions right now " << std::endl;
// 	// 	   	    return false;	    
// 	// 	    }else{

// 	//             collision.initial_state.diatoms.resize(n_diatoms);
// 	// 	    }



// 	// 	double m1 = Units::MASSES[p1->ispecies]*Units::AMU_TO_AUMASS;
// 	// 	double m2 = Units::MASSES[p2->ispecies]*Units::AMU_TO_AUMASS;
		
// 	// 	for (int i =0; i<3;i++){
// 	// 		collision.initial_state.atoms[i].mass = m1;
// 	// 	}
// 	// 	collision.initial_state.diatoms[0].mass = m2;

// 	// 	//Set up indices for placing particles in arrays based on initial arrangement, always 0 for PESs with three identical atoms
// 	// 	collision.initial_arrangement=0;
// 	// 	int a1,a2,a3;
// 	// 	MDIntegrator::getAtomInds(a2,a2,a3,collision.initial_arrangement);

// 	// 	int diatom_ind;


// 	// 	std::cout << "Particle 1 mass: " << collision.initial_state.atoms[a3].mass << std::endl;
// 	//         std::cout << "Particle 2 mass: " << collision.initial_state.diatoms[0].mass << std::endl;


// 	//             std::cout << "Total atoms (including those in diatoms): " << n_atoms << std::endl;
// 	//             std::cout << "Total diatoms: " << n_diatoms << std::endl;
				
// 	// // Determine relative momenta of sparta particles
// 	// 	std::vector<double> rel_vel(3,0.0);
// 	// 	std::vector<double> rel_mom(3,0.0);
// 	// 	std::vector<double> sparta_pair_eelec(2,0.0);
// 	// 	double sparta_erot, sparta_xj; 
// 	// 	double sparta_evib, sparta_xv; 

// 	// 	// Use local sparta particle objects to pass in data
// 	// 	Particle::OnePart* sparta_atom,* sparta_diatom;

// 	// 	if(p1->ispecies == 1 && p2->ispecies != 1){
// 	// 		sparta_atom = p2;
// 	// 		sparta_diatom = p1;
// 	// 	}else if (p1->ispecies != 1 && p2->ispecies == 1){
// 	// 		sparta_atom = p1;
// 	// 		sparta_diatom = p2;
// 	// 	}
			
// 	// 		for (int i=0;i<3;i++){
// 	//      			DEBUG_PRINT( "P1V P2V " << p1->v[i] << " " <<  p2->v[i] );
// 	// 			rel_vel[i]= sparta_atom->v[i] - sparta_diatom->v[i];

// 	//      			DEBUG_PRINT( "REL VEL " << i << " " <<  rel_vel[i] );
// 	// 			rel_vel[i] *= Units::MPS_TO_AUVEL;
// 	//      			DEBUG_PRINT( "REL VEL " << i << " " <<  rel_vel[i] );


// 	// 			sparta_pair_eelec[0] = sparta_atom->eelec;

// 	// 			sparta_erot  = sparta_diatom->erot;
// 	// 			sparta_evib  = sparta_diatom->evib;
// 	// 			sparta_xj  = sparta_diatom->xj;
// 	// 			sparta_xv  = sparta_diatom->xv;
// 	// 			sparta_pair_eelec[1]  = sparta_diatom->eelec;			


// 	// 		}

// 	// 		DEBUG_PRINT( "0.5									 Surface Hopping  : " << collision.integration.surface_hopping_enabled );

// 	// //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// 	// //										  //
// 	// //                                                                                //
// 	// //SETUP COLLISIONDATA STRUCTURE WITH DATA FROM PARTICLES	                  //
// 	// //                                                                                //
// 	// //                                                                                //
// 	// //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// 	// 	try {
// 	// 		//Collision reduced mass
// 	// 		collision.mu_coll = m1*m2/(m1+m2);
// 	// 		collision.mu_dia = collision.initial_state.atoms[a1].mass*collision.initial_state.atoms[a2].mass/(collision.initial_state.atoms[a1].mass+collision.initial_state.atoms[a2].mass);
			
// 	// 		// Store original collision axis in SPARTA Cartesian coordinates
// 	// 		for(int j=0;j<3;j++){
// 	// 			rel_mom[j] = rel_vel[j]*collision.mu_coll;
// 	//      			DEBUG_PRINT( "REL MOM " << j << " " <<  rel_mom[j] );
// 	// 			collision.i_sparta_coll_axis[j] = rel_vel[j]/sqrt(dotProduct(rel_vel,rel_vel));
// 	// 		}

// 	// 	     	DEBUG_PRINT( "REL MOM " <<  sqrt(dotProduct(rel_mom,rel_mom)) );
// 	// 	     	DEBUG_PRINT( "SPARTA COLLISION AXIS " << collision.i_sparta_coll_axis[0] << " " << collision.i_sparta_coll_axis[1]<< " " << collision.i_sparta_coll_axis[2] );
// 	// 		DEBUG_PRINT( "EROT " <<     sparta_erot );
// 	// 	     	DEBUG_PRINT( "EVIB " <<     sparta_evib );
// 	// 	     	DEBUG_PRINT( "EELEC " <<    sparta_pair_eelec[0] << " " <<  sparta_pair_eelec[1]  );
			
// 	// 		collision.etr = dotProduct(rel_mom,rel_mom)/(2*collision.mu_coll);
			
// 	// 		DEBUG_PRINT( "COLL MU " << collision.mu_coll );
			
// 	// 		collision.etrI = collision.etr; 	
				
// 	// 		DEBUG_PRINT( "COLL EREL_TR " <<  collision.etrI );
			
			
// 	// 		// Determine initial electronic state (Surface Number) 
// 	// 		int ies;
// 	// 		MDIntegrator::getEState(ies, collision.initial_state.atoms[a3].eelec , collision.initial_state.diatoms[0].eelec);
// 	// 		collision.initial_state.electronic_state = ies;
		
// 	// 	        collision.initial_state.electronic_state = ies ; 
// 	// 	        collision.final_state.electronic_state = collision.initial_state.electronic_state ; 
			
// 	// 		DEBUG_PRINT( "Initial electronic state is : " << ies );
// 	// 		DEBUG_PRINT( "1									 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
			
			
// 	// 		collision.erot = sparta_erot*Units::J_TO_HARTREE;
// 	// 		collision.erotI = collision.erot;
			
// 	// 		collision.evib = (sparta_evib*Units::J_TO_HARTREE)+ZPE[ies];
// 	// 		//collision.evib = (sparta_evib*Units::J_TO_HARTREE);
// 	// 		collision.evibI = collision.evib;

// 	// 		collision.initial_state.diatoms[0].erot = sparta_erot*Units::J_TO_HARTREE;
// 	// 		collision.initial_state.diatoms[0].evib = sparta_evib*Units::J_TO_HARTREE;

// 	// 		collision.initial_state.diatoms[0].xj = sparta_xj;
// 	// 		collision.initial_state.diatoms[0].xv = sparta_xv;

// 	// 		collision.initial_state.atoms[a3].eelec = sparta_pair_eelec[0]*Units::J_TO_HARTREE;
// 	// 		collision.initial_state.diatoms[0].eelec = sparta_pair_eelec[1]*Units::J_TO_HARTREE;
// 	// 		collision.eelecI =  collision.initial_state.atoms[a3].eelec + collision.initial_state.diatoms[0].eelec;	
			
			
			
// 	// 		collision.integration.max_time_fs = 50000.0;
// 	// 	        collision.impactParameter = 0.0; // 0 Angstrom
// 	// 	        collision.initial_separation = 15.0*Units::ANG_TO_BOHR; // ~15 Angstrom
// 	// 	        collision.bond_check = 12.0*Units::ANG_TO_BOHR; // ~10 Angstrom
				
// 	// 	   	collision.integration.integration_type = integration_type; 
			
// 	// 		collision.n_states=6;
			
			
			
// 	// 		integrator.initializePES();
// 	// 		integrator.setInitialSurface(ies);
		
					
// 	// 	        // Set default parameters for collision data object 
// 	// 	        if(collision.integration.integration_type==10){
// 	// 			collision.integration.surface_hopping_enabled = true;
// 	// 			integrator.enableSurfaceHopping(true);
// 	// 			collision.integration.timestep_fs = 0.005;
// 	// 		DEBUG_PRINT( "2 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
// 	// 		}else{
// 	// 			collision.integration.timestep_fs = 0.05;
// 	// 		}

// 	// 		DEBUG_PRINT( "3 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
				
				
			
			
			
			
// 	// 		// Set up atoms based on collision type
// 	// 	        if (p1->ispecies == 0 && p2->ispecies == 0) {
// 	// 	            // O + O collision
// 	// 	            setupOOCollision(collision, p1, p2);
// 	// 	        } else if ((p1->ispecies == 0 && p2->ispecies == 1) || 
// 	// 	                   (p1->ispecies == 1 && p2->ispecies == 0)) {
// 	// 	            // O2 + O collision
// 	// 	            setupO2OCollision(collision, p1, p2);
// 	// 	        } else if (p1->ispecies == 1 && p2->ispecies == 1) {
// 	// 	            // O2 + O2 collision
// 	// 	            setupO2O2Collision(collision, p1, p2);
// 	// 	        }
				

// 	// 		DEBUG_PRINT( "INITIAL THETA/PHI " << collision.imd_ca_theta << " " << collision.imd_ca_phi );


// 	// 		DEBUG_PRINT( "4 Surface Hopping  : " << collision.integration.surface_hopping_enabled );
// 	// 	        std::cout << "Collision prepared successfully" << std::endl;
				
// 	// #ifndef SPARTA_MD_BUILD	
// 	//     		int goOn;
// 	//     		std::cout << "\n -1 to exit here...";
// 	//     		std::cin >> goOn;
// 	// 		if(goOn==-1){
// 	// 			exit(0);
// 	// 		}
// 	// #endif	
// 	// 		collision.success = true;
// 	// 	        return true;
					
// 	// 	} catch (const std::exception& e) {
// 	// 	            collision.error_msg = std::string("prepare() exception: ") + e.what();
// 	// 	            std::cout << "Error: " << collision.error_msg << std::endl;
// 	// 	            return false;
// 	// 	}
// 	// }






/// COLLIDE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// 	bool CollideMD::collide(MDIntegrator& integrator , CollisionData& collision) {
//         std::cout << "\n=== Colliding() ===" << std::endl;
	
// 		// Create a basic config object based on your new IntegrationControls
// 		CRDSConfig config;
// 		config.max_time_fs = 50000.0;


// 		if(true) {
// 				std::cout << "\n VELOCITY VERLET" << std::endl;
// 			return integrator.integrateTrajectory(collision);
// 		}else if(false){
// 				std::cout << "\n BURLISCH STOER" << std::endl;
// 			return integrator.integrateTrajectoryBS(collision);
// 			}	
// 		return false;
// }


bool CollideMD::process(MDIntegrator& integrator , CollisionData& collision) {
        DEBUG_PRINT( "\n=== Testing process() ===" );
        
        try {
            // Simple outcome analysis for testing
            //collision.reaction_outcome = analyzeOutcome(collision);
            //std::cout << "Reaction outcome: " << collision.reaction_outcome << std::endl;
            	  
	    // Calculate energy change
            collision.energy_change = collision.final_state.potential_energy + collision.final_state.kinetic_energy -
                                    collision.initial_state.potential_energy - collision.initial_state.kinetic_energy;

            std::cout << "Energy change: " << collision.energy_change << " Hartree " << std::endl;
            
	    std::cout << "Initial total, PE, KE: " << collision.initial_state.potential_energy + collision.initial_state.kinetic_energy << " " << collision.initial_state.potential_energy << " " << collision.initial_state.kinetic_energy << std::endl;
	    std::cout << "Final total, PE, KE: " << collision.final_state.potential_energy + collision.final_state.kinetic_energy << " " << collision.final_state.potential_energy << " " << collision.final_state.kinetic_energy << std::endl;

	    collision.reaction_outcome = integrator.processOutcome(collision);

	    DEBUG_PRINT( " FINAL DIATOM MOM" <<  sqrt(dotProduct(collision.final_state.diatoms[0].mom,collision.final_state.diatoms[0].mom)) );

            return true;
            
        } catch (const std::exception& e) {
            collision.error_msg = std::string("process() exception: ") + e.what();
            std::cout << "Error: " << collision.error_msg << std::endl;
            return false;
        }
}



// void CollideMD::setupOOCollision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2) {
//         std::cout << "Setting up O + O collision" << std::endl;
        
//         // Atom 1 (O)
//         collision.initial_state.atoms[0].pos[0] = -7.5; // -7.5 Angstrom
//         collision.initial_state.atoms[0].pos[1] = 0.0;
//         collision.initial_state.atoms[0].pos[2] = 0.0;
//         collision.initial_state.atoms[0].mom[0] = 10.0; // 1000 m/s
//         collision.initial_state.atoms[0].mom[1] = 0.0;
//         collision.initial_state.atoms[0].mom[2] = 0.0;
//         collision.initial_state.atoms[0].mass = 15.999*Units::AMU_TO_AUMASS; // O atomic mass
//         collision.initial_state.atoms[0].atomic_number = 8;
        
//         // Atom 2 (O)
//         collision.initial_state.atoms[1].pos[0] = 7.5;  // +7.5 Angstrom
//         collision.initial_state.atoms[1].pos[1] = collision.impactParameter;
//         collision.initial_state.atoms[1].pos[2] = 0.0;
//         collision.initial_state.atoms[1].mom[0] = -10.0; // -1000 m/s
//         collision.initial_state.atoms[1].mom[1] = 0.0;
//         collision.initial_state.atoms[1].mom[2] = 0.0;
//         collision.initial_state.atoms[1].mass = 15.999*Units::AMU_TO_AUMASS;
//         collision.initial_state.atoms[1].atomic_number = 8;
//     }
    


// 	void CollideMD::setupO2OCollision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2) {

//  // Set up basics into atom data containers
// 	std::cout << "Setting up O2 + O collision" << std::endl;
		
//         collision.initial_state.atoms[0].mass = 15.994915*Units::AMU_TO_AUMASS;
//         collision.initial_state.atoms[0].atomic_number = 8;

//         collision.initial_state.atoms[1].mass = 15.994915*Units::AMU_TO_AUMASS;
//         collision.initial_state.atoms[1].atomic_number = 8;

//         collision.initial_state.atoms[2].mass = 15.994915*Units::AMU_TO_AUMASS;
//         collision.initial_state.atoms[2].atomic_number = 8;
       
// 	collision.initial_arrangement = 0 ;
// 	int a1ind, a2ind, a3ind; 
// 	MDIntegrator::getAtomInds(a1ind,a2ind,a3ind, collision.initial_arrangement);
	
//         double diatom_mu;
//         diatom_mu = (collision.initial_state.atoms[a1ind].mass*collision.initial_state.atoms[a2ind].mass)/(collision.initial_state.atoms[a1ind].mass + collision.initial_state.atoms[a2ind].mass);
//         std::vector<double> mass_scale(2);
//         mass_scale[0]= diatom_mu/collision.initial_state.atoms[a1ind].mass;
//         mass_scale[1]= diatom_mu/collision.initial_state.atoms[a2ind].mass;


//         int sample_mode;
// 	#ifdef SPARTA_MD_BUILD
// 	// Get info from SPARTA particles
// 	sample_mode = 1; 

// 	#else
// 	// Get sampling mode from user 
// 			std::cout << "\nEnter collision init type \n"<< std::endl;
// 		std::cout << "-2 for test against benchmark \n -1 for init from jv file with contrived orientation \n 0 for normal init from jv file \n 1 for input from SPARTA Particle  ";
// 			std::cin >> sample_mode;
// 		if(sample_mode==-3){
// 			exit(0);
// 		};
// 	#endif  
// 	//Initialize variables for later
// 	#ifdef FIXED_RNG
// 		std::mt19937 ic_rng(7);
// 	#else
// 		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
// 		std::cout << "RANDOM SEED 1 " << seed <<std::endl;
// 		std::mt19937 ic_rng(seed); 
// 	#endif
// 		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		
// 	double rand1 = distribution(ic_rng);
// 	double rand2 = distribution(ic_rng);
// 	double rand3 = distribution(ic_rng);
// 	double rand4 = distribution(ic_rng);
// 	double rand5 = distribution(ic_rng);
// 	double rand6 = distribution(ic_rng); 
// 	double rand7 = distribution(ic_rng);  
// 	double rand8 = distribution(ic_rng);  
	
// 	std::cout << "RANDOMS: " << rand1 << " " << rand2 << " " << rand3 << " " << rand4 << std::endl;
// 	std::cout << "RANDOMS: " << rand5 << " " << rand6 << " " << rand7 << " " << rand8 << std::endl;








// 		std::vector<jvData> jvdata;

// 		int chosen_state; 
// 		double vib_phase; 
// 		double drift_time; 
// 		double bond_length; //;jvdata[chosen_state].rin + 0.5*(jvdata[chosen_state].rout - jvdata[chosen_state].rin);
// 		double tau, r_in, r_out;
		
// 		double bond_theta;
// 		double bond_phi; 

// 		double vib_mom; 
// 		double rot_mom;

// 		double am_theta;
// 		double am_phi;
// 		int state_num = collision.initial_state.electronic_state ;
		
// 		if(sample_mode==-2){

// 			readICsFromFile("ICs.txt", collision);

// 		}else if (sample_mode==-1  || sample_mode==0){  //Choose JV state read from file 
				
			
// 			for  (int e_state = 0; e_state < collision.n_states; e_state++){
// 				std::string JVFile = "JV_levels_O2-"+std::to_string(e_state)+".txt";	
			
// 				std::ifstream JVFileStream(JVFile);

// 					if (JVFileStream.good()) {
// 						std::cout << "The file '" << JVFile  << "' exists." << std::endl;
// 					jvdata = read_jv_levels_from_file( JVFile);
// 						std::cout << "Read " << jvdata.size() << " J,V levels for O2 State " << std::to_string(e_state) << std::endl;
// 					} else {
// 						std::cout << "The file '"  << JVFile << "' does not exist." << std::endl;

// 					jvdata = calculate_all_jv_levels("O2", 300, 70,e_state);
// 						std::cout << "Found " << jvdata.size() << " J,V levels for O2 State " << std::to_string(e_state) << std::endl;
// 					// Write to file
// 						write_jv_levels_to_file(jvdata,"O2", e_state);

// 				}	

// 			}



// 			std::string JVFile = "JV_levels_O2-"+std::to_string(collision.initial_state.electronic_state)+".txt";	
// 			jvdata = read_jv_levels_from_file( JVFile );
// 				std::cout << "Read " << jvdata.size() << " J,V levels for O2 State " << std::to_string(state_num) << std::endl;



// 				std::cout << "\nEnter a JV state index...";
// 				std::cin >> chosen_state;
// 				chosen_state--;
		
			
// 				// Set rotational and vibrational energy  
// 				vib_phase = jvdata[chosen_state].tau; 
// 				bond_length = jvdata[chosen_state].rin; //+ 0.5*(jvdata[chosen_state].rout - jvdata[chosen_state].rin);
				
// 				vib_mom = 0.0; //jvdata[chosen_state].rin + rand*(jvdata[chosen_state].rout - jvdata[chosen_state].rin);
// 				double pmag; 
// 				double XJ = static_cast<double>(jvdata[chosen_state].j);
				
// 				calcXJ(rot_mom ,XJ ,-1);
// 				rot_mom = sqrt(rot_mom);
				
// 				//rot_mom =  sqrt(jvdata[chosen_state].erot* 2.0 * diatom_mu * pow(bond_length,2));
	
// 				std::cout << "Chosen State " << chosen_state << std::endl; 
// 				std::cout << "J " << jvdata[chosen_state].j << " V " << jvdata[chosen_state].v << std::endl; 
// 				std::cout << "edia " << jvdata[chosen_state].energy << std::endl;
// 				std::cout << "erot " << jvdata[chosen_state].erot << " evib " << jvdata[chosen_state].evib << std::endl; 
				
// 				collision.initial_state.diatoms[0].xj = static_cast<double>(jvdata[chosen_state].j);
// 				collision.initial_state.diatoms[0].xv = static_cast<double>(jvdata[chosen_state].v);


// 				//collision.initial_state.diatoms[0].erot = jvdata[chosen_state].erot;
// 				//collision.initial_state.diatoms[0].evib = ;


// 				auto [edia_final, r_in, r_out] = calc_edia(collision.initial_state.diatoms[0].xj, collision.initial_state.diatoms[0].xv, collision.initial_state.electronic_state);
			
// 					//auto [rmin, rmax, emin, emax] = MM_dia(collision.initial_state.diatoms[0].xj, collision.initial_state.electronic_state);
// 					//auto [r_in, r_out] = turning_points(jvdata[chosen_state].energy, collision.initial_state.diatoms[0].xj, rmin, emin, rmax, emax, collision.initial_state.electronic_state);
	




// 				//auto [dumv, r_in, r_out] = calc_v(jvdata[chosen_state].energy, collision.initial_state.diatoms[0].xj, 0);
// 				//std::cout << "dum v " << dumv << std::endl;
// 				std::cout<< " rin, rin2 " << jvdata[chosen_state].rin << " " << r_in << std::endl; 
// 				DEBUG_PRINT( "rout, rout2 " <<  jvdata[chosen_state].rout << " " << r_out );
				
// 				tau = calc_tau(edia_final, r_in, r_out, collision.initial_state.diatoms[0].xj, collision.initial_state.electronic_state);	
				
// 				DEBUG_PRINT( "Edia, Edia2 " << jvdata[chosen_state].energy << " " << edia_final ); 
// 				DEBUG_PRINT( "TAU1, Tau2 " << jvdata[chosen_state].tau << " " << tau ); 
// 				std::cout<< " rin, rin2 " << jvdata[chosen_state].rin << " " << r_in << std::endl; 
// 				DEBUG_PRINT( "rout, rout2 " <<  jvdata[chosen_state].rout << " " << r_out );
				
// 				//rand5 = distribution(ic_rng);
// 				DEBUG_PRINT( "RAND 5: " << rand5 ); 
// 				double test_time = rand5*tau*Units::AUTIME_TO_FS;			
// 				if(test_time >= 0.50*tau){
// 					bond_length = r_out;
// 					drift_time = test_time - 0.50*tau;
// 				}else if (test_time < 0.50*tau){
// 					bond_length = r_in;
// 					drift_time = test_time;
// 				}
// 				std::cout << "DRIFT TIME " << drift_time << std::endl;
			


// 		}else if (sample_mode==1){
// 			double edia = collision.initial_state.diatoms[0].evib + collision.initial_state.diatoms[0].erot;

// 			DEBUG_PRINT( "                                       XJ,XV IN PREPARE: " << collision.initial_state.diatoms[0].xj << " " << collision.initial_state.diatoms[0].xv );
// 			auto [edia_jv, r_in, r_out] = calc_edia(collision.initial_state.diatoms[0].xj, collision.initial_state.diatoms[0].xv, collision.initial_state.electronic_state);

// 			DEBUG_PRINT( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PREPARE EDIA from J,V: " <<  edia_jv );
// 			DEBUG_PRINT( "EDIA from Erot Evib " << edia );
			
		
			
// 			//tau = calc_tau(edia, r_in, r_out, collision.initial_state.diatoms[0].xj, collision.initial_state.electronic_state);	
// 			tau = calc_tau(edia_jv, r_in, r_out, collision.initial_state.diatoms[0].xj, collision.initial_state.electronic_state);	


// 				//rand5 = distribution(ic_rng);
// 			DEBUG_PRINT( "RAND 5: " << rand5 ); 
// 				double test_time = rand5*tau*Units::AUTIME_TO_FS;			
// 			if(test_time >= 0.50*tau){
// 				bond_length = r_out;
// 				drift_time = test_time - 0.50*tau;
// 			}else if (test_time < 0.50*tau){
// 				bond_length = r_in;
// 				drift_time = test_time;
// 			}
// 			vib_mom = 0.0;
// 			DEBUG_PRINT( "								DRIFT TIME :" <<  drift_time );
// 			//double xj = jFromErot(collision.initial_state.diatoms[0].erot);
// 			//double xv = vFromEvib(collision.initial_state.diatoms[0].evib);

// 			//bond_length = getBondLength(xj,xv,collision.initial_state.electronic_state);
			
// 			//vib_mom = sqrt(collision.initial_state.diatoms[0].evib * 2. * collision.mu_coll) ; //jvdata[chosen_state].rin + rand*(jvdata[chosen_state].rout - jvdata[chosen_state].rin);
// 			//rot_mom = sqrt(collision.initial_state.diatoms[0].erot * 2. * collision.mu_dia * pow(bond_length,2)) ;//sqrt(jvdata[chosen_state].erot* 2.0 * diatom_mu * pow(bond_length,2));
// 			double jsq;
// 			calcXJ(jsq,collision.initial_state.diatoms[0].xj,-1);
// 			rot_mom = sqrt(jsq);

// 			DEBUG_PRINT( "ROT MOM NUM " << rot_mom );

// 			DEBUG_PRINT( "EDIA, TAU, BL, AND J IN PREPARE: " << edia << " " << tau << " " << bond_length << " " << rot_mom );
// 			DEBUG_PRINT( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PREPARE, PE, Rrel, J " << edia_jv << " " << bond_length << " " << rot_mom);

// 		}
			
// 	if(sample_mode == 0 || sample_mode==1){
// 	/// Randomized orientation sampling
// 	bond_theta = rand1*Units::PI;
// 	bond_phi = rand2*2.0*Units::PI; 

// 	am_theta = rand3*Units::PI;
// 	am_phi = rand4*2.0*Units::PI; 
	
// 	}else if(sample_mode == -1){

// 					// CONTRIVED TESTING CASE: BOND SITS ALONG Y AXIS AND ANGULAR MOMENTUM IS ALONG Z AXIS
				
// 				bond_theta = 0.5*Units::PI;
// 				bond_phi = 0.5*Units::PI; 
		
// 				am_theta = 0.0;
// 				am_phi = 0.0; 
// 	}


// 	if(sample_mode !=-2){
// 	//define unit vector along bond from a1 to a2 
// 	std::vector<double> rel_vec(3, 0.0);
// 	rel_vec[0] = sin(bond_theta)*cos(bond_phi);
// 	rel_vec[1] = sin(bond_theta)*sin(bond_phi);
// 	rel_vec[2] = cos(bond_theta);

// 	std::vector<double> am_vec(3, 0.0);
// 	am_vec[0] = sin(am_theta)*cos(am_phi);
// 	am_vec[1] = sin(am_theta)*sin(am_phi);
// 	am_vec[2] = cos(am_theta);

// 	//define unit vector of velocity for atom 1 that gives desired AM 
// 	std::vector<double> rot_vec(3, 0.0);
// 	rot_vec = crossProduct(rel_vec, am_vec);

// 	double rvMag = sqrt(dotProduct(rot_vec,rot_vec));
// 		//DEBUG_PRINT(  "1: ROTVEC MAG " << rvMag );
// 		for (int i = 0; i<3; i++){
// 		rot_vec[i] = -1.0*rot_vec[i]/rvMag;
// 		};

// 		//DEBUG_PRINT(  "2: ROTVEC MAG " << sqrt(dotProduct(rot_vec,rot_vec)) );
// 	//sample impact parameter 

// 		collision.impactParameter =  0.0   //Units::bMax*sqrt(rand6);
// 		DEBUG_PRINT(  "Impact Parameter: " << collision.impactParameter );



// 	// Give atoms their initial positions and momenta 
// 	for (int i = 0; i<3; i++){
// 		collision.initial_state.atoms[a1ind].pos[i] = bond_length*rel_vec[i]*mass_scale[0];
// 		collision.initial_state.atoms[a2ind].pos[i] = -bond_length*rel_vec[i]*mass_scale[1];
		
// 		DEBUG_PRINT("A1 R " << i+1 << "  " << collision.initial_state.atoms[a1ind].pos[i] );
// 		DEBUG_PRINT("A2 R " << i+1 << "  " <<  collision.initial_state.atoms[a2ind].pos[i] );
	
// 		//collision.initial_state.atoms[a1ind].mom[i] = rot_mom*rot_vec[i]*mass_scale[0]/bond_length;
// 		//collision.initial_state.atoms[a2ind].mom[i] = -rot_mom*rot_vec[i]*mass_scale[1]/bond_length;
// 		collision.initial_state.atoms[a1ind].mom[i] = rot_mom*rot_vec[i]/bond_length;
// 		collision.initial_state.atoms[a2ind].mom[i] = -rot_mom*rot_vec[i]/bond_length;
		
// 		DEBUG_PRINT("A1 RotP " << i+1 << "  " <<  rot_mom*rot_vec[i]*mass_scale[0]/bond_length );
// 		DEBUG_PRINT("A2 RotP " << i+1 << "  " <<  -rot_mom*rot_vec[i]*mass_scale[1]/bond_length );
	
// 		collision.initial_state.atoms[a1ind].mom[i] += vib_mom*rel_vec[i]*mass_scale[0];
// 		collision.initial_state.atoms[a2ind].mom[i] += -vib_mom*rel_vec[i]*mass_scale[1];
		
// 		DEBUG_PRINT("A1 VibP " << i+1 << "  " <<  vib_mom*rel_vec[i]*mass_scale[0] );
// 		DEBUG_PRINT("A2 VibP " << i+1 << "  " <<  -vib_mom*rel_vec[i]*mass_scale[1] );
		
// 		collision.initial_state.atoms[a3ind].pos[i] = 1e10 ;
// 		collision.initial_state.atoms[a3ind].mom[i] = 0.000;

// 	}
	
// 	//DEBUG_PRINT( "J,V " << i_to_jv[chosen_state][1] << " " <<  i_to_jv[chosen_state][2] ); 
// 	DEBUG_PRINT( "Bond Length " << bond_length ); 
	
// 	//DEBUG_PRINT( "Mass 1 " << collision.initial_state.atoms[a1ind].mass ); 
// 	//DEBUG_PRINT( "Mass 2 " << collision.initial_state.atoms[a2ind].mass ); 

// 	//DEBUG_PRINT( "Mass Scale 1 " << mass_scale[0] ); 
// 	//DEBUG_PRINT( "Mass Scale 2 " << mass_scale[1] ); 
	
// 	DEBUG_PRINT( "1: Angular Momentum " << rot_mom ); 
// 	DEBUG_PRINT( "1: Linear Momentum " << rot_mom/bond_length ); 
// 	DEBUG_PRINT( "1: Erot " << collision.initial_state.diatoms[0].erot ); 
// 	DEBUG_PRINT( "1: Evib " << collision.initial_state.diatoms[0].evib ); 

// 	std::vector<double> rrel(3,0.0);
// 	std::vector<double> prel(3,0.0);
// 	for (int i = 0; i<3; i++){
// 	rrel[i]= collision.initial_state.atoms[a1ind].pos[i]-collision.initial_state.atoms[a2ind].pos[i];
// 	prel[i]= collision.initial_state.atoms[a1ind].mom[i]-collision.initial_state.atoms[a2ind].mom[i];
// 	}

// 	double L = sqrt(dotProduct(crossProduct(rrel,prel), crossProduct(rrel,prel)));

// 		}

// 	///////////////// END OF JV ASSIGNMENT LOGIC


// 		DEBUG_PRINT( "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CSURF BEFORE ANALYSIS  " << collision.initial_state.electronic_state  ); 
// 		DEBUG_PRINT( "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~POSITIONS BEFORE ANALYSIS  " ); 
// 	for(int i = 0 ; i<3;i++){
// 		DEBUG_PRINT( "               " << i );
// 		DEBUG_PRINT( collision.initial_state.atoms[i].pos[0] << " " << collision.initial_state.atoms[i].pos[1] << " " << collision.initial_state.atoms[i].pos[2] ); 

// 	}
// 		DEBUG_PRINT( "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MOMENTA BEFORE ANALYSIS  " );
// 	for(int i = 0 ; i<3;i++){
// 		DEBUG_PRINT( "               " << i );
// 		DEBUG_PRINT( collision.initial_state.atoms[i].mom[0] << " " << collision.initial_state.atoms[i].mom[1] << " " << collision.initial_state.atoms[i].mom[2] ); 

// 	}
// 		MDIntegrator integrator;

// 		// Copy initial conditions to working version
// 		collision.final_state.atoms = collision.initial_state.atoms; 
// 		MDIntegrator::setDiatomFromAtoms(collision, 0, true);
// 		collision.final_state.diatoms = collision.initial_state.diatoms; 
// 		collision.final_state.electronic_state = collision.initial_state.electronic_state; 
		
// 		DEBUG_PRINT( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ANALYZING PARTICLES AS PLACED " );
// 		integrator.printCollision(collision);
// 		int iok = integrator.processOutcome(collision);
// 		DEBUG_PRINT( "POST JV ASSIGN OUTCOME CODE: " << iok );


// 		int go_drift;
// 	#ifdef SPARTA_MD_BUILD
// 	go_drift=1;
// 	#else
// 		std::cout << "\n1 to proceed to drift, 0 to skip, -1 to exit...";
// 		std::cin >> go_drift;
// 	#endif

// 		if(go_drift==0){
// 			std::cout << "Skipping drift..." << std::endl;
// 		}else if(go_drift==1){
// 			std::cout << "Drifting for " << drift_time  <<std::endl;
			
// 		integrator.initializePES();
// 		integrator.setInitialSurface(collision.initial_state.electronic_state);



// 		bool drift = integrator.driftDiatom(collision,drift_time);
// 			int iok2 = integrator.processOutcome(collision);
// 			DEBUG_PRINT( "POST JV ASSIGN OUTCOME CODE: " << iok2 );
		
// 		}else if(go_drift=-1){
// 			exit(0);
// 		}


		
// 	std::vector<double> coll_axis(3,0.0);
// 	coll_axis[2]=1.0;

// 	double mag = 0.0;
// 	std::vector<double> coll_relpos(3,0.0);
// 	for (int i=0;i<3;i++){
// 				coll_relpos[i]=collision.initial_state.atoms[a3ind].pos[i];
// 			mag+=pow(collision.initial_state.atoms[a3ind].pos[i],2);
// 	}
// 	mag=sqrt(mag);

// 	for (int i=0;i<3;i++){
// 		coll_relpos[i] /= mag;
// 	}



// 		// Give the two colliding bodies their initial relative momentum
// 	std::vector<double> coll_mass(2,0.0);
// 	coll_mass[0] = collision.initial_state.atoms[a1ind].mass + collision.initial_state.atoms[a2ind].mass;
// 	coll_mass[1] = collision.initial_state.atoms[a3ind].mass;
// 	double mu_coll = (coll_mass[0]*coll_mass[1])/(coll_mass[0]+coll_mass[1]);

// 	if(mu_coll!=collision.mu_coll){
// 		std::cout << "SERIOUS ERROR IN CALCULATING MU COLL " << collision.mu_coll << " " << mu_coll << std::endl;
// 			exit(1);
// 	}

// 	int etr_sample_flag;
// 	double erel_coll=0.0;

// 	if(sample_mode == 1){
// 		etr_sample_flag = 1;   
// 	}else{
// 		etr_sample_flag = 0;
// 	}

// 	if(etr_sample_flag == -1 ){
// 	erel_coll = 0.15;
// 	}else if (etr_sample_flag == 0){
// 	double temp = 20000.0;
// 	double emax_sample = 15.0;
// 	bool sample_success = false;
	
// 	double erel_sample=0.0;
// 	while(!sample_success){  
	
// 	erel_sample = rand7*emax_sample;
// 	double p_sample = exp(-erel_sample/(Units::KB_EV*temp));
// 	//std::cout << " rand7 " << rand7 << std::endl;
	
// 	if(rand8<=p_sample){
// 		erel_coll = erel_sample*Units::EV_TO_HARTREE;
// 		std::cout << " SAMPLE SUCCESS! " << erel_coll << std::endl;
// 		sample_success=true;
// 	}
// 	rand7 = distribution(ic_rng);  
// 	rand8 = distribution(ic_rng);  
// 	}
	
// 	}else{
// 		erel_coll = collision.etrI;
// 	}

// 	std::cout << " INITIAL COLLISIONAL ENERGY " << collision.etrI << std::endl;
	
// 	double prel_coll = sqrt(2*collision.mu_coll*erel_coll);
// 	DEBUG_PRINT( " INITIAL COLLISIONAL MOMENTA " << prel_coll );
// 	//exit(0);

// 		std::vector<double>mass_fac_coll(2,0.0); 
// 			mass_fac_coll[0] = collision.mu_coll/coll_mass[0];
// 			mass_fac_coll[1] = collision.mu_coll/coll_mass[1];

// 	for (int i=0;i<3;i++){
// 	collision.initial_state.atoms[a1ind].mom[i] += prel_coll*coll_axis[i]*0.5;
// 	collision.initial_state.atoms[a2ind].mom[i] += prel_coll*coll_axis[i]*0.5;
// 	collision.initial_state.atoms[a3ind].mom[i] -= prel_coll*coll_axis[i];
	
	
// 	collision.initial_state.atoms[a3ind].pos[i] = collision.initial_separation*coll_axis[i];
// 	}

// 	std::vector<double> relPos = MDIntegrator::getCollRRel(collision,0); 
// 	double relD = sqrt(dotProduct(relPos,relPos));

// 	collision.imd_ca_theta = acos(relPos[2]/relD);
// 	collision.imd_ca_phi = atan2(relPos[1] , relPos[0]); 

// 	//Add impact parameter, +ve y, collaxis Z means collisional angular momentum is along -ve X
// 	collision.initial_state.atoms[a3ind].pos[1] = collision.impactParameter;
// 	std::vector<double> relPos2 = MDIntegrator::getCollRRel(collision,0); 

// 	DEBUG_PRINT( "INIT SEP VEC " <<  relPos2[0] << " " << relPos2[1] << " " << relPos2[2] );
// 	DEBUG_PRINT( "IMPACT PARAMETER " << collision.initial_state.atoms[a3ind].pos[1] );

//     // Copy initial conditions to working version
//     collision.final_state.atoms=collision.initial_state.atoms; 
//     MDIntegrator::setDiatomFromAtoms(collision, 0, true);
//     collision.final_state.diatoms=collision.initial_state.diatoms; 
// }

//  void CollideMD::setupO2O2Collision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2) {
//         std::cout << "Setting up O2 + O2 collision" << std::endl;
        
//         double bond_length = 2.21; // Angstrom
        
//         // First O2 molecule (atoms 0 and 1)
//         collision.initial_state.atoms[0].pos[0] = -5.0 - bond_length/2.0;
//         collision.initial_state.atoms[0].pos[1] = 0.0;
//         collision.initial_state.atoms[0].pos[2] = 0.0;
//         collision.initial_state.atoms[0].mom[0] = 500.0;
//         collision.initial_state.atoms[0].mom[1] = 0.0;
//         collision.initial_state.atoms[0].mom[2] = 0.0;
//         collision.initial_state.atoms[0].mass = 15.999;
//         collision.initial_state.atoms[0].atomic_number = 8;
        
//         collision.initial_state.atoms[1].pos[0] = -5.0 + bond_length/2.0;
//         collision.initial_state.atoms[1].pos[1] = 0.0;
//         collision.initial_state.atoms[1].pos[2] = 0.0;
//         collision.initial_state.atoms[1].mom[0] = 500.0;
//         collision.initial_state.atoms[1].mom[1] = 0.0;
//         collision.initial_state.atoms[1].mom[2] = 0.0;
//         collision.initial_state.atoms[2].pos[0] = 5.0 - bond_length/2.0;
//         collision.initial_state.atoms[2].pos[1] = collision.impactParameter;
//         collision.initial_state.atoms[2].pos[2] = 0.0;
//         collision.initial_state.atoms[2].mom[0] = -500.0;
//         collision.initial_state.atoms[2].mom[1] = 0.0;
//         collision.initial_state.atoms[2].mom[2] = 0.0;
//         collision.initial_state.atoms[2].mass = 15.999;
//         collision.initial_state.atoms[2].atomic_number = 8;
        
//         collision.initial_state.atoms[3].pos[0] = 5.0 + bond_length/2.0;
//         collision.initial_state.atoms[3].pos[1] = collision.impactParameter;
//         collision.initial_state.atoms[3].pos[2] = 0.0;
//         collision.initial_state.atoms[3].mom[0] = -500.0;
//         collision.initial_state.atoms[3].mom[1] = 0.0;
//         collision.initial_state.atoms[3].mom[2] = 0.0;
//         collision.initial_state.atoms[3].mass = 15.999;
//         collision.initial_state.atoms[3].atomic_number = 8;
// }



	////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

int CollideMD::analyzeOutcome(const CollisionData& collision) {
        // Simple analysis based on final atom separations
        double max_separation = 0.0;
        
        for (size_t i = 0; i < collision.final_state.atoms.size(); ++i) {
            for (size_t j = i + 1; j < collision.final_state.atoms.size(); ++j) {
                double dx = collision.final_state.atoms[i].pos[0] - collision.final_state.atoms[j].pos[0];
                double dy = collision.final_state.atoms[i].pos[1] - collision.final_state.atoms[j].pos[1];
                double dz = collision.final_state.atoms[i].pos[2] - collision.final_state.atoms[j].pos[2];
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist > max_separation) {
                    max_separation = dist;
                }
            }
        }
        
        if (max_separation > 10.0) {
            return 2; // Dissociation
        } else if (max_separation > 5.0) {
            return 1; // Inelastic
        } else {
            return 0; // Elastic
        }
}


bool CollideMD::updateColl(CollisionData& collision, Particle::OnePart* sparta_p1, 
                           Particle::OnePart* sparta_p2, Particle::OnePart* sparta_p3) {

    if (!collision.success) {
        std::cerr << "Warning: Updating SPARTA from failed collision" << std::endl;
        return false;
    }

    // 1. Identify which SPARTA pointer is the atom and which is the diatom
    char* p1_name = particle->species[sparta_p1->ispecies].id;
    char* p2_name = particle->species[sparta_p2->ispecies].id;
    bool is_p1_diatom = (strlen(p1_name) >= 2);
    
    Particle::OnePart* sparta_atom = is_p1_diatom ? sparta_p2 : sparta_p1;
    Particle::OnePart* sparta_diatom = is_p1_diatom ? sparta_p1 : sparta_p2;

    int atom_species_idx = sparta_atom->ispecies;

    // 2. Calculate Laboratory Center of Mass (COM) Velocity
    double m_atom_kg = particle->species[sparta_atom->ispecies].molwt;
    double m_diatom_kg = particle->species[sparta_diatom->ispecies].molwt;
    double m_tot = m_atom_kg + m_diatom_kg;
    
    std::vector<double> v_com(3, 0.0);
    for (int i = 0; i < 3; i++) {
        v_com[i] = (m_atom_kg * sparta_atom->v[i] + m_diatom_kg * sparta_diatom->v[i]) / m_tot;
    }

    // 3. Process the collision outcome
    switch(collision.reaction_outcome) {
        
        // =======================================================================
        // CASE 0: ELASTIC COLLISION
        // Internal states remain identical. Only scattering angles/velocities change.
        // =======================================================================
        case 0: 
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~ ELASTIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

            std::vector<double> relMom = MDIntegrator::getCollPRel(collision, 1); 
            double relV = sqrt(dotProduct(relMom, relMom));
            if (relV < 1e-12) relV = 1e-12;

            collision.fmd_ca_theta = acos(relMom[2] / relV);
            collision.fmd_ca_phi = atan2(relMom[1], relMom[0]); 

            std::vector<double> new_sparta_uvrel = transformVelocityAfterMD(
                collision.i_sparta_coll_axis, 1.0, collision.fmd_ca_theta, collision.fmd_ca_phi);

            double suvr_mag = sqrt(dotProduct(new_sparta_uvrel, new_sparta_uvrel));
            for (int i = 0; i < 3; i++) new_sparta_uvrel[i] /= suvr_mag;

            double v_rel_final_au = sqrt(2.0 * collision.final_state.etr / collision.mu_coll);
            double v_rel_final_mps = v_rel_final_au / Units::MPS_TO_AUVEL;

            double v_atom_mag = v_rel_final_mps * (m_diatom_kg / m_tot);
            double v_diatom_mag = v_rel_final_mps * (m_atom_kg / m_tot);

            for (int i = 0; i < 3; i++) {
                sparta_atom->v[i]   = v_com[i] + (v_atom_mag * new_sparta_uvrel[i]);
                sparta_diatom->v[i] = v_com[i] - (v_diatom_mag * new_sparta_uvrel[i]);
            }
            
            // Note: erot, evib, eelec are left completely untouched!
            break;
        }

        // =======================================================================
        // CASE 1: INELASTIC COLLISION
        // Energies transfer between translation, rotation, and vibration.
        // =======================================================================
        case 1: 
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~ INELASTIC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

            // 1. Update Velocities (Same math as elastic)
            std::vector<double> relMom = MDIntegrator::getCollPRel(collision, 1); 
            double relV = sqrt(dotProduct(relMom, relMom));
            if (relV < 1e-12) relV = 1e-12;

            collision.fmd_ca_theta = acos(relMom[2] / relV);
            collision.fmd_ca_phi = atan2(relMom[1], relMom[0]); 

            std::vector<double> new_sparta_uvrel = transformVelocityAfterMD(
                collision.i_sparta_coll_axis, 1.0, collision.fmd_ca_theta, collision.fmd_ca_phi);

            double suvr_mag = sqrt(dotProduct(new_sparta_uvrel, new_sparta_uvrel));
            for (int i = 0; i < 3; i++) new_sparta_uvrel[i] /= suvr_mag;

            double v_rel_final_au = sqrt(2.0 * collision.final_state.etr / collision.mu_coll);
            double v_rel_final_mps = v_rel_final_au / Units::MPS_TO_AUVEL;

            double v_atom_mag = v_rel_final_mps * (m_diatom_kg / m_tot);
            double v_diatom_mag = v_rel_final_mps * (m_atom_kg / m_tot);

            for (int i = 0; i < 3; i++) {
                sparta_atom->v[i]   = v_com[i] + (v_atom_mag * new_sparta_uvrel[i]);
                sparta_diatom->v[i] = v_com[i] - (v_diatom_mag * new_sparta_uvrel[i]);
            }

            // 2. Update Internal Energies
            DiAtomData diatom = collision.final_state.diatoms[0];
            sparta_diatom->erot = diatom.erot * Units::HARTREE_TO_J;
            sparta_diatom->xj = diatom.xj;
            sparta_diatom->evib = diatom.evib * Units::HARTREE_TO_J; 
            sparta_diatom->xv = diatom.xv;
            
            double eelec_atom, eelec_diatom;
            MDIntegrator::getEelec(collision.final_state.electronic_state, eelec_atom, eelec_diatom);
            sparta_atom->eelec = eelec_atom * Units::HARTREE_TO_J;
            sparta_diatom->eelec = eelec_diatom * Units::HARTREE_TO_J;
            
            break;
        }

        // =======================================================================
        // CASE 2: EXCHANGE COLLISION
        // Atoms swap. e.g., N + O2 -> NO + O.
        // =======================================================================
        case 2: 
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~ EXCHANGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

            // Optional future feature: If the atoms exchanged were different species, 
            // the mass and ispecies IDs must be updated here!
            // e.g., sparta_atom->ispecies = new_atom_id;
            //       m_atom_kg = new_atom_mass;
            
            // 1. Update Velocities
            std::vector<double> relMom = MDIntegrator::getCollPRel(collision, 1); 
            double relV = sqrt(dotProduct(relMom, relMom));
            if (relV < 1e-12) relV = 1e-12;

            collision.fmd_ca_theta = acos(relMom[2] / relV);
            collision.fmd_ca_phi = atan2(relMom[1], relMom[0]); 

            std::vector<double> new_sparta_uvrel = transformVelocityAfterMD(
                collision.i_sparta_coll_axis, 1.0, collision.fmd_ca_theta, collision.fmd_ca_phi);

            double suvr_mag = sqrt(dotProduct(new_sparta_uvrel, new_sparta_uvrel));
            for (int i = 0; i < 3; i++) new_sparta_uvrel[i] /= suvr_mag;

            // Use the UPDATED mu_coll in case the species masses changed during exchange
            double v_rel_final_au = sqrt(2.0 * collision.final_state.etr / collision.mu_coll);
            double v_rel_final_mps = v_rel_final_au / Units::MPS_TO_AUVEL;

            double v_atom_mag = v_rel_final_mps * (m_diatom_kg / m_tot);
            double v_diatom_mag = v_rel_final_mps * (m_atom_kg / m_tot);

            for (int i = 0; i < 3; i++) {
                sparta_atom->v[i]   = v_com[i] + (v_atom_mag * new_sparta_uvrel[i]);
                sparta_diatom->v[i] = v_com[i] - (v_diatom_mag * new_sparta_uvrel[i]);
            }

            // 2. Update Internal Energies
            DiAtomData diatom = collision.final_state.diatoms[0];
            sparta_diatom->erot = diatom.erot * Units::HARTREE_TO_J;
            sparta_diatom->xj = diatom.xj;
            sparta_diatom->evib = diatom.evib * Units::HARTREE_TO_J; 
            sparta_diatom->xv = diatom.xv;
            

            // Explicitly initialize to 0.0 to prevent garbage memory leaks!
            double eelec_atom = 0.0;
            double eelec_diatom = 0.0;

            // Also, double-check that the MD library actually populated final_state.electronic_state!
            // If it didn't, fall back to the initial state to be safe.
            int final_es = collision.final_state.electronic_state;
            if (final_es < 0 || final_es > 6) final_es = collision.initial_state.electronic_state;

            MDIntegrator::getEelec(final_es, eelec_atom, eelec_diatom);

            // double eelec_atom, eelec_diatom;
            // MDIntegrator::getEelec(collision.final_state.electronic_state, eelec_atom, eelec_diatom);
            sparta_atom->eelec = eelec_atom * Units::HARTREE_TO_J;
            sparta_diatom->eelec = eelec_diatom * Units::HARTREE_TO_J;
            
            break;
        }

        // =======================================================================
        // CASE 3: DISSOCIATION
        // Diatom breaks apart into two free atoms.
        // =======================================================================
        case 3: 
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~ DISSOCIATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

            // The diatom is destroyed. Internal energies become 0 for all atoms.
            sparta_atom->erot = 0.0; sparta_atom->xj = 0.0;
            sparta_atom->evib = 0.0; sparta_atom->xv = 0.0;
            sparta_atom->eelec = 0.0; 

            AtomData atom1 = collision.final_state.atoms[0];
            AtomData atom2 = collision.final_state.atoms[1];
            AtomData atom3 = collision.final_state.atoms[2];
            
            // Transform the 3 MD momenta back into the SPARTA laboratory orientation
            auto velocities_A = transformThreeParticleVelocities(
                collision.i_sparta_coll_axis, atom1.mom, atom2.mom, atom3.mom);

            std::vector<double> v1(3, 0.0), v2(3, 0.0), v3(3, 0.0);
            for (int i = 0; i < 3; i++) {
                // V = (P / mass) converted to m/s
                v1[i] = (velocities_A[0][i] / atom1.mass) / Units::MPS_TO_AUVEL;
                v2[i] = (velocities_A[1][i] / atom2.mass) / Units::MPS_TO_AUVEL;
                v3[i] = (velocities_A[2][i] / atom3.mass) / Units::MPS_TO_AUVEL;
            }

            // Apply the COM velocity and update the original two particles
            sparta_diatom->ispecies = atom_species_idx; // Convert diatom to an atom species!
            
            for (int i = 0; i < 3; i++) {
                sparta_atom->v[i]   = v_com[i] + v1[i];
                sparta_diatom->v[i] = v_com[i] + v2[i];
            }
            
            sparta_diatom->erot = 0.0; sparta_diatom->xj = 0.0;
            sparta_diatom->evib = 0.0; sparta_diatom->xv = 0.0;
            sparta_diatom->eelec = 0.0;

            // Create the third atom if it hasn't been instantiated yet
            if (sparta_p3 == nullptr) {
                int id = MAXSMALLINT * random->uniform();
                Particle::OnePart *particles = particle->particles;
                double x[3];
                double v[3]; 
                
                // Spawn it at the exact location of the original atom
                memcpy(x, sparta_atom->x, 3 * sizeof(double));
                memcpy(v, sparta_atom->v, 3 * sizeof(double)); // Dummy initial V, overwritten below
                
                int ielectron_flag = (ambiflag && sparta_atom->ispecies == ambispecies);
                int jelectron_flag = (ambiflag && sparta_diatom->ispecies == ambispecies);
                
                int reallocflag = particle->add_particle(id, atom_species_idx, sparta_atom->icell, x, v, 0.0, 0.0, 0.0, 0.0, 0.0);
                
                if (reallocflag) {
                    if (!ielectron_flag) sparta_atom = particle->particles + (sparta_atom - particles);
                    if (!jelectron_flag) sparta_diatom = particle->particles + (sparta_diatom - particles);
                }
                sparta_p3 = &particle->particles[particle->nlocal - 1];
            }

            // Assign the third atom its COM-corrected velocity
            sparta_p3->ispecies = atom_species_idx; 
            sparta_p3->icell = sparta_diatom->icell;
            for (int i = 0; i < 3; i++) {
                sparta_p3->x[i] = 0.0; // Handled by SPARTA position advection later
                sparta_p3->v[i] = v_com[i] + v3[i];
            }
            
            sparta_p3->erot = 0.0; sparta_p3->xj = 0.0;
            sparta_p3->evib = 0.0; sparta_p3->xv = 0.0;
            sparta_p3->eelec = 0.0;
            sparta_p3->flag = 0;
            sparta_p3->dtremain = sparta_diatom->dtremain;
            sparta_p3->weight = sparta_diatom->weight;
            
            break;
        }
    }

    return true;
}


// bool CollideMD::updateColl(CollisionData& collision, Particle::OnePart* sparta_p1, 
//                            Particle::OnePart* sparta_p2,
//                            Particle::OnePart* sparta_p3) {

//     if (!collision.success) {
//         std::cerr << "Warning: Updating SPARTA from failed collision" << std::endl;
//         return false;
//     }

//     // 1. Identify which SPARTA pointer is the atom and which is the diatom dynamically
//     char* p1_name = particle->species[sparta_p1->ispecies].id;
//     char* p2_name = particle->species[sparta_p2->ispecies].id;
//     bool is_p1_diatom = (strlen(p1_name) >= 2);
//     bool is_p2_diatom = (strlen(p2_name) >= 2);

//     Particle::OnePart* sparta_atom = nullptr;
//     Particle::OnePart* sparta_diatom = nullptr;
//     Particle::OnePart* sparta_atom2 = nullptr;
//     Particle::OnePart* sparta_atom3 = nullptr;

//     if (is_p1_diatom && !is_p2_diatom) {
//         sparta_diatom = sparta_p1;
//         sparta_atom = sparta_p2;
//     } else if (!is_p1_diatom && is_p2_diatom) {
//         sparta_atom = sparta_p1;
//         sparta_diatom = sparta_p2;
//     } else {
//         std::cerr << "Warning: updateColl currently only supports atom-diatom." << std::endl;
//         return false;
//     }

//     // 2. Kinematics setup
//     std::vector<double> relMom = MDIntegrator::getCollPRel(collision, 1); 
//     double relV = sqrt(dotProduct(relMom, relMom));

//     collision.fmd_ca_theta = acos(relMom[2] / relV);
//     collision.fmd_ca_phi = atan2(relMom[1], relMom[0]); 

//     DEBUG_PRINT("FINAL THETA/PHI :      " << collision.fmd_ca_theta << " " << collision.fmd_ca_phi);

//     std::vector<double> new_sparta_uvrel = transformVelocityAfterMD(collision.i_sparta_coll_axis, 1.0, 
//                                                                     collision.fmd_ca_theta, collision.fmd_ca_phi);

//     double suvr_mag = sqrt(dotProduct(new_sparta_uvrel, new_sparta_uvrel));
//     for (int i = 0; i < 3; i++) {
//         new_sparta_uvrel[i] /= suvr_mag;
//     }

//     DEBUG_PRINT("NEW SUVR " << new_sparta_uvrel[0] << " " << new_sparta_uvrel[1] << " " << new_sparta_uvrel[2]);
//     DEBUG_PRINT("COLLISION OUTCOME : " << collision.reaction_outcome);

//     // 3. Print Energy Changes
//     std::cout << "ENERGY CHANGE: EROTF-EROTI " << collision.final_state.erot - collision.initial_state.erot << std::endl;	
//     std::cout << "ENERGY CHANGE: EVIBF-EVIBI " << collision.final_state.evib - collision.initial_state.evib << std::endl;	
//     if (collision.final_state.diatoms.size() > 0 && collision.initial_state.diatoms.size() > 0) {
//         std::cout << "ENERGY CHANGE: XJF , XJI " << collision.final_state.diatoms[0].xj << " , " << collision.initial_state.diatoms[0].xj << std::endl;	
//         std::cout << "ENERGY CHANGE: XVF, XVI " << collision.final_state.diatoms[0].xv << " , " << collision.initial_state.diatoms[0].xv << std::endl;	
//     }
//     std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

//     // 4. Update SPARTA based on reaction outcome
//     std::vector<double> atomMom(3, 0.0);
//     std::vector<double> diAtomMom(3, 0.0);
//     double atomMomMag, diAtomMomMag;

//     switch(collision.reaction_outcome) {
//         case 0:  
//             std::cout << "~~~~~~~~~~~~~~~~~~~~~ ELASTIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
//             break;
//         case 1:
//         case 2: {
//             if (collision.reaction_outcome == 1) std::cout << "~~~~~~~~~~~~~~~~~~~~~ INELASTIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
//             if (collision.reaction_outcome == 2) std::cout << "~~~~~~~~~~~~~~~~~~~~~ EXCHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            
//             int a1 = collision.final_state.diatoms[0].atom1_index;
//             int a2 = collision.final_state.diatoms[0].atom2_index;
//             int a3 = MDIntegrator::getFreeAtomIndex(collision.final_state.diatoms[0]);

//             AtomData atom = collision.final_state.atoms[a3];
//             DiAtomData diatom = collision.final_state.diatoms[0];

//             atomMomMag = sqrt(dotProduct(atom.mom, atom.mom));
//             diAtomMomMag = sqrt(dotProduct(diatom.mom, diatom.mom));

//             for (int i = 0; i < 3; i++){
//                 atomMom[i] = atomMomMag * new_sparta_uvrel[i];
//                 diAtomMom[i] = -diAtomMomMag * new_sparta_uvrel[i];

//                 sparta_atom->v[i] = (atomMom[i] / atom.mass) / Units::MPS_TO_AUVEL;	
//                 sparta_diatom->v[i] = (diAtomMom[i] / diatom.mass) / Units::MPS_TO_AUVEL;	
//             }
        
//             sparta_atom->erot = 0.0;
//             sparta_atom->xj = 0.0;
//             sparta_atom->evib = 0.0;
//             sparta_atom->xv = 0.0;
            
//             sparta_diatom->erot = diatom.erot * Units::HARTREE_TO_J; 
//             sparta_diatom->xj = diatom.xj; 
            
//             // Note: If you don't have a global ZPE array anymore, just remove it or use a getZPE() function
//             sparta_diatom->evib = diatom.evib * Units::HARTREE_TO_J; 
//             sparta_diatom->xv = diatom.xv;
            
//             // Map electronic state back to particle energies
//             double eelec_atom, eelec_diatom;
//             MDIntegrator::getEelec(collision.final_state.electronic_state, eelec_atom, eelec_diatom);
//             sparta_atom->eelec = eelec_atom * Units::HARTREE_TO_J;
//             sparta_diatom->eelec = eelec_diatom * Units::HARTREE_TO_J;
            
//             break;	
//         }
//         case 3: {
//             std::cout << "~~~~~~~~~~~~~~~~~~~~~ DISSOCIATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

//             sparta_atom->erot = 0.0;
//             sparta_atom->xj = 0.0;
//             sparta_atom->evib = 0.0;
//             sparta_atom->xv = 0.0;
//             sparta_atom->eelec = 0.0; // Update with proper electronic energy if needed

//             sparta_atom2 = sparta_p2;
//             sparta_atom3 = sparta_p3;

//             AtomData atom1 = collision.final_state.atoms[0];
//             AtomData atom2 = collision.final_state.atoms[1];
//             AtomData atom3 = collision.final_state.atoms[2];
            
//             auto velocities_A = transformThreeParticleVelocities(collision.i_sparta_coll_axis, atom1.mom, atom2.mom, atom3.mom);

//             std::vector<double> v1(3, 0.0), v2(3, 0.0), v3(3, 0.0);
//             for (int i = 0; i < 3; i++){
//                 v1[i] = velocities_A[0][i] / atom1.mass;
//                 v2[i] = velocities_A[1][i] / atom2.mass;
//                 v3[i] = velocities_A[2][i] / atom3.mass;
//             }

//             // Update original diatom to become an atom
//             sparta_diatom->ispecies = 0; /* species index for atom */
//             for (int i = 0; i < 3; i++) {
//                 sparta_atom->v[i] = v1[i] / Units::MPS_TO_AUVEL;
//                 sparta_diatom->v[i] = v2[i] / Units::MPS_TO_AUVEL;
//             }
//             sparta_diatom->erot = 0.0;
//             sparta_diatom->xj = 0.0;
//             sparta_diatom->evib = 0.0;
//             sparta_diatom->xv = 0.0;
//             sparta_diatom->eelec = 0.0;

//             if (sparta_p3 == nullptr) {
//                 int id = MAXSMALLINT * random->uniform();
//                 Particle::OnePart *particles = particle->particles;
//                 double x[3];
//                 double v[3]; 
//                 memcpy(x, sparta_atom->x, 3 * sizeof(double));
//                 memcpy(v, sparta_atom->v, 3 * sizeof(double));
                
//                 int ielectron_flag = (ambiflag && sparta_atom->ispecies == ambispecies);
//                 int jelectron_flag = (ambiflag && sparta_diatom->ispecies == ambispecies);
                
//                 int reallocflag = particle->add_particle(id, 0, sparta_atom->icell, x, v, 0.0, 0.0, 0.0, 0.0, 0.0);
//                 if (reallocflag) {
//                     if (!ielectron_flag) sparta_atom = particle->particles + (sparta_atom - particles);
//                     if (!jelectron_flag) sparta_diatom = particle->particles + (sparta_diatom - particles);
//                 }
//                 sparta_p3 = &particle->particles[particle->nlocal - 1];

//             }

//             sparta_p3->ispecies = 0; /* species index for second atom */
//             sparta_p3->icell = sparta_diatom->icell;
//             for (int i = 0; i < 3; i++) {
//                 sparta_p3->x[i] = 0.0;
//                 sparta_p3->v[i] = v3[i] / Units::MPS_TO_AUVEL;
//             }
//             sparta_p3->erot = 0.0;
//             sparta_p3->xj = 0.0;
//             sparta_p3->evib = 0.0;
//             sparta_p3->xv = 0.0;
//             sparta_p3->eelec = 0.0;
//             sparta_p3->flag = 0;
//             sparta_p3->dtremain = sparta_diatom->dtremain;
//             sparta_p3->weight = sparta_diatom->weight;
            
//             break;
//         }
//     }

//     return true;
// }


// bool CollideMD::updateColl(CollisionData& collision, Particle::OnePart* sparta_p1, 
//                        Particle::OnePart* sparta_p2,
//                        Particle::OnePart* sparta_p3) {



// 	Particle::OnePart* sparta_atom;
// 	Particle::OnePart* sparta_atom2;
// 	Particle::OnePart* sparta_atom3;
//         Particle::OnePart* sparta_diatom;


// 	if(collision.p1Type==0 && collision.p2Type==1){
// 		sparta_atom = sparta_p1;
// 		sparta_diatom = sparta_p2;
// 	}else if(collision.p1Type==1 && collision.p2Type==0){
// 		sparta_atom = sparta_p2;
// 		sparta_diatom = sparta_p1;

// 	}

// 	if (!collision.success) {
//         std::cerr << "Warning: Updating SPARTA from failed collision" << std::endl;
//         return false;
//     }

// 	std::vector<double> relMom = MDIntegrator::getCollPRel(collision,1); 
// 	double relV = sqrt(dotProduct(relMom,relMom));

// 	collision.fmd_ca_theta = acos(relMom[2]/relV);
// 	collision.fmd_ca_phi = atan2(relMom[1] , relMom[0]); 

// 	DEBUG_PRINT( "FINAL THETA/PHI :      " << collision.fmd_ca_theta << " " << collision.fmd_ca_phi );


// 	std::vector<double> new_sparta_uvrel(3,0.0);
// 	//new_sparta_uvrel = rotateCollAxis(collision.i_sparta_coll_axis, collision.fmd_ca_theta , collision.fmd_ca_phi);


// 	DEBUG_PRINT( " SIZE OF COLL I SPARTA AXIS " << collision.i_sparta_coll_axis.size() );
// 	DEBUG_PRINT( " INITIAL SPARTA COLL AXIS : " << collision.i_sparta_coll_axis[0] << " " << collision.i_sparta_coll_axis[1] << " " << collision.i_sparta_coll_axis[2] );

// 		new_sparta_uvrel = transformVelocityAfterMD(collision.i_sparta_coll_axis,
// 																1.0,
// 																collision.fmd_ca_theta,
// 																collision.fmd_ca_phi);
		

// 		DEBUG_PRINT( "Initial relative velocity: (" << collision.i_sparta_coll_axis[0] << ", " 
// 				<< collision.i_sparta_coll_axis[1] << ", " << collision.i_sparta_coll_axis[2] );





// 		DEBUG_PRINT( collision.fmd_ca_theta << " " << collision.fmd_ca_phi );

// 		double suvr_mag = sqrt(dotProduct(new_sparta_uvrel,new_sparta_uvrel));

// 		for (int i=0;i<3;i++){
// 		new_sparta_uvrel[i] /= suvr_mag;
// 	}


// 		int a1,a2,a3;
// 		int a1 = collision.final_state.diatoms[0].atom1_index;
// 		int a2 = collision.final_state.diatoms[0].atom2_index;
// 		int a3 = MDIntegrator::getFreeAtomIndex(collision.final_state.diatoms[0]);
// 		//MDIntegrator::getAtomInds(a1,a2,a3, collision.final_arrangement);


		
// 		DEBUG_PRINT( " NEW SUVR " << new_sparta_uvrel[0] << " " << new_sparta_uvrel[1] << " " << new_sparta_uvrel[2] );
// 		DEBUG_PRINT( "MAG : " << sqrt(dotProduct(new_sparta_uvrel,new_sparta_uvrel)) );
// 		DEBUG_PRINT( "COLLISION OUTCOME : " << collision.reaction_outcome );


// 	//CASES: 
// 	//0 ELASTIC: No Update?
// 	//1 INELASTIC: In at least one DoF, update all of them
// 	//2 REACTIVE:Exchange: update all
// 	//3 Reactive: Dissociation, update all and make a new sparta particle with new values
// 	//Need to register react tally for sparta somehow


// 			std::vector<double> atomMom(3,0.0);
// 			std::vector<double> diAtomMom(3,0.0);
// 			double atomMomMag;
// 			double diAtomMomMag;

// 	CollisionData::AtomData atom;
// 	CollisionData::AtomData atom2;
// 	CollisionData::AtomData atom3;

// 	CollisionData::DiAtomData diatom;
// 	switch(collision.reaction_outcome){
// 			case 0:  
// 			std::cout << "~~~~~~~~~~~~~~~~~~~~~ ELASTIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
// 				break;
// 		case 1:
// 			std::cout << "~~~~~~~~~~~~~~~~~~~~~ INELASTIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
// 				break;
// 		case 2:
// 			std::cout << "~~~~~~~~~~~~~~~~~~~~~ EXCHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
// 				break;
// 		case 3:
// 			std::cout << "~~~~~~~~~~~~~~~~~~~~~ DISSOCIATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

// 			if (sparta_p3 != nullptr) {
// 				#ifdef SPARTA_MD_BUILD
// 				// SPARTA version - DON'T delete, SPARTA owns this memory
// 				sparta_p3 = nullptr;  // Just set to nullptr
// 				#else
// 				// Standalone version - you created it with 'new', so you delete it
// 					delete sparta_p3;
// 					sparta_p3 = nullptr;
// 				#endif
// 			}
		
// 			break;
// 	}

// 	std::cout << "ENERGY CHANGE: EROTF-EROTI " << collision.erotF-collision.erotI << std::endl;	
// 	std::cout << "ENERGY CHANGE: EVIBF-EVIBI " << collision.evibF-collision.evibI << std::endl;	
// 	std::cout << "ENERGY CHANGE: XJF , XJI " << collision.final_state.diatoms[0].xj << " , " << collision.initial_state.diatoms[0].xj << std::endl;	
// 	std::cout << "ENERGY CHANGE: XVF, XVI " << collision.final_state.diatoms[0].xv << " , " << collision.initial_state.diatoms[0].xv << std::endl;	
// 	std::cout << "ENERGY CHANGE: EELECF-EELECI " << collision.eelecF-collision.eelecI << std::endl;	
// 	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

// 	switch(collision.reaction_outcome){
// 		case 0: 
// 			break;
// 		case 1:
// 		case 2:
		
		
// 			atom = collision.final_state.atoms[a3];
// 			diatom = collision.final_state.diatoms[0];

// 			atomMomMag  = sqrt(dotProduct(atom.mom,atom.mom));
// 			diAtomMomMag  = sqrt(dotProduct(diatom.mom,diatom.mom));


// 			//Update atom particle as atom, diatom particle as diatom
// 			for (int i=0; i<3; i++){

// 				atomMom[i] = atomMomMag*new_sparta_uvrel[i];
// 				diAtomMom[i] = -diAtomMomMag*new_sparta_uvrel[i];

// 				sparta_atom->v[i] = (atomMom[i]/atom.mass)/Units::MPS_TO_AUVEL ;	
			
// 				sparta_diatom->v[i] = (diAtomMom[i]/diatom.mass)/Units::MPS_TO_AUVEL;	
// 				//DEBUG_PRINT( "SPARTA V " << sparta_p2->v[i]  );
// 				//DEBUG_PRINT( "MOM/MASS " << diatom.mom[i] << " " << diatom.mass );
// 			}
		
// 			sparta_atom->erot = 0.0*Units::HARTREE_TO_J;
// 			sparta_atom->xj = 0.0;
// 			sparta_atom->evib = 0.0*Units::HARTREE_TO_J;
// 			sparta_atom->xv = 0.0;
// 			sparta_atom->eelec = 0.0*Units::HARTREE_TO_J;

// 			sparta_diatom->erot = diatom.erot*Units::HARTREE_TO_J; 
// 			sparta_diatom->xj = diatom.xj; 
// 			sparta_diatom->evib = (diatom.evib-ZPE[collision.final_state.electronic_state])*Units::HARTREE_TO_J;
// 			//sparta_diatom->evib = (diatom.evib)*Units::HARTREE_TO_J;
// 			sparta_diatom->xv = diatom.xv;
// 			sparta_diatom->eelec = diatom.eelec*Units::HARTREE_TO_J;
			
// 			break;	

// 		case(3):

// 			sparta_atom->erot = 0.0*Units::HARTREE_TO_J;
// 			sparta_atom->xj = 0.0;
// 			sparta_atom->evib = 0.0*Units::HARTREE_TO_J;
// 			sparta_atom->xv = 0.0;
// 			sparta_atom->eelec = 0.0*Units::HARTREE_TO_J;


// 			sparta_atom2 = sparta_p2;
// 			sparta_atom3 = sparta_p3;

// 			atom = collision.final_state.atoms[a1];
// 			atom2 = collision.final_state.atoms[a2];
// 			atom3 = collision.final_state.atoms[a3];
// 			auto velocities_A = transformThreeParticleVelocities(collision.i_sparta_coll_axis, atom.mom, atom2.mom, atom3.mom);

// 			std::vector<double> a1_vf = velocities_A[0];
// 			std::vector<double> a2_vf = velocities_A[1]; 
// 			std::vector<double> a3_vf = velocities_A[2];
			
// 			std::vector<double> v1(3,0.0);
// 			std::vector<double> v2(3,0.0);
// 			std::vector<double> v3(3,0.0);
			
			
// 			for (int i=0; i<3; i++){
// 				v1[i] = a1_vf[i]/atom.mass;
// 				v2[i] = a2_vf[i]/atom2.mass;
// 				v3[i] = a3_vf[i]/atom3.mass;
// 			}


// 				// Update original diatom to become an atom
// 				sparta_diatom->ispecies = 0/* species index for atom */;
// 				for (int i = 0; i < 3; i++) {
// 					sparta_atom->v[i] = v1[i]/Units::MPS_TO_AUVEL;
// 					sparta_diatom->v[i] = v2[i]/Units::MPS_TO_AUVEL;
// 				}
// 				sparta_diatom->erot = 0.0*Units::HARTREE_TO_J;
// 				sparta_diatom->xj = 0.0;
// 				sparta_diatom->evib = 0.0*Units::HARTREE_TO_J;
// 				sparta_diatom->xv = 0.0;
// 				sparta_diatom->eelec = 0.0*Units::HARTREE_TO_J;


// 			DEBUG_PRINT( "3RD PARTICLE POINTER: " << sparta_p3 );
				
			
			
// 			if (sparta_p3 == nullptr) {
// 	#ifdef SPARTA_MD_BUILD
				
// 			DEBUG_PRINT( "DEBUG: Before particle creation - grid=" << grid );
				
// 			if (grid) {
// 					DEBUG_PRINT( "DEBUG: grid->cells=" << grid->cells );
// 				}else{
// 				DEBUG_PRINT( " NO GRID OBJECT " ) ; 
// 			}	
// 			// add 3rd K particle if reaction created it
// 			// index of new K particle = nlocal-1
// 			// if add_particle() performs a realloc:
// 			//   make copy of x,v, then repoint ip,jp to new particles data struct
// 			//   unless electron

// 			int id = MAXSMALLINT*random->uniform();
		
			
// 			Particle::OnePart *particles = particle->particles;
// 			double x[3];
// 			double v[3]; 
// 				memcpy(x,sparta_atom->x,3*sizeof(double));
// 					memcpy(v,sparta_atom->v,3*sizeof(double));
				
			
			
// 			int ielectron_flag = (ambiflag && sparta_atom->ispecies == ambispecies);
// 			int jelectron_flag = (ambiflag && sparta_diatom->ispecies == ambispecies);
			
// 			int reallocflag =
// 				particle->add_particle(id,0,sparta_atom->icell,x,v,0.0,0.0,0.0,0.0,0.0);
// 			if (reallocflag) {
// 				if (!ielectron_flag)
// 				sparta_atom = particle->particles + (sparta_atom - particles);
// 				if (!jelectron_flag)
// 				sparta_diatom = particle->particles + (sparta_diatom - particles);
// 			}

// 			sparta_p3 = &particle->particles[particle->nlocal-1];
				
// 			DEBUG_PRINT( "DEBUG: After particle creation - grid=" << grid );
// 					if (grid) {
// 						DEBUG_PRINT( "DEBUG: grid->cells=" << grid->cells );
// 					}else{
// 					DEBUG_PRINT( " NO GRID OBJECT " ) ; 
// 				}
							
				
// 	#else
// 				// Standalone version - create a local particle or handle differently
// 				// Option 1: Create a local OnePart object
// 				//Particle::OnePart standalone_particle;
// 				//sparta_p3 = &standalone_particle;
				
// 				// Generate a simple ID for standalone
// 				int particle_id_counter = 2;
// 				//parta_p3->id = particle_id_counter++;
				
// 				// Option 2: Or allocate dynamically if you need persistence
// 				sparta_p3 = new Particle::OnePart();
// 				sparta_p3->id = particle_id_counter++;
				
// 	#endif

// 			}


// 				// Set up the new particle (sparta_p3)
// 				sparta_p3->ispecies = 0   /* species index for second atom */;
// 				sparta_p3->icell = sparta_diatom->icell;
// 				for (int i = 0; i < 3; i++) {
// 					sparta_p3->x[i] = 0.0;
// 					sparta_p3->v[i] = v3[i]/Units::MPS_TO_AUVEL;
// 				}
// 				sparta_p3->erot = 0.0*Units::HARTREE_TO_J;
// 				sparta_p3->xj = 0.0;
// 				sparta_p3->evib = 0.0*Units::HARTREE_TO_J;
// 				sparta_p3->xv = 0.0;
// 				sparta_p3->eelec = 0.0*Units::HARTREE_TO_J;
// 				sparta_p3->flag = 0;
// 				sparta_p3->dtremain = sparta_diatom->dtremain;
// 				sparta_p3->weight = sparta_diatom->weight;
			
// 			break;


// 	}

// 	return true;
// }



// Function to run a single collision fully
bool CollideMD::runCollision( 
	               int                 integration_type,
		       Particle::OnePart*  sparta_p1,
		       Particle::OnePart*  sparta_p2,
		       Particle::OnePart*  sparta_p3
			) {
    
    std::cout << "\nRunning collision " << std::endl;
    CollisionData collision;
    //MDIntegrator integrator;

    bool success = true;
    
    //Prepare MD particles from SPARTA ones 
    success &= prepare(config_ ,collision, sparta_p1, sparta_p2, integration_type);
    
	if (success) {
        integrator_.setInitialSurface(collision.initial_state.electronic_state);
    }

    // Collide MD particles
    if (success) success &= QCT::collide(config_,integrator_, collision);

    // Process results of MD particles
    if (success) success &= QCT::process(integrator_, collision);

    if (success) success &= QCT::finalize(collision);
   
    // Update SPARTA particles from MD results
    if (success) success &= this->updateColl(collision, sparta_p1, sparta_p2, sparta_p3);


    std::cout << "\n=== Results =======" << std::endl;
    collision.printSummary();
    
    if (!success) {
        std::cout << "MD Collision FAILED: "  << std::endl;
    } else {
        std::cout << "MD Collision FINISHED: "  << std::endl;
    }

    //exit(0);
    return success;
}



void CollideMD::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp) {
	//bool success = CollideMD::prepare(CollisionData& collision, ip, jp, 0); 
}

int CollideMD::perform_collision(Particle::OnePart *&ip, Particle::OnePart *&jp, Particle::OnePart *&kp) {

	if ((ip->ispecies == 1 && jp->ispecies == 1) || (ip->ispecies == 0 && jp->ispecies == 0)) {

		//std::cout << "Skipping MD collision for pair of species: " << ip->ispecies << " " << jp->ispecies << std::endl; 
		return true; 

	}else {      	
  		std::cout << "Setting up MD collision" << std::endl; 


#ifdef NONADIA_MD_BUILD
#pragma message "COMPILING NONADIABATIC MD"
		
		return runCollision(10,ip, jp, kp);
#else
		
		return runCollision(0,ip, jp, kp);
#endif
	}

}


#ifdef SPARTA_MD_BUILD
#pragma message "COMPILING SPARTA VERSION - SPARTA_MD_BUILD is defined"

double CollideMD::vremax_init(int igroup, int jgroup) {

 // parent has set mixture ptr

  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double cxs = params[isp][jsp].diam*params[isp][jsp].diam*MY_PI;
      //prefactor[isp][jsp] = cxs * pow(2.0*update->boltz*params[isp][jsp].tref/
        //params[isp][jsp].mr,params[isp][jsp].omega-0.5) /
        //tgamma(2.5-params[isp][jsp].omega);
      double beta = MAX(vscale[isp],vscale[jsp]);
      double vrm = 2.0 * cxs * beta;
      vrmgroup = MAX(vrmgroup,vrm);
    }
  }

  return vrmgroup;    

}
double CollideMD::attempt_collision(int icell, int np, double volume) 
{
  double fnum = update->fnum;
  double dt = update->dt;

  double nattempt;

  if (remainflag) {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
    remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
  } else {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + random->uniform();
  }

  return nattempt;
  }

    
  

double CollideMD::attempt_collision(int icell, int igroup, int jgroup, double volume ) {

double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 // return 2x the value for igroup != jgroup, since no J,I pairing

 double npairs;
 if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 else npairs = ngroup[igroup] * (ngroup[jgroup]);
 //else npairs = 0.5 * ngroup[igroup] * (ngroup[jgroup]);

 nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

 if (remainflag) {
   nattempt += remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else nattempt += random->uniform();

 return nattempt;    

}
int CollideMD::test_collision(int ispecies, int jspecies, int igroup, 
                              Particle::OnePart *ip, Particle::OnePart *jp) {
    // Return 1 if collision should happen, 0 otherwise
    // This is where you'd implement your collision probability logic
    return 1;  // placeholder - always collide for now
}

  void CollideMD::read_param_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open VSS parameter file %s",fname);
    error->one(FLERR,str);
  }

  // set all species diameters to -1, so can detect if not read
  // set all cross-species parameters to -1 to catch no-reads, as
  // well as user-selected average

  for (int i = 0; i < nparams; i++) {
    params[i][i].diam = -1.0;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].diam = params[i][j].omega = params[i][j].tref = -1.0;
      params[i][j].alpha = params[i][j].rotc1 = params[i][j].rotc2 = -1.0;
      params[i][j].rotc3 = params[i][j].vibc1 = params[i][j].vibc2 = -1.0;
      params[i][j].elecc1 = params[i][j].elecc2 = -1.0;
    }
  }

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least REQWORDS, which depends on VARIABLE flag

  int REQWORDS = 5;
  if (relaxflag == VARIABLE) REQWORDS = 11;
  char **words = new char*[REQWORDS+1]; // one extra word in cross-species lines
  char line[MAXLINE];
  int isp,jsp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    int nwords = wordparse(REQWORDS+1,line,words);
    if (nwords < REQWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    jsp = particle->find_species(words[1]);

    // if we don't match a species with second word, but it's not a number,
    // skip the line (it involves a species we aren't using)
    if ( jsp < 0 &&  !(atof(words[1]) > 0) ) continue;

    if (jsp < 0 ) {
      params[isp][isp].diam = atof(words[1]);
      params[isp][isp].omega = atof(words[2]);
      params[isp][isp].tref = atof(words[3]);
      params[isp][isp].alpha = atof(words[4]);
      if (relaxflag == VARIABLE) {
        params[isp][isp].rotc1 = atof(words[5]);
        params[isp][isp].rotc2 = atof(words[6]);
        params[isp][isp].rotc3 = (MY_PI+MY_PI2*MY_PI2)*params[isp][isp].rotc2;
        params[isp][isp].rotc2 = (MY_PI*MY_PIS/2.)*sqrt(params[isp][isp].rotc2);
        params[isp][isp].vibc1 = atof(words[7]);
        params[isp][isp].vibc2 = atof(words[8]);
        params[isp][isp].elecc1 = atof(words[9]);
        params[isp][isp].elecc2 = atof(words[10]);
      }
    }else {
      if (nwords < REQWORDS+1)  // one extra word in cross-species lines
        error->one(FLERR,"Incorrect line format in VSS parameter file");
      params[isp][jsp].diam = params[jsp][isp].diam = atof(words[2]);
      params[isp][jsp].omega = params[jsp][isp].omega = atof(words[3]);
      params[isp][jsp].tref = params[jsp][isp].tref = atof(words[4]);
      params[isp][jsp].alpha = params[jsp][isp].alpha = atof(words[5]);
      if (relaxflag == VARIABLE) {
        params[isp][jsp].rotc1 = params[jsp][isp].rotc1 = atof(words[6]);
        params[isp][jsp].rotc2 = atof(words[7]);
        params[isp][jsp].rotc3 = params[jsp][isp].rotc3 =
                        (MY_PI+MY_PI2*MY_PI2)*params[isp][jsp].rotc2;
        if(params[isp][jsp].rotc2 > 0)
                params[isp][jsp].rotc2 = params[jsp][isp].rotc2 =
                                (MY_PI*MY_PIS/2.)*sqrt(params[isp][jsp].rotc2);
        params[isp][jsp].vibc1 = params[jsp][isp].vibc1= atof(words[8]);
        params[isp][jsp].vibc2 = params[jsp][isp].vibc2= atof(words[9]);
        params[isp][jsp].elecc1 = params[jsp][isp].elecc1= atof(words[10]);
        params[isp][jsp].elecc2 = params[jsp][isp].elecc2= atof(words[11]);
      }
    }
  }

  delete [] words;
  fclose(fp);

  // check that params were read for all species
  for (int i = 0; i < nparams; i++) {

    if (params[i][i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
              particle->species[i].id);
      error->one(FLERR,str);
    }
  }

  for ( int i = 0; i<nparams; i++) {
    params[i][i].mr = particle->species[i].mass / 2;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].mr = params[j][i].mr = particle->species[i].mass *
        particle->species[j].mass / (particle->species[i].mass + particle->species[j].mass);

      if(params[i][j].diam < 0) params[i][j].diam = params[j][i].diam =
                                  0.5*(params[i][i].diam + params[j][j].diam);
      if(params[i][j].omega < 0) params[i][j].omega = params[j][i].omega =
                                   0.5*(params[i][i].omega + params[j][j].omega);
      if(params[i][j].tref < 0) params[i][j].tref = params[j][i].tref =
                                  0.5*(params[i][i].tref + params[j][j].tref);
      if(params[i][j].alpha < 0) params[i][j].alpha = params[j][i].alpha =
                                   0.5*(params[i][i].alpha + params[j][j].alpha);

      if (relaxflag == VARIABLE) {
        if(params[i][j].rotc1 < 0) params[i][j].rotc1 = params[j][i].rotc1 =
                                     0.5*(params[i][i].rotc1 + params[j][j].rotc1);
        if(params[i][j].rotc2 < 0) params[i][j].rotc2 = params[j][i].rotc2 =
                                     0.5*(params[i][i].rotc2 + params[j][j].rotc2);
        if(params[i][j].rotc3 < 0) params[i][j].rotc3 = params[j][i].rotc3 =
                                     0.5*(params[i][i].rotc3 + params[j][j].rotc3);
        if(params[i][j].vibc1 < 0) params[i][j].vibc1 = params[j][i].vibc1 =
                                     0.5*(params[i][i].vibc1 + params[j][j].vibc1);
        if(params[i][j].vibc2 < 0) params[i][j].vibc2 = params[j][i].vibc2 =
                                     0.5*(params[i][i].vibc2 + params[j][j].vibc2);
        if(params[i][j].elecc1 < 0) params[i][j].elecc1 = params[j][i].elecc1 =
                                     0.5*(params[i][i].elecc1 + params[j][j].elecc1);
        if(params[i][j].elecc2 < 0) params[i][j].elecc2 = params[j][i].elecc2 =
                                     0.5*(params[i][i].elecc2 + params[j][j].elecc2);
      }
    }
  }
}




/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
------------------------------------------------------------------------- */

int CollideMD::wordparse(int maxwords, char *line, char **words)
{
  int nwords = 1;
  char * word;

  words[0] = strtok(line," \t\n");
  while ((word = strtok(NULL," \t\n")) != NULL && nwords < maxwords) {
    words[nwords++] = word;
  }
  return nwords;
}

#endif

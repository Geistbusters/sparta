#ifdef COLLIDE_CLASS

CollideStyle(md,CollideMD)

#else



#ifndef SPARTA_COLLIDE_MD_H
#define SPARTA_COLLIDE_MD_H

#include "particle.h"
#include "collide.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <cstring>
#include "collision_data.h"  // Include collision data structure
#include "md_integrator.h"   // Include MD integrator
//#include "debug_utils.h"
#include <string>

	

namespace SPARTA_NS {

#ifdef SPARTA_MD_BUILD
class CollideMD : public Collide {
#else
class CollideMD {  // Standalone version doesn't inherit
#endif

	
public:



//namespace SPARTA_NS {


//class CollideMD : public CollideVSS {
// public:
#ifdef SPARTA_MD_BUILD
  CollideMD(class SPARTA *, int, char **);  // Full SPARTA constructor
  virtual void init();

  struct State {      // two-particle state
    double vr2;
    double vr;
    double imass,jmass;
    double ave_rotdof;
    double ave_vibdof;
    double ave_elecdof;
    double ave_dof;
    double etrans;
    double erot;
    double evib;
    double eelec; 
    double eexchange;
    double eint;
    double etotal;
    double ucmf;
    double vcmf;
    double wcmf;
  };

  struct Params {             // VSS model parameters
    double diam;
    double omega;
    double tref;
    double alpha;
    double rotc1;
    double rotc2;
    double rotc3;
    double vibc1;
    double vibc2;
    double elecc1;
    double elecc2;
    double mr;
  };

// Required pure virtual functions from base Collide class:
  virtual double vremax_init(int, int);
  virtual double attempt_collision(int, int, double);
  virtual double attempt_collision(int, int, int, double);  // overloaded version
  virtual int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *);

  void setup_collision(Particle::OnePart *, Particle::OnePart *) override;

  /**
 * @brief Prepare collision from existing particle data (SPARTA integration mode)
 *
 * Uses predetermined parameters from DSMC particle data:
 * - Relative collision velocity (from particle velocities)
 * - Rovibrational states J, v (from particle internal energy)
 * - Electronic states (from particle electronic energy)
 *
 * Only samples:
 * - Impact parameter (uniform within collision cross section)
 * - Molecular orientations (isotropic)
 *
 * This mode allows CRDS to function as a collision model within SPARTA,
 * providing quantum-state-resolved outcomes for DSMC simulations.
 *
 * @param[in,out] config Configuration parameters
 * @param[out] collision Collision data structure to be populated
 * @param[in] p1 First colliding particle from SPARTA
 * @param[in] p2 Second colliding particle from SPARTA
 * @return true if preparation successful, false if invalid particle data
 *
 * @note Compiled only when SPARTA_MD_BUILD is defined
 * @see prepareFromSampling() for pure QCT mode
 */
bool prepareFromParticles(CRDSConfig& config, CollisionData& collision,
                          Particle::OnePart* p1,
                          Particle::OnePart* p2);



  int perform_collision(Particle::OnePart *&, Particle::OnePart *&, Particle::OnePart *&) override;


    bool updateColl(CollisionData& collision, Particle::OnePart* sparta_p1, 
                       Particle::OnePart* sparta_p2,
                       Particle::OnePart* sparta_p3 = nullptr ) ;

#else  

  void setup_collision(Particle::OnePart *, Particle::OnePart *) ;
  
  int perform_collision(Particle::OnePart *&, Particle::OnePart *&, Particle::OnePart *&) ;
   bool updateColl(CollisionData& collision, Particle::OnePart* sparta_p1, 
                       Particle::OnePart* sparta_p2,
                       Particle::OnePart* sparta_p3 = nullptr ) ;

#endif 
  CollideMD();

virtual	~CollideMD();

  //SPARTA side 
  // Override pure virtual functions from base Collide class
  //double vremax_init(int, int) override;
  //double attempt_collision(int, int, double) override;
  //double attempt_collision(int, int, int, double) override;
  //int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *) override;
  
 // void setup_collision(Particle::OnePart *, Particle::OnePart *) override;
  
 // int perform_collision(Particle::OnePart *&, Particle::OnePart *&, Particle::OnePart *&) override;

  bool runCollision( 
	               int                 integration_type,
		       Particle::OnePart*  sparta_p1,
		       Particle::OnePart*  sparta_p2,
		       Particle::OnePart*  sparta_p3 = nullptr 
			);

  // MD Side of things  
  void setVerbosity(int level);
    void setPESPath(const std::string& path);
    bool initializePES();
    void configureIntegrator(bool surface_hopping_enabled, int initial_surface);

  bool prepare(CRDSConfig& config,CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2, int integration_type);
//    static bool prepare(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2, int integration_type);
  bool prepare(CollisionData& collision,  int integration_type);
  bool collide(MDIntegrator& integrator, CollisionData& collision);
  bool process(MDIntegrator& integrator, CollisionData& collision);

  static void setupOOCollision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2);
  static void setupO2OCollision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2);
  static void setupO2O2Collision(CollisionData& collision, Particle::OnePart* p1, Particle::OnePart* p2);	
  static void testSetupO2OCollision(CollisionData& collision);   
  static int analyzeOutcome(const CollisionData& collision);

  double sample_impact_parameter();
  
 private:
    MDIntegrator integrator_;
    CRDSConfig config_;
    int verbosity_;
    std::string pes_path_;
    int integration_type_;



////////// VSS Stuff





#ifdef SPARTA_MD_BUILD
 protected:
  int relaxflag,eng_exchange;
  double vr_indice;
  double **prefactor; // static portion of collision attempt frequency

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  Params **params;             // VSS params for each species
  int nparams;                // # of per-species params read in


  double sample_bl(RanKnuth *, double, double);
  double rotrel (int, double);
  double vibrel (int, double);
  double elecrel (int, double);

  void read_param_file(char *);
  int wordparse(int, char *, char **);
#endif  

};

}

//#ifdef COLLIDE_CLASS
//CollideStyle(md,CollideMD)
#endif

#endif

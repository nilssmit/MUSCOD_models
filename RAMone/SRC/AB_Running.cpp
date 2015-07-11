#include <cmath>
#include "def_usrmod.hpp"
#include <AB_Constraints.h>


#define  NMOS   3  /* Number of phases (MOdel Stages) */
/* There are three different stages for the symmetrically running prismatic biped:
 * 0) Flight until touchdown left (FlowMap)
 * 1) Touchdown collision of the left leg (JumpMap)
 * 2) Stance left (FlowMap)
 */
#define  NXA    0  /* Number of algebraic states */
#define  NPR    0  /* Number of local parameters */

/** \brief Entry point for the muscod application */
extern "C" void def_model(void);
void def_model(void)
{
	/* Define problem dimensions */
	def_mdims(NMOS, NPFree, rcfcn, rcfcne);
	/* Define the first phase */
	/* def_mstage(I, NXD, NXA, NU, mfcn, lfcn, jacmlo, jacmup, astruc, afcn, ffcn, gfcn, rwh, iwh)
	 * Call to define a model stage with index I, where
	 * NXD is the differential state dimension,
	 * NXA the algebraic state dimension, and
	 * NU is the control dimension.
	 * mfcn is a pointer to a Mayer term function (or NULL) to be evaluated at the end of the stage, and
	 * lfcn a pointer to a Lagrange term (or NULL).
	 * For documentation of the left-hand side matrix function afcn, and of the integers jacmlo,
	 * jacmup, and astruc that provide structural matrix information please consult the
	 * DAESOL-manual [BBS99]; setting the integers to zero is equivalent to not defining
	 * any structural information.
	 * ffcn is a pointer to the differential right hand side function,
	 * gfcn the pointer to the algebraic right hand side function (or NULL).
	 * rwh, iwh are real and integer work arrays which can be used to pass a common workspace to the stage functions.
	 *
	 */
	def_mstage( 0, // 0) Flight until touchdown left (FlowMap)
				NY, NXA, NU,
				NULL, NULL,
				0, 0, 0, NULL, ffcn_flight, NULL,
				NULL, NULL
				);
	def_mstage( 1, // 1) Touchdown collision of the left leg (JumpMap)
				NY, NXA, NU,
				NULL, NULL,
				0, 0, 0, NULL, ffcn_collisionLEFT, NULL,
				NULL, NULL
				);
	def_mstage( 2, // 2) Stance left (FlowMap)
				NY, NXA, NU,
				mfcn_COT, NULL,
				0, 0, 0, NULL, ffcn_stanceLEFT, NULL,
				0, NULL
				);

	/* Define constraints at the start point
	 * i.e., the coupled constraints for periodicity and the conditions for liftoff of the right
	 * leg (after the left leg was already in the air)	 */
	def_mpc(0, "Start Point", NPR, rdfcn_singleLeg_n, rdfcn_singleLeg_ne, rdfcn_liftoffRIGHT_lR, rcfcn_beginning);
	/* Define constraints throughout the flight phase:
	 * i.e., the decoupled constraints for ground clearance and actuator speed and force	 */
	def_mpc(0, "Interior Point", NPR, rdfcn_neConstraints_n, rdfcn_neConstraints_ne, rdfcn_neConstraints_lr, NULL);
	/* Define constraints at the start of phase 1,
	 * i.e., the conditions for touchdown of the left leg (with both legs being in the air)	 */
	def_mpc(1,"Start Point", NPR, rdfcn_singleLeg_n, rdfcn_singleLeg_ne, rdfcn_touchdownLEFT_lr, NULL);
	/* Define constraints throughout the touchdown collision:
	 * i.e., the decoupled constraints for ground clearance and actuator speed and force	 */
	def_mpc(1, "Interior Point", NPR, rdfcn_neConstraints_n, rdfcn_neConstraints_ne, rdfcn_neConstraints_lr, NULL);
	/* Define constraints at the beginning of the left stance phase:
	 * i.e., the decoupled constraints for ground clearance and actuator speed and force	 */
	def_mpc(2, "Start Point", NPR, rdfcn_neConstraints_n, rdfcn_neConstraints_ne, rdfcn_neConstraints_Lr, NULL);
	/* Define constraints throughout the left stance phase:
	 * i.e., the decoupled constraints for ground clearance and actuator speed and force	 */
	def_mpc(2, "Interior Point", NPR, rdfcn_neConstraints_n, rdfcn_neConstraints_ne, rdfcn_neConstraints_Lr, NULL);
	/* Define constraints at the end point of phase 2
	 * i.e., the coupled constraints for periodicity and the conditions for the average velocity	 */
	def_mpc(2, "End Point", NPR, rdfcn_avgSpeed_n, rdfcn_avgSpeed_ne, rdfcn_avgSpeed, rcfcn_endSymmetric);

	// Create output to a .mot and a .plt file
	def_mio (NULL , motion_output, plot_output);
}

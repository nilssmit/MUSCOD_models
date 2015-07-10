/*
 * PB_Constraints.h
 *
 * Common Constraint Header for a Prismatic Biped.
 *
 *  Created on: April 14, 2013
 *      Author: cdremy@umich.edu
 */

#include "PB_Dynamics.h"
#include "def_usrmod.hpp"

#ifndef PB_DDH_CONSTRAINTS_H_
#define PB_DDH_CONSTRAINTS_H_

// Create the instance of the model object
static PB_Dynamics syst;

/*********************************
 * FlowMaps
 *********************************/
// The Flow Maps are used as differential right hand side functions
// Both legs in the air
void ffcn_flight(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);
// Left leg in the air, right on the ground
void ffcn_stanceRIGHT(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);
// Right leg in the air, left on the ground
void ffcn_stanceLEFT(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);
// Both legs on the ground
void ffcn_stanceBOTH(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);

/*********************************
 * JumpMaps
 *********************************/
// The Jump Maps are used as differential right hand side functions
// Touchdown collision of the right leg
void ffcn_collisionRIGHT(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);
// Touchdown collision of the left leg
void ffcn_collisionLEFT(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);
// Touchdown collision of both legs (the outcome of the collisions doesn't depend on the previous contact configuration, so we just need one function)
void ffcn_collisionBOTH(double *t, double *xd, double *xa, double *u, double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info);

/*********************************
 * JumpSets
 *********************************/
// The Jummp Sets are used as start and end point constraints
// In addition to the constraints for lift-off and touchdown,
// constraints on ground clearance (2), actuator torque(8) and
// actuator speed (8) are enforced.  They can be turned on and
// off via the free parameter const_sel, that takes binary flags
// as follows:
// 1 ground clearance constraint on
// 2 actuator torque constraint on
// 4 actuator velocity constraint on
static int rdfcn_singleLeg_n  = 19;
static int rdfcn_singleLeg_ne = 1;
// Since the detection of touchdown and liftoff depends on the current
// contact configuration, every function is provided with two versions:
// capital means, the leg is on the ground, small it is in the air
void rdfcn_liftoffRIGHT_LR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_liftoffRIGHT_lR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_liftoffLEFT_LR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_liftoffLEFT_Lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_touchdownRIGHT_lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_touchdownRIGHT_Lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_touchdownLEFT_lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_touchdownLEFT_lR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
static int rdfcn_bothLegs_n  = 20;
static int rdfcn_bothLegs_ne = 2;
// Liftoff of both legs at once
void rdfcn_liftoffBOTH(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
// Touchdown of both legs at once
void rdfcn_touchdownBOTH(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
// Liftoff of the left leg, while the right leg touches the ground
void rdfcn_switchRIGHTLEFT(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
// Liftoff of the right leg, while the left leg touches the ground
void rdfcn_switchLEFTRIGHT(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);

/*********************************
 * Periodicity Constraints
 *********************************/
static int rcfcn  = NY + NU - 1;  /* Number of coupled constraints for periodicity */
static int rcfcne = NY + NU - 1;  /* Number of coupled equality constraints for periodicity  */
enum constNames { dx_const, y_const, dy_const, phi_const, dphi_const,
                  alphaL_const, dalphaL_const, ualphaL_const, dualphaL_const,
                  alphaR_const, dalphaR_const, ualphaR_const, dualphaR_const,
                  lL_const, dlL_const, ulL_const, dulL_const,
                  lR_const, dlR_const, ulR_const, dulR_const,
                  TalphaL_const, TalphaR_const,
                  FlL_const, FlR_const
};

// Periodicity constraint at the beginning of the stride
void rcfcn_beginning(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
// Periodicity constraint at the middle of the stride that switches left and right leg
void rcfcn_endSymmetric(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
// Periodicity constraint at the end of the stride
void rcfcn_endPeriodic(double *ts, double *sd, double *sa, double *u, double *p, double *pr, double *res, long *dpnd, InfoPtr *info);

/*********************************
 * End point Constraints
 *********************************/
// Enforce the desired average speed:
static int rdfcn_avgSpeed_n  = 1;
static int rdfcn_avgSpeed_ne = 1;
void rdfcn_avgSpeed(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);

/*********************************
 * Interior point Constraints
 *********************************/
// Constraints on ground clearance (2), actuator torque(8) and
// actuator speed (8) are enforced.  They can be turned on and
// off via the free parameter const_sel, that takes binary flags
// as follows:
// 1 ground clearance constraint on
// 2 actuator torque constraint on
// 4 actuator velocity constraint on
static int rdfcn_neConstraints_n  = 18;
static int rdfcn_neConstraints_ne = 0;
void rdfcn_neConstraints_LR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_neConstraints_lR(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_neConstraints_Lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);
void rdfcn_neConstraints_lr(double *ts, double *sd, double *sa, double *u,  double *p, double *pr, double *res, long *dpnd, InfoPtr *info);



/*********************************
 * Objective Function
 *********************************/
void mfcn_COT(double *ts, double *sd, double *sa, double *u, double *p, double *mval, long *dpnd, InfoPtr *info);

/*********************************
 * Output Functions
 *********************************/
static const char *figureview_file_name = "PB_DDH_MOTION.fv";
void figureview_output
(
		long   *imos,      ///< index of model stage (I)
		long   *imsn,      ///< index of m.s. node on current model stage (I)
		double *ts,        ///< time at m.s. node (I)
		double *te,        ///< time at end of m.s. interval (I)
		double *sd,        ///< differential states at m.s. node (I)
		double *sa,        ///< algebraic states at m.s. node (I)
		double *u,         ///< controls at m.s. node (I)
		double *udot,      ///< control slopes at m.s. node (I)
		double *ue,        ///< controls at end of m.s. interval (I)
		double *uedot,     ///< control slopes at end of m.s. interval (I)
		double *p,         ///< global model parameters (I)
		double *pr,        ///< local i.p.c. parameters (I)
		double *ccxd,
		double *mul_ccxd,  ///< multipliers of continuity conditions (I)
#if defined(PRSQP) || defined(EXTPRSQP)
		double *ares,
		double *mul_ares,
#endif
		double *rd,
		double *mul_rd,    ///< multipliers of decoupled i.p.c. (I)
		double *rc,
		double *mul_rc,    ///< multipliers of coupled i.p.c. (I)
		double *obj,
		double *rwh,       ///< real work array (I)
		long   *iwh        ///< integer work array (I)
);
static const char *motion_file_name = "PB_MOTION.txt";
static const char *plot_file_name = "PB_PLOT.txt";
void motion_output
(
		long   *imos,      ///< index of model stage (I)
		long   *imsn,      ///< index of m.s. node on current model stage (I)
		double *ts,        ///< time at m.s. node (I)
		double *te,        ///< time at end of m.s. interval (I)
		double *sd,        ///< differential states at m.s. node (I)
		double *sa,        ///< algebraic states at m.s. node (I)
		double *u,         ///< controls at m.s. node (I)
		double *udot,      ///< control slopes at m.s. node (I)
		double *ue,        ///< controls at end of m.s. interval (I)
		double *uedot,     ///< control slopes at end of m.s. interval (I)
		double *p,         ///< global model parameters (I)
		double *pr,        ///< local i.p.c. parameters (I)
		double *ccxd,
		double *mul_ccxd,  ///< multipliers of continuity conditions (I)
#if defined(PRSQP) || defined(EXTPRSQP)
		double *ares,
		double *mul_ares,
#endif
		double *rd,
		double *mul_rd,    ///< multipliers of decoupled i.p.c. (I)
		double *rc,
		double *mul_rc,    ///< multipliers of coupled i.p.c. (I)
		double *obj,
		double *rwh,       ///< real work array (I)
		long   *iwh        ///< integer work array (I)
);
void plot_output(
  double *t,      ///< time (I)
  double *xd,     ///< differential states (I)
  double *xa,     ///< algebraic states (I)
  double *u,      ///< controls (I)
  double *p,      ///< global model parameters (I)
  double *rwh,    ///< real work array (I)
  long   *iwh,     ///< integer work array (I)
  InfoPtr *info
);

/*********************************
 * Misc Functionality
 *********************************/
// PFree
#define  NPFree 9 /* Number of adjustable system parameters */
enum pFreeNames { hipJointType_free, legJointType_free,
	              cost_fct_sel_free, v_avg_free,
	              const_sel_free,
				  kalpha_free, du_max_alpha_free,
		          kl_free, du_max_l_free
};
void convertParameters(double *p_free, double *p);
void computeConstraints(int const_selFLAG, int offset, double *res);

#endif /* PB_DDH_CONSTRAINTS_H_ */

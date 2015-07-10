/*
 * AB_Constraints.cpp
 *
 * Common Constraint code for an Articulating Biped.
 *
 *  Created on: May 15, 2015
 *      Author: nilssmit@umich.edu
 */

#include "AB_Constraints.h"
#include <stdio.h>
#include "model.hpp"

/*********************************
 * FlowMaps
 *********************************/
// Both legs in the air
void ffcn_flight(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Flow Map:
	syst.FlowMap(rhs);
}
// Left leg in the air, right on the ground
void ffcn_stanceRIGHT(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Flow Map:
	syst.FlowMap(rhs);
}
// Right leg in the air, left on the ground
void ffcn_stanceLEFT(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Flow Map:
	syst.FlowMap(rhs);
}
// Both legs on the ground
void ffcn_stanceBOTH(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Flow Map:
	syst.FlowMap(rhs);
}

/*********************************
 * JumpMaps
 *********************************/
// The outcome of the collisions doesn't depend on the previous
// contact configuration, so we just need one function per touchdown event
// Touchdown collision of the right leg
void ffcn_collisionRIGHT(double *t, double *xd, double *xa, double *u,
		double *p_free, double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Jump Map:
	syst.JumpMap(rhs);
}
// Touchdown collision of the left leg
void ffcn_collisionLEFT(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Jump Map:
	syst.JumpMap(rhs);
}
// Touchdown collision of both legs
void ffcn_collisionBOTH(double *t, double *xd, double *xa, double *u, double *p_free,
		double *rhs, double *rwh, long *iwh, InfoPtr *info) {
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(xd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	// Compute Jump Map:
	syst.JumpMap(rhs);
}

/*********************************
 * JumpSets
 *********************************/
// Since the detection of touchdown and liftoff depends on the current
// contact configuration, every function is provided with two versions:
// capital means, the leg is on the ground, small it is in the air
void rdfcn_liftoffRIGHT_LR(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[liftoffR];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_liftoffRIGHT_lR(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[liftoffL];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_liftoffLEFT_LR(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[liftoffL];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_liftoffLEFT_Lr(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[liftoffL];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_touchdownRIGHT_lr(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[touchdownR];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_touchdownRIGHT_Lr(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[touchdownR];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_touchdownLEFT_lr(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[touchdownL];
	computeConstraints(p_free[const_sel_free], 1, res);
}
void rdfcn_touchdownLEFT_lR(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[touchdownL];
	computeConstraints(p_free[const_sel_free], 1, res);
}
// Liftoff of both legs at once
void rdfcn_liftoffBOTH(double *ts, double *sd, double *sa, double *u, double *p_free,
		double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[liftoffL];
	res[1] = eventVal[liftoffR];
	computeConstraints(p_free[const_sel_free], 2, res);
}
// Touchdown of both legs at once
void rdfcn_touchdownBOTH(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = flight;
	// Update States:
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	double eventVal[NEV];
	syst.JumpSet(eventVal);
	res[0] = eventVal[touchdownL];
	res[1] = eventVal[touchdownR];
	computeConstraints(p_free[const_sel_free], 2, res);
}

/*********************************
 * Periodicity Constraints
 *********************************/
// Periodicity constraint at the beginning of the stride
void rcfcn_beginning(double *ts, double *sd, double *sa, double *u, double *p_free,
		double *pr, double *res, long *dpnd, InfoPtr *info) {
	// Sum all periodic values in positive direction
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
		return;
	}
	res[dx_const]       = sd[dx];
	res[y_const]        = sd[y];
	res[dy_const]       = sd[dy];
	res[phi_const]      = sd[phi];
	res[dphi_const]     = sd[dphi];
	res[alphaL_const]   = sd[alphaL];
	res[dalphaL_const]  = sd[dalphaL];
	res[ualphaL_const]  = sd[ualphaL];
	res[dualphaL_const] = sd[dualphaL];
	res[alphaR_const]   = sd[alphaR];
	res[dalphaR_const]  = sd[dalphaR];
	res[ualphaR_const]  = sd[ualphaR];
	res[dualphaR_const] = sd[dualphaR];
	res[betaL_const]       = sd[betaL];
	res[dbetaL_const]      = sd[dbetaL];
	res[ubetaL_const]      = sd[ubetaL];
	res[dubetaL_const]     = sd[dubetaL];
	res[betaR_const]       = sd[betaR];
	res[dbetaR_const]      = sd[dbetaR];
	res[ubetaR_const]      = sd[ubetaR];
	res[dubetaR_const]     = sd[dubetaR];
	res[TalphaL_const]  = u[TalphaL];
	res[TalphaR_const]  = u[TalphaR];
	res[TbetaL_const]      = u[TbetaL];
	res[TbetaR_const]      = u[TbetaR];
}
// Periodicity constraint at the middle of the stride that switches left and right leg
void rcfcn_endSymmetric(double *ts, double *sd, double *sa, double *u,
		double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
		return;
	}
	// Sum all periodic values in negative direction
	// in the constraints, we have to flip left and right:
	res[dx_const]       = -sd[dx];
	res[y_const]        = -sd[y];
	res[dy_const]       = -sd[dy];
	res[phi_const]      = -sd[phi];
	res[dphi_const]     = -sd[dphi];
	res[alphaL_const]   = -sd[alphaR];
	res[dalphaL_const]  = -sd[dalphaR];
	res[ualphaL_const]  = -sd[ualphaR];
	res[dualphaL_const] = -sd[dualphaR];
	res[alphaR_const]   = -sd[alphaL];
	res[dalphaR_const]  = -sd[dalphaL];
	res[ualphaR_const]  = -sd[ualphaL];
	res[dualphaR_const] = -sd[dualphaL];
	res[betaL_const]       = -sd[betaR];
	res[dbetaL_const]      = -sd[dbetaR];
	res[ubetaL_const]      = -sd[ubetaR];
	res[dubetaL_const]     = -sd[dubetaR];
	res[betaR_const]       = -sd[betaL];
	res[dbetaR_const]      = -sd[dbetaL];
	res[ubetaR_const]      = -sd[ubetaL];
	res[dubetaR_const]     = -sd[dubetaL];
	res[TalphaL_const]  = -u[TalphaR];
	res[TalphaR_const]  = -u[TalphaL];
	res[TbetaL_const]      = -u[TbetaR];
	res[TbetaR_const]      = -u[TbetaL];
}
// Periodicity constraint at the end of the stride
void rcfcn_endPeriodic(double *ts, double *sd, double *sa, double *u, double *p_free,
		double *pr, double *res, long *dpnd, InfoPtr *info) {
	// Sum all periodic values in negative direction
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
		return;
	}
	res[dx_const]       = -sd[dx];
	res[y_const]        = -sd[y];
	res[dy_const]       = -sd[dy];
	res[phi_const]      = -sd[phi];
	res[dphi_const]     = -sd[dphi];
	res[alphaL_const]   = -sd[alphaL];
	res[dalphaL_const]  = -sd[dalphaL];
	res[ualphaL_const]  = -sd[alphaL];
	res[dualphaL_const] = -sd[dalphaL];
	res[alphaR_const]   = -sd[alphaR];
	res[dalphaR_const]  = -sd[ualphaR];
	res[ualphaR_const]  = -sd[ualphaR];
	res[dualphaR_const] = -sd[dualphaR];
	res[betaL_const]       = -sd[betaL];
	res[dbetaL_const]      = -sd[dbetaL];
	res[ubetaL_const]      = -sd[ubetaL];
	res[dubetaL_const]     = -sd[dubetaL];
	res[betaR_const]       = -sd[betaR];
	res[dbetaR_const]      = -sd[dbetaR];
	res[ubetaR_const]      = -sd[ubetaR];
	res[dubetaR_const]     = -sd[dubetaR];
	res[TalphaL_const]  = -u[TalphaL];
	res[TalphaR_const]  = -u[TalphaR];
	res[TbetaL_const]      = -u[TbetaL];
	res[TbetaR_const]      = -u[TbetaR];
}

/*********************************
 * Other Constraints
 *********************************/
// Enforce the desired average speed:
void rdfcn_avgSpeed(double *ts, double *sd, double *sa, double *u, double *p_free,
		double *pr, double *res, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = RFCN_DPND(*ts, *sd, 0, 0, *p_free, 0);
		return;
	}
	res[0] = sd[x] - p_free[v_avg_free] * (*ts);
}

/*********************************
 * Interior point Constraints
 *********************************/
// Ground clearances:
// Since the computation of the constraints depends on the current contact configuration,
// four different configurations exist:
// capital means, the leg is on the ground, small it is in the air

void rdfcn_neConstraints_LR(double *ts, double *sd, double *sa, double *u,  double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info){
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = stance;
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	computeConstraints(p_free[const_sel_free], 0, res);
}
void rdfcn_neConstraints_lR(double *ts, double *sd, double *sa, double *u,  double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info){
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = stance;
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	computeConstraints(p_free[const_sel_free], 0, res);
}

void rdfcn_neConstraints_Lr(double *ts, double *sd, double *sa, double *u,  double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info){
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = stance;
	zVecPtr[phaseR] = flight;
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	computeConstraints(p_free[const_sel_free], 0, res);
}
void rdfcn_neConstraints_lr(double *ts, double *sd, double *sa, double *u,  double *p_free, double *pr, double *res, long *dpnd, InfoPtr *info){
	if (*dpnd) {
		*dpnd = RFCN_DPND(0, *sd, 0, *u, *p_free, 0);
		return;
	}
	// Set phases:
	double zVecPtr[NZ];
	zVecPtr[phaseL] = flight;
	zVecPtr[phaseR] = flight;
	syst.setContState(sd);
	syst.setDiscState(zVecPtr);
	double p[NP];
	convertParameters(p_free,p);
	syst.setSystParam(p);
	syst.setExctState(u);
	computeConstraints(p_free[const_sel_free], 0, res);
}




/*********************************
 * Objective Functions
 *********************************/
void mfcn_COT(double *ts, double *sd, double *sa, double *p_free, double *pr,
		double *mval, long *dpnd, InfoPtr *info) {
	if (*dpnd) {
		*dpnd = MFCN_DPND(0, *sd, 0, *p_free, 0);
		return;
	}
	switch(int(p_free[cost_fct_sel_free])){
	case 1: // Actuator Cost Of Transport
		if (sd[x] > 0.0001)
			*mval = sd[posActWork] / sd[x];
		else
			*mval = sd[posActWork];
		break;
	case 2: // Electrical Cost Of Transport
		if (sd[x] > 0.0001)
			*mval = sd[posElWork] / sd[x];
		else
			*mval = sd[posElWork];
		break;
	case 3: // Electrical Loss Cost Of Transport
		if (sd[x] > 0.0001)
			*mval = sd[totElLoss] / sd[x];
		else
			*mval = sd[totElLoss];
		break;
	default:
		cout << "Unknown parameter value for ""cost_fct_sel""." << endl;
		cout << "Must be 1 for posActWork" << endl;
		cout << "        2 for posElWork" << endl;
		cout << "        3 for totElLoss" << endl;
		break;
	}
}


/*********************************
 * Output Functions
 *********************************/
void figureview_output(long *imos, ///< index of model stage (I)
		long *imsn, ///< index of m.s. node on current model stage (I)
		double *ts, ///< time at m.s. node (I)
		double *te, ///< time at end of m.s. interval (I)
		double *sd, ///< differential states at m.s. node (I)
		double *sa, ///< algebraic states at m.s. node (I)
		double *u, ///< controls at m.s. node (I)
		double *udot, ///< control slopes at m.s. node (I)
		double *ue, ///< controls at end of m.s. interval (I)
		double *uedot, ///< control slopes at end of m.s. interval (I)
		double *p_free, ///< global model parameters (I)
		double *pr, ///< local i.p.c. parameters (I)
		double *ccxd, double *mul_ccxd, ///< multipliers of continuity conditions (I)
#if defined(PRSQP) || defined(EXTPRSQP)
		double *ares, double *mul_ares,
#endif
		double *rd, double *mul_rd, ///< multipliers of decoupled i.p.c. (I)
		double *rc, double *mul_rc, ///< multipliers of coupled i.p.c. (I)
		double *obj, double *rwh, ///< real work array (I)
		long *iwh ///< integer work array (I)
		) {
	FILE *fp;
	stringstream output_filename("");
	const char* fv_header_state_names = "# Time\n"
			"# LowerTrunk_TX\n"
			"# LowerTrunk_TY\n"
			"# LowerTrunk_RZ\n"
			"# UpperLeg_R_RZ\n"
			"# UpperLeg_R_TY\n"
			"# UpperLeg_L_RZ\n"
			"# UpperLeg_L_TY\n"
			"# LowerLeg_R_RZ\n"
			"# LowerLeg_R_TY\n"
			"# LowerLeg_L_RZ\n"
			"# LowerLeg_L_TY\n";
	// Filename
	output_filename << "./RES/" << figureview_file_name;

	/* Create and reset file when we're in the first MS-knot.
	 * otherwise just append data
	 */
	if ((*imsn == 0) && (*imos == 0)) {
		fp = fopen(output_filename.str().c_str(), "w");

		if (!fp) {
			fprintf(stderr, "Error: Could not open file '%s'!\n",
					output_filename.str().c_str());
			return;
		}

		// Write the header of the figureview file

		fprintf(
				fp,
				"# JAFV Format: DOF_0 DOF_1 ... DOF_N TAG0_x TAG0_y ... TAGN_z\n");
		fprintf(fp, "# Format for cal3dview:\n");

		fprintf(fp, "# MODEL hopper\n");
		fprintf(fp, "# SYMMETRIC_DATA\n");
		fprintf(fp, "# BEGIN_FORMAT\n");
		// the degrees of freedom
		fprintf(fp, "%s", fv_header_state_names);
		fprintf(fp, "# END_FORMAT\n");
		// More header information
		fprintf(
				fp,
				"#\n#\n# Articulating biped hopper\n#\n# Variables per frame: \n#\n# %d\n#\n",
				1 + NQ + 2);
	} else {
		fp = fopen(output_filename.str().c_str(), "a");

		if (!fp) {
			fprintf(stderr, "Error: Could not open file '%s'!\n",
					output_filename.str().c_str());
			return;
		}
	}

	// Write a line of data:
	fprintf(fp, "%e\t", *ts);
	fprintf(fp, "%e\t", sd[x]);
	fprintf(fp, "%e\t", sd[y]);
	fprintf(fp, "%e\t", sd[phi]);
	fprintf(fp, "%e\t", sd[alphaR]);
	fprintf(fp, "0\t");
	fprintf(fp, "%e\t", sd[alphaL]);
	fprintf(fp, "0\t");
	fprintf(fp, "0\t");
	fprintf(fp, "%e\t", 1 - sd[betaR]);
	fprintf(fp, "0\t");
	fprintf(fp, "%e\t", 1 - sd[betaL]);
	fprintf(fp, "\n");
	// Close file
	fclose(fp);
}

void motion_output(long *imos, ///< index of model stage (I)
			long *imsn, ///< index of m.s. node on current model stage (I)
			double *ts, ///< time at m.s. node (I)
			double *te, ///< time at end of m.s. interval (I)
			double *sd, ///< differential states at m.s. node (I)
			double *sa, ///< algebraic states at m.s. node (I)
			double *u, ///< controls at m.s. node (I)
			double *udot, ///< control slopes at m.s. node (I)
			double *ue, ///< controls at end of m.s. interval (I)
			double *uedot, ///< control slopes at end of m.s. interval (I)
			double *p_free, ///< global model parameters (I)
			double *pr, ///< local i.p.c. parameters (I)
			double *ccxd, double *mul_ccxd, ///< multipliers of continuity conditions (I)
	#if defined(PRSQP) || defined(EXTPRSQP)
			double *ares, double *mul_ares,
	#endif
			double *rd, double *mul_rd, ///< multipliers of decoupled i.p.c. (I)
			double *rc, double *mul_rc, ///< multipliers of coupled i.p.c. (I)
			double *obj, double *rwh, ///< real work array (I)
			long *iwh ///< integer work array (I)
			) {
	FILE *fp;
	stringstream output_filename("");
	// Filename
	char problemname[50];
	get_pname(problemname);
	output_filename << "./RES/" << problemname <<".mot";
	if ((*imsn == 0) && (*imos == 0))
				fp = fopen(output_filename.str().c_str(), "w");
			else
				fp = fopen(output_filename.str().c_str(), "a");
			if (!fp) {
				fprintf(stderr, "Error: Could not open file '%s'!\n",
						output_filename.str().c_str());
				return;
			}
		/* Create and reset file when we're in the first MS-knot.
		 * otherwise just append data
		 */
		if ((*imsn == 0) && (*imos == 0)) {
			// write header:
			/// Get names of model quantities.
			char **xd_name[1];  /* differential state names (O) */
			char **xa_name[1];  /* algebraic state names (O) */
			char **u_name[1];   /* control names (O) */
			char **h_name[1];   /* model stage duration names (O) */
			char **p_name[1];   /* global model parameter names (O) */
			char *of_name[1];   /* objective name (O) */
			get_names(xd_name, xa_name, u_name, h_name, p_name, of_name);
			// Objective:
			fprintf(fp, "COT: %e\n", *obj);
			// Print parameters first:
			fprintf(fp, "//");
			for (int i = 0; i<NPFree; i++){
				fprintf(fp, "%s\t",p_name[0][i]);
			}
			fprintf(fp, "\n");
			for (int i = 0; i<NPFree; i++){
				fprintf(fp, "%e\t",p_free[i]);
			}
			fprintf(fp, "\n");
			// States and timing:
			fprintf(fp, "//imos\t");
			fprintf(fp, "imsn\t");
			fprintf(fp, "time\t");
			for (int i = 0; i<NY; i++){
				fprintf(fp, "%s\t",xd_name[0][i]);
			}
			for (int i = 0; i<NU; i++){
				fprintf(fp, "%s\t",u_name[0][i]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "%li\t", *imos);
		fprintf(fp, "%li\t", *imsn);
		fprintf(fp, "%e\t", *ts);
		for (int i = 0; i<NY; i++){
			fprintf(fp, "%e\t",sd[i]);
		}
		for (int i = 0; i<NU; i++){
			fprintf(fp, "%e\t",u[i]);
		}
		fprintf(fp, "\n");
		// Close file
		fclose(fp);
}

void plot_output(
  double *t,      ///< time (I)
  double *xd,     ///< differential states (I)
  double *xa,     ///< algebraic states (I)
  double *u,      ///< controls (I)
  double *p_free,      ///< global model parameters (I)
  double *rwh,    ///< real work array (I)
  long   *iwh,     ///< integer work array (I)
  InfoPtr *info
){
	FILE *fp;
	stringstream output_filename("");
	// Filename
	char problemname[50];
	get_pname(problemname);
	output_filename << "./RES/" << problemname <<".plt";
	if (*t == 0)
		fp = fopen(output_filename.str().c_str(), "w");
	else
		fp = fopen(output_filename.str().c_str(), "a");
	if (!fp) {
		fprintf(stderr, "Error: Could not open file '%s'!\n",
				output_filename.str().c_str());
		return;
	}
	if (*t == 0) {
		// write header:
		/// Get names of model quantities.
		char **xd_name[1];  /* differential state names (O) */
		char **xa_name[1];  /* algebraic state names (O) */
		char **u_name[1];   /* control names (O) */
		char **h_name[1];   /* model stage duration names (O) */
		char **p_name[1];   /* global model parameter names (O) */
		char *of_name[1];   /* objective name (O) */
		get_names(xd_name, xa_name, u_name, h_name, p_name, of_name);
		fprintf(fp, "//time\t");
		for (int i = 0; i<NY; i++){
			fprintf(fp, "%s\t",xd_name[0][i]);
		}
		for (int i = 0; i<NU; i++){
			fprintf(fp, "%s\t",u_name[0][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "%e\t", *t);
	for (int i = 0; i<NY; i++){
		fprintf(fp, "%e\t",xd[i]);
	}
	for (int i = 0; i<NU; i++){
		fprintf(fp, "%e\t",u[i]);
	}

	fprintf(fp, "\n");
	// Close file
	fclose(fp);
}


/*********************************
 * Misc Functions
 *********************************/
void convertParameters(double *p_free, double *p){
	/* this is a quite tricky function, since it maps a set of free
	 * parameters (which can be changed in MUSCOD) to the full parameter
	 * vector of the dynamic representation (or the other way around). We
	 * hence have to be very careful with the two different parameter enums
	 */
	p[cost_fct_sel]   = p_free[cost_fct_sel_free];
	p[v_avg]          = p_free[v_avg_free];
    p[sigma]          = p_free[sigma_free];
	p[const_sel_free] = p_free[const_sel];
	
    p[hip_jnt_type] = 0;
    p[leg_jnt_type] = 0;
	p[g]     		= 9.81;
	p[m1]    		= 6.54;
	p[m2]    		= 1.56;
	p[m3]    		= 0.42;
	p[lH]    		= 0.137675;
	p[lL1]   		= 0.2;
	p[lL2]   		= 0.2385;
	p[l2]    		= 0.016;
	p[l3]    		= 0.1585;
	p[rFoot] 		= 0.0282;
	p[j1_]   		= 0.04;
	p[j2]    		= 0.0085;
	p[j3]    		= 0.0035;
	p[kalpha] 		= 70;
	p[balpha] 		= 7;
	p[kbeta1] 		= 80;
	p[kbeta2] 		= 36;
	p[bbeta1] 	    = 10;
	p[bbeta2] 	    = 1;
	p[betaSmall] 	= 0.05;
	p[P_max_alpha]  = 1860;
	p[du_max_alpha] = 7.84;
	p[c_lim_alpha]  = 0.0707;
	p[j_mot_alpha]  = 0.3025;
	p[P_max_beta]   = 1860;
	p[du_max_beta]  = 7.01;
	p[c_lim_beta]   = 0.0707;
	p[j_mot_beta]   = 0.378;
}

// Computes the inequality constraints on ground clearance (2),
// actuator torque(8) and actuator speed (8).  They can be turned on
// and off via the free parameter const_selFLAG, that takes binary flags
// as follows:
// 1 ground clearance constraint on
// 2 actuator torque constraint on
// 4 actuator velocity constraint on
// 'offset' defines, how many equality constraints are used in the res-vector
void computeConstraints(int const_selFLAG, int offset, double *res){
	double gcL;
	double gcR;
	double du_max[8];
	double F_max[8];
	syst.Constraints(&gcL, &gcR, du_max, F_max);
	// map to residuals:
	if (const_selFLAG & 1){
		res[0+offset] = gcL;
		res[1+offset] = gcR;
	} else {
		res[0+offset] = 1;
		res[1+offset] = 1;
	}
	if (const_selFLAG & 2){
		for (int i = 0; i<8; i++){
			res[i+2+offset] = F_max[i];
		}
	} else {
		for (int i = 0; i<8; i++){
			res[i+2+offset] = 1;
		}

	}
	if (const_selFLAG & 4){
		for (int i = 0; i<8; i++){
			res[i+10+offset] = du_max[i];
		}
	} else {
		for (int i = 0; i<8; i++){
			res[i+10+offset] = 1;
		}
	}
}

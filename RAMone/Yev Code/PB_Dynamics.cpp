/*
 * PB_Dynamics.cpp
 *
 * Common Dynamic code for a Prismatic Biped.
 *
 *  Created on: April 13, 2013
 *      Author: cdremy@umich.edu
 */

#include "PB_Dynamics.h"
#include <algorithm>

PB_Dynamics::PB_Dynamics() {}  // Constructor (Empty)
PB_Dynamics::~PB_Dynamics() {} // Destructor (Empty)

/* Hybrid Dynamic Equations */
/* The FlowMap computes the right hand side of the governing differential
 * equations.  Prior to it's call, the values for continuous states y,
 * discrete states z, excitation states u, and system parameters p must
 * be set to the newest values.  The function updates the dynamics (i.e.,
 * computes mass matrix, differentiable forces, etc. automatically.  It
 * returns an array (size NY) of the right hand side.
 */
void PB_Dynamics::FlowMap(double* dydtPtr){
	updateDynamics();

	VectorYd dydt;

	/**********************
	 * Actuator dynamics:
	 **********************/
	// Map velocities to position derivatives:
	dydt(ualphaL)  = yVec(dualphaL);
	dydt(ualphaR)  = yVec(dualphaR);
	dydt(ulL)      = yVec(dulL);
	dydt(ulR)      = yVec(dulR);
	if (pVec(hip_jnt_type) == PEA){ // No actuator dynamics in the parallel actuated joints
		dydt(dualphaL) = 0;
		dydt(dualphaR) = 0;
	} else { // Actuators are driven by the difference between motor and joint force
		dydt(dualphaL) = (uVec(TalphaL) - tau(qalphaL))/j_mot_alpha;
		dydt(dualphaR) = (uVec(TalphaR) - tau(qalphaR))/j_mot_alpha;
	}
	if (pVec(leg_jnt_type) == PEA){ // No actuator dynamics in the parallel actuated joints
		dydt(dulL)     = 0;
		dydt(dulR)     = 0;
	} else { // Actuators are driven by the difference between motor and joint force
		dydt(dulL)     = (uVec(FlL) - tau(qlL))/j_mot_l;
		dydt(dulR)     = (uVec(FlR) - tau(qlR))/j_mot_l;
	}

	/**********************
	 * Mechanical dynamics:
	 **********************/
	// Map velocities to position derivatives:
	dydt(x)      = yVec(dx);
	dydt(y)      = yVec(dy);
	dydt(phi)    = yVec(dphi);
	dydt(alphaL) = yVec(dalphaL);
	dydt(alphaR) = yVec(dalphaR);
	dydt(lL)     = yVec(dlL);
	dydt(lR)     = yVec(dlR);

	// Compute contact forces:
	VectorQd Jlambda;
	if (nContact>0) // check if some ground contact is present
	{
		MatrixXd projM(nContact*2, nContact*2);
		VectorXd projF(nContact*2, 1);
		VectorXd lambda(nContact*2, 1);
		projM = J * Mdecomp.solve(J.transpose());
		projF = J * Mdecomp.solve(tau + h);
		lambda = projM.partialPivLu().solve(-projF - dJdtTimesdqdt);
		// Project these forces back into the generalized coordinate space
		Jlambda = J.transpose()*lambda;
	} else {
		Jlambda = VectorQd::Zero();
	}

	// Evaluate the EQM:
	VectorQd dd_q;
	dd_q = Mdecomp.solve(tau + h + Jlambda);

	// Map the generalized accelerations back into continuous state derivatives:
	dydt(dx)      = dd_q(qx);
	dydt(dy)      = dd_q(qy);
	dydt(dphi)    = dd_q(qphi);
	dydt(dalphaL) = dd_q(qalphaL);
	dydt(dalphaR) = dd_q(qalphaR);
	dydt(dlL)     = dd_q(qlL);
	dydt(dlR)     = dd_q(qlR);

	/*****************
	 * Cost Functions:
	 *****************/
	// Compute the mechanical power performed by the actuator (that is actuator velocity times actuator torque/force):
	double P_act_alphaL, P_act_alphaR, P_act_lL, P_act_lR;
	if (pVec(hip_jnt_type) == PEA){ // Actuator velocity is joint velocity
		P_act_alphaL = yVec(dalphaL)*uVec(TalphaL);
		P_act_alphaR = yVec(dalphaR)*uVec(TalphaR);
	} else { // Actuator velocity is different from joint velocity
		P_act_alphaL = yVec(dualphaL)*uVec(TalphaL);
		P_act_alphaR = yVec(dualphaR)*uVec(TalphaR);
	}
	if (pVec(leg_jnt_type) == PEA){ // Actuator velocity is joint velocity
		P_act_lL = yVec(dlL)*uVec(FlL);
		P_act_lR = yVec(dlR)*uVec(FlR);
	} else { // Actuator velocity is different from joint velocity
		P_act_lL = yVec(dulL)*uVec(FlL);
		P_act_lR = yVec(dulR)*uVec(FlR);
	}
	// Compute the electrical losses:
	double P_loss_el_alphaL = uVec(TalphaL)*uVec(TalphaL) * pVec(du_max_alpha)*pVec(du_max_alpha) / P_max_alpha;
	double P_loss_el_alphaR = uVec(TalphaR)*uVec(TalphaR) * pVec(du_max_alpha)*pVec(du_max_alpha) / P_max_alpha;
	double P_loss_el_lL = uVec(FlL)*uVec(FlL) * pVec(du_max_l)*pVec(du_max_l) / P_max_l;
	double P_loss_el_lR = uVec(FlR)*uVec(FlR) * pVec(du_max_l)*pVec(du_max_l) / P_max_l;

	// Compute the accumulated positive actuator work (assuming that negative work at EACH actuator can not be recovered):
	dydt(posActWork)  = max(0.0, P_act_alphaL) +
                        max(0.0, P_act_alphaR) +
                        max(0.0, P_act_lL) +
                        max(0.0, P_act_lR);
    // Compute the total electrical Losses
	dydt(totElLoss)  = P_loss_el_alphaL + P_loss_el_alphaR + P_loss_el_lL + P_loss_el_lR;
	// Compute the positive electrical work (which allows using negative work in one actuator to perform positive work in another):
    dydt(posElWork)  = max(0.0, P_act_alphaL + P_act_alphaR + P_act_lL + P_act_lR + double(dydt(totElLoss)));


	/*******************
	 * Write to pointer:
	 *******************/
    for (int i = 0; i<NY; i++)
    	dydtPtr[i] = dydt(i);
}
/* The JumpMap computes the changes in velocities, that are necessary to
 * bring the contact points to a velocity of zero.  Prior to it's call,
 * the values for continuous states y, discrete states z, excitation
 * states u, and system parameters p must be set to the newest values.
 * **NOTE** the values of z need to represent the phases AFTER the
 * collision.  The function updates the dynamics (i.e., computes mass
 * matrix, differentiable forces, etc. automatically.  It returns an
 * array (size NY) of the state after the collision
 */
void PB_Dynamics::JumpMap(double* yPlusPtr){
	updateDynamics();

	// As most things remain unchanged we copy the current states
	for (int i = 0; i<NY; i++)
		yPlusPtr[i] = yVec(i);

	MatrixXd projM(nContact*2, nContact*2);
	MatrixQd projMinv;
	VectorQd dqMINUS;
	VectorQd dqPLUS;
	// At touchdown, the contact point comes to a complete rest, a fully plastic collision is computed:
	// Velocities before collision:
	dqMINUS(qx)      = yVec(dx);
	dqMINUS(qy)      = yVec(dy);
	dqMINUS(qphi)    = yVec(dphi);
	dqMINUS(qalphaL) = yVec(dalphaL);
	dqMINUS(qalphaR) = yVec(dalphaR);
	dqMINUS(qlL)     = yVec(dlL);
	dqMINUS(qlR)     = yVec(dlR);
	// Project EoM into the contact space:
	// Rest after contact: J*qPlus = 0
	// with: M*(qPlus-qMinus) = I_cont_q % with the contact impulse I_cont
	// qPlus = inv(M)*I_cont_q + qMinus
	// -> J*inv(M)*I_cont_q + qMinus = 0
	// -> J*inv(M)*J'*I_cont_x + qMinus = 0
	// -> I_cont_x = -inv(J*inv(M)*J')*qMinus
	// -> I_cont_q = -J'*inv(J*inv(M)*J')*J*qMinus
	// qPlus = -inv(M)*J'*inv(J*inv(M)*J')*J*qMinus + qMinus
	projM = J * Mdecomp.solve(J.transpose());
	projMinv = J.transpose() * projM.partialPivLu().solve(J);
	// Velocities after
	dqPLUS = Mdecomp.solve(M*dqMINUS - projMinv*dqMINUS);
	yPlusPtr[dx]       = dqPLUS(qx);
	yPlusPtr[dy]       = dqPLUS(qy);
	yPlusPtr[dphi]     = dqPLUS(qphi);
	yPlusPtr[dalphaL]  = dqPLUS(qalphaL);
	yPlusPtr[dalphaR]  = dqPLUS(qalphaR);
	yPlusPtr[dlL]      = dqPLUS(qlL);
	yPlusPtr[dlR]      = dqPLUS(qlR);
}
/* The JumpSet tracks the various events that can happen throughout a
 * stride. Each event value will have a zero crossing in positive
 * direction, when the corresponding event happens  Prior to it's call,
 * the values for continuous states y, discrete states z, excitation
 * states u, and system parameters p must be set to the newest values.
 * **NOTE** the values of z need to represent the phases BEFORE the
 * event.  If one seeks to identify multiple events simultaneously,
 * all involved event values must be zero.  They are returned in an
 * array (size NEV).
 */
void PB_Dynamics::JumpSet(double* eventValPtr){
	updateDynamics();

	// Compute the contact forces:
	VectorXd lambda(nContact*2, 1);
	if (nContact>0)
	{  // check if some ground contact is present
		MatrixXd projM(nContact*2, nContact*2);
		VectorXd projF(nContact*2, 1);
		projM = J * Mdecomp.solve(J.transpose());
		projF = J * Mdecomp.solve(tau + h);
		lambda = projM.partialPivLu().solve(-projF - dJdtTimesdqdt);
	} else {
		// Do nothing, lambda is an empty vector anyway
	}
	int counter = 0;
	if (zVec(phaseL) == stance) {	//  Event is detected if the vertical contact force of the left leg becomes negative.
		eventValPtr[liftoffL] = -lambda (counter+1);
		counter = counter + 2;
	} else { // But only in flight
		eventValPtr[liftoffL] = -1;
	}
	if (zVec(phaseR) == stance) {	//  Event is detected if the vertical contact force of the right leg becomes negative.
		eventValPtr[liftoffR] = -lambda (counter+1);
		counter = counter + 2;
	} else { // But only in flight
		eventValPtr[liftoffR] = -1;
	}
	if (zVec(phaseL) == flight) {	// Event is detected if the left foot goes below the ground during flight
		eventValPtr[touchdownL] = -posL(1);
	} else { // But only in flight
		eventValPtr[touchdownL] = -1;
	}
	if (zVec(phaseR) == flight) {	// Event is detected if the right foot goes below the ground during flight
		eventValPtr[touchdownR] = -posR(1);
	} else { // But only in flight
		eventValPtr[touchdownR] = -1;
	}
	// Detect apex transit
	eventValPtr[apexTransit] = yVec(dy);
}

/* Compute a number of constraint values for
 * a) The position of the left and right foot above the ground
 * b) The velocity du of each actuator below du_max*(1-c_lim)
 * c) The force Fu of each actuator below F_max = c_lim*P_max/du_max
 * I.e., each constraint is met as soon it is larger than 0.
 * the vectors for the actuator speeds and velocities are given in the order:
 * [alphaL, alphaR, lL, lR]
 */
void PB_Dynamics::Constraints(double* gcL, double* gcR, double* du_max, double* F_max){
	updateDynamics();

	// Ground clearance:
	if (zVec(phaseL) == flight) {
		*gcL = posL(1);
	} else {
		*gcL = 1; // if the foot is on the ground, the constraint is always met
	}

	if (zVec(phaseR) == flight) {
		*gcR = posR(1);
	} else {
		*gcR = 1; // If the foot is on the ground, the constraint is always met
	}

	if (pVec(hip_jnt_type) == PEA){ // Actuator velocity is joint velocity
		// Compute how much the actual joint velocity is below its upper limit:
		du_max[0] = pVec(du_max_alpha)*(1-c_lim_alpha)-yVec(dalphaL);
		du_max[1] = pVec(du_max_alpha)*(1-c_lim_alpha)-yVec(dalphaR);
		// Compute how much the actual joint velocity is above its lower limit:
		du_max[4] = pVec(du_max_alpha)*(1-c_lim_alpha)+yVec(dalphaL);
		du_max[5] = pVec(du_max_alpha)*(1-c_lim_alpha)+yVec(dalphaR);
	} else { // Actuator velocity is different from joint velocity
		// Compute how much the actual actuator velocity is below its upper limit:
		du_max[0] = pVec(du_max_alpha)*(1-c_lim_alpha)-yVec(dualphaL);
		du_max[1] = pVec(du_max_alpha)*(1-c_lim_alpha)-yVec(dualphaR);
		// Compute how much the actual actuator velocity is above its lower limit:
		du_max[4] = pVec(du_max_alpha)*(1-c_lim_alpha)+yVec(dualphaL);
		du_max[5] = pVec(du_max_alpha)*(1-c_lim_alpha)+yVec(dualphaR);
	}
	if (pVec(leg_jnt_type) == PEA){ // Actuator velocity is joint velocity
		// Compute how much the actual joint velocity is below its upper limit:
		du_max[2] = pVec(du_max_l)*(1-c_lim_l)-yVec(dlL);
		du_max[3] = pVec(du_max_l)*(1-c_lim_l)-yVec(dlR);
		// Compute how much the actual joint velocity is above its lower limit:
		du_max[6] = pVec(du_max_l)*(1-c_lim_l)+yVec(dlL);
		du_max[7] = pVec(du_max_l)*(1-c_lim_l)+yVec(dlR);
	} else { // Actuator velocity is different from joint velocity
		// Compute how much the actual actuator velocity is below its upper limit:
		du_max[2] = pVec(du_max_l)*(1-c_lim_l)-yVec(dulL);
		du_max[3] = pVec(du_max_l)*(1-c_lim_l)-yVec(dulR);
		// Compute how much the actual actuator velocity is above its lower limit:
		du_max[6] = pVec(du_max_l)*(1-c_lim_l)+yVec(dulL);
		du_max[7] = pVec(du_max_l)*(1-c_lim_l)+yVec(dulR);
	}

	// Compute how much the actual actuator torque is below its upper limit:
	F_max[0] = c_lim_alpha*P_max_alpha/pVec(du_max_alpha) - uVec(TalphaL);
	F_max[1] = c_lim_alpha*P_max_alpha/pVec(du_max_alpha) - uVec(TalphaR);
	F_max[2] = c_lim_l*P_max_l/pVec(du_max_l) - uVec(FlL);
	F_max[3] = c_lim_l*P_max_l/pVec(du_max_l) - uVec(FlR);
	// Compute how much the actual actuator torque is above its lower limit:
	F_max[4] = c_lim_alpha*P_max_alpha/pVec(du_max_alpha) + uVec(TalphaL);
	F_max[5] = c_lim_alpha*P_max_alpha/pVec(du_max_alpha) + uVec(TalphaR);
	F_max[6] = c_lim_l*P_max_l/pVec(du_max_l) + uVec(FlL);
	F_max[7] = c_lim_l*P_max_l/pVec(du_max_l) + uVec(FlR);
}

// Update current states
void PB_Dynamics::setContState(double* yPtr){yVec = Map<VectorYd>(yPtr);}
void PB_Dynamics::setDiscState(double* zPtr){zVec = Map<VectorZd>(zPtr);}
void PB_Dynamics::setSystParam(double* pPtr){pVec = Map<VectorPd>(pPtr);}
void PB_Dynamics::setExctState(double* uPtr){uVec = Map<VectorUd>(uPtr);}
// Retrieve the current states
void PB_Dynamics::getContState(double* yPtr){
	for (int i = 0; i<NY; i++)
		yPtr[i] = yVec(i);
}
void PB_Dynamics::getDiscState(double* zPtr){
	for (int i = 0; i<NZ; i++)
		zPtr[i] = zVec(i);
}
void PB_Dynamics::getSystParam(double* pPtr){
	for (int i = 0; i<NP; i++)
		pPtr[i] = pVec(i);
}
void PB_Dynamics::getExctState(double* uPtr){
	for (int i = 0; i<NU; i++)
		uPtr[i] = uVec(i);
}

// Prepare all the components for the equations of motion:
void PB_Dynamics::updateDynamics(){
	// Before this function is called, we assume that all states and parameters have been updated.
	// Update all parameters (inertia values, damping coefficients, ...)
	ComputeDependentParameters();
	// We compute new values for all the 'dynamics'-components:
	ComputeMassMatrix();  // M & Mdecomp
	ComputeDiffForces();  // h
	ComputeJointForces(); // f (this also computes T_spring_alpha and F_spring_l)
	nContact = 0;
	// And for the contact processing:
	if (zVec(phaseL) == flight) {
		ComputeContactPointL();     // posL
	} else {
		nContact++;
		ComputeContactJacobianL();  //JL
		ComputeContactJacobianDtTIMESdqdtL();  // dJLdtTimesdqdt
	}
	if (zVec(phaseR) == flight) {
		ComputeContactPointR();     // posL
	} else {
		nContact++;
		ComputeContactJacobianR();  //JR
		ComputeContactJacobianDtTIMESdqdtR();  // dJRdtTimesdqdt
	}
	//cout << "updateDynamics"  << nContact << endl;
	// And build the compound Jacobian (which starts with L, then R)
	J.resize(nContact*2, NQ);
	dJdtTimesdqdt.resize(nContact*2, 1);
	int counter = 0;
	if (zVec(phaseL) == stance) {
		J.row(counter)   = JL.row(0);
		J.row(counter+1) = JL.row(1);
		dJdtTimesdqdt(counter)   = dJLdtTimesdqdt(0);
		dJdtTimesdqdt(counter+1) = dJLdtTimesdqdt(1);
		counter = counter + 2;
	}
	if (zVec(phaseR) == stance) {
		J.row(counter)   = JR.row(0);
		J.row(counter+1) = JR.row(1);
		dJdtTimesdqdt(counter)   = dJRdtTimesdqdt(0);
		dJdtTimesdqdt(counter+1) = dJRdtTimesdqdt(1);
		counter = counter + 2;
	}
}

// Update dependent parameters
void PB_Dynamics::ComputeDependentParameters(){
	// Compute maximal power and the unscaled actuator inertia according to the scaling laws and
	// the values provided for rho_alpha and rho_l:
	P_max_alpha  = pVec(P_max)  * pow(double(pVec(rho_alpha)), +1.35);
	P_max_l      = pVec(P_max)  * pow(double(pVec(rho_l)),     +1.35);
	// Percentage of the maximal power that can be used
	c_lim_alpha  = pVec(c_lim)  * pow(double(pVec(rho_alpha)), -0.37);
	c_lim_l      = pVec(c_lim)  * pow(double(pVec(rho_l)),     -0.37);
	// Compute the actual actuator inertias as a function of the gear ratios (maximal speeds)
	double j_unsc_alpha = pVec(j_unsc) * pow(double(pVec(rho_alpha)), +1.26);
	double j_unsc_l     = pVec(j_unsc) * pow(double(pVec(rho_l)),     +1.26);
	j_mot_alpha = j_unsc_alpha/(pVec(du_max_alpha)*pVec(du_max_alpha));
	j_mot_l     = j_unsc_l/(pVec(du_max_l)*pVec(du_max_l));
	// Compute the viscous damping coefficient of the springs, according to the desired damping ratio:
	double j_leg   = pVec(j3) + pow(double(pVec(l_0) - pVec(l3)),2)*pVec(m3) + pVec(j2) + pow(double(pVec(l2)),2)*pVec(m2); // total leg inertia wrt the hip
	b_alpha  = pVec(balphaRat)*2*sqrt(double(pVec(kalpha)*j_leg));
	b_l      = pVec(blRat)*2*sqrt(double(pVec(kl)*pVec(m3)));
	// To prevent instability damping values are forced to be non-negative:
	b_alpha = max(0.0, b_alpha);
	b_l     = max(0.0, b_l);
}

// Components for the Equations of Motion
void PB_Dynamics::ComputeMassMatrix(){
	M = MatrixQd::Zero();
	double t293 = yVec(alphaL)+yVec(phi);
	double t294 = cos(t293);
	double t295 = yVec(alphaR)+yVec(phi);
	double t296 = cos(t295);
	double t297 = pVec(l2)*pVec(m2)*t294;
	double t298 = yVec(lL)*pVec(m3)*t294;
	double t299 = pVec(l2)*pVec(m2)*t296;
	double t300 = yVec(lR)*pVec(m3)*t296;
	double t301 = pVec(m2)*2.0;
	double t302 = pVec(m3)*2.0;
	double t303 = pVec(m1)+t301+t302;
	double t304 = sin(t293);
	double t305 = sin(t295);
	double t306 = pVec(l2)*pVec(m2)*t304;
	double t307 = yVec(lL)*pVec(m3)*t304;
	double t308 = pVec(l2)*pVec(m2)*t305;
	double t309 = yVec(lR)*pVec(m3)*t305;
	double t320 = pVec(l3)*pVec(m3)*t294;
	double t325 = pVec(l3)*pVec(m3)*t296;
	double t310 = t297+t298+t299+t300-t320-t325;
	double t321 = pVec(l3)*pVec(m3)*t304;
	double t326 = pVec(l3)*pVec(m3)*t305;
	double t311 = t306+t307+t308+t309-t321-t326;
	double t312 = pVec(l2)*pVec(l2);
	double t313 = pVec(l3)*pVec(l3);
	double t314 = yVec(lL)*yVec(lL);
	double t315 = pVec(m3)*t314;
	double t316 = pVec(m2)*t312;
	double t317 = pVec(m3)*t313;
	double t318 = yVec(lR)*yVec(lR);
	double t319 = pVec(m3)*t318;
	double t323 = pVec(l3)*yVec(lL)*pVec(m3)*2.0;
	double t322 = pVec(j2)+pVec(j3)+t315+t316+t317-t323;
	double t324 = pVec(m3)*t304;
	double t328 = pVec(l3)*yVec(lR)*pVec(m3)*2.0;
	double t327 = pVec(j2)+pVec(j3)+t316+t317+t319-t328;
	double t329 = pVec(m3)*t305;
	M(0,0) = t303;
	M(0,2) = t310;
	M(0,3) = t297+t298-pVec(l3)*pVec(m3)*t294;
	M(0,4) = t324;
	M(0,5) = t299+t300-pVec(l3)*pVec(m3)*t296;
	M(0,6) = t329;
	M(1,1) = t303;
	M(1,2) = t311;
	M(1,3) = t306+t307-pVec(l3)*pVec(m3)*t304;
	M(1,4) = -pVec(m3)*t294;
	M(1,5) = t308+t309-pVec(l3)*pVec(m3)*t305;
	M(1,6) = -pVec(m3)*t296;
	M(2,0) = t310;
	M(2,1) = t311;
	M(2,2) = pVec(j1_)+pVec(j2)*2.0+pVec(j3)*2.0+t315+t319+pVec(m2)*t312*2.0+pVec(m3)*t313*2.0-pVec(l3)*yVec(lL)*pVec(m3)*2.0-pVec(l3)*yVec(lR)*pVec(m3)*2.0;
	M(2,3) = t322;
	M(2,5) = t327;
	M(3,0) = t297+t298-t320;
	M(3,1) = t306+t307-t321;
	M(3,2) = t322;
	M(3,3) = t322;
	M(4,0) = t324;
	M(4,1) = -pVec(m3)*t294;
	M(4,4) = pVec(m3);
	M(5,0) = t299+t300-t325;
	M(5,1) = t308+t309-t326;
	M(5,2) = t327;
	M(5,5) = t327;
	M(6,0) = t329;
	M(6,1) = -pVec(m3)*t296;
	M(6,6) = pVec(m3);
	// Increase inertia matrix values for DOFs with parallel elastic actuation:
	// HIP
	if (pVec(hip_jnt_type) == PEA){
		M(qalphaL,qalphaL) = M(qalphaL,qalphaL) + j_mot_alpha;
		M(qalphaR,qalphaR) = M(qalphaR,qalphaR) + j_mot_alpha;
	}
	// LEG
	if (pVec(leg_jnt_type) == PEA){
		M(qlL,qlL) = M(qlL,qlL) + j_mot_l;
		M(qlR,qlR) = M(qlR,qlR) + j_mot_l;
	}
	// directly decompose the Mass-matrix:
	Mdecomp = M.ldlt();
}
void PB_Dynamics::ComputeDiffForces(){
	h = VectorQd::Zero();
	double t331 = yVec(alphaL)+yVec(phi);
	double t332 = sin(t331);
	double t333 = yVec(alphaR)+yVec(phi);
	double t334 = sin(t333);
	double t335 = cos(t331);
	double t336 = cos(t333);
	double t337 = yVec(dalphaL)*pVec(l2)*pVec(m2)*t332;
	double t338 = yVec(dalphaR)*pVec(l2)*pVec(m2)*t334;
	double t339 = yVec(dphi)*pVec(l2)*pVec(m2)*t332;
	double t340 = yVec(dphi)*pVec(l2)*pVec(m2)*t334;
	double t341 = yVec(dalphaL)*yVec(lL)*pVec(m3)*t332;
	double t342 = yVec(dalphaR)*yVec(lR)*pVec(m3)*t334;
	double t343 = yVec(dphi)*yVec(lL)*pVec(m3)*t332;
	double t344 = yVec(dphi)*yVec(lR)*pVec(m3)*t334;
	double t345 = yVec(dlL)*pVec(m3)*t332;
	double t346 = yVec(dlR)*pVec(m3)*t334;
	double t347 = yVec(dalphaL)*pVec(l2)*pVec(m2)*t335;
	double t348 = yVec(dalphaR)*pVec(l2)*pVec(m2)*t336;
	double t349 = yVec(dphi)*pVec(l2)*pVec(m2)*t335;
	double t350 = yVec(dphi)*pVec(l2)*pVec(m2)*t336;
	double t351 = yVec(dalphaL)*yVec(lL)*pVec(m3)*t335;
	double t352 = yVec(dalphaR)*yVec(lR)*pVec(m3)*t336;
	double t353 = yVec(dphi)*yVec(lL)*pVec(m3)*t335;
	double t354 = yVec(dphi)*yVec(lR)*pVec(m3)*t336;
	double t355 = yVec(dy)*pVec(l3)*pVec(m3)*t335;
	double t356 = yVec(dy)*pVec(l3)*pVec(m3)*t336;
	double t357 = yVec(dx)*pVec(l2)*pVec(m2)*t332;
	double t358 = yVec(dx)*pVec(l2)*pVec(m2)*t334;
	double t359 = yVec(dx)*yVec(lL)*pVec(m3)*t332;
	double t360 = yVec(dx)*yVec(lR)*pVec(m3)*t334;
	double t361 = yVec(dalphaL)*yVec(lL)*pVec(m3)*2.0;
	double t362 = yVec(dphi)*yVec(lL)*pVec(m3)*2.0;
	double t363 = yVec(dx)*pVec(m3)*t335;
	double t364 = yVec(dy)*pVec(m3)*t332;
	double t383 = yVec(dphi)*pVec(l3)*pVec(m3)*2.0;
	double t365 = t361+t362+t363+t364-t383-yVec(dalphaL)*pVec(l3)*pVec(m3)*2.0;
	double t367 = yVec(dy)*pVec(l2)*pVec(m2)*t335;
	double t368 = yVec(dy)*yVec(lL)*pVec(m3)*t335;
	double t369 = yVec(dx)*pVec(l3)*pVec(m3)*t332;
	double t366 = yVec(dalphaL)*(t355+t357+t359-t367-t368-t369);
	double t370 = pVec(l2)*pVec(m2)*t332;
	double t371 = yVec(lL)*pVec(m3)*t332;
	double t372 = yVec(dlL)*yVec(dx)*pVec(m3)*t335;
	double t373 = yVec(dlL)*yVec(dy)*pVec(m3)*t332;
	double t374 = yVec(dalphaL)*yVec(dy)*pVec(l2)*pVec(m2)*t335;
	double t375 = yVec(dphi)*yVec(dy)*pVec(l2)*pVec(m2)*t335;
	double t376 = yVec(dalphaL)*yVec(dy)*yVec(lL)*pVec(m3)*t335;
	double t377 = yVec(dphi)*yVec(dy)*yVec(lL)*pVec(m3)*t335;
	double t378 = yVec(dalphaL)*yVec(dx)*pVec(l3)*pVec(m3)*t332;
	double t379 = yVec(dphi)*yVec(dx)*pVec(l3)*pVec(m3)*t332;
	double t380 = t363+t364;
	double t381 = yVec(dalphaL)*yVec(dalphaL);
	double t382 = yVec(dphi)*yVec(dphi);
	double t384 = yVec(dalphaR)*yVec(lR)*pVec(m3)*2.0;
	double t385 = yVec(dphi)*yVec(lR)*pVec(m3)*2.0;
	double t386 = yVec(dx)*pVec(m3)*t336;
	double t387 = yVec(dy)*pVec(m3)*t334;
	double t389 = yVec(dy)*pVec(l2)*pVec(m2)*t336;
	double t390 = yVec(dy)*yVec(lR)*pVec(m3)*t336;
	double t391 = yVec(dx)*pVec(l3)*pVec(m3)*t334;
	double t388 = yVec(dalphaR)*(t356+t358+t360-t389-t390-t391);
	double t392 = pVec(l2)*pVec(m2)*t334;
	double t393 = yVec(lR)*pVec(m3)*t334;
	double t394 = yVec(dlR)*yVec(dx)*pVec(m3)*t336;
	double t395 = yVec(dlR)*yVec(dy)*pVec(m3)*t334;
	double t396 = yVec(dalphaR)*yVec(dy)*pVec(l2)*pVec(m2)*t336;
	double t397 = yVec(dphi)*yVec(dy)*pVec(l2)*pVec(m2)*t336;
	double t398 = yVec(dalphaR)*yVec(dy)*yVec(lR)*pVec(m3)*t336;
	double t399 = yVec(dphi)*yVec(dy)*yVec(lR)*pVec(m3)*t336;
	double t400 = yVec(dalphaR)*yVec(dx)*pVec(l3)*pVec(m3)*t334;
	double t401 = yVec(dphi)*yVec(dx)*pVec(l3)*pVec(m3)*t334;
	double t402 = t386+t387;
	double t403 = yVec(dalphaR)*yVec(dalphaR);
	h(0,0) = yVec(dphi)*(t337+t338+t339+t340+t341+t342+t343+t344-yVec(dlL)*pVec(m3)*t335-yVec(dlR)*pVec(m3)*t336-yVec(dalphaL)*pVec(l3)*pVec(m3)*t332-yVec(dalphaR)*pVec(l3)*pVec(m3)*t334-yVec(dphi)*pVec(l3)*pVec(m3)*t332-yVec(dphi)*pVec(l3)*pVec(m3)*t334)-yVec(dlL)*(yVec(dalphaL)*pVec(m3)*t335+yVec(dphi)*pVec(m3)*t335)-yVec(dlR)*(yVec(dalphaR)*pVec(m3)*t336+yVec(dphi)*pVec(m3)*t336)+yVec(dalphaL)*(t337+t339+t341+t343-yVec(dlL)*pVec(m3)*t335-yVec(dalphaL)*pVec(l3)*pVec(m3)*t332-yVec(dphi)*pVec(l3)*pVec(m3)*t332)+yVec(dalphaR)*(t338+t340+t342+t344-yVec(dlR)*pVec(m3)*t336-yVec(dalphaR)*pVec(l3)*pVec(m3)*t334-yVec(dphi)*pVec(l3)*pVec(m3)*t334);
	h(1,0) = -pVec(g)*(pVec(m1)+pVec(m2)*2.0+pVec(m3)*2.0)-yVec(dlL)*(yVec(dalphaL)*pVec(m3)*t332+yVec(dphi)*pVec(m3)*t332)-yVec(dlR)*(yVec(dalphaR)*pVec(m3)*t334+yVec(dphi)*pVec(m3)*t334)-yVec(dalphaL)*(t345+t347+t349+t351+t353-yVec(dalphaL)*pVec(l3)*pVec(m3)*t335-yVec(dphi)*pVec(l3)*pVec(m3)*t335)-yVec(dalphaR)*(t346+t348+t350+t352+t354-yVec(dalphaR)*pVec(l3)*pVec(m3)*t336-yVec(dphi)*pVec(l3)*pVec(m3)*t336)-yVec(dphi)*(t345+t346+t347+t348+t349+t350+t351+t352+t353+t354-yVec(dalphaL)*pVec(l3)*pVec(m3)*t335-yVec(dalphaR)*pVec(l3)*pVec(m3)*t336-yVec(dphi)*pVec(l3)*pVec(m3)*t335-yVec(dphi)*pVec(l3)*pVec(m3)*t336);
	h(2,0) = t366+t372+t373+t374+t375+t376+t377+t378+t379+t388+t394+t395+t396+t397+t398+t399+t400+t401-yVec(dlR)*(t384+t385+t386+t387-yVec(dalphaR)*pVec(l3)*pVec(m3)*2.0-yVec(dphi)*pVec(l3)*pVec(m3)*2.0)-yVec(dlL)*t365+yVec(dphi)*(t355+t356+t357+t358+t359+t360-yVec(dy)*pVec(l2)*pVec(m2)*t335-yVec(dy)*pVec(l2)*pVec(m2)*t336-yVec(dx)*pVec(l3)*pVec(m3)*t332-yVec(dx)*pVec(l3)*pVec(m3)*t334-yVec(dy)*yVec(lL)*pVec(m3)*t335-yVec(dy)*yVec(lR)*pVec(m3)*t336)-pVec(g)*(t370+t371+t392+t393-pVec(l3)*pVec(m3)*t332-pVec(l3)*pVec(m3)*t334)-yVec(dalphaL)*yVec(dx)*pVec(l2)*pVec(m2)*t332-yVec(dalphaR)*yVec(dx)*pVec(l2)*pVec(m2)*t334-yVec(dalphaL)*yVec(dy)*pVec(l3)*pVec(m3)*t335-yVec(dalphaR)*yVec(dy)*pVec(l3)*pVec(m3)*t336-yVec(dphi)*yVec(dx)*pVec(l2)*pVec(m2)*t332-yVec(dphi)*yVec(dx)*pVec(l2)*pVec(m2)*t334-yVec(dphi)*yVec(dy)*pVec(l3)*pVec(m3)*t335-yVec(dphi)*yVec(dy)*pVec(l3)*pVec(m3)*t336-yVec(dalphaL)*yVec(dx)*yVec(lL)*pVec(m3)*t332-yVec(dalphaR)*yVec(dx)*yVec(lR)*pVec(m3)*t334-yVec(dphi)*yVec(dx)*yVec(lL)*pVec(m3)*t332-yVec(dphi)*yVec(dx)*yVec(lR)*pVec(m3)*t334;
	h(3,0) = t366+t372+t373+t374+t375+t376+t377+t378+t379-yVec(dlL)*t365-pVec(g)*(t370+t371-pVec(l3)*pVec(m3)*t332)+yVec(dphi)*(t355+t357+t359-t367-t368-t369)-yVec(dalphaL)*yVec(dx)*pVec(l2)*pVec(m2)*t332-yVec(dalphaL)*yVec(dy)*pVec(l3)*pVec(m3)*t335-yVec(dphi)*yVec(dx)*pVec(l2)*pVec(m2)*t332-yVec(dphi)*yVec(dy)*pVec(l3)*pVec(m3)*t335-yVec(dalphaL)*yVec(dx)*yVec(lL)*pVec(m3)*t332-yVec(dphi)*yVec(dx)*yVec(lL)*pVec(m3)*t332;
	h(4,0) = -yVec(dalphaL)*t380-yVec(dphi)*t380+pVec(g)*pVec(m3)*t335-pVec(l3)*pVec(m3)*t381-pVec(l3)*pVec(m3)*t382+yVec(lL)*pVec(m3)*t381+yVec(lL)*pVec(m3)*t382-yVec(dalphaL)*yVec(dphi)*pVec(l3)*pVec(m3)*2.0+yVec(dalphaL)*yVec(dphi)*yVec(lL)*pVec(m3)*2.0+yVec(dalphaL)*yVec(dx)*pVec(m3)*t335+yVec(dalphaL)*yVec(dy)*pVec(m3)*t332+yVec(dphi)*yVec(dx)*pVec(m3)*t335+yVec(dphi)*yVec(dy)*pVec(m3)*t332;
	h(5,0) = t388+t394+t395+t396+t397+t398+t399+t400+t401-yVec(dlR)*(-t383+t384+t385+t386+t387-yVec(dalphaR)*pVec(l3)*pVec(m3)*2.0)-pVec(g)*(t392+t393-pVec(l3)*pVec(m3)*t334)+yVec(dphi)*(t356+t358+t360-t389-t390-t391)-yVec(dalphaR)*yVec(dx)*pVec(l2)*pVec(m2)*t334-yVec(dalphaR)*yVec(dy)*pVec(l3)*pVec(m3)*t336-yVec(dphi)*yVec(dx)*pVec(l2)*pVec(m2)*t334-yVec(dphi)*yVec(dy)*pVec(l3)*pVec(m3)*t336-yVec(dalphaR)*yVec(dx)*yVec(lR)*pVec(m3)*t334-yVec(dphi)*yVec(dx)*yVec(lR)*pVec(m3)*t334;
	h(6,0) = -yVec(dalphaR)*t402-yVec(dphi)*t402+pVec(g)*pVec(m3)*t336-pVec(l3)*pVec(m3)*t382-pVec(l3)*pVec(m3)*t403+yVec(lR)*pVec(m3)*t382+yVec(lR)*pVec(m3)*t403-yVec(dalphaR)*yVec(dphi)*pVec(l3)*pVec(m3)*2.0+yVec(dalphaR)*yVec(dphi)*yVec(lR)*pVec(m3)*2.0+yVec(dalphaR)*yVec(dx)*pVec(m3)*t336+yVec(dalphaR)*yVec(dy)*pVec(m3)*t334+yVec(dphi)*yVec(dx)*pVec(m3)*t336+yVec(dphi)*yVec(dy)*pVec(m3)*t334;}
void  PB_Dynamics::ComputeContactPointL(){
	posL =  Vector2d::Zero();
	double t424 = yVec(alphaL)+yVec(phi);
	posL(0,0) = yVec(x)+yVec(lL)*sin(t424);
	posL(1,0) = -pVec(rFoot)+yVec(y)-yVec(lL)*cos(t424);
}
void  PB_Dynamics::ComputeContactPointR(){
	posR =  Vector2d::Zero();
	double t405 = yVec(alphaR)+yVec(phi);
	posR(0,0) = yVec(x)+yVec(lR)*sin(t405);
	posR(1,0) = -pVec(rFoot)+yVec(y)-yVec(lR)*cos(t405);
}
void PB_Dynamics::ComputeContactJacobianL(){
	JL =  Matrix2Qd::Zero();
	double t426 = yVec(alphaL)+yVec(phi);
	double t427 = cos(t426);
	double t428 = yVec(lL)*t427;
	double t429 = pVec(rFoot)+t428;
	double t430 = sin(t426);
	double t431 = yVec(lL)*t430;
	JL(0,0) = 1.0;
	JL(0,2) = t429;
	JL(0,3) = t429;
	JL(0,4) = t430;
	JL(1,1) = 1.0;
	JL(1,2) = t431;
	JL(1,3) = t431;
	JL(1,4) = -t427;
}
void PB_Dynamics::ComputeContactJacobianR(){
	JR =  Matrix2Qd::Zero();
	double t407 = yVec(alphaR)+yVec(phi);
	double t408 = cos(t407);
	double t409 = yVec(lR)*t408;
	double t410 = sin(t407);
	double t411 = yVec(lR)*t410;
	JR(0,0) = 1.0;
	JR(0,2) = pVec(rFoot)+t409;
	JR(0,5) = pVec(rFoot)+t409;
	JR(0,6) = t410;
	JR(1,1) = 1.0;
	JR(1,2) = t411;
	JR(1,5) = t411;
	JR(1,6) = -t408;
}
void  PB_Dynamics::ComputeContactJacobianDtTIMESdqdtL(){
	dJLdtTimesdqdt =  Vector2d::Zero();
	double t433 = yVec(alphaL)+yVec(phi);
	double t434 = sin(t433);
	double t435 = cos(t433);
	double t436 = yVec(dalphaL)*yVec(lL)*t434;
	double t437 = yVec(dphi)*yVec(lL)*t434;
	double t438 = t436+t437-yVec(dlL)*t435;
	double t439 = yVec(dlL)*t434;
	double t440 = yVec(dalphaL)*yVec(lL)*t435;
	double t441 = yVec(dphi)*yVec(lL)*t435;
	double t442 = t439+t440+t441;
	dJLdtTimesdqdt(0,0) = -yVec(dalphaL)*t438-yVec(dphi)*t438+yVec(dlL)*(yVec(dalphaL)*t435+yVec(dphi)*t435);
	dJLdtTimesdqdt(1,0) = yVec(dalphaL)*t442+yVec(dphi)*t442+yVec(dlL)*(yVec(dalphaL)*t434+yVec(dphi)*t434);
}
void  PB_Dynamics::ComputeContactJacobianDtTIMESdqdtR(){
	dJRdtTimesdqdt =  Vector2d::Zero();
	double t413 = yVec(alphaR)+yVec(phi);
	double t414 = sin(t413);
	double t415 = cos(t413);
	double t416 = yVec(dalphaR)*yVec(lR)*t414;
	double t417 = yVec(dphi)*yVec(lR)*t414;
	double t418 = t416+t417-yVec(dlR)*t415;
	double t419 = yVec(dlR)*t414;
	double t420 = yVec(dalphaR)*yVec(lR)*t415;
	double t421 = yVec(dphi)*yVec(lR)*t415;
	double t422 = t419+t420+t421;
	dJRdtTimesdqdt(0,0) = -yVec(dalphaR)*t418-yVec(dphi)*t418+yVec(dlR)*(yVec(dalphaR)*t415+yVec(dphi)*t415);
	dJRdtTimesdqdt(1,0) = yVec(dalphaR)*t422+yVec(dphi)*t422+yVec(dlR)*(yVec(dalphaR)*t414+yVec(dphi)*t414);
}
void PB_Dynamics::ComputeJointForces(){
	// Joint forces:
	tau = VectorQd::Zero();
	// Compute spring and damping forces:
	double T_spring_alphaL, T_spring_alphaR,
	       F_spring_lL, F_spring_lR;
	// HIP
	if (pVec(hip_jnt_type) == PEA){ // Parallel spring in the hip
		T_spring_alphaL = pVec(kalpha)*(0 - yVec(alphaL)) +
				              b_alpha *(0 - yVec(dalphaL));
		T_spring_alphaR = pVec(kalpha)*(0 - yVec(alphaR)) +
				              b_alpha *(0 - yVec(dalphaR));
		tau(qalphaL) = T_spring_alphaL + uVec(TalphaL);
		tau(qalphaR) = T_spring_alphaR + uVec(TalphaR);
	} else { // Serial spring in the hip
		T_spring_alphaL = pVec(kalpha)*(yVec(ualphaL)  - yVec(alphaL)) +
				              b_alpha *(yVec(dualphaL) - yVec(dalphaL));
		T_spring_alphaR = pVec(kalpha)*(yVec(ualphaR)  - yVec(alphaR)) +
				              b_alpha *(yVec(dualphaR) - yVec(dalphaR));
		tau(qalphaL) = T_spring_alphaL + 0;
		tau(qalphaR) = T_spring_alphaR + 0;
	}
	// LEG
	if (pVec(leg_jnt_type) == PEA){ // Parallel spring in the leg
		F_spring_lL = pVec(kl)*(pVec(l_0) + 0 - yVec(lL)) +
			              b_l *(          + 0 - yVec(dlL));
		F_spring_lR = pVec(kl)*(pVec(l_0) + 0 - yVec(lR)) +
			              b_l *(          + 0 - yVec(dlR));
		tau(qlL) = F_spring_lL + uVec(FlL);
		tau(qlR) = F_spring_lR + uVec(FlR);
	} else { // Serial spring in the leg
		F_spring_lL = pVec(kl)*(pVec(l_0) + yVec(ulL)  - yVec(lL)) +
			              b_l *(          + yVec(dulL) - yVec(dlL));
		F_spring_lR = pVec(kl)*(pVec(l_0) + yVec(ulR)  - yVec(lR)) +
			              b_l *(          + yVec(dulR) - yVec(dlR));
		tau(qlL) = F_spring_lL + 0;
		tau(qlR) = F_spring_lR + 0;
	}
}


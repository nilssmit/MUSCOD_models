/*
 * AB_Dynamics.cpp
 *
 * Common Dynamic code for an Articulating Biped.
 *
 *  Created on: May 15, 2015
 *      Author: nilssmit@umich.edu
 */

#include "AB_Dynamics.h"
#include <algorithm>

AB_Dynamics::AB_Dynamics() {}  // Constructor (Empty)
AB_Dynamics::~AB_Dynamics() {} // Destructor (Empty)

/* Hybrid Dynamic Equations */
/* The FlowMap computes the right hand side of the governing differential
 * equations.  Prior to it's call, the values for continuous states y,
 * discrete states z, excitation states u, and system parameters p must
 * be set to the newest values.  The function updates the dynamics (i.e.,
 * computes mass matrix, differentiable forces, etc. automatically.  It
 * returns an array (size NY) of the right hand side.
 */
void AB_Dynamics::FlowMap(double* dydtPtr){
	updateDynamics();

	VectorYd dydt;

	/**********************
	 * Actuator dynamics:
	 **********************/
	// Map velocities to position derivatives:
	dydt(ualphaL)  = yVec(dualphaL);
	dydt(ualphaR)  = yVec(dualphaR);
	dydt(ubetaL)      = yVec(dubetaL);
	dydt(ubetaR)      = yVec(dubetaR);
	if (pVec(hip_jnt_type) == PEA){ // No actuator dynamics in the parallel actuated joints
		dydt(dualphaL) = 0;
		dydt(dualphaR) = 0;
	} else { // Actuators are driven by the difference between motor and joint force
		dydt(dualphaL) = (uVec(TalphaL) - tau(qalphaL))/pVec(j_mot_alpha);
		dydt(dualphaR) = (uVec(TalphaR) - tau(qalphaR))/pVec(j_mot_alpha);
	}
	if (pVec(knee_jnt_type) == PEA){ // No actuator dynamics in the parallel actuated joints
		dydt(dubetaL)     = 0;
		dydt(dubetaR)     = 0;
	} else { // Actuators are driven by the difference between motor and joint force
		dydt(dubetaL) = (uVec(TbetaL) - tau(qbetaL))/pVec(j_mot_beta);
		dydt(dubetaR) = (uVec(TbetaR) - tau(qbetaR))/pVec(j_mot_beta);
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
	dydt(betaL)     = yVec(dbetaL);
	dydt(betaR)     = yVec(dbetaR);

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
	dydt(dbetaL)     = dd_q(qbetaL);
	dydt(dbetaR)     = dd_q(qbetaR);

	/*****************
	 * Cost Functions:
	 *****************/
	// Compute the mechanical power performed by the actuator (that is actuator velocity times actuator torque/force):
	double P_act_alphaL, P_act_alphaR, P_act_betaL, P_act_betaR;
	if (pVec(hip_jnt_type) == PEA){ // Actuator velocity is joint velocity
		P_act_alphaL = yVec(dalphaL)*uVec(TalphaL);
		P_act_alphaR = yVec(dalphaR)*uVec(TalphaR);
	} else { // Actuator velocity is different from joint velocity
		P_act_alphaL = yVec(dualphaL)*uVec(TalphaL);
		P_act_alphaR = yVec(dualphaR)*uVec(TalphaR);
	}
	if (pVec(knee_jnt_type) == PEA){ // Actuator velocity is joint velocity
		P_act_betaL = yVec(dbetaL)*uVec(TbetaL);
		P_act_betaR = yVec(dbetaR)*uVec(TbetaR);
	} else { // Actuator velocity is different from joint velocity
		P_act_betaL = yVec(dubetaL)*uVec(TbetaL);
		P_act_betaR = yVec(dubetaR)*uVec(TbetaR);
	}
	// Compute the electrical losses:
	double P_loss_el_alphaL = uVec(TalphaL)*uVec(TalphaL) * pVec(du_max_alpha)*pVec(du_max_alpha) / pVec(P_max_alpha);
	double P_loss_el_alphaR = uVec(TalphaR)*uVec(TalphaR) * pVec(du_max_alpha)*pVec(du_max_alpha) / pVec(P_max_alpha);
	double P_loss_el_betaL = uVec(TbetaL)*uVec(TbetaL) * pVec(du_max_beta)*pVec(du_max_beta) / pVec(P_max_beta);
	double P_loss_el_betaR = uVec(TbetaR)*uVec(TbetaR) * pVec(du_max_beta)*pVec(du_max_beta) / pVec(P_max_beta);

	// Compute the accumulated positive actuator work (assuming that negative work at EACH actuator can not be recovered):
	dydt(posActWork)  = max(0.0, P_act_alphaL) +
                        max(0.0, P_act_alphaR) +
                        max(0.0, P_act_betaL) +
                        max(0.0, P_act_betaR);
    // Compute the total electrical Losses
	dydt(totElLoss)  = P_loss_el_alphaL + P_loss_el_alphaR + P_loss_el_betaL + P_loss_el_betaR;
	// Compute the positive electrical work (which allows using negative work in one actuator to perform positive work in another):
    dydt(posElWork)  = max(0.0, P_act_alphaL + P_act_alphaR + P_act_betaL + P_act_betaR + double(dydt(totElLoss)));


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
void AB_Dynamics::JumpMap(double* yPlusPtr){
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
	dqMINUS(qbetaL)     = yVec(dbetaL);
	dqMINUS(qbetaR)     = yVec(dbetaR);
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
	yPlusPtr[dbetaL]      = dqPLUS(qbetaL);
	yPlusPtr[dbetaR]      = dqPLUS(qbetaR);
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
void AB_Dynamics::JumpSet(double* eventValPtr){
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
 * [alphaL, alphaR, betaL, betaR]
 */
void AB_Dynamics::Constraints(double* gcL, double* gcR, double* du_max, double* F_max){
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
		du_max[0] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))-yVec(dalphaL);
		du_max[1] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))-yVec(dalphaR);
		// Compute how much the actual joint velocity is above its lower limit:
		du_max[4] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))+yVec(dalphaL);
		du_max[5] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))+yVec(dalphaR);
	} else { // Actuator velocity is different from joint velocity
		// Compute how much the actual actuator velocity is below its upper limit:
		du_max[0] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))-yVec(dualphaL);
		du_max[1] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))-yVec(dualphaR);
		// Compute how much the actual actuator velocity is above its lower limit:
		du_max[4] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))+yVec(dualphaL);
		du_max[5] = pVec(du_max_alpha)*(1-pVec(c_lim_alpha))+yVec(dualphaR);
	}
	if (pVec(knee_jnt_type) == PEA){ // Actuator velocity is joint velocity
		// Compute how much the actual joint velocity is below its upper limit:
		du_max[2] = pVec(du_max_beta)*(1-pVec(c_lim_beta))-yVec(dbetaL);
		du_max[3] = pVec(du_max_beta)*(1-pVec(c_lim_beta))-yVec(dbetaR);
		// Compute how much the actual joint velocity is above its lower limit:
		du_max[6] = pVec(du_max_beta)*(1-pVec(c_lim_beta))+yVec(dbetaL);
		du_max[7] = pVec(du_max_beta)*(1-pVec(c_lim_beta))+yVec(dbetaR);
	} else { // Actuator velocity is different from joint velocity
		// Compute how much the actual actuator velocity is below its upper limit:
		du_max[2] = pVec(du_max_beta)*(1-pVec(c_lim_beta))-yVec(dubetaL);
		du_max[3] = pVec(du_max_beta)*(1-pVec(c_lim_beta))-yVec(dubetaR);
		// Compute how much the actual actuator velocity is above its lower limit:
		du_max[6] = pVec(du_max_beta)*(1-pVec(c_lim_beta))+yVec(dubetaL);
		du_max[7] = pVec(du_max_beta)*(1-pVec(c_lim_beta))+yVec(dubetaR);
	}

	// Compute how much the actual actuator torque is below its upper limit:
	F_max[0] = pVec(c_lim_alpha)*pVec(P_max_alpha)/pVec(du_max_alpha) - uVec(TalphaL);
	F_max[1] = pVec(c_lim_alpha)*pVec(P_max_alpha)/pVec(du_max_alpha) - uVec(TalphaR);
	F_max[2] = pVec(c_lim_beta)*pVec(P_max_beta)/pVec(du_max_beta) - uVec(TbetaL);
	F_max[3] = pVec(c_lim_beta)*pVec(P_max_beta)/pVec(du_max_beta) - uVec(TbetaR);
	// Compute how much the actual actuator torque is above its lower limit:
	F_max[4] = pVec(c_lim_alpha)*pVec(P_max_alpha)/pVec(du_max_alpha) + uVec(TalphaL);
	F_max[5] = pVec(c_lim_alpha)*pVec(P_max_alpha)/pVec(du_max_alpha) + uVec(TalphaR);
	F_max[6] = pVec(c_lim_beta)*pVec(P_max_beta)/pVec(du_max_beta) + uVec(TbetaL);
	F_max[7] = pVec(c_lim_beta)*pVec(P_max_beta)/pVec(du_max_beta) + uVec(TbetaR);
}

// Update current states
void AB_Dynamics::setContState(double* yPtr){yVec = Map<VectorYd>(yPtr);}
void AB_Dynamics::setDiscState(double* zPtr){zVec = Map<VectorZd>(zPtr);}
void AB_Dynamics::setSystParam(double* pPtr){pVec = Map<VectorPd>(pPtr);}
void AB_Dynamics::setExctState(double* uPtr){uVec = Map<VectorUd>(uPtr);}
// Retrieve the current states
void AB_Dynamics::getContState(double* yPtr){
	for (int i = 0; i<NY; i++)
		yPtr[i] = yVec(i);
}
void AB_Dynamics::getDiscState(double* zPtr){
	for (int i = 0; i<NZ; i++)
		zPtr[i] = zVec(i);
}
void AB_Dynamics::getSystParam(double* pPtr){
	for (int i = 0; i<NP; i++)
		pPtr[i] = pVec(i);
}
void AB_Dynamics::getExctState(double* uPtr){
	for (int i = 0; i<NU; i++)
		uPtr[i] = uVec(i);
}

// Prepare all the components for the equations of motion:
void AB_Dynamics::updateDynamics(){
	// Before this function is called, we assume that all states and parameters have been updated.
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

// Components for the Equations of Motion
void AB_Dynamics::ComputeMassMatrix(){
	M = MatrixQd::Zero();
	double t2 = cos(yVec(phi));
	double t3 = yVec(alphaL)+yVec(betaL)+yVec(phi);
	double t4 = cos(t3);
	double t5 = yVec(alphaR)+yVec(betaR)+yVec(phi);
	double t6 = cos(t5);
	double t7 = yVec(alphaL)+yVec(phi);
	double t8 = cos(t7);
	double t9 = yVec(alphaR)+yVec(phi);
	double t10 = cos(t9);
	double t11 = pVec(lL2)*pVec(m3)*t4;
	double t12 = pVec(l2)*pVec(m2)*t8;
	double t13 = pVec(lL1)*pVec(m3)*t8;
	double t14 = pVec(lL2)*pVec(m3)*t6;
	double t15 = pVec(l2)*pVec(m2)*t10;
	double t16 = pVec(lL1)*pVec(m3)*t10;
	double t17 = pVec(l3)-pVec(lL2);
	double t18 = pVec(m2)*2.0;
	double t19 = pVec(m3)*2.0;
	double t20 = pVec(m1)+t18+t19;
	double t21 = sin(yVec(phi));
	double t22 = sin(t3);
	double t23 = sin(t5);
	double t24 = sin(t7);
	double t25 = sin(t9);
	double t26 = pVec(lL2)*pVec(m3)*t22;
	double t27 = pVec(l2)*pVec(m2)*t24;
	double t28 = pVec(lL1)*pVec(m3)*t24;
	double t29 = pVec(lL2)*pVec(m3)*t23;
	double t30 = pVec(l2)*pVec(m2)*t25;
	double t31 = pVec(lL1)*pVec(m3)*t25;
	double t32 = pVec(lH)*pVec(m2)*t2*2.0;
	double t33 = pVec(lH)*pVec(m3)*t2*2.0;
	double t59 = pVec(l3)*pVec(m3)*t4;
	double t69 = pVec(l3)*pVec(m3)*t6;
	double t34 = t11+t12+t13+t14+t15+t16+t32+t33-t59-t69;
	double t35 = pVec(lH)*pVec(m2)*t21*2.0;
	double t36 = pVec(lH)*pVec(m3)*t21*2.0;
	double t60 = pVec(l3)*pVec(m3)*t22;
	double t70 = pVec(l3)*pVec(m3)*t23;
	double t37 = t26+t27+t28+t29+t30+t31+t35+t36-t60-t70;
	double t38 = pVec(lH)*pVec(lH);
	double t39 = yVec(alphaL)+yVec(betaL);
	double t40 = cos(t39);
	double t41 = yVec(alphaR)+yVec(betaR);
	double t42 = cos(t41);
	double t43 = cos(yVec(alphaL));
	double t44 = cos(yVec(alphaR));
	double t45 = cos(yVec(betaL));
	double t46 = cos(yVec(betaR));
	double t47 = pVec(l2)*pVec(l2);
	double t48 = pVec(l3)*pVec(l3);
	double t49 = pVec(lL1)*pVec(lL1);
	double t50 = pVec(lL2)*pVec(lL2);
	double t51 = pVec(lL1)*pVec(lL2)*pVec(m3)*t45*2.0;
	double t52 = pVec(m3)*t48;
	double t53 = pVec(m3)*t50;
	double t54 = pVec(lH)*pVec(lL2)*pVec(m3)*t40;
	double t55 = pVec(m2)*t47;
	double t56 = pVec(m3)*t49;
	double t57 = pVec(lL1)*pVec(lL2)*pVec(m3)*t46*2.0;
	double t58 = pVec(lH)*pVec(lL2)*pVec(m3)*t42;
	double t61 = pVec(l2)*pVec(lH)*pVec(m2)*t43;
	double t62 = pVec(lH)*pVec(lL1)*pVec(m3)*t43;
	double t64 = pVec(l3)*pVec(lL2)*pVec(m3)*2.0;
	double t65 = pVec(l3)*pVec(lL1)*pVec(m3)*t45*2.0;
	double t67 = pVec(l3)*pVec(lH)*pVec(m3)*t40;
	double t63 = pVec(j2)+pVec(j3)+t51+t52+t53+t54+t55+t56+t61+t62-t64-t65-t67;
	double t66 = pVec(lL1)*pVec(lL2)*pVec(m3)*t45;
	double t68 = pVec(j3)+t52+t53-t64+t66-pVec(l3)*pVec(lL1)*pVec(m3)*t45;
	double t71 = pVec(l2)*pVec(lH)*pVec(m2)*t44;
	double t72 = pVec(lH)*pVec(lL1)*pVec(m3)*t44;
	double t73 = pVec(lL1)*pVec(lL2)*pVec(m3)*t46;
	double t74 = pVec(j3)+t52+t53-t64+t73-pVec(l3)*pVec(lL1)*pVec(m3)*t46;
	double t75 = pVec(j3)+t52+t53-t64;
	M(0,0) = t20;
	M(0,2) = t34;
	M(0,3) = t11+t12+t13-pVec(l3)*pVec(m3)*t4;
	M(0,4) = -pVec(m3)*t4*t17;
	M(0,5) = t14+t15+t16-pVec(l3)*pVec(m3)*t6;
	M(0,6) = -pVec(m3)*t6*t17;
	M(1,1) = t20;
	M(1,2) = t37;
	M(1,3) = t26+t27+t28-pVec(l3)*pVec(m3)*t22;
	M(1,4) = -pVec(m3)*t17*t22;
	M(1,5) = t29+t30+t31-pVec(l3)*pVec(m3)*t23;
	M(1,6) = -pVec(m3)*t17*t23;
	M(2,0) = t34;
	M(2,1) = t37;
	M(2,2) = pVec(j1_)+pVec(j2)*2.0+pVec(j3)*2.0+t51+t57+pVec(m2)*t38*2.0+pVec(m3)*t38*2.0+pVec(m2)*t47*2.0+pVec(m3)*t48*2.0+pVec(m3)*t49*2.0+pVec(m3)*t50*2.0-pVec(l3)*pVec(lL2)*pVec(m3)*4.0-pVec(l3)*pVec(lH)*pVec(m3)*t40*2.0+pVec(l2)*pVec(lH)*pVec(m2)*t43*2.0+pVec(l2)*pVec(lH)*pVec(m2)*t44*2.0-pVec(l3)*pVec(lH)*pVec(m3)*t42*2.0-pVec(l3)*pVec(lL1)*pVec(m3)*t45*2.0-pVec(l3)*pVec(lL1)*pVec(m3)*t46*2.0+pVec(lH)*pVec(lL2)*pVec(m3)*t40*2.0+pVec(lH)*pVec(lL1)*pVec(m3)*t43*2.0+pVec(lH)*pVec(lL2)*pVec(m3)*t42*2.0+pVec(lH)*pVec(lL1)*pVec(m3)*t44*2.0;
	M(2,3) = t63;
	M(2,4) = pVec(j3)+t52+t53+t54+t66-pVec(l3)*pVec(lL2)*pVec(m3)*2.0-pVec(l3)*pVec(lH)*pVec(m3)*t40-pVec(l3)*pVec(lL1)*pVec(m3)*t45;
	M(2,5) = pVec(j2)+pVec(j3)+t52+t53+t55+t56+t57+t58+t71+t72-pVec(l3)*pVec(lL2)*pVec(m3)*2.0-pVec(l3)*pVec(lH)*pVec(m3)*t42-pVec(l3)*pVec(lL1)*pVec(m3)*t46*2.0;
	M(2,6) = pVec(j3)+t52+t53+t58+t73-pVec(l3)*pVec(lL2)*pVec(m3)*2.0-pVec(l3)*pVec(lH)*pVec(m3)*t42-pVec(l3)*pVec(lL1)*pVec(m3)*t46;
	M(3,0) = t11+t12+t13-t59;
	M(3,1) = t26+t27+t28-t60;
	M(3,2) = t63;
	M(3,3) = pVec(j2)+pVec(j3)+t51+t52+t53+t55+t56-t64-t65;
	M(3,4) = t68;
	M(4,0) = -pVec(m3)*t4*t17;
	M(4,1) = -pVec(m3)*t17*t22;
	M(4,2) = pVec(j3)+t52+t53+t54-t64+t66-t67-pVec(l3)*pVec(lL1)*pVec(m3)*t45;
	M(4,3) = t68;
	M(4,4) = t75;
	M(5,0) = t14+t15+t16-t69;
	M(5,1) = t29+t30+t31-t70;
	M(5,2) = pVec(j2)+pVec(j3)+t52+t53+t55+t56+t57+t58-t64+t71+t72-pVec(l3)*pVec(lH)*pVec(m3)*t42-pVec(l3)*pVec(lL1)*pVec(m3)*t46*2.0;
	M(5,5) = pVec(j2)+pVec(j3)+t52+t53+t55+t56+t57-t64-pVec(l3)*pVec(lL1)*pVec(m3)*t46*2.0;
	M(5,6) = t74;
	M(6,0) = -pVec(m3)*t6*t17;
	M(6,1) = -pVec(m3)*t17*t23;
	M(6,2) = pVec(j3)+t52+t53+t58-t64+t73-pVec(l3)*pVec(lH)*pVec(m3)*t42-pVec(l3)*pVec(lL1)*pVec(m3)*t46;
	M(6,5) = t74;
	M(6,6) = t75;
	// Increase inertia matrix values for DOFs with parallel elastic actuation:
	// HIP
	if (pVec(hip_jnt_type) == PEA){
		M(qalphaL,qalphaL) = M(qalphaL,qalphaL) + pVec(j_mot_alpha);
		M(qalphaR,qalphaR) = M(qalphaR,qalphaR) + pVec(j_mot_alpha);
	}
	// LEG
	if (pVec(knee_jnt_type) == PEA){
		M(qbetaL,qbetaL) = M(qbetaL,qbetaL) + pVec(j_mot_beta);
		M(qbetaR,qbetaR) = M(qbetaR,qbetaR) + pVec(j_mot_beta);
	}
	// directly decompose the Mass-matrix:
	Mdecomp = M.ldlt();
}
void AB_Dynamics::ComputeDiffForces(){
	h = VectorQd::Zero();
	double t2 = yVec(alphaL)+yVec(betaL)+yVec(phi);
	double t3 = sin(t2);
	double t4 = yVec(alphaR)+yVec(betaR)+yVec(phi);
	double t5 = sin(t4);
	double t6 = yVec(dphi)*yVec(dphi);
	double t7 = yVec(dalphaL)*yVec(dalphaL);
	double t8 = yVec(dalphaR)*yVec(dalphaR);
	double t9 = yVec(dbetaL)*yVec(dbetaL);
	double t10 = yVec(dbetaR)*yVec(dbetaR);
	double t11 = yVec(alphaL)+yVec(phi);
	double t12 = sin(t11);
	double t13 = yVec(alphaR)+yVec(phi);
	double t14 = sin(t13);
	double t15 = sin(yVec(phi));
	double t16 = cos(t11);
	double t17 = cos(t13);
	double t18 = cos(yVec(phi));
	double t19 = cos(t2);
	double t20 = cos(t4);
	double t21 = sin(yVec(alphaL));
	double t22 = sin(yVec(alphaR));
	double t23 = sin(yVec(betaL));
	double t24 = sin(yVec(betaR));
	double t25 = yVec(alphaL)+yVec(betaL);
	double t26 = sin(t25);
	double t27 = yVec(alphaR)+yVec(betaR);
	double t28 = sin(t27);
	double t29 = pVec(g)*pVec(l3)*pVec(m3)*t3;
	double t30 = pVec(lL1)*pVec(lL2)*pVec(m3)*t9*t23;
	double t31 = yVec(dalphaL)*yVec(dbetaL)*pVec(lL1)*pVec(lL2)*pVec(m3)*t23*2.0;
	double t32 = yVec(dbetaL)*yVec(dphi)*pVec(lL1)*pVec(lL2)*pVec(m3)*t23*2.0;
	double t33 = pVec(g)*pVec(l3)*pVec(m3)*t5;
	double t34 = pVec(lL1)*pVec(lL2)*pVec(m3)*t10*t24;
	double t35 = yVec(dalphaR)*yVec(dbetaR)*pVec(lL1)*pVec(lL2)*pVec(m3)*t24*2.0;
	double t36 = yVec(dbetaR)*yVec(dphi)*pVec(lL1)*pVec(lL2)*pVec(m3)*t24*2.0;
	double t37 = pVec(l3)-pVec(lL2);
	h(0,0) = -pVec(l3)*pVec(m3)*t3*t6-pVec(l3)*pVec(m3)*t3*t7-pVec(l3)*pVec(m3)*t5*t6-pVec(l3)*pVec(m3)*t3*t9-pVec(l3)*pVec(m3)*t5*t8-pVec(l3)*pVec(m3)*t5*t10+pVec(l2)*pVec(m2)*t6*t12+pVec(l2)*pVec(m2)*t7*t12+pVec(l2)*pVec(m2)*t6*t14+pVec(l2)*pVec(m2)*t8*t14+pVec(lH)*pVec(m2)*t6*t15*2.0+pVec(lH)*pVec(m3)*t6*t15*2.0+pVec(lL2)*pVec(m3)*t3*t6+pVec(lL2)*pVec(m3)*t3*t7+pVec(lL2)*pVec(m3)*t5*t6+pVec(lL2)*pVec(m3)*t3*t9+pVec(lL2)*pVec(m3)*t5*t8+pVec(lL2)*pVec(m3)*t5*t10+pVec(lL1)*pVec(m3)*t6*t12+pVec(lL1)*pVec(m3)*t7*t12+pVec(lL1)*pVec(m3)*t6*t14+pVec(lL1)*pVec(m3)*t8*t14-yVec(dalphaL)*yVec(dbetaL)*pVec(l3)*pVec(m3)*t3*2.0-yVec(dalphaR)*yVec(dbetaR)*pVec(l3)*pVec(m3)*t5*2.0-yVec(dalphaL)*yVec(dphi)*pVec(l3)*pVec(m3)*t3*2.0+yVec(dalphaL)*yVec(dphi)*pVec(l2)*pVec(m2)*t12*2.0-yVec(dalphaR)*yVec(dphi)*pVec(l3)*pVec(m3)*t5*2.0+yVec(dalphaR)*yVec(dphi)*pVec(l2)*pVec(m2)*t14*2.0-yVec(dbetaL)*yVec(dphi)*pVec(l3)*pVec(m3)*t3*2.0-yVec(dbetaR)*yVec(dphi)*pVec(l3)*pVec(m3)*t5*2.0+yVec(dalphaL)*yVec(dbetaL)*pVec(lL2)*pVec(m3)*t3*2.0+yVec(dalphaR)*yVec(dbetaR)*pVec(lL2)*pVec(m3)*t5*2.0+yVec(dalphaL)*yVec(dphi)*pVec(lL2)*pVec(m3)*t3*2.0+yVec(dalphaL)*yVec(dphi)*pVec(lL1)*pVec(m3)*t12*2.0+yVec(dalphaR)*yVec(dphi)*pVec(lL2)*pVec(m3)*t5*2.0+yVec(dalphaR)*yVec(dphi)*pVec(lL1)*pVec(m3)*t14*2.0+yVec(dbetaL)*yVec(dphi)*pVec(lL2)*pVec(m3)*t3*2.0+yVec(dbetaR)*yVec(dphi)*pVec(lL2)*pVec(m3)*t5*2.0;
	h(1,0) = -pVec(g)*pVec(m1)-pVec(g)*pVec(m2)*2.0-pVec(g)*pVec(m3)*2.0-pVec(l2)*pVec(m2)*t6*t16-pVec(l2)*pVec(m2)*t6*t17-pVec(l2)*pVec(m2)*t7*t16-pVec(l2)*pVec(m2)*t8*t17+pVec(l3)*pVec(m3)*t6*t19+pVec(l3)*pVec(m3)*t6*t20+pVec(l3)*pVec(m3)*t7*t19+pVec(l3)*pVec(m3)*t8*t20+pVec(l3)*pVec(m3)*t9*t19+pVec(l3)*pVec(m3)*t10*t20-pVec(lH)*pVec(m2)*t6*t18*2.0-pVec(lH)*pVec(m3)*t6*t18*2.0-pVec(lL1)*pVec(m3)*t6*t16-pVec(lL1)*pVec(m3)*t6*t17-pVec(lL1)*pVec(m3)*t7*t16-pVec(lL1)*pVec(m3)*t8*t17-pVec(lL2)*pVec(m3)*t6*t19-pVec(lL2)*pVec(m3)*t6*t20-pVec(lL2)*pVec(m3)*t7*t19-pVec(lL2)*pVec(m3)*t8*t20-pVec(lL2)*pVec(m3)*t9*t19-pVec(lL2)*pVec(m3)*t10*t20+yVec(dalphaL)*yVec(dbetaL)*pVec(l3)*pVec(m3)*t19*2.0+yVec(dalphaR)*yVec(dbetaR)*pVec(l3)*pVec(m3)*t20*2.0-yVec(dalphaL)*yVec(dphi)*pVec(l2)*pVec(m2)*t16*2.0+yVec(dalphaL)*yVec(dphi)*pVec(l3)*pVec(m3)*t19*2.0-yVec(dalphaR)*yVec(dphi)*pVec(l2)*pVec(m2)*t17*2.0+yVec(dalphaR)*yVec(dphi)*pVec(l3)*pVec(m3)*t20*2.0+yVec(dbetaL)*yVec(dphi)*pVec(l3)*pVec(m3)*t19*2.0+yVec(dbetaR)*yVec(dphi)*pVec(l3)*pVec(m3)*t20*2.0-yVec(dalphaL)*yVec(dbetaL)*pVec(lL2)*pVec(m3)*t19*2.0-yVec(dalphaR)*yVec(dbetaR)*pVec(lL2)*pVec(m3)*t20*2.0-yVec(dalphaL)*yVec(dphi)*pVec(lL1)*pVec(m3)*t16*2.0-yVec(dalphaL)*yVec(dphi)*pVec(lL2)*pVec(m3)*t19*2.0-yVec(dalphaR)*yVec(dphi)*pVec(lL1)*pVec(m3)*t17*2.0-yVec(dalphaR)*yVec(dphi)*pVec(lL2)*pVec(m3)*t20*2.0-yVec(dbetaL)*yVec(dphi)*pVec(lL2)*pVec(m3)*t19*2.0-yVec(dbetaR)*yVec(dphi)*pVec(lL2)*pVec(m3)*t20*2.0;
	h(2,0) = t29+t30+t31+t32+t33+t34+t35+t36-pVec(g)*pVec(l2)*pVec(m2)*t12-pVec(g)*pVec(l2)*pVec(m2)*t14-pVec(g)*pVec(lH)*pVec(m2)*t15*2.0-pVec(g)*pVec(lH)*pVec(m3)*t15*2.0-pVec(g)*pVec(lL2)*pVec(m3)*t3-pVec(g)*pVec(lL2)*pVec(m3)*t5-pVec(g)*pVec(lL1)*pVec(m3)*t12-pVec(g)*pVec(lL1)*pVec(m3)*t14+pVec(l2)*pVec(lH)*pVec(m2)*t7*t21+pVec(l2)*pVec(lH)*pVec(m2)*t8*t22-pVec(l3)*pVec(lH)*pVec(m3)*t7*t26-pVec(l3)*pVec(lH)*pVec(m3)*t9*t26-pVec(l3)*pVec(lH)*pVec(m3)*t8*t28-pVec(l3)*pVec(lH)*pVec(m3)*t10*t28-pVec(l3)*pVec(lL1)*pVec(m3)*t9*t23-pVec(l3)*pVec(lL1)*pVec(m3)*t10*t24+pVec(lH)*pVec(lL1)*pVec(m3)*t7*t21+pVec(lH)*pVec(lL1)*pVec(m3)*t8*t22+pVec(lH)*pVec(lL2)*pVec(m3)*t7*t26+pVec(lH)*pVec(lL2)*pVec(m3)*t9*t26+pVec(lH)*pVec(lL2)*pVec(m3)*t8*t28+pVec(lH)*pVec(lL2)*pVec(m3)*t10*t28-yVec(dalphaL)*yVec(dbetaL)*pVec(l3)*pVec(lH)*pVec(m3)*t26*2.0-yVec(dalphaR)*yVec(dbetaR)*pVec(l3)*pVec(lH)*pVec(m3)*t28*2.0-yVec(dalphaL)*yVec(dbetaL)*pVec(l3)*pVec(lL1)*pVec(m3)*t23*2.0-yVec(dalphaR)*yVec(dbetaR)*pVec(l3)*pVec(lL1)*pVec(m3)*t24*2.0+yVec(dalphaL)*yVec(dphi)*pVec(l2)*pVec(lH)*pVec(m2)*t21*2.0-yVec(dalphaL)*yVec(dphi)*pVec(l3)*pVec(lH)*pVec(m3)*t26*2.0+yVec(dalphaR)*yVec(dphi)*pVec(l2)*pVec(lH)*pVec(m2)*t22*2.0-yVec(dalphaR)*yVec(dphi)*pVec(l3)*pVec(lH)*pVec(m3)*t28*2.0-yVec(dbetaL)*yVec(dphi)*pVec(l3)*pVec(lH)*pVec(m3)*t26*2.0-yVec(dbetaR)*yVec(dphi)*pVec(l3)*pVec(lH)*pVec(m3)*t28*2.0-yVec(dbetaL)*yVec(dphi)*pVec(l3)*pVec(lL1)*pVec(m3)*t23*2.0-yVec(dbetaR)*yVec(dphi)*pVec(l3)*pVec(lL1)*pVec(m3)*t24*2.0+yVec(dalphaL)*yVec(dbetaL)*pVec(lH)*pVec(lL2)*pVec(m3)*t26*2.0+yVec(dalphaR)*yVec(dbetaR)*pVec(lH)*pVec(lL2)*pVec(m3)*t28*2.0+yVec(dalphaL)*yVec(dphi)*pVec(lH)*pVec(lL1)*pVec(m3)*t21*2.0+yVec(dalphaL)*yVec(dphi)*pVec(lH)*pVec(lL2)*pVec(m3)*t26*2.0+yVec(dalphaR)*yVec(dphi)*pVec(lH)*pVec(lL1)*pVec(m3)*t22*2.0+yVec(dalphaR)*yVec(dphi)*pVec(lH)*pVec(lL2)*pVec(m3)*t28*2.0+yVec(dbetaL)*yVec(dphi)*pVec(lH)*pVec(lL2)*pVec(m3)*t26*2.0+yVec(dbetaR)*yVec(dphi)*pVec(lH)*pVec(lL2)*pVec(m3)*t28*2.0;
	h(3,0) = t29+t30+t31+t32-pVec(g)*pVec(l2)*pVec(m2)*t12-pVec(g)*pVec(lL2)*pVec(m3)*t3-pVec(g)*pVec(lL1)*pVec(m3)*t12-pVec(l2)*pVec(lH)*pVec(m2)*t6*t21+pVec(l3)*pVec(lH)*pVec(m3)*t6*t26-pVec(l3)*pVec(lL1)*pVec(m3)*t9*t23-pVec(lH)*pVec(lL1)*pVec(m3)*t6*t21-pVec(lH)*pVec(lL2)*pVec(m3)*t6*t26-yVec(dalphaL)*yVec(dbetaL)*pVec(l3)*pVec(lL1)*pVec(m3)*t23*2.0-yVec(dbetaL)*yVec(dphi)*pVec(l3)*pVec(lL1)*pVec(m3)*t23*2.0;
	h(4,0) = pVec(m3)*t37*(pVec(g)*t3+pVec(lH)*t6*t26+pVec(lL1)*t6*t23+pVec(lL1)*t7*t23+yVec(dalphaL)*yVec(dphi)*pVec(lL1)*t23*2.0);
	h(5,0) = t33+t34+t35+t36-pVec(g)*pVec(l2)*pVec(m2)*t14-pVec(g)*pVec(lL2)*pVec(m3)*t5-pVec(g)*pVec(lL1)*pVec(m3)*t14-pVec(l2)*pVec(lH)*pVec(m2)*t6*t22+pVec(l3)*pVec(lH)*pVec(m3)*t6*t28-pVec(l3)*pVec(lL1)*pVec(m3)*t10*t24-pVec(lH)*pVec(lL1)*pVec(m3)*t6*t22-pVec(lH)*pVec(lL2)*pVec(m3)*t6*t28-yVec(dalphaR)*yVec(dbetaR)*pVec(l3)*pVec(lL1)*pVec(m3)*t24*2.0-yVec(dbetaR)*yVec(dphi)*pVec(l3)*pVec(lL1)*pVec(m3)*t24*2.0;
	h(6,0) = pVec(m3)*t37*(pVec(g)*t5+pVec(lH)*t6*t28+pVec(lL1)*t6*t24+pVec(lL1)*t8*t24+yVec(dalphaR)*yVec(dphi)*pVec(lL1)*t24*2.0);
}
void  AB_Dynamics::ComputeContactPointL(){
	posL =  Vector2d::Zero();
	double t2 = yVec(alphaL)+yVec(phi);
	double t3 = yVec(alphaL)+yVec(betaL)+yVec(phi);
	posL(0,0) = yVec(x)+pVec(lH)*sin(yVec(phi))+pVec(lL1)*sin(t2)+pVec(lL2)*sin(t3);
	posL(1,0) = -pVec(rFoot)+yVec(y)-pVec(lH)*cos(yVec(phi))-pVec(lL1)*cos(t2)-pVec(lL2)*cos(t3);
}
void  AB_Dynamics::ComputeContactPointR(){
	posR =  Vector2d::Zero();
	double t2 = yVec(alphaR)+yVec(phi);
	double t3 = yVec(alphaR)+yVec(betaR)+yVec(phi);
	posR(0,0) = yVec(x)+pVec(lH)*sin(yVec(phi))+pVec(lL1)*sin(t2)+pVec(lL2)*sin(t3);
	posR(1,0) = -pVec(rFoot)+yVec(y)-pVec(lH)*cos(yVec(phi))-pVec(lL1)*cos(t2)-pVec(lL2)*cos(t3);
}
void AB_Dynamics::ComputeContactJacobianL(){
	double t2 = yVec(alphaL)+yVec(phi);
	double t3 = cos(t2);
	double t4 = pVec(lL1)*t3;
	double t5 = yVec(alphaL)+yVec(betaL)+yVec(phi);
	double t6 = cos(t5);
	double t7 = pVec(lL2)*t6;
	double t8 = sin(t2);
	double t9 = pVec(lL1)*t8;
	double t10 = sin(t5);
	double t11 = pVec(lL2)*t10;
	JL(0,0) = 1.0;
	JL(0,2) = pVec(rFoot)+t4+t7+pVec(lH)*cos(yVec(phi));
	JL(0,3) = pVec(rFoot)+t4+t7;
	JL(0,4) = pVec(rFoot)+t7;
	JL(1,1) = 1.0;
	JL(1,2) = t9+t11+pVec(lH)*sin(yVec(phi));
	JL(1,3) = t9+t11;
	JL(1,4) = t11;
}
void AB_Dynamics::ComputeContactJacobianR(){
	double t2 = yVec(alphaR)+yVec(phi);
	double t3 = cos(t2);
	double t4 = pVec(lL1)*t3;
	double t5 = yVec(alphaR)+yVec(betaR)+yVec(phi);
	double t6 = cos(t5);
	double t7 = pVec(lL2)*t6;
	double t8 = sin(t2);
	double t9 = pVec(lL1)*t8;
	double t10 = sin(t5);
	double t11 = pVec(lL2)*t10;
	JR(0,0) = 1.0;
	JR(0,2) = pVec(rFoot)+t4+t7+pVec(lH)*cos(yVec(phi));
	JR(0,5) = pVec(rFoot)+t4+t7;
	JR(0,6) = pVec(rFoot)+t7;
	JR(1,1) = 1.0;
	JR(1,2) = t9+t11+pVec(lH)*sin(yVec(phi));
	JR(1,5) = t9+t11;
	JR(1,6) = t11;
}
void  AB_Dynamics::ComputeContactJacobianDtTIMESdqdtL(){
	dJLdtTimesdqdt =  Vector2d::Zero();
	double t2 = yVec(alphaL)+yVec(phi);
	double t3 = sin(t2);
	double t4 = pVec(lL1)*t3;
	double t5 = yVec(alphaL)+yVec(betaL)+yVec(phi);
	double t6 = sin(t5);
	double t7 = pVec(lL2)*t6;
	double t8 = t4+t7;
	double t9 = yVec(dalphaL)*t8;
	double t10 = yVec(dbetaL)*pVec(lL2)*t6;
	double t11 = cos(t2);
	double t12 = pVec(lL1)*t11;
	double t13 = cos(t5);
	double t14 = pVec(lL2)*t13;
	double t15 = t12+t14;
	double t16 = yVec(dalphaL)*t15;
	double t17 = yVec(dbetaL)*pVec(lL2)*t13;
	double t18 = yVec(dalphaL)+yVec(dbetaL)+yVec(dphi);
	dJLdtTimesdqdt(0,0) = -yVec(dalphaL)*(t9+t10+yVec(dphi)*t8)-yVec(dphi)*(t9+t10+yVec(dphi)*(t4+t7+pVec(lH)*sin(yVec(phi))))-yVec(dbetaL)*pVec(lL2)*t6*t18;
	dJLdtTimesdqdt(1,0) = yVec(dalphaL)*(t16+t17+yVec(dphi)*t15)+yVec(dphi)*(t16+t17+yVec(dphi)*(t12+t14+pVec(lH)*cos(yVec(phi))))+yVec(dbetaL)*pVec(lL2)*t13*t18;
}
void  AB_Dynamics::ComputeContactJacobianDtTIMESdqdtR(){
	dJRdtTimesdqdt =  Vector2d::Zero();
	double t2 = yVec(alphaR)+yVec(phi);
	double t3 = sin(t2);
	double t4 = pVec(lL1)*t3;
	double t5 = yVec(alphaR)+yVec(betaR)+yVec(phi);
	double t6 = sin(t5);
	double t7 = pVec(lL2)*t6;
	double t8 = t4+t7;
	double t9 = yVec(dalphaR)*t8;
	double t10 = yVec(dbetaR)*pVec(lL2)*t6;
	double t11 = cos(t2);
	double t12 = pVec(lL1)*t11;
	double t13 = cos(t5);
	double t14 = pVec(lL2)*t13;
	double t15 = t12+t14;
	double t16 = yVec(dalphaR)*t15;
	double t17 = yVec(dbetaR)*pVec(lL2)*t13;
	double t18 = yVec(dalphaR)+yVec(dbetaR)+yVec(dphi);
	dJRdtTimesdqdt(0,0) = -yVec(dalphaR)*(t9+t10+yVec(dphi)*t8)-yVec(dphi)*(t9+t10+yVec(dphi)*(t4+t7+pVec(lH)*sin(yVec(phi))))-yVec(dbetaR)*pVec(lL2)*t6*t18;
	dJRdtTimesdqdt(1,0) = yVec(dalphaR)*(t16+t17+yVec(dphi)*t15)+yVec(dphi)*(t16+t17+yVec(dphi)*(t12+t14+pVec(lH)*cos(yVec(phi))))+yVec(dbetaR)*pVec(lL2)*t13*t18;
}
void AB_Dynamics::ComputeJointForces(){
	// Joint forces:
	tau = VectorQd::Zero();
	// Compute spring and damping forces:
	double T_spring_alphaL, T_spring_alphaR,
	       T_spring_betaL, T_spring_betaR;
	// Compute nonlinear damping coefficient in the knees
	double b_betaL, b_betaR;
	// HIP
	if (pVec(hip_jnt_type) == PEA){ // Parallel spring in the hip
		T_spring_alphaL = pVec(kalpha)*(0 - yVec(alphaL)) +
				              pVec(balpha) *(0 - yVec(dalphaL));
		T_spring_alphaR = pVec(kalpha)*(0 - yVec(alphaR)) +
				              pVec(balpha) *(0 - yVec(dalphaR));
		tau(qalphaL) = T_spring_alphaL + uVec(TalphaL);
		tau(qalphaR) = T_spring_alphaR + uVec(TalphaR);
	} else { // Serial spring in the hip
		T_spring_alphaL = pVec(kalpha)*(yVec(ualphaL)  - yVec(alphaL)) +
				              pVec(balpha) *(yVec(dualphaL) - yVec(dalphaL));
		T_spring_alphaR = pVec(kalpha)*(yVec(ualphaR)  - yVec(alphaR)) +
				              pVec(balpha) *(yVec(dualphaR) - yVec(dalphaR));
		tau(qalphaL) = T_spring_alphaL + 0;
		tau(qalphaR) = T_spring_alphaR + 0;
	}
	// LEG
	if (pVec(knee_jnt_type) == PEA){ // Parallel spring in the leg
		b_betaL = Sigmoid(-yVec(dbetaL),pVec(sigma))*((pVec(bbeta1) - pVec(bbeta2)) *
				  Sigmoid(pVec(betaSmall)-yVec(betaL),pVec(sigma)) + pVec(bbeta2)) + 
				  Sigmoid(yVec(dbetaL),pVec(sigma))*pVec(bbeta2);
		b_betaL = Sigmoid(-yVec(dbetaR),pVec(sigma))*((pVec(bbeta1) - pVec(bbeta2)) *
				  Sigmoid(pVec(betaSmall)-yVec(betaR),pVec(sigma)) + pVec(bbeta2)) + 
				  Sigmoid(yVec(dbetaR),pVec(sigma))*pVec(bbeta2);
		
		T_spring_betaL = pVec(kbeta1)*(0 - yVec(betaL)) + Logistic(yVec(betaL) - pVec(betaSmall),pVec(sigma))*(pVec(kbeta1) - pVec(kbeta2)) +
			              b_betaL *(0 - yVec(dbetaL));
		T_spring_betaR = pVec(kbeta1)*(0 - yVec(betaR)) + Logistic(yVec(betaR) - pVec(betaSmall),pVec(sigma))*(pVec(kbeta1) - pVec(kbeta2)) +
			              b_betaR *(0 - yVec(dbetaR));
		tau(qbetaL) = T_spring_betaL + uVec(TbetaL);
		tau(qbetaR) = T_spring_betaR + uVec(TbetaR);
	} else { // Serial spring in the leg
		b_betaL = Sigmoid(yVec(dubetaL)-yVec(dbetaL),pVec(sigma))*((pVec(bbeta1) - pVec(bbeta2)) *
				  Sigmoid(pVec(betaSmall)+yVec(dubetaL)-yVec(betaL),pVec(sigma)) + pVec(bbeta2)) + 
				  Sigmoid(yVec(dbetaL)-yVec(dubetaL),pVec(sigma))*pVec(bbeta2);
		b_betaR = Sigmoid(yVec(dubetaR)-yVec(dbetaR),pVec(sigma))*((pVec(bbeta1) - pVec(bbeta2)) *
				  Sigmoid(pVec(betaSmall)+yVec(dubetaR)-yVec(betaR),pVec(sigma)) + pVec(bbeta2)) + 
				  Sigmoid(yVec(dbetaR)-yVec(dubetaR),pVec(sigma))*pVec(bbeta2);
	
		T_spring_betaL = pVec(kbeta1)*(yVec(ubetaL) - yVec(betaL)) + Logistic(yVec(betaL) - yVec(ubetaL) - pVec(betaSmall),pVec(sigma))*(pVec(kbeta1) - pVec(kbeta2)) +
			              b_betaL *(yVec(dubetaL) - yVec(dbetaL));
		T_spring_betaR = pVec(kbeta1)*(yVec(ubetaR) - yVec(betaR)) + Logistic(yVec(betaR) - yVec(ubetaR) - pVec(betaSmall),pVec(sigma))*(pVec(kbeta1) - pVec(kbeta2)) +
			              b_betaR *(yVec(dubetaR) - yVec(dbetaR));
		tau(qbetaL) = T_spring_betaL + 0;
		tau(qbetaR) = T_spring_betaR + 0;
	}
}
double AB_Dynamics::Logistic(double x, double sigma){
	double y = x/sigma;
	
	if (abs(y) > 10){ // In linear domain
		return max(0.0,x);
	} else { // In logistic domain
		return sigma*log(1+exp(y));
	}
}
double AB_Dynamics::Sigmoid(double x, double sigma){
	double y = x/sigma;
	
	if (y > 10){ // In linear domain
		return 1;
	}
	else if (y < -10){
		return 0;
	}
	else {
		return 1/(1+exp(-y));
	}
}
		
		


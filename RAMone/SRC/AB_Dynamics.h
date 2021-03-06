/*
 * AB_Dynamics.h
 *
 * Common Dynamic Header for an Articulating Biped.
 *
 *  Created on: May 15, 2015
 *      Author: nilssmit@umich.edu
 */

#ifndef AB_DYNAMICS_H_
#define AB_DYNAMICS_H_

#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;


/* Coordinates, states, parameters, and events */
//Q
#define  NQ 7  /* Number of generalized coordinates */
/* WITH THESE COORDINATES IT IS IMPORTANT THAT THEY MATCH THE
 * DEFINITIONS IN THE CREATING M-FILE.  OTHERWISE THE AUTOGENERATED
 * FUNCTIONS WILL NOT WORK! */
enum qNames { qx, qy, qphi, qalphaL, qbetaL, qalphaR, qbetaR
};
//Y
#define  NY 25 /* Number of continuous states */
enum yNames { x, dx, y, dy, phi, dphi,
	          alphaL, dalphaL,ualphaL, dualphaL,
	          alphaR, dalphaR,ualphaR, dualphaR,
	          betaL, dbetaL, ubetaL, dubetaL,
	          betaR, dbetaR, ubetaR, dubetaR,
	          posActWork, posElWork, totElLoss
};
// Z
#define  NZ 2  /* Number of discrete states */
enum zNames { phaseL, phaseR };
enum phaseNames { flight, stance };
// P
#define  NP 34 /* Number of system parameters */
enum pNames { hip_jnt_type, knee_jnt_type,
              g, m1, m2, m3,
			  lH, lL1, lL2,
		      l2, l3, rFoot,
		      j1_, j2, j3,
			  kalpha, balpha,
			  kbeta1, kbeta2, bbeta1, bbeta2, betaSmall,
		      P_max_alpha, du_max_alpha, c_lim_alpha, j_mot_alpha,
			  P_max_beta, du_max_beta, c_lim_beta, j_mot_beta,
              cost_fct_sel, v_avg, sigma, const_sel              
};

// U
#define  NU 4 /* Number of excitation states*/
enum uNames { TalphaL, TalphaR,
	          TbetaL, TbetaR
};

// EVENTS
#define  NEV 5 /* Number of events */
enum evNames { touchdownL, touchdownR,
			   liftoffL, liftoffR,
	           apexTransit
};

// CONFIGURATION
enum config { SEA, PEA };

/* Definitions for the Eigen class */
// States and parameters
typedef Matrix<double, NY, 1> VectorYd; // for continuous states
typedef Matrix<double, NZ, 1> VectorZd; // for discrete states
typedef Matrix<double, NP, 1> VectorPd; // for system parameters
typedef Matrix<double, NU, 1> VectorUd; // for excitation states
// Coordinates (Mass matrix, Jacobians)
typedef Matrix<double, NQ, 1> VectorQd;   // for generalized coordinates
typedef Matrix<double, NQ, NQ> MatrixQd;  // for the mass matrix
typedef LDLT<MatrixQd> LDLTQd;            // for the decomposed mass matrix
typedef Matrix<double, 2, NQ> Matrix2Qd;  // for the individual Jacobians
// already defined in Eigen typedef Matrix<double, 2, 1>  Vector2d;   // for the individual jacobian derivatives

typedef Matrix<double, Dynamic, NQ> MatrixXQd; // for the compound jacobian
// already defined in Eigen typedef Matrix<double, Dynamic, 1>  VectorXd;  // for the compound jacobian derivative

// Specific implementations will have to inherent from this class
class AB_Dynamics{
public:
	AB_Dynamics();
	virtual ~AB_Dynamics();
	// Set states
	void setContState(double* yPtr);
	void setDiscState(double* zPtr);
	void setSystParam(double* pPtr);
	void setExctState(double* uPtr);
	// Get states
	void getContState(double* yPtr);
	void getDiscState(double* zPtr);
	void getSystParam(double* pPtr);
	void getExctState(double* uPtr);
	// Dynamic Functions
	void FlowMap(double* dydtPtr);
	void JumpMap(double* yPlusPtr);
	void JumpSet(double* eventValPtr);
	// Constraints
	void Constraints(double* gcL, double* gcR, double* du_max, double* F_max);

private:

	// States and parameters:
	VectorYd yVec;
	VectorZd zVec;
	VectorPd pVec;
	VectorUd uVec;

	// Parameter processing:
	void ComputeDependentParameters();
	// Components of Dynamic Equations
	void ComputeMassMatrix();
	void ComputeDiffForces();
	// Contact handling:
	void ComputeContactPointL();
	void ComputeContactPointR();
	void ComputeContactJacobianL();
	void ComputeContactJacobianR();
	void ComputeContactJacobianDtTIMESdqdtL();
	void ComputeContactJacobianDtTIMESdqdtR();
	// Forces in the Joints & Springs
	void ComputeJointForces();
	// Full update of the dynamics
	void updateDynamics();
	// Implement logistic and sigmoid functions
	double Logistic(double x, double sigma);
	double Sigmoid(double x, double sigma);

	// Variables for dynamic equations:
	MatrixQd  M;   // Mass matrix
	LDLTQd	  Mdecomp;  // ... and its decomposition
	VectorQd  h;   // differentiable forces
	VectorQd  tau; // generalized input forces
	// Contact processing
	Vector2d  posL; // contact point position left
	Vector2d  posR; // contact point position right
	Matrix2Qd JL;   // contact Jacobian left
	Matrix2Qd JR;   // contact Jacobian right
	Vector2d  dJLdtTimesdqdt;  // derivative of Jacobian times velocity left
	Vector2d  dJRdtTimesdqdt;  // derivative of Jacobian times velocity right
	// Contact compounds:
	int nContact;			   // The number of feet that are in contact
	MatrixXQd J;   			   // Compound contact Jacobian
	VectorXd  dJdtTimesdqdt;   // Compound derivative of Jacobian times velocity
};

#endif /* AB_DYNAMICS_H_ */

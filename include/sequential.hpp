#ifndef _SEQUENTIAL_H
#define _SEQUENTIAL_H

#include "../LinearAlgebra/Matrix.hpp"
#include "../LinearAlgebra/LinearSolvers.hpp"

#include "DynamicsModel.hpp"
#include "ControlModel.hpp"
#include "ProcessNoiseModel.hpp"
#include "MeasurementModel.hpp"


namespace SequentialFilter{

// The Classical Kalman Filter (linear)
class Kalman{
public:

	Kalman(){};

	void propagate(Vector & x, Matrix & Pxx, double dt, Vector z, 
				   const MeasurementModel & msmt,
				   const DynamicsModel & dyn,
				   const ControlModel & ctrl,
				   const ProcessNoiseModel & pnoise){

		// propagate using the dynamics model, control and process noise
		Matrix STM = dyn.get_stm(dt, x);
		Vector u = ctrl.get_control(dt, x);
		Matrix G = ctrl.get_gain(dt, x);
		Vector vbar = pnoise.get_mean(dt, x);
		Matrix Q = pnoise.get_covariance(dt, x);
		Matrix PNTM = pnoise.get_pntm(dt, x);

		// time update
		m_prediction = STM*x + G*u + PNTM*vbar;
		Matrix Pbar = STM*Pxx*(~STM) + PNTM*Q*(~PNTM);

		// get data from measurement model
		Matrix H = msmt.get_jacobian(dt, m_prediction);
		Matrix R = msmt.get_covariance(dt, m_prediction);

		// measurement update
		Matrix K = Pbar*(~H)*inv(H*Pbar*(~H)+R);
		m_postfit = K*(z-H*m_prediction)+m_prediction;
		x = m_postfit;
		Pxx = Pbar - K*H*Pbar;
		return;
	}

	Vector postfit_resid(const Vector & z, const MeasurementModel & msmt) const;
	Vector prediction_resid(const Vector & z, const MeasurementModel & msmt) const;

private:

	Vector m_prediction;
	Vector m_postfit;

};

// The first-order nonlinear perturbation to the
// Classical Kalman Filter
class KalmanNonlin{
public:

	KalmanNonlin(){};

	void propagate(Vector & x, Matrix & Pxx, double dt, Vector z, 
				   const MeasurementModel & msmt,
				   const DynamicsModel & dyn,
				   const ControlModel & ctrl,
				   const ProcessNoiseModel & pnoise){

		Vector dx(x); dx.fill(0);

		// propagate using the dynamics model, control and process noise
		Vector u = ctrl.get_control(dt, x);
		Matrix G = ctrl.get_gain(dt, x);
		Vector vbar = pnoise.get_mean(dt, x);
		Matrix Q = pnoise.get_covariance(dt, x);
		Matrix PNTM = pnoise.get_pntm(dt, x);
		Matrix STM = dyn.get_stm(dt, x);
		m_prediction = dyn.propagate(dt, x, u, vbar);
		
		// time update
		Vector dxbar = STM*dx + G*u + PNTM*vbar;
		Matrix Pbar = STM*Pxx*(~STM) + PNTM*Q*(~PNTM);

		// get data from measurement model
		Vector h = msmt.get_value(dt, m_prediction);
		Matrix H = msmt.get_jacobian(dt, m_prediction);
		Matrix R = msmt.get_covariance(dt, m_prediction);
		Vector dz = z-h;

		// measurement update
		Matrix K = Pbar*(~H)*inv(H*Pbar*(~H)+R);
		dx = K*(dz-H*dxbar)+dxbar;
		x = m_prediction+dx;
		Pxx = Pbar - K*H*Pbar;

		return;
	}

	Vector postfit_resid(const Vector & z, const MeasurementModel & msmt) const;
	Vector prediction_resid(const Vector & z, const MeasurementModel & msmt) const;

private:

	Vector m_prediction;
	Vector m_postfit;

};


// The Extended Kalman Filter (Nonlinear)
class ExtendedKalman{
public:
	ExtendedKalman(){};

	void propagate(Vector & x, Matrix & Pxx, double dt, Vector z, 
				   const MeasurementModel & msmt,
				   const DynamicsModel & dyn,
				   const ControlModel & ctrl,
				   const ProcessNoiseModel & pnoise){

		// propagate using the dynamics model, control and process noise
		Vector u = ctrl.get_control(dt, x);
		Matrix G = ctrl.get_gain(dt, x);
		Vector vbar = pnoise.get_mean(dt, x);
		Matrix Q = pnoise.get_covariance(dt, x);
		Matrix PNTM = pnoise.get_pntm(dt, x);
		Matrix STM = dyn.get_stm(dt, x);
		m_prediction = dyn.propagate(dt, x, u, vbar);
		
		// time update
		Matrix Pbar = STM*Pxx*(~STM) + PNTM*Q*(~PNTM);

		// get data from measurement model
		Matrix H = msmt.get_jacobian(dt, m_prediction);
		Matrix R = msmt.get_covariance(dt, m_prediction);
		Vector dz = msmt.get_residual(m_prediction, z);

		// measurement update
		Matrix K = Pbar*(~H)*inv(H*Pbar*(~H)+R);
		x = m_prediction + K*dz;
		Pxx = Pbar - K*H*Pbar;

		return;
	}

	Vector postfit_resid(const Vector & z, const MeasurementModel & msmt) const;
	Vector prediction_resid(const Vector & z, const MeasurementModel & msmt) const;

private:

	Vector m_prediction;
	Vector m_postfit;

};

/*
// The Unscented Kalman Filter (Nonlinear)
class UnscentedKalman{

};

// The Information Filter
class Information{

};

// The Square-Root Information Filter
class SqrtInformation{

};
*/

}


#endif
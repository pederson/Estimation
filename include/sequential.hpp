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


// The Unscented Kalman Filter (Nonlinear)
class UnscentedKalman{
public:
	UnscentedKalman(){};

	void propagate(Vector & x, Matrix & Pxx, double dt, Vector z, 
				   const MeasurementModel & msmt,
				   const DynamicsModel & dyn,
				   const ControlModel & ctrl,
				   const ProcessNoiseModel & pnoise){

		// get requisite data from process noise and control
		Vector v = pnoise.get_mean(dt, x);
		Matrix Q = pnoise.get_covariance(dt, x);
		Vector u = ctrl.get_control(dt, x);
		Matrix R = msmt.get_covariance(dt, x);

		unsigned int n = x.length();
		unsigned int d = v.length();
		unsigned int m = z.length();

		// generate sigma points
		Matrix P(n+d+m,n+d+m); P.fill(0);
		P(0,n-1,0,n-1) = Pxx;
		P(n,n+d-1, n,n+d-1) = Q;
		P(n+d,n+d+m-1, n+d,n+d+m-1) = R;
		Vector xaug(x.length()+v.length()+z.length()); xaug.fill(0);
		xaug.subcol(0,0,x.length()-1) = x;
		xaug.subcol(0,x.length(),x.length()+v.length()-1) = v;
		Matrix xi; double w0_m, w0_c, wi;
		SigmaPoints(xaug, P, 1.0, 3.0, 3.0-(n+d+m), xi, w0_m, w0_c, wi);
		unsigned int nstate = xi.rows();
		unsigned int nsig = xi.cols();

		// do time update from dynamics
		for (auto i=0; i<nsig; i++) xi.subcol(i, 0, n-1) = dyn.propagate(dt, Vector(xi.subcol(i, 0, n-1)), u, Vector(xi.subcol(i, n, n+d-1)));

		// calculate xbar and Pxxbar
		Vector xbar = w0_m*xi.subcol(0,0,n-1);
		for (auto i=1; i<nsig; i++) xbar += wi*xi.subcol(i, 0, n-1);
		Matrix Pxxbar = w0_c*(xi.subcol(0, 0,n-1)-xbar)*~(xi.subcol(0,0,n-1)-xbar);
		for (auto i=1; i<nsig; i++) Pxxbar += wi*(xi.subcol(i, 0,n-1)-xbar)*~(xi.subcol(i,0,n-1)-xbar);

		// generate Zi using measurement model
		Matrix Zi(z.length(), nsig); Zi.fill(0);
		for (auto i=0; i<nsig; i++) Zi.col(i) = msmt.get_value(dt, xi.subcol(i,0,n-1))+xi.subcol(i,n+d,n+d+m-1);
		
		// calculate zbar and Pzzbar
		Vector zbar = w0_m*Zi.col(0);
		for (auto i=1; i<nsig; i++) zbar += wi*Zi.col(i);
		Matrix Pzzbar = w0_c*(Zi.col(0)-zbar)*~(Zi.col(0)-zbar);
		for (auto i=1; i<nsig; i++) Pzzbar += wi*(Zi.col(i)-zbar)*~(Zi.col(i)-zbar);

		// calculate cross covariance Pxzbar
		Matrix Pxzbar = w0_c*(xi.subcol(0,0,n-1)-xbar)*~(Zi.col(0)-zbar);
		for (auto i=1; i<nsig; i++) Pxzbar += wi*(xi.subcol(i,0,n-1)-xbar)*~(Zi.col(i)-zbar);
		
		// measurement update
		Matrix K = Pxzbar*inv(Pzzbar);
		x = xbar + K*(z-zbar);
		Pxx = Pxxbar - K*Pzzbar*~K;

		return;
	}

	Vector postfit_resid(const Vector & z, const MeasurementModel & msmt) const;
	Vector prediction_resid(const Vector & z, const MeasurementModel & msmt) const;

private:

	void SigmaPoints(const Vector & x, const Vector & Pxx, double alpha, double beta, double kappa,
					 Matrix & points, double & w0_m, double & w0_c, double & wi){
		unsigned int n=x.length();
		double lambda = alpha*alpha*(n+kappa)-n;
		Matrix L;
		cholesky(Pxx, L);

		// form 2n+1 sigma points
		Matrix xi(n, 2*n+1);
		xi.col(0) = x;
		for (auto i=0; i<n; i++) xi.col(i+1) = x + sqrt(n+lambda)*L.col(i);
		for (auto i=0; i<n; i++) xi.col(i+n+1) = x - sqrt(n+lambda)*L.col(i);

		// fill in the weights
		w0_m = lambda/(n+lambda);
		w0_c = w0_m + (1-alpha*alpha+beta);
		wi = 0.5/(n+lambda);
		swap(xi,points);
		return;
	}
};

/*
// The Information Filter
class Information{

};

// The Square-Root Information Filter
class SqrtInformation{

};
*/

}


#endif
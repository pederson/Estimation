#include "include/sequential.hpp"

#include <iostream>

using namespace std;


class DynamicsOscillator : public DynamicsModel{
public:
	Matrix get_stm(double dt, const Vector & state) const{
		Matrix STM(2,2);
		STM(0,0) = cos(omega*dt);
		STM(0,1) = 1.0/omega*sin(omega*dt);
		STM(1,0) = -omega*sin(omega*dt);
		STM(1,1) = cos(omega*dt);
		return STM;
	}
	Vector propagate(double dt, const Vector & state, const Vector & control, const Vector & process_noise) const{
		Matrix STM = get_stm(dt, state);
		return STM*state;
	}

private:
	double omega=2.0;
};


class ProcessNoiseOscillator : public ProcessNoiseModel{
public:
	Vector get_mean(double dt, const Vector & state) const{
		Vector out(1);
		out(0) = 0;
		return out;
	}
	Matrix get_covariance(double dt, const Vector & state) const{
		Matrix cov(1,1);
		cov(0,0) = 1;
		return cov;
	}
	Matrix get_pntm(double dt, const Vector & state) const{
		Matrix PNTM(2,1);
		PNTM(0,0) = 1.0/(omega*omega)*(1.0-cos(omega*dt));
		PNTM(1,0) = 1.0/omega*sin(omega*dt);
		return PNTM;
	}

private:
	double omega=2.0;
};

class ZeroNoise : public ProcessNoiseModel{
public:
	Vector get_mean(double dt, const Vector & state) const{
		Vector out(1);
		out(0) = 0;
		return out;
	}
	Matrix get_covariance(double dt, const Vector & state) const{
		Matrix cov(1,1);
		cov(0,0) = 0;
		return cov;
	}
	Matrix get_pntm(double dt, const Vector & state) const{
		Matrix PNTM(2,1);
		PNTM(0,0) = 0;
		PNTM(1,0) = 0;
		return PNTM;
	}

private:
	double omega=2.0;
};


class MeasurementOscillator : public MeasurementModel{
public:
	Vector get_value(double dt, const Vector & state) const {
		double rho = sqrt(state(0)*state(0)+height*height);
		Vector h(2);
		h(0) = rho;
		h(1) = state(0)*state(1)/rho;
		return h;
	}
	Matrix get_jacobian(double dt, const Vector & state) const{
		Matrix H(2,2);
		double rho = sqrt(state(0)*state(0)+height*height);
		H(0,0) = state(0)/rho;
		H(0,1) = 0;
		H(1,0) = state(1)/rho*(1-(state(0)*state(0))/(rho*rho));
		H(1,1) = state(0)/rho;
		return H;
	}
	Matrix get_covariance(double dt, const Vector & state) const{
		Matrix R(2,2);
		R.fill(0);
		R(0,0) = 0.25*0.25;
		R(1,1) = 0.1*0.1;
		return R;
	}
	Vector get_residual(const Vector & state, const Vector & measurement) const{
		Vector h = get_value(0, state);
		return measurement - h;
	}
private:
	double height=5.4;
};


class ZeroControl : public ControlModel{
public:
	Vector get_control(double dt, const Vector & state) const{
		Vector out(state);
		out.fill(0);
		return out;
	}
	Matrix get_gain(double dt, const Vector & state) const{
		return eye(state.length(), state.length());
	}
};

int main(int argc, char * argv[]){

	// read data
	Matrix oscill = dlmread("data/oscillator.dat");
	// cout << oscill << endl;
	Vector tdat = oscill.col(0);
	Vector rdat = oscill.col(1);
	Vector rrdat = oscill.col(2);
	Matrix zdat(tdat.length(), 2); 
	zdat.col(0) = rdat;
	zdat.col(1) = rrdat;

	// build linear model and filter
	SequentialFilter::Kalman ckf;
	MeasurementOscillator msmt;
	DynamicsOscillator dyn;
	ProcessNoiseOscillator pnoise;
	ZeroNoise nonoise;
	ZeroControl ctrl;
	Vector x0(2); x0(0) = 4.5; x0(1) = 0.15;
	Matrix Pxx0(2,2); Pxx0.fill(0); Pxx0(0,0) = 1000; Pxx0(1,1) = 100;

	Vector x = x0;
	Matrix Pxx = Pxx0;
	double tprev=0.0;
	cout << "CKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	for (auto k=0; k<tdat.length(); k++){
		ckf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
		cout << "CKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
		tprev = tdat(k);
	}


	SequentialFilter::KalmanNonlin nkf;
	x = x0;
	Pxx = Pxx0;
	tprev=0.0;
	cout << "NKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	for (auto k=0; k<tdat.length(); k++){
		nkf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
		cout << "NKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
		tprev = tdat(k);
	}


	SequentialFilter::ExtendedKalman ekf;
	x = x0; x(0) = 4.9; x(1) = 0.25;
	Pxx = Pxx0;
	tprev=0.0;
	cout << "EKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	for (auto k=0; k<tdat.length(); k++){
		ekf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
		cout << "EKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
		tprev = tdat(k);
	}



	return 0;
}
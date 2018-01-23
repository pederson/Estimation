#include "include/sequential.hpp"

#include <iostream>

using namespace std;


struct MatrixInverter{
	static Matrix invert(const Matrix & m){
		return inv(m);
	}
};



class DynamicsOscillatorPulse : public DynamicsModel{
public:
	Matrix get_stm(double dt, const Vector & state) const{
		Matrix STM(3,3);
		STM.fill(0);
		STM(0,0) = cos(2*omega*dt);
		STM(0,1) = 1.0/(2*omega)*sin(2*omega*dt);
		STM(0,2) = 0.5*(1.0-cos(2*omega*dt));
		STM(1,0) = -2.0*omega*sin(2*omega*dt);
		STM(1,1) = cos(2*omega*dt);
		STM(1,2) = omega*sin(2*omega*dt);
		STM(2,2) = 1.0;
		return STM;
	}
	Vector propagate(double dt, const Vector & state, const Vector & control, const Vector & process_noise) const{
		Matrix STM = get_stm(dt, state);
		return STM*state;
	}

private:
	double omega=6.2832e9;
};


class ProcessNoiseOscillatorPulse : public ProcessNoiseModel{
public:
	Vector get_mean(double dt, const Vector & state) const{
		Vector out(1);
		out(0) = 0;
		return out;
	}
	Matrix get_covariance(double dt, const Vector & state) const{
		Matrix cov(1,1);
		cov(0,0) = 0.1*fabs(state(2)) + 1.0e-12;
		return cov;
	}
	Matrix get_pntm(double dt, const Vector & state) const{
		Matrix PNTM(3,1);
		PNTM.fill(0);
		PNTM(2,0) = 1.0;
		return PNTM;
	}

private:
	double omega=6.2832e9;
};

// class ZeroNoise : public ProcessNoiseModel{
// public:
// 	Vector get_mean(double dt, const Vector & state) const{
// 		Vector out(1);
// 		out(0) = 0;
// 		return out;
// 	}
// 	Matrix get_covariance(double dt, const Vector & state) const{
// 		Matrix cov(1,1);
// 		cov(0,0) = 0;
// 		return cov;
// 	}
// 	Matrix get_pntm(double dt, const Vector & state) const{
// 		Matrix PNTM(2,1);
// 		PNTM(0,0) = 0;
// 		PNTM(1,0) = 0;
// 		return PNTM;
// 	}

// private:
// 	double omega=6.2832e9;
// };


class MeasurementOscillatorPulse : public MeasurementModel{
public:
	Vector get_value(double dt, const Vector & state) const {
		// double rho = sqrt(state(0)*state(0)+height*height);
		Vector h(1);
		h(0) = state(0);
		return h;
	}
	Matrix get_jacobian(double dt, const Vector & state) const{
		Matrix H(1,3);
		H.fill(0);
		H(0,0) = 1;


		return H;
	}
	Matrix get_covariance(double dt, const Vector & state) const{
		Matrix R(1,1);
		R(0,0) = 1.0e-17;
		return R;
	}
	Vector get_residual(const Vector & state, const Vector & measurement) const{
		Vector h = get_value(0, state);
		return measurement - h;
	}
private:
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
	Matrix oscill = dlmread("data/oscillator_pulse.dat");
	// cout << oscill << endl;
	Vector tdat = oscill.col(0);
	Vector rdat = oscill.col(1);



	// build linear model and filter
	SequentialFilter::Kalman<MatrixInverter> ckf;
	MeasurementOscillatorPulse msmt;
	DynamicsOscillatorPulse dyn;
	ProcessNoiseOscillatorPulse pnoise;
	// ZeroNoise nonoise;
	ZeroControl ctrl;
	Vector x0(3); x0(0) = 0.0; x0(1) = 0.0; x0(2) = 0.01;
	Matrix Pxx0(3,3); Pxx0.fill(0); Pxx0(0,0) = 0.05*0.05; Pxx0(1,1) = 0.05*0.05; Pxx0(2,2) = 0.05*0.05;

	Vector x = x0;
	Matrix Pxx = Pxx0;
	double tprev=0.0;
	cout << "CKF:, t, x1, x2, x3, Pxx1, Pxx2, Pxx3" << endl;
	for (auto k=0; k<tdat.length(); k++){
		Vector ms(1);
		ms(0) = rdat(k)*rdat(k);
		ckf.propagate(x, Pxx, tdat(k)-tprev, ~ ms, msmt, dyn, ctrl, pnoise);
		cout << "CKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(2) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << ", " << Pxx(2,2) << endl;
		tprev = tdat(k);
	}


	// SequentialFilter::KalmanNonlin nkf;
	// x = x0;
	// Pxx = Pxx0;
	// tprev=0.0;
	// cout << "NKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	// for (auto k=0; k<tdat.length(); k++){
	// 	nkf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
	// 	cout << "NKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
	// 	tprev = tdat(k);
	// }


	// SequentialFilter::ExtendedKalman ekf;
	// x = x0; x(0) = 4.9; x(1) = 0.25;
	// Pxx = Pxx0;
	// tprev=0.0;
	// cout << "EKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	// for (auto k=0; k<tdat.length(); k++){
	// 	ekf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
	// 	cout << "EKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
	// 	tprev = tdat(k);
	// }


	// SequentialFilter::UnscentedKalman ukf;
	// x = x0; x(0) = 4.9; x(1) = 0.25;
	// Pxx = Pxx0;
	// tprev=0.0;
	// cout << "UKF:, t, x1, x2, err1, err2, Pxx1, Pxx2" << endl;
	// for (auto k=0; k<tdat.length(); k++){
	// 	ekf.propagate(x, Pxx, tdat(k)-tprev, ~ Vector(zdat.row(k)), msmt, dyn, ctrl, pnoise);
	// 	cout << "UKF:, " << tdat(k) << ", " << x(0) << ", " << x(1) << ", " << x(0)-oscill(k,3) << ", " << x(1)-oscill(k,4) << ", " << Pxx(0,0) << ", " << Pxx(1,1) << endl;
	// 	tprev = tdat(k);
	// }



	return 0;
}
#include "kalman_filter.h"
#include <iostream>
const float PI2 = 2*M_PI;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_in;
  R_laser_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_ ;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
}


VectorXd KalmanFilter::Calculate_z_pred() {
    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    float phi;
    float rho;
    float rhodot;

    rho = sqrt(px*px+py*py);
    phi = atan2(py,px);
    if(rho<0.000001)
        rho = 0.000001;
    rhodot = (px*vx+py*vy)/rho; 
    VectorXd z_pred = VectorXd(3);
    z_pred << rho,phi,rhodot;
    return z_pred;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    VectorXd z_pred = VectorXd(3);
    z_pred  = Calculate_z_pred();
	VectorXd y = z - z_pred;
    while(y(1)>M_PI){
        y(1)-=PI2;
    }
    while(y(1)<-M_PI){
        y(1)+=PI2;
    }
	MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_rader_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K*Hj_)*P_;
}

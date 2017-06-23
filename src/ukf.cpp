#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
  time_us_ = 0;

  // CRTV model used
  n_x_ = 5;
  n_aug_ = 7;

  // UKF free parameter
  lambda_ = 3 - n_x_;

  // Nb sigma points: 2*n_aug_+1
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  weights_ = VectorXd(2*n_aug_+1);

  //set weights
  weights_.fill(1.0);
  weights_ /= (2*(lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Based on ground truth values analysis: cf data/get_sigma.m
  std_a_ = 0.0714;
  std_yawdd_ = 0.097739;

  //std_a_ = 0.1;
  //std_yawdd_ = 0.1;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
    */
    // first measurement
    cout << "UKF: " << endl;

    P_ << 1000, 0, 0, 0, 0,
			    0, 1000, 0, 0, 0,
			    0, 0, 1000, 0, 0,
			    0, 0, 0, 1000, 0,
			    0, 0, 0, 0, 1000;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
		  x_ <<  rho * cos(phi), rho * sin(phi), 0, 0, 0;
      P_ << 1, 0, 0, 0, 0,
			      0, 1, 0, 0, 0,
			      0, 0, 10000, 0, 0,
			      0, 0, 0, 100, 0,
			      0, 0, 0, 0, 100;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		  //set the state with the initial location and zero velocity
		  x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      P_ << 1, 0, 0, 0, 0,
			      0, 1, 0, 0, 0,
			      0, 0, 10000, 0, 0,
			      0, 0, 0, 100, 0,
			      0, 0, 0, 0, 100;
    }

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	//compute the time elapsed between the current and previous measurements
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //------------------------------------------------------
  // 1) Create Augmented State, Augmented Covariance matrix 
  //------------------------------------------------------

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  // mean of last 2 elements, noise values, is 0

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  //------------------------------------------------------
  // 2) Select Augmented Sigma points
  //------------------------------------------------------

  // number of sigma points
  int n_sig = 2*n_aug_+1;
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig);

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  //create augmented sigma points
  float c1 = sqrt(lambda_ + n_aug_);
  Xsig_aug.middleCols(1, n_aug_) = c1 * A_aug;
  Xsig_aug.middleCols(n_aug_+1, n_aug_) = -c1 *A_aug;
  Xsig_aug += x_aug.replicate(1, n_sig);

  //-------------------------------------------------------------
  // 3) Apply non linear CRTV process motion model: dt dependant
  //-------------------------------------------------------------

  // This implements the CTRV process model
  //predict sigma points
  for (int n = 0; n < n_sig; n++) {
    float px = Xsig_aug(0, n);
    float py = Xsig_aug(1, n);
    float v = Xsig_aug(2, n);
    float psi = Xsig_aug(3, n);
    float psi_dot = Xsig_aug(4, n);

    float nu_a = Xsig_aug(5, n);
    float nu_yawdd = Xsig_aug(6, n);

    //avoid division by zero
    //write predicted sigma points into right column

    if (fabs(psi_dot) > 0.001) {
      Xsig_pred_(0, n) = px + (v/psi_dot)*(sin(psi+psi_dot*dt)-sin(psi));
      Xsig_pred_(1, n) = py + (v/psi_dot)*(-cos(psi+psi_dot*dt)+cos(psi));
    } else {
      Xsig_pred_(0, n) = px + v*cos(psi)*dt; 
      Xsig_pred_(1, n) = py + v*sin(psi)*dt; 
    }

    Xsig_pred_(2, n) = v + 0;
    Xsig_pred_(3, n) = psi + psi_dot*dt;
    Xsig_pred_(4, n) = psi_dot + 0;

    Xsig_pred_(0, n) += (0.5*dt*dt*cos(psi)*nu_a);
    Xsig_pred_(1, n) += (0.5*dt*dt*sin(psi)*nu_a);
    Xsig_pred_(2, n) += (dt*nu_a);
    Xsig_pred_(3, n) += (0.5*dt*dt*nu_yawdd);
    Xsig_pred_(4, n) += (dt*nu_yawdd);
  }  

  //----------------------------------------------
  // 4) Predict new Mean and new Covariance
  //----------------------------------------------
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x.fill(0.0);
  for (int i = 0; i < n_sig; i++) {
    x = x + weights_(i) * Xsig_pred_.col(i); 
  }

  //predict state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < n_sig; i++) {
    VectorXd xi = Xsig_pred_.col(i) - x; 

    //angle normalization
    while (xi(3)> M_PI) xi(3)-=2.*M_PI;
    while (xi(3)<-M_PI) xi(3)+=2.*M_PI;

    P = P + weights_(i) * xi * xi.transpose();
  }

  x_ = x;
  P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

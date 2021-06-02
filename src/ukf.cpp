#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    std_a_ = 9;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    /**
    * DO NOT MODIFY measurement noise values below.
    * These are provided by the sensor manufacturer.
    */

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
    * End DO NOT MODIFY section for measurement noise values
    */

    /**
    * TODO: Complete the initialization. See ukf.h for other member properties.
    * Hint: one or more values initialized above might be wildly off...
    */

    is_initialized_ = false;

    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;

    Xsig_pred_ = MatrixXd(n_aug_, 2*n_aug_+1);

    weights = VectorXd(2*n_aug_+1);
    weights.fill(0.5/(lambda_+n_aug_));
    weights(0) = lambda_/(lambda_+n_aug_);

    Rlidar_ = MatrixXd(2, 2);
    Rlidar_ << std_laspx_, 0, 0, std_laspy_;

    Rradar_ = MatrixXd(3, 3);
    Rradar_.fill(0.0);
    Rradar(0, 0) = std_radr_;
    Rradar(1, 1) = std_radphi_;
    Rradar(2, 2) = std_radrd_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Make sure you switch between lidar and radar
    * measurements.
    */
    if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            x_ << meas_package.raw_measurements_(0),
                  meas_package.raw_measurements_(1),
                  0,
                  0,
                  0;
            P_ = MatrixXd::identity(5, 5);
            P_(0, 0) = std_laspx_;
            P_(1, 1) = std_laspy_;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double rho_dot = meas_package.raw_measurements_(2);

            double p_x = rho * cos(phi);
            double p_y = rho * sin(phi);

            x_ << p_x, p_y, rho_dot, phi, 0;
        }
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
}

void UKF::Prediction(double delta_t) {
    /**
    * TODO: Complete this function! Estimate the object's location.
    * Modify the state vector, x_. Predict sigma points, the state,
    * and the state covariance matrix.
    */

    // Augment state and covariance matrix
    MatrixXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug(5) = std_a_;
    x_aug(6) = std_yawdd_;

    MatrixXd P_aug = MatrixXd::identity(n_aug_, n_aug_);
    P_aug.fill(0.0);
    for (int i = 0; i < n_x_; i++) {
        P_aug(i, i) = P_(i, i);
    }
    P_aug(5, 5) = std_a_*std_a;
    P_aug(6, 6) = std_yawdd_*std_yawdd_;

    // Generate Sigma Points
    VectorXd Xsig = VextorXd(n_aug_);
    MatrixXd A = P_aug.llt().matrixL();
    MatrixXd c1 = sqrt(lambda_ + n_x_);
    Xsig_pred_.col(0) = x_;
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Xsig(i+1) = x_aug + c1*A;
        Xsig(i+1+n_aug_) = x_aug - c1*A;
    }

    // Sigma points prediction
    double px = Xsig(0);
    double py = Xsig(1);
    double v = Xsig(2);
    double psi = Xsig(3);
    double psi_dot = Xsig(4);
    double mu_a = Xsig(5);
    double mu_psi_ddot = Xsig(6);

    VectorXd x1 = VectorXd(n_aug_);
    if (psi_dot < 0.0001) {
        x1 << v*cos(psi)*delta_t,
              v*sin(psi)*delta_t,
              0,
              psi_dot*delta_t,
              0;
    } else {
        x1 << V/psi_dot*(sin(psi+psi_dot*delta_t) - sin(psi)),
              V/psi_dot*(-cos(psi+psi_dot*delta_t) + cos(psi)),
              0,
              psi_dot*delta_t,
              0;
    }
    VectorXd x2 = VectorXd(n_aug_);
    x2 << 0.5*delta_t*delta_t*cos(psi)*mu_a,
          0.5*delta_t*delta_t*sin(psi)*mu_a,
          delta_t*mu_a,
          0.5*delta_t*delta_t*mu_psi_ddot,
          delta_t*mu_psi_ddot;

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        for (int j = 0; j < 5; j++) {
            Xsig_pred_(j, i) = px + x1(j) + x2(j);
        }
    }

    // Predict mean and covariance
    VectorXd x = VectorXd(n_x_);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        x += weights(i)*Xsig_pred_;
    }
    x_ = x;

    MatrixXd P = MatrixXd(n_x_, n_x_);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        P += weights(i)*(Xsig_pred_ - x)*(Xsig_pred_ - x).transpose();
    }
    P_ = P;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use lidar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the lidar NIS, if desired.
    */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use radar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the radar NIS, if desired.
    */
}

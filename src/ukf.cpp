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
    std_a_ = 2;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 1;

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

    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

    weights_ = VectorXd(2*n_aug_+1);
    weights_.fill(0.5/(lambda_+n_aug_));
    weights_(0) = lambda_/(lambda_+n_aug_);

    Rlidar_ = MatrixXd(2, 2);
    Rlidar_(0, 0) = std_laspx_*std_laspx_;
    Rlidar_(1, 1) = std_laspy_*std_laspy_;

    Rradar_ = MatrixXd(3, 3);
    Rradar_(0, 0) = std_radr_*std_radr_;
    Rradar_(1, 1) = std_radphi_*std_radphi_;
    Rradar_(2, 2) = std_radrd_*std_radrd_;
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

            P_.fill(0.0);
            P_(0, 0) = std_laspx_*std_laspx_;
            P_(1, 1) = std_laspy_*std_laspy_;
            P_(2, 2) = 1;
            P_(3, 3) = 1;
            P_(4, 4) = 1;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double rho_dot = meas_package.raw_measurements_(2);

            double p_x = rho * cos(phi);
            double p_y = rho * sin(phi);

            x_ << p_x, p_y, rho_dot, phi, 0;

            P_.fill(0.0);
            P_(0, 0) = std_radr_*std_radr_;
            P_(1, 1) = std_radphi_*std_radphi_;
            P_(2, 2) = std_radrd_*std_radrd_;
            P_(3, 3) = 1;
            P_(4, 4) = 1;
        }
        is_initialized_ = true;
        time_us_ = meas_package.timestamp_;
        return;
    }

    double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;

    Prediction(dt);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
        UpdateLidar(meas_package);
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
        UpdateRadar(meas_package);
}

void UKF::Prediction(double delta_t) {
    /**
    * TODO: Complete this function! Estimate the object's location.
    * Modify the state vector, x_. Predict sigma points, the state,
    * and the state covariance matrix.
    */

    // Augment state and covariance matrix
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_*std_a_;
    P_aug(6, 6) = std_yawdd_*std_yawdd_;

    // Generate Sigma Points
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    MatrixXd A = P_aug.llt().matrixL();
    double c1 = sqrt(lambda_ + n_aug_);
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i+1) = x_aug + c1*A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - c1*A.col(i);
    }

    // Sigma points prediction
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
        // extract values for better readability
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        } else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        // add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        // write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;
      }


    // Predict mean and covariance
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        x_ += weights_(i)*Xsig_pred_.col(i);
    }

    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        while (x_diff(3) > M_PI) x_diff(3) -= 2*M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2*M_PI;

        P_ += weights_(i)*x_diff*x_diff.transpose();
    }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use lidar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the lidar NIS, if desired.
    */

    int n_z = 2;
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0),
         meas_package.raw_measurements_(1);

    // Transfer predicted state into measurement state space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Zsig(0, i) = Xsig_pred_(0, i);
        Zsig(1, i) = Xsig_pred_(1, i);
    }

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred += weights_(i)*Zsig.col(i);
    }

    // Predict covatance matrix
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_(i)*z_diff*z_diff.transpose();
    }

    S += Rlidar_;

    // Cross-correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc += weights_(i)*x_diff*z_diff.transpose();
    }

    // Kalman gain
    MatrixXd K = Tc*S.inverse();

    // Update state
    x_ += K*(z - z_pred);
    P_ -= K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    * TODO: Complete this function! Use radar data to update the belief
    * about the object's position. Modify the state vector, x_, and
    * covariance, P_.
    * You can also calculate the radar NIS, if desired.
    */

    int n_z = 3;
    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_(0),
         meas_package.raw_measurements_(1),
         meas_package.raw_measurements_(2);

    // Transfer predicted state into measurement state space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
      // extract values for better readability
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v  = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);

      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;

      // measurement model
      Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                     // r
      Zsig(1, i) = atan2(p_y, p_x);                             // phi
      Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); // r_dot
    }

    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred += weights_(i)*Zsig.col(i);
    }

    // Predicted covariance matrix
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
        S += weights_(i)*z_diff*z_diff.transpose();
    }

    S += Rradar_;

    // Cross-correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;

        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while (x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;

        Tc += weights_(i)*x_diff*z_diff.transpose();
    }

    // Kalman gain
    MatrixXd K = Tc*S.inverse();

    // Update state
    VectorXd z_diff = z - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
    x_ += K*z_diff;
    P_ -= K*S*K.transpose();
}

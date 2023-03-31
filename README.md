# Unscented Kalman Filter

UKF highway project

In this project, I implemented an Unscented Kalman Filter to estimate the state of multiple cars on a highway using noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower than the tolerance outlined in the project rubric.

The main program can be built and ran by doing the following from the project top directory.
```
mkdir build
cd build
cmake ..
make
./ukf_highway
```

```main.cpp``` is using ```highway.h``` to create a straight 3 lanes highway environment with 3 traffic cars and the main ego car at the center. The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car has it's own UKF object generated for it, and will update each individual one during every time step.

The red spheres above cars represent the (x, y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.

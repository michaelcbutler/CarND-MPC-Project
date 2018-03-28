## Model Predictive Control Project

The goal of this project is to implement and configure an MPC algorithm to successfully guide the simulated car around a lap of the test track.

The implementation and configuration are largely based on the solution presented in lesson 20.9 "Solution: Mind the Line." My modifications are discussed in the sections below.

### N & dt

I moved the definition of `N` and `dt` to within the `MPC` class. This adds scoping to the parameters while allowing visibility within `main.cpp`.

```
class MPC {
 public:
  MPC();

  virtual ~MPC();

  // timestep length and duration
  static const size_t N = 15;
  static constexpr double dt = 0.1;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};
```
### Fitting a polynomial to the waypoints

The waypoints are first transformed into the vehicle's frame of reference:
```
// transform (x, y) points to vehicle coordinate system 
for (int i = 0; i < ptsx.size(); ++i) {
  double dx = ptsx[i] - px;
  double dy = ptsy[i] - py;
  ptsx[i] = dx*cos(psi) + dy*sin(psi);
  ptsy[i] = dy*cos(psi) - dx*sin(psi);
}
```
Then the points are converted to the type `Eigen::VectorXd` required by the fitting algorithm:
```
Eigen::VectorXd xvals = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptsx.data(), ptsx.size());
Eigen::VectorXd yvals = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptsy.data(), ptsy.size());
```
And coefficents for a cubic polynomial fit are computed:
```
auto coeffs = polyfit(xvals, yvals, 3);
```

### Cost function

The components of the cost function duplicate the components from lesson 20.9. However, different weighting coefficients were considered:
```
// Minimize CTE, yaw error, velocity error
for (int t = 0; t < MPC::N; t++) {
  fg[0] += wt_cte * CppAD::pow(vars[cte_start + t] - ref_cte, 2);
  fg[0] += wt_epsi * CppAD::pow(vars[epsi_start + t] - ref_epsi, 2);
  fg[0] += wt_v * CppAD::pow(vars[v_start + t] - ref_v, 2);
}

// Minimize the magnitude of steering and throttle inputs
for (int t = 0; t < MPC::N - 1; t++) {
  fg[0] += wt_delta * CppAD::pow(vars[delta_start + t], 2);
  fg[0] += wt_a * CppAD::pow(vars[a_start + t], 2);
}

// Minimize the change of (smooth) steering and throttle inputs
for (int t = 0; t < MPC::N - 2; t++) {
  fg[0] += wt_delta_dot * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
  fg[0] += wt_a_dot * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
}
```

### Initial constraints

We initialize the model to the initial state as lesson 20.9.
```
fg[1 + x_start] = vars[x_start];
fg[1 + y_start] = vars[y_start];
fg[1 + psi_start] = vars[psi_start];
fg[1 + v_start] = vars[v_start];
fg[1 + cte_start] = vars[cte_start];
fg[1 + epsi_start] = vars[epsi_start];
```

### Model constraints

The model constraints differ from lesson 20.9 in three respects. First, the loop indexing begins with `t = 0` rather than `t = 1`. This gives a more natural mapping of the variables at `t` or `t + 1`. Second, the sign of delta is negated for consistency with the simulation code. Third, the computations of `f0` and `psides0` are modified for a cubic polynomial fit.
```
for (int t = 0; t < MPC::N - 1; t++) {
  // The state at time t+1 .
  AD<double> x1 = vars[x_start + t + 1];
  AD<double> y1 = vars[y_start + t + 1];
  AD<double> psi1 = vars[psi_start + t + 1];
  AD<double> v1 = vars[v_start + t + 1];
  AD<double> cte1 = vars[cte_start + t + 1];
  AD<double> epsi1 = vars[epsi_start + t + 1];

  // The state at time t.
  AD<double> x0 = vars[x_start + t];
  AD<double> y0 = vars[y_start + t];
  AD<double> psi0 = vars[psi_start + t];
  AD<double> v0 = vars[v_start + t];
  AD<double> cte0 = vars[cte_start + t];
  AD<double> epsi0 = vars[epsi_start + t];

  // Only consider the actuation at time t.
  AD<double> delta0 = vars[delta_start + t];
  AD<double> a0 = vars[a_start + t];

  AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
  AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0);

  // Here's `x` to get you started.
  // The idea here is to constraint this value to be 0.
  //
  // Recall the equations for the model:
  // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
  // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
  // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
  // v_[t] = v[t-1] + a[t-1] * dt
  // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
  // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
  fg[2 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * MPC::dt);
  fg[2 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * MPC::dt);
  fg[2 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * MPC::dt); // delta sign issue correction
  fg[2 + v_start + t] = v1 - (v0 + a0 * MPC::dt);
  fg[2 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * MPC::dt));
  fg[2 + epsi_start + t] = epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * MPC::dt); // delta sign issue correction
}
```

### Initialization

For the first `mpc.Solve()`, the `state` vector is initialized as:
```
state << 0.0, 0.0, 0.0, v, cte, epsi;
```
where initial `x`, `y`, and `psi` are zero. `cte` and `epsi` evaluate at x = 0 as:
```
double cte = polyeval(coeffs, 0.0);
```
and 
```
double epsi = -atan(coeffs[1]); 
```
### Parameter tuning and latency adjustments

The initial parameter settings were taken from lesson 20.10 "Tuning MPC" with `N` = 25, `dt` - 0.05, `wt_delta` = 100, and `wt_delta_dot` = 500; All other weights were set to unity. With the latency set to 100 msec and reference velocity set to 40, the car would start well but swerve wildly at the left turn after the bridge. Removing the latency postponed the wild swerving by one corner.

Reducing `N` to 15 and increasing `dt` to 0.1 allowed the car to successfully complete the lap. However, at the previously troublesome turns the MPC trajectory would briefly "freak out." 

To force the MPC path to more tightly track the way points, I increased `wt_cte` to 2000 and `wt_epsi` to 4000. Again, the simulation was able to complete the lap, but this time without the "glitches."

But when I re-applied the 100 msec latency, the car was unstable again. I addressed the latency issue with two adjustments. First, I adjusted the vehicle position accounting for latency before transforming the way points into the vehicle coordinate system:
```
  // adjust px, py, psi for latency (acceleration is unknown)
  const double latency = 0.1; // 100 msec
  px += v*cos(psi)*latency;
  py += v*sin(psi)*latency;
  psi -= steer_value*v/MPC::Lf*latency;
```
I ignored acceleration because the relationship to throttle position is ambiguous and since the simulated car moves at near-constant speed I assume acceleration is negligible.

The second adjustment increased `dt` to 0.15 (greater than the latency). I reduced `N` to 10 to keep the total time constant.

With `ref_v` set to 40, the car would successfully complete the lap at approximately 30 mph. With `ref_v` increased to 75, the car would successfully complete the lap at approximately 55 mph.

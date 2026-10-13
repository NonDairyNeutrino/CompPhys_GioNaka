# Exercise 1.1: Constant Acceleration

## Introduction

Consider an object positioned some distance above Earth's surface such that its acceleration due to gravity _g_ is approximately 9.8 m/s^2 toward the Earth's surface.
If the object is released from rest _v\_0 = 0 m/s_, and allowed to fall for _t\_f = 10_ seconds without the influence of air resistance, our model of gravity tells us that the object should continue to accelerate and increase its speed toward Earth for the entire duration.
Becuase the acceleration is constant, the motion of the object can be found analytically via kinematic _x(t) = x\_0 + v\_0 t + 1/2 a t^2_.
Though, this motion can also be simulated by applying the Euler method to the equation of motion: _d^2 y / dt^2 = -g_.

## Methods: The Euler Method

As it is one of the physicist's favorite tools, we consider the Taylor expansion centered at zero _y(t) = y(0) + dy/dt t + 1/2 d^2 y / dt^2 t^2 + ..._.
If we assume the time step _dt_ is small enough, we consider it sufficiently accurate to only consider the terms up to first order in time as the higher powers of the time step are increasingly negligible; this leaves the approximation: _y(t) = y(0) + dy/dt t_.
If the time is discretized, then this becomes _the Euler Method_: _y\_n = y\_(n - 1) + dy/dt dt_, where _y\_n = y(n dt)_.
While this form comes with simplicity, the error per step (i.e. the _local_ error) comes from where we truncate the Taylor expansion, namely dt^2, while the total or _global_ error comes from the product of the number of time steps _t/dt_ and the local error to yield _dt_; thus the error between the simulated and true behavior will e.g. halve for each halving of the time step.

## Results: Linearly Increasing Velocity



## Discussion: Equivalence of Numeric and Analytic solutions



## Conclusion

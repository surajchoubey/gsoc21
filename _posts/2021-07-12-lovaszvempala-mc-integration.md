---
layout: post
title:  "Lovasz-Vempala Integration"
date:   2021-07-13 00:00:00 +0530
categories: jekyll update
permalink: /lv-mc-integration/
---

Before going through this, it is strongly recommended you have the knowledge of Monte Carlo Integration process. You can study it from the [wikipedia page](https://en.wikipedia.org/wiki/Monte_Carlo_integration#Overview). Also, it is recommended that you have gone through the previous post on [Simple MC Integration]({{site.baseurl}}/simple-mc-integration/).

## What is Lovasz Vempala Integration?
Lovasz Vempala Integration algorithm is used from page 7 of this [research paper](https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf). This integration algorithm can be used to integrate logconcave density functions particularly $$g(x)$$, which are expressed in the form of $$f(x) = e^{-g(x)}$$. The logconcave function $$f(x)$$ can be integrated over a polytope using [**Hamiltonian Monte Carlo Sampling**](https://arxiv.org/abs/1701.02434). The integration algorithm uses this algorithm to run the random walks used for sampling for the Monte Carlo integration process.

## Why this Integration Algorithm over the Simple Monte Carlo Integration?
Assuming that the have the knowledge of Simple Monte Carlo Integration, for integration of functions over $$n$$ dimensional polytopes, from the previous post on [Simple Monte Carlo Integration]({{site.baseurl}}/simple-mc-integration/).

In Simple MC Integration, the accuracy decreases as the integral answers deviate proportionately as we increase the dimensions and causes more and more estimation errors. Ideally, more the number of random sample points taken into consideration from the convex body it does help in mixing the random points(which are taken into account to perform the Monte Carlo Integration Algorithm) more uniformly througout the convex body. 

However, it cannot be trusted for anticipating very small estimation error. If the user prefers to go using trustable sampling technique and expect a very small estimation error, Lovasz Vempala Integration algorithm has an edge over Simple MC Integration.

## The Theory Behind Lovasz Vempala Integration Algorithm
The algorithm can be little tricky to explain but I am going to make it easy and simple to understand. The algorithm can be found on page 7 of this [research paper](https://www.cc.gatech.edu/~vempala/acg/www/www/papers/integration.pdf). Now, that we are aware of the caveats of Simple Monte Carlo Integration process, lets jump onto this.

Let $$f(x) = e^{-g(x)}$$ be desired logconcave function to integrate over a convex body. Some variables you need to know about: <br>
* $$\beta$$ is used-defined such that $$0 < \beta < 1$$.
* $$n$$ is the dimension of convex body ($$n > 0$$).
* $$\varepsilon$$ is volume and variable calculation error set by user.
<br><br>
* $$
\begin{equation*}
  B = 2n + 2 \cdot ln \left( \dfrac{1}{\varepsilon}\right) + n\cdot ln\left(\dfrac{1}{\beta}\right)
\end{equation*}
$$
<br><br>
* $$
  \begin{equation*}
    \displaystyle m = \Bigl\lceil \sqrt{n} \cdot ln(B) \Bigr\rceil
  \end{equation*}
$$
<br><br>
* $$
  \begin{equation*}
    \displaystyle k = \Bigl\lceil \dfrac{512}{\varepsilon^{2}} \cdot \sqrt{n} \cdot ln(B) \Bigr\rceil
  \end{equation*}
$$
<br><br>
* $$
  \begin{equation*}
    \displaystyle a_{i} = \dfrac{1}{B}\left(1 + \dfrac{1}{\sqrt{n}} \right)^{i}
  \end{equation*}
$$
<br><br>
* $$
  \begin{equation*}
    \displaystyle f_{i}(x) = f(x)^{a_{i}}
  \end{equation*}
$$
<br><br>
* $$
  \begin{equation*}
    f_{m}(x) = f(x)
  \end{equation*}
$$
<br><br>

#### Algorithm step-by-step:

1. We run warmstart samples from a choosen point say $$x_{0}$$ for $$k$$ times to ensure proper mixing using uniform random walks mechanism around the $$n$$ dimensional convex body $$K$$ ($$n > 0$$). The point is chosen such than $$f(x_{0}) \ge \beta^{n} \cdot f(x_{max})$$, where $$x_{max}$$ is the global maxima of the function within the convex body $$K$$ choosen itself and $$\beta$$ is a parameter).

2. Let $$W_{0} = volume(K)$$ be the estimation of the volume of convex body $$K$$.

3. For $$i = 1,2,...,m$$, do the following: <br>
    * Run the samples $$k$$ times with target density proportional to $$f_{i-1}$$ and starting points $$ \displaystyle X_{i-1}^{1}, X_{i-1}^{2}, ... , X_{i-1}^{k} $$ to get independent random points $$ \displaystyle X_{i}^{1}, X_{i}^{2}, ... , X_{i}^{k} $$

    * Using these points we compute <br>

    {:refdef: style="text-align: center;"}
    $$
        \begin{equation*} 
          \displaystyle W_{i} = \dfrac{1}{k} {\sum}\limits_{j = 1}^{k} f(X_{i}^{j})^{a_{i} - a_{i-1}}
        \end{equation*}
    $$
    {: refdef}

4. Return $$W = W_{0} W_{1} ... W_{m}$$.

$$W$$ is the desired integral value for our logconcave function $$f(x)$$ over the desired subspace $$K$$. By multiplying $$W_{1} W_{2}... W_{m}$$ the telescopic series in the exponents gets canceled out with the subsequent terms. And we are left with the desired value of the the required integral $$W$$.

## The Project
As discussed the theory above in this process, it involves to build an integration function ```lovasz_vempala_integrate()``` using Lovasz Vempala Integration algorithm. Regarding building the function we take in several parameters. Listing them down, and explaining each of them before proceeding to the next part:

### Structure of the integration function<br>
```
template
<
	typename EvaluationFunctor,
	typename GradientFunctor,
	typename Parameters,
	typename WalkType,
	typename Polytope,
	typename Point,
	typename NT
>
NT lovasz_vempala_integrate(EvaluationFunctor &g,
                            GradientFunctor &grad_g,
                            Parameters &params,
                            Polytope &P,
                            Point x0,
                            NT beta = 1.0,
                            volumetype voltype = SOB,
                            unsigned int walk_length = 10,
                            NT epsilon = 0.1)
{
  // code goes here
}                            
```
<br>
* `EvaluationFunctor` is the type of function expressed in terms of $$g(x)$$ which is meant to be integrated in the form of $$f(x) = e^{-g(x)}$$ around the provided convex body.

* `GradientFunctor` is a type of gradient function which returns the gradient of a the `EvaluationFunctor g` as discussed above to return<br>
  {:refdef: style="text-align: center;"}
  $$
    \begin{equation*}
      \left[ \dfrac{\partial f}{\partial x_{1}}, \dfrac{\partial f}{\partial x_{2}}, ... , \dfrac{\partial f}{\partial x_{n}} \right]
    \end{equation*}
  $$
  {:refdef }
  the term right above is the gradient for $$f(x)$$ which is $$\nabla f(x)$$. Also, $$\dfrac{\partial f}{\partial y}$$ represents partial differentiation of a function $$f$$ with respect to $$y$$.

* `Parameters` are the additional variables that let you decide more variables to customize your function and your gradient function. (Ideally you can look here in [oracle_functors.hpp](https://github.com/GeomScale/volume_approximation/blob/develop/include/ode_solvers/oracle_functors.hpp) from [GeomScale/volume_approximation](https://github.com/GeomScale/volume_approximation) repository and have a close look o how gradient functor, evaluation functor and parameters are constructed under one `struct`)

* `Point` is a user-defined datatype to store $$n$$-dimensional points in the $$n$$-dimensional space. $$x_{0}$$ is a point chosen such that $$f(x_{0}) \ge \beta^{n} \cdot f(x_{max})$$ satisfies.

* `NT beta` is a parameter in integration to help decide other $$\beta$$-dependent parameters.

* `volumetype` is an `enum` used to specify the volume algorithm in this [directory](https://github.com/GeomScale/volume_approximation/tree/develop/include/volume) from [GeomScale/volume_approximation](https://github.com/GeomScale/volume_approximation) to calculate the volume of the given convex body. Available options from `enum` ones are `CB`,`CG` and `SOB`. Namely, cooling balls, cooling gaussians and sequence of balls algorithm.

* `walk_length` is the walk length which is going to be used for length of the walks which is supposed to be used in running warmstart samples, volume calculation algorithms and HMC algorithm.

* `NT epsilon` i.e. $$\varepsilon$$ is the permissible error which can be set by the user to increase accuracy for the volume algorithm. $$\varepsilon$$ is also dependent on the $$\beta$$ dependent parameters too.

* `typename WalkType` is type of random walk which will be used to run hit-and-run sampling for the warmstart samples and ensure proper mixing. Available ones are **BallWalk**, **CDHRWalk**, **RDHRWalk**, **BilliardWalk**, **AcceleratedBilliardWalk** and so on.

### How this algorithm is applied to the project?
After supplying the above quantities and objects to the integration function, we get down to the main business.

1. We use an `OptimizationFunctor` which wraps the supplied `EvaluationFunctor` and `GradientFunctor` and samples from the `EvaluationFunctor` proportionally to the variance $$a_{i}$$. `OptimizationFunctor` has its own parameters which helps us setup the variance accordingly.

2. Hamiltonian Monte Carlo Walk which is already implemented in [here](https://github.com/GeomScale/volume_approximation/blob/develop/include/random_walks/hamiltonian_monte_carlo_walk.hpp) from [GeomScale/volume_approximation](https://github.com/GeomScale/volume_approximation) is used to sample from logconcave density functions, which is $$f(x) = e^{-g(x)}$$ here in context.

3. From the supplied `WalkType`, $$k$$ warmstart samples are run to ensure that the HMC walks got properly mixed and follow memoryless property to choose the starting point based on random decision.

4. Proceeding to the main part of the Lovasz Vempala Algorithm, a 2-level nested loop is run. In the outer loop, we keep running the loop until variance $$a_{i} $$ crosses unity. $$W_{i}$$ for $$i = 1,2,...,m$$ is calculated and multiplied subsequently to $$W_{0}$$ which was previous calculated as $$volume(K)$$.

5. For the inner loop HMC walks are executed inside the convex body which samples points properly. $$W_{i}$$ is calculated for $$g(x)$$ expressed as a logconcave function where $$g(x) \propto a_{i} - a_{i-1}$$ for $$k$$ Hamiltonian Monte Carlo walks with `walk_length` set from the function declaration itself. 

6. Return $$W = W_{0} W_{1} ... W_{m}$$.

## Testing 
Testing has been done using **Monte Carlo Integration** functions from [torchquad](https://www.esa.int/gsp/ACT/open_source/torchquad/). You can find the tests on this [github repository](https://github.com/surajchoubey/torchquad-checks). You can go through the `README.md` and it will explain all the required information. 

## Usage
For testing out the `lovasz_vempala_integrate()` it is recommended that you have the knowledge of [oracle_functors.hpp](https://github.com/GeomScale/volume_approximation/blob/develop/include/ode_solvers/oracle_functors.hpp)  in [GeomScale/volume_approximation](https://github.com/GeomScale/volume_approximation) github repository to get basic understanding of how `EvaluationFunctor` and `GradientFunctor` are implemented. 

### Example Code
As explained above, the variables are passed so that they satisfy the conditions of the algorithm. For a more detailed look on the code you can visit to [lv_mc_integration.cpp](https://github.com/surajchoubey/volume_approximation/blob/lovasz-vempala/test/lv_mc_integration.cpp) in this PR.

```ruby


//# 1. g is the EvaluationFunctor
//# 2. grad_g is the GradientFunctor
//# 3. HP is the polytope passed
//# 4. AcceleratedBilliardWalk is the walktype passed for running warmstart samples
//# 5. HPOLYTOPE specifies notation of the polytope( also other examples 
//#    like VPOLYTOPE not tested particularly on it yet)
//# 6. Point is a user defined datatype to represent n-dimensional points in n-dimensional    space
//# 7. NT is the number type(it can be float/double)
//# 8. beta is an integration parameters whose value can be between (0,1]
//# 8. x0 is the point that is the initial point chosen to run the 
//#    warmstart samples so it properly mixes before running the HMC walks
//#    It should be chosen such that it satisfies f(x0) >= beta ^ n * max(f)

HPOLYTOPE HP;
NT integral value;

HP = generate_cube <HPOLYTOPE> (5, false);
integral_value = lovasz_vempala_integrate 
  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
  (g, grad_g, HP, x0, beta, SOB, 5, 0.1);

std::cout << "Integral value for a 5D polytope = " << integral_value << std::endl;

HP = generate_cube <HPOLYTOPE> (10, false);
integral_value = lovasz_vempala_integrate 
  <EvaluationFunctor, GradientFunctor, AcceleratedBilliardWalk, HPOLYTOPE, Point, NT>
  (g, grad_g, HP, x0, beta, SOB, 5, 0.1);

std::cout << "Integral value for a 10D polytope = " << integral_value << std::endl;

```

### Output Code

```
Integral value for a 5D polytope = 7.44508
Integral value for a 10D polytope = 55.7236

```


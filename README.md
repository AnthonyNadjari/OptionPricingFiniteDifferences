
Finite Difference Option Pricer
# PDE Option Pricing - Black Scholes Model  - using Finite Difference method

The python file allows the user to compute the price of any payoff using the Finite Difference method. There are two main purposes: <br>

* **Via the Python File**: the user shall create a class (using the main parameters of an option), that contains the desired payoff. For instance, to create a payoff vvector, user should first define the boundary conditions vector, pass it thourgh the exponential function and finally, use this vector as the argument of a payoff function.

* **Main**: compute the price of a call option and its greeks via finite difference (PDE) and via closed form with the test_mesh function. The user is more than welcome to change the inputs of the function (in the *main()*). It takes approximately 20 seconds with 3 months maturity to get the price and the greeks. The user can uncomment rows 48-49 to see the mesh. <br> <br>
*Remark: The gamma computation seems to have an issue as we end up with a factor 2 between our value of the gamma and the one of the closed form gamma. We were not able to find the source of the issue.* 

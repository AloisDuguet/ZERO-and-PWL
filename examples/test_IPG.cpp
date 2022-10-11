#include <zero.h>


int main(int argc, char **argv) {

  GRBEnv GurobiEnv;
  try {

     Models::IPG::IPGInstance IPG_Instance;
     int                      numItems   = 2;
     int                      numPlayers = 2;

     // Objective terms
     arma::vec    c(numItems);
     arma::sp_mat C(numPlayers * (numItems - 1), numItems);
     // Constraints
     arma::sp_mat a(1, numItems);
     arma::vec    b(1);
     // The index of the integer variables
     arma::vec IntegerIndexes(numItems);
     // Implicit bounds on variables. Could be omitted
     VariableBounds VarBounds = {{0, 1}, {0, 1}};

     for (unsigned int i = 0; i < numItems; ++i)
        IntegerIndexes.at(i) = i;

     C(0, 0) = 2;
     C(1, 1) = 3;

     a(0, 0) = 3;
     a(0, 1) = 4;
     b(0)    = 5;

     // The standard is minimization.
     c(0) = -1;
     c(1) = -2;
     std::cout << "Ciao Khalid \n";
     // Create a parametrized Integer Program
     MathOpt::IP_Param PlayerOne(C, a, b, c, IntegerIndexes, VarBounds, &GurobiEnv);

     // Let's create another parametrized Integer Program for the second player.
     // We'll reuse the previously created matrices and vectors

     C(0, 0) = 5;
     C(1, 1) = 4;

     a(0, 0) = 2;
     a(0, 1) = 5;

     c(0) = -3;
     c(1) = -5;

     MathOpt::IP_Param PlayerTwo(C, a, b, c, IntegerIndexes, VarBounds, &GurobiEnv);


     // Add the players to the instance. We also specify a file path to write the instance
     IPG_Instance.addIPParam(PlayerOne, "A_Parametrized_KnapsackProblem1");
     IPG_Instance.addIPParam(PlayerTwo, "A_Parametrized_KnapsackProblem2");
     // Save the instance with the standardize format
     IPG_Instance.save("A_Knapsack_Game");
     // Create a model from the instance
     Models::IPG::IPG IPG_Model(&GurobiEnv, IPG_Instance);
     // Select the algorithm to compute a Nash Equilibrium
     IPG_Model.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
     // Extra parameters
     IPG_Model.setDeviationTolerance(3e-4);
     IPG_Model.setNumThreads(4);
     IPG_Model.setLCPAlgorithm(Data::LCP::Algorithms::MIP);
     IPG_Model.setGameObjective(Data::IPG::Objectives::Feasibility);
     IPG_Model.setTimeLimit(600);
     // Lock the model
     IPG_Model.finalize();
     // Run!
     IPG_Model.findNashEq();

     // Print the solution

     IPG_Model.getX().at(0).print("Player 1:");
     IPG_Model.getX().at(1).print("\n Player 2:");

  } catch (ZEROException &e) {
     throw ZEROException(e);
  }
}

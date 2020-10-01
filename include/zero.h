/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#pragma once

#include <armadillo>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "support/codes.h"
#include "support/utils.h"
#include "support/version.h"

using perps          = std::vector<std::pair<unsigned int, unsigned int>>;
using DoubleAttrPair = std::vector<std::pair<unsigned int, double>>;
std::ostream &                            operator<<(std::ostream &ost, perps C);
template <class T> std::ostream &         operator<<(std::ostream &ost, std::vector<T> v);
template <class T, class S> std::ostream &operator<<(std::ostream &ost, std::pair<T, S> p);
using spmat_Vec = std::vector<std::unique_ptr<arma::sp_mat>>;
using vec_Vec   = std::vector<std::unique_ptr<arma::vec>>;
template <typename T1, typename T2, typename T3> using triple = std::tuple<T1, T2, T3>;

// Forward declarations
class ZEROException;
enum class ZEROErrorCode;

namespace MathOpt {
  /**
	* @brief This namespace contains the definition of the support structures for the mathematical
	*optimization
	**/
  struct QP_Objective;
  struct QP_Constraints;
  class MP_Param;
  class QP_Param;
  class IP_Param;
  class LCP;
} // namespace MathOpt

namespace Game {
  /**
	* @brief This namespace contains the definitions of the games
	**/
  class NashGame;
  class PolyLCP;
  class OuterLCP;
  class EPEC;
  class IPG;
} // namespace Game
namespace Algorithms {
  /**
	* @brief This namespace contains the definition of the algorithms. For each
	* game, there is a lower-level namespace under which the algorithms are
	* nested.
	*/
  class AbstractAlgorithm; ///< Abstact type for other algorithms
  namespace EPEC {
	 class Algorithm;
	 class PolyBase;
	 class FullEnumeration;
	 class InnerApproximation;
	 class CombinatorialPNE;
	 class OuterApproximation;
  } // namespace EPEC
  namespace IPG {
	 class Algorithm;
	 class Oracle;
  } // namespace IPG
} // namespace Algorithms

class ZEROAlgorithmData;


/**
 * @brief This namespace contains the Data support structures for the
 * algorithms. For any game for which an algorithm is available, the class DataObject contains the
 * control parameters. Other objects (enum classes, etc) are control parameters for their relative
 * objects.
 */
namespace Data {

  namespace EPEC {
	 class DataObject;
  }
  namespace IPG {
	 class DataObject;
  }
  namespace LCP {
	 enum class PolyhedraStrategy;
  }

} // namespace Data


#include "games/games.h"
#include "mathopt/mathopt.h"

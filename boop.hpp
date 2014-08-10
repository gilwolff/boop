/*
 * boop = bond orientational order parameters
 * 
 * This library computes Steinhardt's Ql and Wl order parameters.
 * 
 * The library depends on boost's special functions
 * 
 */

#ifndef BOOP_HPP
#define BOOP_HPP

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <iterator>

#ifndef PI
#define PI (4*atan(1))
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace boop{
#include "periodic.hpp"
#include "particle.hpp"
#include "bond.hpp"
#include "Q.hpp"
#include "box.hpp"
}
#endif

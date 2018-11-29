//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/checks.h"

#if !defined(KRATOS_STATIC_CHECKS_H_INCLUDED )
#define  KRATOS_STATIC_CHECKS_H_INCLUDED

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Static functions
///@{

static inline void KRATOS_STATIC_CHECK(const bool Check)
{
    KRATOS_CHECK(Check);
}

static inline void KRATOS_STATIC_CHECK_IS_FALSE(const bool Check)
{
    KRATOS_CHECK_IS_FALSE(Check);
}

template<class TTypea, class TTypeb>
static inline void KRATOS_STATIC_CHECK_EQUAL(const TTypea& a, const TTypeb& b)
{
    KRATOS_CHECK_EQUAL(a, b);
}

template<class TTypea, class TTypeb>
static inline void KRATOS_STATIC_CHECK_NOT_EQUAL(const TTypea& a, const TTypeb& b)
{
    KRATOS_CHECK_NOT_EQUAL(a, b);
}

static inline void KRATOS_STATIC_CHECK_LESS(const double a, const double b)
{
    KRATOS_CHECK_LESS(a, b);
}

static inline void KRATOS_STATIC_CHECK_LESS_EQUAL(const double a, const double b)
{
    KRATOS_CHECK_LESS_EQUAL(a, b);
}

static inline void KRATOS_STATIC_CHECK_GREATER(const double a, const double b)
{
    KRATOS_CHECK_GREATER(a, b);
}

static inline void KRATOS_STATIC_CHECK_GREATER_EQUAL(const double a, const double b)
{
    KRATOS_CHECK_GREATER_EQUAL(a, b);
}

static inline void KRATOS_STATIC_CHECK_STRING_CONTAIN_SUB_STRING(const std::string& a, const std::string& b)
{
    KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(a, b);
}

static inline void KRATOS_STATIC_CHECK_NEAR(const double a, const double b, const double tolerance)
{
    KRATOS_CHECK_NEAR(a, b, tolerance);
}

static inline void KRATOS_STATIC_CHECK_DOUBLE_EQUAL(const double a, const double b)
{
    KRATOS_CHECK_DOUBLE_EQUAL(a, b);
}
} /// Namespace Kratos

///@}

///@} addtogroup block

#endif // KRATOS_STATIC_CHECKS_H_INCLUDED  defined
